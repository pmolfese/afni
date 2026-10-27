use crate::model::{Bounds, Content, Figure, Plot};
use anyhow::{Context, Result, bail, ensure};
use resvg::{tiny_skia, usvg};
use std::{fmt::Write as _, io::Write, path::Path};

pub const FONT: &[u8] = include_bytes!("../assets/Ubuntu-Light.ttf");
const INK: &str = "#213448";
const MUTED: &str = "#70818f";
const GRID: &str = "#e8eef1";
#[derive(Clone, Copy, Debug)]
pub struct Area {
    pub x: f64,
    pub y: f64,
    pub w: f64,
    pub h: f64,
    pub bounds: Bounds,
}
#[cfg_attr(not(feature = "gui"), allow(dead_code))]
impl Area {
    pub fn point(&self, x: f64, y: f64) -> [f64; 2] {
        [
            self.x + (x - self.bounds.x[0]) / (self.bounds.x[1] - self.bounds.x[0]) * self.w,
            self.y + self.h
                - (y - self.bounds.y[0]) / (self.bounds.y[1] - self.bounds.y[0]) * self.h,
        ]
    }
    pub fn data(&self, x: f64, y: f64) -> [f64; 2] {
        [
            self.bounds.x[0] + (x - self.x) / self.w * (self.bounds.x[1] - self.bounds.x[0]),
            self.bounds.y[1] - (y - self.y) / self.h * (self.bounds.y[1] - self.bounds.y[0]),
        ]
    }
    pub fn contains(&self, x: f64, y: f64) -> bool {
        x >= self.x && x <= self.x + self.w && y >= self.y && y <= self.y + self.h
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TextTarget {
    FigureTitle,
    FigureCaption,
    PanelTitle(usize),
    PanelCaption(usize),
    XLabel(usize),
    YLabel(usize),
}
#[cfg_attr(not(feature = "gui"), allow(dead_code))]
#[derive(Clone, Copy, Debug)]
pub struct TextRegion {
    pub target: TextTarget,
    pub rect: [f64; 4],
}
#[cfg_attr(not(feature = "gui"), allow(dead_code))]
pub struct Scene {
    pub tree: usvg::Tree,
    pub areas: Vec<Area>,
    pub text_regions: Vec<TextRegion>,
    pub logical_height: f64,
}
pub struct Renderer {
    options: usvg::Options<'static>,
}
impl Renderer {
    pub fn new() -> Self {
        let mut options = usvg::Options::default();
        options.fontdb_mut().load_font_data(FONT.to_vec());
        options.font_family = "Ubuntu".into();
        Self { options }
    }
    pub fn scene(&self, fig: &Figure, views: &[Bounds], width: u32, height: u32) -> Result<Scene> {
        ensure!(
            !fig.plots.is_empty(),
            "add a graph before exporting an image"
        );
        ensure!(
            views.len() == fig.plots.len(),
            "view count does not match panels"
        );
        ensure!(views.iter().all(|b| b.valid()), "invalid plot view");
        let h = 1200. * height as f64 / width as f64;
        if fig.plots.len() == 1 && fig.plots[0].qc.margin_off {
            let a = Area {
                x: 0.,
                y: 0.,
                w: 1200.,
                h,
                bounds: views[0],
            };
            let mut svg = format!(
                r#"<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 1200 {h}">"#
            );
            draw_plot(&mut svg, &fig.plots[0], &a, 0);
            svg.push_str("</svg>");
            return Ok(Scene {
                tree: usvg::Tree::from_str(&svg, &self.options)?,
                areas: vec![a],
                text_regions: vec![],
                logical_height: h,
            });
        }
        let columns = fig.columns;
        let rows = fig.plots.len().div_ceil(columns);
        let gap = 20.;
        let margin = 28.;
        let figure_title = wrap(&fig.title, 1144., 30.);
        let title_extra = figure_title.len().saturating_sub(1) as f64 * 36.;
        let figure_caption = wrap(&fig.caption, 1144., 13.);
        let caption_height = if fig.caption.is_empty() {
            0.
        } else {
            figure_caption.len() as f64 * 18. + 15.
        };
        let top = 108. + title_extra;
        let bottom = 37. + caption_height;
        let pw = (1200. - margin * 2. - gap * (columns - 1) as f64) / columns as f64;
        let ph = (h - top - bottom - gap * (rows - 1) as f64) / rows as f64;
        ensure!(
            ph >= 210.,
            "figure is too short for {} panels; increase -height or reduce -width",
            fig.plots.len()
        );
        let mut svg = format!(
            r##"<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 1200 {h}"><rect width="1200" height="{h}" fill="#f3f6f7"/>"##
        );
        text(
            &mut svg,
            margin,
            31.,
            11.,
            "#147d92",
            &fig.subtitle,
            "start",
        );
        for (line, title) in figure_title.iter().enumerate() {
            text(
                &mut svg,
                margin,
                70. + line as f64 * 36.,
                30.,
                INK,
                title,
                "start",
            );
        }
        let mut text_regions = vec![
            TextRegion {
                target: TextTarget::FigureTitle,
                rect: [margin, 38., 1144., 38. + title_extra],
            },
            TextRegion {
                target: TextTarget::FigureCaption,
                rect: [margin, h - bottom + 4., 1144., caption_height.max(18.)],
            },
        ];
        if !fig.caption.is_empty() {
            for (line, caption) in figure_caption.iter().enumerate() {
                text(
                    &mut svg,
                    margin,
                    h - bottom + 17. + line as f64 * 18.,
                    13.,
                    INK,
                    caption,
                    "start",
                );
            }
        }
        let mut areas = vec![];
        for (i, (p, b)) in fig.plots.iter().zip(views).enumerate() {
            let x = margin + (i % columns) as f64 * (pw + gap);
            let y = top + (i / columns) as f64 * (ph + gap);
            rect(&mut svg, x, y, pw, ph, "#ffffff", Some("#dee7eb"), 9.);
            text(
                &mut svg,
                x + 20.,
                y + 26.,
                10.,
                "#147d92",
                &format!("{:02}  /  {}", i + 1, p.kind()),
                "start",
            );
            let titles = wrap(&p.title, pw - 40., 19.);
            let title_extra = titles.len().saturating_sub(1) as f64 * 24.;
            for (line, title) in titles.iter().enumerate() {
                text(
                    &mut svg,
                    x + 20.,
                    y + 52. + line as f64 * 24.,
                    19.,
                    INK,
                    title,
                    "start",
                );
            }
            text(
                &mut svg,
                x + 20.,
                y + 73. + title_extra,
                11.,
                MUTED,
                &short(&p.note, ((pw - 40.) / 5.4) as usize),
                "start",
            );
            let captions = wrap(&p.caption, pw - 40., 13.);
            let caption_height = if p.caption.is_empty() {
                0.
            } else {
                captions.len() as f64 * 18. + 12.
            };
            text_regions.push(TextRegion {
                target: TextTarget::PanelTitle(i),
                rect: [x + 20., y + 33., pw - 40., 24. + title_extra],
            });
            if !p.caption.is_empty() {
                text_regions.push(TextRegion {
                    target: TextTarget::PanelCaption(i),
                    rect: [
                        x + 20.,
                        y + ph - caption_height.max(18.) - 5.,
                        pw - 40.,
                        caption_height.max(18.),
                    ],
                });
            }
            if !p.caption.is_empty() {
                for (line, caption) in captions.iter().enumerate() {
                    text(
                        &mut svg,
                        x + 20.,
                        y + ph - caption_height + 12. + line as f64 * 18.,
                        13.,
                        INK,
                        caption,
                        "start",
                    );
                }
            }
            let heat = matches!(p.content, Content::Heat { .. });
            let mut a = Area {
                x: x + if heat { 93. } else { 73. },
                y: y + 103. + title_extra,
                w: pw - if heat { 172. } else { 100. },
                h: ph - if heat { 179. } else { 162. } - title_extra - caption_height,
                bounds: *b,
            };
            if !p.qc.boxes.is_empty() {
                a.w -= 100.;
            }
            let slot = a;
            let mut rows = legend_rows(p, a.w);
            loop {
                a = slot;
                let extra = rows.saturating_sub(1) as f64 * 16.;
                a.y += extra;
                a.h -= extra;
                ensure!(
                    a.w > 30. && a.h > 30.,
                    "plot area is too small for labels; increase height or reduce columns"
                );
                if matches!(
                    p.content,
                    Content::Heat {
                        square_cells: true,
                        ..
                    }
                ) || p.diagonal
                {
                    let ratio = (b.x[1] - b.x[0]) / (b.y[1] - b.y[0]);
                    let w = a.w.min(a.h * ratio);
                    let ht = w / ratio;
                    a.x += (a.w - w) / 2.;
                    a.y += (a.h - ht) / 2.;
                    a.w = w;
                    a.h = ht;
                }
                ensure!(a.w > 30. && a.h > 30., "plot area is too small");
                let needed = legend_rows(p, a.w);
                if needed <= rows {
                    break;
                }
                rows = needed;
            }
            if let Content::Heat {
                matrix,
                square_cells: false,
                ..
            } = &p.content
            {
                ensure!(
                    a.h >= matrix.len() as f64 * 11.,
                    "effect heatmap labels need more room; increase -height or select fewer -rsa-rois"
                );
            }
            text_regions.push(TextRegion {
                target: TextTarget::XLabel(i),
                rect: [a.x, a.y + a.h + if heat { 52. } else { 31. }, a.w, 18.],
            });
            text_regions.push(TextRegion {
                target: TextTarget::YLabel(i),
                rect: [a.x - if heat { 94. } else { 61. }, a.y, 20., a.h],
            });
            draw_plot(&mut svg, p, &a, i);
            draw_boxes(&mut svg, p, &a, i);
            areas.push(a);
        }
        text(
            &mut svg,
            margin,
            h - 14.,
            10.,
            MUTED,
            if fig.synthetic {
                "Includes synthetic demonstration data"
            } else {
                "1dPetePlot"
            },
            "start",
        );
        text(
            &mut svg,
            1172.,
            h - 14.,
            10.,
            MUTED,
            "1dPetePlot  /  AFNI",
            "end",
        );
        svg.push_str("</svg>");
        let tree =
            usvg::Tree::from_str(&svg, &self.options).context("preparing shared SVG scene")?;
        Ok(Scene {
            tree,
            areas,
            text_regions,
            logical_height: h,
        })
    }
    pub fn raster(&self, scene: &Scene) -> Result<tiny_skia::Pixmap> {
        let size = scene.tree.size().to_int_size();
        let mut pix =
            tiny_skia::Pixmap::new(size.width(), size.height()).context("allocating image")?;
        resvg::render(
            &scene.tree,
            tiny_skia::Transform::identity(),
            &mut pix.as_mut(),
        );
        Ok(pix)
    }
    pub fn export(
        &self,
        fig: &Figure,
        views: &[Bounds],
        width: u32,
        height: u32,
        path: &Path,
        overwrite: bool,
    ) -> Result<()> {
        let ext = path
            .extension()
            .and_then(|s| s.to_str())
            .unwrap_or("")
            .to_ascii_lowercase();
        ensure!(
            ext == "png" || ext == "svg",
            "unsupported output extension '.{ext}'; use .png or .svg"
        );
        ensure!(
            overwrite || !path.try_exists()?,
            "{} already exists; use -overwrite",
            path.display()
        );
        let scene = self.scene(fig, views, width, height)?;
        let data = match ext.as_str() {
            "png" => self.raster(&scene)?.encode_png()?,
            "svg" => scene
                .tree
                .to_string(&usvg::WriteOptions::default())
                .into_bytes(),
            _ => unreachable!(),
        };
        let parent = path
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        let mut temp = tempfile::NamedTempFile::new_in(parent)
            .with_context(|| format!("creating output in {}", parent.display()))?;
        temp.write_all(&data)?;
        temp.flush()?;
        let result = if overwrite {
            temp.persist(path)
        } else {
            temp.persist_noclobber(path)
        };
        if let Err(e) = result {
            bail!(
                "could not save {}: {} (existing output preserved)",
                path.display(),
                e.error
            );
        }
        Ok(())
    }
}
fn xml(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}
fn text_width(text: &str, size: f64) -> f64 {
    static FACE: std::sync::OnceLock<ttf_parser::Face<'static>> = std::sync::OnceLock::new();
    let face = FACE.get_or_init(|| ttf_parser::Face::parse(FONT, 0).expect("bundled font"));
    text.chars()
        .map(|c| {
            face.glyph_index(c)
                .and_then(|id| face.glyph_hor_advance(id))
                .unwrap_or(face.units_per_em()) as f64
                / face.units_per_em() as f64
                * size
        })
        .sum::<f64>()
}
fn wrap(s: &str, width: f64, size: f64) -> Vec<String> {
    let mut lines = vec![];
    for paragraph in s.split('\n') {
        let mut line = String::new();
        for word in paragraph.split_whitespace() {
            let proposed = if line.is_empty() {
                word.to_string()
            } else {
                format!("{line} {word}")
            };
            if !line.is_empty() && text_width(&proposed, size) > width - 2. {
                lines.push(std::mem::take(&mut line));
            }
            if !line.is_empty() {
                line.push(' ');
            }
            for c in word.chars() {
                if !line.is_empty() && text_width(&format!("{line}{c}"), size) > width - 2. {
                    lines.push(std::mem::take(&mut line));
                }
                line.push(c);
            }
        }
        if !line.is_empty() || paragraph.is_empty() {
            lines.push(line);
        }
    }
    if lines.is_empty() {
        lines.push(String::new());
    }
    lines
}
fn short(s: &str, n: usize) -> String {
    if s.chars().count() > n {
        format!(
            "{}…",
            s.chars().take(n.saturating_sub(1)).collect::<String>()
        )
    } else {
        s.to_owned()
    }
}
fn text(s: &mut String, x: f64, y: f64, size: f64, color: &str, t: &str, anchor: &str) {
    let _ = write!(
        s,
        r#"<text x="{x:.3}" y="{y:.3}" font-family="Ubuntu" font-weight="300" font-size="{size:.3}" fill="{color}" text-anchor="{anchor}">{}</text>"#,
        scientific_text(t)
    );
}
// Braces make scientific formatting explicit; file/column names keep literal underscores.
fn scientific_text(t: &str) -> String {
    let chars: Vec<_> = t.chars().collect();
    let mut out = String::new();
    let mut i = 0;
    while i < chars.len() {
        if ['^', '_'].contains(&chars[i])
            && chars.get(i + 1) == Some(&'{')
            && let Some(end) = chars[i + 2..].iter().position(|c| *c == '}')
        {
            let content: String = chars[i + 2..i + 2 + end].iter().collect();
            let shift = if chars[i] == '^' { "super" } else { "sub" };
            let _ = write!(
                out,
                r#"<tspan baseline-shift="{shift}" font-size="70%">{}</tspan>"#,
                xml(&content)
            );
            i += end + 3;
        } else {
            out.push_str(&xml(&chars[i].to_string()));
            i += 1;
        }
    }
    out
}
// SVG rectangle primitive: coordinates, paint, and corner radius.
#[allow(clippy::too_many_arguments)]
fn rect(
    s: &mut String,
    x: f64,
    y: f64,
    w: f64,
    h: f64,
    fill: &str,
    stroke: Option<&str>,
    radius: f64,
) {
    let _ = write!(
        s,
        r#"<rect x="{x:.3}" y="{y:.3}" width="{w:.3}" height="{h:.3}" rx="{radius}" fill="{fill}" stroke="{}"/>"#,
        stroke.unwrap_or("none")
    );
}
fn line(s: &mut String, p: [f64; 2], q: [f64; 2], color: &str, width: f64, dash: bool) {
    let _ = write!(
        s,
        r#"<path d="M {:.3} {:.3} L {:.3} {:.3}" fill="none" stroke="{color}" stroke-width="{width}" {} />"#,
        p[0],
        p[1],
        q[0],
        q[1],
        if dash {
            r#"stroke-dasharray="5 5""#
        } else {
            ""
        }
    );
}
fn ticks(range: [f64; 2]) -> Vec<f64> {
    let raw = (range[1] - range[0]) / 4.;
    let power = 10f64.powf(raw.log10().floor());
    let step = [1., 2., 2.5, 5., 10.]
        .into_iter()
        .find(|v| *v * power >= raw)
        .unwrap_or(10.)
        * power;
    let start = (range[0] / step).ceil() * step;
    (0..10)
        .map(|i| start + i as f64 * step)
        .take_while(|v| *v <= range[1] + step * 1e-9)
        .collect()
}
fn num(v: f64) -> String {
    if v.abs() < 1e-10 {
        "0".into()
    } else if v.abs() >= 10000. || v.abs() < 0.001 {
        format!("{v:.1e}")
    } else {
        format!("{v:.3}")
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_owned()
    }
}
fn palette(t: f64, cmap: &str) -> String {
    let t = t.clamp(0., 1.);
    let stops: &[[f64; 3]] = if cmap == "viridis" {
        &[
            [68., 1., 84.],
            [59., 82., 139.],
            [33., 145., 140.],
            [94., 201., 98.],
            [253., 231., 37.],
        ]
    } else if cmap == "Reds_and_Blues_Inv" {
        &[
            [33., 102., 172.],
            [146., 197., 222.],
            [247., 247., 247.],
            [244., 165., 130.],
            [178., 24., 43.],
        ]
    } else {
        &[
            [178., 24., 43.],
            [244., 165., 130.],
            [247., 247., 247.],
            [146., 197., 222.],
            [33., 102., 172.],
        ]
    };
    let z = t * (stops.len() - 1) as f64;
    let i = (z.floor() as usize).min(stops.len() - 2);
    let f = z - i as f64;
    let c = [0, 1, 2].map(|k| (stops[i][k] * (1. - f) + stops[i + 1][k] * f).round() as u8);
    format!("#{:02x}{:02x}{:02x}", c[0], c[1], c[2])
}
fn draw_plot(svg: &mut String, p: &Plot, a: &Area, id: usize) {
    rect(svg, a.x, a.y, a.w, a.h, "#fbfcfd", None, 0.);
    let heat = matches!(p.content, Content::Heat { .. });
    let _ = write!(
        svg,
        r#"<defs><clipPath id="plot{id}"><rect x="{}" y="{}" width="{}" height="{}"/></clipPath></defs>"#,
        a.x, a.y, a.w, a.h
    );
    let _ = write!(svg, r#"<g clip-path="url(#plot{id})">"#);
    // Independent annotations first; QC censor/run overlays retain visibility.
    for span in &p.qc.annotations {
        let lo = span.x[0].max(a.bounds.x[0]);
        let hi = span.x[1].min(a.bounds.x[1]);
        if hi <= lo {
            continue;
        }
        let l = a.point(lo, a.bounds.y[1])[0];
        let r = a.point(hi, a.bounds.y[1])[0];
        let _ = write!(
            svg,
            r#"<rect x="{l}" y="{}" width="{}" height="{}" fill="{}" fill-opacity="{}"/>"#,
            a.y,
            r - l,
            a.h,
            span.color,
            span.alpha
        );
    }
    for span in &p.qc.spans {
        let l = a.point(span.x[0], a.bounds.y[1])[0];
        let r = a.point(span.x[1], a.bounds.y[1])[0];
        let _ = write!(
            svg,
            r#"<rect x="{l}" y="{}" width="{}" height="{}" fill="{}" fill-opacity="{}"/>"#,
            a.y,
            r - l,
            a.h,
            span.color,
            span.alpha
        );
    }
    svg.push_str("</g>");
    if !heat {
        for t in ticks(a.bounds.y) {
            let q = a.point(a.bounds.x[0], t);
            line(svg, q, [a.x + a.w, q[1]], GRID, 1., false);
            if !p.qc.margin_off {
                text(svg, a.x - 11., q[1] + 4., 11., MUTED, &num(t), "end");
            }
        }
        for t in ticks(a.bounds.x) {
            let q = a.point(t, a.bounds.y[0]);
            line(svg, q, [q[0], a.y], GRID, 1., false);
            if !p.qc.margin_off {
                text(svg, q[0], a.y + a.h + 20., 11., MUTED, &num(t), "middle");
            }
        }
    }
    let _ = write!(svg, r#"<g clip-path="url(#plot{id})">"#);
    match &p.content {
        Content::Lines(series) | Content::Scatter(series) => {
            for s in series {
                if let Some(band) = &s.band {
                    let mut start = 0;
                    while start < s.points.len() {
                        let valid = |i: usize| {
                            s.points[i][1].is_finite() && band[i].iter().all(|v| v.is_finite())
                        };
                        if !valid(start) {
                            start += 1;
                            continue;
                        }
                        let mut end = start + 1;
                        while end < s.points.len() && valid(end) {
                            end += 1;
                        }
                        let mut path = String::new();
                        for (j, i) in (start..end).enumerate() {
                            let q = a.point(s.points[i][0], band[i][1]);
                            let _ = write!(
                                path,
                                "{} {:.3} {:.3} ",
                                if j == 0 { "M" } else { "L" },
                                q[0],
                                q[1]
                            );
                        }
                        for i in (start..end).rev() {
                            let q = a.point(s.points[i][0], band[i][0]);
                            let _ = write!(path, "L {:.3} {:.3} ", q[0], q[1]);
                        }
                        let _ = write!(
                            svg,
                            r#"<path d="{path} Z" fill="{}" fill-opacity="0.16"/>"#,
                            s.color
                        );
                        start = end;
                    }
                }
                if matches!(p.content, Content::Lines(_)) {
                    let mut path = String::new();
                    let mut pen = false;
                    for p in &s.points {
                        if !p[1].is_finite() {
                            pen = false;
                            continue;
                        }
                        let q = a.point(p[0], p[1]);
                        let _ = write!(
                            path,
                            "{} {:.3} {:.3} ",
                            if pen { "L" } else { "M" },
                            q[0],
                            q[1]
                        );
                        pen = true;
                    }
                    let _ = write!(
                        svg,
                        r#"<path d="{path}" fill="none" stroke="{}" stroke-width="2.1" stroke-linejoin="round" stroke-linecap="round"/>"#,
                        s.color
                    );
                    for (i, pt) in s.points.iter().enumerate() {
                        if !pt[1].is_finite()
                            || (i > 0 && s.points[i - 1][1].is_finite())
                            || (i + 1 < s.points.len() && s.points[i + 1][1].is_finite())
                        {
                            continue;
                        }
                        let q = a.point(pt[0], pt[1]);
                        let _ = write!(
                            svg,
                            r#"<circle cx="{}" cy="{}" r="3" fill="{}"/>"#,
                            q[0], q[1], s.color
                        );
                    }
                } else {
                    for pt in s.points.iter().filter(|pt| pt[1].is_finite()) {
                        let q = a.point(pt[0], pt[1]);
                        let _ = write!(
                            svg,
                            r#"<circle cx="{:.3}" cy="{:.3}" r="3.4" fill="{}" fill-opacity="0.72" stroke="white" stroke-width="0.7"/>"#,
                            q[0], q[1], s.color
                        );
                    }
                }
            }
        }
        Content::Histogram { series, edges, .. } => {
            for (s, e) in series.iter().zip(edges) {
                for (i, pt) in s.points.iter().enumerate() {
                    let l = a.point(e[i], pt[1]);
                    let r = a.point(e[i + 1], 0.);
                    let _ = write!(
                        svg,
                        r#"<rect x="{:.3}" y="{:.3}" width="{:.3}" height="{:.3}" fill="{}" fill-opacity="0.78" stroke="white" stroke-width="1"/>"#,
                        l[0],
                        l[1],
                        (r[0] - l[0]).max(0.),
                        (r[1] - l[1]).max(0.),
                        s.color
                    );
                }
            }
        }
        Content::Heat {
            matrix, cmap, clim, ..
        } => {
            let n = matrix.len();
            for (r, row) in matrix.iter().enumerate() {
                for (c, v) in row.iter().enumerate() {
                    let q = a.point(c as f64, (n - r) as f64);
                    let w = a.w / (a.bounds.x[1] - a.bounds.x[0]);
                    let h = a.h / (a.bounds.y[1] - a.bounds.y[0]);
                    rect(
                        svg,
                        q[0],
                        q[1],
                        w + 0.05,
                        h + 0.05,
                        &if v.is_finite() {
                            palette((v - clim[0]) / (clim[1] - clim[0]), cmap)
                        } else {
                            "#dce3e8".into()
                        },
                        None,
                        0.,
                    );
                }
            }
        }
    }
    for &h in &p.qc.hlines {
        line(
            svg,
            a.point(a.bounds.x[0], h),
            a.point(a.bounds.x[1], h),
            "#00a5b8",
            1.6,
            true,
        );
    }
    if p.diagonal {
        let lo = a.bounds.x[0].min(a.bounds.y[0]);
        let hi = a.bounds.x[1].max(a.bounds.y[1]);
        line(svg, a.point(lo, lo), a.point(hi, hi), "#8e9aa5", 1.3, true);
    }
    for (v, label) in &p.refs {
        line(
            svg,
            a.point(*v, a.bounds.y[0]),
            a.point(*v, a.bounds.y[1]),
            if label.is_some() {
                "#c45b53"
            } else {
                "#667985"
            },
            if label.is_some() { 2. } else { 1.2 },
            label.is_none(),
        );
    }
    svg.push_str("</g>");
    if p.qc.margin_off {
        return;
    }
    if let Content::Heat {
        matrix,
        labels,
        row_labels,
        square_cells,
        cmap,
        clim,
    } = &p.content
    {
        let n = matrix.len();
        let skip = labels.len().div_ceil(12);
        for i in (0..labels.len()).step_by(skip) {
            let q = a.point(i as f64 + 0.5, 0.);
            if q[0] >= a.x && q[0] <= a.x + a.w {
                let _ = write!(
                    svg,
                    r#"<g transform="translate({:.3} {:.3}) rotate(-40)">"#,
                    q[0],
                    a.y + a.h + 15.
                );
                text(svg, 0., 0., 9.5, MUTED, &short(&labels[i], 16), "end");
                svg.push_str("</g>");
            }
        }
        let skip = if *square_cells { n.div_ceil(12) } else { 1 };
        for i in (0..n).step_by(skip) {
            let q = a.point(0., (n - i) as f64 - 0.5);
            if q[1] >= a.y && q[1] <= a.y + a.h {
                text(
                    svg,
                    a.x - 9.,
                    q[1] + 3.,
                    9.5,
                    MUTED,
                    &short(&row_labels[i], 16),
                    "end",
                );
            }
        }
        let bx = a.x + a.w + 19.;
        for i in 0..100 {
            let t = i as f64 / 99.;
            rect(
                svg,
                bx,
                a.y + a.h * (1. - t),
                9.,
                a.h / 99. + 0.1,
                &palette(t, cmap),
                None,
                0.,
            );
        }
        for v in [clim[0], (clim[0] + clim[1]) / 2., clim[1]] {
            let y = a.y + a.h * (1. - (v - clim[0]) / (clim[1] - clim[0]));
            text(svg, bx + 15., y + 3., 10., MUTED, &num(v), "start");
        }
    } else {
        line(
            svg,
            [a.x, a.y + a.h],
            [a.x + a.w, a.y + a.h],
            "#bdcbd3",
            1.,
            false,
        );
        line(svg, [a.x, a.y], [a.x, a.y + a.h], "#bdcbd3", 1., false);
    }
    text(
        svg,
        a.x + a.w / 2.,
        a.y + a.h + if heat { 65. } else { 44. },
        12f64.min(12. * a.w / text_width(&p.xlabel, 12.).max(1.)),
        INK,
        &p.xlabel,
        "middle",
    );
    let _ = write!(
        svg,
        r#"<g transform="translate({:.3} {:.3}) rotate(-90)">"#,
        a.x - if heat { 85. } else { 52. },
        a.y + a.h / 2.
    );
    let label = if let Some(n) = p.qc.label_wrap {
        p.ylabel
            .chars()
            .collect::<Vec<_>>()
            .chunks(n)
            .map(|c| c.iter().collect::<String>())
            .collect::<Vec<_>>()
    } else {
        vec![p.ylabel.clone()]
    };
    for (i, l) in label.iter().enumerate() {
        text(
            svg,
            0.,
            (i as f64 - (label.len() - 1) as f64 / 2.) * 13.,
            12f64.min(12. * a.h / text_width(l, 12.).max(1.)),
            INK,
            l,
            "middle",
        );
    }
    svg.push_str("</g>");
    let mut lx = a.x;
    let mut ly = a.y - 15. - legend_rows(p, a.w).saturating_sub(1) as f64 * 16.;
    for (c, l) in legend_entries(p) {
        let label = short(l, ((a.w - 33.) / 5.5).clamp(1., 24.) as usize);
        let w = label.chars().count() as f64 * 5.5 + 33.;
        if lx > a.x && lx + w > a.x + a.w {
            lx = a.x;
            ly += 16.;
        }
        line(svg, [lx, ly], [lx + 15., ly], c, 2.5, false);
        text(svg, lx + 21., ly + 4., 10., MUTED, &label, "start");
        lx += w;
    }
}

fn legend_entries(p: &Plot) -> Vec<(&str, &str)> {
    let mut legend = vec![];
    for span in &p.qc.annotations {
        if let Some(label) = &span.label {
            let entry = (span.color.as_str(), label.as_str());
            if !legend.contains(&entry) {
                legend.push(entry);
            }
        }
    }
    if p.legend {
        match &p.content {
            Content::Lines(s) | Content::Scatter(s) | Content::Histogram { series: s, .. } => {
                for v in s {
                    legend.push((v.color.as_str(), v.label.as_str()));
                }
            }
            _ => {}
        }
    }
    for (_, l) in &p.refs {
        if let Some(l) = l {
            legend.push(("#c45b53", l));
        }
    }
    legend
}
fn legend_rows(p: &Plot, width: f64) -> usize {
    let mut rows = 1;
    let mut used = 0.;
    for (_, label) in legend_entries(p) {
        let label = short(label, ((width - 33.) / 5.5).clamp(1., 24.) as usize);
        let w = label.chars().count() as f64 * 5.5 + 33.;
        if used > 0. && used + w > width {
            rows += 1;
            used = 0.;
        }
        used += w;
    }
    rows
}

fn draw_boxes(svg: &mut String, p: &Plot, a: &Area, id: usize) {
    if p.qc.boxes.is_empty() {
        return;
    }
    let left = a.x + a.w + 22.;
    text(svg, left + 31., a.y - 15., 10., MUTED, "SUMMARY", "middle");
    let _ = write!(
        svg,
        r#"<defs><clipPath id="boxes{id}"><rect x="{left}" y="{}" width="74" height="{}"/></clipPath></defs>"#,
        a.y, a.h
    );
    for (i, b) in p.qc.boxes.iter().enumerate() {
        let x = left + 17. + i as f64 * 38.;
        let color = if b.label == "AC" {
            "#147d92"
        } else {
            "#89929e"
        };
        let _ = write!(svg, r#"<g clip-path="url(#boxes{id})">"#);
        if b.n > 0 {
            let y = b.five.map(|v| a.point(a.bounds.x[0], v)[1]);
            line(svg, [x, y[0]], [x, y[4]], color, 1.3, false);
            for j in [0, 4] {
                line(svg, [x - 7., y[j]], [x + 7., y[j]], color, 1.3, false);
            }
            rect(
                svg,
                x - 10.,
                y[3],
                20.,
                (y[1] - y[3]).max(1.),
                color,
                None,
                2.,
            );
            line(svg, [x - 10., y[2]], [x + 10., y[2]], "white", 2., false);
            for &v in &b.outliers {
                let y = a.point(a.bounds.x[0], v)[1];
                let _ = write!(
                    svg,
                    r#"<circle cx="{x}" cy="{y}" r="2.3" stroke="{color}" fill="none"/>"#
                );
            }
        }
        svg.push_str("</g>");
        text(svg, x, a.y + a.h + 20., 10., INK, &b.label, "middle");
        text(
            svg,
            x,
            a.y + a.h + 34.,
            9.,
            MUTED,
            &format!("n={}", b.n),
            "middle",
        );
        if b.n == 0 {
            text(svg, x, a.y + a.h / 2., 9., MUTED, "empty", "middle");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model;
    #[test]
    fn captions_wrap_using_font_metrics_without_omitting_text_or_covering_axes() {
        let source = "WWWWWWWWWWWWWWWWWWWWWWWWW narrow words and αβγ\nA second paragraph.";
        let lines = wrap(source, 180., 13.);
        assert!(lines.iter().all(|l| text_width(l, 13.) <= 180.));
        assert_eq!(
            lines
                .concat()
                .chars()
                .filter(|c| !c.is_whitespace())
                .collect::<String>(),
            source
                .chars()
                .filter(|c| !c.is_whitespace())
                .collect::<String>()
        );
        let mut fig = crate::model::demo();
        fig.caption = source.into();
        fig.plots[0].caption = source.into();
        let renderer = Renderer::new();
        let views: Vec<_> = fig.plots.iter().map(|p| p.bounds).collect();
        let scene = renderer.scene(&fig, &views, 1500, 1400).unwrap();
        let caption = scene
            .text_regions
            .iter()
            .find(|r| r.target == TextTarget::PanelCaption(0))
            .unwrap();
        let axis = scene
            .text_regions
            .iter()
            .find(|r| r.target == TextTarget::XLabel(0))
            .unwrap();
        assert!(caption.rect[1] >= axis.rect[1] + axis.rect[3]);
        assert!(
            scene
                .text_regions
                .iter()
                .any(|r| r.target == TextTarget::FigureCaption)
        );
    }
    #[test]
    fn output_dimensions_and_atomic_replacement() {
        let d = tempfile::tempdir().unwrap();
        let f = model::demo();
        let v = f.plots.iter().map(|p| p.bounds).collect::<Vec<_>>();
        let r = Renderer::new();
        let p = d.path().join("test.png");
        r.export(&f, &v, 1200, 1000, &p, false).unwrap();
        let bytes = std::fs::read(&p).unwrap();
        assert_eq!(u32::from_be_bytes(bytes[16..20].try_into().unwrap()), 1200);
        assert_eq!(u32::from_be_bytes(bytes[20..24].try_into().unwrap()), 1000);
        assert!(r.export(&f, &v, 1200, 1000, &p, false).is_err());
        assert_eq!(std::fs::read(&p).unwrap(), bytes);
        assert!(r.export(&f, &v, 1200, 100, &p, true).is_err());
        assert_eq!(std::fs::read(&p).unwrap(), bytes);
        r.export(&f, &v, 1200, 1100, &p, true).unwrap();
        assert_ne!(std::fs::read(&p).unwrap(), bytes);
        let p = d.path().join("test.svg");
        r.export(&f, &v, 1200, 1000, &p, false).unwrap();
        let s = std::fs::read_to_string(p).unwrap();
        assert!(s.contains("<path"));
        assert!(!s.contains("<text"));
        assert!(!s.contains("<image"));
    }
    #[test]
    fn mapping_roundtrip() {
        let a = Area {
            x: 70.,
            y: 50.,
            w: 500.,
            h: 200.,
            bounds: Bounds {
                x: [-2., 3.],
                y: [0., 1.],
            },
        };
        let p = a.point(1., 0.7);
        let d = a.data(p[0], p[1]);
        assert!((d[0] - 1.).abs() < 1e-12 && (d[1] - 0.7).abs() < 1e-12);
    }
    #[test]
    fn saved_vector_scene_matches_cpu_preview() {
        let r = Renderer::new();
        let f = model::demo();
        let views = f.plots.iter().map(|p| p.bounds).collect::<Vec<_>>();
        let scene = r.scene(&f, &views, 1200, 1000).unwrap();
        let preview = r.raster(&scene).unwrap();
        let svg = scene.tree.to_string(&usvg::WriteOptions::default());
        // Re-open with an empty font database: SVG must carry all label geometry.
        let tree = usvg::Tree::from_str(&svg, &usvg::Options::default()).unwrap();
        let saved = r
            .raster(&Scene {
                tree,
                areas: vec![],
                text_regions: vec![],
                logical_height: 1000.,
            })
            .unwrap();
        let mismatches = preview
            .data()
            .iter()
            .zip(saved.data())
            .filter(|(a, b)| a.abs_diff(**b) > 3)
            .count();
        assert!(
            mismatches < preview.data().len() / 1000,
            "vector export changed too many raster channels: {mismatches}"
        );
    }
    #[test]
    fn missing_samples_break_both_paths_and_bands_and_annotations_clip() {
        let args = [
            "-infiles",
            "1D: 0 1 NaN 3 4",
            "-missing",
            "gap",
            "-yfiles_pm",
            "1D: 5@0.2",
            "-censor_trs",
            "1",
            "-boxplot_on",
        ];
        let o = crate::cli::parse(&args.map(String::from)).unwrap();
        let f = model::from_options(&o).unwrap();
        let a = Area {
            x: 70.,
            y: 50.,
            w: 500.,
            h: 200.,
            bounds: f.plots[0].bounds,
        };
        let mut svg = String::new();
        draw_plot(&mut svg, &f.plots[0], &a, 0);
        assert!(!svg.contains("NaN"));
        assert!(!svg.contains("inf"));
        assert_eq!(
            svg.matches("fill-opacity=\"0.16\"").count(),
            2,
            "uncertainty must split at the missing sample"
        );
        let path = svg
            .split("stroke-width=\"2.1\"")
            .next()
            .unwrap()
            .rsplit("<path")
            .next()
            .unwrap();
        assert_eq!(
            path.matches("M ").count(),
            2,
            "line must not bridge the missing sample"
        );
        assert!(svg.contains("clip-path=\"url(#plot0)\""));
    }
    #[test]
    fn scientific_labels_preserve_literal_underscores_and_escape_xml() {
        assert_eq!(scientific_text("model_a & b"), "model_a &amp; b");
        let t = scientific_text("β_{model} r^{2}");
        assert!(t.contains("baseline-shift=\"sub\""));
        assert!(t.contains("baseline-shift=\"super\""));
        assert!(scientific_text("r^{<x>}").contains("&lt;x&gt;"));
    }
    #[test]
    fn shaded_intervals_clip_and_wrapped_legends_keep_equal_scatter_axes() {
        let mut args = [
            "-scatter",
            "-diagonal",
            "-xfile",
            "1D: 0 1 2",
            "-infiles",
            "1D: 0 1 2",
        ]
        .map(String::from)
        .to_vec();
        for i in 0..8 {
            args.extend([
                "-shade".into(),
                "-10:10".into(),
                "-shade-label".into(),
                format!("Long condition {i}"),
            ]);
        }
        let o = crate::cli::parse(&args).unwrap();
        let f = model::from_options(&o).unwrap();
        let r = Renderer::new();
        let scene = r.scene(&f, &[f.plots[0].bounds], 1200, 1000).unwrap();
        let a = scene.areas[0];
        assert!((a.w / a.h - 1.).abs() < 1e-10);
        assert!(legend_rows(&f.plots[0], a.w) > 1);
        assert_eq!(legend_entries(&f.plots[0]).len(), 8);
        let mut svg = String::new();
        draw_plot(&mut svg, &f.plots[0], &a, 0);
        assert_eq!(svg.matches("fill-opacity=\"0.18\"").count(), 8);
        assert!(!svg.contains("NaN"));
    }
}
