use crate::{cli::Options, input};
use anyhow::{Result, ensure};

pub const COLORS: [&str; 6] = [
    "#147d92", "#e18a47", "#7064b5", "#ca5773", "#5b9162", "#557cac",
];
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Bounds {
    pub x: [f64; 2],
    pub y: [f64; 2],
}
#[cfg_attr(not(feature = "gui"), allow(dead_code))]
impl Bounds {
    pub fn valid(self) -> bool {
        self.x.iter().chain(self.y.iter()).all(|v| v.is_finite())
            && self.x[1] > self.x[0]
            && self.y[1] > self.y[0]
            && (self.x[1] - self.x[0]).is_finite()
            && (self.y[1] - self.y[0]).is_finite()
    }
    #[cfg(test)]
    pub fn zoom(self, x: f64, y: f64, factor: f64) -> Self {
        Self {
            x: self.x.map(|v| x + (v - x) * factor),
            y: self.y.map(|v| y + (v - y) * factor),
        }
    }
    pub fn pan(self, dx: f64, dy: f64) -> Self {
        Self {
            x: self.x.map(|v| v + dx),
            y: self.y.map(|v| v + dy),
        }
    }
}
#[derive(Clone, Debug)]
pub struct Series {
    pub label: String,
    pub color: String,
    pub points: Vec<[f64; 2]>,
    pub band: Option<Vec<[f64; 2]>>,
    pub point_labels: Vec<String>,
}
#[derive(Clone, Debug)]
pub struct IntervalRow {
    pub label: String,
    pub detail: String,
    /// Estimate, lower endpoint, upper endpoint. None means an incomplete triplet.
    pub values: Option<[f64; 3]>,
}
#[derive(Clone, Debug)]
pub enum Content {
    Lines(Vec<Series>),
    Scatter(Vec<Series>),
    Intervals {
        rows: Vec<IntervalRow>,
        color: String,
        method: String,
    },
    #[cfg_attr(not(feature = "gui"), allow(dead_code))]
    Histogram {
        series: Vec<Series>,
        edges: Vec<Vec<f64>>,
        density: bool,
    },
    Heat {
        matrix: Vec<Vec<f64>>,
        labels: Vec<String>,
        row_labels: Vec<String>,
        square_cells: bool,
        named_conditions: bool,
        cmap: String,
        clim: [f64; 2],
    },
}
#[derive(Clone, Debug)]
pub struct Plot {
    pub style: crate::style::PlotStyle,
    pub id: usize,
    pub caption: String,
    pub title: String,
    pub note: String,
    pub xlabel: String,
    pub ylabel: String,
    pub content: Content,
    pub bounds: Bounds,
    pub diagonal: bool,
    pub refs: Vec<(f64, Option<String>)>,
    pub legend: bool,
    pub qc: crate::qc::Decoration,
}
#[derive(Clone, Debug)]
pub struct Figure {
    pub caption: String,
    pub title: String,
    pub subtitle: String,
    pub plots: Vec<Plot>,
    pub columns: usize,
    pub synthetic: bool,
}
impl Figure {
    pub fn empty() -> Self {
        Self {
            title: "My figure".into(),
            caption: String::new(),
            subtitle: "Imported data".into(),
            plots: vec![],
            columns: 1,
            synthetic: false,
        }
    }
}
impl Plot {
    pub fn kind(&self) -> &str {
        match self.content {
            Content::Lines(_) => "TIME SERIES",
            Content::Scatter(_) if self.qc.xticks.is_some() => "RUN SUMMARIES",
            Content::Scatter(_) => "SCATTER",
            Content::Histogram { .. } => "DISTRIBUTION",
            Content::Heat { .. } => "HEATMAP",
            Content::Intervals { .. } => "CONFIDENCE INTERVALS",
        }
    }
    #[cfg_attr(not(feature = "gui"), allow(dead_code))]
    pub fn hover(&self, x: f64, y: f64) -> String {
        match &self.content {
            Content::Intervals { rows, method, .. } => {
                if y < 0. || y >= rows.len() as f64 {
                    return String::new();
                }
                let row = &rows[rows.len() - 1 - y.floor() as usize];
                match row.values {
                    Some([estimate, lo, hi]) => format!(
                        "{} · {} · estimate {estimate:.5} · CI [{lo:.5}, {hi:.5}] · {method}",
                        row.detail, self.xlabel
                    ),
                    None => format!("{} · incomplete effect/CI · {method}", row.detail),
                }
            }
            Content::Heat {
                matrix,
                labels,
                row_labels,
                ..
            } => {
                let n = matrix.len();
                if x < 0. || x >= labels.len() as f64 || y < 0. || y >= n as f64 {
                    return String::new();
                }
                let col = x.floor() as usize;
                let row = n - 1 - y.floor() as usize;
                format!(
                    "{} × {}   ·   {:.5}",
                    row_labels[row], labels[col], matrix[row][col]
                )
            }
            Content::Histogram {
                series,
                edges,
                density,
            } => {
                let mut out = vec![];
                for (si, (s, e)) in series.iter().zip(edges).enumerate() {
                    if let Some(i) = e.windows(2).position(|v| x >= v[0] && x <= v[1]) {
                        out.push(format!(
                            "{}: [{:.3}, {:.3}] {} = {:.4}",
                            self.style.label(si, &s.label),
                            e[i],
                            e[i + 1],
                            if *density { "density" } else { "count" },
                            s.points[i][1]
                        ));
                    }
                }
                out.join("  |  ")
            }
            Content::Lines(ss) | Content::Scatter(ss) => {
                let dx = self.bounds.x[1] - self.bounds.x[0];
                let dy = self.bounds.y[1] - self.bounds.y[0];
                let best = ss
                    .iter()
                    .enumerate()
                    .flat_map(|(si, s)| {
                        s.points
                            .iter()
                            .enumerate()
                            .filter(|(_, p)| p[1].is_finite())
                            .map(move |(i, p)| (s, i, p, si))
                    })
                    .min_by(|a, b| {
                        let d =
                            |p: &[f64; 2]| ((p[0] - x) / dx).powi(2) + ((p[1] - y) / dy).powi(2);
                        d(a.2).total_cmp(&d(b.2))
                    });
                best.map(|(s, i, p, si)| {
                    format!(
                        "{} · {}   x = {:.4}   y = {:.4}",
                        self.style.label(si, &s.label),
                        s.point_labels
                            .get(i)
                            .cloned()
                            .unwrap_or_else(|| format!("sample {i}")),
                        p[0],
                        p[1]
                    ) + if self.qc.censored.get(i) == Some(&true) {
                        " · CENSORED"
                    } else {
                        ""
                    }
                })
                .unwrap_or_default()
            }
        }
    }
}
pub fn padded(lo: f64, hi: f64) -> [f64; 2] {
    let gap = if lo == hi {
        lo.abs().max(1.) * 0.1
    } else {
        (hi - lo) * 0.07
    };
    [lo - gap, hi + gap]
}
pub(crate) fn bounds(ss: &[Series], include_band: bool) -> Bounds {
    let mut x = [f64::INFINITY, f64::NEG_INFINITY];
    let mut y = x;
    for s in ss {
        for p in &s.points {
            x[0] = x[0].min(p[0]);
            x[1] = x[1].max(p[0]);
            y[0] = y[0].min(p[1]);
            y[1] = y[1].max(p[1]);
        }
        if include_band && let Some(b) = &s.band {
            for v in b {
                y[0] = y[0].min(v[0]);
                y[1] = y[1].max(v[1]);
            }
        }
    }
    Bounds {
        x: padded(x[0], x[1]),
        y: padded(y[0], y[1]),
    }
}
pub fn histogram(
    values: &[f64],
    nbins: Option<usize>,
    width: Option<f64>,
    density: bool,
) -> Result<(Vec<f64>, Vec<[f64; 2]>)> {
    ensure!(
        !values.is_empty() && values.iter().all(|v| v.is_finite()),
        "histogram requires finite observations"
    );
    let mut lo = values.iter().copied().fold(f64::INFINITY, f64::min);
    let mut hi = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if lo == hi {
        let half = width.unwrap_or(1.) / 2.;
        lo -= half;
        hi += half;
    }
    let n = if let Some(w) = width {
        ensure!(w.is_finite() && w > 0., "bin width must be positive");
        ((hi - lo) / w).ceil() as usize
    } else {
        nbins
            .unwrap_or((values.len() as f64).sqrt() as usize)
            .max(if nbins.is_some() { 1 } else { 8 })
    };
    ensure!(
        (1..=10000).contains(&n),
        "histogram requires 1..10000 bins; adjust -histwidth"
    );
    let step = width.unwrap_or((hi - lo) / n as f64);
    ensure!(
        step.is_finite() && step > 0. && (lo + step) > lo,
        "bin width is too small or data range too large"
    );
    let edges: Vec<f64> = (0..=n).map(|i| lo + i as f64 * step).collect();
    ensure!(
        edges.iter().all(|v| v.is_finite()),
        "histogram edges overflow"
    );
    let mut count = vec![0.; n];
    for &v in values {
        let i = (((v - lo) / step).floor() as usize).min(n - 1);
        count[i] += 1.;
    }
    let pts = count
        .iter()
        .enumerate()
        .map(|(i, c)| {
            [
                (edges[i] + edges[i + 1]) / 2.,
                if density {
                    c / (values.len() as f64 * step)
                } else {
                    *c
                },
            ]
        })
        .collect();
    Ok((edges, pts))
}
pub fn color(value: &str) -> Result<String> {
    if let Ok(v) = value.parse::<f64>() {
        ensure!(
            v.is_finite() && (0.0..=1.0).contains(&v),
            "grayscale must be 0..1"
        );
        let v = (v * 255.).round() as u8;
        return Ok(format!("#{v:02x}{v:02x}{v:02x}"));
    }
    if value.len() == 4
        && value.starts_with('#')
        && value[1..].chars().all(|c| c.is_ascii_hexdigit())
    {
        return Ok(format!(
            "#{}",
            value[1..].chars().flat_map(|c| [c, c]).collect::<String>()
        ));
    }
    let c = match value {
        "black" => "#243344",
        "red" => "#c44755",
        "blue" => "#3c70c9",
        "green" => "#4a9563",
        "orange" => "#e18a47",
        "purple" => "#7064b5",
        "gray" | "grey" => "#89929e",
        "white" => "#ffffff",
        "cyan" => "#00b7c7",
        "yellow" => "#ffd766",
        "pink" => "#ffb3b3",
        "violet" => "#a67bd4",
        x => x,
    };
    ensure!(
        c.len() == 7 && c.starts_with('#') && c[1..].chars().all(|x| x.is_ascii_hexdigit()),
        "unsupported color '{value}': use #RRGGBB or a common color name"
    );
    Ok(c.into())
}
pub fn from_options(o: &Options) -> Result<Figure> {
    if !o.rdm.files.is_empty() {
        return crate::rdm::figure(o);
    }
    if o.rsa.file.is_some() {
        return crate::rsa::figure(o);
    }
    if let Some(name) = &o.demo {
        let mut f = demo();
        if name != "gallery" {
            let i = match name.as_str() {
                "line" => 0,
                "histogram" => 1,
                "scatter" => 2,
                _ => 3,
            };
            f.plots = vec![f.plots.remove(i)];
            f.columns = 1;
            f.title = f.plots[0].title.clone();
        }
        if let Some(t) = &o.title {
            f.title = t.clone();
        }
        return presentation(f, o);
    }
    let stdin_count = o
        .files
        .iter()
        .chain(o.xfile.iter())
        .chain(o.reorder.iter())
        .chain(o.qc.pm.iter())
        .chain(o.qc.lower.iter())
        .chain(o.qc.upper.iter())
        .chain(o.qc.censor_files.iter())
        .filter(|s| s.trim_end_matches('\'').split(['[', '{']).next() == Some("-"))
        .count();
    ensure!(stdin_count <= 1, "stdin can only be consumed once");
    let mut plots = vec![];
    if o.kind == "heat" {
        ensure!(o.files.len() == 1, "heatmap needs exactly one input matrix");
        let mut m = input::read_with(&o.files[0], o.missing)?;
        let n = m.len();
        ensure!(
            m.iter().all(|r| r.len() == n),
            "heatmap matrix must be square"
        );
        ensure!(
            n <= 1000,
            "prototype heatmaps support at most 1000 conditions"
        );
        ensure!(
            ["viridis", "Reds_and_Blues", "Reds_and_Blues_Inv"].contains(&o.cmap.as_str()),
            "unsupported color map '{}'",
            o.cmap
        );
        let mut labels = (0..n).map(|i| i.to_string()).collect::<Vec<_>>();
        if let Some(file) = &o.reorder {
            let order = input::column(file)?;
            ensure!(
                order.len() == n
                    && order
                        .iter()
                        .all(|v| v.fract() == 0. && *v >= 0. && *v < n as f64),
                "reorder must contain N integer indices in 0..N-1"
            );
            let order: Vec<usize> = order.iter().map(|v| *v as usize).collect();
            let mut sorted = order.clone();
            sorted.sort_unstable();
            ensure!(
                sorted == (0..n).collect::<Vec<_>>(),
                "reorder must be a permutation without duplicates"
            );
            m = order
                .iter()
                .map(|&r| order.iter().map(|&c| m[r][c]).collect())
                .collect();
            labels = order.iter().map(|i| i.to_string()).collect();
        }
        let lo = m.iter().flatten().copied().fold(f64::INFINITY, f64::min);
        let hi = m
            .iter()
            .flatten()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        ensure!(
            lo.is_finite() && hi.is_finite(),
            "heatmap has no finite values"
        );
        let clim = o.clim.unwrap_or_else(|| {
            if o.zerocenter {
                let a = lo.abs().max(hi.abs()).max(f64::EPSILON);
                [-a, a]
            } else if lo == hi {
                padded(lo, hi)
            } else {
                [lo, hi]
            }
        });
        ensure!(o.ylabels.len() <= 1, "heatmap accepts one y-axis label");
        plots.push(Plot {
            style: Default::default(),
            id: plots.len(),
            caption: String::new(),
            title: "Representational dissimilarity".into(),
            note: format!("{n} × {n} matrix"),
            xlabel: if o.xlabel.is_empty() {
                "Condition".into()
            } else {
                o.xlabel.clone()
            },
            ylabel: o.ylabels.first().cloned().unwrap_or("Condition".into()),
            content: Content::Heat {
                matrix: m,
                row_labels: labels.clone(),
                labels,
                square_cells: true,
                named_conditions: false,
                cmap: o.cmap.clone(),
                clim,
            },
            bounds: Bounds {
                x: [0., n as f64],
                y: [0., n as f64],
            },
            diagonal: false,
            refs: vec![],
            legend: false,
            qc: Default::default(),
        });
    } else {
        let x = if let Some(file) = &o.xfile {
            Some(input::column(file)?)
        } else {
            None
        };
        let mut ss = vec![];
        for file in &o.files {
            let m = input::read_with(file, o.missing)?;
            for c in 0..m[0].len() {
                let y = m.iter().map(|r| r[c]).collect::<Vec<_>>();
                if let Some(x) = &x {
                    ensure!(x.len() == y.len(), "x/y length mismatch for {file}");
                }
                ensure!(
                    y.iter().any(|v| v.is_finite()),
                    "{file}: selected column has no finite values"
                );
                let i = ss.len();
                let label = std::path::Path::new(file)
                    .file_name()
                    .unwrap_or_default()
                    .to_string_lossy();
                ss.push(Series {
                    label: if m[0].len() > 1 {
                        format!("{label} [{c}]")
                    } else {
                        label.into()
                    },
                    color: if o.colors.is_empty() {
                        COLORS[i % COLORS.len()].into()
                    } else {
                        color(&o.colors[i % o.colors.len()])?
                    },
                    points: y
                        .iter()
                        .enumerate()
                        .map(|(j, v)| {
                            [
                                x.as_ref()
                                    .map(|x| x[j])
                                    .unwrap_or(o.xzero + j as f64 * o.dx),
                                *v,
                            ]
                        })
                        .collect(),
                    band: None,
                    point_labels: vec![],
                });
            }
        }
        ensure!(
            o.ylabels.is_empty() || o.ylabels.len() == ss.len(),
            "-ylabels must match the number of selected columns"
        );
        ensure!(
            o.legends.is_empty() || o.legends.len() == ss.len(),
            "-legend_labels must match the number of selected columns"
        );
        for (i, s) in ss.iter_mut().enumerate() {
            if let Some(l) = o.legends.get(i) {
                s.label = l.clone();
            }
        }
        let decorations = if o.kind == "line" {
            crate::qc::prepare(&mut ss, &o.qc, o.missing)?
        } else {
            vec![crate::qc::Decoration::default(); ss.len()]
        };
        if let Some(column) = &o.group_by {
            ss = crate::grouping::split(ss, &o.files[0], column)?;
            if !o.colors.is_empty() {
                for (i, s) in ss.iter_mut().enumerate() {
                    s.color = color(&o.colors[i % o.colors.len()])?;
                }
            }
        }
        let common = bounds(&ss, o.qc.use_pm);
        ensure!(
            o.qc.ylimits.is_empty()
                || o.qc.ylimits.len() == 1
                || (!o.one && o.qc.ylimits.len() == ss.len()),
            "yaxis needs one range or one per separate panel"
        );
        let groups: Vec<Vec<Series>> = if o.one {
            vec![ss]
        } else {
            ss.into_iter().map(|s| vec![s]).collect()
        };
        for (i, mut series) in groups.into_iter().enumerate() {
            let mut decoration = decorations[i].clone();
            if o.one {
                decoration.hlines = decorations
                    .iter()
                    .flat_map(|d| d.hlines.iter().copied())
                    .collect();
            }

            let title = if o.one {
                "All series".into()
            } else {
                series[0].label.clone()
            };
            let mut b = common;
            if o.sepscl {
                b.y = bounds(&series, o.qc.use_pm).y;
            }
            let (content, ylabel, xlabel) = if o.kind == "histogram" {
                let mut edges = vec![];
                let shared_edges = if o.group_by.is_some() {
                    let all: Vec<_> = series
                        .iter()
                        .flat_map(|s| s.points.iter().map(|p| p[1]))
                        .filter(|v| v.is_finite())
                        .collect();
                    Some(histogram(&all, o.bins, o.binwidth, false)?.0)
                } else {
                    None
                };
                for s in &mut series {
                    let values = s
                        .points
                        .iter()
                        .map(|p| p[1])
                        .filter(|v| v.is_finite())
                        .collect::<Vec<_>>();
                    let (e, p) = if let Some(edges) = &shared_edges {
                        let mut counts = vec![0.; edges.len() - 1];
                        for v in &values {
                            let index = edges
                                .partition_point(|e| e <= v)
                                .saturating_sub(1)
                                .min(counts.len() - 1);
                            counts[index] += 1.;
                        }
                        let points = counts
                            .iter()
                            .enumerate()
                            .map(|(i, count)| {
                                [
                                    (edges[i] + edges[i + 1]) / 2.,
                                    if o.density {
                                        count / (values.len() as f64 * (edges[i + 1] - edges[i]))
                                    } else {
                                        *count
                                    },
                                ]
                            })
                            .collect();
                        (edges.clone(), points)
                    } else {
                        histogram(&values, o.bins, o.binwidth, o.density)?
                    };
                    edges.push(e);
                    s.points = p;
                }
                b = bounds(&series, o.qc.use_pm);
                b.x = [
                    edges.iter().map(|e| e[0]).fold(f64::INFINITY, f64::min),
                    edges
                        .iter()
                        .map(|e| *e.last().unwrap())
                        .fold(f64::NEG_INFINITY, f64::max),
                ];
                b.x = padded(b.x[0], b.x[1]);
                b.y[0] = 0.;
                (
                    Content::Histogram {
                        series,
                        edges,
                        density: o.density,
                    },
                    if o.density {
                        "Density".into()
                    } else {
                        "Count".into()
                    },
                    if o.xlabel.is_empty() {
                        "Value".into()
                    } else {
                        o.xlabel.clone()
                    },
                )
            } else {
                let ylabel = o.ylabels.get(i).cloned().unwrap_or("Value".into());
                let xlabel = if o.xlabel.is_empty() {
                    "Sample".into()
                } else {
                    o.xlabel.clone()
                };
                (
                    if o.kind == "scatter" {
                        Content::Scatter(series)
                    } else {
                        Content::Lines(series)
                    },
                    ylabel,
                    xlabel,
                )
            };
            for (x, _) in &o.vlines {
                b.x[0] = b.x[0].min(*x);
                b.x[1] = b.x[1].max(*x);
            }
            if let Some(x) = o.xlim {
                b.x = x;
            }
            if let Some(y) =
                o.qc.ylimits
                    .get(if o.qc.ylimits.len() == 1 { 0 } else { i })
            {
                b.y = [y[0].unwrap_or(b.y[0]), y[1].unwrap_or(b.y[1])];
            }
            if o.diagonal {
                ensure!(
                    o.xlim.is_none() && o.ylim.is_none() || o.xlim == o.ylim,
                    "-diagonal needs identical explicit x/y ranges (or neither)"
                );
                let r = [b.x[0].min(b.y[0]), b.x[1].max(b.y[1])];
                b = Bounds { x: r, y: r };
            }
            plots.push(Plot {
                style: Default::default(),
                id: plots.len(),
                caption: String::new(),
                title,
                note: if decoration.censored.iter().any(|v| *v) {
                    format!(
                        "{} / {} samples censored · shaded intervals",
                        decoration.censored.iter().filter(|v| **v).count(),
                        decoration.censored.len()
                    )
                } else if o.diagonal {
                    "Identity line · equal axes".into()
                } else {
                    String::new()
                },
                xlabel,
                ylabel,
                content,
                bounds: b,
                diagonal: o.diagonal,
                refs: o.vlines.clone(),
                legend: o.legend_on || o.one,
                qc: decoration,
            });
        }
    }
    ensure!(
        plots.len() <= 64,
        "supports at most 64 panels; select fewer columns"
    );
    ensure!(
        !o.qc.margin_off || plots.len() == 1,
        "margin_off needs a single panel; use -one_graph"
    );
    let mut annotations = o.shades.clone();
    if o.shade_design.file.is_some() {
        let Content::Lines(series) = &plots[0].content else {
            unreachable!()
        };
        let x = series[0].points.iter().map(|p| p[0]).collect::<Vec<_>>();
        annotations.extend(crate::shading::design_spans(&o.shade_design, &x)?);
    }
    for p in &mut plots {
        p.qc.annotations = annotations.clone();
        let mut notes = vec![];
        if !p.note.is_empty() {
            notes.push(p.note.clone());
        }
        if !o.qc.lower.is_empty() {
            notes.push("Supplied lower/upper bounds".into());
        } else if !o.qc.pm.is_empty() {
            notes.push("Supplied ± widths".into());
        }
        if o.group_by.is_some() {
            notes.push("Grouped series · no aggregation".into());
        }
        p.note = notes.join(" · ");
        if matches!(p.content, Content::Heat { .. }) {
            if let Some(x) = o.xlim {
                p.bounds.x = x;
            }
            ensure!(o.qc.ylimits.len() <= 1, "heatmap accepts one yaxis range");
            if let Some(y) = o.qc.ylimits.first() {
                p.bounds.y = [y[0].unwrap_or(p.bounds.y[0]), y[1].unwrap_or(p.bounds.y[1])];
            }
        }
        ensure!(p.bounds.valid(), "plot range overflow or invalid bounds");
    }
    presentation(
        Figure {
            caption: String::new(),
            title: o.title.clone().unwrap_or("1dPetePlot".into()),
            subtitle: format!(
                "{} · {} input file{}",
                o.kind,
                o.files.len(),
                if o.files.len() == 1 { "" } else { "s" }
            ),
            plots,
            columns: 1,
            synthetic: false,
        },
        o,
    )
}

// Titles refer to original panels, before applying the display permutation.
pub(crate) fn presentation(mut f: Figure, o: &Options) -> Result<Figure> {
    f.caption = o.caption.clone();
    for (id, p) in f.plots.iter_mut().enumerate() {
        p.id = id;
    }
    for (n, caption) in &o.panel_captions {
        ensure!(
            *n > 0 && *n <= f.plots.len(),
            "caption panel {n} does not exist"
        );
        f.plots[*n - 1].caption = caption.clone();
    }
    for (n, title) in &o.panel_titles {
        ensure!(
            *n > 0 && *n <= f.plots.len(),
            "panel {n} does not exist ({} panels)",
            f.plots.len()
        );
        f.plots[*n - 1].title = title.clone();
    }
    if o.reverse {
        f.plots.reverse();
    }
    if let Some(order) = &o.panel_order {
        let mut sorted = order.clone();
        sorted.sort_unstable();
        ensure!(
            sorted == (1..=f.plots.len()).collect::<Vec<_>>(),
            "-panel-order must list every panel exactly once, using 1..{}",
            f.plots.len()
        );
        f.plots = order.iter().map(|n| f.plots[*n - 1].clone()).collect();
    }
    if let Some(n) = o.columns {
        f.columns = n.min(f.plots.len());
    }
    Ok(f)
}

pub fn demo() -> Figure {
    let n = 120;
    let mut a = Series {
        label: "Model A".into(),
        color: COLORS[0].into(),
        points: vec![],
        band: Some(vec![]),
        point_labels: vec![],
    };
    let mut b = Series {
        label: "Model B".into(),
        color: COLORS[1].into(),
        points: vec![],
        band: None,
        point_labels: vec![],
    };
    for i in 0..n {
        let t = i as f64;
        let y = 0.12 + 0.5 * (-((t - 60.) / 26.).powi(2)).exp() + 0.035 * (t * 0.23).sin();
        a.points.push([t, y]);
        a.band.as_mut().unwrap().push([y - 0.06, y + 0.06]);
        b.points.push([
            t,
            0.1 + 0.34 * (-((t - 72.) / 23.).powi(2)).exp() + 0.03 * (t * 0.19).cos(),
        ]);
    }
    let line = Plot {
        style: Default::default(),
        id: 0,
        caption: String::new(),
        title: "Model evidence through time".into(),
        note: "Mean ± illustrative uncertainty".into(),
        xlabel: "Time (s)".into(),
        ylabel: "Association (r)".into(),
        bounds: Bounds {
            x: [0., 119.],
            y: [0., 0.8],
        },
        content: Content::Lines(vec![a, b]),
        diagonal: false,
        refs: vec![],
        legend: true,
        qc: Default::default(),
    };
    // Deterministic inverse-CDF style samples; synthetic values, not real RSA results.
    let vals = (0..240)
        .map(|i| {
            let u = (i as f64 + 0.5) / 240.;
            0.24 + 0.085 * (u / (1. - u)).ln() + 0.025 * (i as f64 * 2.4).sin()
        })
        .collect::<Vec<_>>();
    let (e, p) = histogram(&vals, Some(22), None, false).unwrap();
    let ymax = p.iter().map(|p| p[1]).fold(0., f64::max) * 1.18;
    let hist = Plot {
        style: Default::default(),
        id: 1,
        caption: String::new(),
        title: "Effects across regions".into(),
        note: "240 synthetic regions · 22 bins".into(),
        xlabel: "RSA association (r)".into(),
        ylabel: "Regions".into(),
        bounds: Bounds {
            x: [-0.35, 0.85],
            y: [0., ymax],
        },
        content: Content::Histogram {
            series: vec![Series {
                label: "ROI effects".into(),
                color: COLORS[0].into(),
                points: p,
                band: None,
                point_labels: vec![],
            }],
            edges: vec![e],
            density: false,
        },
        diagonal: false,
        refs: vec![(0., None), (0.24, Some("Mean 0.24".into()))],
        legend: false,
        qc: Default::default(),
    };
    let s = Series {
        label: "Regions".into(),
        color: COLORS[2].into(),
        points: (0..90)
            .map(|i| {
                let x = 0.03 + i as f64 / 120.;
                [x, 0.06 + 0.83 * x + 0.075 * (i as f64 * 2.3).sin()]
            })
            .collect(),
        band: None,
        point_labels: vec![],
    };
    let scatter = Plot {
        style: Default::default(),
        id: 2,
        caption: String::new(),
        title: "Compare model associations".into(),
        note: "90 synthetic regions · identity reference".into(),
        xlabel: "Model A (r)".into(),
        ylabel: "Model B (r)".into(),
        bounds: Bounds {
            x: [-0.05, 0.85],
            y: [-0.05, 0.85],
        },
        content: Content::Scatter(vec![s]),
        diagonal: true,
        refs: vec![],
        legend: false,
        qc: Default::default(),
    };
    let labels: Vec<String> = [
        "Face A", "Face B", "Face C", "Place A", "Place B", "Place C", "Object A", "Object B",
        "Object C",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    let matrix = (0..9)
        .map(|r| {
            (0..9)
                .map(|c| {
                    if r == c {
                        0.
                    } else {
                        let d = (r as f64 - c as f64).abs();
                        if r / 3 == c / 3 {
                            0.08 + 0.045 * d
                        } else {
                            0.48 + 0.12 * (d * 0.7).sin().abs()
                        }
                    }
                })
                .collect()
        })
        .collect();
    let heat = Plot {
        style: Default::default(),
        id: 3,
        caption: String::new(),
        title: "Representational structure".into(),
        note: "9 conditions · shared color scale".into(),
        xlabel: String::new(),
        ylabel: String::new(),
        bounds: Bounds {
            x: [0., 9.],
            y: [0., 9.],
        },
        content: Content::Heat {
            matrix,
            row_labels: labels.clone(),
            labels,
            square_cells: true,
            named_conditions: false,
            cmap: "viridis".into(),
            clim: [0., 0.7],
        },
        diagonal: false,
        refs: vec![],
        legend: false,
        qc: Default::default(),
    };
    Figure {
        caption: String::new(),
        title: "A closer look at representational structure".into(),
        subtitle: "RSA EXPLORER  /  Synthetic demonstration data".into(),
        plots: vec![line, hist, scatter, heat],
        columns: 2,
        synthetic: true,
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn panel_titles_bind_before_arrangement() {
        let options = crate::cli::parse(
            &[
                "-demo",
                "-panel-title",
                "1",
                "First source",
                "-panel-order",
                "4,1,2,3",
                "-columns",
                "3",
            ]
            .map(String::from),
        )
        .unwrap();
        let original = super::demo();
        let f = super::from_options(&options).unwrap();
        assert_eq!(f.columns, 3);
        assert_eq!(f.plots[0].title, original.plots[3].title);
        assert_eq!(f.plots[1].title, "First source");
        assert_eq!(f.plots[1].bounds, original.plots[0].bounds);
    }

    use super::*;
    #[test]
    fn heat_reordering_preserves_orientation_and_labels() {
        let dir = tempfile::tempdir().unwrap();
        let matrix = dir.path().join("matrix.1D");
        let order = dir.path().join("order.1D");
        std::fs::write(&matrix, "10 20\n30 40\n").unwrap();
        std::fs::write(&order, "1\n0\n").unwrap();
        let options = Options {
            files: vec![matrix.to_str().unwrap().into()],
            kind: "heat".into(),
            reorder: Some(order.to_str().unwrap().into()),
            zerocenter: true,
            ..Options::default()
        };
        let figure = from_options(&options).unwrap();
        match &figure.plots[0].content {
            Content::Heat {
                matrix,
                labels,
                clim,
                ..
            } => {
                assert_eq!(*matrix, vec![vec![40., 30.], vec![20., 10.]]);
                assert_eq!(*labels, vec!["1", "0"]);
                assert_eq!(*clim, [-40., 40.]);
            }
            _ => panic!("expected heatmap"),
        }
    }
    #[test]
    fn histogram_boundaries_and_density() {
        let (e, p) = histogram(&[0., 1., 2.], Some(2), None, false).unwrap();
        assert_eq!(e, vec![0., 1., 2.]);
        assert_eq!(p.iter().map(|p| p[1]).collect::<Vec<_>>(), vec![1., 2.]);
        let (e, p) = histogram(&[-2., 0., 2.], None, Some(1.5), true).unwrap();
        let integral: f64 = p
            .iter()
            .enumerate()
            .map(|(i, p)| p[1] * (e[i + 1] - e[i]))
            .sum();
        assert!((integral - 1.).abs() < 1e-12);
    }
    #[test]
    fn constant_histogram() {
        let (_, p) = histogram(&[3.; 9], Some(8), None, false).unwrap();
        assert_eq!(p.iter().map(|p| p[1]).sum::<f64>(), 9.);
    }
    #[test]
    fn navigation_preserves_ranges() {
        let b = Bounds {
            x: [0., 10.],
            y: [-1., 1.],
        };
        assert_eq!(
            b.zoom(5., 0., 0.5),
            Bounds {
                x: [2.5, 7.5],
                y: [-0.5, 0.5]
            }
        );
        assert_eq!(
            b.pan(1., 2.),
            Bounds {
                x: [1., 11.],
                y: [1., 3.]
            }
        );
    }
    #[test]
    fn uncertainty_autorange_is_opt_in_and_panel_options_follow_reordering() {
        let args = [
            "-infiles",
            "1D: 1 2 3 | 3 2 1",
            "-yfiles_pm",
            "1D: 3@10 | 3@20",
            "-censor_trs",
            "1",
            "-boxplot_on",
            "-yaxis",
            "0:5",
            ":10",
            "-panel-order",
            "2,1",
        ];
        let o = crate::cli::parse(&args.map(String::from)).unwrap();
        let f = from_options(&o).unwrap();
        assert_eq!(f.plots[1].bounds.y, [0., 5.]);
        assert_eq!(f.plots[0].bounds.y[1], 10.);
        assert_eq!(f.plots[0].qc.boxes[1].n, 2);
        assert!(f.plots[0].hover(1., 2.).contains("CENSORED"));
        let args = ["-infiles", "1D: 1 2 3", "-yfiles_pm", "1D: 3@10"];
        let mut o = crate::cli::parse(&args.map(String::from)).unwrap();
        assert!(from_options(&o).unwrap().plots[0].bounds.y[1] < 4.);
        o.qc.use_pm = true;
        assert!(from_options(&o).unwrap().plots[0].bounds.y[1] > 13.);
    }
    #[test]
    fn independent_shading_does_not_censor_or_change_limits_or_boxplots() {
        let args = ["-infiles", "1D: 0 1 2 3", "-boxplot_on", "-censor_trs", "1"];
        let o = crate::cli::parse(&args.map(String::from)).unwrap();
        let base = from_options(&o).unwrap();
        let mut args = args.map(String::from).to_vec();
        args.extend(
            [
                "-shade",
                "-10:20",
                "-shade-color",
                "blue",
                "-shade-label",
                "Task A",
                "-shade",
                "1:2",
                "-shade-color",
                "orange",
                "-shade-label",
                "Task B",
            ]
            .map(String::from),
        );
        let o = crate::cli::parse(&args).unwrap();
        let f = from_options(&o).unwrap();
        assert_eq!(f.plots[0].bounds, base.plots[0].bounds);
        assert_eq!(f.plots[0].qc.censored, base.plots[0].qc.censored);
        assert_eq!(f.plots[0].qc.boxes[1].n, 3);
        assert_eq!(f.plots[0].qc.boxes[1].five, base.plots[0].qc.boxes[1].five);
        assert_eq!(f.plots[0].qc.annotations.len(), 2);
        assert_eq!(f.plots[0].qc.annotations[0].color, color("blue").unwrap());
    }
}
