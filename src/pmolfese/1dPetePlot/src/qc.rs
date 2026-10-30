//! QC annotations and statistics. Censor indices refer to the selected sample rows.
use crate::{
    input,
    model::{Series, color},
};
use anyhow::{Context, Result, ensure};

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Options {
    pub pm: Vec<String>,
    pub lower: Vec<String>,
    pub upper: Vec<String>,
    pub use_pm: bool,
    pub censor_files: Vec<String>,
    pub censor_trs: Vec<String>,
    pub censor_color: Vec<String>,
    pub hlines: Vec<Option<f64>>,
    pub patches: Vec<usize>,
    pub boxplot: bool,
    pub box_view: String,
    pub scales: Vec<String>,
    pub ylimits: Vec<[Option<f64>; 2]>,
    pub margin_off: bool,
    pub label_wrap: Option<usize>,
}
#[derive(Clone, Debug)]
pub struct Span {
    pub x: [f64; 2],
    pub color: String,
    pub alpha: f64,
}
#[derive(Clone, Debug)]
pub struct BoxStats {
    pub label: String,
    pub n: usize,
    pub five: [f64; 5], // lower whisker, Q1, median, Q3, upper whisker
    pub outliers: Vec<f64>,
}
#[derive(Clone, Debug, Default)]
pub struct Decoration {
    pub xticks: Option<Vec<(f64, String)>>,
    pub spans: Vec<Span>,
    pub annotations: Vec<crate::shading::Shade>,
    pub censored: Vec<bool>,
    pub hlines: Vec<f64>,
    pub boxes: Vec<BoxStats>,
    pub margin_off: bool,
    pub label_wrap: Option<usize>,
}
pub fn box_stats(label: &str, values: impl Iterator<Item = f64>) -> BoxStats {
    let mut v: Vec<_> = values.filter(|v| v.is_finite()).collect();
    v.sort_by(f64::total_cmp);
    let n = v.len();
    if n == 0 {
        return BoxStats {
            label: label.into(),
            n,
            five: [0.; 5],
            outliers: vec![],
        };
    }
    let q = |f: f64| {
        let x = (n - 1) as f64 * f;
        let i = x.floor() as usize;
        v[i] + (v[x.ceil() as usize] - v[i]) * x.fract()
    };
    let (q1, med, q3) = (q(0.25), q(0.5), q(0.75));
    let low = q1 - 1.5 * (q3 - q1);
    let high = q3 + 1.5 * (q3 - q1);
    BoxStats {
        label: label.into(),
        n,
        five: [
            *v.iter().find(|&&x| x >= low).unwrap(),
            q1,
            med,
            q3,
            *v.iter().rev().find(|&&x| x <= high).unwrap(),
        ],
        outliers: v.into_iter().filter(|&x| x < low || x > high).collect(),
    }
}
fn rgba(parts: &[String]) -> Result<(String, f64)> {
    if parts.is_empty() {
        return Ok(("#ffb3b3".into(), 0.45));
    }
    let tokens: Vec<_> = parts.iter().flat_map(|s| s.split_whitespace()).collect();
    if tokens.len() <= 2 {
        let alpha = if tokens.len() == 2 {
            tokens[1].parse::<f64>()?
        } else {
            0.45
        };
        ensure!(
            alpha.is_finite() && (0.0..=1.0).contains(&alpha),
            "censor opacity must be 0..1"
        );
        return Ok((color(tokens[0])?, alpha));
    }
    ensure!(
        tokens.len() == 3 || tokens.len() == 4,
        "censor_RGB needs a color [alpha] or R G B [A]"
    );
    let v = tokens
        .iter()
        .map(|s| s.parse::<f64>())
        .collect::<std::result::Result<Vec<_>, _>>()?;
    ensure!(
        v.iter().all(|x| x.is_finite() && (0.0..=1.0).contains(x)),
        "censor RGB(A) values must be 0..1"
    );
    Ok((
        format!(
            "#{:02x}{:02x}{:02x}",
            (v[0] * 255.).round() as u8,
            (v[1] * 255.).round() as u8,
            (v[2] * 255.).round() as u8
        ),
        *v.get(3).unwrap_or(&0.45),
    ))
}
pub fn prepare(ss: &mut [Series], o: &Options, missing: input::Missing) -> Result<Vec<Decoration>> {
    let count = ss.len();
    let n = ss[0].points.len();
    ensure!(
        ss.iter().all(|s| s.points.len() == n),
        "time-series inputs must have equal selected row counts"
    );
    ensure!(
        o.hlines.is_empty() || o.hlines.len() == 1 || o.hlines.len() == count,
        "censor_hline needs one value or one per series"
    );
    ensure!(
        o.scales.is_empty() || o.scales.len() == 1 || o.scales.len() == count,
        "scale needs one value or one per series"
    );
    ensure!(
        o.ylimits.is_empty() || o.ylimits.len() == 1 || o.ylimits.len() == count,
        "yaxis needs one range or one per series"
    );
    let mut mask = vec![false; n];
    for file in &o.censor_files {
        let m = input::column(file)?;
        ensure!(
            m.len() == n,
            "censor mask must match the selected sample count ({n})"
        );
        ensure!(
            m.iter().all(|v| *v == 0. || *v == 1.),
            "censor mask must contain 0 (censored) or 1 (retained)"
        );
        for (i, v) in m.iter().enumerate() {
            mask[i] |= *v == 0.;
        }
    }
    for spec in &o.censor_trs {
        for i in input::indices(spec, n)? {
            mask[i] = true;
        }
    }
    let censor_on = !o.censor_files.is_empty() || !o.censor_trs.is_empty();
    let mut spans = vec![];
    if censor_on || !o.patches.is_empty() {
        let x: Vec<_> = ss[0].points.iter().map(|p| p[0]).collect();
        ensure!(
            x.windows(2).all(|v| v[1] > v[0]),
            "censoring and run patches require increasing x coordinates"
        );
        let mut edges = vec![if n == 1 {
            x[0] - 0.5
        } else {
            x[0] - (x[1] - x[0]) / 2.
        }];
        edges.extend(x.windows(2).map(|v| v[0] + (v[1] - v[0]) / 2.));
        edges.push(if n == 1 {
            x[0] + 0.5
        } else {
            x[n - 1] + (x[n - 1] - x[n - 2]) / 2.
        });
        ensure!(
            edges.iter().all(|v| v.is_finite()),
            "sample boundaries overflow"
        );
        if !o.patches.is_empty() {
            let total = o
                .patches
                .iter()
                .try_fold(0usize, |a, &b| a.checked_add(b))
                .context("run lengths overflow")?;
            ensure!(
                total == n,
                "patch run lengths must sum to selected sample count ({n})"
            );
            let mut start = 0;
            for (i, &length) in o.patches.iter().enumerate() {
                if i % 2 == 1 && length > 0 {
                    spans.push(Span {
                        x: [edges[start], edges[start + length]],
                        color: "#dce3e8".into(),
                        alpha: 0.5,
                    });
                }
                start += length;
            }
        }
        let (color, alpha) = rgba(&o.censor_color)?;
        let mut i = 0;
        while i < n {
            if !mask[i] {
                i += 1;
                continue;
            }
            let start = i;
            while i < n && mask[i] {
                i += 1;
            }
            spans.push(Span {
                x: [edges[start], edges[i]],
                color: color.clone(),
                alpha,
            });
        }
    } else {
        rgba(&o.censor_color)?;
    }
    let mut pm = vec![];
    for file in &o.pm {
        let m = input::read_with(file, missing)?;
        ensure!(
            m.len() == n,
            "yfiles_pm must match the selected sample count ({n})"
        );
        for c in 0..m[0].len() {
            pm.push(m.iter().map(|r| r[c].abs()).collect::<Vec<_>>());
        }
    }
    ensure!(
        pm.is_empty() || pm.len() == count,
        "yfiles_pm must match the number of selected y columns"
    );
    let endpoints = |files: &[String]| -> Result<Vec<Vec<f64>>> {
        let mut columns = vec![];
        for file in files {
            let m = input::read_with(file, missing)?;
            ensure!(
                m.len() == n,
                "uncertainty bounds must match the selected sample count ({n})"
            );
            for c in 0..m[0].len() {
                columns.push(m.iter().map(|r| r[c]).collect());
            }
        }
        ensure!(
            columns.is_empty() || columns.len() == count,
            "uncertainty bounds must match selected Y columns"
        );
        Ok(columns)
    };
    let lower = endpoints(&o.lower)?;
    let upper = endpoints(&o.upper)?;
    ensure!(
        lower.len() == upper.len(),
        "supply matching lower and upper bounds"
    );
    let mut out = vec![];
    for (i, s) in ss.iter_mut().enumerate() {
        let threshold = o
            .hlines
            .get(if o.hlines.len() == 1 { 0 } else { i })
            .copied()
            .flatten();
        let scale = if let Some(scale) = o.scales.get(if o.scales.len() == 1 { 0 } else { i }) {
            match scale.as_str() {
                "SCALE_TO_HLINE" => {
                    let h = threshold
                        .context("SCALE_TO_HLINE requires a threshold for every series")?;
                    ensure!(h != 0., "cannot scale to a zero threshold");
                    1. / h
                }
                "SCALE_TO_MAX" => {
                    let max = s
                        .points
                        .iter()
                        .map(|p| p[1])
                        .filter(|v| v.is_finite())
                        .fold(f64::NEG_INFINITY, f64::max);
                    ensure!(
                        max.is_finite() && max != 0.,
                        "cannot scale to a zero or missing maximum"
                    );
                    1. / max
                }
                _ => scale
                    .parse::<f64>()
                    .context("scale must be numeric, SCALE_TO_HLINE or SCALE_TO_MAX")?,
            }
        } else {
            1.
        };
        ensure!(scale.is_finite(), "scale must be finite");
        let supplied = if !lower.is_empty() {
            let mut band = vec![];
            for (row, p) in s.points.iter().enumerate() {
                let (lo, hi) = (lower[i][row], upper[i][row]);
                ensure!(
                    !(lo.is_finite() && hi.is_finite()) || lo <= hi,
                    "lower uncertainty endpoint exceeds upper at sample {row}"
                );
                ensure!(
                    !p[1].is_finite()
                        || ((!lo.is_finite() || lo <= p[1]) && (!hi.is_finite() || p[1] <= hi)),
                    "uncertainty bounds do not contain the estimate at sample {row}"
                );
                let pair = if !lo.is_finite() || !hi.is_finite() {
                    [f64::NAN; 2]
                } else if scale >= 0. {
                    [lo * scale, hi * scale]
                } else {
                    [hi * scale, lo * scale]
                };
                ensure!(
                    pair.iter().all(|v| !v.is_infinite()),
                    "scaled uncertainty bounds overflow"
                );
                band.push(pair);
            }
            Some(band)
        } else {
            None
        };
        if supplied.is_some() {
            s.band = supplied;
        }
        for p in &mut s.points {
            let finite = p[1].is_finite();
            p[1] *= scale;
            ensure!(!finite || p[1].is_finite(), "scaled value overflow");
        }
        if let Some(widths) = pm.get(i) {
            let b: Vec<_> = s
                .points
                .iter()
                .zip(widths)
                .map(|(p, w)| [p[1] - w * scale.abs(), p[1] + w * scale.abs()])
                .collect();
            ensure!(
                b.iter().flatten().all(|v| !v.is_infinite()),
                "uncertainty band overflow"
            );
            s.band = Some(b);
        }
        let mut d = Decoration {
            spans: spans.clone(),
            censored: mask.clone(),
            margin_off: o.margin_off,
            label_wrap: o.label_wrap,
            ..Default::default()
        };
        if let Some(h) = threshold {
            let h = h * scale;
            ensure!(h.is_finite(), "scaled threshold overflow");
            d.hlines.push(h);
        }
        if o.boxplot {
            if !censor_on || o.box_view != "AC_ONLY" {
                d.boxes.push(box_stats(
                    if censor_on { "BC" } else { "All" },
                    s.points.iter().map(|p| p[1]),
                ));
            }
            if censor_on && o.box_view != "BC_ONLY" {
                d.boxes.push(box_stats(
                    "AC",
                    s.points
                        .iter()
                        .enumerate()
                        .filter(|(j, _)| !mask[*j])
                        .map(|(_, p)| p[1]),
                ));
            }
        }
        out.push(d);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn series(y: &[f64]) -> Series {
        Series {
            label: "test".into(),
            color: "#123456".into(),
            points: y
                .iter()
                .enumerate()
                .map(|(i, v)| [i as f64 * 2., *v])
                .collect(),
            band: None,
            point_labels: vec![],
        }
    }
    #[test]
    fn box_statistics_use_linear_percentiles_and_observed_whiskers() {
        let b = box_stats("BC", [0., 1., 2., 3., 100., f64::NAN].into_iter());
        assert_eq!(b.n, 5);
        assert_eq!(b.five, [0., 1., 2., 3., 3.]);
        assert_eq!(b.outliers, [100.]);
        let b = box_stats("AC", [0., 1., 2., 3.].into_iter());
        assert_eq!(b.five, [0., 0.75, 1.5, 2.25, 3.]);
        assert_eq!(box_stats("AC", [].into_iter()).n, 0);
    }
    #[test]
    fn masks_union_without_deleting_samples_and_follow_x_spacing() {
        let d = tempfile::tempdir().unwrap();
        let f = d.path().join("mask.1D");
        std::fs::write(&f, "1\n0\n1\n1\n").unwrap();
        let o = Options {
            censor_files: vec![f.to_str().unwrap().into()],
            censor_trs: vec!["2".into()],
            patches: vec![2, 2],
            boxplot: true,
            ..Default::default()
        };
        let mut ss = vec![series(&[0., 1., 2., 100.])];
        let ds = prepare(&mut ss, &o, input::Missing::Error).unwrap();
        let d = &ds[0];
        assert_eq!(ss[0].points.len(), 4);
        assert_eq!(d.censored, [false, true, true, false]);
        assert_eq!(d.spans[0].x, [3., 7.]);
        assert_eq!(d.spans[1].x, [1., 5.]);
        assert_eq!(d.boxes[0].n, 4);
        assert_eq!(d.boxes[1].n, 2);
        assert_eq!(d.boxes[1].five[2], 50.);
    }
    #[test]
    fn scales_apply_to_curve_band_and_threshold_but_not_axis_limits() {
        let o = Options {
            pm: vec!["1D: 4@0.5".into()],
            scales: vec!["-2".into()],
            hlines: vec![Some(1.)],
            ylimits: vec![[Some(0.), Some(4.)]],
            ..Default::default()
        };
        let mut ss = vec![series(&[0., 1., 2., 3.])];
        let ds = prepare(&mut ss, &o, input::Missing::Error).unwrap();
        assert_eq!(ss[0].points[1][1], -2.);
        assert_eq!(ss[0].band.as_ref().unwrap()[1], [-3., -1.]);
        assert_eq!(ds[0].hlines, [-2.]);
        assert_eq!(o.ylimits[0], [Some(0.), Some(4.)]);
        let o = Options {
            scales: vec!["SCALE_TO_MAX".into()],
            ..Default::default()
        };
        let mut ss = vec![series(&[-3., -1.])];
        prepare(&mut ss, &o, input::Missing::Error).unwrap();
        assert_eq!(ss[0].points[0][1], 3.); // maximum, not maximum absolute value
    }
    #[test]
    fn invalid_masks_and_zero_denominators_fail() {
        for o in [
            Options {
                censor_files: vec!["1D: 1 0".into()],
                ..Default::default()
            },
            Options {
                censor_files: vec!["1D: 1 2 1".into()],
                ..Default::default()
            },
            Options {
                patches: vec![2, 2],
                ..Default::default()
            },
            Options {
                scales: vec!["SCALE_TO_HLINE".into()],
                hlines: vec![Some(0.)],
                ..Default::default()
            },
            Options {
                pm: vec!["1D: 1 1".into()],
                ..Default::default()
            },
        ] {
            assert!(prepare(&mut [series(&[0., 1., 2.])], &o, input::Missing::Error).is_err());
        }
        let o = Options {
            censor_trs: vec!["0..$".into()],
            boxplot: true,
            box_view: "AC_ONLY".into(),
            ..Default::default()
        };
        let d = prepare(&mut [series(&[0., 1., 2.])], &o, input::Missing::Error).unwrap();
        assert_eq!(d[0].boxes.len(), 1);
        assert_eq!(d[0].boxes[0].n, 0);
    }
    #[test]
    fn supplied_bounds_preserve_asymmetry_and_reverse_correctly_under_negative_scaling() {
        let options = crate::cli::parse(
            &[
                "-infiles",
                "1D: 1 3 5",
                "-yfiles_lo",
                "1D: 0 1 4",
                "-yfiles_hi",
                "1D: 4 5 6",
                "-scale",
                "-2",
                "-ylim_use_pm",
            ]
            .map(str::to_owned),
        )
        .unwrap();
        let fig = crate::model::from_options(&options).unwrap();
        let crate::model::Content::Lines(s) = &fig.plots[0].content else {
            panic!()
        };
        assert_eq!(
            s[0].band.as_ref().unwrap(),
            &vec![[-8., 0.], [-10., -2.], [-12., -8.]]
        );
        assert_eq!(s[0].points[1][1], -6.);
        assert!(fig.plots[0].bounds.y[0] < -12.);
    }
}
