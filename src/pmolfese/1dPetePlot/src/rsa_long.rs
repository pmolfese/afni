//! Identity-based readers for the two long schemas emitted by 3dRSA.
//! Time labels are annotations, never inferred timestamps; run summaries are categorical.
use crate::{
    cli::Options,
    input::Missing,
    model::{self, Bounds, Content, Figure, IntervalRow, Series},
    rsa::{self, Roi, Table},
};
use anyhow::{Context, Result, ensure};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone, Debug)]
pub struct LongTable {
    pub time: bool,
    /// Index strings for time tables; literal summary names for run tables.
    pub steps: Vec<String>,
    pub labels: BTreeMap<String, String>,
    pub rows: BTreeMap<(i32, String, String), Vec<f64>>,
    pub preview: Vec<Vec<String>>,
}
pub fn detect(text: &str) -> bool {
    text.lines()
        .filter_map(|l| l.trim().strip_prefix('#'))
        .any(|l| {
            let f: Vec<_> = l.split_whitespace().collect();
            f.first() == Some(&"ROI")
                && f.iter()
                    .any(|s| ["summary", "time_index", "time_label"].contains(s))
        })
}
fn token(s: &str) -> bool {
    !s.is_empty() && s.len() <= 128 && !s.contains(',') && !s.starts_with('-')
}
fn number(s: &str) -> Result<f64> {
    let v = if ["NA", "N/A", "NAN"].contains(&s.to_ascii_uppercase().as_str()) {
        f64::NAN
    } else {
        s.parse().context("invalid numeric RSA value")?
    };
    ensure!(!v.is_infinite(), "infinite RSA value");
    Ok(v)
}
pub fn read(text: &str) -> Result<Table> {
    let mut table = Table {
        long: None,
        columns: vec![],
        models: vec![],
        effects: BTreeMap::new(),
        rows: BTreeMap::new(),
        metadata: vec![],
    };
    let mut long = LongTable {
        time: false,
        steps: vec![],
        labels: BTreeMap::new(),
        rows: BTreeMap::new(),
        preview: vec![],
    };
    for (line_no, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(comment) = line.strip_prefix('#') {
            let f: Vec<_> = comment.split_whitespace().collect();
            if f.first() != Some(&"ROI") {
                if table.metadata.len() < 100 {
                    table.metadata.push(comment.trim().into());
                }
                continue;
            }
            ensure!(table.columns.is_empty(), "multiple RSA column headers");
            ensure!(
                f.len() >= 9 && f[..3] == ["ROI", "label", "nvox"],
                "incomplete long RSA header"
            );
            long.time = f[3..5] == ["time_index", "time_label"];
            ensure!(
                long.time || f[3..5] == ["summary", "model"],
                "unknown long RSA identity columns"
            );
            ensure!(
                f.len() <= 512
                    && f.iter().all(|s| s.len() <= 128)
                    && f.iter().collect::<BTreeSet<_>>().len() == f.len(),
                "invalid/duplicate long RSA columns"
            );
            let required = if long.time {
                vec!["effect", "stat", "p", "q", "z"]
            } else if f.contains(&"beta") {
                vec!["beta", "partial_r", "p", "q"]
            } else {
                vec!["effect", "p", "q"]
            };
            ensure!(
                required.iter().all(|s| f.contains(s)),
                "incomplete long RSA schema"
            );
            if !long.time {
                ensure!(
                    f.contains(&"z") != f.contains(&"uncalibrated"),
                    "run table requires exactly one z/uncalibrated column"
                );
            }
            for c in &f[5..] {
                let allowed = if long.time {
                    [
                        "effect",
                        "stat",
                        "p",
                        "q",
                        "pfwe",
                        "z",
                        "zfwe",
                        "bootLo",
                        "bootHi",
                        "cbootLo",
                        "cbootHi",
                        "dualLo",
                        "dualHi",
                        "reliability",
                        "nc_low",
                        "nc_high",
                    ]
                    .contains(c)
                } else {
                    [
                        "effect",
                        "beta",
                        "partial_r",
                        "p",
                        "q",
                        "pfwe",
                        "z",
                        "zfwe",
                        "uncalibrated",
                    ]
                    .contains(c)
                };
                ensure!(allowed, "unsupported long RSA column '{c}'");
            }
            ensure!(
                long.time
                    || (f.contains(&"effect") != f.contains(&"beta"))
                        && (f.contains(&"beta") == f.contains(&"partial_r")),
                "mixed/incomplete run effect roles"
            );
            table.columns = f.iter().map(|s| s.to_string()).collect();
            continue;
        }
        let n = table.columns.len();
        ensure!(n > 0, "RSA data precedes header");
        let fields: Vec<_> = line.split_whitespace().collect();
        ensure!(
            fields.len() >= n,
            "long RSA row {} has too few fields",
            line_no + 1
        );
        // AFNI ROI labels may have spaces; all trailing schema fields have fixed roles.
        let start = fields.len() - (n - 2);
        let id: i32 = fields[0].parse().context("invalid ROI ID")?;
        ensure!(id > 0, "ROI IDs must be positive");
        let label = fields[1..start].join(" ");
        ensure!(!label.is_empty() && label.len() <= 256, "invalid ROI label");
        let nvox: u32 = fields[start].parse().context("invalid nvox")?;
        let (step, step_label, model) = if long.time {
            let index: u32 = fields[start + 1]
                .parse()
                .context("time_index must be a nonnegative integer")?;
            (
                index.to_string(),
                fields[start + 2].to_string(),
                "series".to_string(),
            )
        } else {
            (
                fields[start + 1].to_string(),
                fields[start + 1].to_string(),
                fields[start + 2].to_string(),
            )
        };
        ensure!(
            token(&step)
                && (if long.time {
                    !step_label.is_empty() && step_label.len() <= 511
                } else {
                    token(&step_label)
                })
                && token(&model),
            "model/summary names must be 1..128 bytes without commas or a leading dash; time labels may contain up to 511 bytes"
        );
        if let Some(old) = long.labels.get(&step) {
            ensure!(
                old == &step_label,
                "time_index {step} has conflicting time labels"
            );
        } else {
            ensure!(
                !long.time || !long.labels.values().any(|label| label == &step_label),
                "duplicate time_label across different indices"
            );
            long.steps.push(step.clone());
            long.labels.insert(step.clone(), step_label);
        }
        if !table.models.contains(&model) {
            table.models.push(model.clone());
        }
        let mut values = vec![f64::NAN; 5];
        for (column, field) in table.columns[5..].iter().zip(&fields[start + 3..]) {
            let value = number(field)?;
            if ["p", "q", "pfwe"].contains(&column.as_str()) {
                ensure!(
                    value.is_nan() || (0.0..=1.0).contains(&value),
                    "{column} must be in 0..1 or missing"
                );
            }
            values.push(value);
        }
        if let Some(old) = table.rows.get(&id) {
            ensure!(
                old.label == label && old.nvox == nvox,
                "ROI {id} has conflicting label/nvox metadata"
            );
        } else {
            table.rows.insert(
                id,
                Roi {
                    id,
                    label: label.clone(),
                    nvox,
                    values: values.clone(),
                },
            );
        }
        ensure!(
            long.rows.insert((id, model, step), values).is_none(),
            "duplicate ROI/model/time or ROI/model/summary identity on row {}",
            line_no + 1
        );
        ensure!(
            table.rows.len() <= 10_000 && long.rows.len() * n <= 1_000_000,
            "long RSA table exceeds 10000 ROIs / 1000000 cells"
        );
        if long.preview.len() < 8 {
            let mut row = vec![id.to_string(), label];
            row.extend(fields[start..].iter().map(|s| s.to_string()));
            long.preview.push(row);
        }
    }
    ensure!(!long.rows.is_empty(), "RSA input contains no rows");
    if long.time {
        long.steps.sort_by_key(|s| s.parse::<u32>().unwrap());
    }
    for (i, c) in table.columns.iter().enumerate().skip(5) {
        let kind = match c.as_str() {
            "effect" => Some("effect"),
            "beta" => Some("b"),
            "partial_r" => Some("pr"),
            _ => None,
        };
        if let Some(kind) = kind {
            for m in &table.models {
                table.effects.insert((m.clone(), kind.into()), i);
            }
        }
    }
    for method in ["boot", "cboot", "dual"] {
        let pair = ci_columns(&table, method);
        let lo = table.columns.contains(&format!("{method}Lo"));
        let hi = table.columns.contains(&format!("{method}Hi"));
        ensure!(lo == hi, "incomplete {method} CI bounds");
        if let Some([i, j]) = pair {
            for v in long.rows.values() {
                ensure!(
                    v[i].is_nan() || v[j].is_nan() || v[i] <= v[j],
                    "reversed {method} CI endpoints"
                );
            }
        }
    }
    table.long = Some(Box::new(long));
    Ok(table)
}
fn ci_columns(table: &Table, method: &str) -> Option<[usize; 2]> {
    Some([
        table
            .columns
            .iter()
            .position(|s| s == &format!("{method}Lo"))?,
        table
            .columns
            .iter()
            .position(|s| s == &format!("{method}Hi"))?,
    ])
}
struct Selection {
    models: Vec<String>,
    rois: Vec<i32>,
    steps: Vec<String>,
    effect: String,
    modes: Vec<String>,
    ci: Option<String>,
}
fn select(o: &Options, t: &Table, l: &LongTable) -> Result<Selection> {
    let r = &o.rsa;
    ensure!(
        r.version == 1 && r.pair.is_none(),
        "long tables use preset v1; -rsa-pair is for wide tables"
    );
    ensure!(o.missing != Missing::Mean, "RSA effects are never imputed");
    ensure!(
        o.bins.is_none() && o.binwidth.is_none() && !o.density,
        "long RSA tables do not pool observations into histograms"
    );
    let models = r.models.clone().unwrap_or_else(|| t.models.clone());
    ensure!(
        !models.is_empty() && models.len() <= 8 && models.iter().all(|m| t.models.contains(m)),
        "choose up to 8 models from {}",
        t.models.join(",")
    );
    let rois = r
        .rois
        .clone()
        .unwrap_or_else(|| t.rows.keys().copied().collect());
    ensure!(
        !rois.is_empty()
            && rois.len() <= 40
            && rois.iter().all(|id| t.rows.contains_key(id))
            && rois.iter().collect::<BTreeSet<_>>().len() == rois.len(),
        "choose up to 40 unique existing ROI IDs with -rsa-rois"
    );
    ensure!(
        if l.time {
            r.summaries.is_none()
        } else {
            r.times.is_none()
        },
        "-rsa-times requires a time table; -rsa-summaries requires a run table"
    );
    let steps = if l.time {
        if let Some(times) = &r.times {
            let times: Vec<_> = times
                .iter()
                .map(|s| s.parse::<u32>().map(|i| i.to_string()))
                .collect::<std::result::Result<_, _>>()
                .context("invalid time_index selection")?;
            ensure!(
                times.iter().all(|s| l.steps.contains(s))
                    && times.iter().collect::<BTreeSet<_>>().len() == times.len(),
                "unknown or duplicate time_index"
            );
            l.steps
                .iter()
                .filter(|s| times.contains(s))
                .cloned()
                .collect()
        } else {
            l.steps.clone()
        }
    } else {
        r.summaries.clone().unwrap_or_else(|| l.steps.clone())
    };
    ensure!(
        !steps.is_empty()
            && steps.len() <= 2000
            && steps.iter().all(|s| l.steps.contains(s))
            && steps.iter().collect::<BTreeSet<_>>().len() == steps.len(),
        "choose 1..2000 existing time indices or summaries"
    );
    ensure!(
        steps.len() * rois.len() * models.len() <= 50_000,
        "select at most 50000 long RSA observations"
    );
    let effect = if r.effect == "auto" {
        if t.columns.iter().any(|c| c == "beta") {
            "b"
        } else {
            "effect"
        }
    } else {
        &r.effect
    }
    .to_string();
    ensure!(
        models
            .iter()
            .all(|m| t.effects.contains_key(&(m.clone(), effect.clone()))),
        "unavailable effect kind; long tables use effect, or b/pr for joint run regression (stat is not an effect)"
    );
    let modes = r
        .plots
        .clone()
        .unwrap_or_else(|| vec![if l.time { "curves" } else { "points" }.into()]);
    ensure!(
        !modes.is_empty()
            && modes.iter().all(|s| [
                "heatmap",
                "intervals",
                if l.time { "curves" } else { "points" }
            ]
            .contains(&s.as_str())),
        "time tables support curves,heatmap,intervals; run tables support points,heatmap"
    );
    ensure!(
        modes.len() * rois.len() <= 64,
        "select at most 64 plot panels"
    );
    ensure!(
        !modes.iter().any(|s| s == "heatmap") || steps.len() <= 40,
        "heatmaps support at most 40 time/summary columns; use a selection"
    );
    let wants_ci = modes.iter().any(|s| s == "intervals") || r.ci != "auto";
    let ci = if wants_ci {
        ensure!(
            l.time && effect == "effect",
            "supplied CI endpoints are available only for time-series effect tables"
        );
        let methods: Vec<_> = ["boot", "cboot", "dual"]
            .into_iter()
            .filter(|m| ci_columns(t, m).is_some())
            .collect();
        let method = if r.ci == "auto" {
            ensure!(
                methods.len() == 1,
                "choose -rsa-ci boot|cboot|dual: automatic CI selection requires exactly one supplied method"
            );
            methods[0]
        } else {
            ensure!(
                methods.contains(&r.ci.as_str()),
                "selected CI method is unavailable"
            );
            &r.ci
        };
        ensure!(
            !modes.iter().any(|s| s == "intervals") || steps.len() * models.len() <= 40,
            "interval plots support at most 40 rows per ROI"
        );
        ensure!(
            modes
                .iter()
                .any(|s| ["curves", "intervals"].contains(&s.as_str())),
            "-rsa-ci requires curves or intervals"
        );
        Some(method.into())
    } else {
        None
    };
    Ok(Selection {
        models,
        rois,
        steps,
        effect,
        modes,
        ci,
    })
}
pub fn freeze(o: &Options, t: &Table, l: &LongTable) -> Result<Vec<String>> {
    let s = select(o, t, l)?;
    let mut args = o.arguments.clone();
    let mut flags = vec![
        ("-rsa-version", "1".into()),
        ("-rsa-models", s.models.join(",")),
        (
            "-rsa-rois",
            s.rois
                .iter()
                .map(i32::to_string)
                .collect::<Vec<_>>()
                .join(","),
        ),
        ("-rsa-effect", s.effect),
        ("-rsa-plots", s.modes.join(",")),
        (
            if l.time {
                "-rsa-times"
            } else {
                "-rsa-summaries"
            },
            s.steps.join(","),
        ),
    ];
    if let Some(ci) = s.ci {
        flags.push(("-rsa-ci", ci));
    }
    for (flag, value) in flags {
        if let Some(i) = args.iter().position(|s| s == flag) {
            args[i + 1] = value;
        } else {
            args.extend([flag.into(), value]);
        }
    }
    Ok(args)
}
pub fn figure(o: &Options, t: &Table, l: &LongTable) -> Result<Figure> {
    let s = select(o, t, l)?;
    let ci = s.ci.as_deref().and_then(|m| ci_columns(t, m));
    let ci_label =
        s.ci.as_deref()
            .map(|m| rsa::interval_label(t, m))
            .transpose()?;
    let kind = match s.effect.as_str() {
        "b" => "Standardized beta",
        "pr" => "Partial r",
        _ => "RSA effect (source scale)",
    };
    let ticks: Vec<_> = s
        .steps
        .iter()
        .enumerate()
        .map(|(i, step)| {
            (
                if l.time {
                    step.parse::<u32>().unwrap() as f64
                } else {
                    i as f64
                },
                if l.time {
                    format!("{step}: {}", l.labels[step])
                } else {
                    step.clone()
                },
            )
        })
        .collect();
    let mut plots = vec![];
    for id in &s.rois {
        let roi = &t.rows[id];
        let title = format!("{}: {}", id, roi.label);
        let mut series = vec![];
        for (i, m) in s.models.iter().enumerate() {
            let c = t.effects[&(m.clone(), s.effect.clone())];
            let mut points = vec![];
            let mut band = vec![];
            let mut labels = vec![];
            for (j, step) in s.steps.iter().enumerate() {
                let row = l.rows.get(&(*id, m.clone(), step.clone()));
                let effect = row.map(|v| v[c]).unwrap_or(f64::NAN);
                let bounds =
                    ci.map(|[lo, hi]| row.map(|v| [v[lo], v[hi]]).unwrap_or([f64::NAN; 2]));
                ensure!(
                    o.missing != Missing::Error
                        || effect.is_finite()
                            && bounds.is_none_or(|b| b.iter().all(|v| v.is_finite())),
                    "ROI {id}, model {m}, {}: missing observation/bounds; use -missing gap",
                    ticks[j].1
                );
                points.push([ticks[j].0, effect]);
                if let Some(b) = bounds {
                    band.push(if effect.is_finite() { b } else { [f64::NAN; 2] });
                }
                let mut detail = format!(
                    "ROI {id}: {} · nvox {} · {m} · {}",
                    roi.label, roi.nvox, ticks[j].1
                );
                if let Some(row) = row {
                    for (k, col) in t.columns.iter().enumerate().skip(5) {
                        if k != c {
                            detail.push_str(&format!(" · {col}={}", row[k]));
                        }
                    }
                } else {
                    detail.push_str(" · absent row");
                }
                labels.push(detail);
            }
            series.push(Series {
                label: m.clone(),
                color: model::color(
                    o.colors
                        .get(i % o.colors.len().max(1))
                        .map(String::as_str)
                        .unwrap_or(model::COLORS[i % model::COLORS.len()]),
                )?,
                points,
                band: ci.map(|_| band),
                point_labels: labels,
            });
        }
        ensure!(
            series
                .iter()
                .flat_map(|s| &s.points)
                .any(|p| p[1].is_finite()),
            "ROI {id}: no finite selected effects"
        );
        for mode in &s.modes {
            let axis = if l.time {
                "Time index: source label (not seconds)"
            } else {
                "Run / summary (categorical)"
            };
            let mut p = match mode.as_str() {
                "curves" | "points" => {
                    let bounds = model::bounds(&series, ci.is_some());
                    let content = if mode == "curves" {
                        Content::Lines(series.clone())
                    } else {
                        Content::Scatter(series.clone())
                    };
                    let mut p = rsa::plot(content, title.clone(), axis.into(), kind.into(), bounds);
                    p.qc.xticks = Some(ticks.clone());
                    p.legend = !l.time;
                    p
                }
                "heatmap" => {
                    let matrix: Vec<_> = series
                        .iter()
                        .map(|s| s.points.iter().map(|p| p[1]).collect::<Vec<_>>())
                        .collect();
                    let max = matrix
                        .iter()
                        .flatten()
                        .filter(|v| v.is_finite())
                        .map(|v| v.abs())
                        .fold(0., f64::max)
                        .max(1e-12);
                    rsa::plot(
                        Content::Heat {
                            matrix,
                            labels: ticks.iter().map(|(_, s)| s.clone()).collect(),
                            row_labels: s.models.clone(),
                            square_cells: false,
                            named_conditions: false,
                            cmap: "Reds_and_Blues_Inv".into(),
                            clim: [-max, max],
                        },
                        title.clone(),
                        axis.into(),
                        "Model".into(),
                        Bounds {
                            x: [0., s.steps.len() as f64],
                            y: [0., s.models.len() as f64],
                        },
                    )
                }
                "intervals" => {
                    let mut rows = vec![];
                    for ss in &series {
                        for (j, point) in ss.points.iter().enumerate() {
                            let b = ss.band.as_ref().unwrap()[j];
                            let v = [point[1], b[0], b[1]];
                            rows.push(IntervalRow {
                                label: ticks[j].1.clone(),
                                detail: ss.point_labels[j].clone(),
                                values: v.iter().all(|v| v.is_finite()).then_some(v),
                            });
                        }
                    }
                    let finite: Vec<_> = rows.iter().filter_map(|r| r.values).flatten().collect();
                    ensure!(
                        !finite.is_empty(),
                        "ROI {id}: no complete supplied CI triplets"
                    );
                    let range = model::padded(
                        finite.iter().copied().fold(0., f64::min),
                        finite.iter().copied().fold(0., f64::max),
                    );
                    let bounds = Bounds {
                        x: range,
                        y: [0., rows.len() as f64],
                    };
                    rsa::plot(
                        Content::Intervals {
                            rows,
                            color: series[0].color.clone(),
                            method: ci_label.clone().unwrap(),
                        },
                        title.clone(),
                        kind.into(),
                        "Time index: source label".into(),
                        bounds,
                    )
                }
                _ => unreachable!(),
            };
            p.note = if let Some(label) = &ci_label {
                label.clone()
            } else if l.time {
                "Source time indices · full labels and statistics on hover".into()
            } else {
                "MEAN, runs and contrasts are distinct summaries; no connecting lines".into()
            };
            ensure!(p.bounds.valid(), "RSA plot bounds overflow");
            plots.push(p);
        }
    }
    let mut f = model::presentation(
        Figure {
            title: o
                .title
                .clone()
                .unwrap_or_else(|| "RSA · long results".into()),
            subtitle: format!(
                "RSA preset v1 · {kind} · {} ROIs · {} {}",
                s.rois.len(),
                s.steps.len(),
                if l.time { "time points" } else { "summaries" }
            ),
            caption: String::new(),
            columns: 1,
            plots,
            synthetic: false,
        },
        o,
    )?;
    if f.caption.is_empty() {
        f.caption = if l.time { "Effects retain source time_index order and scale. Labels are not interpreted as elapsed time. Supplied CI endpoints only when requested; no inference, filtering or imputation." } else { "Supplied run effects and named summaries. MEAN is 3dRSA's equal-run estimate; contrasts are separate estimands. No averaging, significance filtering or imputation." }.into();
        for meta in t.metadata.iter().filter(|m| {
            m.starts_with("metric")
                || m.starts_with("run contrast:")
                || m.starts_with("model_series")
        }) {
            f.caption.push_str(&format!(" {meta}."));
        }
    }
    Ok(f)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;
    fn fixture(name: &str) -> String {
        format!("{}/examples/rsa/{name}.rsa.1D", env!("CARGO_MANIFEST_DIR"))
    }
    fn options(file: &str, flags: &[&str]) -> Options {
        let mut args = vec!["-rsa".into(), file.into()];
        args.extend(flags.iter().map(|s| s.to_string()));
        crate::cli::parse(&args).unwrap()
    }
    #[test]
    fn time_curves_keep_index_label_roi_and_ci_identity_without_treating_stat_as_effect() {
        let file = fixture("series");
        let t = rsa::read(Path::new(&file)).unwrap();
        assert_eq!(t.effects.len(), 1);
        assert_eq!(t.effects[&("series".into(), "effect".into())], 5);
        let o = options(
            &file,
            &[
                "-rsa-rois",
                "12,3",
                "-rsa-times",
                "10,2,0",
                "-rsa-ci",
                "boot",
            ],
        );
        let f = rsa::figure(&o).unwrap();
        assert_eq!(f.plots.len(), 2);
        assert!(f.plots[0].title.starts_with("12:"));
        let Content::Lines(ss) = &f.plots[0].content else {
            panic!()
        };
        assert_eq!(
            ss[0].points.iter().map(|p| p[0]).collect::<Vec<_>>(),
            [0., 2., 10.]
        );
        assert_eq!(
            f.plots[0].qc.xticks.as_ref().unwrap()[2],
            (10., "10: phase_10".into())
        );
        assert_eq!(ss[0].points[0][1], 0.129376);
        assert_eq!(ss[0].band.as_ref().unwrap()[0], [0.039376, 0.249376]);
        assert!(f.plots[0].hover(0., 0.129376).contains("ROI 12: Occipital"));
        assert!(f.plots[0].hover(0., 0.129376).contains("stat="));
        assert!(f.plots[0].note.contains("95%"));
        let frozen = rsa::freeze(&o).unwrap();
        assert_eq!(
            frozen[frozen.iter().position(|s| s == "-rsa-times").unwrap() + 1],
            "0,2,10"
        );
        for flags in [
            &["-rsa-effect", "stat"][..],
            &["-rsa-models", "missing"],
            &["-rsa-times", "2,999"],
            &["-rsa-summaries", "MEAN"],
            &["-rsa-plots", "histogram"],
        ] {
            assert!(rsa::figure(&options(&file, flags)).is_err());
        }
    }
    #[test]
    fn verbatim_negative_time_labels_are_annotations_and_empty_intervals_fail() {
        let source = "#ROI label nvox time_index time_label effect stat p q z bootLo bootHi\n1 A 2 0 -200,ms 0.2 4 1 1 0 NA NA\n1 A 2 1 100,ms 0.3 5 1 1 0 NA NA\n";
        let t = read(source).unwrap();
        assert_eq!(t.long.as_ref().unwrap().labels["0"], "-200,ms");
        assert!(read(&source.replace("100,ms", "-200,ms")).is_err());
        let d = tempfile::tempdir().unwrap();
        let p = d.path().join("empty.1D");
        std::fs::write(&p, source).unwrap();
        assert!(
            rsa::figure(&options(
                p.to_str().unwrap(),
                &["-rsa-plots", "intervals", "-missing", "gap"]
            ))
            .unwrap_err()
            .to_string()
            .contains("no complete")
        );
    }
    #[test]
    fn runs_are_unconnected_categories_with_distinct_models_beta_and_partial_r() {
        let file = fixture("runs-joint");
        let f = rsa::figure(&options(
            &file,
            &[
                "-rsa-models",
                "visual,semantic",
                "-rsa-rois",
                "3",
                "-rsa-summaries",
                "task,rest,task_minus_rest",
                "-rsa-effect",
                "pr",
            ],
        ))
        .unwrap();
        let Content::Scatter(ss) = &f.plots[0].content else {
            panic!()
        };
        assert_eq!(ss[0].label, "visual");
        assert_eq!(ss[0].points, [[0., 0.55], [1., 0.41], [2., 0.14]]);
        assert!(ss[0].point_labels[2].contains("task_minus_rest"));
        assert!(f.caption.contains("fixed-run estimand"));
        let heat = rsa::figure(&options(
            &file,
            &[
                "-rsa-rois",
                "3",
                "-rsa-models",
                "visual",
                "-rsa-summaries",
                "task,rest",
                "-rsa-plots",
                "heatmap",
            ],
        ))
        .unwrap();
        let Content::Heat {
            matrix,
            labels,
            row_labels,
            ..
        } = &heat.plots[0].content
        else {
            panic!()
        };
        assert_eq!(labels, &["task", "rest"]);
        assert_eq!(row_labels, &["visual"]);
        assert_eq!(matrix, &[vec![1.1, 0.82]]);
        for flags in [
            &["-rsa-ci", "boot"][..],
            &["-rsa-times", "0"],
            &["-rsa-plots", "curves"],
            &["-rsa-effect", "r"],
            &["-rsa-summaries", "missing"],
        ] {
            assert!(rsa::figure(&options(&file, flags)).is_err());
        }
    }
    #[test]
    fn shuffled_rows_join_by_identity_missing_combinations_gap_and_conflicts_fail() {
        let source = std::fs::read_to_string(fixture("series")).unwrap();
        let mut lines: Vec<_> = source
            .lines()
            .filter(|s| !s.starts_with("12 Occipital 210 1 "))
            .collect();
        let header = lines.iter().position(|l| l.starts_with("#ROI")).unwrap();
        lines[header + 1..].reverse();
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("series.1D");
        std::fs::write(&path, lines.join("\n")).unwrap();
        let file = path.to_str().unwrap();
        assert!(rsa::figure(&options(file, &[])).is_err());
        let f = rsa::figure(&options(
            file,
            &["-rsa-rois", "12", "-missing", "gap", "-rsa-ci", "boot"],
        ))
        .unwrap();
        let Content::Lines(ss) = &f.plots[0].content else {
            panic!()
        };
        assert!(ss[0].points[1][1].is_nan());
        assert!(ss[0].band.as_ref().unwrap()[1][0].is_nan());
        assert_eq!(ss[0].points[2][0], 2.);
        let first = source.lines().find(|l| !l.starts_with('#')).unwrap();
        for bad in [
            format!("{source}{first}\n"),
            source.replacen("phase_00", "wrong_label", 1),
            source.replacen("12 Occipital 210", "12 Other 210", 1),
            source.replacen("0.039376", "9.0", 1),
            source.replace("bootHi", "unknown"),
        ] {
            assert!(read(&bad).is_err());
        }
    }
    #[test]
    fn ci_method_ambiguity_and_partial_triplets_require_explicit_choices() {
        let source = "# bootstrap: 100 percentile CI: 90%\n# cond_bootstrap: 100 percentile CI: 95%\n#ROI label nvox time_index time_label effect stat p q z bootLo bootHi cbootLo cbootHi\n1 A 2 0 zero 0.2 5 1 1 0 0.1 0.3 0.05 0.4\n1 A 2 1 one 0.3 6 1 1 0 NA 0.4 0.2 0.5\n";
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("ci.1D");
        std::fs::write(&p, source).unwrap();
        let file = p.to_str().unwrap();
        assert!(rsa::figure(&options(file, &["-rsa-plots", "intervals"])).is_err());
        assert!(rsa::figure(&options(file, &["-rsa-ci", "boot"])).is_err());
        let f = rsa::figure(&options(
            file,
            &[
                "-rsa-plots",
                "intervals",
                "-rsa-ci",
                "boot",
                "-missing",
                "gap",
            ],
        ))
        .unwrap();
        let Content::Intervals { rows, method, .. } = &f.plots[0].content else {
            panic!()
        };
        assert_eq!(rows[0].values, Some([0.2, 0.1, 0.3]));
        assert!(rows[1].values.is_none());
        assert!(method.contains("90%"));
        let f = rsa::figure(&options(file, &["-rsa-ci", "boot", "-missing", "gap"])).unwrap();
        let Content::Lines(ss) = &f.plots[0].content else {
            panic!()
        };
        assert_eq!(ss[0].points[1][1], 0.3);
        assert!(ss[0].band.as_ref().unwrap()[1][0].is_nan());
    }
    #[test]
    fn long_projects_freeze_selections_relocate_and_replay_current_view_pixels() {
        use crate::document::{self, Origin, Presentation, Snapshot};
        for name in ["series", "runs"] {
            let d = tempfile::tempdir().unwrap();
            let old = d.path().join("old");
            std::fs::create_dir(&old).unwrap();
            let p = old.join("source.1D");
            std::fs::copy(fixture(name), &p).unwrap();
            let o = options(p.to_str().unwrap(), &["-rsa-rois", "3"]);
            let f = rsa::figure(&o).unwrap();
            let origin = Origin::capture(&o).unwrap();
            let state = Snapshot {
                presentation: Presentation::capture(&f),
                views: vec![f.plots[0].bounds.pan(0.2, 0.)],
                width: 1000,
                height: 650,
                full_plot: false,
            };
            document::save(&old.join("figure.json"), &origin, &state, false).unwrap();
            let moved = d.path().join("moved");
            std::fs::rename(old, &moved).unwrap();
            let loaded = document::load(&moved.join("figure.json")).unwrap();
            let r = crate::render::Renderer::new();
            let pixels = |f: &Figure| {
                r.raster(&r.scene(f, &state.views, 1000, 650).unwrap())
                    .unwrap()
            };
            assert_eq!(pixels(&f).data(), pixels(&loaded.figure).data());
            std::fs::write(moved.join("source.1D"), "changed").unwrap();
            assert!(document::load(&moved.join("figure.json")).is_err());
        }
    }
}
