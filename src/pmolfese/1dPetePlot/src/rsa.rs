//! Versioned, deterministic presets for ordinary ROI-wide 3dRSA results.
//! No inference, significance filtering, RDM reconstruction, or cross-file row-position joins.
use crate::{
    cli::Options,
    input::Missing,
    model::{self, Bounds, Content, Figure, Plot, Series},
};
use anyhow::{Context, Result, bail, ensure};
use std::{
    collections::{BTreeMap, BTreeSet},
    path::Path,
};

#[derive(Clone, Debug)]
pub struct Request {
    pub file: Option<String>,
    pub models: Option<Vec<String>>,
    pub rois: Option<Vec<i32>>,
    pub plots: Option<Vec<String>>,
    pub effect: String,
    pub pair: Option<Vec<String>>,
    pub version: u32,
}
impl Default for Request {
    fn default() -> Self {
        Self {
            file: None,
            models: None,
            rois: None,
            plots: None,
            effect: "auto".into(),
            pair: None,
            version: 1,
        }
    }
}
pub fn names(s: &str) -> Result<Vec<String>> {
    let values: Vec<String> = s.split(',').map(|s| s.trim().to_owned()).collect();
    ensure!(
        !values.is_empty()
            && values.iter().all(|s| !s.is_empty())
            && values.iter().collect::<BTreeSet<_>>().len() == values.len(),
        "use a comma-separated list of unique, nonempty values"
    );
    Ok(values)
}
#[derive(Clone, Debug)]
pub struct Roi {
    pub id: i32,
    pub label: String,
    pub nvox: u32,
    pub values: Vec<f64>,
}
#[derive(Clone, Debug)]
pub struct Table {
    pub columns: Vec<String>,
    pub models: Vec<String>,
    pub effects: BTreeMap<(String, String), usize>,
    pub rows: BTreeMap<i32, Roi>,
    pub metadata: Vec<String>,
}
// Extra scientific columns are recognized but never selected as primary effects.
fn auxiliary(name: &str) -> bool {
    ["reliability", "nc_low", "nc_high"].contains(&name)
        || [
            "_p",
            "_q",
            "_pfwe",
            "_bootLo",
            "_bootHi",
            "_cbootLo",
            "_cbootHi",
            "_dualLo",
            "_dualHi",
            "_cvR",
            "_cvP",
            "_cvQ",
            "_cvPfwe",
            "_cvDiff",
            "_looR",
            "_looP",
            "_looQ",
            "_looPfwe",
            "_looBootLo",
            "_looBootHi",
            "_rDiff",
            "_zDiff",
            "_diff",
        ]
        .iter()
        .any(|suffix| name.ends_with(suffix))
}
pub fn read(path: &Path) -> Result<Table> {
    let stat =
        std::fs::metadata(path).with_context(|| format!("opening RSA table {}", path.display()))?;
    ensure!(
        stat.is_file() && stat.len() <= 16 * 1024 * 1024,
        "RSA input must be a regular file up to 16 MiB"
    );
    use std::io::Read;
    let mut text = String::new();
    std::fs::File::open(path)?
        .take(16 * 1024 * 1024 + 1)
        .read_to_string(&mut text)?;
    ensure!(
        text.len() <= 16 * 1024 * 1024,
        "RSA input grew beyond 16 MiB"
    );
    let mut table = Table {
        columns: vec![],
        models: vec![],
        effects: BTreeMap::new(),
        rows: BTreeMap::new(),
        metadata: vec![],
    };
    for (line_no, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(comment) = line.strip_prefix('#') {
            let fields: Vec<_> = comment.split_whitespace().collect();
            if fields.first() == Some(&"ROI") {
                ensure!(table.columns.is_empty(), "multiple RSA column headers");
                ensure!(
                    fields.len() >= 6 && fields[..3] == ["ROI", "label", "nvox"],
                    "expected ordinary ROI-wide RSA header: #ROI label nvox MODEL_r ..."
                );
                ensure!(
                    !fields
                        .iter()
                        .any(|s| ["summary", "time_index", "time_label", "model"].contains(s)),
                    "long series/run-resolved RSA tables are not supported yet; use ordinary ROI-wide results"
                );
                ensure!(
                    fields.len() <= 512 && fields.iter().all(|s| s.len() <= 128),
                    "RSA header exceeds 512 columns / 128 bytes per name"
                );
                ensure!(
                    fields.iter().collect::<BTreeSet<_>>().len() == fields.len(),
                    "duplicate RSA column names"
                );
                table.columns = fields.iter().map(|s| s.to_string()).collect();
            } else if table.metadata.len() < 100 {
                table.metadata.push(comment.trim().to_owned());
            }
            continue;
        }
        ensure!(
            !table.columns.is_empty(),
            "RSA data requires a #ROI label nvox ... header; numeric RDMs use -heat instead"
        );
        let fields: Vec<_> = line.split_whitespace().collect();
        let n = table.columns.len();
        ensure!(
            fields.len() >= n,
            "RSA row {} has too few fields: expected {n}",
            line_no + 1
        );
        // The label is the only text column, and AFNI can emit labels containing spaces.
        let numeric_start = fields.len() - (n - 2);
        let id: i32 = fields[0]
            .parse()
            .with_context(|| format!("invalid ROI ID on row {}", line_no + 1))?;
        ensure!(id > 0, "ROI IDs must be positive");
        let label = fields[1..numeric_start].join(" ");
        ensure!(
            !label.is_empty() && label.len() <= 256,
            "RSA label must be 1..256 bytes"
        );
        let nvox = fields[numeric_start]
            .parse::<u32>()
            .context("nvox must be a nonnegative integer")?;
        let mut values = vec![f64::NAN; 3];
        for (column, field) in table.columns[3..].iter().zip(&fields[numeric_start + 1..]) {
            let value = if ["NA", "N/A", "NAN"].contains(&field.to_ascii_uppercase().as_str()) {
                f64::NAN
            } else {
                field
                    .parse::<f64>()
                    .with_context(|| format!("ROI {id}: invalid numeric {column}: {field}"))?
            };
            ensure!(!value.is_infinite(), "ROI {id}: {column} is infinite");
            if [
                "_p", "_q", "_pfwe", "_cvP", "_cvQ", "_cvPfwe", "_looP", "_looQ", "_looPfwe",
            ]
            .iter()
            .any(|s| column.ends_with(s))
            {
                ensure!(
                    value.is_nan() || (0.0..=1.0).contains(&value),
                    "ROI {id}: {column} must be a probability in 0..1 or missing"
                );
            }
            values.push(value);
        }
        ensure!(
            table
                .rows
                .insert(
                    id,
                    Roi {
                        id,
                        label,
                        nvox,
                        values
                    }
                )
                .is_none(),
            "duplicate ROI ID {id}; long/duplicated tables cannot be paired by row position"
        );
        ensure!(
            table.rows.len() <= 10_000 && table.rows.len() * n <= 1_000_000,
            "RSA input exceeds 10000 ROIs / 1000000 cells; subset the table explicitly"
        );
    }
    ensure!(!table.rows.is_empty(), "RSA input contains no ROI rows");
    let weight_prefixes: Vec<_> = table
        .columns
        .iter()
        .filter_map(|c| c.strip_suffix("_cvR").map(|fit| format!("{fit}_w_")))
        .collect();
    for (i, name) in table.columns.iter().enumerate().skip(3) {
        if weight_prefixes
            .iter()
            .any(|prefix| name.starts_with(prefix))
        {
            continue;
        }
        let effect = ["r", "b", "pr"].iter().find_map(|kind| {
            name.strip_suffix(&format!("_{kind}"))
                .map(|model| (model, *kind))
        });
        if let Some((model, kind)) = effect {
            ensure!(
                !model.is_empty() && !model.contains(',') && !model.starts_with('-'),
                "RSA model names must be nonempty and cannot contain commas"
            );
            if !table.models.iter().any(|s| s == model) {
                table.models.push(model.into());
            }
            table.effects.insert((model.into(), kind.into()), i);
        } else {
            ensure!(
                auxiliary(name),
                "unsupported RSA column '{name}'; this preset supports ordinary model r/b/pr effects (not commonality or unknown schemas)"
            );
        }
    }
    ensure!(
        !table.models.is_empty(),
        "RSA table has no supported MODEL_r or MODEL_b/MODEL_pr effects"
    );
    for model in &table.models {
        for suffix in ["p", "q"] {
            ensure!(
                table.columns.contains(&format!("{model}_{suffix}")),
                "model {model} is missing its _{suffix} column; unsupported RSA schema"
            );
        }
        let has = |kind: &str| table.effects.contains_key(&(model.clone(), kind.into()));
        ensure!(
            (has("r") && !has("b") && !has("pr")) || (!has("r") && has("b") && has("pr")),
            "model {model}: expected r or paired b/pr effects, not mixed/incomplete roles"
        );
    }
    for (lo, hi) in [
        ("_bootLo", "_bootHi"),
        ("_cbootLo", "_cbootHi"),
        ("_dualLo", "_dualHi"),
        ("_looBootLo", "_looBootHi"),
    ] {
        for (i, column) in table.columns.iter().enumerate() {
            if let Some(base) = column.strip_suffix(lo) {
                let j = table
                    .columns
                    .iter()
                    .position(|s| s == &format!("{base}{hi}"))
                    .context("RSA confidence interval is missing its upper bound")?;
                for row in table.rows.values() {
                    ensure!(
                        row.values[i].is_nan()
                            || row.values[j].is_nan()
                            || row.values[i] <= row.values[j],
                        "ROI {}: reversed interval bounds for {base}",
                        row.id
                    );
                }
            } else if let Some(base) = column.strip_suffix(hi) {
                ensure!(
                    table.columns.contains(&format!("{base}{lo}")),
                    "RSA confidence interval is missing its lower bound"
                );
            }
        }
    }
    Ok(table)
}
#[derive(Clone, Debug)]
struct Selection {
    models: Vec<String>,
    rois: Vec<i32>,
    plots: Vec<String>,
    effect: String,
    pair: Vec<String>,
    heat_omitted: bool,
}
fn resolve(table: &Table, request: &Request) -> Result<Selection> {
    ensure!(request.version == 1, "unsupported RSA preset version");
    let models = request
        .models
        .clone()
        .unwrap_or_else(|| table.models.clone());
    ensure!(
        !models.is_empty() && models.len() <= 8,
        "choose at most 8 models with -rsa-models NAME,NAME (available: {})",
        table.models.join(",")
    );
    for m in &models {
        ensure!(
            table.models.contains(m),
            "unknown RSA model '{m}'; available: {}",
            table.models.join(",")
        );
    }
    let rois = request
        .rois
        .clone()
        .unwrap_or_else(|| table.rows.keys().copied().collect());
    ensure!(
        !rois.is_empty() && rois.iter().collect::<BTreeSet<_>>().len() == rois.len(),
        "select unique ROI IDs"
    );
    for id in &rois {
        ensure!(
            table.rows.contains_key(id),
            "requested ROI {id} is absent from the table"
        );
    }
    let effect = if request.effect == "auto" {
        if models.iter().all(|m| table.effects.contains_key(&(m.clone(), "r".into()))) { "r" }
        else if models.iter().all(|m| table.effects.contains_key(&(m.clone(), "b".into()))) { "b" }
        else { bail!("selected models have different effect kinds; choose compatible models and -rsa-effect r|b|pr"); }
    } else { request.effect.as_str() }.to_owned();
    ensure!(
        ["r", "b", "pr"].contains(&effect.as_str()),
        "RSA effect must be auto, r, b, or pr"
    );
    for m in &models {
        ensure!(
            table.effects.contains_key(&(m.clone(), effect.clone())),
            "model {m} has no {effect} effect"
        );
    }
    let heat_omitted = request.plots.is_none() && rois.len() > 40;
    let plots = request.plots.clone().unwrap_or_else(|| {
        let mut p = vec!["histogram".into()];
        if models.len() > 1 {
            p.push("scatter".into());
        }
        if rois.len() <= 40 {
            p.push("heatmap".into());
        }
        p
    });
    ensure!(
        !plots.is_empty()
            && plots
                .iter()
                .all(|s| ["histogram", "scatter", "heatmap"].contains(&s.as_str())),
        "RSA plots must be histogram,scatter,heatmap"
    );
    ensure!(
        !plots.iter().any(|s| s == "heatmap") || rois.len() <= 40,
        "effect heatmaps support up to 40 ROIs; select -rsa-rois IDs or -rsa-plots histogram,scatter"
    );
    let pair = request
        .pair
        .clone()
        .unwrap_or_else(|| models.iter().take(2).cloned().collect());
    if plots.iter().any(|s| s == "scatter") {
        ensure!(
            pair.len() == 2 && pair[0] != pair[1] && pair.iter().all(|m| models.contains(m)),
            "scatter needs two selected models; use -rsa-pair XMODEL,YMODEL"
        );
    } else {
        ensure!(
            request.pair.is_none(),
            "-rsa-pair requires scatter in -rsa-plots"
        );
    }
    Ok(Selection {
        models,
        rois,
        plots,
        effect,
        pair,
        heat_omitted,
    })
}
/// Persist resolved names, ROI IDs, modes and version rather than rerunning future defaults.
pub fn freeze(o: &Options) -> Result<Vec<String>> {
    let Some(file) = &o.rsa.file else {
        return Ok(o.arguments.clone());
    };
    let table = read(Path::new(file))?;
    let selected = resolve(&table, &o.rsa)?;
    let mut args = o.arguments.clone();
    for (flag, value) in [
        ("-rsa-version", "1".into()),
        ("-rsa-models", selected.models.join(",")),
        (
            "-rsa-rois",
            selected
                .rois
                .iter()
                .map(i32::to_string)
                .collect::<Vec<_>>()
                .join(","),
        ),
        ("-rsa-effect", selected.effect),
        ("-rsa-plots", selected.plots.join(",")),
    ] {
        if !args.iter().any(|a| a == flag) {
            args.extend([flag.into(), value]);
        }
    }
    // Explicit auto must also freeze to the resolved effect.
    if let Some(i) = args.iter().position(|a| a == "-rsa-effect") {
        args[i + 1] = resolve(&table, &o.rsa)?.effect;
    }
    if selected.plots.iter().any(|s| s == "scatter") && o.rsa.pair.is_none() {
        args.extend(["-rsa-pair".into(), selected.pair.join(",")]);
    }
    Ok(args)
}
fn plot(content: Content, title: String, xlabel: String, ylabel: String, bounds: Bounds) -> Plot {
    Plot {
        id: 0,
        caption: String::new(),
        title,
        note: String::new(),
        xlabel,
        ylabel,
        content,
        bounds,
        diagonal: false,
        refs: vec![],
        legend: false,
        qc: Default::default(),
    }
}
pub fn figure(o: &Options) -> Result<Figure> {
    let path = Path::new(o.rsa.file.as_ref().context("missing -rsa FILE")?);
    let table = read(path)?;
    let s = resolve(&table, &o.rsa)?;
    ensure!(
        s.plots.iter().any(|p| p == "histogram")
            || (o.bins.is_none() && o.binwidth.is_none() && !o.density),
        "histogram settings require histogram in -rsa-plots"
    );
    let kind = match s.effect.as_str() {
        "b" => "beta (b)",
        "pr" => "partial r",
        _ => "RSA comparison (r)",
    };
    let column = |m: &str| table.effects[&(m.into(), s.effect.clone())];
    let selected: Vec<_> = s.rois.iter().map(|id| &table.rows[id]).collect();
    if o.missing == Missing::Error {
        for row in &selected {
            for m in &s.models {
                ensure!(
                    row.values[column(m)].is_finite(),
                    "ROI {}: {m}_{} is missing; use -missing gap to omit missing values explicitly",
                    row.id,
                    s.effect
                );
            }
        }
    }
    ensure!(
        o.missing != Missing::Mean,
        "RSA presets do not impute effects"
    );
    let color = |i: usize| {
        model::color(
            o.colors
                .get(i % o.colors.len().max(1))
                .map(String::as_str)
                .unwrap_or(model::COLORS[i % model::COLORS.len()]),
        )
    };
    let mut plots = vec![];
    for mode in &s.plots {
        match mode.as_str() {
            "histogram" => {
                for (i, m) in s.models.iter().enumerate() {
                    let values: Vec<_> = selected
                        .iter()
                        .map(|r| r.values[column(m)])
                        .filter(|v| v.is_finite())
                        .collect();
                    ensure!(
                        !values.is_empty(),
                        "model {m} has no finite effects in the selected ROIs"
                    );
                    let (edges, points) = model::histogram(&values, o.bins, o.binwidth, o.density)?;
                    let max = points.iter().map(|p| p[1]).fold(0., f64::max);
                    let bounds = Bounds {
                        x: model::padded(edges[0], *edges.last().unwrap()),
                        y: [0., max * 1.12],
                    };
                    let mut p = plot(
                        Content::Histogram {
                            series: vec![Series {
                                label: m.clone(),
                                color: color(i)?,
                                points,
                                band: None,
                                point_labels: vec![],
                            }],
                            edges: vec![edges],
                            density: o.density,
                        },
                        format!("{m} · {kind}"),
                        kind.into(),
                        if o.density { "Density" } else { "ROI count" }.into(),
                        bounds,
                    );
                    p.refs = vec![
                        (0., None),
                        (
                            values.iter().sum::<f64>() / values.len() as f64,
                            Some("Mean".into()),
                        ),
                    ];
                    p.note = format!(
                        "{} finite / {} selected ROIs · no significance filter",
                        values.len(),
                        selected.len()
                    );
                    plots.push(p);
                }
            }
            "scatter" => {
                let x = &s.pair[0];
                let y = &s.pair[1];
                // Both lookups use the same unique ROI identity; gaps omit the pair together.
                let pairs: Vec<_> = selected
                    .iter()
                    .filter_map(|r| {
                        let p = [r.values[column(x)], r.values[column(y)]];
                        p.iter().all(|v| v.is_finite()).then_some((*r, p))
                    })
                    .collect();
                ensure!(
                    !pairs.is_empty(),
                    "no paired finite ROI effects for {x} vs {y}"
                );
                let lo = pairs
                    .iter()
                    .flat_map(|(_, p)| p.iter())
                    .copied()
                    .fold(f64::INFINITY, f64::min);
                let hi = pairs
                    .iter()
                    .flat_map(|(_, p)| p.iter())
                    .copied()
                    .fold(f64::NEG_INFINITY, f64::max);
                let range = model::padded(lo, hi);
                let mut p = plot(
                    Content::Scatter(vec![Series {
                        label: format!("{x} / {y}"),
                        color: color(0)?,
                        points: pairs.iter().map(|(_, p)| *p).collect(),
                        band: None,
                        point_labels: pairs
                            .iter()
                            .map(|(r, _)| format!("ROI {} · {} · nvox {}", r.id, r.label, r.nvox))
                            .collect(),
                    }]),
                    format!("{x} vs {y}"),
                    format!("{x} · {kind}"),
                    format!("{y} · {kind}"),
                    Bounds { x: range, y: range },
                );
                p.diagonal = true;
                p.note = format!(
                    "{} paired / {} selected ROIs · identity line",
                    pairs.len(),
                    selected.len()
                );
                plots.push(p);
            }
            "heatmap" => {
                let matrix: Vec<Vec<f64>> = selected
                    .iter()
                    .map(|r| s.models.iter().map(|m| r.values[column(m)]).collect())
                    .collect();
                let finite: Vec<_> = matrix
                    .iter()
                    .flatten()
                    .copied()
                    .filter(|v| v.is_finite())
                    .collect();
                ensure!(!finite.is_empty(), "no finite effects for the heatmap");
                let max = finite.iter().map(|v| v.abs()).fold(0., f64::max).max(1e-12);
                let limit = if s.effect == "b" { max } else { max.max(1.) };
                let mut p = plot(
                    Content::Heat {
                        matrix,
                        labels: s.models.clone(),
                        row_labels: selected
                            .iter()
                            .map(|r| format!("{}: {}", r.id, r.label))
                            .collect(),
                        square_cells: false,
                        cmap: "Reds_and_Blues_Inv".into(),
                        clim: [-limit, limit],
                    },
                    format!("ROI × model · {kind} effects"),
                    "Model".into(),
                    "ROI".into(),
                    Bounds {
                        x: [0., s.models.len() as f64],
                        y: [0., selected.len() as f64],
                    },
                );
                p.note = format!(
                    "{} ROIs × {} models · gray = missing · not an RDM",
                    selected.len(),
                    s.models.len()
                );
                plots.push(p);
            }
            _ => unreachable!(),
        }
    }
    let mut figure = model::presentation(
        Figure {
            title: o.title.clone().unwrap_or_else(|| {
                format!(
                    "RSA · {}",
                    path.file_name().unwrap_or_default().to_string_lossy()
                )
            }),
            caption: String::new(),
            subtitle: format!(
                "RSA preset v1 · {} · {} ROIs · {} models",
                kind,
                selected.len(),
                s.models.len()
            ),
            columns: 2.min(plots.len()),
            plots,
            synthetic: false,
        },
        o,
    )?;
    let mut provenance = "Supplied effects only; p/q/FWE, intervals, fits, contrasts and reliability are not plotted or recomputed.".to_string();
    if s.heat_omitted {
        provenance.push_str(" Heatmap omitted above 40 ROIs; select -rsa-rois to include it.");
    }
    if let Some(metric) = table.metadata.iter().find(|m| m.starts_with("metric")) {
        provenance.push_str(&format!(" Source {metric}."));
    }
    if figure.caption.is_empty() {
        figure.caption = provenance;
    }
    for p in &figure.plots {
        ensure!(p.bounds.valid(), "RSA plot bounds overflow");
    }
    Ok(figure)
}

/// Leave enough vertical room for every ROI label in an effect heatmap.
pub fn minimum_height(fig: &Figure, width: u32) -> u32 {
    let panel = fig
        .plots
        .iter()
        .filter_map(|p| match &p.content {
            Content::Heat {
                matrix,
                square_cells: false,
                ..
            } => Some(matrix.len() as f64 * 13. + 260.),
            _ => None,
        })
        .fold(0., f64::max);
    if panel == 0. {
        0
    } else {
        ((200. + fig.plots.len().div_ceil(fig.columns) as f64 * panel) * width as f64 / 1200.)
            .ceil() as u32
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn fixture(name: &str) -> String {
        format!("{}/examples/rsa/{name}.rsa.1D", env!("CARGO_MANIFEST_DIR"))
    }
    fn options(file: &str, extra: &[&str]) -> Options {
        let mut args = vec!["-rsa".to_owned(), file.to_owned()];
        args.extend(extra.iter().map(|s| s.to_string()));
        crate::cli::parse(&args).unwrap()
    }
    #[test]
    fn ordinary_roles_pair_by_roi_and_rectangular_heatmap_matches_source() {
        let file = fixture("ordinary");
        let table = read(Path::new(&file)).unwrap();
        assert_eq!(table.models, ["semantic", "visual"]);
        assert_eq!(table.rows[&3].label, "Left temporal");
        assert_eq!(table.effects.len(), 2); // Fit, p/q, intervals and nvox never become effects.
        let f = figure(&options(&file, &[])).unwrap();
        assert_eq!(f.plots.len(), 4);
        let Content::Histogram { series, .. } = &f.plots[0].content else {
            panic!()
        };
        assert_eq!(series[0].points.iter().map(|p| p[1]).sum::<f64>(), 8.);
        let Content::Scatter(series) = &f.plots[2].content else {
            panic!()
        };
        assert_eq!(series[0].points[0], [0.68, 0.21]); // ROI 3, not the first file row.
        assert!(
            f.plots[2]
                .hover(0.68, 0.21)
                .contains("ROI 3 · Left temporal · nvox 185")
        );
        let Content::Heat {
            matrix,
            labels,
            row_labels,
            square_cells,
            ..
        } = &f.plots[3].content
        else {
            panic!()
        };
        assert_eq!(matrix.len(), 8);
        assert_eq!(matrix[0], [0.68, 0.21]);
        assert_eq!(labels, &["semantic", "visual"]);
        assert_eq!(row_labels[0], "3: Left temporal");
        assert!(!square_cells);
        assert!(f.plots[3].hover(1.5, 7.5).contains("0.21000"));
        assert_eq!(f.plots[3].hover(2.0, 7.5), "");
    }
    #[test]
    fn regression_and_selections_keep_effect_roles_and_requested_order() {
        let file = fixture("regression");
        for (effect, expected) in [("b", 1.4), ("pr", 0.6)] {
            let f = figure(&options(
                &file,
                &[
                    "-rsa-rois",
                    "9,2",
                    "-rsa-models",
                    "B,A",
                    "-rsa-effect",
                    effect,
                    "-rsa-plots",
                    "heatmap",
                ],
            ))
            .unwrap();
            let Content::Heat { matrix, labels, .. } = &f.plots[0].content else {
                panic!()
            };
            assert_eq!(labels, &["B", "A"]);
            assert_eq!(matrix[0][1], expected);
        }
        let f = figure(&options(&file, &[])).unwrap();
        assert!(f.subtitle.contains("beta"));
        for args in [
            vec!["-rsa-models", "missing"],
            vec!["-rsa-rois", "99"],
            vec!["-rsa-effect", "r"],
            vec!["-rsa-pair", "A,A"],
            vec!["-rsa-plots", "unknown"],
        ] {
            let mut command = vec!["-rsa".into(), file.clone()];
            command.extend(args.iter().map(|s| s.to_string()));
            assert!(
                crate::cli::parse(&command)
                    .and_then(|o| figure(&o))
                    .is_err()
            );
        }
    }
    #[test]
    fn missing_values_omit_complete_pairs_and_do_not_impute_or_filter_by_significance() {
        let d = tempfile::tempdir().unwrap();
        let path = d.path().join("missing.rsa.1D");
        std::fs::write(&path, "#ROI label nvox A_r A_p A_q B_r B_p B_q\n2 two 10 0.8 1 1 NaN 1 1\n1 one 20 0.2 1 1 0.3 1 1\n3 three 30 NaN NaN NaN 0.5 0.01 0.01\n").unwrap();
        let file = path.to_str().unwrap();
        assert!(
            figure(&options(file, &[]))
                .unwrap_err()
                .to_string()
                .contains("missing")
        );
        let f = figure(&options(file, &["-missing", "gap"])).unwrap();
        let Content::Scatter(series) = &f.plots[2].content else {
            panic!()
        };
        assert_eq!(series[0].points, [[0.2, 0.3]]);
        assert!(series[0].point_labels[0].contains("ROI 1"));
        let Content::Heat { matrix, .. } = &f.plots[3].content else {
            panic!()
        };
        assert!(matrix[1][1].is_nan());
        assert!(matrix[2][0].is_nan());
        assert!(
            crate::cli::parse(&["-rsa".into(), file.into(), "-missing".into(), "mean".into()])
                .is_err()
        );
    }
    #[test]
    fn unsupported_schemas_and_malformed_scientific_metadata_fail() {
        let d = tempfile::tempdir().unwrap();
        let path = d.path().join("bad.1D");
        for text in [
            "1 2 3\n",
            "#ROI label nvox summary model effect p q\n1 A 2 MEAN M 0.1 0.2 0.3\n",
            "#ROI label nvox time_index time_label effect stat p q\n1 A 2 0 t0 0.1 0.2 0.3 0.4\n",
            "#ROI label nvox M_r M_p M_q\n1 A 2 .1 .2 .3\n1 B 3 .2 .3 .4\n",
            "#ROI label nvox M_r M_p M_q mystery\n1 A 2 .1 .2 .3 .4\n",
            "#ROI label nvox M_r M_p M_q\n1 A 2 .1 2 .3\n",
            "#ROI label nvox M_r M_p M_q M_bootLo M_bootHi\n1 A 2 .1 .2 .3 .8 .2\n",
            "#ROI label nvox M_b M_p M_q\n1 A 2 .1 .2 .3\n",
            "#ROI label nvox M_r M_p M_q\n1 A 2 inf .2 .3\n",
            "#ROI label nvox M_r M_p M_q M_bootLo\n1 A 2 .1 .2 .3 .8\n",
        ] {
            std::fs::write(&path, text).unwrap();
            assert!(read(&path).is_err(), "accepted {text}");
        }
    }
    #[test]
    fn bounded_defaults_never_truncate_distribution_observations() {
        let d = tempfile::tempdir().unwrap();
        let path = d.path().join("large.1D");
        let mut text = "#ROI label nvox M_r M_p M_q\n".to_string();
        for i in 1..=50 {
            text.push_str(&format!("{i} roi_{i} 10 {} 1 1\n", i as f64 / 100.));
        }
        std::fs::write(&path, text).unwrap();
        let o = options(path.to_str().unwrap(), &[]);
        let f = figure(&o).unwrap();
        assert_eq!(f.plots.len(), 1);
        let Content::Histogram { series, .. } = &f.plots[0].content else {
            panic!()
        };
        assert_eq!(series[0].points.iter().map(|p| p[1]).sum::<f64>(), 50.);
        assert!(f.caption.contains("Heatmap omitted"));
        assert!(figure(&options(path.to_str().unwrap(), &["-rsa-plots", "heatmap"])).is_err());
    }
    #[test]
    fn project_freezes_preset_and_replays_relocated_sources_with_identical_pixels() {
        let d = tempfile::tempdir().unwrap();
        let path = d.path().join("results.rsa.1D");
        std::fs::copy(fixture("ordinary"), &path).unwrap();
        let o = options(path.to_str().unwrap(), &["-rsa-effect", "auto"]);
        let f = figure(&o).unwrap();
        let origin = crate::document::Origin::capture(&o).unwrap();
        let state = crate::document::Snapshot {
            presentation: crate::document::Presentation::capture(&f),
            views: f.plots.iter().map(|p| p.bounds).collect(),
            width: 1500,
            height: 1300,
            full_plot: true,
        };
        let recipe = d.path().join("figure.json");
        crate::document::save(&recipe, &origin, &state, false).unwrap();
        let json = std::fs::read_to_string(&recipe).unwrap();
        assert!(json.contains("-rsa-version"));
        assert!(!json.contains("\"auto\""));
        let moved = d.path().join("moved");
        std::fs::create_dir(&moved).unwrap();
        std::fs::rename(&path, moved.join("results.rsa.1D")).unwrap();
        std::fs::rename(&recipe, moved.join("figure.json")).unwrap();
        let loaded = crate::document::load(&moved.join("figure.json")).unwrap();
        let renderer = crate::render::Renderer::new();
        let pixels = |fig: &Figure| {
            renderer
                .raster(
                    &renderer
                        .scene(fig, &state.views, state.width, state.height)
                        .unwrap(),
                )
                .unwrap()
        };
        assert_eq!(pixels(&f).data(), pixels(&loaded.figure).data());
        std::fs::write(moved.join("results.rsa.1D"), "changed").unwrap();
        assert!(crate::document::load(&moved.join("figure.json")).is_err());
    }
}
