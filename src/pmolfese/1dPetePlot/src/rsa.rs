//! Versioned, deterministic presets for 3dRSA results; long schemas dispatch to rsa_long.
//! No inference, significance filtering, RDM reconstruction, or cross-file row-position joins.
use crate::{
    cli::Options,
    input::Missing,
    model::{self, Bounds, Content, Figure, IntervalRow, Plot, Series},
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
    pub ci: String,
    pub times: Option<Vec<String>>,
    pub summaries: Option<Vec<String>>,
    pub family: String,
    pub measures: Option<Vec<String>>,
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
            ci: "auto".into(),
            times: None,
            summaries: None,
            family: "effects".into(),
            measures: None,
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
    pub long: Option<Box<crate::rsa_long::LongTable>>,
    pub columns: Vec<String>,
    pub models: Vec<String>,
    pub effects: BTreeMap<(String, String), usize>,
    /// ROI identities; long-table observations live in `long.rows`.
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
    if crate::rsa_long::detect(&text) {
        return crate::rsa_long::read(&text);
    }
    let mut table = Table {
        long: None,
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
                    "unsupported mixed RSA identity schema; expected a complete wide or long header"
                );
                ensure!(
                    fields.len() <= 512 && fields.iter().all(|s| s.len() <= 256),
                    "RSA header exceeds 512 columns / 256 bytes per name"
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
        if crate::rsa_extra::component(&table, name) {
            continue;
        }
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
                auxiliary(name) || crate::rsa_extra::component(&table, name),
                "unsupported RSA column '{name}'; unknown scientific schema"
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
    crate::rsa_extra::validate(&table)?;
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
    ci: Option<String>,
}
fn interval_columns(table: &Table, model: &str, method: &str) -> Option<[usize; 2]> {
    Some([
        table
            .columns
            .iter()
            .position(|c| c == &format!("{model}_{method}Lo"))?,
        table
            .columns
            .iter()
            .position(|c| c == &format!("{model}_{method}Hi"))?,
    ])
}
pub(crate) fn interval_label(table: &Table, method: &str) -> Result<String> {
    let (prefix, label) = match method {
        "boot" => ("bootstrap:", "Subject bootstrap"),
        "cboot" => ("cond_bootstrap:", "Condition bootstrap"),
        "dual" => ("dual_bootstrap:", "Dual bootstrap"),
        _ => unreachable!(),
    };
    let mut level = None;
    let mut construction = None;
    for meta in table.metadata.iter().filter(|m| m.starts_with(prefix)) {
        if let Some((_, tail)) = meta.split_once("CI:") {
            let recorded = if meta.contains("corrected variance t CI:") {
                Some("corrected-variance t")
            } else if meta.contains("percentile CI:") {
                Some("percentile")
            } else {
                None
            };
            if let Some(recorded) = recorded {
                ensure!(
                    construction.is_none_or(|old| old == recorded),
                    "conflicting RSA CI construction metadata"
                );
                construction = Some(recorded);
            }
            let token = tail
                .trim()
                .split_once('%')
                .context("invalid RSA confidence-level metadata")?
                .0
                .trim();
            let value: f64 = token
                .parse()
                .context("invalid RSA confidence-level metadata")?;
            ensure!(
                value.is_finite() && value > 0. && value < 100.,
                "confidence level must be between 0 and 100 percent"
            );
            ensure!(
                level.is_none_or(|old| old == value),
                "conflicting RSA confidence-level metadata"
            );
            level = Some(value);
        }
    }
    let mut label = match level {
        Some(level) => format!("{label} · {level}% CI"),
        None => format!("{label} · CI level not recorded"),
    };
    if let Some(construction) = construction {
        label.push_str(&format!(" ({construction})"));
    }
    Ok(label)
}
fn resolve(table: &Table, request: &Request) -> Result<Selection> {
    ensure!(request.version == 1, "unsupported RSA preset version");
    ensure!(
        request.times.is_none() && request.summaries.is_none(),
        "-rsa-times/-rsa-summaries require a long RSA table"
    );
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
                .all(|s| ["histogram", "scatter", "heatmap", "intervals"].contains(&s.as_str())),
        "RSA plots must be histogram,scatter,heatmap,intervals"
    );
    ensure!(
        ["auto", "boot", "cboot", "dual"].contains(&request.ci.as_str()),
        "-rsa-ci must be auto, boot, cboot, or dual"
    );
    let ci = if plots.iter().any(|s| s == "intervals") {
        ensure!(
            rois.len() <= 40,
            "interval plots support up to 40 ROIs; select -rsa-rois IDs"
        );
        ensure!(
            effect != "pr",
            "3dRSA regression CI bounds belong to beta (b), not partial r; use -rsa-effect b"
        );
        let available: Vec<_> = ["boot", "cboot", "dual"]
            .into_iter()
            .filter(|method| {
                models
                    .iter()
                    .all(|m| interval_columns(table, m, method).is_some())
            })
            .collect();
        let method = if request.ci == "auto" {
            ensure!(
                !available.is_empty(),
                "no common supplied CI bounds for selected models; select models with matching interval columns"
            );
            ensure!(
                available.len() == 1,
                "multiple CI methods available ({}); choose -rsa-ci boot|cboot|dual",
                available.join(", ")
            );
            available[0]
        } else {
            ensure!(
                available.contains(&request.ci.as_str()),
                "selected CI method {} has no bounds for every selected model",
                request.ci
            );
            &request.ci
        };
        Some(method.to_owned())
    } else {
        ensure!(
            request.ci == "auto",
            "-rsa-ci requires intervals in -rsa-plots"
        );
        None
    };
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
        ci,
    })
}
/// Persist resolved names, ROI IDs, modes and version rather than rerunning future defaults.
pub fn freeze(o: &Options) -> Result<Vec<String>> {
    let Some(file) = &o.rsa.file else {
        return Ok(o.arguments.clone());
    };
    let table = read(Path::new(file))?;
    if o.rsa.family != "effects" {
        return crate::rsa_extra::freeze(o, &table);
    }
    ensure!(
        o.rsa.measures.is_none(),
        "-rsa-measures requires a fit/commonality family"
    );
    if let Some(long) = &table.long {
        return crate::rsa_long::freeze(o, &table, long);
    }
    freeze_table(o, &table)
}
pub(crate) fn freeze_table(o: &Options, table: &Table) -> Result<Vec<String>> {
    let selected = resolve(table, &o.rsa)?;
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
        args[i + 1] = resolve(table, &o.rsa)?.effect;
    }
    if selected.plots.iter().any(|s| s == "scatter") && o.rsa.pair.is_none() {
        args.extend(["-rsa-pair".into(), selected.pair.join(",")]);
    }
    if let Some(method) = selected.ci {
        if let Some(i) = args.iter().position(|a| a == "-rsa-ci") {
            args[i + 1] = method;
        } else {
            args.extend(["-rsa-ci".into(), method]);
        }
    }
    Ok(args)
}
pub(crate) fn plot(
    content: Content,
    title: String,
    xlabel: String,
    ylabel: String,
    bounds: Bounds,
) -> Plot {
    Plot {
        style: Default::default(),
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
    if o.rsa.family != "effects" {
        return crate::rsa_extra::figure(o, &table);
    }
    ensure!(
        o.rsa.measures.is_none(),
        "-rsa-measures requires a fit/commonality family"
    );
    if let Some(long) = &table.long {
        return crate::rsa_long::figure(o, &table, long);
    }
    figure_table(o, &table, None)
}
pub(crate) fn figure_table(o: &Options, table: &Table, label: Option<&str>) -> Result<Figure> {
    let path = Path::new(o.rsa.file.as_ref().context("missing RSA file")?);
    let s = resolve(table, &o.rsa)?;
    ensure!(
        s.plots.iter().any(|p| p == "histogram")
            || (o.bins.is_none() && o.binwidth.is_none() && !o.density),
        "histogram settings require histogram in -rsa-plots"
    );
    let kind = label.unwrap_or(match s.effect.as_str() {
        "b" => "beta (b)",
        "pr" => "partial r",
        _ => "RSA comparison (r)",
    });
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
            "intervals" => {
                let method = s.ci.as_deref().unwrap();
                let method_label = interval_label(table, method)?;
                let mut interval_plots = vec![];
                let (mut lo, mut hi) = (0_f64, 0_f64);
                for (i, m) in s.models.iter().enumerate() {
                    let [lower, upper] = interval_columns(table, m, method).unwrap();
                    let mut rows = vec![];
                    for r in &selected {
                        let values = [r.values[column(m)], r.values[lower], r.values[upper]];
                        let complete = values.iter().all(|v| v.is_finite());
                        ensure!(
                            complete || o.missing == Missing::Gap,
                            "ROI {}: {m} effect/CI is incomplete; use -missing gap to retain its row marked missing",
                            r.id
                        );
                        if complete {
                            for v in values {
                                lo = lo.min(v);
                                hi = hi.max(v);
                            }
                        }
                        rows.push(IntervalRow {
                            label: format!("{}: {}", r.id, r.label),
                            detail: format!("ROI {} · {} · nvox {}", r.id, r.label, r.nvox),
                            values: complete.then_some(values),
                        });
                    }
                    let count = rows.iter().filter(|r| r.values.is_some()).count();
                    ensure!(
                        count > 0,
                        "model {m} has no complete effect/CI triplets in the selected ROIs"
                    );
                    let mut p = plot(
                        Content::Intervals {
                            rows,
                            color: color(i)?,
                            method: method_label.clone(),
                        },
                        format!("{m} · {kind}"),
                        kind.into(),
                        "ROI".into(),
                        Bounds {
                            x: [0., 1.],
                            y: [0., selected.len() as f64],
                        },
                    );
                    p.note = method_label.clone();
                    p.caption = format!(
                        "{count}/{} complete ROIs · supplied endpoints; no significance filter.",
                        selected.len()
                    );
                    p.refs.push((0., None));
                    interval_plots.push(p);
                }
                // Common effect scale across the models, including every endpoint and zero.
                let range = model::padded(lo, hi);
                for p in &mut interval_plots {
                    p.bounds.x = range;
                }
                plots.extend(interval_plots);
            }
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
                        named_conditions: false,
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
    let mut provenance = if s.ci.is_some() {
        "Supplied effects and CI endpoints; no recalculation or scale conversion. Missing triplets retain a labeled row without a point or interval. p/q/FWE, fits, contrasts and reliability are not plotted."
    } else {
        "Supplied effects only; p/q/FWE, intervals, fits, contrasts and reliability are not plotted or recomputed."
    }.to_string();
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

/// Leave enough vertical room for every ROI label in an effect heatmap or interval plot.
pub fn minimum_height(fig: &Figure, width: u32) -> u32 {
    let panel = fig
        .plots
        .iter()
        .filter_map(|p| match &p.content {
            Content::Intervals { rows, .. } => Some(rows.len() as f64 * 17. + 290.),
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
    fn intervals_preserve_asymmetric_endpoints_roi_order_and_effect_scale() {
        let f = figure(&options(
            &fixture("ordinary"),
            &["-rsa-plots", "intervals", "-rsa-rois", "12,3"],
        ))
        .unwrap();
        assert_eq!(f.plots.len(), 2);
        let Content::Intervals { rows, method, .. } = &f.plots[0].content else {
            panic!()
        };
        assert_eq!(rows[0].label, "12: Occipital");
        assert_eq!(rows[0].values, Some([0.18, -0.04, 0.34]));
        assert_eq!(rows[1].values, Some([0.68, 0.51, 0.79]));
        assert_eq!(method, "Subject bootstrap · 95% CI (percentile)");
        assert_eq!(f.plots[0].bounds.x, f.plots[1].bounds.x);
        assert!(f.plots[0].bounds.x[0] < -0.04 && f.plots[0].bounds.x[1] > 0.82);
        let hover = f.plots[0].hover(0.18, 1.5);
        assert!(
            hover.contains("ROI 12")
                && hover.contains("[-0.04000, 0.34000]")
                && hover.contains("95% CI")
        );
        assert!(f.plots[0].hover(0., 2.).is_empty());
        let file = fixture("regression");
        let f = figure(&options(
            &file,
            &["-rsa-plots", "intervals", "-rsa-rois", "9,2"],
        ))
        .unwrap();
        let Content::Intervals { rows, method, .. } = &f.plots[0].content else {
            panic!()
        };
        assert_eq!(rows[0].values, Some([1.4, 0.8, 1.8])); // No tanh or clipping to correlation range.
        assert_eq!(method, "Dual bootstrap · 90% CI (corrected-variance t)");
        assert!(f.plots[0].xlabel.contains("beta"));
        assert!(
            figure(&options(
                &file,
                &["-rsa-plots", "intervals", "-rsa-effect", "pr"]
            ))
            .unwrap_err()
            .to_string()
            .contains("not partial r")
        );
    }
    #[test]
    fn interval_method_is_explicit_when_ambiguous_and_level_is_never_assumed() {
        let d = tempfile::tempdir().unwrap();
        let path = d.path().join("methods.1D");
        let data = "#ROI label nvox M_r M_p M_q M_bootLo M_bootHi M_cbootLo M_cbootHi\n1 A 20 .4 1 1 .1 .5 -.2 .8\n";
        std::fs::write(&path, format!("# bootstrap: 1000 percentile CI: 95%\n# cond_bootstrap: 1000 groups: 6 usable: 999 percentile CI: 80%\n{data}")).unwrap();
        let file = path.to_str().unwrap();
        assert!(
            figure(&options(file, &["-rsa-plots", "intervals"]))
                .unwrap_err()
                .to_string()
                .contains("multiple CI methods")
        );
        for (method, bounds, level) in [
            ("boot", [0.4, 0.1, 0.5], "95%"),
            ("cboot", [0.4, -0.2, 0.8], "80%"),
        ] {
            let f = figure(&options(
                file,
                &["-rsa-plots", "intervals", "-rsa-ci", method],
            ))
            .unwrap();
            let Content::Intervals { rows, method, .. } = &f.plots[0].content else {
                panic!()
            };
            assert_eq!(rows[0].values, Some(bounds));
            assert!(method.contains(level));
        }
        for extra in [
            vec!["-rsa-plots", "intervals", "-rsa-ci", "dual"],
            vec!["-rsa-ci", "boot"],
            vec!["-rsa-plots", "intervals", "-rsa-ci", "unknown"],
        ] {
            assert!(figure(&options(file, &extra)).is_err());
        }
        std::fs::write(
            &path,
            format!("# loo bootstrap: 1000 percentile CI: 99%\n{data}"),
        )
        .unwrap();
        let f = figure(&options(
            file,
            &["-rsa-plots", "intervals", "-rsa-ci", "boot"],
        ))
        .unwrap();
        assert!(f.plots[0].note.contains("level not recorded"));
        for metadata in [
            "# bootstrap: 10 percentile CI: 100%",
            "# bootstrap: 10 percentile CI: NaN%",
            "# bootstrap: 10 percentile CI: bad",
            "# bootstrap: 10 percentile CI: 95%\n# bootstrap: 10 percentile CI: 90%",
        ] {
            std::fs::write(&path, format!("{metadata}\n{data}")).unwrap();
            assert!(
                figure(&options(
                    file,
                    &["-rsa-plots", "intervals", "-rsa-ci", "boot"]
                ))
                .is_err()
            );
        }
    }
    #[test]
    fn missing_ci_triplets_retain_identity_without_fabricating_points() {
        let d = tempfile::tempdir().unwrap();
        let path = d.path().join("missing.1D");
        let data = "#ROI label nvox M_r M_p M_q M_bootLo M_bootHi\n1 good 20 .4 1 1 .1 .5\n2 no_lo 20 .3 1 1 NA .7\n3 no_effect 20 NA 1 1 -.2 .8\n4 no_hi 20 .1 1 1 -.1 NA\n";
        std::fs::write(&path, data).unwrap();
        let file = path.to_str().unwrap();
        assert!(
            figure(&options(
                file,
                &["-rsa-plots", "intervals", "-rsa-rois", "1,2"]
            ))
            .unwrap_err()
            .to_string()
            .contains("incomplete")
        );
        let f = figure(&options(
            file,
            &["-rsa-plots", "intervals", "-missing", "gap"],
        ))
        .unwrap();
        let Content::Intervals { rows, .. } = &f.plots[0].content else {
            panic!()
        };
        assert_eq!(rows.len(), 4);
        assert_eq!(rows[0].values, Some([0.4, 0.1, 0.5]));
        assert!(rows[1..].iter().all(|r| r.values.is_none()));
        assert!(f.plots[0].hover(0., 2.5).contains("ROI 2 · no_lo"));
        assert!(
            figure(&options(
                file,
                &[
                    "-rsa-plots",
                    "intervals",
                    "-missing",
                    "gap",
                    "-rsa-rois",
                    "2,3,4"
                ]
            ))
            .unwrap_err()
            .to_string()
            .contains("no complete")
        );
    }
    #[test]
    fn forty_interval_rows_render_at_automatic_height_and_larger_selections_fail() {
        let d = tempfile::tempdir().unwrap();
        let path = d.path().join("many.1D");
        let mut text = "#ROI label nvox M_r M_p M_q M_bootLo M_bootHi\n".to_owned();
        for id in 1..=41 {
            text.push_str(&format!("{id} ROI_{id} 10 .4 1 1 .1 .5\n"));
        }
        std::fs::write(&path, text).unwrap();
        let file = path.to_str().unwrap();
        assert!(
            figure(&options(file, &["-rsa-plots", "intervals"]))
                .unwrap_err()
                .to_string()
                .contains("40 ROIs")
        );
        let ids = (1..=40)
            .map(|i| i.to_string())
            .collect::<Vec<_>>()
            .join(",");
        let f = figure(&options(
            file,
            &["-rsa-plots", "intervals", "-rsa-rois", &ids],
        ))
        .unwrap();
        let renderer = crate::render::Renderer::new();
        let h = minimum_height(&f, 1200);
        let scene = renderer.scene(&f, &[f.plots[0].bounds], 1200, h).unwrap();
        assert!(scene.areas[0].h / 40. >= 13.);
        assert!(renderer.scene(&f, &[f.plots[0].bounds], 1200, 500).is_err());
    }
    #[test]
    fn interval_recipe_freezes_auto_method_and_replays_current_view_pixels() {
        let file = fixture("ordinary");
        let o = options(&file, &["-rsa-plots", "intervals", "-rsa-ci", "auto"]);
        let f = figure(&o).unwrap();
        let origin = crate::document::Origin::capture(&o).unwrap();
        let mut views: Vec<_> = f.plots.iter().map(|p| p.bounds).collect();
        views[0] = Bounds {
            x: [0.1, 0.6],
            y: [1., 6.],
        };
        let state = crate::document::Snapshot {
            presentation: crate::document::Presentation::capture(&f),
            views,
            width: 1500,
            height: 1000,
            full_plot: false,
        };
        let d = tempfile::tempdir().unwrap();
        let recipe = d.path().join("intervals.json");
        crate::document::save(&recipe, &origin, &state, false).unwrap();
        let json = std::fs::read_to_string(&recipe).unwrap();
        assert!(
            json.contains("-rsa-ci") && json.contains("\"boot\"") && !json.contains("\"auto\"")
        );
        let loaded = crate::document::load(&recipe).unwrap();
        let r = crate::render::Renderer::new();
        let pixels = |fig: &Figure| {
            r.raster(
                &r.scene(fig, &state.views, state.width, state.height)
                    .unwrap(),
            )
            .unwrap()
        };
        assert_eq!(pixels(&f).data(), pixels(&loaded.figure).data());
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
