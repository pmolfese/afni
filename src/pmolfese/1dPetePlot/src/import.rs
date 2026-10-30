//! Bounded local table inspection and deterministic plot recipes. No inference or sampling.
use crate::{
    cli,
    document::{self, Origin, Roles},
    model::{self, Figure},
    render::Renderer,
};
use anyhow::{Context, Result, bail, ensure};
use std::{
    io::{BufRead, BufReader, Read},
    path::{Path, PathBuf},
    sync::{
        Arc,
        atomic::{AtomicU64, Ordering},
        mpsc,
    },
};

pub const MAX_BYTES: u64 = 16 * 1024 * 1024;
pub const MAX_CELLS: usize = 1_000_000;
#[derive(Clone, Debug)]
pub struct Column {
    pub name: String,
    pub numeric: bool,
    pub missing: usize,
    pub invalid: bool,
    pub role: String,
}
#[derive(Clone, Debug)]
pub struct Table {
    pub path: PathBuf,
    pub fingerprint: String,
    pub columns: Vec<Column>,
    pub rows: usize,
    pub preview: Vec<Vec<String>>,
    pub rsa: Option<crate::rsa::Table>,
}
impl Table {
    pub fn long(&self) -> Option<&crate::rsa_long::LongTable> {
        self.rsa.as_ref().and_then(|t| t.long.as_deref())
    }
    pub fn numeric(&self) -> Vec<usize> {
        self.columns
            .iter()
            .enumerate()
            .filter(|(_, c)| c.numeric && !c.invalid)
            .map(|(i, _)| i)
            .collect()
    }
    pub fn square(&self) -> bool {
        if self.rsa.is_some() {
            self.rows <= 40
        } else {
            self.rows == self.numeric().len() && self.rows <= 256
        }
    }
}
fn cancelled(check: &impl Fn() -> bool) -> Result<()> {
    ensure!(!check(), "Import cancelled");
    Ok(())
}
pub fn inspect(path: &Path, cancel: impl Fn() -> bool) -> Result<Table> {
    cancelled(&cancel)?;
    let metadata =
        std::fs::metadata(path).with_context(|| format!("opening {}", path.display()))?;
    ensure!(
        metadata.is_file(),
        "choose a regular .1D or headered TSV file"
    );
    ensure!(
        metadata.len() <= MAX_BYTES,
        "interactive import supports files up to 16 MiB; select a smaller table or use the CLI"
    );
    let name = path.to_str().context("input path must be UTF-8")?;
    ensure!(
        !name.ends_with(['\'', ']', '}']),
        "rename this file: its name ends with AFNI selector syntax"
    );
    let tsv = name.to_ascii_lowercase().ends_with(".tsv");
    ensure!(
        tsv || name.to_ascii_lowercase().ends_with(".1d"),
        "supported imports are numeric .1D and headered .tsv files"
    );
    let fingerprint = document::fingerprint(path)?;
    if !tsv {
        // Detect the semantic header before the generic numeric .1D reader.
        let mut rsa_header = false;
        for line in BufReader::new(std::fs::File::open(path)?.take(MAX_BYTES + 1)).lines() {
            cancelled(&cancel)?;
            let line = line?;
            let line = line.trim();
            if line
                .strip_prefix('#')
                .is_some_and(|c| c.split_whitespace().next() == Some("ROI"))
            {
                rsa_header = true;
                break;
            }
            if !line.is_empty() && !line.starts_with('#') {
                break;
            }
        }
        if rsa_header {
            let rsa = crate::rsa::read(path)?;
            ensure!(
                rsa.columns.len() <= 256,
                "interactive RSA import supports at most 256 columns; use -rsa in the CLI"
            );
            cancelled(&cancel)?;
            ensure!(
                document::fingerprint(path)? == fingerprint,
                "source changed while reading; retry import"
            );
            let columns = rsa
                .columns
                .iter()
                .enumerate()
                .map(|(i, name)| {
                    let effect = rsa
                        .effects
                        .iter()
                        .find(|(_, index)| **index == i)
                        .map(|((_, kind), _)| kind.as_str())
                        .or_else(|| crate::rsa_extra::family(&rsa, i));
                    Column {
                        name: name.clone(),
                        numeric: effect.is_some(),
                        invalid: false,
                        missing: if let Some(long) = &rsa.long {
                            if effect.is_some() {
                                long.rows.values().filter(|r| !r[i].is_finite()).count()
                            } else {
                                0
                            }
                        } else if i < 3 {
                            0
                        } else {
                            rsa.rows
                                .values()
                                .filter(|r| !r.values[i].is_finite())
                                .count()
                        },
                        role: effect
                            .map(|kind| format!("effect: {kind}"))
                            .unwrap_or_else(|| {
                                if rsa.long.is_some() && i < 5 {
                                    "scientific identity"
                                } else if i < 3 {
                                    "ROI metadata"
                                } else {
                                    "auxiliary statistic"
                                }
                                .into()
                            }),
                    }
                })
                .collect();
            let preview = rsa
                .rows
                .values()
                .take(8)
                .map(|r| {
                    let mut values = vec![r.id.to_string(), r.label.clone(), r.nvox.to_string()];
                    values.extend(r.values[3..].iter().map(f64::to_string));
                    values
                })
                .collect();
            let preview = rsa
                .long
                .as_ref()
                .map(|l| l.preview.clone())
                .unwrap_or(preview);
            let parent = path
                .parent()
                .filter(|p| !p.as_os_str().is_empty())
                .unwrap_or(Path::new("."));
            return Ok(Table {
                path: std::fs::canonicalize(parent)?.join(path.file_name().unwrap()),
                fingerprint,
                columns,
                rows: rsa
                    .long
                    .as_ref()
                    .map(|l| l.rows.len())
                    .unwrap_or(rsa.rows.len()),
                preview,
                rsa: Some(rsa),
            });
        }
    }
    let file = std::fs::File::open(path)?;
    let mut columns: Vec<Column> = vec![];
    let mut preview = vec![];
    let mut rows = 0usize;
    let mut bytes = 0usize;
    for (line_no, line) in BufReader::new(file.take(MAX_BYTES + 1)).lines().enumerate() {
        cancelled(&cancel)?;
        let line = line?;
        bytes += line.len() + 1;
        ensure!(
            bytes <= MAX_BYTES as usize,
            "file grew beyond the 16 MiB import limit"
        );
        let data = line.split('#').next().unwrap_or("");
        if data.trim().is_empty() {
            continue;
        }
        let cells: Vec<_> = if tsv {
            data.trim_end_matches('\r')
                .split('\t')
                .take(257)
                .map(str::trim)
                .collect()
        } else {
            data.split_whitespace().take(257).collect()
        };
        ensure!(
            cells.len() <= 256,
            "interactive import supports at most 256 columns"
        );
        if columns.is_empty() {
            let names: Vec<_> = if tsv {
                cells.iter().map(|s| s.to_string()).collect()
            } else {
                (0..cells.len())
                    .map(|c| format!("Column {}", c + 1))
                    .collect()
            };
            let unique: std::collections::BTreeSet<_> = names.iter().collect();
            ensure!(
                unique.len() == names.len()
                    && names.iter().all(|s| !s.is_empty() && s.len() <= 512),
                "TSV column names must be nonempty, unique, and at most 512 bytes"
            );
            columns = names
                .into_iter()
                .map(|name| Column {
                    name,
                    numeric: true,
                    missing: 0,
                    invalid: false,
                    role: String::new(),
                })
                .collect();
            if tsv {
                continue;
            }
        }
        ensure!(
            cells.len() == columns.len(),
            "row {} has {} fields; expected {}",
            line_no + 1,
            cells.len(),
            columns.len()
        );
        rows += 1;
        ensure!(
            rows <= 100_000 && rows * columns.len() <= MAX_CELLS,
            "interactive import supports up to 100000 rows / 1000000 cells; use a smaller table or the CLI"
        );
        for (column, cell) in columns.iter_mut().zip(&cells) {
            if ["", "na", "n/a", "nan"].contains(&cell.to_ascii_lowercase().as_str()) {
                column.missing += 1;
            } else {
                match cell.parse::<f64>() {
                    Ok(v) if v.is_infinite() => column.invalid = true,
                    Ok(v) if v.is_nan() => column.missing += 1,
                    Ok(_) => {}
                    Err(_) if tsv => column.numeric = false,
                    Err(_) => bail!(
                        "row {}: .1D files must be numeric; use headered TSV for text columns",
                        line_no + 1
                    ),
                }
            }
        }
        if preview.len() < 8 {
            preview.push(
                cells
                    .iter()
                    .map(|s| s.chars().take(160).collect())
                    .collect(),
            );
        }
    }
    ensure!(
        rows > 0 && !columns.is_empty(),
        "input contains no data rows"
    );
    cancelled(&cancel)?;
    ensure!(
        document::fingerprint(path)? == fingerprint,
        "source changed while reading; retry import"
    );
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    Ok(Table {
        path: std::fs::canonicalize(parent)?.join(path.file_name().unwrap()),
        fingerprint,
        columns,
        rows,
        preview,
        rsa: None,
    })
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Line,
    Histogram,
    Scatter,
    Heat,
    Intervals,
    Rdm,
}
impl Kind {
    pub const ALL: [Self; 6] = [
        Self::Line,
        Self::Histogram,
        Self::Scatter,
        Self::Heat,
        Self::Intervals,
        Self::Rdm,
    ];
    pub fn name(self) -> &'static str {
        match self {
            Self::Line => "Lines",
            Self::Histogram => "Distribution",
            Self::Scatter => "Scatter",
            Self::Heat => "Heatmap",
            Self::Intervals => "Confidence intervals",
            Self::Rdm => "RDM / named matrix",
        }
    }
    pub fn explanation(self) -> &'static str {
        match self {
            Self::Line => "Numeric columns as curves; choose x or use row index.",
            Self::Histogram => {
                "Distribution of selected values; counts use all finite observations."
            }
            Self::Scatter => {
                "Pair selected y columns with an explicit x column from the same rows."
            }
            Self::Heat => "Square numeric matrix; shape alone does not imply an RDM.",
            Self::Rdm => {
                "A symmetric pairwise matrix with an explicit condition-label file; values are preserved."
            }
            Self::Intervals => {
                "RSA effects with supplied asymmetric confidence bounds, one row per ROI."
            }
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct RdmMapping {
    pub labels: String,
    pub order: String,
    pub sense: String,
    pub clim: String,
}
impl Default for RdmMapping {
    fn default() -> Self {
        Self {
            labels: String::new(),
            order: String::new(),
            sense: "unspecified".into(),
            clim: String::new(),
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct Mapping {
    pub kind: Kind,
    pub x: Option<usize>,
    pub y: Vec<usize>,
    pub missing: String,
    pub roi_filter: String,
    pub ci: String,
    pub rsa_models: String,
    pub rsa_family: String,
    pub rsa_steps: String,
    pub rsa_band: bool,
    pub presentation: Option<document::Panel>,
    pub reset_style: bool,
    pub rdm: RdmMapping,
    pub attachments: crate::attachments::Attachments,
}
impl Mapping {
    pub fn for_kind(table: &Table, kind: Kind) -> Self {
        let mut n = table.numeric();
        if let Some(rsa) = &table.rsa {
            let primary = if rsa.effects.keys().any(|(_, kind)| kind == "b") {
                "b"
            } else if rsa.long.is_some() {
                "effect"
            } else {
                "r"
            };
            n.retain(|i| {
                rsa.effects
                    .iter()
                    .any(|((_, kind), index)| kind == primary && index == i)
            });
        }
        Self {
            kind,
            x: if kind == Kind::Scatter && table.long().is_none() {
                n.first().copied()
            } else {
                None
            },
            y: if table.long().is_some() {
                n.into_iter().take(1).collect()
            } else {
                match kind {
                    Kind::Heat | Kind::Rdm => n,
                    Kind::Scatter => n.into_iter().skip(1).take(1).collect(),
                    _ => n.into_iter().take(1).collect(),
                }
            },
            missing: "error".into(),
            roi_filter: table
                .rsa
                .as_ref()
                .filter(|t| t.long.is_some())
                .and_then(|t| t.rows.keys().next())
                .map(i32::to_string)
                .unwrap_or_default(),
            ci: "auto".into(),
            rsa_models: String::new(),
            rsa_family: "effects".into(),
            rsa_steps: String::new(),
            rsa_band: false,
            presentation: None,
            reset_style: false,
            rdm: Default::default(),
            attachments: crate::attachments::Attachments {
                include_bounds: true,
                threshold_percent: 5,
                ..Default::default()
            },
        }
    }
    pub fn allows(&self, table: &Table, i: usize) -> bool {
        table.rsa.as_ref().is_none_or(|t| {
            if self.rsa_family == "effects" {
                t.effects.values().any(|v| *v == i)
            } else {
                crate::rsa_extra::family(t, i) == Some(self.rsa_family.as_str())
            }
        })
    }
    pub fn set_rsa_family(&mut self, table: &Table, family: &str) {
        self.rsa_family = family.into();
        if family == "effects" {
            let original = Self::for_kind(table, self.kind);
            self.x = original.x;
            self.y = original.y;
            return;
        }
        let n: Vec<_> = table
            .numeric()
            .into_iter()
            .filter(|i| {
                table.rsa.as_ref().is_some_and(|t| {
                    if family == "effects" {
                        t.effects.values().any(|v| v == i)
                    } else {
                        crate::rsa_extra::family(t, *i) == Some(family)
                    }
                })
            })
            .collect();
        self.x = if self.kind == Kind::Scatter {
            n.first().copied()
        } else {
            None
        };
        self.y = match self.kind {
            Kind::Heat => n,
            Kind::Scatter => n.into_iter().skip(1).take(1).collect(),
            _ => n.into_iter().take(1).collect(),
        };
    }
    pub fn restore(
        table: &Table,
        o: &cli::Options,
        roles: &Roles,
        presentation: document::Panel,
    ) -> Result<Self> {
        ensure!(
            roles.columns
                == table
                    .columns
                    .iter()
                    .map(|c| c.name.clone())
                    .collect::<Vec<_>>(),
            "source columns changed; reload/relink or reimport before editing roles"
        );
        let kind = Kind::ALL
            .into_iter()
            .find(|k| k.name() == roles.kind)
            .context("unknown saved graph type")?;
        ensure!(
            table.rsa.is_none() && matches!(kind, Kind::Line | Kind::Scatter | Kind::Histogram),
            "attachment editing requires an ordinary line, scatter or distribution graph"
        );
        let mut m = Self::for_kind(table, kind);
        m.x = roles.x;
        m.y = roles.y.clone();
        m.missing = roles.missing.clone();
        m.presentation = Some(presentation);
        let column = |s: &str| -> Result<usize> {
            s.parse::<usize>()
                .ok()
                .filter(|i| *i < table.columns.len())
                .or_else(|| table.columns.iter().position(|c| c.name == s))
                .context("saved column role no longer exists")
        };
        let selected = |spec: &str| -> Result<Option<Vec<usize>>> {
            let (name, suffix) = document::split_spec(spec);
            if Path::new(name) != table.path || suffix.is_empty() {
                return Ok(None);
            }
            ensure!(
                suffix.starts_with('[') && suffix.ends_with(']'),
                "unsupported saved attachment selection"
            );
            Ok(Some(
                suffix[1..suffix.len() - 1]
                    .split(',')
                    .map(column)
                    .collect::<Result<_>>()?,
            ))
        };
        m.attachments.group = o.group_by.as_deref().map(column).transpose()?;
        m.attachments.include_bounds = o.qc.use_pm;
        ensure!(
            o.qc.pm.len() <= 1
                && o.qc.lower.len() <= 1
                && o.qc.upper.len() <= 1
                && o.qc.censor_files.len() <= 1,
            "edit one attachment file per role"
        );
        if let Some(spec) = o.qc.pm.first() {
            if let Some(cols) = selected(spec)? {
                ensure!(cols.len() == m.y.len(), "width/Y mapping mismatch");
                m.attachments.uncertainty = "widths".into();
                for (y, c) in m.y.iter().zip(cols) {
                    m.attachments.bands.insert(*y, [Some(c), None]);
                }
            } else {
                m.attachments.uncertainty = "file".into();
                m.attachments.width_file = spec.clone();
            }
        } else if let Some(spec) = o.qc.lower.first() {
            let lo =
                selected(spec)?.context("editing external asymmetric bounds requires reimport")?;
            let hi = selected(o.qc.upper.first().context("missing upper bounds")?)?
                .context("upper bounds must be columns in the graph's source")?;
            ensure!(
                lo.len() == m.y.len() && hi.len() == m.y.len(),
                "bound/Y mapping mismatch"
            );
            m.attachments.uncertainty = "bounds".into();
            for ((y, l), h) in m.y.iter().zip(lo).zip(hi) {
                m.attachments.bands.insert(*y, [Some(l), Some(h)]);
            }
        }
        if let Some(spec) = o.qc.censor_files.first() {
            if let Some(cols) = selected(spec)? {
                ensure!(cols.len() == 1, "censor needs one column");
                m.attachments.censor_column = Some(cols[0]);
            } else {
                m.attachments.censor_file = spec.clone();
            }
        }
        m.attachments.xmat = o.shade_design.file.clone().unwrap_or_default();
        m.attachments.palette = o.shade_design.colors.join(" ");
        m.attachments.threshold_percent = (o.shade_design.threshold * 100.).round() as u32;
        // Opening an editor is not renewed evidence that sidecar rows align.
        m.attachments.aligned = false;
        Ok(m)
    }
    pub fn roles(&self, table: &Table) -> Roles {
        Roles {
            columns: table.columns.iter().map(|c| c.name.clone()).collect(),
            x: self.x,
            y: self.y.clone(),
            kind: self.kind.name().into(),
            missing: self.missing.clone(),
        }
    }
    fn rsa_arguments(&self, table: &Table, rsa: &crate::rsa::Table) -> Result<Vec<String>> {
        if self.rsa_family != "effects" {
            ensure!(
                rsa.long.is_none(),
                "choose primary effects for a long table"
            );
            ensure!(!self.y.is_empty(), "choose a scientific quantity");
            let mode = match self.kind {
                Kind::Histogram => "histogram",
                Kind::Scatter => "scatter",
                Kind::Heat => "heatmap",
                Kind::Intervals => "intervals",
                _ => anyhow::bail!("choose distribution, scatter, heatmap or intervals"),
            };
            let mut indices = self.y.clone();
            if self.kind == Kind::Scatter {
                ensure!(indices.len() == 1, "select one Y quantity");
                indices.insert(0, self.x.context("choose an X quantity")?);
            }
            if matches!(self.kind, Kind::Histogram | Kind::Intervals) {
                ensure!(
                    indices.len() == 1,
                    "select one quantity per graph; add other quantities separately"
                );
            }
            ensure!(
                indices
                    .iter()
                    .all(|i| crate::rsa_extra::family(rsa, *i) == Some(self.rsa_family.as_str())),
                "quantities must share one scientific family; probabilities and weights cannot substitute for effects"
            );
            let mut args = vec![
                "-rsa".into(),
                table.path.to_str().context("path must be UTF-8")?.into(),
                "-rsa-family".into(),
                self.rsa_family.clone(),
                "-rsa-measures".into(),
                indices
                    .iter()
                    .map(|i| rsa.columns[*i].clone())
                    .collect::<Vec<_>>()
                    .join(","),
                "-rsa-plots".into(),
                mode.into(),
                "-missing".into(),
                self.missing.clone(),
            ];
            if !self.roi_filter.trim().is_empty() {
                args.extend(["-rsa-rois".into(), self.roi_filter.trim().into()]);
            }
            if self.kind == Kind::Intervals {
                args.extend(["-rsa-ci".into(), self.ci.clone()]);
            }
            return Ok(args);
        }
        if let Some(long) = &rsa.long {
            ensure!(
                self.y.len() == 1,
                "choose one effect kind for this long-table graph"
            );
            let effect = rsa
                .effects
                .iter()
                .find(|(_, i)| **i == self.y[0])
                .map(|((_, kind), _)| kind)
                .context("select an effect, not a statistic")?;
            let mode = match self.kind {
                Kind::Line if long.time => "curves",
                Kind::Scatter if !long.time => "points",
                Kind::Heat => "heatmap",
                Kind::Intervals if long.time => "intervals",
                _ => anyhow::bail!("choose a plot compatible with this long RSA schema"),
            };
            let rois = crate::rsa::names(self.roi_filter.trim()).context(
                "choose one ROI ID per imported graph; the CLI can make multiple panels",
            )?;
            ensure!(
                rois.len() == 1,
                "choose one ROI ID per imported graph; add other ROIs separately"
            );
            let mut args = vec![
                "-rsa".into(),
                table.path.to_str().context("path must be UTF-8")?.into(),
                "-rsa-effect".into(),
                effect.clone(),
                "-rsa-plots".into(),
                mode.into(),
                "-rsa-rois".into(),
                rois[0].clone(),
                "-missing".into(),
                self.missing.clone(),
            ];
            if !self.rsa_models.trim().is_empty() {
                args.extend(["-rsa-models".into(), self.rsa_models.trim().into()]);
            }
            if !self.rsa_steps.trim().is_empty() {
                args.extend([
                    if long.time {
                        "-rsa-times"
                    } else {
                        "-rsa-summaries"
                    }
                    .into(),
                    self.rsa_steps.trim().into(),
                ]);
            }
            if self.kind == Kind::Intervals || self.kind == Kind::Line && self.rsa_band {
                let method = if self.ci == "auto" && self.kind == Kind::Line {
                    let methods: Vec<_> = ["boot", "cboot", "dual"]
                        .into_iter()
                        .filter(|m| rsa.columns.contains(&format!("{m}Lo")))
                        .collect();
                    ensure!(
                        methods.len() == 1,
                        "choose a supplied CI method; auto requires exactly one"
                    );
                    methods[0].to_string()
                } else {
                    self.ci.clone()
                };
                args.extend(["-rsa-ci".into(), method]);
            }
            return Ok(args);
        }
        ensure!(
            !matches!(self.kind, Kind::Line | Kind::Rdm),
            "RSA effects use distribution, paired scatter, or ROI-by-model heatmap"
        );
        ensure!(
            self.missing != "mean",
            "RSA effects are not imputed; choose reject or gap"
        );
        ensure!(!self.y.is_empty(), "select RSA effect columns");
        let mut indices = vec![];
        if self.kind == Kind::Scatter {
            indices.push(self.x.context("scatter needs an X effect")?);
            ensure!(
                self.y.len() == 1 && Some(self.y[0]) != self.x,
                "select one Y effect different from X"
            );
        } else if matches!(self.kind, Kind::Histogram | Kind::Intervals) {
            ensure!(
                self.y.len() == 1,
                "select one model effect per graph; add other models separately"
            );
        }
        indices.extend(&self.y);
        let roles: Vec<_> = indices
            .iter()
            .map(|index| {
                rsa.effects
                    .iter()
                    .find(|(_, i)| *i == index)
                    .map(|(role, _)| role)
                    .context("select an r, b, or pr effect, not metadata/statistics")
            })
            .collect::<Result<_>>()?;
        let effect = &roles[0].1;
        ensure!(
            roles.iter().all(|(_, kind)| kind == effect),
            "selected RSA effects must share one kind: r, b, or pr"
        );
        let models = roles
            .iter()
            .map(|(model, _)| model.as_str())
            .collect::<Vec<_>>()
            .join(",");
        let mut args = vec![
            "-rsa".into(),
            table.path.to_str().context("path must be UTF-8")?.into(),
            "-rsa-models".into(),
            models,
            "-rsa-effect".into(),
            effect.clone(),
            "-rsa-plots".into(),
            match self.kind {
                Kind::Histogram => "histogram",
                Kind::Scatter => "scatter",
                Kind::Heat => "heatmap",
                Kind::Intervals => "intervals",
                _ => unreachable!(),
            }
            .into(),
            "-missing".into(),
            self.missing.clone(),
        ];
        if !self.roi_filter.trim().is_empty() {
            args.extend(["-rsa-rois".into(), self.roi_filter.trim().into()]);
        }
        if self.kind == Kind::Intervals {
            args.extend(["-rsa-ci".into(), self.ci.clone()]);
        }
        Ok(args)
    }
    pub fn arguments(&self, table: &Table) -> Result<Vec<String>> {
        let attachments =
            self.attachments
                .arguments(table, self.kind, self.x, &self.y, &self.missing)?;
        if self.kind == Kind::Rdm {
            ensure!(
                table.rsa.is_none(),
                "RSA effect tables cannot reconstruct an RDM"
            );
            ensure!(
                !self.rdm.labels.trim().is_empty(),
                "choose a condition-label file, one label per matrix row/column"
            );
            let mut args = vec![
                "-rdm".into(),
                table.path.to_str().context("path must be UTF-8")?.into(),
                "-rdm-labels".into(),
                self.rdm.labels.clone(),
                "-rdm-kind".into(),
                self.rdm.sense.clone(),
                "-missing".into(),
                self.missing.clone(),
            ];
            if !self.rdm.order.trim().is_empty() {
                args.extend(["-rdm-order".into(), self.rdm.order.clone()]);
            }
            if !self.rdm.clim.trim().is_empty() {
                args.extend(["-clim".into(), self.rdm.clim.trim().into()]);
            }
            return Ok(args);
        }
        if let Some(rsa) = &table.rsa {
            return self.rsa_arguments(table, rsa);
        }
        ensure!(
            self.kind != Kind::Intervals,
            "confidence intervals require an ordinary 3dRSA ROI table with supplied bounds"
        );
        ensure!(!self.y.is_empty(), "select at least one numeric y column");
        ensure!(
            self.y.len() <= if self.kind == Kind::Heat { 256 } else { 16 },
            "select at most 16 series per graph"
        );
        let mut sorted = self.y.clone();
        sorted.sort_unstable();
        sorted.dedup();
        ensure!(sorted.len() == self.y.len(), "duplicate y columns");
        ensure!(
            self.y.iter().chain(self.x.iter()).all(|&c| table
                .columns
                .get(c)
                .is_some_and(|c| c.numeric && !c.invalid)),
            "selected columns must be numeric without infinity"
        );
        ensure!(
            ["error", "gap", "mean"].contains(&self.missing.as_str()),
            "choose a missing-value policy"
        );
        if let Some(x) = self.x {
            ensure!(!self.y.contains(&x), "choose different x and y columns");
            ensure!(
                table.columns[x].missing == 0,
                "x coordinates must be finite; choose another column or row index"
            );
        }
        if self.kind == Kind::Scatter {
            ensure!(self.x.is_some(), "scatter needs an x column");
        }
        if self.kind == Kind::Heat {
            ensure!(
                self.x.is_none() && self.y.len() == table.rows && table.rows <= 256,
                "heatmap needs a square selection of at most 256 rows/columns"
            );
        }
        if matches!(self.kind, Kind::Line | Kind::Scatter) {
            ensure!(
                table.rows * self.y.len() <= 50_000,
                "interactive line/scatter previews support 50000 points; select fewer columns, use a histogram, or use the CLI"
            );
        }
        let file = table.path.to_str().context("path must be UTF-8")?;
        let selector = self
            .y
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join(",");
        let basename = table.path.file_name().unwrap().to_string_lossy();
        let title = format!("{} · {}", basename, self.kind.name());
        let mut args = vec![
            "-infiles".into(),
            format!("{file}[{selector}]"),
            "-missing".into(),
            self.missing.clone(),
            "-title".into(),
            title.clone(),
            "-panel-title".into(),
            "1".into(),
            title,
        ];
        match self.kind {
            Kind::Line => {}
            Kind::Histogram => args.push("-histogram".into()),
            Kind::Scatter => args.push("-scatter".into()),
            Kind::Heat => args.push("-heat".into()),
            Kind::Intervals | Kind::Rdm => unreachable!(),
        }
        if self.kind != Kind::Heat {
            args.extend(["-one".into(), "-legend_on".into(), "-legend_labels".into()]);
            // CLI list values starting with '-' need a harmless display prefix.
            args.extend(self.y.iter().map(|&c| label(&table.columns[c].name)));
        }
        if let Some(x) = self.x {
            args.extend([
                "-xfile".into(),
                format!("{file}[{x}]"),
                "-xlabel".into(),
                label(&table.columns[x].name),
            ]);
        } else if self.kind == Kind::Heat {
            args.extend([
                "-xlabel".into(),
                "Column index".into(),
                "-ylabels".into(),
                "Row index".into(),
            ]);
        } else if self.kind == Kind::Histogram && self.y.len() == 1 {
            args.extend(["-xlabel".into(), label(&table.columns[self.y[0]].name)]);
        } else if self.kind == Kind::Line {
            args.extend(["-xlabel".into(), "Row index (0-based)".into()]);
        }
        args.extend(attachments);
        Ok(args)
    }
}
fn label(s: &str) -> String {
    if s.starts_with('-') {
        format!("Column {s}")
    } else {
        s.into()
    }
}
pub struct Prepared {
    pub figure: Figure,
    pub origin: Origin,
    pub roles: Roles,
    pub rgba: Vec<u8>,
    pub size: [usize; 2],
    pub mapping: Mapping,
}
pub fn prepare(table: &Table, mapping: &Mapping, cancel: impl Fn() -> bool) -> Result<Prepared> {
    cancelled(&cancel)?;
    ensure!(
        document::fingerprint(&table.path)? == table.fingerprint,
        "source changed since inspection; retry import"
    );
    mapping.attachments.validate_files()?;
    let options = cli::parse(&mapping.arguments(table)?)?;
    let origin = Origin::capture(&options)?;
    let mut figure = model::from_options(&options)?;
    if let Some(p) = &mapping.presentation {
        ensure!(
            figure.plots.len() == 1,
            "an attachment edit must produce one graph"
        );
        let plot = &mut figure.plots[0];
        if !mapping.reset_style {
            p.style
                .validate(plot)
                .context("style no longer fits the selected series; enable Reset styling")?;
            plot.style = p.style.clone();
        }
        plot.title = p.title.clone();
        plot.caption = p.caption.clone();
        plot.xlabel = p.xlabel.clone();
        plot.ylabel = p.ylabel.clone();
    }
    if matches!(mapping.kind, Kind::Line | Kind::Scatter) {
        let points: usize = figure
            .plots
            .iter()
            .map(|p| match &p.content {
                model::Content::Lines(s) | model::Content::Scatter(s) => {
                    s.iter().map(|s| s.points.len()).sum()
                }
                _ => 0,
            })
            .sum();
        ensure!(
            points <= 50_000,
            "grouped previews support 50000 stored points; select fewer Y columns/groups or use the CLI"
        );
    }
    origin.verify()?;
    cancelled(&cancel)?;
    ensure!(
        document::fingerprint(&table.path)? == table.fingerprint,
        "source changed while preparing preview; retry import"
    );
    let renderer = Renderer::new();
    let height = 600.max(crate::rsa::minimum_height(&figure, 900));
    let scene = renderer.scene(&figure, &[figure.plots[0].bounds], 900, height)?;
    let pix = renderer.raster(&scene)?;
    cancelled(&cancel)?;
    Ok(Prepared {
        figure,
        origin,
        roles: mapping.roles(table),
        rgba: pix.data().to_vec(),
        size: [900, height as usize],
        mapping: mapping.clone(),
    })
}

pub enum Work {
    Inspect(PathBuf),
    Preview(Arc<Table>, Box<Mapping>),
    Confirm(Box<Prepared>),
}
pub enum Answer {
    Table(Arc<Table>),
    Preview(Box<Prepared>),
    Accepted(Box<Prepared>),
}
struct Job {
    generation: u64,
    work: Work,
}
/// One worker, one queued job, one replaceable pending job. Stale results cannot reach the UI.
pub struct Worker {
    sender: mpsc::SyncSender<Job>,
    receiver: mpsc::Receiver<(u64, Result<Answer, String>)>,
    generation: Arc<AtomicU64>,
    pending: Option<Job>,
}
impl Worker {
    pub fn new(wake: impl Fn() + Send + 'static) -> Self {
        let (sender, jobs) = mpsc::sync_channel::<Job>(1);
        let (answers, receiver) = mpsc::channel();
        let generation = Arc::new(AtomicU64::new(0));
        let token = generation.clone();
        std::thread::spawn(move || {
            while let Ok(job) = jobs.recv() {
                let cancel = || token.load(Ordering::Relaxed) != job.generation;
                if cancel() {
                    continue;
                }
                let result = match job.work {
                    Work::Inspect(path) => {
                        inspect(&path, cancel).map(|t| Answer::Table(Arc::new(t)))
                    }
                    Work::Preview(table, mapping) => {
                        prepare(&table, &mapping, cancel).map(|p| Answer::Preview(Box::new(p)))
                    }
                    Work::Confirm(prepared) => {
                        prepared.origin.verify().map(|_| Answer::Accepted(prepared))
                    }
                }
                .map_err(|e| format!("{e:#}"));
                if !cancel() && answers.send((job.generation, result)).is_ok() {
                    wake();
                }
            }
        });
        Self {
            sender,
            receiver,
            generation,
            pending: None,
        }
    }
    pub fn start(&mut self, work: Work) {
        let generation = self.generation.fetch_add(1, Ordering::Relaxed) + 1;
        self.pending = Some(Job { generation, work });
        self.flush();
    }
    fn flush(&mut self) {
        if let Some(job) = self.pending.take()
            && let Err(mpsc::TrySendError::Full(job)) = self.sender.try_send(job)
        {
            self.pending = Some(job);
        }
    }
    pub fn cancel(&mut self) {
        self.generation.fetch_add(1, Ordering::Relaxed);
        self.pending = None;
    }
    pub fn poll(&mut self) -> Option<Result<Answer, String>> {
        self.flush();
        loop {
            match self.receiver.try_recv() {
                Ok((generation, result)) => {
                    if generation == self.generation.load(Ordering::Relaxed) {
                        return Some(result);
                    }
                }
                Err(mpsc::TryRecvError::Empty) => return None,
                Err(mpsc::TryRecvError::Disconnected) => {
                    return Some(Err(
                        "Import worker stopped; reopen the application to retry.".into(),
                    ));
                }
            }
        }
    }
}
impl Drop for Worker {
    fn drop(&mut self) {
        self.cancel();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn data(dir: &Path, name: &str, content: &str) -> PathBuf {
        let path = dir.join(name);
        std::fs::write(&path, content).unwrap();
        path
    }
    #[test]
    fn typed_preview_preserves_text_and_requires_explicit_missing_policy() {
        let d = tempfile::tempdir().unwrap();
        let path = data(
            d.path(),
            "table.tsv",
            "roi\ttime\teffect\nA\t0\t0.2\nB\t1\tNA\nC\t2\t0.6\n",
        );
        let table = inspect(&path, || false).unwrap();
        assert_eq!(table.numeric(), [1, 2]);
        assert_eq!(table.preview[0][0], "A");
        assert_eq!(table.columns[2].missing, 1);
        let mut mapping = Mapping {
            kind: Kind::Line,
            x: Some(1),
            y: vec![2],
            missing: "error".into(),
            roi_filter: String::new(),
            ci: "auto".into(),
            rsa_models: String::new(),
            rsa_family: "effects".into(),
            rsa_steps: String::new(),
            rsa_band: false,
            presentation: None,
            reset_style: false,
            rdm: Default::default(),
            attachments: crate::attachments::Attachments {
                include_bounds: true,
                threshold_percent: 5,
                ..Default::default()
            },
        };
        assert!(prepare(&table, &mapping, || false).is_err());
        mapping.missing = "gap".into();
        let preview = prepare(&table, &mapping, || false).unwrap();
        let model::Content::Lines(series) = &preview.figure.plots[0].content else {
            panic!()
        };
        assert!(series[0].points[1][1].is_nan());
        assert_eq!(series[0].label, "effect");
        assert_eq!(preview.figure.plots[0].xlabel, "time");
        mapping.missing = "mean".into();
        let preview = prepare(&table, &mapping, || false).unwrap();
        let model::Content::Lines(series) = &preview.figure.plots[0].content else {
            panic!()
        };
        assert!((series[0].points[1][1] - 0.4).abs() < 1e-12);
        assert_eq!(
            std::fs::read_to_string(&path).unwrap(),
            "roi\ttime\teffect\nA\t0\t0.2\nB\t1\tNA\nC\t2\t0.6\n"
        );
        mapping.y = vec![0];
        assert!(mapping.arguments(&table).is_err());
    }
    #[test]
    fn four_candidates_use_selected_full_data_and_validate_pairing() {
        let d = tempfile::tempdir().unwrap();
        let path = data(d.path(), "matrix.1D", "0 1 2\n1 0 3\n2 3 0\n");
        let table = inspect(&path, || false).unwrap();
        assert!(table.square());
        for kind in Kind::ALL
            .into_iter()
            .filter(|k| !matches!(k, Kind::Intervals | Kind::Rdm))
        {
            let mapping = Mapping::for_kind(&table, kind);
            let p = prepare(&table, &mapping, || false).unwrap();
            assert_eq!(p.figure.plots.len(), 1);
            assert_eq!(p.size, [900, 600]);
            match &p.figure.plots[0].content {
                model::Content::Intervals { .. } => panic!("generic data has no CI roles"),
                model::Content::Histogram { series, .. } => {
                    assert_eq!(series[0].points.iter().map(|p| p[1]).sum::<f64>(), 3.)
                }
                model::Content::Heat { matrix, .. } => {
                    assert_eq!(matrix[0], [0., 1., 2.]);
                    assert_eq!(p.figure.plots[0].kind(), "HEATMAP");
                }
                model::Content::Scatter(series) => {
                    assert_eq!(series[0].points, [[0., 1.], [1., 0.], [2., 3.]])
                }
                model::Content::Lines(series) => {
                    assert_eq!(series[0].points, [[0., 0.], [1., 1.], [2., 2.]])
                }
            }
        }
        let bad = Mapping {
            kind: Kind::Scatter,
            x: Some(1),
            y: vec![1],
            missing: "error".into(),
            roi_filter: String::new(),
            ci: "auto".into(),
            rsa_models: String::new(),
            rsa_family: "effects".into(),
            rsa_steps: String::new(),
            rsa_band: false,
            presentation: None,
            reset_style: false,
            rdm: Default::default(),
            attachments: crate::attachments::Attachments {
                include_bounds: true,
                threshold_percent: 5,
                ..Default::default()
            },
        };
        assert!(bad.arguments(&table).is_err());
        let non_square = Mapping {
            kind: Kind::Heat,
            x: None,
            y: vec![0, 1],
            missing: "error".into(),
            roi_filter: String::new(),
            ci: "auto".into(),
            rsa_models: String::new(),
            rsa_family: "effects".into(),
            rsa_steps: String::new(),
            rsa_band: false,
            presentation: None,
            reset_style: false,
            rdm: Default::default(),
            attachments: crate::attachments::Attachments {
                include_bounds: true,
                threshold_percent: 5,
                ..Default::default()
            },
        };
        assert!(non_square.arguments(&table).is_err());
    }
    #[test]
    fn malformed_changed_cancelled_and_large_inputs_fail_without_mutation() {
        let d = tempfile::tempdir().unwrap();
        for (name, text) in [
            ("ragged.1D", "1 2\n3\n"),
            ("bad.1D", "1 nope\n2 3\n"),
            ("bad.tsv", "a\ta\n1\t2\n"),
            ("empty.1D", "# empty\n"),
            ("ignored.csv", "1,2\n"),
        ] {
            assert!(inspect(&data(d.path(), name, text), || false).is_err());
        }
        let path = data(d.path(), "valid.1D", "1\n2\n3\n");
        assert!(inspect(&path, || true).is_err());
        let table = inspect(&path, || false).unwrap();
        let mapping = Mapping::for_kind(&table, Kind::Line);
        std::fs::write(&path, "2\n3\n4\n").unwrap();
        assert!(
            prepare(&table, &mapping, || false)
                .err()
                .unwrap()
                .to_string()
                .contains("changed")
        );
        let large = d.path().join("huge.1D");
        let file = std::fs::File::create(&large).unwrap();
        file.set_len(MAX_BYTES + 1).unwrap();
        assert!(
            inspect(&large, || false)
                .unwrap_err()
                .to_string()
                .contains("16 MiB")
        );
        assert!(inspect(d.path(), || false).is_err());
    }
    #[test]
    fn bounded_large_preview_uses_every_point_and_rejects_oversized_selection() {
        let d = tempfile::tempdir().unwrap();
        let text: String = (0..25_001)
            .map(|i| format!("{} {}\n", i % 100, i % 31))
            .collect();
        let path = data(d.path(), "large.1D", &text);
        let start = std::time::Instant::now();
        let table = inspect(&path, || false).unwrap();
        let mut mapping = Mapping::for_kind(&table, Kind::Line);
        let prepared = prepare(&table, &mapping, || false).unwrap();
        let model::Content::Lines(series) = &prepared.figure.plots[0].content else {
            panic!()
        };
        assert_eq!(series[0].points.len(), 25_001);
        eprintln!(
            "Import inspection + 25001-point PNG preview: {:?}",
            start.elapsed()
        );
        mapping.y = vec![0, 1];
        assert!(
            mapping
                .arguments(&table)
                .unwrap_err()
                .to_string()
                .contains("50000")
        );
        mapping.kind = Kind::Histogram;
        let prepared = prepare(&table, &mapping, || false).unwrap();
        let model::Content::Histogram { series, .. } = &prepared.figure.plots[0].content else {
            panic!()
        };
        assert_eq!(series.len(), 2);
        for s in series {
            assert_eq!(s.points.iter().map(|p| p[1]).sum::<f64>(), 25_001.);
        }
    }
    #[test]
    fn worker_cancellation_and_replacement_deliver_only_the_current_request() {
        let d = tempfile::tempdir().unwrap();
        let first = data(d.path(), "first.1D", "0\n1\n2\n");
        let last = data(d.path(), "last.1D", "3\n4\n5\n");
        let mut worker = Worker::new(|| {});
        worker.start(Work::Inspect(first.clone()));
        worker.cancel();
        worker.start(Work::Inspect(first));
        worker.start(Work::Inspect(last.clone()));
        let deadline = std::time::Instant::now() + std::time::Duration::from_secs(10);
        loop {
            if let Some(result) = worker.poll() {
                let Answer::Table(table) = result.unwrap() else {
                    panic!()
                };
                assert_eq!(table.path.file_name(), last.file_name());
                worker.start(Work::Preview(
                    table.clone(),
                    Box::new(Mapping::for_kind(&table, Kind::Line)),
                ));
                worker.cancel();
                break;
            }
            assert!(std::time::Instant::now() < deadline, "worker timed out");
            std::thread::sleep(std::time::Duration::from_millis(5));
        }
        assert!(worker.poll().is_none());
    }
}
