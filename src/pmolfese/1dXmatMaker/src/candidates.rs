// Developed at the NIH
//
// Part of AFNI. See the repository LICENSE.txt for licensing terms and
// disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
// Third-party and adapted portions retain their existing attribution and license
// terms; see THIRD_PARTY.md.

//! Generate aligned trial/ONOFF candidates and Picard's versioned identities.
//!
//! Preparation owns a small collection of HRF curves. Publication holds only
//! one full design at a time and stages its artifacts before moving on; it
//! never allocates candidates x rows x columns. All paths are preflighted as
//! one rollback operation. Read load_trials, load_library, then generate.
//!
//! Picard performs selection/fitting. This producer only emits ordinary AFNI
//! matrices plus explicit stable IDs. No afni-core, afni-io, fitter or shell
//! command is called. Source waveform provenance never stands in for an
//! independently validated library: the built-in GLMsingle data are pinned.
use crate::{
    basis::Basis,
    design::{Design, EventMode, EventStimulus, Nuisance, Source, build_with_events},
    error::{Error, Result},
    input,
    model::Acquisition,
    publication,
    sampled::{self, Curve},
    timing::{Timing, TimingMode},
    xmat::{RowSelection, write_matrix},
};
use std::{
    collections::HashSet,
    path::{Path, PathBuf},
};

/// One explicit trial. IDs are nonnegative; runs are zero-based acquisition IDs.
#[derive(Debug, Clone)]
pub struct Trial {
    /// Stable event/trial ID, independent of file row and matrix column.
    pub id: usize,
    /// Repeated-condition ID used by Picard tuning and downstream adapters.
    pub condition: usize,
    /// Original run ID, unaffected by censoring.
    pub run: usize,
    /// Seconds relative to the run start.
    pub onset: f64,
    /// Event duration in seconds; all candidates use this same duration.
    pub duration: f64,
    /// User item identifier or '-' for absent; never guessed from a label.
    pub item: String,
}
/// One prepared HRF with the provenance required by Picard's HRF policy.
#[derive(Debug)]
pub struct Candidate {
    /// Stable candidate ID, exactly representable in an AFNI float ID map.
    pub id: usize,
    /// Prepared waveform, already convolved for the declared duration.
    pub curve: Curve,
    /// Original source sampling interval in seconds, before resampling.
    pub source_dt: f64,
    /// Library or user source identifier.
    pub source: String,
    /// Pinned source revision or user library version.
    pub version: String,
    /// Common units/normalization convention across the collection.
    pub normalization: String,
}

/// Read the explicit TSV event contract; preserve file order for IM columns.
/// Identical onsets are legal distinct trials. Duplicate IDs and malformed
/// missing values are errors. Duration varies only across separate campaigns.
pub fn load_trials(path: &Path, a: &Acquisition) -> Result<Vec<Trial>> {
    let rows = tsv(path, "trial\tcondition\trun\tonset\tduration\titem")?;
    let mut trials = Vec::new();
    let mut ids = HashSet::new();
    for row in rows {
        if row.len() != 6 {
            return Err(Error::input("events need six TSV fields"));
        }
        let id = identifier(&row[0])?;
        let condition = identifier(&row[1])?;
        let run = identifier(&row[2])?;
        let onset = number(&row[3])?;
        let duration = number(&row[4])?;
        if !ids.insert(id)
            || run >= a.run_lengths().len()
            || onset < 0.0
            || onset > (a.run_lengths()[run] - 1) as f64 * a.tr()
            || !(0.0..=999.0).contains(&duration)
        {
            return Err(Error::input(
                "invalid/duplicate trial identity, onset or duration",
            ));
        }
        trials.push(Trial {
            id,
            condition,
            run,
            onset,
            duration,
            item: row[5].clone(),
        });
    }
    if trials.is_empty() || trials.iter().any(|t| t.duration != trials[0].duration) {
        return Err(Error::input(
            "candidate events need a nonempty common-duration trial table",
        ));
    }
    Ok(trials)
}
/// Build the embedded GLMsingle collection or prepare a versioned custom TSV.
/// Custom kind=impulse applies one continuous box integral; kind=convolved
/// requires matching duration and applies NO convolution. This explicit tag
/// prevents a plausible-looking but incorrect double-duration response.
pub fn load_library(
    spec: &str,
    duration: f64,
    tr: f64,
    normalization: &str,
) -> Result<Vec<Candidate>> {
    if !matches!(normalization, "peak" | "none" | "library") {
        return Err(Error::input("normalization must be peak, none or library"));
    }
    if spec == "glmsingle" {
        if normalization == "none" {
            return Err(Error::input(
                "GLMsingle supports peak or library normalization",
            ));
        }
        return sampled::glmsingle(duration, tr, normalization == "peak")?
            .into_iter()
            .enumerate()
            .map(|(id, curve)| {
                Ok(Candidate {
                    id,
                    curve,
                    source_dt: 0.1,
                    source: "GLMsingle".into(),
                    version: sampled::GLMSINGLE_REVISION.into(),
                    normalization: if normalization == "peak" {
                        "glmsingle-fit-peak1"
                    } else {
                        "glmsingle-library-globalpeak1"
                    }
                    .into(),
                })
            })
            .collect();
    }
    if normalization == "library" {
        return Err(Error::input(
            "library normalization is specific to GLMsingle",
        ));
    }
    let path = Path::new(spec);
    let parent = path.parent().unwrap_or(Path::new("."));
    let rows = tsv(path, "id\twaveform\tdt\tkind\tduration\tsource\tversion")?;
    let mut result = Vec::new();
    let mut ids = HashSet::new();
    for row in rows {
        if row.len() != 7 {
            return Err(Error::input("library needs seven TSV fields"));
        }
        let id = identifier(&row[0])?;
        if !ids.insert(id) {
            return Err(Error::input("duplicate candidate ID"));
        }
        let tokens = input::rows(&parent.join(&row[1]), false)?;
        if tokens.iter().any(|r| r.len() != 1) {
            return Err(Error::input("custom waveform requires one numeric column"));
        }
        let values = tokens
            .iter()
            .map(|r| number(&r[0]))
            .collect::<Result<Vec<_>>>()?;
        let dt = number(&row[2])?;
        let declared = number(&row[4])?;
        let curve = Curve::new(values, dt, declared)?;
        let mut curve = match row[3].as_str() {
            "impulse" if declared == 0.0 => curve.convolve_box(duration)?,
            "convolved" if declared == duration => curve,
            _ => {
                return Err(Error::input(
                    "impulse needs declared duration 0; convolved duration must match events",
                ));
            }
        };
        if normalization == "peak" {
            curve = curve.peak_normalized()?;
        }
        result.push(Candidate {
            id,
            curve,
            source_dt: dt,
            source: row[5].clone(),
            version: row[6].clone(),
            normalization: normalization.into(),
        });
    }
    Ok(result)
}
/// Parse finite numeric metadata without permitting partial text conversion.
fn number(s: &str) -> Result<f64> {
    s.parse::<f64>()
        .ok()
        .filter(|v| v.is_finite())
        .ok_or_else(|| Error::input("invalid finite number in candidate input"))
}
/// Bound IDs to the integer range exactly represented by AFNI float maps.
fn identifier(s: &str) -> Result<usize> {
    s.parse::<usize>()
        .ok()
        .filter(|&n| n <= 16_777_215)
        .ok_or_else(|| Error::input("candidate/trial IDs must be in 0..16777215"))
}
/// Strict TAB grammar. Fields retain spaces; tabs/newlines cannot be data.
fn tsv(path: &Path, header: &str) -> Result<Vec<Vec<String>>> {
    let text = std::fs::read_to_string(path).map_err(|e| Error::io("cannot read", path, e))?;
    let mut lines = text
        .lines()
        .filter(|s| !s.trim().is_empty() && !s.starts_with('#'));
    if lines.next() != Some(header) {
        return Err(Error::input(format!(
            "{}: expected TSV header {header}",
            path.display()
        )));
    }
    lines
        .map(|line| {
            if line.len() > 65536 || line.chars().any(|c| c.is_control() && c != '\t') {
                return Err(Error::input("invalid TSV control/line length"));
            }
            let fields = line.split('\t').map(str::to_owned).collect::<Vec<_>>();
            if fields.iter().any(String::is_empty) {
                return Err(Error::input(
                    "empty TSV field; use explicit '-' where allowed",
                ));
            }
            Ok(fields)
        })
        .collect()
}

/// Collection inputs, moved into generation so output construction has one owner.
/// Global nuisance columns are allowed, but Picard repeated-condition tuning
/// additionally requires run-local nuisance; the consumer enforces that policy.
pub struct Request {
    /// Validated acquisition layout.
    pub acquisition: Acquisition,
    /// Explicit trial identities in source order.
    pub trials: Vec<Trial>,
    /// Two or more prepared and aligned HRF candidates.
    pub library: Vec<Candidate>,
    /// Supplied nuisance tables appended after the event source.
    pub nuisance: Vec<Nuisance>,
    /// Original timeline keep mask; None keeps all rows.
    pub keep: Option<Vec<bool>>,
    /// Output basename; its parent directory must already exist.
    pub prefix: PathBuf,
}
/// Staged artifact type. Carrying the type beside its path avoids fragile
/// arithmetic offsets when adding an optional output in the future.
enum Artifact {
    Matrix(usize, bool, bool),
    Columns(bool),
    Trials,
    Events,
    Waveform(usize),
    Catalog(bool),
}

/// Generate trial and condition-collapsed ONOFF designs as ordinary xmat files.
/// One cache contains one candidate/mode design; all stages remain unpublished
/// until every design, identity and waveform has been written successfully.
pub fn generate(request: Request, command: &str) -> Result<()> {
    let Request {
        acquisition,
        trials,
        library,
        nuisance,
        keep,
        prefix,
    } = request;
    if !(2..=4096).contains(&library.len()) || trials.is_empty() {
        return Err(Error::input(
            "need 2..4096 candidates and at least one trial",
        ));
    }
    // Request fields are public library inputs; do not assume that callers
    // used load_trials. Validate again before indexing a run or writing a TSV.
    let mut trial_ids = HashSet::new();
    for t in &trials {
        if !trial_ids.insert(t.id)
            || t.id > 16_777_215
            || t.condition > 16_777_215
            || t.run >= acquisition.run_lengths().len()
            || !t.onset.is_finite()
            || t.onset < 0.0
            || t.onset > (acquisition.run_lengths()[t.run] - 1) as f64 * acquisition.tr()
            || !t.duration.is_finite()
            || !(0.0..=999.0).contains(&t.duration)
            || t.duration != trials[0].duration
            || t.item.is_empty()
            || t.item.chars().any(char::is_control)
        {
            return Err(Error::input(
                "invalid candidate trial identity, coordinates or item",
            ));
        }
    }
    let mut ids = HashSet::new();
    for candidate in &library {
        if !ids.insert(candidate.id)
            || candidate.id > 16_777_215
            || !candidate.source_dt.is_finite()
            || candidate.source_dt <= 0.0
            || candidate.curve.duration() != trials[0].duration
            || candidate.normalization != library[0].normalization
        {
            return Err(Error::input(
                "candidate IDs, duration or normalization disagree",
            ));
        }
        for text in [
            &candidate.source,
            &candidate.version,
            &candidate.normalization,
        ] {
            if text.is_empty() || text == "-" || text.chars().any(char::is_control) {
                return Err(Error::input(
                    "candidate provenance must be explicit single-line text",
                ));
            }
        }
    }
    let basename = prefix
        .file_name()
        .and_then(|s| s.to_str())
        .ok_or_else(|| Error::input("candidate prefix must be a UTF-8 basename"))?;
    if basename.chars().any(char::is_control) {
        return Err(Error::input("invalid candidate prefix"));
    }
    let parent = prefix
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let mut artifacts: Vec<(PathBuf, Artifact)> = Vec::new();
    let mut add =
        |suffix: String, kind| artifacts.push((parent.join(format!("{basename}{suffix}")), kind));
    for (index, c) in library.iter().enumerate() {
        for onoff in [false, true] {
            for full in [false, true] {
                add(
                    format!(
                        ".id{}.{}{}.xmat.1D",
                        c.id,
                        if onoff { "onoff" } else { "trials" },
                        if full { ".full" } else { "" }
                    ),
                    Artifact::Matrix(index, onoff, full),
                );
            }
        }
        add(format!(".id{}.hrf.1D", c.id), Artifact::Waveform(index));
    }
    add(".trials.columns.tsv".into(), Artifact::Columns(false));
    add(".onoff.columns.tsv".into(), Artifact::Columns(true));
    add(".trial_info.1D".into(), Artifact::Trials);
    add(".events.tsv".into(), Artifact::Events);
    add(".candidates.tsv".into(), Artifact::Catalog(false));
    add(".onoff.candidates.tsv".into(), Artifact::Catalog(true));
    let paths = artifacts
        .iter()
        .map(|(p, _)| p.as_path())
        .collect::<Vec<_>>();
    // All input events stay in source order. Global timing lets IM preserve
    // that order even when the table interleaves runs; explicit run IDs remain
    // in the identity table instead of being inferred later from nonzero rows.
    let timing = trials
        .iter()
        .map(|t| {
            (t.onset + acquisition.run_starts()[t.run].0 as f64 * acquisition.tr()).to_string()
        })
        .collect::<Vec<_>>()
        .join(" ");
    let build = |index: usize, onoff: bool| -> Result<Design> {
        let candidate = &library[index];
        let basis = Basis::sampled(
            candidate.curve.clone(),
            format!(
                "SAMPLED(id={},dt={},duration={})",
                candidate.id,
                candidate.curve.dt(),
                candidate.curve.duration()
            ),
        )?;
        let parsed = Timing::parse(
            &timing,
            "candidate events",
            TimingMode::Global,
            &acquisition,
        )?;
        let events = EventStimulus::new(
            if onoff { "ONOFF" } else { "trial" },
            basis,
            parsed,
            if onoff {
                EventMode::Ordinary
            } else {
                EventMode::Individual
            },
        )?;
        let mut sources = vec![Source::Events(events)];
        sources.extend(nuisance.iter().cloned().map(Source::Nuisance));
        let design = build_with_events(acquisition.clone(), sources, keep.as_deref())?;
        if design.stimuli().len() != 1 {
            return Err(Error::input(
                "candidate nuisance sources must be baseline columns",
            ));
        }
        Ok(design)
    };
    let mut cache: Option<(usize, bool, Design)> = None;
    publication::publish(&paths, |index, writer| {
        let wanted = match artifacts[index].1 {
            Artifact::Matrix(c, o, _) => Some((c, o)),
            Artifact::Columns(o) => Some((0, o)),
            Artifact::Trials => Some((0, false)),
            _ => None,
        };
        if let Some((c, o)) = wanted
            && cache
                .as_ref()
                .is_none_or(|(previous, mode, _)| (*previous, *mode) != (c, o))
        {
            cache = Some((c, o, build(c, o)?));
        }
        if let Artifact::Matrix(_, _, full) = artifacts[index].1 {
            return write_matrix(
                &cache
                    .as_ref()
                    .ok_or_else(|| Error::resource("missing candidate cache"))?
                    .2,
                if full {
                    RowSelection::Full
                } else {
                    RowSelection::Retained
                },
                command,
                writer,
            );
        }
        // Non-matrix output errors pass through one I/O conversion. The
        // borrowed writer belongs to publication; this closure never closes it.
        let result = (|| -> std::io::Result<()> {
            match artifacts[index].1 {
                Artifact::Waveform(c) => {
                    writeln!(
                        writer,
                        "# prepared HRF: dt={} duration={} normalization={}",
                        library[c].curve.dt(),
                        library[c].curve.duration(),
                        library[c].normalization
                    )?;
                    for v in library[c].curve.values() {
                        writeln!(writer, "{v:.16e}")?;
                    }
                }
                Artifact::Catalog(onoff) => {
                    writeln!(
                        writer,
                        "# 3dPicard candidates 1\nid\tmatrix\tcolumns\tsource\tversion\tsampling\tduration\tnormalization"
                    )?;
                    let mode = if onoff { "onoff" } else { "trials" };
                    for c in &library {
                        writeln!(
                            writer,
                            "{}\t{basename}.id{}.{}.xmat.1D\t{basename}.{mode}.columns.tsv\t{}\t{}\t{}\tprepared_seconds={}\t{}",
                            c.id,
                            c.id,
                            mode,
                            c.source,
                            c.version,
                            c.source_dt,
                            c.curve.duration(),
                            c.normalization
                        )?;
                    }
                }
                Artifact::Events => {
                    writeln!(writer, "trial\tcondition\trun\tonset\tduration\titem")?;
                    for t in &trials {
                        writeln!(
                            writer,
                            "{}\t{}\t{}\t{}\t{}\t{}",
                            t.id, t.condition, t.run, t.onset, t.duration, t.item
                        )?;
                    }
                }
                Artifact::Columns(_) | Artifact::Trials => {
                    let onoff = matches!(artifacts[index].1, Artifact::Columns(true));
                    let design = &cache
                        .as_ref()
                        .ok_or_else(|| std::io::Error::other("missing identity cache"))?
                        .2;
                    let first = design.stimuli()[0].first.0;
                    if matches!(artifacts[index].1, Artifact::Trials) {
                        writeln!(writer, "# column trial condition run")?;
                        for (i, t) in trials.iter().enumerate() {
                            writeln!(writer, "{} {} {} {}", first + i, t.id, t.condition, t.run)?;
                        }
                    } else {
                        writeln!(
                            writer,
                            "column\tlabel\trole\ttrial\tcondition\trun\tbasis\titem"
                        )?;
                        for (column, meta) in design.columns().iter().enumerate() {
                            if meta.group() > 0 && !onoff {
                                let t = &trials[column - first];
                                writeln!(
                                    writer,
                                    "{column}\t{}\tpenalized\t{}\t{}\t{}\t0\t{}",
                                    meta.label(),
                                    t.id,
                                    t.condition,
                                    t.run,
                                    t.item
                                )?;
                            } else {
                                writeln!(
                                    writer,
                                    "{column}\t{}\t{}\t-\t-\t-\t0\t-",
                                    meta.label(),
                                    if meta.group() > 0 {
                                        "penalized"
                                    } else {
                                        "nuisance"
                                    }
                                )?;
                            }
                        }
                    }
                }
                Artifact::Matrix(..) => unreachable!("matrix handled above"),
            }
            Ok(())
        })();
        result.map_err(|e| Error::resource(format!("cannot write candidate artifact: {e}")))
    })
}
