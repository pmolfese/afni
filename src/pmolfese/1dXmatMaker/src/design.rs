// Developed at the NIH
//
// Part of AFNI. See the repository LICENSE.txt for licensing terms and
// disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
// Third-party and adapted portions retain their existing attribution and license
// terms; see THIRD_PARTY.md.

//! Assemble baseline, nuisance, and event columns on the original timeline.
//!
//! Read `build_with_events` first; `build_design` is the nuisance-only adapter.
//! The builder validates identities, allocates one full matrix, fills run-local
//! polynomials, then processes indexed nuisance/events and appended ortvecs. The
//! retained-row map is stored separately and never changes these calculations.
//! This independently implemented algorithm follows the documented behavior
//! of AFNI Deconvolve.c:init_indep_var_matrix at the pinned M0 revision.
//!
//! Rust reading note: `Source` is a tagged union (enum). A match arm can only
//! access the fields belonging to that source's actual kind. `&sources` borrows
//! the collection during validation; later consuming it transfers ownership.
//! `&mut values[start..end]` grants a temporary exclusive view of one column,
//! without another allocation. The end of a Rust range is excluded.
//! No afni-core or afni-io call occurs: this module constructs numbers in Rust;
//! consumer integration happens through the serialized xmat, not native handles.

use crate::{
    basis::Basis,
    error::{Error, Result, reserve},
    manifest::{BasisIndex, EventIndex, EventMapping, StimulusIndex},
    model::{Acquisition, ColumnIndex, FullRow, RetainedRow, normalize_label},
    table::NumericTable,
    timing::Timing,
};
use std::collections::HashSet;

/// How a supplied nuisance source is labelled and ordered in the matrix.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NuisanceKind {
    /// A sampled -stim_file, optionally expanded into lag columns.
    StimFile,
    /// One -ortvec array, whose labels include [column] before #0.
    Ortvec,
}

/// An owned, validated nuisance table and its normalized source label.
#[derive(Debug, Clone)]
pub struct Nuisance {
    label: String,
    kind: NuisanceKind,
    table: NumericTable,
    baseline: bool,
    first_lag: usize,
}

impl Nuisance {
    /// Take ownership of a table and validate its label and source shape.
    ///
    /// StimFile requires one column; Ortvec accepts one or more. Expansion
    /// must fit AFNI's 60-byte pre-#0 label boundary without truncation.
    pub fn new(label: &str, kind: NuisanceKind, table: NumericTable) -> Result<Self> {
        let label = normalize_label(label)?;
        if kind == NuisanceKind::StimFile && table.ncols() != 1 {
            return Err(Error::input(format!(
                "'{label}': -stim_file requires one numeric input column"
            )));
        }
        if kind == NuisanceKind::Ortvec {
            normalize_label(&format!("{label}[{}]", table.ncols() - 1))?;
        }
        Ok(Self {
            label,
            kind,
            table,
            baseline: true,
            first_lag: 0,
        })
    }
    /// Expand a sampled stimulus with AFNI's historical global lag indexing.
    /// Input has nrows*nptr samples; output column k reads row*nptr-(minlag+k),
    /// with zeros before time zero. Unlike event responses, lags can cross a
    /// run boundary. No input padding/truncation or per-run reset is performed.
    pub fn sampled(
        label: &str,
        table: NumericTable,
        nrows: usize,
        nptr: usize,
        minlag: usize,
        maxlag: usize,
        baseline: bool,
    ) -> Result<Self> {
        if nptr == 0
            || minlag > maxlag
            || table.ncols() != 1
            || nrows.checked_mul(nptr) != Some(table.nrows())
        {
            return Err(Error::input(
                "invalid sampled stimulus shape, lags or sampling rate",
            ));
        }
        let count = maxlag
            .checked_sub(minlag)
            .and_then(|n| n.checked_add(1))
            .ok_or_else(|| Error::input("lag count overflow"))?;
        let mut values = Vec::new();
        reserve(
            &mut values,
            nrows
                .checked_mul(count)
                .ok_or_else(|| Error::input("lag expansion overflow"))?,
        )?;
        let input = table
            .column(0)
            .ok_or_else(|| Error::input("missing sampled column"))?;
        for lag in minlag..=maxlag {
            for row in 0..nrows {
                values.push((row * nptr).checked_sub(lag).map_or(0.0, |i| input[i]));
            }
        }
        Ok(Self {
            label: normalize_label(label)?,
            kind: NuisanceKind::StimFile,
            table: NumericTable::new(nrows, count, values)?,
            baseline,
            first_lag: minlag,
        })
    }
}

/// Ordinary events share basis columns; IM allocates components for each event.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EventMode {
    /// Sum responses from every event in this condition.
    Ordinary,
    /// Preserve each source event as its own set of basis columns.
    Individual,
    /// Sum supplied amplitudes into one response set.
    Am1,
    /// Mean response plus one centered response set per amplitude.
    Am2,
    /// Historical AMx spelling; same numerical model as AM2.
    Amx,
    /// FSL onset/duration/amplitude rows, equivalent to AM1.
    Fsl,
    /// FSL timing with amplitudes forced to one.
    Fs1,
}
impl EventMode {
    pub(crate) fn paired(self) -> bool {
        matches!(self, Self::Am2 | Self::Amx)
    }
    /// Exact AFNI option spelling used in basis metadata.
    pub fn option(self) -> &'static str {
        match self {
            Self::Ordinary => "-stim_times",
            Self::Individual => "-stim_times_IM",
            Self::Am1 => "-stim_times_AM1",
            Self::Am2 => "-stim_times_AM2",
            Self::Amx => "-stim_times_AMx",
            Self::Fsl => "-stim_times_FSL",
            Self::Fs1 => "-stim_times_FS1",
        }
    }
}

/// Owned response basis and validated onset list for one indexed stimulus.
#[derive(Debug)]
pub struct EventStimulus {
    label: String,
    original_label: String,
    basis: Basis,
    timing: Timing,
    mode: EventMode,
    centers: Vec<f32>,
    factors: Vec<f32>,
}
impl EventStimulus {
    /// Validate the label and coefficient count, then own the basis/timing.
    pub fn new(label: &str, basis: Basis, timing: Timing, mode: EventMode) -> Result<Self> {
        if basis.requires_duration() != timing.has_duration {
            return Err(Error::input(
                "basis and timing disagree on per-event duration",
            ));
        }
        if matches!(mode, EventMode::Ordinary | EventMode::Individual) && timing.modulators != 0 {
            return Err(Error::input(
                "ordinary/IM events do not accept amplitude modulators",
            ));
        }
        if mode == EventMode::Ordinary && timing.has_duration {
            return Err(Error::input("duration timing requires IM, AM or FSL mode"));
        }
        if mode == EventMode::Individual && timing.events().is_empty() {
            return Err(Error::input("IM requires at least one event"));
        }
        if matches!(mode, EventMode::Am1 | EventMode::Am2 | EventMode::Amx)
            && timing.modulators == 0
            && !timing.has_duration
        {
            return Err(Error::input(
                "AM timing requires an amplitude or duration field",
            ));
        }
        let mut centers = vec![0.0_f32; timing.modulators];
        // AFNI accumulates means in float across ALL events/runs, before
        // censoring. A supplied override is applied afterward by centered().
        for event in timing.events() {
            for (sum, &a) in centers.iter_mut().zip(&event.amplitudes) {
                *sum += a;
            }
        }
        if !timing.events().is_empty() {
            for sum in &mut centers {
                *sum /= timing.events().len() as f32;
            }
        }
        let factors = vec![1.0; basis.components()];
        let result = Self {
            label: normalize_label(label)?,
            original_label: label.to_owned(),
            basis,
            timing,
            mode,
            centers,
            factors,
        };
        result.ncols()?;
        Ok(result)
    }
    /// Override AM2 centering, or explicitly disable it. `None` entries retain
    /// the event-wide mean. The override vector must cover every modulator;
    /// raw mode and explicit overrides are mutually exclusive.
    pub fn centered(mut self, raw: bool, overrides: Option<&[Option<f32>]>) -> Result<Self> {
        if overrides.is_some() && (!self.mode.paired() || raw) {
            return Err(Error::input(
                "centering override requires AM2/AMx without -rawAM2",
            ));
        }
        if raw {
            self.centers.fill(0.0);
        }
        if let Some(values) = overrides {
            if values.len() != self.centers.len() {
                return Err(Error::input(
                    "centering override count must match amplitudes",
                ));
            }
            for (center, value) in self.centers.iter_mut().zip(values) {
                if let Some(v) = value {
                    if !v.is_finite() {
                        return Err(Error::input("nonfinite centering override"));
                    }
                    *center = *v;
                }
            }
        }
        Ok(self)
    }
    /// Normalize each basis component by AFNI's bounded peak-search grid.
    /// dm families use the reference parser's initial duration of ONE second;
    /// normalization is not recomputed for each event. EXPR random evaluation
    /// consumes its local stream here just as ordinary sampling would.
    pub fn normalized(mut self, target: f32, tr: f32) -> Result<Self> {
        if !target.is_finite() || target <= 0.0 || !tr.is_normal() || tr <= 0.0 {
            return Err(Error::input("invalid basis normalization/TR"));
        }
        let (bot, top) = self.basis.support();
        let mut dt = (top - bot) / 1999.0;
        let mut count = 1999;
        if dt < 0.01 * tr {
            dt = 0.01 * tr;
            count = 1 + ((top - bot) / dt) as usize;
        } else if dt > 0.1 * tr {
            dt = 0.1 * tr;
            count = 1 + ((top - bot) / dt) as usize;
        }
        if count > 1_000_000 {
            return Err(Error::input("normalization grid exceeds resource limit"));
        }
        for c in 0..self.factors.len() {
            let mut peak = 0.0_f32;
            for i in 0..=count {
                let t = bot + i as f32 * dt;
                let v = if self.basis.requires_duration() {
                    self.basis.sample_with_duration(t as f64, 1.0)
                } else {
                    self.basis.sample(c, t as f64)
                }
                .ok_or_else(|| Error::input("nonfinite normalization sample"))?;
                peak = peak.max((v as f32).abs());
            }
            if peak > 0.0 {
                self.factors[c] = target / peak;
            }
        }
        Ok(self)
    }
    fn ncols(&self) -> Result<usize> {
        let count = if self.mode == EventMode::Individual {
            self.timing.events().len()
        } else if self.mode.paired() {
            self.timing.modulators + 1
        } else {
            1
        };
        self.basis
            .components()
            .checked_mul(count)
            .filter(|&n| n <= i32::MAX as usize)
            .ok_or_else(|| Error::input("event coefficient count exceeds AFNI's index range"))
    }
}

/// One input source in indexed stimulus order; Ortvec arrays move to the end.
#[derive(Debug)]
pub enum Source {
    /// Supplied baseline vectors, including an ortvec array.
    Nuisance(Nuisance),
    /// Ordinary or IM timing and basis.
    Events(EventStimulus),
}
impl Source {
    fn label(&self) -> &str {
        match self {
            Self::Nuisance(n) => &n.label,
            Self::Events(e) => &e.label,
        }
    }
    fn ncols(&self) -> Result<usize> {
        match self {
            Self::Nuisance(n) => Ok(n.table.ncols()),
            Self::Events(e) => e.ncols(),
        }
    }
    fn is_ortvec(&self) -> bool {
        matches!(self, Self::Nuisance(n) if n.kind == NuisanceKind::Ortvec)
    }
}

/// AFNI metadata for one signal stimulus; nuisance-only sources are excluded.
#[derive(Debug)]
pub struct StimulusMetadata {
    /// Original zero-based stimulus identity, including intervening nuisance.
    pub index: StimulusIndex,
    /// Normalized stimulus label for StimLabels.
    pub label: String,
    /// Original label for BasisName (AFNI preserves whitespace here).
    pub basis_name: String,
    /// Original supported basis expression.
    pub formula: String,
    /// Ordinary or IM expansion.
    pub mode: EventMode,
    /// First matrix column, inclusive.
    pub first: ColumnIndex,
    /// Last matrix column, inclusive.
    pub last: ColumnIndex,
}

/// Meaning of a matrix column, including its AFNI stimulus group.
#[derive(Debug)]
pub struct Column {
    label: String,
    group: i32,
}

impl Column {
    /// Exact normalized label written in ColumnLabels.
    pub fn label(&self) -> &str {
        &self.label
    }
    /// AFNI group: -1 polynomial, 0 supplied nuisance, positive original stimulus.
    pub fn group(&self) -> i32 {
        self.group
    }
}

/// Immutable generated matrix plus acquisition and censor identities.
#[derive(Debug)]
pub struct Design {
    acquisition: Acquisition,
    columns: Vec<Column>,
    values: Vec<f64>,
    retained: Vec<FullRow>,
    stimuli: Vec<StimulusMetadata>,
    events: Vec<EventMapping>,
    basis_nstim: usize,
    // Symbolic GLT map: (source name, first physical column, minimum user
    // subscript, coefficient count). The minimum is zero for events/Ort but
    // equals minlag for sampled sources. Neither number is a stimulus ID.
    pub(crate) ranges: Vec<(String, usize, usize, usize)>,
    pub(crate) contrasts: Vec<crate::contrast::Contrast>,
}

impl Design {
    /// Signal ranges and basis metadata in original stimulus order.
    pub fn stimuli(&self) -> &[StimulusMetadata] {
        &self.stimuli
    }
    /// Event identities/column relations in stable event then component order.
    pub fn events(&self) -> &[EventMapping] {
        &self.events
    }
    /// AFNI's extended source count, including individual ortvec columns.
    pub fn basis_nstim(&self) -> usize {
        self.basis_nstim
    }

    /// Borrow acquisition metadata; no caller can alter the validated layout.
    pub fn acquisition(&self) -> &Acquisition {
        &self.acquisition
    }
    /// Borrow labels/groups in the same order as numeric columns.
    pub fn columns(&self) -> &[Column] {
        &self.columns
    }
    /// Map retained positions to full rows. Never use this as a run-start map.
    pub fn retained_rows(&self) -> &[FullRow] {
        &self.retained
    }
    /// Resolve a retained-row identity, or None when out of range.
    pub fn full_row(&self, row: RetainedRow) -> Option<FullRow> {
        self.retained.get(row.0).copied()
    }
    /// Borrow a full-timeline column; None means the requested column is absent.
    pub fn column(&self, column: ColumnIndex) -> Option<&[f64]> {
        if column.0 >= self.columns.len() {
            return None;
        }
        let start = column.0 * self.acquisition.nrows();
        Some(&self.values[start..start + self.acquisition.nrows()])
    }
    /// Read a full-timeline value with checked row and column identities.
    pub fn value(&self, row: FullRow, column: ColumnIndex) -> Option<f64> {
        self.column(column)?.get(row.0).copied()
    }
}

/// Create a full design and an ordered retained-row map.
///
/// Nuisance input order is significant: StimFile sources retain their relative
/// order (CLI stimulus-number order), followed by Ortvec arrays in their
/// relative order. The tables are borrowed during copying, then dropped when
/// this function returns; the resulting Design owns all required numeric data.
/// `keep=None` retains every row. Reject empty retained sets, duplicate source
/// or final labels, nonfinite generated values, and incompatible dimensions.
/// Zero/duplicate numeric columns and ncolumns>=nrows are deliberately legal.
pub fn build_design(
    acquisition: Acquisition,
    nuisance: Vec<Nuisance>,
    keep: Option<&[bool]>,
) -> Result<Design> {
    build_with_events(
        acquisition,
        nuisance.into_iter().map(Source::Nuisance).collect(),
        keep,
    )
}

/// Build baseline, nuisance, and event columns in AFNI source order.
///
/// Sources are owned; numerical buffers are copied once into the full design.
/// Indexed sources preserve order, and ortvec arrays follow them. Timing must
/// have been validated against this acquisition. Event columns are never
/// demeaned, and all generation precedes censoring. No support/rank filtering
/// changes column IDs, even for duplicate events or fully censored responses.
pub fn build_with_events(
    acquisition: Acquisition,
    sources: Vec<Source>,
    keep: Option<&[bool]>,
) -> Result<Design> {
    let nrows = acquisition.nrows();
    // Degree p needs p+1 coefficients: P0 through Pp. With p=-1, this
    // expression is exactly zero. Widen before addition to avoid i32 overflow.
    let per_run = (i64::from(acquisition.polort()) + 1) as usize;
    let polynomial_count = per_run
        .checked_mul(acquisition.run_lengths().len())
        .ok_or_else(|| Error::input("polynomial column count overflow"))?;
    let mut ncols = polynomial_count;
    let mut names = HashSet::new();
    for source in &sources {
        match source {
            Source::Nuisance(n) if n.table.nrows() != nrows => {
                return Err(Error::input(format!(
                    "'{}': {} rows; expected {nrows} full rows",
                    n.label,
                    n.table.nrows()
                )));
            }
            Source::Events(e) => e.timing.validate_for(&acquisition)?,
            _ => {}
        }
        if !names.insert(source.label()) {
            return Err(Error::input(format!(
                "duplicate normalized source label '{}'",
                source.label()
            )));
        }
        ncols = ncols
            .checked_add(source.ncols()?)
            .ok_or_else(|| Error::input("column count overflow"))?;
    }
    if ncols == 0 || ncols > i32::MAX as usize {
        return Err(Error::input("design must have 1..2147483647 columns"));
    }
    let elements = nrows
        .checked_mul(ncols)
        .filter(|&n| n <= isize::MAX as usize / size_of::<f64>())
        .ok_or_else(|| Error::input("design dimensions exceed addressable f64 storage"))?;
    if let Some(mask) = keep {
        if mask.len() != nrows {
            return Err(Error::input("censor length must equal full row count"));
        }
        if !mask.iter().any(|&v| v) {
            return Err(Error::input("censor removes every row"));
        }
    }

    let mut retained = Vec::new();
    reserve(
        &mut retained,
        keep.map_or(nrows, |m| m.iter().filter(|&&v| v).count()),
    )?;
    for row in 0..nrows {
        if keep.is_none_or(|mask| mask[row]) {
            retained.push(FullRow(row));
        }
    }
    let mut columns = Vec::new();
    reserve(&mut columns, ncols)?;
    let mut values = Vec::new();
    reserve(&mut values, elements)?;
    values.resize(elements, 0.0);

    // Every column is contiguous: values[column*nrows + full_row]. Initial
    // zeros outside each run are structural zeros, never altered by demeaning.
    for (run, &length) in acquisition.run_lengths().iter().enumerate() {
        let first = acquisition.run_starts()[run].0;
        for degree in 0..per_run {
            let column = run * per_run + degree;
            let begin = column * nrows + first;
            let target = &mut values[begin..begin + length];
            for (local_row, value) in target.iter_mut().enumerate() {
                // A one-sample run is allowed only for P0=1 (or no baseline).
                let x = if length == 1 {
                    0.0
                } else {
                    2.0 * local_row as f64 / (length - 1) as f64 - 1.0
                };
                *value = if acquisition.power_baseline() {
                    (local_row as f64).powi(degree as i32)
                } else {
                    legendre(x, degree)
                };
            }
            if degree > 0 && acquisition.demean() {
                demean(target)?;
            }
            columns.push(Column {
                label: format!("Run#{}Pol#{degree}", run + 1),
                group: -1,
            });
        }
    }

    let mut stimuli = Vec::new();
    let mut events = Vec::new();
    let mut event_count = 0usize;
    let mut basis_nstim = 0usize;
    let mut ranges = Vec::new();
    if polynomial_count > 0 {
        ranges.push(("Ort".to_owned(), 0, 0, polynomial_count));
    }
    // Two passes preserve indexed order while putting ortvec arrays last.
    for ortvec_pass in [false, true] {
        for source in sources.iter().filter(|s| s.is_ortvec() == ortvec_pass) {
            match source {
                Source::Nuisance(source) => {
                    let first = columns.len();
                    // Only indexed sources have unambiguous symbolic names.
                    // AFNI expands ortvec into labels containing brackets, which
                    // collide with GLT selector syntax; use numeric GLTs for them.
                    if source.kind == NuisanceKind::StimFile {
                        ranges.push((
                            source.label.clone(),
                            first,
                            source.first_lag,
                            source.table.ncols(),
                        ));
                    }
                    for local_column in 0..source.table.ncols() {
                        let label = match source.kind {
                            NuisanceKind::StimFile => format!("{}#{local_column}", source.label),
                            NuisanceKind::Ortvec => format!("{}[{local_column}]#0", source.label),
                        };
                        let start = columns.len() * nrows;
                        let target = &mut values[start..start + nrows];
                        let supplied = source
                            .table
                            .column(local_column)
                            .ok_or_else(|| Error::resource("internal nuisance column mismatch"))?;
                        target.copy_from_slice(supplied);
                        if source.baseline && acquisition.demean() {
                            demean(target)?;
                        }
                        if label.len() > 63 {
                            return Err(Error::input(
                                "expanded sampled column label exceeds AFNI's 63-byte limit",
                            ));
                        }
                        columns.push(Column {
                            label,
                            group: if source.baseline {
                                0
                            } else {
                                (basis_nstim + 1) as i32
                            },
                        });
                    }
                    if !source.baseline {
                        stimuli.push(StimulusMetadata {
                            index: StimulusIndex(basis_nstim),
                            label: source.label.clone(),
                            basis_name: source.label.clone(),
                            formula: String::new(),
                            mode: EventMode::Ordinary,
                            first: ColumnIndex(first),
                            last: ColumnIndex(columns.len() - 1),
                        });
                    }
                    basis_nstim += if source.kind == NuisanceKind::Ortvec {
                        source.table.ncols()
                    } else {
                        1
                    };
                }
                Source::Events(source) => {
                    let first = columns.len();
                    let count = source.ncols()?;
                    ranges.push((source.label.clone(), first, 0, count));
                    for column in 0..count {
                        // AFNI's 64-byte column buffer must hold the decimal
                        // coefficient suffix too, not only its 60-byte prefix.
                        let label = format!("{}#{column}", source.label);
                        if label.len() > 63 {
                            return Err(Error::input(
                                "expanded event column label exceeds AFNI's 63-byte limit",
                            ));
                        }
                        columns.push(Column {
                            label,
                            group: (basis_nstim + 1) as i32,
                        });
                    }
                    let index = StimulusIndex(basis_nstim);
                    fill_events(
                        source,
                        &acquisition,
                        first,
                        index,
                        event_count,
                        &mut values,
                        &mut events,
                    )?;
                    event_count = event_count
                        .checked_add(source.timing.events().len())
                        .ok_or_else(|| Error::input("event identity overflow"))?;
                    stimuli.push(StimulusMetadata {
                        index,
                        label: source.label.clone(),
                        basis_name: source.original_label.clone(),
                        formula: source.basis.formula().to_owned(),
                        mode: source.mode,
                        first: ColumnIndex(first),
                        last: ColumnIndex(first + count - 1),
                    });
                    basis_nstim += 1;
                }
            }
        }
    }
    // Stable IDs follow run/source order, while global IM coefficient order
    // remains source-token order. Sorting relations never reorders the matrix.
    events.sort_unstable_by_key(|m| (m.event.0, m.basis.0));
    let mut final_labels = HashSet::new();
    for column in &columns {
        if !final_labels.insert(column.label()) {
            return Err(Error::input(format!(
                "duplicate expanded column label '{}'",
                column.label()
            )));
        }
    }
    if values.iter().any(|v| !v.is_finite()) {
        return Err(Error::input("generation produced a nonfinite matrix value"));
    }
    Ok(Design {
        acquisition,
        columns,
        values,
        retained,
        stimuli,
        events,
        basis_nstim,
        ranges,
        contrasts: Vec::new(),
    })
}

/// Sample each event only inside its original run, preserving AFNI's f32
/// timeline arithmetic and event-wise float accumulation. The final storage is
/// f64. myceil/myfloor's 0.00005 index tolerance and the separate 0.001*TR
/// support check mirror 3dDeconvolve.c's sampling loop at the pinned revision.
fn fill_events(
    source: &EventStimulus,
    acquisition: &Acquisition,
    first: usize,
    stimulus: StimulusIndex,
    first_event: usize,
    values: &mut [f64],
    mappings: &mut Vec<EventMapping>,
) -> Result<()> {
    let nrows = acquisition.nrows();
    let nf = source.basis.components();
    let tr = acquisition.tr() as f32;
    let (bot, top) = source.basis.support();
    let input = source.timing.events();
    let mut order = Vec::new();
    reserve(&mut order, input.len())?;
    order.extend(0..input.len());
    // Sort INDICES, not events. This gives stable run-based IDs while the
    // original vector still determines IM coefficient order. The key includes
    // the original index, so distinct events never have equal sorting keys.
    order.sort_unstable_by_key(|&i| (input[i].run.0, i));
    let mut ids = Vec::new();
    reserve(&mut ids, input.len())?;
    ids.resize(input.len(), 0);
    for (id, &original) in order.iter().enumerate() {
        ids[original] = first_event
            .checked_add(id)
            .ok_or_else(|| Error::input("event identity overflow"))?;
    }
    let sets = if source.mode.paired() {
        source.timing.modulators + 1
    } else {
        1
    };
    reserve(
        mappings,
        input
            .len()
            .checked_mul(nf)
            .and_then(|n| n.checked_mul(sets))
            .ok_or_else(|| Error::input("event manifest size overflow"))?,
    )?;
    for (ordinal, event) in input.iter().enumerate() {
        let base = first
            + if source.mode == EventMode::Individual {
                ordinal * nf
            } else {
                0
            };
        let run = acquisition
            .run_range(event.run)
            .ok_or_else(|| Error::input("event run is out of bounds"))?;
        // Clamp floating bounds BEFORE conversion to usize. This avoids
        // negative-index casts and needless work outside a response's support.
        let lower = ((event.timedex + bot / tr) as f64 - 0.00005)
            .ceil()
            .max(run.start as f64);
        let upper = ((event.timedex + top / tr) as f64 + 0.00005)
            .floor()
            .min((run.end - 1) as f64);
        let sets = if source.mode.paired() {
            source.timing.modulators + 1
        } else {
            1
        };
        // A modulation set is a complete copy of the response basis. AM2 set
        // zero is the mean response; sets 1.. use amplitude minus its center.
        // Flattening set*nf+component keeps metadata and manifest ordering in
        // the same coefficient coordinates used by symbolic contrasts.
        for set in 0..sets {
            let amplitude = if source.mode.paired() && set > 0 {
                event.amplitudes[set - 1] - source.centers[set - 1]
            } else if matches!(
                source.mode,
                EventMode::Am1 | EventMode::Fsl | EventMode::Fs1
            ) && !event.amplitudes.is_empty()
            {
                event.amplitudes.iter().fold(0.0_f32, |sum, &a| sum + a)
            } else {
                1.0
            };
            for component in 0..nf {
                mappings.push(EventMapping {
                    event: EventIndex(ids[ordinal]),
                    stimulus,
                    condition: stimulus,
                    run: event.run,
                    onset: event.onset,
                    duration: event.duration.map_or(source.basis.duration(), f64::from),
                    basis: BasisIndex(set * nf + component),
                    column: ColumnIndex(base + set * nf + component),
                });
                if lower > upper {
                    continue;
                }
                for row in lower as usize..=upper as usize {
                    let seconds = tr * (row as f32 - event.timedex);
                    if seconds + 0.001 * tr < bot || seconds - 0.001 * tr > top {
                        continue;
                    }
                    // Zero amplitude does not evaluate the waveform. This also
                    // keeps the documented local EXPR stream from consuming a draw.
                    if amplitude == 0.0 {
                        continue;
                    }
                    let sample = if let Some(duration)=event.duration {
                    source.basis.sample_with_duration(seconds as f64,duration as f64)
                } else {source.basis.sample(component,seconds as f64)}
                    .ok_or_else(|| Error::input(format!(
                        "stimulus {} ('{}'), run {}, basis {}, full row {row}: nonfinite response",
                        stimulus.0+1, source.label, event.run.0+1, component
                    )))?;
                    let location = (base + set * nf + component) * nrows + row;
                    values[location] = ((values[location] as f32)
                        + amplitude * ((sample as f32) * source.factors[component]))
                        as f64;
                }
            }
        }
    }
    Ok(())
}

/// Evaluate Legendre P_n(x) using the standard three-term recurrence.
/// P0=1, P1=x, Pn=((2n-1)x P(n-1)-(n-1)P(n-2))/n.
/// Keeping two prior values avoids recursion and does not require a solver.
fn legendre(x: f64, degree: usize) -> f64 {
    if degree == 0 {
        return 1.0;
    }
    let mut previous = 1.0;
    let mut current = x;
    for n in 2..=degree {
        let next = ((2 * n - 1) as f64 * x * current - (n - 1) as f64 * previous) / n as f64;
        previous = current;
        current = next;
    }
    current
}

/// Subtract the full vector's mean, using f64 rather than AFNI's float sum.
///
/// Scale before summing to avoid overflowing a finite mean. Compensation
/// reduces loss of small terms. The frozen oracle uses float inputs/sums;
/// small roundoff differences are checked with documented parity tolerances.
/// Extreme inputs whose centered differences overflow are rejected explicitly.
fn demean(values: &mut [f64]) -> Result<()> {
    let count = values.len() as f64;
    let mut mean = 0.0;
    let mut correction = 0.0;
    for &value in values.iter() {
        let term = value / count - correction;
        let next = mean + term;
        correction = (next - mean) - term;
        mean = next;
    }
    for value in values {
        *value -= mean;
        if !value.is_finite() {
            return Err(Error::input("demeaning overflows finite f64 values"));
        }
    }
    Ok(())
}
