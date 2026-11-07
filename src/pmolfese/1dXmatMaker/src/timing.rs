// Developed at the NIH
//
// Part of AFNI. See the repository LICENSE.txt for licensing terms and
// disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
// Third-party and adapted portions retain their existing attribution and license
// terms; see THIRD_PARTY.md.

//! Explicit local/global timing and event identities, before censoring.
//!
//! Local noncomment rows correspond to runs; '*' alone represents an empty
//! run. Global tokens keep their input order for AFNI IM column compatibility.
//! No sorting, amplitude repair, or discarded out-of-range event
//! is allowed. Width-based inference requires explicit TimingMode::Guess. A Timing object remembers its acquisition layout so it cannot
//! accidentally be reused with a different run/TR model.
//!
//! A timing event has two coordinates on purpose: onset preserves the input
//! seconds for the manifest, while timedex preserves AFNI's float sampling
//! arithmetic. Reconstructing one from the other can move an off-grid sample.
//! `fs::read_to_string` is standard-library UTF-8 I/O, not afni-io. The parser
//! borrows that temporary text, but copies validated numbers into owned Events;
//! those events remain valid after the text buffer leaves scope.

use crate::{
    error::{Error, Result, reserve},
    model::{Acquisition, RunIndex},
};
use std::path::Path;

/// Explicit interpretation of onset seconds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TimingMode {
    /// One logical row per run; onsets are relative to that run's start.
    Local,
    /// Tokens on any noncomment row are seconds on the concatenated timeline.
    Global,
    /// AFNI inference: one token per row means global; otherwise local.
    Guess,
}

/// A single onset, before expansion into basis components.
#[derive(Debug)]
pub struct Event {
    pub(crate) run: RunIndex,
    pub(crate) onset: f64,
    pub(crate) timedex: f32,
    pub(crate) amplitudes: Vec<f32>,
    pub(crate) duration: Option<f32>,
}
impl Event {
    /// Zero-based run identity.
    pub fn run(&self) -> RunIndex {
        self.run
    }
    /// Original local onset in seconds, retained for manifest provenance.
    pub fn onset(&self) -> f64 {
        self.onset
    }
}

/// Validated events in input order, with the acquisition contract attached.
#[derive(Debug)]
pub struct Timing {
    events: Vec<Event>,
    lengths: Vec<usize>,
    tr: f64,
    pub(crate) modulators: usize,
    pub(crate) has_duration: bool,
}
/// Input layout, independent of the later column-expansion mode.
#[derive(Debug, Clone, Copy, Default)]
pub enum TimingFormat {
    /// Plain onset tokens, with no auxiliary fields.
    #[default]
    Plain,
    /// Married amplitudes, optionally followed by a nonlinear duration.
    Married {
        /// Require one duration per event for a dm response.
        duration: bool,
    },
    /// FSL onset/duration/amplitude rows separated into runs by '*'.
    Fsl {
        /// FS1 ignores the supplied amplitude and uses one.
        unit_amplitude: bool,
    },
}

/// Explicit preprocessing. Defaults preserve the M2 plain-seconds API.
#[derive(Debug, Clone, Copy, Default)]
pub struct TimingOptions {
    /// Input syntax to decode before assigning events to runs.
    pub format: TimingFormat,
    /// Convert onset milliseconds to seconds; durations remain seconds.
    pub milliseconds: bool,
    /// Subtract this many seconds AFTER conversion.
    pub subtract: f32,
    /// Permit empty conditions and explicit -1 placeholders (not bad onsets).
    pub allow_empty: bool,
}

impl Timing {
    /// Parse plain onset text; errors name source and physical line/token.
    ///
    /// Local onsets must be within 0..(run_length-1)*TR, matching AFNI's last
    /// sampled time. Global onsets are accepted through the last global sample.
    /// Negative, nonfinite, married amplitude/duration, and mixed '*' rows fail.
    /// Empty individual runs are legal; an entirely empty condition is rejected
    /// in this prototype rather than inventing an IM coefficient.
    pub fn parse(
        text: &str,
        source: &str,
        mode: TimingMode,
        acquisition: &Acquisition,
    ) -> Result<Self> {
        Self::parse_with_options(text, source, mode, acquisition, TimingOptions::default())
    }

    /// Decode auxiliaries and transform onsets BEFORE run assignment. Values
    /// are owned, so the returned events never borrow the input text. Amplitude
    /// centering is deferred to EventStimulus, which knows AM1 versus AM2.
    pub fn parse_with_options(
        text: &str,
        source: &str,
        mode: TimingMode,
        acquisition: &Acquisition,
        options: TimingOptions,
    ) -> Result<Self> {
        let fail = |message: &str| Error::input(format!("{source}: {message}"));
        let tr = acquisition.tr() as f32;
        if !tr.is_normal()
            || acquisition.nrows() > 1 << 24
            || !(acquisition.nrows() as f32 * tr).is_finite()
            || !options.subtract.is_finite()
        {
            return Err(fail(
                "event timeline/offset must be representable in finite AFNI float coordinates",
            ));
        }
        let mut rows: Vec<Vec<String>> = text
            .lines()
            .filter_map(|line| {
                let line = line.split('#').next().unwrap_or("").trim();
                if line.is_empty() {
                    None
                } else {
                    Some(line.split_whitespace().map(str::to_owned).collect())
                }
            })
            .collect();
        let has_duration = matches!(
            options.format,
            TimingFormat::Married { duration: true } | TimingFormat::Fsl { .. }
        );
        if let TimingFormat::Fsl { unit_amplitude } = options.format {
            // Convert FSL rows to the same married-token representation. '*'
            // separates runs; two consecutive separators preserve an empty run.
            let mut runs = vec![Vec::new()];
            for row in rows {
                if row == ["*"] {
                    runs.push(Vec::new());
                    continue;
                }
                if row.is_empty() || row.len() > 3 {
                    return Err(fail("FSL requires 1..3 fields per event"));
                }
                let duration = row.get(1).map(String::as_str).unwrap_or("0.1");
                let amplitude = if unit_amplitude {
                    "1"
                } else {
                    row.get(2).map(String::as_str).unwrap_or("1")
                };
                runs.last_mut()
                    .unwrap()
                    .push(format!("{}*{amplitude}:{duration}", row[0]));
            }
            for row in &mut runs {
                if row.is_empty() {
                    row.push("*".into());
                }
            }
            rows = runs;
        }
        let mode = if mode == TimingMode::Guess {
            if rows.iter().all(|r| r.len() <= 1) {
                TimingMode::Global
            } else {
                TimingMode::Local
            }
        } else {
            mode
        };
        if mode == TimingMode::Local && rows.len() != acquisition.run_lengths().len() {
            return Err(fail(
                "local timing requires exactly one row per run, including '*' for empty runs",
            ));
        }
        let mut events = Vec::new();
        let mut width = None;
        for (local_run, row) in rows.iter().enumerate() {
            if row == &["*"] {
                continue;
            }
            if row.iter().any(|s| s == "*") {
                return Err(fail("'*' must be alone on an empty timing row"));
            }
            for token in row {
                let number = |s: &str| -> Result<f32> {
                    let v = s
                        .parse::<f32>()
                        .map_err(|_| fail("invalid timing/auxiliary number"))?;
                    if !v.is_finite() {
                        return Err(fail("nonfinite timing/auxiliary value"));
                    }
                    Ok(v)
                };
                let (front, duration) = if let Some((a, d)) = token.split_once(':') {
                    (a, Some(number(d)?))
                } else {
                    (token.as_str(), None)
                };
                let (onset, amps) = if let Some((t, a)) = front.split_once('*') {
                    (t, a.split(',').map(number).collect::<Result<Vec<_>>>()?)
                } else {
                    (front, Vec::new())
                };
                if matches!(options.format, TimingFormat::Plain)
                    && (!amps.is_empty() || duration.is_some())
                {
                    return Err(fail(
                        "plain timing does not accept amplitude/duration fields",
                    ));
                }
                if duration.is_some() != has_duration
                    || duration.is_some_and(|d| !(0.0..=999.0).contains(&d))
                {
                    return Err(fail(
                        "duration must be present exactly for dm models and lie in 0..999 seconds",
                    ));
                }
                if amps.len() + usize::from(has_duration) > 21
                    || width.is_some_and(|n| n != amps.len())
                {
                    return Err(fail(
                        "inconsistent auxiliary count or more than 21 auxiliaries",
                    ));
                }
                width = Some(amps.len());
                let raw = number(onset)?;
                if raw == -1.0 && options.allow_empty {
                    continue;
                }
                let time = raw * if options.milliseconds { 0.001 } else { 1.0 } - options.subtract;
                if raw < 0.0 || time < 0.0 || !time.is_finite() {
                    return Err(fail(
                        "onset must remain finite and nonnegative after transformation",
                    ));
                }
                // Preserve original decimal onset for untransformed manifests;
                // scientific sampling always uses AFNI float coordinates.
                let seconds = if !options.milliseconds && options.subtract == 0.0 {
                    onset.parse::<f64>().map_err(|_| fail("invalid onset"))?
                } else {
                    f64::from(time)
                };
                let (run, onset, timedex) = if mode == TimingMode::Local {
                    if time > (acquisition.run_lengths()[local_run] - 1) as f32 * tr {
                        return Err(fail("local onset exceeds the last sampled time"));
                    }
                    (
                        local_run,
                        seconds,
                        time / tr + acquisition.run_starts()[local_run].0 as f32,
                    )
                } else {
                    if time > (acquisition.nrows() - 1) as f32 * tr {
                        return Err(fail("global onset exceeds the last sampled time"));
                    }
                    let index = time / tr;
                    let run = acquisition
                        .run_starts()
                        .partition_point(|s| s.0 as f32 <= index)
                        - 1;
                    (
                        run,
                        seconds - acquisition.run_starts()[run].0 as f64 * acquisition.tr(),
                        index,
                    )
                };
                let range = acquisition
                    .run_range(RunIndex(run))
                    .ok_or_else(|| fail("invalid event run"))?;
                if onset < 0.0 || timedex < range.start as f32 || timedex >= range.end as f32 {
                    return Err(fail("onset rounds across a run boundary"));
                }
                reserve(&mut events, 1)?;
                events.push(Event {
                    run: RunIndex(run),
                    onset,
                    timedex,
                    amplitudes: amps,
                    duration,
                });
            }
        }
        if events.is_empty() && !options.allow_empty {
            return Err(fail("event condition is entirely empty"));
        }
        Ok(Self {
            events,
            lengths: acquisition.run_lengths().to_vec(),
            tr: acquisition.tr(),
            modulators: width.unwrap_or(0),
            has_duration,
        })
    }
    /// Borrow source-order events, never sorted by onset.
    pub fn events(&self) -> &[Event] {
        &self.events
    }
    /// Check the stored acquisition contract before any event-column indexing.
    pub(crate) fn validate_for(&self, acquisition: &Acquisition) -> Result<()> {
        if self.lengths != acquisition.run_lengths() || self.tr != acquisition.tr() {
            return Err(Error::input(
                "timing was validated for a different acquisition",
            ));
        }
        Ok(())
    }
}

/// Read plain onset timing after input.rs applies selectors or an inline literal.
pub fn read_timing(path: &Path, mode: TimingMode, acquisition: &Acquisition) -> Result<Timing> {
    read_timing_with_options(path, mode, acquisition, TimingOptions::default())
}

/// Shared input adapter: selectors and inline literals are decoded in input.rs,
/// then timing validation owns the scientific interpretation. No subprocess or
/// afni-io reader is called; the returned Timing owns all parsed numbers.
pub fn read_timing_with_options(
    path: &Path,
    mode: TimingMode,
    acquisition: &Acquisition,
    options: TimingOptions,
) -> Result<Timing> {
    let rows = crate::input::rows(path, true)?;
    let text = rows
        .iter()
        .map(|r| r.join(" "))
        .collect::<Vec<_>>()
        .join("\n");
    Timing::parse_with_options(&text, &path.to_string_lossy(), mode, acquisition, options)
}
