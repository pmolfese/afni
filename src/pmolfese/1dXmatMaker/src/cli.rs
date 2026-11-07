// Developed at the NIH
//
// Part of AFNI. See the repository LICENSE.txt for licensing terms and
// disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
// Third-party and adapted portions retain their existing attribution and license
// terms; see THIRD_PARTY.md.

//! AFNI-style single-dash CLI for ordinary event/baseline/nuisance generation.
//!
//! Read parse() for the accepted grammar, then execute() for file loading.
//! PendingStim stores incomplete options until every declared indexed stimulus
//! has a source, label, and compatible role. No library object is built from these
//! incomplete records. `Option<T>` means a setting has not yet been supplied;
//! unlike a null pointer, it must be handled explicitly before accessing T.
//!
//! `use afni_xmat_maker::{...}` imports THIS package's library, not an
//! external AFNI crate. In execute(), read_table/read_timing load owned local
//! types; build_with_events performs the science; write_outputs_with_manifest
//! publishes the results. There are no afni-io/afni-core calls to explain yet.
//! An import only brings a name into scope; it does not execute that function.
//! Keep file operations here or in dedicated I/O modules, not inside bases.

use afni_xmat_maker::{
    VERSION,
    basis::Basis,
    design::{EventMode, EventStimulus, Nuisance, NuisanceKind, Source, build_with_events},
    error::{Error, Result},
    model::Acquisition,
    table::{read_censor, read_table},
    timing::{TimingFormat, TimingMode, TimingOptions, read_timing_with_options},
    xmat::write_all_outputs,
};
use std::{env, path::PathBuf};

// A raw string preserves shell backslashes and quotes in the worked examples.
// Keep input sizes, run indexing and output descriptions beside each command.
const HELP: &str = r#"1dXmatMaker: independent AFNI design-matrix generator (0.8 standalone release)

Seven worked shell examples follow the option reference below.

Acquisition and output (required):
  -run_lengths N,N,... -TR SECONDS   Original sample counts and common TR.
  -nodata NT TR [-concat FILE]      Alternative spelling; concat lists starts.
  -polort P                        Integer >=-1 or A (automatic order).
  -num_stimts N                    Indexed source count; 0 is allowed.
  -x1D FILE                        Retained-row matrix; never overwritten.

Indexed sources (K starts at 1; each requires -stim_label K LABEL):
  -stim_file K FILE                One sampled signal column.
  -stim_base K                     Classify sampled source as baseline nuisance.
  -stim_minlag K N -stim_maxlag K N -stim_nptr K N
                                  Nonnegative lag limits, samples per TR.
  -stim_times K FILE BASIS         Ordinary condition response.
  -stim_times_IM K FILE BASIS      One response set per event; dm allowed.
  -stim_times_AM1 K FILE BASIS     Sum amplitude-weighted responses.
  -stim_times_AM2 K FILE BASIS     Mean plus centered modulator responses.
  -stim_times_AMx K FILE BASIS     AM2 alias.
  -stim_times_FSL K FILE BASIS     FSL onset/duration/amplitude rows; dm required.
  -stim_times_FS1 K FILE BASIS     FSL with amplitudes forced to one.
  AM2/AMx optionally take :CENTER:x:... after BASIS; x keeps that mean.
  Married syntax: onset*amp1,amp2:duration (duration only for dm bases).

Timing settings (before sources):
  -local_times / -global_times / -guess_times   Explicit mode required for events.
  -stim_times_millisec             Convert ONSETS only; durations stay seconds.
  -stim_times_subtract SECONDS     Subtract after unit conversion.
  -rawAM2                         Disable AM2 centering explicitly.
  -basis_normall PEAK              Normalize every following response component.
  -allzero_OK                     Allow empty ordinary/AM conditions; not empty IM.
  -expr_seed N                    Explicit local EXPR seed (default 305441741).

Baselines and row selection:
  -ortvec FILE LABEL / -inati FILE LABEL   Repeatable numeric nuisance tables.
  -legendre / -nolegendre          Legendre default or legacy sample-index powers.
                                  polort A forces Legendre, matching AFNI.
  -dmbase / -nodmbase              Mean removal; polort=-1 disables it.
  -censor FILE                    One full-row 0/1 column (1 keeps).
  -CENSORTR '1:0-2 2:4,10'         Run:local or global indices; *: selects all runs.
  -nfirst N -nlast N               Inclusive per-run limits. Default nfirst is
                                  maximum ceil(maxlag/nptr); explicit 0 overrides.

Contrasts and companions:
  -num_glt N                      Optional exact declared number of tests.
  -gltsym 'SYM: +event[0] -event[1]' / -gltsym FILE
  -glt NROWS FILE                 Numeric rows covering ALL design columns.
  -glt_label K LABEL              Test label (K starts at 1).
  Selectors: [0..$], [0,2], [*]; [[0..2]] expands rows. | separates rows.
  Ort addresses polynomial columns; indexed source names include baseline inputs.
  Use numeric GLTs to address ortvec columns. Repeated symbolic terms assign.
  -x1D_uncensored FILE            Full-row matrix without GLT metadata.
  -x1D_regcensored FILE            Full rows plus one nuisance spike per removed row.
  -event_manifest FILE            Optional version-1 event-to-column TSV.
  -help / -h / -version           Use alone.

Put acquisition, timing and normalization settings before sources. Unknown,
repeated singleton and incompatible options fail. nfirst/nlast/CENSORTR/censor
combine AFTER full matrix construction, demeaning and amplitude centering.
Event responses stop at run edges; sampled lags can cross those edges.

Text input supports # comments, [column] and {row} selectors, trailing transpose,
and single-column 1D: literals with n@value repetition. Selectors support ranges,
positive steps and $. No padding/truncation or stdin. Out-of-range events fail.
Local timing has one row per run; '*' alone is empty. FSL '*' separates runs.
All-zero/redundant/wide matrices are retained without OLS rank rejection.
Output parents must exist. Files are staged together and never overwritten;
ordinary publication errors roll back new files (not a multi-file atomic commit).

Bases (seconds): GAM, GAMpw, TWOGAM, TWOGAMpw; BLOCK/IGFUN/UBLOCK[4/5];
dmBLOCK/dmUBLOCK[4/5]; TENT[zero], CSPLIN[zero], POLY, SIN, TRIG;
SPMG[1/2/3], MION/MIONN, WAV, EXPR. Quote parameterized basis arguments.
See README.md for basis parameters, EXPR policy and compatibility notes.
Candidate mode (separate from indexed stimulus options):
  -candidate_library glmsingle|LIBRARY.tsv -events EVENTS.tsv
  -run_lengths N,N -TR SECONDS -polort INTEGER -candidate_prefix PREFIX
  [-candidate_ids 0,1,19] [-hrf_normalization peak|library|none]
  [-censor FILE] [-ortvec FILE LABEL ...]
Events header: trial<TAB>condition<TAB>run<TAB>onset<TAB>duration<TAB>item.
Run IDs are zero-based. Fixed duration must agree across events. Outputs include
ordinary trial/ONOFF matrices, prepared waveforms, Picard catalogs, column IDs
and trial_info. See README.md for custom library grammar and provenance.

No silent waveform substitutions. No dataset loading or fitting. No external
AFNI commands/libraries or afni-io/afni-core calls. AFNI environment settings and
.afnirc are ignored. The real GLMsingle library is embedded in the executable.
Exit codes: 0 success, 2 invalid input, 1 I/O/resource failure.

Worked examples (POSIX shell)
=============================

Run these in a new writable directory with 1dXmatMaker on PATH. The printf
commands create small illustrative input files; use your study's actual timing
and acquisition settings for analysis. Output names must not already exist.
Examples 4 and 6 note their additional inputs or indexing conventions below.
All onsets and durations are in seconds; -run_lengths counts samples, not seconds.
These commands generate matrices only; they do not fit image data.

1) Baseline-only design for two runs of unequal length

  TR is 2 seconds; each run gets a constant, linear trend and quadratic trend.
  The result has 220 rows and six columns. No timing or image file is needed.

    1dXmatMaker -run_lengths 120,100 -TR 2 -polort 2 \
      -num_stimts 0 -x1D ex1.baseline.xmat.1D

2) Two conditions, local timing, and an A-minus-B contrast

  Each timing file has one line per run. Onsets restart at zero for run 2.
  GAM supplies one response column per condition; polort 1 adds two baseline
  columns per run. The 240-row, six-column matrix stores the contrast as a GLT.

    printf '10 40 70\n12 42 72\n' > ex2.A.1D
    printf '25 55 85\n27 57 87\n' > ex2.B.1D
    1dXmatMaker -run_lengths 120,120 -TR 1 -polort 1 \
      -num_stimts 2 -local_times \
      -stim_times 1 ex2.A.1D GAM -stim_label 1 A \
      -stim_times 2 ex2.B.1D GAM -stim_label 2 B \
      -gltsym 'SYM: +A -B' -glt_label 1 A_minus_B \
      -x1D ex2.conditions.xmat.1D

3) Individual trials (IM), censoring, and a full-row companion

  Four trials per run produce eight trial columns plus four baseline columns.
  -CENSORTR uses ONE-based run numbers and ZERO-based sample numbers: remove
  samples 0, 1 and 2 from each run. The retained matrix has 234 rows; the full
  companion has all 240 rows. Censoring does not shift event times or identities.
  The event manifest maps events to design columns; it is not a beta manifest.

    printf '10 35 60 85\n12 37 62 87\n' > ex3.trials.1D
    1dXmatMaker -run_lengths 120,120 -TR 1 -polort 1 \
      -num_stimts 1 -local_times -CENSORTR '1:0-2 2:0-2' \
      -stim_times_IM 1 ex3.trials.1D GAM -stim_label 1 trial \
      -x1D ex3.trials.xmat.1D -x1D_uncensored ex3.trials.full.xmat.1D \
      -event_manifest ex3.events.tsv

4) Add motion nuisance regressors and a censor file to example 2

  First run example 2 to create ex2.A.1D and ex2.B.1D. Supply motion.1D with
  240 rows and six numeric motion columns, and keep.1D with 240 rows of one
  0/1 value each (1 keeps the sample). Both inputs cover the ORIGINAL timeline,
  even when rows will be censored. Runs are concatenated in acquisition order.
  -ortvec classifies motion as nuisance; it does not increase -num_stimts.
  The result has 12 columns and as many rows as keep.1D has ones.

    1dXmatMaker -run_lengths 120,120 -TR 1 -polort 1 \
      -num_stimts 2 -local_times -ortvec motion.1D motion -censor keep.1D \
      -stim_times 1 ex2.A.1D GAM -stim_label 1 A \
      -stim_times 2 ex2.B.1D GAM -stim_label 2 B \
      -x1D ex4.motion.xmat.1D -x1D_uncensored ex4.motion.full.xmat.1D

5) Amplitude modulation with variable event durations

  Married timing tokens below are onset*rating:duration. AM2 produces a mean
  response and a centered rating-modulation response, plus the run baselines.
  dmBLOCK uses each event's duration. Quote its parameterized basis expression
  so the shell passes the parentheses unchanged. The result has six columns.

    printf '10*1:2 40*3:4 75*2:1\n12*2:2 42*4:4 77*1:1\n' > ex5.ratings.1D
    1dXmatMaker -run_lengths 120,120 -TR 1 -polort 1 \
      -num_stimts 1 -local_times \
      -stim_times_AM2 1 ex5.ratings.1D 'dmBLOCK(1)' -stim_label 1 rating \
      -x1D ex5.ratings.xmat.1D -event_manifest ex5.events.tsv

6) Generate a small GLMsingle HRF candidate collection

  printf writes real TAB separators for the explicit event table. Its run IDs
  are ZERO-based (unlike -CENSORTR above). Trial IDs are unique across runs;
  condition and item IDs may repeat. All six events here last one second.

    printf 'trial\tcondition\trun\tonset\tduration\titem\n' > ex6.events.tsv
    printf '100\t0\t0\t10\t1\tface\n' >> ex6.events.tsv
    printf '101\t1\t0\t40\t1\thouse\n' >> ex6.events.tsv
    printf '102\t2\t0\t70\t1\ttool\n' >> ex6.events.tsv
    printf '200\t0\t1\t12\t1\tface\n' >> ex6.events.tsv
    printf '201\t1\t1\t42\t1\thouse\n' >> ex6.events.tsv
    printf '202\t2\t1\t72\t1\ttool\n' >> ex6.events.tsv
    1dXmatMaker -candidate_library glmsingle -candidate_ids 0,10,19 \
      -events ex6.events.tsv -run_lengths 120,120 -TR 1 -polort 1 \
      -candidate_prefix ex6.hrf

  This writes trial and condition-collapsed ONOFF matrices for each candidate,
  their full-row companions, prepared waveforms, and shared identity tables.
  For example, ex6.hrf.id10.trials.xmat.1D has six trial and four baseline columns;
  ex6.hrf.id10.onoff.xmat.1D has one ONOFF and four baseline columns.
  ex6.hrf.candidates.tsv is the trial-design catalog for Picard's HRF policy.
  Omit -candidate_ids to generate all 20 embedded HRFs. Default normalization
  makes each prepared HRF's positive peak one. Candidate generation requires
  a common event duration; the ordinary dmBLOCK example above allows it to vary.

7) Two conditions and two runs: 1dXmatMaker -> 3dPicard

  Assume two preprocessed runs, run1.nii.gz and run2.nii.gz, each with
  120 volumes and TR = 2 seconds. Each timing file has one line per run;
  onsets are in seconds relative to that run's start. Use new output names.

    # Create timing files: line 1 = run 1, line 2 = run 2.
    printf '10 30 50 70 90 110 130 150 170 190\n12 32 52 72 92 112 132 152 172 192\n' \
      > frequent.1D
    printf '60 120 180\n62 122 182\n' > target.1D

    # Generate condition regressors and separate polynomial baselines per run.
    1dXmatMaker -run_lengths 120,120 -TR 2 -polort 2 \
      -num_stimts 2 -local_times \
      -stim_times 1 frequent.1D GAM -stim_label 1 frequent \
      -stim_times 2 target.1D GAM -stim_label 2 target \
      -x1D frequent_target.xmat.1D

    # Quote the input list so AFNI concatenates the runs in the stated order.
    3dPicard -input 'run1.nii.gz run2.nii.gz' \
      -matrix frequent_target.xmat.1D \
      -fraction 0.6 \
      -prefix frequent_target.picard.nii.gz \
      -fitts frequent_target.fitts.nii.gz \
      -errts frequent_target.errts.nii.gz

  The matrix has eight columns: six unpenalized baseline columns (constant,
  linear and quadratic for each run), plus the frequent and target columns.
  Picard estimates one coefficient per condition shared across both runs;
  it uses the baseline columns already in the matrix rather than adding polort.
  The input order must match the timing-file rows and -run_lengths. Run
  boundaries come from the matrix. Quoted-list concatenation uses the standard
  AFNI reader settings (AFNI_PATH_SPACES_OK must not be YES).
  The fraction 0.6 is an illustrative fixed choice, not a tuned value.

For fitting these matrices with LSS, REML or Picard, see README.md.
For custom sampled HRFs and candidate provenance, see README.md.
"#;

#[derive(Default)]
struct PendingStim {
    source: Option<PendingSource>,
    label: Option<String>,
    baseline: bool,
    minlag: Option<usize>,
    maxlag: Option<usize>,
    nptr: Option<usize>,
}

/// A source is either a sampled nuisance file or timing plus a parsed basis.
/// Keeping this as an enum makes competing -stim_file/-stim_times definitions
/// a duplicate, rather than leaving two inconsistent optional fields.
enum PendingSource {
    File(PathBuf),
    Events {
        path: PathBuf,
        basis: Basis,
        mode: EventMode,
        centers: Option<Vec<Option<f32>>>,
    },
}

/// Pending GLTs are resolved only after all matrix columns exist.
enum PendingGlt {
    Symbolic(String),
    Numeric(usize, PathBuf),
}

// Parsing owns paths and settings, but does not yet allocate the design.
// This boundary lets execute() validate every source before publishing files.
struct Options {
    regcen: Option<PathBuf>,
    raw: bool,
    milliseconds: bool,
    subtract: f32,
    allow_empty: bool,
    normalization: Option<f32>,
    first: usize,
    last: Option<usize>,
    censortr: Vec<String>,
    glts: Vec<PendingGlt>,
    glt_labels: Vec<(usize, String)>,
    acquisition: Acquisition,
    stimuli: Vec<PendingStim>,
    ortvec: Vec<(PathBuf, String)>,
    censor: Option<PathBuf>,
    primary: PathBuf,
    full: Option<PathBuf>,
    manifest: Option<PathBuf>,
    timing_mode: Option<TimingMode>,
}

/// Borrow the next token from the argument vector; no copied argument buffers.
fn argument<'a>(args: &'a [String], cursor: &mut usize, option: &str) -> Result<&'a str> {
    let value = args
        .get(*cursor)
        .ok_or_else(|| Error::input(format!("{option} needs another argument")))?;
    *cursor += 1;
    Ok(value)
}

/// Set a singleton option once. Repeating it is an error, even if values agree.
fn set_once<T>(slot: &mut Option<T>, value: T, option: &str) -> Result<()> {
    if slot.is_some() {
        return Err(Error::input(format!("duplicate {option}")));
    }
    *slot = Some(value);
    Ok(())
}

/// Require an explicitly supplied value rather than inventing a scientific default.
fn required<T>(value: Option<T>, option: &str) -> Result<T> {
    value.ok_or_else(|| Error::input(format!("missing required {option}")))
}

/// Convert a user stimulus number into an in-range zero-based vector position.
fn stimulus_index(
    args: &[String],
    cursor: &mut usize,
    option: &str,
    stimuli: &Option<Vec<PendingStim>>,
) -> Result<usize> {
    let value = argument(args, cursor, option)?;
    let index = value.parse::<usize>().ok().and_then(|k| k.checked_sub(1));
    let entries = stimuli
        .as_ref()
        .ok_or_else(|| Error::input("-num_stimts must precede indexed stimulus options"))?;
    index.filter(|&i| i < entries.len()).ok_or_else(|| {
        Error::input(format!(
            "{option}: stimulus '{value}' must be in 1..{}",
            entries.len()
        ))
    })
}

/// Parse the ordinary matrix grammar, preserving index/option context in errors.
fn parse(args: &[String]) -> Result<Options> {
    let (mut lengths, mut tr, mut polort, mut demean) = (None, None, None, None);
    let (mut stimuli, mut censor, mut primary, mut full) = (None, None, None, None);
    let (mut manifest, mut timing_mode) = (None, None);
    let (mut nodata, mut concat) = (None, None);
    let (mut raw, mut milliseconds, mut allow_empty) = (false, false, false);
    let (mut subtract, mut normalization, mut first, mut last) = (None, None, None, None);
    let (mut regcen, mut num_glt) = (None, None);
    let (mut censortr, mut glts, mut glt_labels) = (Vec::new(), Vec::new(), Vec::new());
    let mut power = false;
    let mut expression_seed: Option<u32> = None;
    let mut ortvec = Vec::new();
    let mut started_sources = false;
    let mut legendre_seen = false;
    let mut cursor = 0;
    while cursor < args.len() {
        let option = argument(args, &mut cursor, "option")?;
        if started_sources
            && matches!(
                option,
                "-run_lengths"
                    | "-TR"
                    | "-polort"
                    | "-num_stimts"
                    | "-legendre"
                    | "-dmbase"
                    | "-nodmbase"
                    | "-local_times"
                    | "-global_times"
                    | "-expr_seed"
                    | "-nodata"
                    | "-concat"
                    | "-nolegendre"
                    | "-guess_times"
                    | "-rawAM2"
                    | "-stim_times_millisec"
                    | "-stim_times_subtract"
                    | "-basis_normall"
                    | "-allzero_OK"
            )
        {
            return Err(Error::input(format!(
                "{option} must precede stimulus definitions"
            )));
        }
        match option {
            "-nodata" => {
                let n = argument(args, &mut cursor, option)?
                    .parse::<usize>()
                    .map_err(|_| Error::input("invalid nodata length"))?;
                let dt = argument(args, &mut cursor, option)?
                    .parse::<f64>()
                    .map_err(|_| Error::input("invalid nodata TR"))?;
                set_once(&mut nodata, n, option)?;
                set_once(&mut tr, dt, "-TR/-nodata")?;
            }
            "-concat" => set_once(
                &mut concat,
                PathBuf::from(argument(args, &mut cursor, option)?),
                option,
            )?,
            "-rawAM2" | "-stim_times_millisec" | "-allzero_OK" => {
                let flag = match option {
                    "-rawAM2" => &mut raw,
                    "-allzero_OK" => &mut allow_empty,
                    _ => &mut milliseconds,
                };
                if *flag {
                    return Err(Error::input(format!("duplicate {option}")));
                }
                *flag = true;
            }
            "-stim_times_subtract" | "-basis_normall" => {
                let v = argument(args, &mut cursor, option)?
                    .parse::<f32>()
                    .map_err(|_| Error::input("invalid timing/normalization value"))?;
                if !v.is_finite() || (option == "-basis_normall" && v <= 0.0) {
                    return Err(Error::input("invalid timing/normalization value"));
                }
                set_once(
                    if option == "-basis_normall" {
                        &mut normalization
                    } else {
                        &mut subtract
                    },
                    v,
                    option,
                )?;
            }
            "-nfirst" | "-nlast" => {
                let v = argument(args, &mut cursor, option)?
                    .parse::<usize>()
                    .map_err(|_| Error::input("invalid first/last sample"))?;
                set_once(
                    if option == "-nfirst" {
                        &mut first
                    } else {
                        &mut last
                    },
                    v,
                    option,
                )?;
            }
            "-CENSORTR" => censortr.push(argument(args, &mut cursor, option)?.to_owned()),
            "-x1D_regcensored" => set_once(
                &mut regcen,
                PathBuf::from(argument(args, &mut cursor, option)?),
                option,
            )?,
            "-num_glt" => set_once(
                &mut num_glt,
                argument(args, &mut cursor, option)?
                    .parse::<usize>()
                    .map_err(|_| Error::input("invalid GLT count"))?,
                option,
            )?,
            "-gltsym" => glts.push(PendingGlt::Symbolic(
                argument(args, &mut cursor, option)?.to_owned(),
            )),
            "-glt" => {
                let n = argument(args, &mut cursor, option)?
                    .parse::<usize>()
                    .map_err(|_| Error::input("invalid GLT rows"))?;
                glts.push(PendingGlt::Numeric(
                    n,
                    PathBuf::from(argument(args, &mut cursor, option)?),
                ));
            }
            "-glt_label" => {
                let n = argument(args, &mut cursor, option)?
                    .parse::<usize>()
                    .ok()
                    .and_then(|n| n.checked_sub(1))
                    .ok_or_else(|| Error::input("GLT labels start at 1"))?;
                if glt_labels.iter().any(|(i, _)| *i == n) {
                    return Err(Error::input("duplicate GLT label index"));
                }
                glt_labels.push((n, argument(args, &mut cursor, option)?.to_owned()));
            }
            "-expr_seed" => {
                let value = argument(args, &mut cursor, option)?;
                let seed = value.parse::<u32>().map_err(|_| {
                    Error::input("-expr_seed requires a decimal integer in 0..4294967295")
                })?;
                set_once(&mut expression_seed, seed, option)?;
            }
            "-run_lengths" => {
                let value = argument(args, &mut cursor, option)?;
                let parsed = value
                    .split(',')
                    .map(|part| {
                        part.parse::<usize>()
                            .map_err(|_| Error::input(format!("invalid run length '{part}'")))
                    })
                    .collect::<Result<Vec<_>>>()?;
                set_once(&mut lengths, parsed, option)?;
            }
            "-TR" => {
                let value = argument(args, &mut cursor, option)?;
                let parsed = value
                    .parse::<f64>()
                    .map_err(|_| Error::input(format!("invalid TR '{value}'")))?;
                set_once(&mut tr, parsed, option)?;
            }
            "-polort" => {
                let value = argument(args, &mut cursor, option)?;
                let parsed = value.to_owned();
                set_once(&mut polort, parsed, option)?;
            }
            "-dmbase" | "-nodmbase" => {
                set_once(&mut demean, option == "-dmbase", "-dmbase/-nodmbase")?
            }
            "-legendre" | "-nolegendre" => {
                if legendre_seen {
                    return Err(Error::input("duplicate -legendre"));
                }
                legendre_seen = true;
                power = option == "-nolegendre";
            }
            "-num_stimts" => {
                let value = argument(args, &mut cursor, option)?;
                let count = value
                    .parse::<usize>()
                    .ok()
                    .filter(|&n| n <= i32::MAX as usize)
                    .ok_or_else(|| {
                        Error::input("-num_stimts must be a nonnegative 32-bit integer")
                    })?;
                let mut entries = Vec::new();
                entries.try_reserve_exact(count).map_err(|e| {
                    Error::resource(format!("cannot allocate stimulus records: {e}"))
                })?;
                entries.resize_with(count, PendingStim::default);
                set_once(&mut stimuli, entries, option)?;
            }
            "-local_times" | "-global_times" | "-guess_times" => {
                set_once(
                    &mut timing_mode,
                    if option == "-local_times" {
                        TimingMode::Local
                    } else if option == "-guess_times" {
                        TimingMode::Guess
                    } else {
                        TimingMode::Global
                    },
                    "timing mode",
                )?;
            }
            "-event_manifest" => set_once(
                &mut manifest,
                PathBuf::from(argument(args, &mut cursor, option)?),
                option,
            )?,
            "-stim_file" | "-stim_label" | "-stim_base" | "-stim_times" | "-stim_times_IM"
            | "-stim_times_AM1" | "-stim_times_AM2" | "-stim_times_AMx" | "-stim_times_FSL"
            | "-stim_times_FS1" | "-stim_minlag" | "-stim_maxlag" | "-stim_nptr" => {
                started_sources = true;
                let index = stimulus_index(args, &mut cursor, option, &stimuli)?;
                // The index check above also established Some(stimuli). Match
                // explicitly to keep user-input handling free of unwrap panics.
                let entries = stimuli
                    .as_mut()
                    .ok_or_else(|| Error::input("missing -num_stimts"))?;
                let entry = &mut entries[index];
                match option {
                    "-stim_minlag" | "-stim_maxlag" | "-stim_nptr" => {
                        let n = argument(args, &mut cursor, option)?
                            .parse::<usize>()
                            .map_err(|_| Error::input("invalid lag/sampling rate"))?;
                        let slot = match option {
                            "-stim_minlag" => &mut entry.minlag,
                            "-stim_maxlag" => &mut entry.maxlag,
                            _ => &mut entry.nptr,
                        };
                        set_once(slot, n, option)?;
                    }
                    "-stim_file" => set_once(
                        &mut entry.source,
                        PendingSource::File(PathBuf::from(argument(args, &mut cursor, option)?)),
                        option,
                    )?,
                    "-stim_times" | "-stim_times_IM" | "-stim_times_AM1" | "-stim_times_AM2"
                    | "-stim_times_AMx" | "-stim_times_FSL" | "-stim_times_FS1" => {
                        let path = PathBuf::from(argument(args, &mut cursor, option)?);
                        let spec = argument(args, &mut cursor, option)?;
                        // Colons after an AM2 basis select explicit modulator centers.
                        // EXPR can contain arbitrary text, so only paired modes split.
                        let (basis_text, mut centers) =
                            if matches!(option, "-stim_times_AM2" | "-stim_times_AMx") {
                                if let Some((b, c)) = spec.split_once(':') {
                                    let values = c
                                        .split(':')
                                        .map(|v| {
                                            if v == "x" {
                                                Ok(None)
                                            } else {
                                                v.parse::<f32>()
                                                    .map(Some)
                                                    .map_err(|_| Error::input("invalid AM2 center"))
                                            }
                                        })
                                        .collect::<Result<Vec<_>>>()?;
                                    (b, Some(values))
                                } else {
                                    (spec, None)
                                }
                            } else {
                                (spec, None)
                            };
                        // AFNI's documented spelling uses a separate :center token.
                        // Accept the attached spelling too, but never both.
                        if matches!(option, "-stim_times_AM2" | "-stim_times_AMx")
                            && args.get(cursor).is_some_and(|v| v.starts_with(':'))
                        {
                            if centers.is_some() {
                                return Err(Error::input("duplicate AM2 centers"));
                            }
                            let text = argument(args, &mut cursor, option)?;
                            centers = Some(
                                text[1..]
                                    .split(':')
                                    .map(|v| {
                                        if v == "x" {
                                            Ok(None)
                                        } else {
                                            v.parse::<f32>()
                                                .map(Some)
                                                .map_err(|_| Error::input("invalid AM2 center"))
                                        }
                                    })
                                    .collect::<Result<Vec<_>>>()?,
                            );
                        }
                        let mut basis = Basis::parse(basis_text)?;
                        // Offset by original stimulus index, not option order.
                        // This keeps random columns reproducible when options
                        // are reordered and avoids identical source streams.
                        basis.set_expression_seed(
                            expression_seed
                                .unwrap_or(0x1234_abcd)
                                .wrapping_add((index as u32).wrapping_mul(1024)),
                        );
                        let mode = match option {
                            "-stim_times_IM" => EventMode::Individual,
                            "-stim_times_AM1" => EventMode::Am1,
                            "-stim_times_AM2" => EventMode::Am2,
                            "-stim_times_AMx" => EventMode::Amx,
                            "-stim_times_FSL" => EventMode::Fsl,
                            "-stim_times_FS1" => EventMode::Fs1,
                            _ => EventMode::Ordinary,
                        };
                        set_once(
                            &mut entry.source,
                            PendingSource::Events {
                                path,
                                basis,
                                mode,
                                centers,
                            },
                            "stimulus source",
                        )?;
                    }
                    "-stim_label" => set_once(
                        &mut entry.label,
                        argument(args, &mut cursor, option)?.to_owned(),
                        option,
                    )?,
                    _ => {
                        if entry.baseline {
                            return Err(Error::input(format!(
                                "duplicate -stim_base {}",
                                index + 1
                            )));
                        }
                        entry.baseline = true;
                    }
                }
            }
            "-ortvec" | "-inati" => {
                started_sources = true;
                if stimuli.is_none() {
                    return Err(Error::input(
                        "-num_stimts must precede -ortvec (use 0 when appropriate)",
                    ));
                }
                let path = PathBuf::from(argument(args, &mut cursor, option)?);
                let label = argument(args, &mut cursor, option)?.to_owned();
                ortvec.push((path, label));
            }
            "-censor" => set_once(
                &mut censor,
                PathBuf::from(argument(args, &mut cursor, option)?),
                option,
            )?,
            "-x1D" => set_once(
                &mut primary,
                PathBuf::from(argument(args, &mut cursor, option)?),
                option,
            )?,
            "-x1D_uncensored" => set_once(
                &mut full,
                PathBuf::from(argument(args, &mut cursor, option)?),
                option,
            )?,
            _ => {
                return Err(Error::input(format!(
                    "unsupported option '{option}'; see -help and README.md"
                )));
            }
        }
    }
    let tr = required(tr, "-TR or -nodata")?;
    if let Some(n) = nodata {
        if lengths.is_some() {
            return Err(Error::input("-nodata conflicts with -run_lengths"));
        }
        let starts = if let Some(path) = concat.as_deref() {
            afni_xmat_maker::input::rows(path, false)?
                .into_iter()
                .flatten()
                .map(|v| {
                    v.parse::<usize>()
                        .map_err(|_| Error::input("concat starts must be integers"))
                })
                .collect::<Result<Vec<_>>>()?
        } else {
            vec![0]
        };
        if starts.first() != Some(&0)
            || starts.windows(2).any(|w| w[0] >= w[1])
            || starts.last().is_none_or(|&v| v >= n)
        {
            return Err(Error::input(
                "concat starts must increase from zero within nodata length",
            ));
        }
        let ends = starts.iter().skip(1).copied().chain(std::iter::once(n));
        lengths = Some(starts.iter().zip(ends).map(|(&s, e)| e - s).collect());
    } else if concat.is_some() {
        return Err(Error::input("-concat requires -nodata"));
    }
    let lengths = required(lengths, "-run_lengths or -nodata")?;
    let polort = required(polort, "-polort")?;
    // The reference explicitly resets power mode for automatic order. State
    // this visible compatibility rule rather than silently using a new model.
    if polort == "A" && power {
        eprintln!(
            "++ -polort A selects Legendre baselines, including with -nolegendre (AFNI compatibility)"
        );
        power = false;
    }
    let polort = if polort == "A" {
        // AFNI computes duration in float, then myfloor adds 0.00005 before
        // truncation. That tolerance changes the degree just below 150-second
        // boundaries, so ordinary f64 floor would choose a different model.
        let duration = (lengths.iter().copied().max().unwrap_or(0) as f32) * (tr as f32);
        let order = (f64::from(duration) / 150.0 + 0.00005).floor();
        if !order.is_finite() || order < 0.0 || order >= (i32::MAX - 1) as f64 {
            return Err(Error::input("automatic polort exceeds finite index limits"));
        }
        1 + order as i32
    } else {
        polort
            .parse::<i32>()
            .map_err(|_| Error::input("polort requires integer or A"))?
    };
    let acquisition = Acquisition::new(lengths, tr, polort, demean)?.with_power_baseline(power);
    if num_glt.is_some_and(|n| n != glts.len()) || glt_labels.iter().any(|(i, _)| *i >= glts.len())
    {
        return Err(Error::input("GLT count/label index mismatch"));
    }
    let stimuli = required(stimuli, "-num_stimts")?;
    for (index, entry) in stimuli.iter().enumerate() {
        if entry.source.is_none() || entry.label.is_none() {
            return Err(Error::input(format!(
                "stimulus {} requires one source and -stim_label",
                index + 1
            )));
        }
        match &entry.source {
            Some(PendingSource::Events { .. }) if entry.baseline => {
                return Err(Error::input(
                    "-stim_base is illegal with event timing, matching AFNI",
                ));
            }
            Some(PendingSource::Events { .. }) if timing_mode.is_none() => {
                return Err(Error::input(
                    "events require explicit -local_times or -global_times",
                ));
            }
            Some(PendingSource::Events { .. })
                if entry.minlag.is_some() || entry.maxlag.is_some() || entry.nptr.is_some() =>
            {
                return Err(Error::input("lag/sampling options require stim_file"));
            }
            _ => {}
        }
    }
    // Legacy lagged files default to dropping the initial lag interval in each
    // run; an explicit nfirst=0 intentionally retains those rows.
    let first = first.unwrap_or(
        stimuli
            .iter()
            .map(|s| s.maxlag.unwrap_or(0).div_ceil(s.nptr.unwrap_or(1).max(1)))
            .max()
            .unwrap_or(0),
    );
    Ok(Options {
        regcen,
        raw,
        milliseconds,
        subtract: subtract.unwrap_or(0.0),
        allow_empty,
        normalization,
        first,
        last,
        censortr,
        glts,
        glt_labels,
        acquisition,
        stimuli,
        ortvec,
        censor,
        primary: required(primary, "-x1D")?,
        full,
        manifest,
        timing_mode,
    })
}

/// Load validated inputs, build the numerical design, then publish requested files.
/// Options are consumed (moved) so the Acquisition can become part of Design
/// without copying its vectors. `as_deref()` borrows an optional slice/path.
fn execute(options: Options, command_line: &str) -> Result<()> {
    let nrows = options.acquisition.nrows();
    let mut sources = Vec::new();
    for entry in options.stimuli {
        let label = required(entry.label, "-stim_label")?;
        let source = match required(entry.source, "stimulus source")? {
            PendingSource::File(path) => {
                let nptr = entry.nptr.unwrap_or(1);
                let count = nrows
                    .checked_mul(nptr)
                    .ok_or_else(|| Error::input("sampled input size overflow"))?;
                Source::Nuisance(Nuisance::sampled(
                    &label,
                    read_table(&path, count)?,
                    nrows,
                    nptr,
                    entry.minlag.unwrap_or(0),
                    entry.maxlag.unwrap_or(0),
                    entry.baseline,
                )?)
            }
            PendingSource::Events {
                path,
                basis,
                mode,
                centers,
            } => {
                let format = match mode {
                    EventMode::Fsl => TimingFormat::Fsl {
                        unit_amplitude: false,
                    },
                    EventMode::Fs1 => TimingFormat::Fsl {
                        unit_amplitude: true,
                    },
                    EventMode::Ordinary => TimingFormat::Plain,
                    _ => TimingFormat::Married {
                        duration: basis.requires_duration(),
                    },
                };
                let timing = read_timing_with_options(
                    &path,
                    required(options.timing_mode, "timing mode")?,
                    &options.acquisition,
                    TimingOptions {
                        format,
                        milliseconds: options.milliseconds,
                        subtract: options.subtract,
                        allow_empty: options.allow_empty,
                    },
                )?;
                let mut stimulus = EventStimulus::new(&label, basis, timing, mode)?
                    .centered(options.raw, centers.as_deref())?;
                if let Some(target) = options.normalization {
                    stimulus = stimulus.normalized(target, options.acquisition.tr() as f32)?;
                }
                Source::Events(stimulus)
            }
        };
        sources.push(source);
    }
    for (path, label) in options.ortvec {
        sources.push(Source::Nuisance(Nuisance::new(
            &label,
            NuisanceKind::Ortvec,
            read_table(&path, nrows)?,
        )?));
    }
    // map invokes read_censor only for Some(path). The result would be
    // Option<Result<Vec<bool>>>; transpose makes it Result<Option<Vec<bool>>,
    // so ? propagates an I/O error without mistaking it for "no censor".
    let keep = options
        .censor
        .as_deref()
        .map(|p| read_censor(p, nrows))
        .transpose()?;
    // This local library call composes an original-timeline mask. It does not
    // shorten any input buffer; build_with_events retains the full matrix so
    // means, lag coordinates and event identities cannot depend on censoring.
    let keep = afni_xmat_maker::censor::combine(
        &options.acquisition,
        keep.unwrap_or_else(|| vec![true; nrows]),
        options.first,
        options.last,
        &options.censortr,
    )?;
    let mut design = build_with_events(options.acquisition, sources, Some(&keep))?;
    // GLTs are deliberately last: a stimulus may expand to many AM/IM/lag
    // columns, and resolving names earlier would use incomplete offsets.
    for (index, glt) in options.glts.into_iter().enumerate() {
        let label = options
            .glt_labels
            .iter()
            .find(|(i, _)| *i == index)
            .map(|(_, s)| s.clone())
            .unwrap_or_else(|| format!("GLT#{}_", index + 1));
        match glt {
            PendingGlt::Numeric(n, path) => design.add_numeric_contrast(&label, &path, n)?,
            PendingGlt::Symbolic(spec) => {
                let text = if spec.starts_with("SYM:") {
                    spec
                } else {
                    std::fs::read_to_string(&spec)
                        .map_err(|e| Error::resource(format!("cannot read GLT {spec}: {e}")))?
                };
                design.add_symbolic_contrast(&label, &text)?;
            }
        }
    }
    write_all_outputs(
        &design,
        &options.primary,
        options.full.as_deref(),
        options.manifest.as_deref(),
        options.regcen.as_deref(),
        command_line,
    )?;
    eprintln!(
        "++ wrote {} retained rows, {} full rows, {} columns to {}",
        design.retained_rows().len(),
        design.acquisition().nrows(),
        design.columns().len(),
        options.primary.display()
    );
    Ok(())
}

/// Quote arguments for human-readable provenance only; no shell is invoked.
fn command_line(args: &[String]) -> String {
    let mut words = vec!["1dXmatMaker".to_owned()];
    for arg in args {
        if !arg.is_empty()
            && arg
                .chars()
                .all(|c| c.is_ascii_alphanumeric() || "_./,:=-".contains(c))
        {
            words.push(arg.clone());
        } else {
            words.push(format!("'{}'", arg.replace('\'', "'\\''")));
        }
    }
    words.join(" ")
}

/// Process entry point. UTF-8 is required for portable labels/provenance.
pub(crate) fn run() -> Result<()> {
    let args = env::args_os()
        .skip(1)
        .map(|arg| {
            arg.into_string()
                .map_err(|_| Error::input("Command arguments and paths must be UTF-8"))
        })
        .collect::<Result<Vec<_>>>()?;
    if args.len() == 1 && matches!(args[0].as_str(), "-help" | "-h") {
        print!("{HELP}");
        return Ok(());
    }
    if args.len() == 1 && args[0] == "-version" {
        println!("1dXmatMaker {VERSION} (0.8 standalone release)");
        return Ok(());
    }
    if args.iter().any(|a| a == "-candidate_library") {
        return crate::candidate_cli::run(&args, &command_line(&args));
    }
    let options = parse(&args)?;
    execute(options, &command_line(&args))
}
