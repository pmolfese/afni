<!--
Developed at the NIH

Part of AFNI. See the repository LICENSE.txt for licensing terms and
disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
Third-party and adapted portions retain their existing attribution and license
terms; see THIRD_PARTY.md.
-->

# 1dXmatMaker

A standalone Rust program for generating AFNI design matrices without loading
images or running a model fit. Use it to prepare ordinary condition regressors,
single-trial designs, or alternative HRF designs for tools such as `3dPicard`,
`3dREMLfit`, and `3dLSS`.

Version 0.8.0 uses only the Rust standard library. Generation requires no AFNI
installation, external Rust crates, Python, afni-io, or afni-core. The executable
contains the GLMsingle waveform bank. AFNI is needed when fitting the matrices.

## Features

- Separate polynomial baselines for each run; sampled signals, lags, motion and
  other nuisance regressors.
- Ordinary, individual-event (IM), amplitude-modulated (AM1/AM2), FSL and
  duration-modulated timing, with the basis families listed below.
- Censoring, retained/full matrices, nuisance-spike companions, symbolic or
  numeric contrasts, and optional event-to-column identity tables.
- Twenty GLMsingle HRF candidates or custom sampled libraries, producing aligned
  trial and ONOFF matrices with the metadata needed by Picard.
- Picard beta identities export to MVPA/RSA through `1dXmatTool -beta_manifest`.

## Build and install

Run from this directory with **Rust 1.99.1** (the pinned toolchain) and a native
linker/SDK installed. There are no registry dependencies to download.

```sh
cargo build --release --frozen
./target/release/1dXmatMaker -help
cargo test --frozen
```

Make provides build, test and install wrappers:

```sh
make all check install PREFIX="$HOME/.local"
```

Alternatively, use standalone CMake:

```sh
cmake -S . -B target/cmake -DCMAKE_INSTALL_PREFIX="$HOME/.local"
cmake --build target/cmake
ctest --test-dir target/cmake --output-on-failure
cmake --install target/cmake
```

Both wrappers install `bin/1dXmatMaker` and
`share/1dXmatMaker/` containing this README, example inputs, attribution, licenses
and HRF provenance. Add the chosen `bin` directory to PATH. These wrappers build
independently of the parent AFNI build. Matrix generation requires no Python.

Builds and scientific checks have been exercised on native ARM64 macOS. Linux,
Intel macOS and older macOS versions still need platform qualification.
`--frozen` prevents dependency resolution and lockfile changes; the compiler
and system SDK must already be available. Everything under `target/` is generated
and can be deleted before a fresh build.

## Quick start

This creates two condition regressors with separate quadratic baselines for two
120-volume runs at TR = 2 seconds. Each timing-file row describes one run, with
onsets in seconds relative to that run's start.

```sh
printf '10 30 50 70 90 110 130 150 170 190\n12 32 52 72 92 112 132 152 172 192\n' > frequent.1D
printf '60 120 180\n62 122 182\n' > target.1D

1dXmatMaker -run_lengths 120,120 -TR 2 -polort 2 \
  -num_stimts 2 -local_times \
  -stim_times 1 frequent.1D GAM -stim_label 1 frequent \
  -stim_times 2 target.1D GAM -stim_label 2 target \
  -x1D frequent_target.xmat.1D

3dPicard -input 'run1.nii.gz run2.nii.gz' \
  -matrix frequent_target.xmat.1D -fraction 0.6 \
  -prefix frequent_target.picard.nii.gz \
  -fitts frequent_target.fitts.nii.gz \
  -errts frequent_target.errts.nii.gz
```

Supply two preprocessed image runs with those lengths and TR. The matrix has
240 rows and eight columns: six unpenalized baseline columns and two condition
columns. Picard uses the supplied baselines and estimates one coefficient per
condition shared across runs. Its quoted input list concatenates the images in
order using the AFNI reader (`AFNI_PATH_SPACES_OK` must not be `YES`). Input order
must match the timing rows and run lengths. Fraction 0.6 is illustrative, not tuned.

Use new output names: existing files are never overwritten, and parent directories
must exist. Run `1dXmatMaker -help` for the option reference and seven worked examples,
including contrasts, motion nuisance, IM censoring and amplitude/duration modulation.

## Timing, outputs and compatibility

Put acquisition, timing and normalization settings before the sources. Choose
`-local_times`, `-global_times`, or explicitly opt into `-guess_times`. Local timing
needs one noncomment row per run; `*` alone means no events in that run. Negative
or out-of-acquisition onsets are rejected. Timing units are seconds;
`-stim_times_millisec` converts onsets only, not durations.

Married timing uses `onset*amp1,amp2:duration`; duration is required for dm bases.
AM1 sums supplied amplitudes to weight a response. AM2 adds a mean response and
one response per centered amplitude; centering uses all events across all runs
of that stimulus before censoring. `-rawAM2` disables centering explicitly.
FSL input requires a dm basis; a `*` separator starts a new run.

`-polort 2` creates constant, linear and quadratic columns per run. `-polort A`
chooses order from the longest run duration and uses Legendre polynomials.
`-polort -1` removes polynomial baselines and disables baseline demeaning.
Mark sampled nuisance with `-stim_base`, or use `-ortvec FILE LABEL` for a table.
Event responses stop at run edges; historical sampled-source lags can cross them.

The complete design is generated and normalized before row selection. Censor
files contain one 0/1 value per original volume (1 keeps). Censoring, `-CENSORTR`,
`-nfirst` and `-nlast` intersect. CENSORTR run numbers are one-based, local/global
sample indices zero-based; commas separate independent tokens. With lagged
sampled sources, default nfirst removes unsupported initial rows; `-nfirst 0`
explicitly keeps them. Nlast is inclusive within each run.

| Output option | Contents |
|---|---|
| `-x1D FILE` | Required retained-row matrix with AFNI metadata and requested GLTs. |
| `-x1D_uncensored FILE` | Full original timeline, same columns, without GLT metadata. |
| `-x1D_regcensored FILE` | Full timeline plus one nuisance spike per excluded row. |
| `-event_manifest FILE` | Versioned TSV mapping events and basis components to columns. |

Keep the full original image time series when fitting a retained-row matrix;
its `GoodList` identifies included samples and `RunStart` preserves run boundaries.
Text inputs support comments, column/row selectors, transpose and `1D:` literals.
No implicit padding, truncation, stdin input, dataset loading or option abbreviation
is provided. AFNI environment settings and `.afnirc` do not affect generation.

This is a matrix generator, not a complete emulation of `3dDeconvolve`. It preserves
zero/redundant columns and does not impose an OLS rank test. Validate for the intended
consumer, for example `1dXmatTool -infile X.xmat.1D -validate -for reml`.
LSS requires one IM stimulus with one basis component per trial. A zero-support
trial may prevent run splitting because its run cannot be inferred from nonzero
matrix entries. Rank-deficient designs need appropriate fitting choices;
regularization does not make unidentified effects identifiable. Nuisance-spike
and row-censoring representations need not give the same REML autocorrelation fit.

Symbolic GLTs resolve after column ordering. `Ort` addresses polynomial columns;
use numeric GLTs for ortvec columns. Numeric contrasts must cover all columns,
including baseline zeros. Repeated symbolic terms assign rather than accumulate.
Matrices and companions are staged without replacement; ordinary failures roll
back new outputs, but publication is not a crash-atomic multi-file transaction.
Exit codes are 0 for success, 2 for invalid input, and 1 for I/O/resource failures.

## Basis functions and expressions

Quote parameterized bases in the shell, for example `'BLOCK(2,1)'` or
`'TENT(0,20,11)'`. Times are seconds; names are case-sensitive except expression
functions/variables. Invalid parameters and unknown families are errors.

| Family | Forms and parameter rules |
|---|---|
| GAM | Bare `GAM` uses p=8.6, q=0.547 and cutoff 11.1 s. `GAM(p,q)` requires p,q>0 and uses p*q+9.9*sqrt(p)*q. `GAM(p,q,d)` accepts d>=0, replaces nonpositive p/q with the reference defaults, and uses the reference duration support/table path. |
| GAMpw | Bare `GAMpw` aliases default GAM. `GAMpw(peak,FWHM[,d])` requires positive peak/FWHM; the AFNI peak/width search converts these to gamma parameters. |
| TWOGAM, TWOGAMpw | Five parameters: primary p,q, signed undershoot ratio, secondary p,q; pw substitutes peak/FWHM pairs. Both pairs must be positive and abs(ratio)<=1. Optional sixth parameter d>=0; d=0 retains the impulse path. Absolute peak normalization preserves signed lobes. |
| BLOCK, BLOCK4, BLOCK5 | `(d[,peak])`, d>0, peak>=0. Bare BLOCK means power four. Omitted/zero peak preserves natural integrated amplitude; positive peak sets the peak. |
| IGFUN, IGFUN4, IGFUN5 | Same fixed BLOCK waveform and peak convention, with explicit d<=999. |
| UBLOCK, UBLOCK4, UBLOCK5 | `(d[,peak])`, 0<d<=999, peak>=0. Zero/omitted peak normalizes to a long-duration reference; positive peak sets the response peak. |
| dmBLOCK[4/5], dmUBLOCK[4/5] | Optional `(peak)`; duration comes from married timing, in [0,999]. Zero duration gives zero. Negative dmBLOCK peak retains natural amplitude; negative dmUBLOCK peak means the reference duration whose amplitude is one (at most 999 s). |
| TENT, TENTzero | `(start,end,n)`, start<end, integer n>=2 (>=3 for zero). Negative start is allowed. Zero drops the first and last hats, not their sampling times. |
| CSPLIN, CSPLINzero | `(start,end,n)`, start<end, integer n>=4. Cardinal cubic interpolation with reference endpoint corrections; zero drops endpoint columns. |
| POLY | `(start,end,n)`, start<end, 1<=n<=20. Produces n columns of degree 0 through n-1, not n+1. |
| SIN | `(start,end,n)`, start<end, integer n>=1. Sine frequencies k*pi/(end-start), k=1..n. |
| TRIG | `(start,end,n)`, start<end, integer n>=3. Constant, alternating sine/cosine pairs using the reference ordering. |
| SPMG, SPMG1/2/3 | Optional `(d)`, d>=0. Bare SPMG aliases SPMG2. Components are canonical response, time derivative and dispersion finite difference; support is 0..25 s before duration convolution. |
| MION, MIONN | Optional `(d)`, d>=0. The reference exponential response and its sign-reversed counterpart, with 60 s impulse support. |
| WAV | Optional `(d,delay,rise,fall,undershoot,restore)`. Defaults 0,2,4,6,0.2,2. Omitted/empty fields retain defaults; all supplied values are nonnegative. Rise/fall/restore are clamped to at least 0.1 s as in AFNI. Undershoot is dimensionless. |
| EXPR | `EXPR(start,end) expression[~expression...]`, start<end. Each expression produces a component; details below. |

Duration-table families are limited to 999 seconds and bounded computational
work. Explicit zero-duration and omitted-duration forms can differ, especially
GAM and SPMG derivatives. BLOCK's stable integral can remain finite where the
historical C expression overflows; bitwise agreement with every C build is not promised.

For example, `'EXPR(-2,12) t~sin(PI*x)~z^2'` creates three components. Whitespace
or `~` separates expressions; use `t+x` without internal spaces for one component.
`t` is relative seconds, `x=(t-start)/(end-start)`, and `z=2*x-1`. Each component
must mention t/x/z, except literal `1`. Powers `^`/`**` associate to the right.
AFNI domain conventions apply: sqrt/log use absolute input; division by zero
returns zero; invalid real powers retain the base. Remaining nonfinite results fail.
Unimplemented AFNI functions AI, BI, I0, I1, K0, K1, DAI, DBI and ABSEXTREME are rejected.

Random EXPR uses a local reproducible stream, not AFNI's clock/PID-seeded sequence.
`-expr_seed N` accepts a u32 (default 305441741); source index s and component c use
`N + 1024*s + c` modulo 2^32. Each component uses a 48-bit drand48-style recurrence;
GRAN uses Box–Muller. Sampling follows event, component, then full-row order;
all function arguments are evaluated, including both IFELSE branches. Normalization
can consume random samples; censoring selects completed rows without redrawing.

## Sampled HRF candidates

From the package directory, or the installed `share/1dXmatMaker` directory, use a
new output prefix in a writable directory:

```sh
1dXmatMaker -candidate_library glmsingle -candidate_ids 0,10,19 \
  -events examples/candidate_events.tsv -run_lengths 144,144,144 \
  -TR 1 -polort 0 -candidate_prefix hrf
```

Omit `-candidate_ids` for all twenty shapes. Original IDs are preserved. This mode
requires 2–4096 candidates, integer polort and one fixed duration across events.
It accepts censoring and ortvec nuisance, but does not combine with indexed
stimulus options, GLTs, AM modulation or automatic polort.

Events use a literal-tab-separated header and rows:

```text
trial	condition	run	onset	duration	item
100	0	0	3	1	item0
101	1	0	23	1	item1
```

Trial IDs are unique; conditions can repeat. IDs are integers in 0..16777215;
run IDs are zero-based. Onsets are run-local seconds through the last acquired
sample; durations are 0..999 seconds. File order sets trial-column order. `item`
is explicit text or `-`. No identity is inferred from a label.

A custom library TSV has this exact header:

```text
id	waveform	dt	kind	duration	source	version
7	triangle.1D	1	impulse	0	custom-example	v1
3	box_triangle.1D	1	convolved	1	custom-example	v1
```

Waveform paths are relative to the library. Supply one finite numeric column,
2–1000000 samples starting at time zero, and positive dt in seconds. An `impulse`
with declared duration 0 is integrated against the event-duration box once;
a `convolved` waveform must declare the event duration and is not convolved again.
Source/version fields are required. Try `examples/custom_library.tsv` with
`-hrf_normalization none`; default `peak` requires a positive waveform maximum.

The embedded GLMsingle bank comes from commit
`082851f9db87ff30c49370f27365ac1112f18ff7`: twenty 501-sample curves at 0.1-second
spacing, representing a 0.1-second stimulus response. Rust reproduces its duration
box summation, PCHIP resampling at TR, and global float32 normalization.
`-hrf_normalization library` retains that scaling; default `peak` additionally
normalizes each curve's peak. `none` is invalid for this bank. All twenty curves
are prepared before selecting IDs. Fractional onsets use linear interpolation
of prepared curves, with zero outside support and clipping at run boundaries.
This does not import GLMsingle's image fitting, denoising or cross-validation.

Each candidate writes trial and collapsed ONOFF matrices, retained/full versions,
and its prepared HRF. Shared outputs include `.candidates.tsv`,
`.onoff.candidates.tsv`, column-identity tables, `.trial_info.1D` and `.events.tsv`.
ONOFF sums all trials across conditions; it has no single-trial identity.

Fit matching images with Picard, for example:

```sh
3dPicard -input signal.nii -candidate_manifest hrf.candidates.tsv \
  -candidate_policy hrf -fractions 1,.6 -trial_info hrf.trial_info.1D \
  -save_candidates -write_manifest -prefix fitted.nii
```

The example events require three 144-volume runs at TR 1. Repeated-condition tuning
needs sufficient repeated conditions and run-local nuisance structure. For fixed
shrinkage replace `-fractions ... -trial_info ...` with `-fraction .6`. Picard performs
selection; Maker does not automatically transfer an ONOFF winner into a trial fit.
Voxelwise selection has no single selected matrix: preserve the candidate catalog
and DesignID map with the betas.

## Downstream tables with 1dXmatTool

Picard's `-write_manifest` records the actual beta volumes. Use the updated
`1dXmatTool` to prepare downstream inputs; no separate Python helper is installed.
From Picard's original working directory:

```sh
1dXmatTool -beta_manifest fitted.betas.tsv -export_manifest \
  -subject s01 -input_base . -manifest_format manifest -prefix mvpa
3dMVPA -mode classify -trial_manifest mvpa.trials.1D \
  -target_column Condition -cv_column Run -mask mask.nii \
  -method ridge -prefix decoded

1dXmatTool -beta_manifest fitted.betas.tsv -export_manifest \
  -subject s01 -input_base . -manifest_format rsa_runwise -prefix rsa
3dRSA -mode RSA -runwiseTable rsa.runwise.1D -mask mask.nii \
  -model_mat category model.1D -classic_null conditions \
  -nperm 720 -seed 717 -prefix rsa_fit
```

Supply an appropriate condition model matrix; the example candidate events use
six conditions (0 through 5). Runwise RSA reads generated per-run Trial/Condition
tables and selectors into the original bucket. It averages repeated conditions
within each run and computes crossnobis distances; the export itself does no
averaging or image fitting. There is no inferred residual whitening.

For ordinary classic RSA, use `-manifest_format rsa`, optionally `-runs 1`, and
pass the resulting `PREFIX.rsa.1D` with `-dataTableFile`, `-condition_column Condition`
and an explicit `-condition_order`. This requires one selected beta per condition;
repeated cells fail. Both RSA exports write `PREFIX.identities.1D` to retain
candidate/column provenance that the primary RSA table does not accept.

Fixed-design sidecars may lack condition, run and trial identities. Supply
`-events FILE` with header `Volume Condition Run Trial` (optional Item/Basis).
Volume is the actual global beta-bucket index, not the source matrix column;
Run in this annotations file is one-based. One complete, unique annotation row
per beta is required. Supplied annotations cannot contradict known identities.
Do not manufacture run identities for condition coefficients shared across runs.

Picard manifest runs are zero-based; exported runs and `-runs` are one-based.
`SourceRun` preserves the input value, or -1 when absent. For a previous LSS
manifest from 1dXmatTool, explicitly use `-manifest_run_base 1`. Relative paths
resolve beside the manifest unless `-input_base` supplies Picard's original
working directory. Candidate catalogs and voxelwise DesignID maps are retained
without inventing a single selected Xmat. Table export does not establish run
independence or prevent HRF/parameter-selection leakage into downstream CV;
a predictive analysis must account for those fitting/selection steps.

## Source layout and testing

`src/` contains extensively commented Rust code. `examples/` contains small user
inputs. `tests/` contains the test runner, numerical probes and release packager.
`fixtures/afni/` groups reference cases by behavior; `fixtures/glmsingle/` holds
the independent HRF references. Keep these frozen scientific results for
regression tests; they are not needed at runtime and are not installed.
`data/glmsingle/` contains compile-time waveform data plus its license and provenance.
`licenses/` and [THIRD_PARTY.md](THIRD_PARTY.md) retain redistribution notices.

For optional scientific checks, see `tests/README.md`. They require NumPy/SciPy
and AFNI consumers; ordinary builds and Rust tests do not. New reports and all
build products belong under `target/`. Historical development reports are not
part of the distribution. See [AFNI licensing terms](licenses/AFNI.txt) and
[GLMsingle attribution](THIRD_PARTY.md#glmsingle-sampled-library).
