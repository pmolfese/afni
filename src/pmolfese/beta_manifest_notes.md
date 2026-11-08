<!--
Developed at the NIH.
Part of AFNI; see the repository LICENSE.txt for licensing terms and disclaimers.
-->

# Beta manifest workflows

Use the updated `1dXmatTool -beta_manifest` route to connect
`1dXmatMaker -> 3dPicard -> 1dXmatTool -> 3dMVPA/3dRSA`.
The former `picard_tables.py` helper is retired. Existing LSS matrix/events
exports remain available through `1dXmatTool -export_manifest`.

These commands were exercised on synthetic data with the Make and CMake builds
on native ARM64 macOS. Use the newly built tools on PATH and new output names.
The examples illustrate file/identity handling, not a prescribed analysis model.

## Inputs and design

Work in an analysis directory containing three preprocessed runs named
`run1.nii.gz`, `run2.nii.gz`, and `run3.nii.gz`. Each has **120 volumes, TR = 2
seconds**, and matching spatial geometry. Supply a matching `mask.nii.gz`.
The six conditions C0–C5 each occur once per run; six conditions meet the current
3dRSA minimum. Adjust both timing and annotations to match your real events.

The quoted Picard input list uses standard AFNI concatenation;
`AFNI_PATH_SPACES_OK` must not be `YES` for this form. Run order must match
the timing rows and `-run_lengths`.

The example uses a fixed GAM response and `-fraction 1`, the unregularized limit,
so it does not introduce data-dependent shrinkage across downstream validation
runs. Fractional-ridge and HRF-candidate fits use the same table export route;
any fitting, preprocessing or selection in a predictive analysis must respect
its held-out data boundaries.

```sh
# Timing rows correspond to runs 1, 2 and 3.
# Within each run the six events are conditions C0, C1, ..., C5.
cat > trials.1D <<'TIMES'
10 40 70 100 130 160
10 40 70 100 130 160
10 40 70 100 130 160
TIMES

# Volume is the fitted beta index; Trial is a unique event ID; Run is one-based.
cat > annotations.1D <<'EVENTS'
Volume Condition Run Trial
0  C0 1 0
1  C1 1 1
2  C2 1 2
3  C3 1 3
4  C4 1 4
5  C5 1 5
6  C0 2 6
7  C1 2 7
8  C2 2 8
9  C3 2 9
10 C4 2 10
11 C5 2 11
12 C0 3 12
13 C1 3 13
14 C2 3 14
15 C3 3 15
16 C4 3 16
17 C5 3 17
EVENTS

1dXmatMaker -run_lengths 120,120,120 -TR 2 -polort 1 \
  -num_stimts 1 -local_times \
  -stim_times_IM 1 trials.1D GAM -stim_label 1 trial \
  -x1D trials.xmat.1D

3dPicard -input 'run1.nii.gz run2.nii.gz run3.nii.gz' \
  -matrix trials.xmat.1D -fraction 1 -write_manifest \
  -prefix picard.nii.gz
```

The matrix has 360 rows and 24 columns: six run-local baseline columns and
18 trial columns. Picard's sidecar maps each actual beta Volume to its source
XColumn. They are different indices: the first trial is beta volume 0 but design
column 6. Fixed-design sidecars do not necessarily carry condition/run/trial
identities, so `annotations.1D` supplies them explicitly.

## MVPA

Export a whole-bucket trial manifest, then classify conditions with runs defining
the validation folds:

```sh
1dXmatTool -beta_manifest picard.betas.tsv -export_manifest \
  -events annotations.1D -subject s01 -input_base . \
  -manifest_format manifest -prefix mvpa_input

3dMVPA -mode classify -trial_manifest mvpa_input.trials.1D \
  -target_column Condition -cv_column Run \
  -mask mask.nii.gz -method ridge -prefix mvpa
```

`mvpa_input.trials.1D` preserves the beta-volume/design-column mapping and adds
the supplied subject and event identities. `-manifest_format selected` is also
available when a direct table with `InputFile[Volume]` selectors is needed.

## RSA model

This illustrative model treats the six conditions as ordered positions. Replace
it with your scientific hypothesis. Rows and columns are C0, C1, C2, C3, C4, C5.

```sh
# Replace this with the model appropriate to your conditions.
cat > model.1D <<'MODEL'
0 1 2 3 4 5
1 0 1 2 3 4
2 1 0 1 2 3
3 2 1 0 1 2
4 3 2 1 0 1
5 4 3 2 1 0
MODEL
```

## Runwise RSA

Use all three runs for unwhitened cross-validated distances:

```sh
1dXmatTool -beta_manifest picard.betas.tsv -export_manifest \
  -events annotations.1D -subject s01 -input_base . \
  -manifest_format rsa_runwise -prefix rsa_input

3dRSA -mode RSA -runwiseTable rsa_input.runwise.1D \
  -mask mask.nii.gz -model_mat ordinal model.1D \
  -classic_null conditions -nperm 720 -seed 717 \
  -prefix rsa_runwise
```

The export creates `rsa_input.runwise.1D`, per-run files such as
`rsa_input.r01.trials.1D`, and `rsa_input.identities.1D`. The runwise table contains
selectors into the original beta bucket; no beta images are copied. Each TrialFile
has `Trial Condition` rows in the exact order of its corresponding selector.

When conditions repeat within a run, 3dRSA averages their trial betas within that
run. 1dXmatTool does no averaging. Export requires explicit unique Trial IDs,
one coefficient per trial, at least two selected runs, and one bucket per run
(the same whole-session bucket can supply every run via different selectors).
Residual files and whitening are not inferred. Valid cross-run independence is
an analysis requirement, not something table formatting establishes.

## Classic RSA using run 1

Select one beta per condition and write the minimal classic-RSA table:

```sh
1dXmatTool -beta_manifest picard.betas.tsv -export_manifest \
  -events annotations.1D -subject s01 -input_base . \
  -runs 1 -manifest_format rsa -prefix rsa_classic_input

3dRSA -mode RSA -dataTableFile rsa_classic_input.rsa.1D \
  -condition_column Condition -condition_order C0,C1,C2,C3,C4,C5 \
  -mask mask.nii.gz -model_mat ordinal model.1D \
  -classic_null conditions -nperm 720 -seed 717 \
  -prefix rsa_classic
```

Classic export rejects duplicate conditions instead of silently averaging them.
It also writes `rsa_classic_input.identities.1D`, preserving metadata that the
primary RSA table cannot accept as varying subject covariates. Both RSA examples
use single-subject condition-label permutations, not population-level subject
inference. Keep the complete annotation table even when selecting a run: all
source mappings are checked before filtering.

## Identity and path conventions

| Field or option | Meaning |
|---|---|
| `Volume` | Zero-based index in the actual fitted beta bucket. |
| `XColumn` / candidate `Column` | Zero-based physical source design column; not a beta index. |
| Picard manifest `Run` | Zero-based original run; -1 denotes missing identity. |
| Exported `Run`, annotation `Run`, `-runs` | One-based original run; selection never renumbers runs. |
| `SourceRun` | Original input manifest's run value, or -1 when absent. |
| `Trial` | Explicit stable event identity, not a condition or matrix-column guess. |
| `-input_base .` | Resolve relative manifest paths from Picard's original working directory. |
| `-manifest_run_base 1` | Use when importing a prior 1dXmatTool LSS manifest instead of Picard's convention. |

Without `-input_base`, relative paths resolve beside the input manifest. Use the
explicit base when Picard's prefix contains a subdirectory. Exports use absolute
paths and refuse to overwrite existing files, including companion tables.

For candidate fits, `CandidateManifest`, `DesignIDFile`, `DesignIDVolume`, and
trial/basis identities are retained. A voxelwise-selected fit has no single
selected Xmat. Candidate sidecars with complete identities do not need the
supplemental `-events` file; supplied annotations must never contradict them.
The `-events` join uses complete unique Volume keys and requires one input bucket.

Implementation: [1dXmatTool, including beta export](../1dXmatTool.c), and
[focused end-to-end checks](../../tests/scripts/check_1dXmatTool_betas.py).
