# 1dXmatTool: design inspection, selection and run extraction

`1dXmatTool` is a C program for manipulating AFNI design matrices. Its
`-split_runs` operation extracts run-local matrices from a single multi-run
`3dDeconvolve` design. The implementation lives in `thd_xmat.[ch]`, so another
program can use the same extraction, mapping and writing functions.

This enables a single multi-run `afni_proc.py` invocation to supply matrices for
separate per-run `3dLSS` fits. `3dLSS` itself is unchanged.

## Usage

From the results directory containing the matrices produced by `afni_proc.py`:

```bash
mkdir -p split
1dXmatTool \
    -infile X.xmat.1D \
    -uncensored X.nocensor.xmat.1D \
    -split_runs \
    -prefix split/X

# Use the matching full-length, scaled EPI for each run. Adjust the processing
# block number and subject identifier to the actual afni_proc.py outputs.
3dLSS \
    -matrix split/X.r01.xmat.1D \
    -input pb03.SUBJ.r01.scale+orig \
    -prefix LSS.r01
```

`3dLSS` applies the local `GoodList`; supply every time point in that run, including
censored time points. `NRowFull` is the full run length, not the number of matrix
rows. Do not supply the entire concatenated multi-run dataset to a run matrix.

`-runs '1,3..5'` selects specific original run numbers. Run selection is one-based;
row, column, stimulus and volume indices in mapping tables are zero-based.
`-zero_tol 0` requires exact nonzero support; the default is `1e-6` relative to
each column's largest absolute value over the support matrix.

The output directory must already exist. Existing files, including symlinks,
are never replaced. All selected runs are extracted and checked before writing;
ordinary write failures remove files created by that invocation. A process
interruption can leave incomplete outputs, which must be removed or given a new
prefix before retrying.

## Outputs and provenance

For each selected original run, such as run 2:

| File | Contents |
| --- | --- |
| `PREFIX.r02.xmat.1D` | NIML-comment matrix readable by existing AFNI consumers, including `3dLSS`; local timeline and remapped stimulus metadata. |
| `PREFIX.r02.rows.1D` | Output row, input row, local full-timeline TR, original full-timeline TR. |
| `PREFIX.r02.columns.1D` | Every input column, its output column or `-1`, original/local stimulus and IM volume indices, and numeric `Kept` status; labels are comments. |
| `PREFIX.runs.1D` | One row per output run: original run number, original start TR, full length, retained row count, column count; paths are comments. |

`SourceRow` indexes rows in the input matrix, which may already be censored.
`SourceTR` indexes its full timeline. These are deliberately separate.
`SourceIMVolume` and `IMVolume` are indices within the associated IM stimulus.
For the single-IM-stimulus design required by `3dLSS`, `IMVolume` is the output
beta-volume index. Original column labels retain their original suffixes; use the
mapping, not those suffixes, to locate local LSS volumes.

Mapping comments escape backslashes, tabs, newlines and carriage returns as `\\`,
`\t`, `\n`, and `\r`. Column status `Kept` is 1 or 0; unavailable local
indices are `-1`.

The matrix header also stores `XmatSourceFile`, `XmatSourceRun`, `XmatSourceStart`,
`XmatSourceRows`, `XmatSourceColumns`, `XmatSourceCommandLine`, `XmatSupportFile`
and `XmatZeroTolerance`. These describe the **immediate input**, not a flattened
history of all previous extractions. Source/support paths are made absolute
when they resolve on disk. The original command is provenance; the
output `CommandLine` identifies the extraction operation.

## Extraction policy

- Requires explicit `RunStart`, `GoodList` and `NRowFull`.
- Retains row and column order and copies selected matrix values. It does not
  refit, normalize, demean or regenerate regressors.
- The existing `THD_xmat` reader uses floats. Written matrices preserve those
  float values exactly on reread; original double text precision beyond float
  precision is not retained.
- Drops columns without significant support in the selected run or its retained
  rows. Shared nuisance and ordinary stimulus regressors can remain in several
  outputs. A column's within-run peak must exceed `zero_tol * global_peak`.
- IM columns must have significant support in exactly one run. Zero or cross-run
  IM support is rejected rather than assigning a trial by guesswork. IM columns
  assigned to the selected run must retain support after censoring. A fully
  censored selected run is an error.
- `-uncensored` supplies the corresponding full-row design for support inference.
  Run boundaries, TR, labels, stimulus ranges/types and retained matrix values
  must agree. Without it, the tool can only infer support in the retained rows;
  censoring could hide an otherwise cross-run regressor.
- Rebuilds `ColumnLabels`, `ColumnGroups`, `GoodList`, `NRowFull`, `RunStart`,
  stimulus ranges and numbered `Basis*` attributes. Empty stimulus families
  disappear and surviving families are renumbered. Other attributes are retained.
- Drops GLT attributes and `ColumnMasks`, recording that removal in the header.
  Projecting contrasts into a smaller model requires a separate explicit policy.
- Does not guarantee matrix rank or downstream estimator suitability. `3dLSS`
  retains its requirements for one IM stimulus, multiple IM columns and enough
  observations. Constant or dependent retained nuisance regressors are not
  silently removed based on rank.

## Library API

- `THD_xmat_extract_run(source, uncensored, run, zero_tol, &map, err, nerr)`
  returns a new matrix. `run` is zero-based. `uncensored` may be NULL.
- `THD_xmat_run_map` holds output-to-input row/column maps, an inverse column
  map (`-1` for dropped columns), and original run/start/length information.
- `THD_xmat_write(matrix, path, err, nerr)` writes the complete matrix header and
  numeric data without replacing an existing path.
- `THD_xmat_write_run_map(source, result, map, rows_path, columns_path, err, nerr)`
  writes both provenance tables, including dropped source columns.
- Free results with `THD_xmat_free` and `THD_xmat_run_map_free`.

The extracted in-memory matrix has no file identity (`fname == NULL`). Writing
and rereading it establishes the new file identity. Row/column maps always refer
to the immediate source object. The writer regenerates core timeline/stimulus
attributes from the object; extraction is responsible for remapping indexed
basis and group metadata.

The CLI is registered in CMake and the classic Makefile build. The CMake binary
is `build/targets_built/1dXmatTool`.

## Verification, 2026-10-25

Automated checks:

```bash
cmake --build build --target 1dXmatTool thd_xmat_unit -j 8
ctest --test-dir build \
    -R '^(thd_xmat:unit|1dXmatTool:split-runs)$' --output-on-failure
```

The library unit test checks direct/inverse maps, float-exact write/read,
metadata, invalid timelines, collision handling and map-write rollback.
`tests/scripts/check_1dXmatTool.py` generates two unequal runs with boundary
censoring, shared/run-local nuisance regressors and IM trials. It compares the
split matrices with fresh single-run `3dDeconvolve` controls, compares existing
`3dLSS` estimators with NumPy pseudoinverse calculations, and fits full-length
synthetic EPI to verify application of local censoring. Its failure cases cover
missing metadata, malformed selections, mismatched reference designs,
fully censored runs, lost/ambiguous IM trials, GLT removal and overwrite safety.

Local Haxby verification uses the existing 12-run, 288-column design under
`haxby_mvpa/afni_proc_lss/preproc/sub-1.results`. It checks both the original
uncensored design and a copy with the 53 frame exclusions from the previously
verified demo. It does not rerun preprocessing or equate the two preprocessing
variants. Results:

- Twelve run matrices, each with 24 columns and eight trials; 1,399 retained
  rows over 1,452 full-timeline frames.
- All selected matrix values exactly preserved at float precision, with all 96
  trial identities mapped back to the original design.
- Existing `3dLSS` fits completed on each full-volume scaled EPI run.
- Across 864 sampled real-data trial betas, maximum absolute difference from
  independent NumPy LSS fits: `3.2701e-5`.
- Maximum difference between saved LSS estimator weights and NumPy: `5.2362e-7`.

Freshly regenerating Haxby designs separately per run produced small stimulus
regressor differences, reaching `8.5664e-4` at a response tail, and estimator
weight differences up to `1.757e-4`. Those regenerated controls are therefore
not numerically identical to extracting the original design. The exact
preservation claim concerns the original matrix; the independent fit check uses
that preserved matrix. Nuisance demeaning can also introduce constant offsets
between separately generated and joint designs.

Local evidence is in
`validation_20261025/xmat_split/real_validation_v3/verification.json`, with the
commands and outputs alongside it. The verification script is
`validation_20261025/xmat_split/verify_haxby.py`. Haxby checks use the existing
`isccomp` environment and no more than eight cores.

## Inspection and validation

```bash
1dXmatTool -infile X.xmat.1D -uncensored X.nocensor.xmat.1D \
    -describe -json design.json > design.1D
1dXmatTool -infile X.xmat.1D -uncensored X.nocensor.xmat.1D \
    -validate -for lss > lss_checks.1D
1dXmatTool -infile X.xmat.1D -validate -for reml
```

Reports include full/per-run dimensions, censor counts (full minus retained),
stimulus ranges, rank, residual degrees of freedom, condition number, zero
columns and proportional/sign-flipped duplicate pairs. Rank and condition use
unit-norm columns and SVD, with relative cutoff `-rank_tol 1e-7`; zero columns
remain zero. A deficient matrix has infinite condition (JSON `null`). Duplicate
pairs satisfy `1-abs(cosine) <= rank_tol`. Description prints diagnostics without
failing solely because a matrix is deficient. Validation returns 2 for failed
checks, 1 for invalid arguments/read/write errors, and 0 on success.

Default validation checks full/per-run rank and residual DF. LSS preflight
checks every run's collapsed trial model (nuisance + other trials + target trial),
requires one IM family, at least two trials, and the existing 3dLSS row-count
restriction. It also checks trial support/censoring through run extraction.
REML preflight checks full/per-run rank, residual DF and the minimum nine rows.
These are design checks; data, geometry and runtime checks remain with the fitter.

## Select columns, stimuli or groups

```bash
1dXmatTool -infile X.xmat.1D -select_columns '0..3,8' -prefix subset
1dXmatTool -infile X.xmat.1D -select_labels 'motion_demean[0],motion_demean[1]' \
    -prefix motion
1dXmatTool -infile X.xmat.1D -select_stimuli blocks \
    -select_groups baseline,nuisance -prefix selected
```

Indices are zero-based. Names are exact, comma-separated labels. Groups accept
`baseline` (-1), `nuisance` (0), `stimuli` (positive groups), or a numeric group ID.
Selection options combine by union and keep source order. Outputs are
`PREFIX.xmat.1D`, `PREFIX.rows.1D` and `PREFIX.columns.1D`; all input rows and
censor/run metadata remain. Stimulus/basis ranges and groups are remapped.
If GLT attributes exist, selection requires explicit `-drop_glts`; it does not
project contrasts. ColumnMasks are dropped and the removal recorded.

## Compare designs

```bash
1dXmatTool -infile X.xmat.1D -compare other.xmat.1D \
    -atol 1e-6 -rtol 1e-6 > comparison.1D
```

Unique column labels identify matching regressors; full-timeline GoodList indices
identify matching rows. Reports distinguish column order from changed values,
column membership, censoring, run boundaries/TR, stimulus label/option/formula,
and column-group changes. Reordering alone is equivalent. Values differ when
`abs(a-b) > atol + rtol*max(abs(a),abs(b))`. Exit codes are 0 for equivalence,
2 for differences, and 1 for errors. This compares design semantics listed above,
not every arbitrary header or provenance attribute. GoodList must be increasing;
an arbitrary row permutation is invalid input rather than an equivalent design.

## Export trial manifests

Supply conditions explicitly; the tool never infers categories from IM labels.
An events file has one unique, complete row per run-local LSS volume:

```text
# Columns: Run Volume Condition Item
1 0 face block01
1 1 house block02
```

`Item` is optional; absent values use original column labels. Run is one-based,
Volume is zero-based within each run, and both refer to the immediate input.
Use the original multi-run design when exporting all run-specific buckets.

```bash
1dXmatTool -infile X.xmat.1D -uncensored X.nocensor.xmat.1D \
    -export_manifest -events events.1D -subject sub-1 \
    -beta_pattern 'lss/sub-1_run-{run}.LSS.nii.gz' -prefix samples
3dMVPA -trial_manifest samples.trials.1D ...
```

Exactly one `{run}` becomes `01`, `02`, etc. Single-run inputs may use `-beta FILE`
instead. Files need not exist at export time. The default manifest contains
Sample, Subj, Run, Item, Condition, Volume, Xmat, XColumn and InputFile; XColumn
maps to the original source column and InputFile names the entire beta bucket.
Use `-manifest_format selected` for a direct `-dataTableFile` input to 3dMVPA,
3dRSA or 3dSRM: InputFile includes `[Volume]`, and provenance columns are named
SourceXmat/SourceColumn to avoid triggering implicit Xmat item expansion.
There must be exactly one IM family; missing, duplicate or extra event keys fail.

## AFNI 1D format and shared APIs

Run/row/column maps contain numeric data separated by spaces. `# Columns:` names
fields; `# ColumnLabel` and `# RunFiles` hold textual metadata. They can be read
by AFNI numeric tools such as `1dcat`. AFNI's numeric reader in `mri_read.c`
accepts both spaces and tabs and ignores `#` comment lines.

Manifests and 3dMVPA result/audit tables also use `.1D` with a `# Columns:` schema.
They contain text labels as well as numbers and use the shared `THD_datatable`
reader, **not** the numeric-only reader. Text containing spaces is quoted;
embedded quotes are doubled. Plain un-commented table headers remain supported.
3dMVPA's automatic table filenames now end in `.1D`; explicitly requested
`-save_*` filenames are honored. Its numeric/spatial dataset formats are unchanged.

The library additionally exposes `THD_xmat_select_columns` and
`THD_xmat_diagnostics`; the latter reports SVD-based diagnostics without modifying
the design. The new feature integration test independently checks these results
against NumPy SVD, exercises selection and comparison failures, reads numeric
maps with 1dcat, and feeds both manifest styles into 3dMVPA. Run it with:

```bash
ctest --test-dir build -R '^1dXmatTool:features$' --output-on-failure
```

The new operations were also exercised on the existing Haxby data. All twelve
24-column run designs passed LSS preflight, with 1,399 retained rows across the
1,452-frame timeline. All 96 exported condition/volume/source-column mappings
matched the previous manifest exactly. Both face/house and eight-category ridge
fits completed and their commented `.1D` prediction tables were read back.
Evidence: `validation_20261025/xmat_features_haxby/verification.json`. Existing
preprocessing and LSS betas were reused; this was not a raw-data rerun.
