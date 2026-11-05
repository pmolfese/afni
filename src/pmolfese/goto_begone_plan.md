# Goto removal checklist and plan

- [ ] `src/thd_srm.c` — 146 occurrences
- [ ] `src/thd_mvpa.c` — 58 occurrences
- [ ] `src/thd_crossval.c` — 29 occurrences
- [ ] `src/thd_patterns.c` — 17 occurrences
- [ ] `src/thd_permute.c` — 14 occurrences
- [ ] `src/thd_rsa.c` — 11 occurrences
- [ ] `src/thd_simmatrix.c` — 10 occurrences
- [x] ~~`src/thd_xmat.c`~~ — 9 → 0 occurrences (2026-11-05)
- [ ] `src/thd_sampleset.c` — 5 occurrences
- [ ] `src/thd_datatable.c` — 3 occurrences
- [x] ~~`src/1dXmatTool.c`~~ — 32 → 0 occurrences (2026-11-05)
- [ ] `src/3dRSA.c` — 23 occurrences
- [x] ~~`src/SUMA/SUMA_SurfClustSim_core.c`~~ — 26 → 0 occurrences (2026-11-05)
- [x] ~~`src/SUMA/SUMA_SurfClustSim.c`~~ — 12 → 0 occurrences (2026-11-05)
- [x] ~~`src/SUMA/SUMA_SurfFWHM.c`~~ — 6 → 0 occurrences (2026-11-05)
- [x] ~~`src/SUMA/SUMA_SurfACF.h`~~ — 2 → 0 occurrences (2026-11-05)
- [ ] `src/pmolfese/tests/test_patterns.c` — 18 occurrences
- [ ] `src/pmolfese/tests/test_permute.c` — 14 occurrences
- [ ] `src/pmolfese/tests/test_rsa.c` — 5 occurrences
- [ ] `src/pmolfese/tests/test_crossval.c` — 1 occurrences
- [ ] `src/pmolfese/tests/test_mvpa.c` — 1 occurrences
- [ ] `src/pmolfese/tests/test_sampleset.c` — 1 occurrences

Cross out each file and mark its checkbox when its jumps have been removed and
focused validation is complete. Keep original counts as the audit baseline.

## Audit scope

Audited current files at commit `a9358ec31` on 2026-11-05, using Peter Molfese's
Git authorship history and source credits. Found 443 written goto occurrences
in 22 production/test files: 302 in 10 thd_ files, 101 in other production files,
and 40 in tests. Comments and strings are excluded. Macro definitions count once,
not once per expansion. SUMA_SurfFWHM.c is an existing file Peter modified; its
six jump lines are attributed to his commits. Git attribution does not establish
which AI agent generated individual statements.

No jumps in current 3dMVPA.c, 3dInSync.c, 3dSRM.c, thd_insync.c,
thd_insync_stats.c, thd_mapinfer.c, thd_phasefft.c, or thd_progress.c.
This inventory concerns the current checkout, not files only on other branches.

## Refactoring approach

Preserve scientific behavior, CLI defaults, error messages/status, numerical
operation order, random-number consumption, metadata, and ownership transfer.
Use if/else, early returns, ordinary loop termination and meaningful helpers.
Do not replace jumps with deeply nested conditionals or artificial one-pass loops.
For multi-resource operations, let an outer function own and release resources
while inner operations return status. Reuse existing destructors. Do not duplicate
long cleanup sequences at every failure. No source AI-attribution comments are
requested; PR attribution remains required by AGENTS.md.

1. Replace ordinary branches in 3dRSA.c: the two have_datasets jumps become
   runset/input-RDM/ordinary-input branches; subject analysis and general analysis
   become explicit alternatives followed by common cleanup. In thd_simmatrix.c,
   condition permutation work on whether a permutation set exists.
2. Refactor small constructors/parsers in thd_sampleset.c, thd_datatable.c and
   thd_xmat.c. Validate before allocation where practical and release partially
   constructed objects with existing destructors.
3. Split shared workflows into meaningful stages. Extract one X-matrix row's
   work in thd_patterns.c; separate validation, construction and processing in
   thd_crossval.c/thd_permute.c. Preserve statistical failure outputs in
   thd_rsa.c/thd_simmatrix.c.
4. Handle MVPA and SRM in focused batches. For MVPA, separate paired-prediction
   validation, scoring, allocation and serialization; remove jumps hidden in
   macros. For SRM, separate parsing from fitting and extract preparation,
   initialization, iteration and result transfer with explicit status returns.
5. Finish 1dXmatTool and SUMA routines, then refactor test cleanup independently.

User-selected first batch: thd_xmat.c and 1dXmatTool.c.

## Validation requirements

- Run relevant existing numerical and integration tests against verified paths
  to the new build; compare unchanged baseline behavior.
- Cover malformed inputs, partial output rollback, existing output preservation,
  invalid metadata, empty/nonfinite/rank-deficient data where applicable, and
  partial allocation failure where practical.
- Build affected targets with both platform Make and separate-directory CMake
  when supported. Record any unavailable build or dependency, without claiming
  parity. Preserve the local Makefile and existing build products.
- Check affected OpenMP behavior where relevant, scan for remaining jumps,
  review the diff, and run git diff --check.

## Completed work and validation

### 2026-11-05: thd_xmat.c and 1dXmatTool.c

Removed all 41 written jumps in this batch, including error macros. The baseline
inventory had 402 occurrences remaining in 20 files after this batch.

- thd_xmat.c: subset construction and run extraction return status to wrappers
  that own scratch storage and partial output objects. Column selection uses
  guarded construction and existing destructors. Numerical loops and provenance
  semantics are preserved.
- 1dXmatTool.c: separate strict parsing, option validation, column/group selection,
  event validation, per-run manifest writing and split output ownership. Early
  returns from processing helpers reach one owning cleanup path. Split rollback
  retains its count of completed outputs and ownership of the summary file.
- The column-selection parser now always frees its temporary translated list,
  including invalid later tokens. Newly touched run-list and event scratch
  allocations have explicit failure checks.
- Added regression coverage to tests/scripts/check_1dXmatTool_features.py for
  existing row/column-map collisions, invalid group selections, later-run missing
  events, unmatched events and invalid beta patterns. Partial outputs disappear;
  pre-existing maps remain intact.

Validation completed on macOS arm64 with the local GCC 14 Make configuration
and the existing separate CMake build at /tmp/afni-pr1030-cmake:

- `make -C src 1dXmatTool`: passed.
- `cmake --build /tmp/afni-pr1030-cmake --target 1dXmatTool thd_xmat_unit --parallel 4`:
  passed. No source/object list changes were needed because no source files or
  targets were added.
- `ctest --test-dir /tmp/afni-pr1030-cmake -R '^thd_xmat:unit$' --output-on-failure`:
  passed (1/1).
- check_1dXmatTool.py and check_1dXmatTool_features.py: passed with both the Make
  executable at src/1dXmatTool and the newly built CMake executable. The maximum
  synthetic LSS beta error was 7.09e-07. Supporting 3dDeconvolve, 3dLSS and 1dcat
  executables were supplied explicitly from src/macos_13_ARM; MVPA ingestion
  used src/3dMVPA. These supporting executables were not rebuilt for this batch.
- Compared 80 CLI calls from those suites against an executable built before
  editing: identical exit codes, stdout, stderr and output-file bytes, including
  the added failure cases. The comparison restored identical input/output state
  before each candidate run.
- Source scan: no goto tokens remain in either target file. Diff reviewed;
  `git diff --check` passed.

Logs and the temporary baseline comparison harness are in
/tmp/goto-xmat-EvRWGC. Tests used /tmp/pr1030-codespell/bin/python, which provides
NumPy and nibabel. The first test invocation used absent support executables in
src/ and was rerun with the explicit existing paths above. Builds retain existing
compiler warnings. Linux and exhaustive allocation-failure injection were not
run. No binaries, build trees or temporary test outputs were added to Git.


### 2026-11-05: SUMA simulation and surface ACF routines

Removed all 46 jumps from SUMA_SurfClustSim_core.c, SUMA_SurfClustSim.c,
SUMA_SurfFWHM.c and SUMA_SurfACF.h. Across the original inventory, 87 have now
been removed, leaving 356 in 16 files.

- Graph construction, adaptive smoothing and ACF calibration now return status
  to functions that own scratch and partial results. Calibration still resets
  its file-static optimizer state after success or failure; this refactor does
  not change its serial-only behavior.
- Surface-mask loading, volume-to-surface field construction and multi-threshold
  table output use guarded stages or borrowed-resource processing helpers.
- SUMA_SurfFWHM.c changes are confined to SUMA_SurfFWHM_ACF: 16 added/15 removed
  lines. The existing ordinary/global/local FWHM estimator is untouched.
- SUMA_SurfACF.h commits an accumulated curve only on success, as before.
  The estimator's formulas, binning, sample order, variance normalization and
  fitting calls remain the same. Simulation seeds, smoothing order, convergence
  criteria and numerical tolerances are unchanged.
- Added SUMA_SurfACF.h prerequisites to the SurfFWHM/SurfClustSim targets in
  SUMA_Makefile_NoDev (included by the platform SUMA_Makefile). CMake tracks
  these included headers automatically; no new production source/target exists.

Numerical verification uses baseline binaries built from the unchanged files,
with each candidate compared to its own build-system/compiler baseline. The
standalone test is tests/scripts/check_surf_goto_equivalence.py. It needs NumPy,
nibabel, explicit executable paths and the test-only seed library compiled from
tests/scripts/surf_fixed_seed.c. The usual pytest environment was unavailable
(pytest, datalad and xvfbwrapper are absent from the available test interpreter),
so this task used that standalone synthetic-data regression instead.

Important baseline finding: repeated runs of the unchanged SurfFWHM binary give
slightly different fitted ACF parameters. powell_int.c seeds its randomized
starts with time()+getpid(), overriding AFNI_RANDOM_SEEDVAL. The empirical curves
remain identical. A test-only srand48 interposer gives baseline and candidate
identical optimizer starts; it is NOT linked into or enabled for production.
No tolerance was loosened and no production random behavior was changed.

Completed checks:

- Make: `make -C src suma_exec 'SUMA_MAKE_LIST=SurfClustSim SurfFWHM'` passed
  with the existing local GCC 14 / macOS arm64 configuration.
- CMake: configured a separate tree at /tmp/goto-suma-41bmOa/cmake and built
  SurfClustSim and SurfFWHM with GCC 14. Used already available NIfTI/GIFTI
  sources and the locally built GTS library to avoid network fetching. Enabled
  REMOVE_BUILD_PARITY_CHECKS only in this temporary configuration because the
  existing 1dXmatTool installation-component entry is missing; no shared build
  settings or packaging files were changed. Existing compiler warnings remain.
- Make baseline comparisons: all 19 cases passed with exact output-file bytes,
  stdout, exit status and generated ACF parameter reports. Cases include fixed
  smoothing, all four sidedness modes, radii -1/-2/30, masked runs at 1/4 threads,
  adaptive block sizes 2/12/60 (including both detrending reference types),
  volume-to-surface mapping, multi-threshold calibration, ACF synthesis, ordinary
  FWHM, masked/detrended ACF, constant/mixed data and too few ACF bins.
  Raw masked maxima also match across thread counts. The existing ACF synthesis
  acceptance bound (abs(generated_a - requested_a) < 0.15) is preserved as a
  separate scientific check; baseline comparisons themselves use no tolerance.
- Temporary instrumented CMake binaries captured full-precision empirical ACF
  doubles and fitted floats: byte-identical for ordinary, masked and mixed data;
  constant data took the same no-fit path. This goes beyond rounded text output.
- CMake baseline comparisons: all 19 cases also passed with exact equality.
- Sparse NIML datasets, with and without an additional mask: exact stdout and
  ACF curve files under both Make and CMake (four additional comparisons).
- No goto tokens remain in the four requested files. Diff reviewed and
  git diff --check passed. Linux and exhaustive allocation-failure injection
  were not tested. No generated binaries or build artifacts were added to Git.

Validation logs, baseline binaries and the raw-array/sparse-data harnesses are
under /tmp/goto-suma-41bmOa. Baseline CMake objects were compiled from untouched
HEAD copies using the target's exact generated compile/link commands and the
same built libraries. The test-only seed library was built with:

```sh
cc -dynamiclib tests/scripts/surf_fixed_seed.c -o /tmp/goto-suma-41bmOa/surf_fixed_seed.dylib
```

Pass that path with --seed-library, together with --sim, --fwhm, --baseline-sim,
--baseline-fwhm, and --ico, to check_surf_goto_equivalence.py. On Linux the helper
uses LD_PRELOAD and can be compiled with `cc -shared -fPIC` instead.
