# 3dPicard: voxelwise penalized regression and the 1dXmatMaker roadmap

Created 2026-10-29; library-boundary review updated 2026-10-30.
Status: design proposal; no implementation authorized by this plan.
The program name is settled: **3dPicard**, a reference to Jean-Luc
Picard and a tribute to Patrick Stewart. Library names should describe their
numerical responsibilities rather than require a PICARD acronym.

**Development order: Picard first using existing 3dDeconvolve matrices, then
1dXmatMaker.** Establish and validate Picard's matrix input, nuisance handling,
fixed-design fractional ridge and an initial tuning workflow before making a
new design generator part of its development path. Existing `.xmat.1D` files
remain supported after Maker arrives; Maker is an additional producer, not a
prerequisite for Picard or a required replacement in existing AFNI workflows.

The [Picard milestones](#suggested-implementation-sequence-and-acceptance-gates)
come first below. The [1dXmatMaker roadmap](#1dxmatmaker-independent-design-generation-roadmap)
defines the independent generator and its later integration. All milestones
are planned; naming and documenting them does not authorize implementation.

The user reports the 3dMVPA libraries complete, and its roadmap records elastic
net (#28), banded ridge (#35), PLS (#30), and SVR (#31) as completed. The
2026-10-30 source review below identifies concrete extraction boundaries in
`thd_mvpa.[ch]`; completion of those estimators does not mean they have already
been split into standalone libraries. This update changes the plan only and
does not rerun or recertify the MVPA validation results.

## Purpose and boundaries

**Hard constraint: do not modify or refactor 3dDeconvolve.** The new
**1dXmatMaker** tool will independently produce standard `.xmat.1D` files
for **3dLSS, 3dREMLfit, and 3dPicard**. 3dDeconvolve remains an untouched
reference for compatibility tests; it will not be migrated onto new libraries
or made dependent on this tool. This supersedes the earlier conversational
suggestion of extracting a library used by both executables.

One command, provisionally `-penalty ridge|lasso|enet`, for Gaussian linear
models with unpenalized nuisance terms. Primary applications are single-trial
beta estimation and voxelwise encoding. Fractional ridge is the first numerical
deliverable; lasso/elastic net and banded ridge follow through shared solvers.

The distinctive addition is the combination of fractional ridge, time-series
design/censor/run handling, nuisance-aware estimation, and appropriate tuning.
AFNI already has lasso in `3dTfitter`/`thd_lasso.c`, and MVPA has ordinary ridge
and completed continuous-estimator support. Avoid creating a third independent
implementation of the same elastic-net mathematics.

A future single-trial orchestrator can add an HRF library, noise-pool PCs,
PC-count selection, and staged outputs. Its name is undecided. Fractional ridge
alone is not GLMsingle equivalence. Classification, general non-Gaussian GLMs,
and conventional penalized-coefficient t/p maps are outside the initial scope.

## Model and numerical contract

Write the model as `Y = Z Gamma + X B + E`: rows are retained observations,
X contains penalized predictors, Z contains explicit unpenalized columns, and
columns of Y are voxel responses. Define intercept handling explicitly; a
low-level fit must not silently add or recenter an intercept already handled by Z.

For one response, use the documented elastic-net objective:

`||y - Z gamma - X beta||^2/(2n) + lambda * [rho*||beta||_1 + (1-rho)*||beta||^2/2]`.

- `rho=1` is lasso; `rho=0` is ridge. Prefer `lambda` and `l1_ratio`/`rho` in
  shared APIs to avoid incompatible meanings of alpha between packages.
- The original proposal's `-alpha` mixing option remains a CLI decision. The
  inspected MVPA code calls overall strength alpha and mixing rho. Translate
  explicitly at interfaces; do not change existing MVPA option semantics.
- Existing `THD_ridge_fit` adds its lambda directly to X'X. For the normalized
  objective above, its diagonal penalty is `n * lambda` at rho=0. Record and
  test conversions; a shared function name does not imply shared units.
- For rho>0, lambda_max is `max(abs(X' y))/(n*rho)` in the specified projected
  and scaled fitting space. Pure ridge has no finite universal all-zero cutoff.
- Fractional ridge uses `||beta_ridge||_2 / ||beta_reference||_2` in a declared
  coefficient space. Use the minimum-norm least-squares reference when rank
  deficient. Specify rank tolerance, the zero-reference case, fraction=1
  (unregularized reference), and fraction=0 (zero penalized coefficients).
- Default candidate policy: preserve supplied design scaling for fractional
  ridge; offer explicit predictor standardization for sparse regression, with
  training-derived scales and coefficients back-transformed to original units.
  Freeze the exact centering/projection/scaling order before implementation.
- Nuisance terms receive no penalty. Internal QR/SVD normalization of nuisance
  columns may improve conditioning without changing their fitted span; this
  is different from changing which coefficients are penalized.
- Compute factorizations, accumulation, coefficients, and convergence checks
  in double precision. Float AFNI inputs/outputs are acceptable with documented
  conversion. The current xmat reader stores floats; converting those values
  to double does not restore precision already lost when parsing.
- Coordinate descent needs convergence diagnostics and full KKT verification,
  including rechecking any screened variables. Warm starts and Gram caching
  are optimizations after the basic numerical contract is validated.
- A support-only OLS refit is a later, explicitly named option. Define rank-
  deficient support behavior and tune/evaluate the resulting estimator itself.
  It is not automatically the complete relaxed-lasso procedure used by glmnet.

For a fixed fit, project both X and Y against the nuisance span, fit the
penalized coefficients, then recover nuisance coefficients from Y-XB. Keep
the information needed to reconstruct predictions in original units. With
whitening, transform Z, X and Y consistently before nuisance projection.

## Input, validation, and output contract

Provisional inputs: `-input`, `-mask`, `-matrix`, explicit penalty/nuisance
column selection, optional `-ortvec`/`-polort`, and run/trial metadata. Input
datasets must share the spatial domain and match the full time axis. Matrix
GoodList, NRowFull and RunStart map retained rows to original runs/time points.
Reject ambiguous censoring or duplicated baseline specifications.

Read ColumnGroups and stimulus metadata through thd_xmat, with explicit role
overrides. Penalty groups are a separate concept from AFNI stimulus groups,
run groups, and CV groups. Do not infer one from another without a documented
mapping. Add a validate/show-design mode before expensive fitting.

There are two distinct tuning problems:

1. **Shared-predictor prediction:** learn feature coefficients on training runs
   and predict held-out runs. Support run-wise folds, then contiguous blocks
   with temporal exclusion gaps, explicit fold manifests, and nested evaluation
   when an independently evaluated predictive score is requested.
2. **Single-trial beta estimation:** a regressor unique to a held-out run has no
   coefficient estimable from other runs. Ordinary leave-one-run-out prediction
   on a concatenated trial-specific design is not valid tuning for those betas.
   Specify a repeated-condition/trial-identity tuning procedure separately,
   informed by GLMsingle. Without sufficient repeats, allow fixed fractions
   rather than inventing an apparently cross-validated selection.

Fit scaling and learned preprocessing within the appropriate training split.
Specify how run-specific nuisance columns are scored on held-out runs: fitting
held-out nuisance coefficients from held-out Y changes the estimand and must
not be silently called raw out-of-sample prediction. Likewise, distinguish
selection loss from outer-test performance. Define the R2 baseline, pooling,
undefined-score behavior, and handling of negative R2 explicitly.

One-SE selection means stronger regularization within the declared uncertainty
rule, not always a sparser model. Preserve fold-level losses needed for the
rule and define weighting for unequal folds. Global selection needs an explicit
voxel mask, response normalization, and aggregation policy; it is not guaranteed
to be faster merely because it produces one selected parameter.

Provisional outputs:

- Labeled penalized betas; optional nuisance betas, fitted series and residuals.
- Selected fraction and actual lambda for ridge; lambda and mixing parameter
  for elastic net; candidate index and group penalties when relevant.
- Clearly labeled tuning scores and, when available, separate outer-test scores.
- Nonzero counts where meaningful and convergence/rank/status diagnostics.
- Optional held-out predictions with fold/row mapping; censored rows must not
  appear to be predictions. Prediction beyond retained rows requires a matching
  full design, not values invented from the censored matrix.
- Provenance with labels, trial/run identities, folds, candidate grids, scaling,
  solver tolerances and versions. Full coefficient paths only on request for
  limited subsets. Emit a beta/trial manifest usable by RSA and MVPA adapters.

## Concrete library boundaries

### Source baseline and dependency direction

The completed numerical implementations are still in `thd_mvpa.c`, with public
interfaces in `thd_mvpa.h`. No `thd_regression`, `thd_linproj`, or `thd_picard`
module exists at this review. These are planned boundaries, not existing APIs.
Use one regression library for the estimator family, not a library per penalty.

```text
3dMVPA.c -> thd_mvpa -> thd_regression
                     -> thd_crossval / existing pattern and inference modules

3dPicard.c -> thd_picard -> thd_regression
                       -> thd_linproj
                       -> thd_crossval
           -> thd_xmat / AFNI dataset and output APIs

1dXmatMaker -> independent basis/design generation -> thd_xmat -> .xmat.1D
                                                                    |
                                       3dLSS / 3dREMLfit / 3dPicard <-+
```

Arrows express direct responsibilities; existing transitive AFNI header/link
requirements are not a mandate for a broad dependency cleanup. Picard must not
call MVPA's CV runner or depend on classification, permutation, saved-model or
searchlight workflow code. RSA remains a downstream consumer of beta patterns.
3dDeconvolve is outside this graph and stays untouched.

| File/module | Owns | Does not own | Planned action |
|---|---|---|---|
| `thd_regression.[ch]` | Ridge and continuous-estimator numerical fits; coefficients, fitted numeric state, diagnostics; new prepared-design/path fits. | Dataset I/O, labels, censor/run interpretation, fold construction, candidate selection, inference. | New module: move the completed `THD_ridge_*` and `THD_regression_*` blocks together, then add Picard numerical capabilities. |
| `thd_linproj.[ch]` | Factorization of supplied unpenalized Z; rank policy, projection of multiple RHS, nuisance coefficient recovery. | Automatic intercepts, column-role guesses, fold selection, scaling policy, whitening estimation. | New small numerical module for Picard. MVPA nuisance migration is a separate optional change, not a prerequisite. |
| `thd_picard.[ch]` | Fitting workflow over numeric designs/responses: roles and maps supplied by caller, training-only preparation, scaling/back-transforms, scores, response-specific selection and final refits. | AFNI dataset handles, command parsing, file writing, HRF generation, RSA statistics. | New module; thin application layer above shared solvers and folds. |
| `thd_crossval.[ch]` | Index-based fold plans, nesting, validation, generic search traversal. | HRF/fraction interpretation, nuisance fits, repeated-condition beta score definitions. | Reuse; add metadata-only fold construction and a separate path/batch traversal when needed. Preserve current scalar-search API. |
| `thd_xmat.[ch]` | Matrix parsing/writing, labels/attributes, retained/full row mappings, run/column transformations, metadata adapters. | Penalized fitting, response-specific preprocessing, HRF candidate selection. | Reuse existing API; extend only general metadata/mapping needs. |
| `thd_mvpa.[ch]` | MVPA preprocessing, classifiers, selection/PCA workflows, CV policies, permutations, audits, activation patterns and frozen-model formats. | Ownership of the extracted regression implementations. | Include `thd_regression.h`, retain workflow APIs and serialized behavior, replace private solver-field access. |
| `3dPicard.c` | CLI, xmat/manifest ingestion, domain checks, mask and response batching, worker scheduling, AFNI outputs, provenance. | Numerical algorithms and candidate scoring logic. | New executable adapter; start with existing 3dDeconvolve matrices. |
| `thd_patterns`, `thd_sampleset`, `thd_datatable` | Existing spatial extraction and table/sample identities. | A mandatory ROI/searchlight representation for voxelwise fitting. | Reuse only where needed; begin with a local bounded response reader in the executable. |
| `thd_rsa`, `thd_simmatrix`, `thd_mapinfer`, `thd_permute` | Their established representational/inference responsibilities. | Picard's initial solver or model-selection engine. | No Picard-driven changes or new inference dependencies. |

### Stage 1: exact estimator extraction and compatibility surface

Move the following declarations from `thd_mvpa.h` into `thd_regression.h`,
re-exported by including that header from `thd_mvpa.h`:

- `THD_RIDGE_AUTO/PRIMAL/DUAL`, opaque `THD_ridge_workspace`, and every
  `THD_ridge_*` allocation/fit/predict/coefficient/free entry point.
- `THD_REG_*`, `THD_regression_config`, `THD_regression_diagnostic`, opaque
  `THD_regression`, and its allocation/fit/predict/coefficient/diagnostic/free
  functions. Keep enum values, configuration meaning and function signatures.

Move the matching implementations and their private helpers/structs into
`thd_regression.c`: `ridge_fit_penalties`, `reg_soft`, `reg_kkt`, `reg_elastic`,
`reg_banded`, `reg_pls`, and `reg_svr`. Keep all four completed continuous
estimators together. PLS and SVR remain supported for MVPA, without becoming
Picard features. Preserve existing float input/output boundaries, double
accumulation, intercept rules, parameter units and error behavior in this step.

The source has real internal coupling that must be resolved before moving it:

| Existing access in `thd_mvpa.c` | Required boundary |
|---|---|
| `ridge_record_patterns` reads ridge predictor means, coefficients and response means. | Read-only fitted-state access retaining double precision and the current centered-score arithmetic; do not round through the float prediction API. |
| `regression_cv_mean` reads ridge/regression response means. | Numeric training-response means accessor, including SVR where mean is distinct from fitted intercept. |
| `regression_cv_record` reads coefficients, predictor means, group IDs and diagnostics. | Read-only model state plus existing diagnostics accessor. Selection remapping and descriptive group variance/covariance stay in MVPA. |
| `ridge_candidate_evaluate` reads regression diagnostics directly. | Use `THD_regression_diagnostics`; keep fold/candidate aggregation in MVPA. |
| `THD_mvpa_ridge_freeze` reads fitted flag, active dimensions, means and coefficients. | Copy state into the existing frozen-model representation through checked accessors. Preserve serialization fields, arithmetic and prediction behavior. |
| `reg_banded` borrows `THD_ridge_workspace` factor/RHS buffers. | Keep both implementations in the same numerical module; no public scratch-buffer accessor. |
| `reg_banded` calls the classifier-local `lda_solve`. | Give the shared numerical module a narrow SPD solve helper, callable by LDA through a module-internal header. Preserve arithmetic and caller-specific failure messages; do not move the LDA model or create a general algebra framework. |

Concrete accessor design: add `THD_ridge_get_state` and
`THD_regression_get_state` returning a borrowed, const numeric view containing
active n/p/q, predictor means, response means, coefficients, and applicable
group membership. A request fails if the model is unfitted. Views expire on
refit/free; callers may copy into owned storage, but cannot mutate solver state.
Existing coefficient-export functions continue to return caller-owned copies
and original-input-space intercepts. Header/type spelling is proposed; these
functions are not implemented yet. New public headers contain no AFNI datasets,
THD_sampleset, THD_foldplan or permutation types. A private
`thd_regression_internal.h` may declare the SPD helper for LDA; it must not
expose workspace structs to MVPA as a shortcut around these accessors.

Leave `THD_response_ridge`, `THD_mvpa_ridge_workspace`, `THD_mvpa_set_regression`,
`THD_mvpa_run_ridge_cv`, tuning/permutation runners, frozen-model code, all
classifiers and `THD_preprocessor` in `thd_mvpa`. `THD_pca_*` and activation-
pattern routines also stay there initially. A later noise-pool PCA requirement
must justify its own extraction; PLS is not a substitute for nuisance PCA.

### Stages 2 and 4: prepared designs and parameter paths

Add these capabilities inside `thd_regression`, alongside legacy wrappers:

| Proposed object/operation | Contract |
|---|---|
| `THD_regression_design` / prepare | Own a double copy of the supplied effective X and immutable decomposition/cache state. X is already row-selected, nuisance-projected and scaled by the caller. No implicit centering/intercept. Backend capabilities are explicit. |
| `THD_regression_workspace` / allocate | Own per-worker response, coefficient, residual and convergence scratch for bounded n/p/q capacities. Borrow the prepared design; allocate before numerical loops. |
| Fit fixed settings / fit path | Accept double response batches and explicit fraction/lambda/mixing/group settings. Return coefficients and diagnostics; do not choose the winning candidate. Cold start is explicit; warm starts only within a compatible response/design path. |
| Path result sink | Deliver one candidate's coefficients and per-response actual parameters/diagnostics to the workflow, then reuse buffers. Do not require a full path x voxel x coefficient tensor. Callbacks must copy anything retained. |
| Fit selected settings | Refit a response batch using per-response selected settings on a design prepared from all eligible training rows. Share compatible numerical work; never substitute a tuning-fold fit for the final fit. |

Stage 2 implements only the SVD/fractional-ridge backend and its endpoint rules.
Stage 4 adapts the existing elastic-net coordinate-descent kernel to prepared
inputs and path warm starts. Current `THD_regression_fit` zeros coefficients
for each candidate; do not describe warm starts or screening as already present.
Stage 6 adapts the existing banded-ridge kernel, preserving the shared numerical
implementation rather than writing a second banded solver. A candidate may
share immutable Gram data but needs mutable factorization scratch per worker.

Existing APIs retain their centering and legacy penalty units. New prepared
fits accept already transformed doubles and never recenter them. Internal
kernels should serve both interfaces; changes to kernel plumbing require
legacy-parity checks. For normalized Picard lambda, diagonal ridge/group
penalties use n_train*lambda; completed MVPA elastic-net alpha is the normalized
strength, whereas ordinary/banded ridge lambdas are unnormalized. Fractional
ridge reports which convention its actual lambda uses. `alpha`, `rho` and
candidate ordering must not silently change meaning for existing MVPA users.

Prepared state is keyed by effective X, row selection, nuisance span, scaling,
whitening and candidate design. Response-specific fractional lambdas do not
invalidate an otherwise shared SVD. Different folds or HRFs generally do.
Do not reuse response-dependent path state between voxels or folds. Existing
PLS joint-response domains remain unchanged; Picard's ridge/EN batching is
independent-response fitting, not multitask estimation.

### Stage 2: projection, workflow and data ownership

`thd_linproj` exposes an opaque `THD_linproj` with prepare, project,
solve/recover, rank-report and free operations. It factors exactly the supplied
Z in double precision, including an intercept only if present. A rank-revealing
QR/SVD policy is explicit; no normal equations for the projection. Reuse one
factorization to project X and each Y batch. Recover Gamma from Y-XB after
coefficient back-transformation. Preserve column mappings when Z is redundant
and distinguish identifiable fitted nuisance signal from nonunique coefficients.

MVPA's `THD_preprocess_set_nuisance` remains predictor-side adjustment with its
current strict rank/intercept behavior. Do not route it through the new module
in the initial extraction. Sharing its underlying mechanism later requires a
separate parity change; the same word "nuisance" does not make the two policies
interchangeable. Whitening is supplied consistently to Z/X/Y by a later workflow
adapter; `thd_linproj` does not estimate temporal noise or know about runs.

`thd_picard` owns two classes of state:

- An analysis/design plan: owned copies of numeric designs, column-role maps,
  retained-to-full row IDs, run/trial identities needed by the selected mode,
  scaling configuration, candidate descriptions, and validated fold plans.
  The executable parses files into these inputs; the workflow never reads them.
- A worker workspace: response-batch scratch, training-only transforms,
  references to prepared designs, scores, candidate diagnostics, selection and
  final-refit outputs. Borrowed immutable plan state outlives every worker.

Proposed workflow calls are `THD_picard_plan_create`,
`THD_picard_workspace_alloc`, `THD_picard_fit_batch`,
`THD_picard_tune_batch`, and matching free operations. These are interface
requirements, not frozen C signatures. Public data remain numeric arrays and
explicit metadata IDs, with bounded errors returned to the caller. Fatal CLI
exits, logging destinations and file ownership belong to `3dPicard.c`.

The executable loads bounded Y blocks in retained-row order from full-timeline
AFNI data, applies dataset scale factors, and converts to double once. The xmat
adapter converts column-major floats to compact row-major doubles once per input
plan, without claiming to recover lost text precision. X is [row][predictor],
Y is [row][response], and B is [predictor][response], using active dimensions.
Z and nuisance coefficients follow the same conventions. Output mapping and
float conversion occur only at the executable boundary. Existing MVPA buffers
and layout are not changed by this decision.

Each worker uses a private workspace and writes disjoint result slots. File
writing is serialized by the executable. No mutable global solver settings;
failed fits invalidate their result and produce an explicit status/error, not
stale coefficients. Bound cached designs and Y/path scratch with a memory
budget; chunk size or worker count must not change the selected estimator.

### Stage 3: folds and selection boundaries

Reuse `THD_foldplan`, `THD_nested_foldplan` and their validators. Add a numeric
constructor from row/run/group IDs so Picard does not fabricate dataset-backed
THD_sampleset objects. Temporal blocks/gaps use original time/run coordinates,
then map to retained rows; gaps are not distances on the compressed censored
axis. Preserve the current non-fixed-plan requirement of one test assignment
per sample. Repeated/partial test coverage needs its own later contract.

`THD_nested_candidate_search` returns one scalar winner and assumes a nested
outer/inner plan. It is reusable for scalar/global selection, but is not already
a per-voxel path selector or a final all-data tuning driver. Do not run it once
per voxel if that repeats design preparation. Add a separate generic
`THD_cv_path_search` traversal when stage 3 needs it: explicit inner plan,
prepare callback once per fold, path-evaluation callback once per fold/response
batch, and finish callback. Results and scoring live in the caller context.
Keep the existing scalar function and its semantics unchanged. Reuse validated
index mapping and lifetime rules; nested traversal may compose the new driver.

`thd_picard` supplies the callbacks and owns fold-by-candidate-by-response losses,
R2 baselines, one-SE summaries, per-response winners and optional global
aggregation. Store fold-level loss/weight information needed by the selected
rule; stream or batch audit output. Final refits use a separately prepared
full-training design. Prediction CV and repeated-condition beta tuning are
separate workflow modes; only generic scheduling is shared. Do not carry over
MVPA's row-permutation inference into temporally dependent Picard data.

### Build and validation boundary

For extraction, update `src/CMakeLists_binaries.txt`'s explicit 3dMVPA source
list and `src/Makefile.INCLUDE`'s MVPA object/header dependencies to compile
`thd_regression.c` exactly once. Register Picard with `thd_picard`,
`thd_linproj`, `thd_regression`, `thd_crossval` and the required xmat/AFNI adapters.
Account for current crossval/xmat transitive table/sample link dependencies;
removing those dependencies is not part of this change. Do not add thd_mvpa
or inference modules to Picard merely to resolve numerical symbols. Leave
3dDeconvolve, REML and LSS target definitions untouched.

Stage-1 acceptance includes the existing tests registered as `3dMVPA:ridge`,
`response-ridge`, `regression`, `regression_null`, `regression-surface`,
`nuisance`, `activation-patterns` and `saved-models` (with the 3dMVPA prefix),
plus affected nested-tuning, classifier/LDA, coefficient and pipeline coverage.
Retain numerical fixtures in `src/pmolfese/tests/test_mvpa.c`; move or add
solver-only coverage in `test_regression.c` as appropriate, without removing
MVPA integration assertions. Validate all four moved estimator branches, not
only those exposed by Picard. Capture a pre-extraction reference build/results.

New focused tests are `test_linproj.c` (span/rank/recovery),
`test_picard.c` (numeric workflow and training isolation), and Picard CLI/oracle
fixtures for censoring, labels, batch/thread invariance and fraction selection.
Names are planned. Test symbols/ownership, failed-fit invalidation, accessor
lifetimes and double precision where activation patterns depend on it. No
scientific tests were run for this documentation-only boundary review.

## Suggested implementation sequence and acceptance gates

These are the Picard milestones; their existing stage IDs are retained.
Stages 0-2 and the first selected tuning milestone (3a or 3b) should use
unchanged 3dDeconvolve-produced matrices and explicit metadata fixtures.
Choose 3a for an encoding-first workflow or 3b for a single-trial-first workflow.
This gives the generator an already tested consumer and separates fitting
errors from design-generation differences.

| Stage | Deliverable | Gate before proceeding |
|---|---|---|
| 0 | Capture the completed MVPA reference build; specify model roles, units, scaling, state accessors and fixture identities. | Reviewed extraction checklist and baseline results; library boundaries above agreed; no dependency on Maker. |
| 1 | Create `thd_regression.[ch]`, move ridge plus all four continuous estimators, replace MVPA private-field access and update build lists. | Behavior-preserving extraction: ridge/EN/banded/PLS/SVR, LDA helper, CV, patterns, coefficients and saved-model coverage passes. No new numerical algorithm in this change. |
| 2 | Add SVD/fractional paths to `thd_regression`, projection/recovery in `thd_linproj`, fixed-design workflow in `thd_picard`, and the `3dPicard.c` adapter. | Existing 3dDeconvolve matrices work; fracridge and direct-model oracles pass, including rank/endpoint cases; nuisance recovery, labels, batch/thread behavior and memory bounds verified. |
| 3a | Extend numeric fold construction/path traversal in `thd_crossval`; implement shared-feature prediction scoring/selection/final refit in `thd_picard`. | Fold-local preprocessing and held-out perturbation checks; defined score baseline, deterministic per-response/global selection and undefined-score handling. |
| 3b | Add separate repeated-condition beta-tuning mode to `thd_picard` using explicit event identities and shared generic traversal where applicable. | Hand-checkable repeated-trial fixture; reject unsupported tuning when repeats are absent; component-level GLMsingle comparison. |
| 4 | Extend the extracted EN kernel with prepared-double fits and warm-start paths; add Picard sparse-penalty configuration. | Matched objective/scaling/folds against glmnet or scikit-learn; full KKT and warm/cold agreement; legacy MVPA EN results preserved. |
| 5 | Add candidate-design plans, identity validation and HRF-selection policy to `thd_picard`; ingest manifests in the executable. | Correct per-voxel design IDs and aligned beta identities; independent evaluation when advertised; ordinary single-matrix fitting remains independent of Maker. |
| 6 | Adapt extracted banded ridge to Picard's prepared fits; add prewhitening/support-refit adapters and measured optimizations as separate subchanges. | Legacy banded results preserved; per-feature numerical/scientific gates, runtime/peak memory and thread-count checks pass; no duplicate banded solver. |
| 7 | Single-trial orchestrator: noise-pool PCs, PC-count selection, staged outputs. | Pinned-reference component and end-to-end comparisons; documented differences from GLMsingle. |

Stages 3a and 3b are separate deliverables, not interchangeable implementations.
For a single-trial-first release, prioritize 3b and stage 5 over sparse penalties.
MVPA's numerical prerequisites are now complete. Stage 1 is the planned first
code change, with stage 2 building on the extracted library. Do not copy EN or
banded implementations into Picard or broaden the extraction into a rewrite of
MVPA's preprocessing/inference system.

The recommended handoff to Maker implementation is after stage 2 and one of
3a/3b pass their acceptance gates. Maker requirements and reference fixtures
can be recorded earlier. Picard need not finish stages 4-7 before Maker begins;
conversely, those stages can continue using externally prepared matrices.
Picard stage 5 and Maker M6 meet at the candidate-matrix/identity contract.
End-to-end HRF-library testing follows when both are ready. Picard must always
accept a single ordinary matrix independently of the candidate-library feature.

Update classic Makefile dependencies and CMake source/target registration when
libraries/programs are actually introduced. Move solver-focused tests with the
extraction, retaining MVPA integration coverage. Add tests for scientific failure
modes and compatibility; avoid tests that merely repeat the new implementation.

Use synthetic known-answer tests before whole-brain comparisons. Compare
predictions/objectives as well as coefficients where collinear lasso solutions
are nonunique. Evaluate LSA/LSS/ridge split-half RDM reliability through 3dRSA
as an end-to-end scientific benchmark, not as proof of solver correctness or
guaranteed improvement on every dataset. Hyperparameter sharing across halves
can couple estimates; disclose and control that in reliability experiments.

## Performance and deferred decisions

- A design decomposition can be shared only for identical effective designs:
  same retained rows/fold, scaling, nuisance projection, candidate HRF and
  whitening. Cache by these properties, not just by matrix filename.
- Voxel-specific ARMA whitening from 3dREMLfit can break one-SVD-per-design
  reuse. Specify import format, noise-estimation splits and covariance handling
  before promising `-prewhite`; never whiten across CV boundaries implicitly.
- Separate immutable prepared designs from per-thread mutable response scratch.
  Bound response batches and caches; do not allocate all voxels x all candidates
  x all coefficients. Include Gram-versus-residual memory thresholds explicitly.
- Repeated `-matrix` candidates must agree on time/sample identities and intended
  coefficient meanings. General candidate selection and GLMsingle HRF replication
  need distinct policies; the latter does not simply select every stage by CV R2.
- Decide release scope, surface support, automatic versus explicit column roles,
  default fraction/path grids, global versus per-voxel tuning, and exact output
  names before freezing the CLI. These choices remain proposals.

## 1dXmatMaker: independent design-generation roadmap

### Scope and compatibility boundary

The name is settled: **1dXmatMaker**. It is a standalone matrix producer
feeding **3dLSS, 3dREMLfit and 3dPicard**. The implementation language remains
undecided. `1dXmatTool` retains inspection, validation, selection and splitting
responsibilities. The producer's interchange contract
is the existing AFNI matrix format plus optional provenance/event/candidate
manifests. New fitting policies or mandatory new header fields must not be
required of existing consumers.

All current 3dDeconvolve response-function families and their relevant design-
generation semantics are the compatibility target. Establish a versioned
inventory from the source/help, including variants and aliases, before coding:
GAM/GAMpw, TWOGAM/TWOGAMpw, BLOCK/UBLOCK and duration-modulated variants,
TENT/CSPLIN and zero-endpoint variants, SPMG families, WAV, MION/MIONN,
POLY/SIN/TRIG, EXPR, and other accepted forms such as IGFUN. This list is a
starting inventory, not a declaration of exhaustive compatibility.

Cover ordinary, IM, AM1/AM2 and supported FSL timing inputs, event durations,
amplitude centering, local/global times, timing units/offsets, fractional-TR
onsets, normalization, run-edge truncation, baseline/nuisance columns, censoring,
labels, stimulus ranges and supported GLTs. Compatibility involves the complete
sampled design and metadata, not just HRF curves. Existing environment-dependent
behavior must be inventoried and made explicit or documented in the new tool.
Implementation may be staged, but full basis coverage remains the target;
unsupported features must fail explicitly rather than use substitutes.

Independent libraries (provisionally thd_basis/thd_design if C is chosen) may
adapt source from a recorded AFNI revision with attribution and applicable
license terms preserved, or implement independently against reference fixtures.
They belong to the new tool and potential new consumers; they do not replace
3dDeconvolve code. No textual inclusion of 3dDeconvolve.c or runtime delegation
to its executable is intended as the production architecture. Copied numerical
behavior is an intentional maintenance tradeoff, with upstream changes reviewed
and parity fixtures retained. Existing shared-library changes must not alter
3dDeconvolve behavior or add dependencies to its build.

For Picard/GLMsingle-oriented work, design creation should additionally support
sampled custom HRFs and versioned HRF libraries, explicit sampling/interpolation
and normalization policies, duration handling, condition-collapsed ONOFF and
single-trial designs, stable trial/condition/run identities, and candidate
matrices with identical row/column meaning. Preserve whether an imported HRF
is an impulse response or already includes stimulus-duration convolution, so
duration is not applied twice. Reproducing the GLMsingle library requires its
actual reference waveforms and generation/resampling convention, not merely
twenty similar analytic HRFs. Record source version and redistribution terms.

The generator creates candidates and accepts supplied nuisance vectors. Voxelwise
HRF selection, noise-pool estimation/PCA, PC-count tuning, fractional ridge and
post-fit beta rescaling remain in Picard or its future orchestrator. A manifest
can describe multiple ordinary xmat files; LSS/REML consume a selected compatible
file individually, while Picard can later consume the candidate collection.
Generation does not remove LSS's one-IM-stimulus and rank requirements or REML's
estimability requirements. Provide consumer-specific validation and per-run
outputs/splitting where appropriate, without rejecting rank-deficient designs
globally when a penalized Picard fit can support them.

With the file interface and untouched 3dDeconvolve, a standalone Rust generator
is more feasible than in the shared-C-refactor proposal. C still offers direct
reuse of existing AFNI basis/expression-parser behavior and xmat infrastructure.
Evaluate language choice against full compatibility effort, build/distribution
support and maintainability; do not assume a speedup from language alone. A Rust
front end with C numerical routines is also possible, with explicit FFI costs.

### Milestones and acceptance gates

M0-M8 are stable IDs for Maker work, separate from Picard's stage numbers.
All are planned. Each milestone should end with reviewable outputs, documented
coverage and a validation report; partial coverage must not be advertised as
complete 3dDeconvolve compatibility.

| ID | Deliverable | Dependencies | Acceptance gate |
|---|---|---|---|
| M0 | Compatibility inventory, pinned reference corpus, language/build decision, and input/output specification. | Picard contracts; fixtures may be collected before its initial release. | Every current basis family, variant/alias, timing mode and relevant option has a named test or an explicit unresolved item. Record AFNI reference revision, environment settings, precision and numerical tolerances. |
| M1 | Independent executable and design-building core; explicit run lengths/TR; baseline and supplied nuisance columns; xmat metadata/writer integration. | M0; recommended start after Picard stage 2 and 3a or 3b. | Baseline/nuisance-only fixtures match reference values and metadata; unequal runs and censor maps round-trip correctly; supported builds require no changes to 3dDeconvolve or its dependencies. |
| M2 | Core event regressors and single-trial expansion: GAM/BLOCK and selected TENT/SPMG cases, ordinary and IM designs, stable event-to-column mappings. | M1. | Match off-grid and boundary-event fixtures, repeated conditions, empty runs where valid, and label/order conventions. Fit matched generated/reference matrices with Picard; validate a suitable one-IM design with 3dLSS. This is an explicitly partial prototype. |
| M3 | Complete basis-family inventory, including variants, duration convolution, expression-defined models and legacy aliases. | M2. | Every M0 basis entry passes curve and sampled-column tests across documented parameter cases. Invalid parameters fail clearly; no approximation or fallback substitutes for an unsupported family. |
| M4 | Complete timing/modulation and design semantics: AM1/AM2, variable duration, FSL forms, units/offsets, local/global times, normalization, baseline/nuisance combinations and supported GLTs. | M3. | Full-design parity corpus passes, including interactions between basis, modulation, run boundaries and censoring. Relevant environment-dependent behavior is explicit and recorded. Metadata and contrast mappings agree with their intended columns. |
| M5 | Consumer integration and compatibility release: Maker -> xmat -> LSS/REML/Picard; Tool inspection/splitting interoperability. | M4; initial Picard fitting gates passed. | Matched matrices yield equivalent downstream fits within specified tolerances. Unsupported consumer designs receive specific diagnostics; valid penalized designs are not rejected by a universal OLS rank rule. |
| M6 | Sampled custom HRFs and versioned libraries; ONOFF/condition/trial design variants; candidate manifest with stable identities. | M5; candidate contract agreed with Picard stage 5. | Match a pinned GLMsingle HRF-generation/convolution reference for selected fixtures; verify sample timing, normalization, duration treatment and IDs across candidates. Test actual library waveforms as well as custom inputs. |
| M7 | Picard candidate-design workflow and downstream beta-manifest integration. | M6 and Picard stage 5. | Picard produces correct candidate IDs and trial-labeled betas; RSA/MVPA adapters preserve identities. LSS/REML can still read an individual compatible candidate without understanding the manifest. |
| M8 | Performance, reproducibility, packaging and workflow documentation. | M5 for AFNI compatibility; M7 for library workflow claims. | Report wall time and peak memory against matched 3dDeconvolve matrix-only runs, including diagnostic settings. Repeat-build/run checks pass; examples cover all three consumers and make no unsupported speed/equivalence claims. |

M3 and M4 complete the coverage promised by "all current basis functions and
design semantics"; M2 alone does not. Coverage follows the pinned reference,
with an explicit review process for later AFNI additions. Cross-product cases
matter: a family passing an isolated waveform test does not establish its
duration-modulated or multi-run behavior. If a combination is rejected by the
reference, record that rather than inventing a supported equivalent.

### Inputs, products, and interface decisions

Keep the first interface familiar to AFNI users: explicit acquisition timing,
event timing/basis specifications, run definitions, nuisance inputs and censor
information. A named-stimulus or declarative interface can be added over the
same internal model; it must not require a second numerical implementation.
Generation should work from timing and acquisition metadata without loading
voxel time series. Any optional dataset input should be used only for needed
acquisition metadata, with consistency checks against explicit settings.

Proposed output contract, with exact option/file names settled in M0:

- An ordinary censored `.xmat.1D` and an explicitly requested or default paired
  full-row matrix, with shared column meanings and consistent GoodList/NRowFull.
  Build regressors on the original run timelines before applying row censoring;
  censoring must not compress event times or change convolution.
- Complete labels, groups, stimulus ranges, basis metadata, run starts, TR,
  provenance and supported contrast metadata needed by existing consumers.
- An event-to-column manifest identifying event/trial, condition, item if
  supplied, run and basis component. Do not equate matrix-column indices with
  output beta-volume indices without the fitter's mapping; multi-basis events
  can generate several columns.
- For M6, one ordinary matrix per candidate plus a versioned manifest recording
  candidate ID, file path, HRF source/version, sampling, duration and normalization
  conventions. Keep row and coefficient identities aligned; varying basis
  dimensionality requires a separate explicit comparison contract.
- A validation/coverage report distinguishing malformed designs from consumer-
  specific estimability limits. Maker must not infer penalty choices or silently
  remove columns to satisfy a particular fitter.

Prefer reusing 1dXmatTool for inspection and per-run splitting over duplicating
its command interface. LSS examples should use compatible run-local designs
with the matching full-length run data. REML and Picard examples should retain
the correct full/censored timeline mapping. Existing fitting programs need no
new Maker-specific metadata to perform their ordinary matrix-based fits.

### Implementation and testing policy

- C is the current recommendation for compatibility effort; Rust and a hybrid
  remain alternatives until M0. Record the decision against complete basis and
  EXPR support, supported platforms, build integration, ownership and maintenance.
- If C is chosen, tentative independent modules are `thd_basis.[ch]` for basis
  parsing/evaluation and `thd_design.[ch]` for event/run/design assembly, using
  `thd_xmat` where its precision and metadata contract fit. These names and the
  split are proposals, not changes required of 3dDeconvolve.
- Pin reference outputs before adapting numerical routines. Retain independent
  analytic/hand-computed fixtures too: copied implementation and copied oracle
  logic are not independent validation. Compare arrays and semantic metadata;
  source-command provenance need not be byte-identical.
- Full compatibility fixtures cover normalization and scaling, fractional onsets,
  run-edge clipping, timing placeholders, modulation centering, duration handling,
  all-zero or censored trial support, redundant nuisance terms, and contrast
  indexing. Errors and valid-but-unestimable matrices must remain distinguishable.
- Check duplicate labels, nonfinite values, invalid event/run mappings and output
  collisions before writing. Define cleanup behavior for partially written
  candidate collections and preserve reproducible IDs on rerun.
- Benchmark generation and optional diagnostics separately. A comparison against
  `3dDeconvolve -nodata ... -x1D_stop` must record whether condition checks and
  full-row outputs are enabled. HRF-library batching/caching is a separate
  benchmark; no language or executable split guarantees a speedup.
- Document the independent maintenance boundary: upstream changes are reviewed
  and ported deliberately, not pulled in automatically. Tests may invoke an
  unchanged 3dDeconvolve reference executable; production Maker does not.

## Reference material

- Original user-supplied design notes: "AFNI idea: 3dRidge / 3dPenFit (ridge +
  lasso + elastic net)". This plan incorporates corrections discussed in chat.
- [Fractional ridge reference implementation](https://github.com/nrdg/fracridge)
- [glmnet objective and solver documentation](https://glmnet.stanford.edu/articles/glmnet.html)
- [GLMsingle methodological FAQ](https://glmsingle.readthedocs.io/en/latest/wiki.html)
- [GLMsingle implementation and HRF-library interface](https://github.com/cvnlab/GLMsingle/blob/main/glmsingle/glmsingle.py)
- [3dDeconvolve basis and timing documentation](https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dDeconvolve.html)
- Related tool documentation: `src/pmolfese/1dXmatTool.md`.
- Local source reviewed: `src/thd_mvpa.[ch]`, `src/thd_crossval.[ch]`,
  `src/thd_xmat.[ch]`, `src/thd_patterns.h`, `src/thd_sampleset.h`,
  `src/thd_rsa.[ch]`, `src/3dRSA.c`, `src/3dMVPA.c`, `src/3dTfitter.c`,
  `src/3dREMLfit.c`, `src/CMakeLists_binaries.txt`, `src/Makefile.INCLUDE`,
  and `tests/CMakeLists.txt`.
- Related roadmap: `src/pmolfese/3dMVPA_PLAN.md`. Its milestone status remains
  authoritative for MVPA; this document does not advance or replace it.
