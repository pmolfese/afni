# 3dPicard: voxelwise penalized regression and the 1dXmatMaker roadmap

Created 2026-10-29. Status: design proposal; no implementation authorized by
this plan. The program name is settled: **3dPicard**, a reference to Jean-Luc
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

This records the initial proposal and a source-informed refactoring assessment.
The user is actively developing elastic-net support for 3dMVPA. References below
describe the working tree inspected on this date, including unfinished changes;
they do not certify those changes as completed or validated. Recheck the settled
MVPA implementation before any extraction. This planning change touches no code.

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
and ongoing continuous-estimator work. Avoid creating a third independent
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

## Refactoring assessment from the current source

### Reuse established layers

| Existing module | Evidence in the current tree | Proposed use by Picard |
|---|---|---|
| `thd_xmat.[ch]` | Reader, raw attributes, run/censor maps, column selection, diagnostics, explicit trial manifests. | Reuse parsing and mappings. Add narrowly scoped role/row helpers only when necessary; keep penalty policy in the caller. |
| `thd_crossval.[ch]` | Explicit train/test lists, grouped and nested plans, callback-based candidate search. | Reuse plan representation and compatible search infrastructure; add time-block/gap construction and per-response/fold-loss support where needed. |
| `thd_datatable.[ch]` / `thd_sampleset.[ch]` | Metadata and canonical sample/dataset mappings used by MVPA. | Reuse for manifest/table interoperability. Do not require a dataset-backed sampleset just to describe a numerical time-series fold. |
| `thd_patterns.[ch]` | Spatial extraction, ROI/searchlight support, surface domain helpers. | Reuse appropriate domain/extraction utilities; measure batching needs before adding a small time-series extraction helper. Voxelwise fitting does not require searchlights. |
| `thd_rsa.[ch]` / `thd_simmatrix.[ch]` | RDM construction and representational geometry. | Keep as downstream consumers for beta reliability validation. They are not penalized-regression engines. |

### Extract the estimator layer from thd_mvpa after active work settles

The strongest candidate is a new `thd_regression.[ch]` containing dataset-
independent numerical fits. At inspection, `thd_mvpa.c` contains:

- `THD_ridge_*`: centered primal/dual Cholesky ridge, positive penalties,
  multiple responses, and sharing of factors for equal response penalties.
- `THD_regression_*`: a continuous-estimator interface with elastic-net,
  banded-ridge, PLS and SVR branches in the working tree. Their presence is
  not an assertion that the corresponding roadmap milestones have shipped.
- `reg_elastic`: residual-update coordinate descent and full KKT checks.
  `THD_regression_fit` currently zeros coefficients for every candidate, so
  a warm-started path would be an extension, not an existing capability.

Extract the stable estimator code mechanically first, retaining function names
and compatibility declarations/includes for existing callers. A small movement
of the whole existing continuous-estimator block may be safer than splitting
its internal shared buffers by penalty family. PLS/SVR can remain available to
MVPA without becoming Picard features. Decide exact extraction size after the
active changes settle, and preserve model serialization and CLI behavior.

Then introduce an explicit prepared-design interface: factor/prepare X once,
fit bounded blocks of Y, and evaluate parameter paths. Keep independent response
penalties and diagnostics separate from shared design state. Add fractional
ridge as an SVD backend; it cannot be implemented simply by relabeling existing
positive-lambda Cholesky fits. A fraction can map to different lambdas for
different responses, so one factorization per distinct lambda is insufficient
for an efficient whole-brain fractional path.

The extracted numerical layer must support caller-controlled intercept and
centering behavior. Existing MVPA entry points should preserve their centered
fit semantics through wrappers; Picard needs fits on its explicitly projected
design, particularly once whitening is supported. Do not silently recenter
a whitened/projected system using an ordinary constant vector.

Choose matrix layout and ownership once at the adapter boundary: xmat is
column-major float, while current MVPA solvers accept row-major float and
compute internally in double. Avoid repeating conversions per voxel or candidate.
Provide a double-precision fitting path for prepared designs without forcing
existing MVPA callers to change storage conventions.

### Separate reusable projection mechanics from nuisance policy

MVPA's `THD_preprocess_set_nuisance` currently adjusts predictor columns before
scaling, adds an intercept, normalizes nuisance columns internally, and requires
a full-rank nuisance design in every training split. Picard's partial penalized
model needs projection of both X and Y and reconstruction of nuisance coefficients.
These are related numerical operations but different scientific policies.

A small `thd_linproj.[ch]` is a candidate if both consumers benefit: reusable
QR/SVD factorization, projection, coefficient recovery, explicit rank policy,
and caller-owned scratch. Keep MVPA's existing strict rank/centering policy in
its wrapper. Define Picard's treatment of redundant or fold-zero nuisance
columns separately. Preserve fitted nuisance span and predictor identifiability;
do not silently delete scientifically meaningful columns.

This extraction is optional for the first fixed-design prototype. Share a
concrete second use rather than creating a general linear-algebra framework.
Use existing AFNI factorization facilities where their precision and contracts
fit; audit `thd_lasso.c` for useful numerical precedent without inheriting its
file-static configuration or assuming it has the desired thread/path interface.

### Extend folds without importing the entire MVPA workflow

`THD_nested_candidate_search` already separates preparation, evaluation, loss
aggregation and tie-breaking. Keep it generic. Picard-specific fraction paths,
trial-repeat tuning, HRF selection and final voxel refitting belong in a thin
Picard workflow layer, provisionally `thd_picard.[ch]`, when needed.

Current grouped-fold constructors use THD_sampleset, and current validation
requires exactly one test assignment per sample for non-fixed plans. Temporal
gaps can exclude rows from training while retaining test coverage, but repeated
CV or arbitrary partial test coverage needs a deliberate extension, not bypassed
validation. A metadata-only constructor may be sufficient for early time-series
work. Trial-repeat beta tuning requires its own scoring contract; a shared fold
struct does not make it the same problem as prediction CV.

Keep MVPA classification, permutations, feature selection, activation patterns,
saved-model workflows and model comparison in their existing domain layer.
Keep RSA's RDM regression and inference there as well. Picard should export
compatible beta datasets and manifests; adding it must not require RSA to link
to the new penalized solver or changing its scientific nulls.

## Suggested implementation sequence and acceptance gates

These are the Picard milestones; their existing stage IDs are retained.
Stages 0-2 and the first selected tuning milestone (3a or 3b) should use
unchanged 3dDeconvolve-produced matrices and explicit metadata fixtures.
Choose 3a for an encoding-first workflow or 3b for a single-trial-first workflow.
This gives the generator an already tested consumer and separates fitting
errors from design-generation differences.

| Stage | Deliverable | Gate before proceeding |
|---|---|---|
| 0 | Settle design roles, objective units, scaling, output identities and tuning modes; review completed MVPA work. | Written contracts and reuse decision; no concurrent extraction of active solver changes. |
| 1 | Mechanically extract stable estimator code if ready; reuse xmat infrastructure. | Existing ridge/EN, preprocessing, CV, coefficient and saved-model checks remain equivalent for the moved code; supported builds link correctly. |
| 2 | Fixed-design fractional ridge with unpenalized nuisance, fixed fraction/path and diagnostics. | Independent fracridge oracle and direct augmented-model checks; zero/rank-deficient/high-dimensional cases; coefficient reconstruction. |
| 3a | Shared-feature run-wise prediction tuning and final refit. | Fold-local preprocessing and held-out perturbation checks; defined score baseline, deterministic selection and undefined-score handling. |
| 3b | Repeated-condition single-trial tuning with explicit event identities. | Hand-checkable repeated-trial fixture; reject unsupported tuning when repeats are absent; component-level GLMsingle comparison. |
| 4 | Lasso/elastic-net paths using the settled shared solver. | Matched objective/scaling/folds against glmnet or scikit-learn; KKT failures visible; warm/cold path agreement. |
| 5 | Multiple compatible design matrices and HRF selection. | Explicit selection protocol; correct per-voxel design IDs and aligned beta identities; independent evaluation when advertised. |
| 6 | Banded ridge, prewhitening, support refits, performance refinements. | Separate numerical and scientific validation per feature, measured runtime/peak memory, thread-count agreement within documented tolerances. |
| 7 | Single-trial orchestrator: noise-pool PCs, PC-count selection, staged outputs. | Pinned-reference component and end-to-end comparisons; documented differences from GLMsingle. |

Stages 3a and 3b are separate deliverables, not interchangeable implementations.
For a single-trial-first release, prioritize 3b and stage 5 over sparse penalties.
If MVPA extraction is not ready, stage 2 can start with new fractional-ridge code
without copying its existing elastic-net implementation.

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
  `src/3dREMLfit.c`, `src/Makefile.INCLUDE`, and `tests/CMakeLists.txt`.
- Related roadmap: `src/pmolfese/3dMVPA_PLAN.md`. Its milestone status remains
  authoritative for MVPA; this document does not advance or replace it.
