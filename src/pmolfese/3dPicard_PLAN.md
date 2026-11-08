# 3dPicard: next-generation single-trial estimation roadmap

Updated 2026-11-08. **Status: planning only; no next-generation implementation
or benchmark results are claimed.** This update authorizes documentation, not
implementation. New milestones use NG identifiers; completed stages 1–7 and
the previous roadmap are preserved below the horizontal rule. This section
supersedes historical proposals where they conflict. Review each milestone's
results before authorizing the next implementation stage.

## Next milestone — NG1

**white/AR(1)/ARMA covariance support with principled gap handling, plus the crossed simulation benchmark**

Here ARMA means **ARMA(1,1)**. White noise and AR(1) are nested special cases;
REML is a covariance-estimation method, not a fourth covariance family.
The objective is better single-trial estimates and preserved information for
MVPA, RSA, and related analyses. A larger solver menu or agreement with a
numerical reference does not establish a scientific advantage. Picard should
remain useful with ordinary AFNI matrices independently of Maker.

### Current foundation and missing capabilities

| Workstream | Implemented foundation | Missing next-generation capability |
|---|---|---|
| Combine HRF adaptation, denoising, and penalized GLS | Per-voxel selection among supplied compatible HRF designs; staged run-local noise PCs and PC-count selection; fractional ridge, elastic net/lasso and banded ridge; fixed scalar AR(1) applied before nuisance projection. | Supplied voxelwise covariance, estimated white/AR(1)/ARMA models, covariance preserved across censor gaps, and validation of the combined adaptive workflow. |
| Support the downstream analysis contract | Native beta labels/provenance, retained/full row maps, optional beta manifests, and the documented Picard -> 1dXmatTool -> MVPA/RSA route. | Machine-readable learning/evaluation boundaries, partition-specific fit orchestration, explicit estimands/scales and uncertainty contracts, and tests of downstream geometry rather than file interchange alone. |
| Make the shrinkage target scientifically explicit | Explicit penalized/nuisance roles, supplied coefficient units, zero-centered ridge, band-specific penalties, and repeated-condition or prediction tuning. | First-class fixed targets or condition/behavior-centered trial deviations, identifiable parameterization, target-aware tuning and exports, and validation that meaningful trial variation survives. |
| Inference | Coefficients and fit/selection diagnostics. | No current conventional t/F inference. Conditional uncertainty is proposed separately; no ordinary t/F labels for penalized or selected estimates. |

Source baseline: `3dPicard.c`, `thd_picard.[ch]`, `thd_regression.[ch]`, and
[beta manifest workflows](beta_manifest_notes.md). The dedicated
[Maker plan](1dXmatMaker_PLAN.md) records M0–M8 complete for its qualified native
ARM64 macOS release; older Maker-pending statements below are historical.
This planning review does not rerun or recertify those validation records.

### NG1 scope and sequence

NG1 contains reviewable substeps. Imported/fixed covariance support and a
reproducible crossed benchmark come first. Native covariance estimation is an
explicit later substep, not implied by accepting imported parameter maps.

| Step | Deliverable | Acceptance gate |
|---|---|---|
| NG1a | Freeze covariance, scale, gaps, estimation provenance and benchmark contracts; independent dense reference solver. | Hand-computable white/AR(1) cases; documented ARMA convention, stationary initialization and parameter domain; benchmark configuration and metrics fixed before comparisons. |
| NG1b | Common covariance preparation/application interface; fixed scalar and supplied voxelwise AR(1)/ARMA parameters; explicit gap policies. | Independent GLS/ridge agreement; identity and ARMA-to-AR(1) limits; censor/run/fold checks; compatibility with existing fixed AR(1) mode. |
| NG1c | Import 3dREMLfit parameter maps and group compatible effective designs for bounded reuse. | Matched-covariance unpenalized agreement with 3dREMLfit under declared grid/truncation/rank settings; geometry and provenance rejection tests; memory/thread measurements. |
| NG1d | Crossed simulation harness and first report using known covariance, supplied estimates and an external REML reference. | Reproducible seeds/configurations, independent numerical oracles, tuning separated from evaluation, paired Monte Carlo uncertainty, all failures reported. External AR(1)-restricted estimation must be explicitly implemented or verified, not inferred from an unrestricted ARMA result. |
| NG1e | Native covariance estimation, initially a documented two-stage procedure; run-shared versus run-specific and pooling policies explicitly scoped. | Likelihood/parameter recovery references, low-information diagnostics, fold-local estimation, end-to-end rerun of the crossed benchmark. If deferred, NG1 is reported as partial, with imported covariance support only. |

Do not begin by jointly optimizing HRF, PC count, covariance and penalty in one
large search. Establish fixed-covariance correctness first. Joint variance-
component estimation of shrinkage and temporal covariance is a later research
extension, particularly for datasets without repeated conditions.

### Covariance and penalized GLS contract

For a voxel, write `y = X beta + Z gamma + epsilon`, with
`Cov(epsilon) = sigma^2 R(theta)`. Normalize R to unit diagonal so covariance
scale does not silently change penalty meaning. For retained rows, factor
`R = L L'` and apply `W = L^-1` to **X, Z and y**. Project the transformed
nuisance span before penalized fitting; recover nuisance coefficients and
fitted/residual values in original units afterward.

The initial ridge objective is
`||W(y-X beta-Z gamma)||^2/(2n) + lambda ||beta||^2/2`, where n is the number
of retained observations. Specify conversions to the diagonal penalty
`n*lambda`, as well as any innovation-variance versus marginal-variance
normalization. The current AR(1) filter has a different overall normalization
from unit-marginal-correlation whitening; preserve its legacy behavior and
explicitly test penalty-unit conversion. Fractional ridge must use the
minimum-norm reference from the same covariance/nuisance/design transformation.

Implementation responsibilities:

- Audit existing AFNI REML covariance routines and ownership/build boundaries
  before choosing reuse or a Picard-local adapter. Do not modify 3dREMLfit or
  3dDeconvolve as an implicit part of this work. Keep shared regression kernels
  independent of AFNI dataset I/O and covariance estimation.
- Prefer extending the existing Picard files over introducing many new files.
  Define immutable prepared covariance/design state and worker-local scratch;
  no explicit dense inverse is required. Dense Cholesky is the small-case
  oracle; a production factorization must have measured memory/runtime costs.
- Accept fixed parameters and voxel maps through a validated schema: model,
  AR/MA sign convention, run scope, source, parameter grid, scale and estimation
  rows/design. Check image geometry, finite values, stationarity/invertibility,
  positive definiteness and mask coverage. Do not silently clip invalid maps.
- Import `3dREMLfit -Rvar` a/b values with verified parameterization. Its noise
  scale is a separate quantity; it is not automatically needed for beta fitting.
  Full-data imported estimates must be labeled as such, not presented as
  training-only estimates during CV.
- Cache only identical effective designs: retained rows, run/gap policy,
  covariance parameters, HRF, nuisance/PC columns, scaling and fold all matter.
  Any quantization of parameters is an explicit approximation to benchmark.
  Bound caches and response batches; avoid voxel-by-candidate coefficient cubes.
- Expose covariance fit status, boundary estimates, likelihood/rank diagnostics,
  model/parameter maps and fallback reasons. An unidentifiable fit must not
  silently become a successful white-noise fit. Default fallback policy remains
  a decision for NG1a, with explicit failure available.
- Native estimation initially profiles a documented Gaussian restricted
  likelihood on a specified pilot design with sufficient residual rank.
  Compare full-trial and compact designs; neither is universally preferred.
  Rank-deficient designs require rank-aware treatment; saturated designs cannot
  support ordinary residual REML estimation. Do not substitute penalized SSE
  into an ordinary REML formula and call it joint penalized REML.

### Principled censor-gap and fold handling

Proposed default for the new covariance mode: preserve temporal covariance
across missing observations within a run. Select the retained-row submatrix
from covariance defined on original acquisition times. For stationary AR(1),
retained samples separated by Delta TRs have correlation `rho^Delta`.
ARMA uses its corresponding autocovariance at the actual lag.

- Censored measurements remain excluded. Do not retain contaminated values,
  compress the clock, or interpolate values and treat them as observations.
- Keep independent run blocks. Offer an explicit reset/segment policy for
  declared process discontinuities. A censor flag alone is not proof of a
  discontinuity. Test isolated missing samples, long gaps, leading/trailing
  gaps, single-observation segments and entirely censored runs.
- Preserve legacy `-prewhite_ar1` reset semantics; new options and help must
  distinguish the new preserve-gap policy. No silent numerical default change
  in existing command lines.
- In CV, form training covariance on training observation IDs without reading
  held-out responses. Temporal dependence across a train/test boundary still
  affects evaluation: prefer independent runs or justified blocked splits and
  buffers. Never whiten the full response series and then split it.
- Optional predictions at censored times are a separate, labeled output and
  must not enter the fit, sample count, validation score or likelihood as data.

### Crossed simulation benchmark and possible paper

The central question is: **when do temporal covariance modeling, shrinkage,
HRF adaptation and denoising improve single-trial recovery and downstream
representations, and when do their estimation costs or biases outweigh gains?**
This is a candidate paper question, not a novelty or superiority claim.

Cross generating noise with fitted covariance:

| True noise | Fit white | Fit AR(1) | Fit ARMA(1,1) |
|---|---|---|---|
| White | Matched baseline | Unnecessary complexity | Unnecessary complexity |
| AR(1) | Ignored correlation | Matched | Extra parameter |
| ARMA(1,1) | Ignored correlation | Misspecified | Matched |

For every cell compare unpenalized and tuned ridge fits. Separate known true
covariance (where applicable), externally estimated covariance and native
estimation. For misspecified families, report the actual estimation procedure;
there is no matching true parameter to supply. Oracle penalties selected using
true betas are upper-bound probes only, separate from deployable tuning.

Specific questions and planned evidence:

| Question | Comparison / outcome |
|---|---|
| Does shrinkage still add value after adequate whitening? | Ridge versus GLS under identical known and estimated covariance; beta MSE decomposed into variance and bias. |
| When is ARMA worth its extra parameter? | Matched AR(1)/ARMA and white truth across run lengths/SNR; recovery, estimation instability and runtime, including near-canceling AR/MA parameters. |
| How much loss comes from covariance estimation rather than fitting? | Known-covariance versus estimated-covariance methods with identical penalties/tuning rules. |
| Can ridge mask covariance misspecification? | White/AR(1) fits under ARMA truth; selected shrinkage, residual structure, amplitude attenuation and trial mixing. |
| Which pilot noise model is reliable? | Full-trial, compact condition model and explicitly scoped pooling; vary trial variability and residual rank, report covariance error as well as beta error. |
| Does preserving correlation across censor gaps help? | Preserve versus reset on identical censor masks; random missingness and artifact-driven censoring analyzed separately. |
| Does HRF mismatch masquerade as temporal noise? | Correct, fixed-mismatched and library-selected HRFs crossed with covariance fitting. |
| Are PCs and whitening complementary? | Neither, PCs only, whitening only and both under shared spatial noise plus temporal noise. |
| Are useful trial differences preserved? | Known behavior/learning/adaptation effects, condition contrasts, amplitude calibration and trial-to-trial mixing. |
| Do better beta metrics yield better downstream analysis? | Multivoxel known geometry, RDM error, independent-run decoding and test-retest reliability; no tuning on evaluation labels. |
| Can an apparent benefit arise through shared selection? | Independent-partition versus shared HRF/PC/penalty/covariance choices; label the latter conditional/transductive where appropriate. |

Stage the experiment rather than running an unmanageable full factorial:

1. Fixed correct HRF, no learned PCs: the 3x3 core, multiple covariance
   strengths, white and near-white limits, positive/negative MA examples.
2. Stress tests varying trial overlap, run length, SNR, trial-amplitude variance,
   effective design rank and censoring. Hold noise marginal variance comparable
   across generating families; specify stationary initialization/burn-in.
3. Selected interactions with HRF mismatch, spatially shared noise and PCs;
   multivoxel ground truth for RSA/MVPA. Include structured trial variation
   that repeated-condition tuning might otherwise classify as noise.
4. Real independent-session/test-retest comparisons with pinned GLMsingle,
   3dREMLfit, LSA/LSS and relevant fmrilss configurations. Match preprocessing,
   design/scaling where possible and disclose intentional policy differences.

Record seeds, generator/estimator versions, design rank, tuning splits,
parameter grids, compiler, thread counts, wall time and peak memory. Use paired
realizations and confidence intervals/Monte Carlo standard errors; choose
replicate counts by precision requirements. Include failures and boundary fits.
Separate deterministic numerical CI tests from expensive scientific benchmarks.
Residual whiteness is diagnostic, not proof of better betas; high reliability
alone can reflect excessive shrinkage. Report amplitude bias and geometry too.

Before a paper claim, perform a systematic prior-art comparison with GLMsingle,
fmrilss/fmriAR, covariance-aware trial estimation such as ITEM, and existing
AFNI methods. Report regimes of failure as well as improvements. The earlier
`sim_whiten_singletrial.py` is a pilot, not this benchmark: its ARMA truth is
fit with AR(1), it does not run 3dREMLfit, and its oracle penalty comparisons do
not establish deployable tuning performance or failure of ARMA REML.

## NG2 — Combine HRF adaptation, denoising, and penalized GLS

**Current position:** the orchestration already selects supplied HRFs per voxel,
locks that choice, builds run-specific PCs, selects PC count, and fits/tunes a
penalty. Fixed scalar AR(1) can participate. This is a substantial foundation;
we are missing estimated covariance and validation of the complete interaction.
HRF candidates preserve coefficient identities and equal projected ranks;
selection among arbitrary basis dimensions is a separate future contract.

Implementation after NG1:

1. Introduce an explicit stage configuration recording which design estimates
   covariance, when HRFs/PCs are selected and which choices are then frozen.
   Start with a declared pilot covariance followed by HRF selection, PC-count
   selection and penalized refit. Compare a limited alternative ordering in
   simulations; do not assert that this order is universally optimal.
2. Apply the selected covariance consistently to every candidate's task and
   nuisance columns and response. Existing PCA is computed in original time
   space after nuisance removal; retain that as the initial explicit policy,
   then whiten PC columns as nuisance during fitting. Whitening before PCA
   would be a different denoising method requiring its own validation.
3. Specify HRF scoring: with covariance fixed across HRF candidates, use a
   declared common weighted-fit criterion. If covariance is re-estimated for
   each candidate, raw residual SSE alone is not a likelihood comparison;
   determinant terms and differing fixed-effect spaces need careful treatment.
   Do not compare REML values across changing mean-model spaces uncritically.
4. Wire covariance-aware plans into fixed, repeat and prediction paths without
   duplicating regression solvers. Keep original-unit beta and output contracts;
   expose intermediate stages and the effective settings of each final voxel.
5. Implement partition-aware selection under NG3 before advertising independent
   downstream prediction performance for the adaptive pipeline.

Acceptance: each ablation matches its isolated reference; the joint pipeline
passes independent multi-HRF/PC/covariance/penalty fixtures, run/censor and thread
checks; held-out responses cannot change trained parameters in inductive mode.
Full GLMsingle equivalence is neither the criterion nor an implied claim.

## NG3 — Support the downstream analysis contract

**Current position:** beta-to-design-column mappings and external event
annotations already support MVPA/RSA interchange via 1dXmatTool. Ordinary
fixed-design sidecars may need annotations to supply condition/run/trial IDs.
This is not yet a complete guarantee of independent downstream evaluation or
of comparable uncertainty across trial estimates.

Implementation:

- Extend/version the manifest contract: subject/session/run/trial/condition,
  basis component, source design column, actual beta sub-brick, output stage,
  units/scaling, estimand, target/penalty groups, and evaluation partition.
  Distinguish total trial beta, deviation, condition mean and modulation slope.
  Preserve existing consumers and explicit annotation workflows.
- Record training identities for HRF, covariance, PC-count, penalty, scales and
  shrinkage target separately. A header saying merely "CV" is insufficient.
  Make conditional selection scores distinguishable from independent test scores.
- Add an outer-partition orchestrator and a frozen-settings representation.
  Learn settings on training runs and estimate held-out trial amplitudes under
  those settings. Held-out amplitudes necessarily use held-out observations;
  they must not update training settings or consume held-out condition labels
  to construct a supervised shrinkage target in a decoding evaluation.
- Resolve run-local PCs explicitly: learning temporal PCs from the held-out
  run's noise pool is test-data adaptation, even without labels. Offer a
  separately labeled transductive protocol or independently supplied nuisance
  terms/validated transferable noise model. Do not pretend temporal PC vectors
  from a training run can simply be copied into another run.
- Test permutation/reordering invariance, missing/duplicate trial rejection,
  multi-basis identities, optional stage outputs and partition isolation through
  Maker -> Picard -> Tool -> MVPA/RSA. Consumer changes need their own scoped
  review; do not silently alter existing RSA/MVPA statistics.
- Add optional conditional uncertainty for fixed covariance/design/penalty.
  With linear estimator `beta_hat = A y`, sampling covariance is `A V A'`;
  posterior covariance is a different object. Fraction-selected lambda and
  learned HRFs/PCs/targets are frozen for this conditional calculation.
  Export diagonals or requested contrasts first; evaluate compact operators
  or blocks before considering full per-voxel trial covariance arrays.
- Label omitted hyperparameter/model-selection uncertainty. Test conditional
  covariance against Monte Carlo for fixed settings; assess empirical coverage
  separately if intervals are offered. General t/F inference remains outside
  this milestone. Any later conventional GLS inference needs a fixed-design,
  estimability and degrees-of-freedom contract of its own.

Acceptance: round-trip identities/scales, partition leakage tests, independent
consumer fixtures, covariance oracles, and bounded output/memory cost. Document
that within-voxel trial covariance is not the full spatial noise covariance
needed by some multivariate distance estimators.

## NG4 — Make the shrinkage target scientifically explicit

**Current position:** ridge shrinks penalized coefficients toward zero in their
supplied units. Nuisance roles and banded penalties provide building blocks,
but there is no first-class condition/behavior-centered trial estimator.
An unpenalized condition-mean regressor plus all unconstrained trial columns
has a redundant decomposition; adding that column alone is not a complete
specification, especially at the unpenalized endpoint.

Implementation in two steps:

1. **Fixed target.** Support a declared target `m` through a penalty on
   `beta-m`. Fit deviations using `y-Xm`, reconstruct total betas and original-
   unit predictions, and record target source/units. Fractional norms refer to
   deviations from m and the matching unpenalized deviation reference, not
   silently to total beta norms. Targets learned from data obey training splits.
2. **Structured target.** Model trial amplitudes as `beta = H alpha + D delta`:
   H encodes protected condition means or prespecified behavioral effects;
   D spans a declared complementary deviation space, with constraints such as
   within-condition zero-sum deviations. Fit `XH alpha` unpenalized and shrink
   deviations. Define weighting and basis normalization so arbitrary coding of
   D does not change the intended penalty on trial-space deviations. Specify
   singleton groups, absent conditions, collinear modulators, all-censored trial
   support and rank-deficient designs. Check estimability after nuisance and
   covariance transformations, not just rank of H.

Export total trial amplitudes, protected effects and deviations as distinct
products with explicit labels. Handle nuisance recovery, fraction endpoints,
zero-reference cases, coefficient units and band penalties in the reparameterized
space. Establish ridge first; sparse penalties on deviations require additional
basis-dependence and scientific interpretation decisions.

IM designs are the primary target. AM1/AM2 designs ordinarily estimate weighted
or mean-plus-modulation effects, not one free amplitude per trial. Do not label
AM slopes as trial betas or apply repeat tuning without a valid identity and
estimand contract. Validate modulation centering, scale changes and protected
slope recovery separately. A hybrid trial-plus-modulation model needs the
explicit decomposition above rather than duplicate unconstrained columns.

Tuning must match the scientific objective. Repeated-condition agreement can
suppress real adaptation, learning or behavior-related trial variability.
Compare tuning for repeat reliability, prediction and preservation of known
structured effects; report their tradeoffs. In decoding, condition labels used
for shrinkage can inject the answer into the features. Require a label-safe
protocol and null-label tests; do not offer condition-informed held-out shrinkage
as an ordinary independent decoding workflow.

Acceptance: independent constrained/augmented-solve references, coding-invariant
ridge predictions under the declared metric, correct total/deviation reconstruction,
zero/unpenalized endpoints, fold-local targets, and known modulation/learning
signals preserved within justified tolerances. Covariance exports from NG3 must
be updated for the joint protected-effect/deviation estimator.

## Delivery order, evidence and boundaries

NG1 is next; NG2 integration, NG3 downstream guarantees and NG4 explicit targets
are subsequent plans. Freeze NG3 identity/split fields early enough that NG1/NG2
provenance need not be discarded. Start with ridge and fixed covariance, then
broaden to existing penalties only with estimator-specific validation. Do not
claim native covariance estimation when only import is complete.

For implementation: read test setup, verify the executable/module under test,
maintain both classic Make and CMake registrations, and validate focused targets
in both where supported. Preserve unrelated files and existing AFNI behavior.
Use independent small numerical oracles, sanitizer checks for changed memory
paths, and documented thread determinism before expensive scientific benchmarks.
A milestone report must state implemented versus deferred capabilities, tests
actually run, failed cases, performance conditions and scientific limitations.
No benchmark result or speedup is asserted by this planning update.

Reference anchors for implementation and paper review:

- [3dREMLfit covariance, censoring and parameter interface](https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dREMLfit.html).
- [GLMsingle methodological FAQ and repeat-tuning limitations](https://glmsingle.readthedocs.io/en/latest/wiki.html).
- [fmrilss prewhitening and residual-model discussion](https://bbuchsbaum.github.io/fmrilss/articles/prewhitening.html).
- [fmrilss OASIS: ridge and optional prewhitening](https://bbuchsbaum.github.io/fmrilss/articles/oasis_method.html).

<hr>

## Historical record — completed first-generation work and earlier proposals

The original plan follows unchanged. It contains completed implementation and
validation records **and** older proposals; placement below this divider does
not mark every historical proposal complete. Historical status statements,
including Maker progress and the November 7 pilot interpretation, must be read
with the corrections and current priorities above. Preserve this record for
traceability rather than treating its chronology as the active roadmap.

# 3dPicard: voxelwise penalized regression and the 1dXmatMaker roadmap

Created 2026-10-29; library-boundary review updated 2026-10-30;
stages 1, 2, 3a, 3b, 4, 5, 6 and 7 completed 2026-11-06.
Status: Picard stages 1, 2, 3a, 3b, 4, 5, 6 and 7 implemented; Maker milestones are separate.
This document alone does not authorize further implementation.
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
defines the independent generator and its later integration. Except for the
stage 1, 2, 3a, 3b, 4, 5, 6 and 7 completions recorded below, milestones remain planned; naming and
documenting them does not authorize implementation.

## Source consolidation and readability — 2026-11-06

User request: keep the Picard implementation in as few files as practical,
audit for `goto` statements, and explain the code for a recent CS graduate.
This is a source-organization/readability revision, not a new numerical phase.
Historical phase records below retain their original filenames; this section
and the current boundary table supersede those file locations.

The previous **18 C/header files are now six**:

| Current file | Contents |
|---|---|
| `3dPicard.c` | AFNI command, options, dataset/mask handling, candidate manifest parsing, staged orchestration, output/provenance and transaction cleanup. |
| `thd_picard.c` | Five numbered sections: nuisance projection; fixed/repeated-condition fits; held-out prediction/global selection; candidate/locked-HRF selection; noise PCs and PC-count scoring. |
| `thd_picard.h` | All Picard/projection APIs in matching sections. A guarded private section provides the plan layout to the implementation and its adapter. Ordinary clients retain opaque plans. |
| `thd_regression.c` | Shared MVPA/Picard solvers, including the thin SVD adapter and fractional-ridge implementation formerly in a separate file. |
| `thd_regression.h` | Shared public estimator/design/workspace contracts, units and result views. |
| `thd_regression_internal.h` | Small private Cholesky/SVD interface shared by regression, Picard and MVPA; keeps these implementation services outside the public estimator API. |

- Removed the separate `thd_linproj.[ch]`, `thd_picard_predict.[ch]`,
  `thd_picard_candidates.[ch]`, `thd_picard_single.[ch]`,
  `thd_picard_internal.h`, `thd_regression_svd.c`, and `picard_manifest.[ch]`.
  This is physical consolidation; public function names and statistical
  contracts stay the same. Consumers of the retired Picard headers now include
  `thd_picard.h`; manifest routines are file-local to the command.
- Private helpers have descriptive section-specific names to avoid collisions
  in the larger translation units. Added reading maps and expanded comments
  for matrix strides, ownership, lifetimes, projection, centering, primal/dual
  solves, warm starts, KKT convergence, PLS deflation, SVR bias, SVD storage,
  fractional-penalty bisection and strict manifest parsing. Expanded compact
  shared solver code into normally spaced, <=80-column C.
- Consolidation makes the shared regression object link to the existing
  EISPACK library even for MVPA; approved CMake registrations now declare it.
  Focused unit tests retain their separate behavioral coverage while linking
  the consolidated implementation and existing fold support. Classic Make
  has one shared regression-object recipe, avoiding duplicate rules.
- An executable-token audit, ignoring comments/formatting and accounting for
  private helper renames and manifest linkage, found no changes to computation
  or command logic. A lexer-aware scan found **zero `goto` statements** in all
  six current files. All six source/header files remain <=80 columns.
- Validation passed: both build systems, all 35 selected Picard/MVPA/shared-module
  tests, classic-binary Picard integration, strict C99 warning checks, and
  ASan/UBSan regression/staged-adapter checks. Evidence and pre-change copies
  are in `/tmp/picard_consolidation_reference` for this session.
- Existing AFNI implementation files, loaders, EISPACK and MVPA source were
  not edited. Changes are confined to the new Picard/shared-regression files,
  focused test includes, this record, and previously approved build/test
  registrations. Existing staged changes remain staged exactly as before;
  this revision is left unstaged for review. No commit was created.

Purpose: fewer files to navigate while retaining explicit numerical boundaries
and test coverage. No next milestone is started by this cleanup.

## Post-phase-5 AFNI I/O and readability revision — 2026-11-06

User authorization: use AFNI's dataset-loading system, leave fixes to the
NIfTI reader and `THD_load` for the user's separate work, and add more
explanatory source comments. This does not start phase 6.

- Both input volumes and masks now use `THD_open_dataset` / `DSET_load`.
  Removed Picard's raw NIfTI payload/byte-swapping/scaling path and separate
  payload ownership. Reader behavior, including conversion and NaN/Inf repair,
  belongs to AFNI; Picard no longer overrides `AFNI_FLOATSCAN`.
- Picard gathers doubles from loaded AFNI bricks and applies their brick
  factors. Numerical finite-value checks remain, but cannot detect original
  file values already repaired by the loader. Help now states this boundary.
- Added a source reading guide and comments explaining ownership, row/voxel
  maps, active batch strides, thread-private workspaces, output sentinels,
  cleanup, nuisance projection/recovery, and selection versus full-data refits.
  Expanded compact nuisance-projection code into readable blocks. Numerical
  library changes are comments/formatting only, with identical C token streams.
- Changes are confined to Picard's new source, numerical modules, tests and
  this plan. Existing AFNI readers and shared build registrations are untouched.

Validation: CMake and classic Make builds; Picard integration against both
binaries; fixed/repeat, prediction, candidate-selection and nuisance-projection
unit tests. Loader fixtures compare direct NIfTI with `3dcopy` HEAD/BRIK for
integer scaling, big-endian compressed float64 conversion, nonzero intercepts,
NaN/Inf repair and masks. A native NaN fixture verifies that `AFNI_FLOATSCAN=YES`
loads repaired data and `NO` allows Picard's finite-value rejection to act.
Edited C sources remain within 80 columns and contain no `goto` statements.

## Phase 1 implementation record — 2026-11-06

User authorization: implement phase 1, recap its purpose, and prompt before
each next phase. Avoid `goto` statements in this program and its libraries.
Work began on branch `MakeItSo` after verifying commit `cd5d7a737`
(`refactored out the goto statements`), no tracked modifications, and no
`goto` tokens in `3dRSA.c` or `thd_mvpa.c`.

Completed:

- Added `thd_regression.[ch]`, moving ordinary/per-response ridge and all four
  continuous estimators (elastic net, banded ridge, PLS, SVR) from MVPA.
  Existing signatures, configuration ownership, intercept rules, parameter
  units, diagnostics, float I/O and double arithmetic remain intact.
- Added checked `THD_ridge_get_state` / `THD_regression_get_state` accessors.
  `THD_regression_state` exposes borrowed const coefficient/centering/response
  mean arrays, active dimensions and applicable feature groups. Views expire
  on every fit attempt or free; failed requests clear the output view.
  For SVR, predictor offsets are zero and the training response mean remains
  distinct from the fitted intercept. Original-space intercepts continue to
  use the existing coefficient-export APIs.
- Updated CV baselines, coefficient/group audits, activation patterns and
  frozen-model copying to use checked state views and public diagnostics.
  No private solver workspace fields remain exposed to MVPA.
- Added the private `thd_regression_internal.h` SPD helper boundary shared by
  LDA and banded ridge. Arithmetic and legacy failure messages are preserved,
  including the existing LDA-worded banded factorization failure message.
- Re-exported the numerical API from `thd_mvpa.h`; updated CMake and classic
  Make source/object/header dependencies. The standalone regression unit
  target links only the new numerical source and libm.
- Added `pmolfese/tests/test_regression.c` for state/ownership/precision,
  invalidation/refit, penalty units, SVR mean/intercept, and SPD contracts.
  Existing MVPA numerical and integration fixtures remain in place.

Validation:

- Captured the pre-extraction executable and sources; all 23 existing MVPA
  CTest cases passed before editing.
- CMake build and `ctest --test-dir build --output-on-failure
  -R '3dMVPA|thd_mvpa|thd_regression' -j 4`: all 24 cases passed afterward.
  Coverage includes all moved estimators, LDA/classifiers, nested tuning,
  nuisance handling, coefficients, response-specific ridge, saved models,
  activation patterns, volume/surface workflows and worker determinism.
- `make -C src -j 4 3dMVPA` passed with the existing platform Makefile;
  regression-oracle and saved-model integration suites also passed against
  the resulting `src/3dMVPA` executable.
- Standalone C99 compilation with `-Wall -Wextra -Werror -pedantic` passed;
  standalone AddressSanitizer/UndefinedBehaviorSanitizer tests passed.
- Direct reference/current CMake comparisons matched 357 regression, 168
  saved-model-suite and 166 activation-pattern text artifacts (these counts
  include fixture tables). Comment lines were excluded from table comparison;
  all 23 saved-model files were compared byte-for-byte, without filtering.
- Solver-body comparison confirmed no changes beyond helper routing, local
  error-helper naming, comments and whitespace. Symbol inspection confirmed
  solver definitions reside in the new object and MVPA references them.
  No `goto` statements were introduced; `git diff --check` passed.

Local reference binaries, logs and the comparison harness are in
`/tmp/picard_phase1_reference` (temporary validation evidence, not a committed
fixture archive). The existing integration scripts retain their independent
numerical oracles.

Purpose: Picard can now share the tested estimator implementation without
depending on MVPA's classifiers, CV workflow, permutation or saved-model code.
Phase 1 introduces no `3dPicard` executable or fractional-ridge algorithm.
The next phase is stage 2: prepared SVD/fractional ridge, explicit nuisance
projection/recovery, a numeric Picard workflow and its executable adapter.
The user subsequently authorized phase 2; its completion is recorded below.
`3dDeconvolve` is unchanged.

The user reports the 3dMVPA libraries complete, and its roadmap records elastic
net (#28), banded ridge (#35), PLS (#30), and SVR (#31) as completed. The
2026-10-30 source review below identifies concrete extraction boundaries in
`thd_mvpa.[ch]`; completion of those estimators does not mean they have already
been split into standalone libraries. This update changes the plan only and
does not rerun or recertify the MVPA validation results.

## Phase 2 implementation record — 2026-11-06

User direction and standing boundaries:

- The scientific roadmap targets GLMsingle-style single-trial estimation;
  “single” in the user's clarification meant **GLMsingle**. Fractional ridge
  alone is still not GLMsingle equivalence. Phase 3b is the natural next
  milestone for this goal, followed by candidate-HRF and noise-pool work.
- Default to AFNI HEAD/BRIK buckets and minimize output files. Keep explicit
  NIfTI export available. Make provenance part of the dataset header.
- **Do not modify existing AFNI files that the user did not author without
  express permission.** The user explicitly approved retaining scoped build/test
  registration edits in `src/CMakeLists_binaries.txt`, `src/Makefile.INCLUDE`,
  and `tests/CMakeLists.txt` (“Allow these scoped build/test edits”). This is
  not permission to change other AFNI files or existing library implementations.
- Avoid `goto` statements in Picard and its new libraries. Recap each completed
  phase and prompt before starting the next one.

Completed implementation:

- Added `thd_regression_svd.c` to the new regression module, with an immutable
  prepared-double design, per-worker response scratch, fixed/per-response
  fractional fits and a streaming path callback. This separate translation
  unit keeps the existing MVPA regression kernels and their link requirements
  unchanged. Both `thd_mvpa.[ch]` and `thd_regression.c` match the end of phase 1.
- Added `thd_linproj.[ch]`: SVD projection of exactly the supplied nuisance
  span, rank reporting, and minimum-norm coefficient recovery. Added
  `thd_picard.[ch]`: owned numeric designs/role/row/run maps, projected fitting,
  nuisance recovery from Y-XB, and original-unit fitted/residual batches.
- Added `3dPicard.c`: existing xmat ingestion, explicit/default column roles,
  validate/show-design mode, full/retained timeline checks, volume/mask input,
  bounded response batches, private OpenMP workers, memory-budget checks,
  labelled AFNI outputs and header provenance. It does not require MVPA,
  cross-validation, inference, or Maker modules.
- Reused the existing EISPACK SVD, xmat reader and dataset interfaces without
  modifying their implementations. The initial raw NIfTI adapter and forced
  `AFNI_FLOATSCAN=NO` setting were removed in the post-phase-5 I/O revision
  below; all formats now use the standard AFNI dataset loader.
- Added `test_linproj.c`, `test_picard.c`, `test_picard.py`, and the pinned
  numeric fixture `picard_fracridge_reference.json`.

Numerical and CLI decisions frozen for this phase:

- Preserve supplied design scaling. Add no intercept, baseline, whitening,
  standardization, or tuning. Default `ColumnGroups>0` to penalized columns
  and `<=0` to nuisance; explicit column lists override roles. Missing or
  malformed default role metadata requires explicit roles. Every column has
  exactly one role; duplicate labels and ambiguous censor maps are rejected.
- Rank retains singular values `s > rank_tol*smax`, default `rank_tol=1e-10`.
  Use a direct thin SVD, including the tall transpose for wide designs, rather
  than normal equations. Validate SVD reconstruction before accepting factors.
  Projection uses a second pass and sets residual vectors to numerical zero
  at `64*DBL_EPSILON*max(n,k)*input_norm`; this policy is explicit in the help.
- Define fraction in supplied penalized-coefficient units relative to the
  minimum-norm projected LS reference. Solve its monotone norm equation by
  bounded bisection to relative fraction tolerance `1e-10`. Fraction 1 is the
  reference; fraction 0 is zero penalized coefficients and infinite penalty.
  A zero reference yields zero betas, lambda=0 and undefined achieved fraction.
- Report normalized lambda for `SSE/(2n)+lambda*||beta||^2/2`, alongside the
  diagonal `PenaltyXtX=n*lambda`. Public diagnostics use infinity/NaN for the
  mathematical endpoint/undefined ratio. Output maps encode those values as
  -1 and distinguish normal (0), zero reference (1), zero endpoint (2), and
  outside-mask (-1) statuses. Failed fits invalidate results and write no maps.
- Require `-fraction` for fitting. `-validate_only` / `-show_design` may run
  without image input. Full input time length must equal NRowFull. Data are
  real scalar volume series; surface data and dataset selectors are deferred.
- Default output is **one** `-prefix` HEAD/BRIK bucket (compressed BRIK follows
  AFNI's normal settings), containing penalized `*_Coef` bricks, optional
  nuisance `*_Coef` bricks (`-save_nuisance`), then seven `Picard_*` diagnostics.
  Explicit `.nii`/`.nii.gz` prefixes request NIfTI. `-fitts PREFIX` and
  `-errts PREFIX` request optional companion buckets; `-save_fitted` and
  `-save_residuals` generate `_fitts` and `_errts` names automatically.
- Store model/settings, labels, column roles, matrix/input/mask sources,
  stimulus metadata, TR, ranks and full/retained row/run identities in
  `PICARD_*` attributes and history. Fitted/residual buckets contain retained
  rows only, labelled by original time; their header mapping avoids inventing
  predictions at censored rows. `-write_manifest` explicitly requests the
  beta TSV and retained/full/run `.rows.1D` sidecars. No text files by default.
- Default worker batch=64, jobs=1, memory budget=2048 MB. The budget checks a
  conservative estimate of resident input/output payload and design/worker
  buffers; it is not a hard process-RSS limit and excludes library bookkeeping
  and compression buffers. No candidate-by-voxel coefficient tensor is kept.

Validation results:

- CMake and platform Make builds of `3dPicard` passed. Both binaries passed the
  same integration/oracle script. Both new C unit tests passed; AddressSanitizer
  and UndefinedBehaviorSanitizer checks passed for the new numerical workflow
  and nuisance projection tests. Existing 24 MVPA/regression tests passed.
- Actual unchanged `3dDeconvolve` output fixtures cover censoring, two runs,
  baseline roles and labels. Independent NumPy/LAPACK augmented least squares
  and SciPy Brent-root oracles check betas, nuisance recovery, predictions,
  residuals, norm ratios and lambda units. C fixtures check streamed paths,
  per-response fractions, copied state, failed-fit invalidation and rank cutoff.
- Tall/wide/rank-deficient designs, redundant nuisance columns, zero references,
  both fraction endpoints, preserved supplied units, masks and malformed inputs
  pass. Tests also cover native AFNI input, integer scaling, big-endian/gzipped
  double NIfTI with nonzero scaling intercept, and nonfinite retained inputs.
- Batch sizes 2/5 and jobs 1/2 produce identical output values. Default output
  is only a HEAD/BRIK pair; coefficient labels, header maps, optional manifests,
  explicit NIfTI output, familiar fitts/errts prefixes, collision protection,
  failed-output cleanup and the memory-budget rejection are checked.
- Compared against [nrdg/fracridge at ada68ca](https://github.com/nrdg/fracridge/blob/ada68ca27c881b000d704961e4cafe085e809271/fracridge/fracridge.py).
  The saved fixture records the full commit, source SHA-256 and parameters.
  That implementation interpolates a log-alpha grid; Picard solves the norm
  equation directly. Reference comparisons allow 2% relative / 0.003 absolute
  difference, while independent direct-model output checks use 3e-5 relative /
  3e-6 absolute tolerances for float maps. Tests require no network access.
- No `goto` statements were introduced. Existing AFNI library implementations
  and 3dDeconvolve remain untouched. `git diff --check` passed. Temporary build
  and validation logs are in `/tmp/picard_phase2_reference`.

Example (native AFNI output):

```sh
3dPicard -input epi+orig -matrix X.xmat.1D -fraction 0.5 \
         -prefix stats.picard -errts errts.picard
```

Build/run the new tests from the repository root:

```sh
cmake --build build --target 3dPicard thd_linproj_unit thd_picard_unit -j 4
ctest --test-dir build --output-on-failure -R '3dPicard|thd_linproj|thd_picard'
```

The Python oracle script requires NumPy, SciPy and nibabel, plus existing
3dDeconvolve/3dinfo/3dcopy fixture tools. Purpose: phase 2 now performs fixed
voxelwise fractional ridge with unpenalized nuisance terms on ordinary AFNI
matrices and produces native AFNI buckets. Stop and prompt before phase 3b;
repeated-condition tuning, HRF selection and noise-pool estimation are not yet
implemented.

## Phase 3b implementation record — 2026-11-06

User authorization: implement Picard phase 3b and improve source readability.
The user clarified that **Maker M0 and M1 are not included**. They remain
planned; completing this phase satisfies the recommended Picard prerequisite
for beginning Maker, not Maker's own acceptance gates.

Completed implementation:

- Extended only the new `thd_picard.[ch]` numerical module with immutable
  repeated-condition plans, private worker buffers, per-fold candidate losses,
  response-specific selection and final full-design refitting. Existing
  regression/SVD/projection implementations and AFNI libraries are unchanged.
- Added `-fractions LIST` with `-trial_info FILE` as an alternative to fixed
  `-fraction F`. The whitespace-delimited table has four nonnegative integer
  fields: `column trial condition run`. Column and run indices are zero-based;
  run indices refer to xmat `RunStart`. Comments begin with `#`. Table order is
  arbitrary, but every penalized column must occur exactly once, trial IDs must
  be unique, and one column per trial is supported. Multi-basis trials are
  explicitly outside this phase.
- Require each trial column to have nonzero retained support confined to its
  declared run. Require nuisance columns to be run-local; reject shared
  cross-run nuisance vectors rather than coupling training and held-out fits.
  All-zero nuisance columns are harmless. Each retained run needs at least one
  condition repeated in another run. Unmatched conditions remain in fits but
  do not contribute to selection. Designs without repeats can still use the
  existing fixed-fraction mode.
- Construct one leave-one-run-out training design per retained run. Remove the
  omitted run's trial/nuisance columns and factor only training rows. Fractions
  use each training fit's own minimum-norm reference, without held-out response
  values in the penalty calculation. Full-design OLS supplies held-out targets:
  the required run-local support makes this OLS problem block-separable.
- Compare **every** training-trial beta with every matching-condition held-out
  OLS beta. Score ordered pairs, not condition averages. Pool squared errors
  across pairs and divide by total pair count, so unequal folds receive pair
  weights. Scores are MSE in supplied coefficient units, without session
  normalization. Users must supply comparable trial scaling/response units.
- Choose a fraction separately per voxel; exact ties choose the largest
  fraction. Zero responses therefore choose the largest grid value. Always
  compute an OLS selection baseline, even if fraction 1 is absent from the
  grid. Refit the complete supplied design at the selected fraction and report
  that fit's lambda, nuisance coefficients, fitted values and residuals. No
  post-fit beta scaling/offset is applied. These are selection losses, not
  independently evaluated prediction scores, R2, or inferential statistics.
- Retain default **one HEAD/BRIK bucket** output. Tuning adds `Picard_CV_MSE`,
  `Picard_CV_OLS_MSE`, and `Picard_CV_Index`; `Picard_Fraction` contains the
  selected fraction. Optional `-save_cv` appends pooled MSE maps named
  `Picard_CV_Grid000`, etc., in the supplied grid order. CV maps use -1 outside
  the mask. No extra score files are created.
- Store grid order, column/trial/condition/run identities, fold run IDs and
  pair counts, scoring units and the fold-local policy in `PICARD_*` header
  attributes. Optional `-write_manifest` adds Trial/Condition/Run columns to
  the existing beta TSV. Include training plans and worker loss buffers in the
  conservative pre-allocation memory estimate; never allocate an all-voxel
  candidate coefficient tensor. Candidate paths are streamed per worker batch.
- Reformatted all of `3dPicard.c` within 80 columns. Added ownership/pipeline
  documentation, comments on every function, descriptive stage names, named
  output slots and expanded control flow. Explained matrix storage conversion,
  censor identities, scaling, threading, borrowed numerical results, and output
  rollback. Preserved the Picard dedication and Trek references in help.
  Also formatted/commented the new numerical tuning API and implementation.

Scientific comparison and validation:

- Pinned GLMsingle scoring reference:
  [`calcbadness.py` at 082851f9db87ff30c49370f27365ac1112f18ff7](https://github.com/cvnlab/GLMsingle/blob/082851f9db87ff30c49370f27365ac1112f18ff7/glmsingle/ssq/calcbadness.py).
  `picard_glmsingle_reference.json` records source SHA-256 values and numeric
  outputs from the unmodified helper plus `zerodiv.py`.
  `make_picard_glmsingle_fixture.py` regenerates it from downloaded pinned
  sources; normal tests require no network or GLMsingle installation.
- This is a **component comparison**, not a GLMsingle-equivalence claim.
  GLMsingle evaluates a full-data fraction path and session-normalizes betas;
  Picard uses fold-local fraction fits and raw coefficient-unit MSE. The pinned
  fixture uses equal-norm orthogonal trials, where both fraction policies agree,
  and one session. Convert Picard MSE by `pair_count/reference_sample_variance`
  to match GLMsingle's normalized SSE. The scored values and winners agree.
  HRF selection, noise-pool PCs and autoscaling are not implemented here.
- Independent NumPy/SciPy augmented-least-squares and Brent-root oracles check
  general correlated designs separately, including training-fold nuisance
  fits, all grid scores, selected fractions, full refit betas and residuals.
  Perturbing one run strongly exercises training isolation. Real unchanged
  `3dDeconvolve -stim_times_IM -local_times` fixtures cover unequal run lengths,
  explicit conditions, multiple repeats within a run and censoring.
- Hand-computable C fixture checks pair counts/MSE, per-response winners,
  endpoint refits, zero-response ties, metadata ownership, smaller subsequent
  batches and result invalidation. Integration checks include rank-deficient
  trial designs, redundant nuisance, grid reordering, missing fraction 1,
  malformed/missing/duplicate identities, insufficient repeats, wrong run
  support, shared nuisance rejection, zero-support trials and memory rejection.
- Both CMake and classic Make builds pass. The three Picard CTest cases pass;
  the integration/oracle suite also passes against the classic executable.
  HEAD/BRIK defaults, explicit NIfTI export, score/identity metadata, masks,
  optional manifests, rollback, and batch/thread equality are covered. Existing
  fixed-fraction tests continue to pass. Strict C99 warning checks and
  AddressSanitizer/UndefinedBehaviorSanitizer checks pass for the numerical core.
- No `goto` statements were introduced. No existing AFNI source, numerical
  library, or build/test registration file needed modification in this phase.
  `thd_crossval`'s table/sample-based search API does not fit this numeric
  beta-pair traversal; no dependency or speculative generic API was added.
  Temporary references and validation logs: `/tmp/picard_phase3b_reference`.

Example table (assuming columns 0/1 are run-local nuisance columns, and
columns 2–5 are the four single-trial regressors):

```text
# column trial condition run
2        100   0         0
3        101   1         0
4        102   0         1
5        103   1         1
```

```sh
3dPicard -input epi+orig -matrix X.xmat.1D \
         -fractions '1,.8,.6,.4,.2,0' -trial_info trials.1D \
         -prefix stats.picard -save_cv
```

Purpose: Picard can now select shrinkage from repeated-condition evidence
instead of requiring one fixed fraction for every voxel, while retaining
AFNI-native outputs and explicit trial identities. Stop and prompt for the
next phase. Candidate HRF work (Picard stage 5) and independent Maker M0/M1
remain separate possible next steps; neither is authorized by this completion.

## Stage 3a implementation record — 2026-11-06

The user requested completing stage 3a if it was missing, and phase 4 only if
3a had already been implemented. Inspection confirmed 3a was missing. This
change completes 3a; phase 4 and Maker M0/M1 remain unimplemented.

Completed implementation:

- Added numeric time-series fold construction, deep copying and generic batch
  path traversal to the user-authored `thd_crossval.[ch]` (authorship checked
  against the file history). Existing `.c` implementation bytes, scalar search,
  nested-fold semantics and validators are preserved. The `THD_foldplan` struct
  now has a tag for forward declarations; its layout is unchanged.
- `THD_foldplan_time_series` uses explicit original row/run identities. Default
  folds leave retained runs out; optional blocks start at each original
  `RunStart`. An explicit assignment vector can define other full-coverage
  partitions. Temporal gaps exclude nearby training rows **within the same
  run**, measured on original acquisition coordinates before censoring. Every
  retained row is tested exactly once; empty training/test folds fail. Copying
  a plan preserves excluded rows, not merely its test assignments.
- `THD_cv_path_search` calls response preparation and path evaluation once per
  fold/batch, followed by caller-owned pooling/selection. It performs no fit-time
  allocation. Prepared factors are immutable, created before response data are
  seen, and shared across workers. Existing scalar/nested search is unchanged.
- Added `thd_picard_predict.[ch]` for shared-feature prediction: train-only
  nuisance projection, fractional reference norms, candidate fitting,
  test-row scoring, deterministic selection and separately prepared full-design
  refitting. `thd_picard_internal.h` holds the unchanged private plan layout so
  workflows can share numeric storage without exposing it to the application.
  Every penalized column needs nonzero training support in every fold;
  trial-specific columns must use the separate repeated-condition mode.
- Keep AFNI I/O, masks and scheduling in `3dPicard.c`. Global normalization and
  accumulation also live in the numeric workflow, through an explicit bounded
  accumulator API. The existing crossval/sample-table headers/link dependencies
  are reused without a broad library refactor; prediction never reads datasets
  or constructs sample tables. Existing AFNI numerical libraries are unchanged.

Scientific contract and CLI:

- `-cv_predict -fractions LIST` enables prediction tuning and is mutually
  exclusive with `-trial_info` and fixed `-fraction`. Default folds leave one
  run out. `-cv_block_size N` requests run-local blocks; `-cv_folds FILE` instead
  supplies two nonnegative integer fields per retained row: original full-row
  index and fold ID. File order is immaterial; folds follow first retained-row
  occurrence of each ID. `#` comments are allowed. `-cv_gap N` defaults to zero.
- `-cv_nuisance train` (default) predicts `X_test beta + Z_test gamma` with
  training coefficients only. An SVD row-space check rejects test nuisance
  vectors unidentifiable from training (such as an unseen run intercept),
  while accepting redundant shared nuisance columns.
- `-cv_nuisance project` explicitly scores `P(Z_test)Y_test` against
  `P(Z_test)X_test beta`. This is **conditional signal prediction**: held-out
  responses are used to remove a supplied nuisance span, but never to learn
  beta, lambda, scaling or new nuisance regressors. Reject a nuisance span
  that removes every test degree of freedom. No implicit intercept, centering,
  predictor standardization, whitening or learned noise-pool terms are added.
- Pool SSE by test-row count; MSE is pooled SSE divided by retained rows.
  The baseline is each training response mean, passed through the same test
  projection. Report `R2 = 1 - SSE/baseline_SSE`; retain negative R2. Exactly
  zero baseline SSE makes R2 undefined (public API NaN; dataset R2=0 with
  `Picard_PredStatus=1`). MSE remains defined and selectable. Normal status is
  0 and outside-mask status is -1. These are **selection scores**, not outer-test
  performance or inferential statistics. Nested evaluation is not advertised
  or implemented in this stage.
- `-cv_rule min` selects the lowest pooled MSE, resolving exact ties toward the
  largest fraction. `-cv_rule 1se` selects the smallest fraction within the best
  candidate's pooled mean plus its weighted-fold standard error. With
  `w_f=ntest_f/n`, use
  `SE^2 = sum(w_f*(MSE_f-mean)^2)*sum(w_f^2)/(1-sum(w_f^2))`.
  This is an explicitly defined selection heuristic, not a confidence interval
  for overlapping training folds. Full-data refitting uses the winning fraction
  and its own full-design norm/lambda, not coefficients borrowed from a CV fold.
- `-cv_global` requires an explicit `-mask`. Normalize each included voxel's
  fold losses by its pooled baseline MSE, then average across voxels with
  positive baselines. Exclude zero-baseline voxels; if all are excluded, choose
  the largest fraction. Use the same minimum/one-SE selector. A serial prepass
  accumulates voxels in spatial order across batches; the final parallel pass
  refits at the common fraction. This is bounded in memory and deterministic
  across batch/job choices, and does not promise faster global fitting.
- Preserve the default single AFNI HEAD/BRIK bucket. Prediction mode adds
  `Picard_PredMSE`, `Picard_PredBaseline`, `Picard_PredR2`, `Picard_PredIndex`, and
  `Picard_PredStatus`, instead of repeat-tuning diagnostics. Optional `-save_cv`
  appends candidate MSE maps to the same bucket. `PICARD_PRED_*` attributes
  record exact train/test original-row membership (including gaps), scoring
  policy, rule and any global selection/normalization. The existing optional
  fitted/residual buckets remain **full-fit retained-row** outputs, not held-out
  predictions. No extra score sidecars are created.
- Expanded help, function comments and memory estimates; retained 80-column
  source/help readability and the Picard dedication. Added scoped registration
  changes only in the three previously approved CMake/Make/test files.

Validation:

- Both CMake and classic Make builds pass. The combined 29 selected CTest cases
  (Picard, crossval and MVPA/regression) pass. After the global-accumulator
  boundary refinement, the five targeted Picard/crossval cases pass again.
  The complete integration/oracle suite also passes against the classic binary.
- New `test_picard_predict.c` has hand-computable MSE, R2 and full-refit values;
  held-out perturbation verifies the training fractional fit remains unchanged.
  It checks numeric block/gap coordinates, copied fold ownership, generic path
  abort behavior, unseen nuisance rejection, explicit conditional projection,
  weighted one-SE selection, zero-baseline handling, forced/global selection,
  accumulator normalization and failed-fit/accumulator invalidation.
- Independent NumPy/SciPy direct-fit oracles check every candidate's prediction
  MSE, training-mean baseline, selected fraction, R2 and full refit. Cases cover
  unequal runs, censored rows, large held-out perturbations, correlated and
  rank-deficient predictors, redundant nuisance, original-time blocks/gaps,
  shuffled explicit fold tables, multiscale voxel responses, minimum/one-SE and
  per-voxel/global selection, all-undefined global fallback and masks. Real
  unchanged `3dDeconvolve` matrices exercise run-local baseline projection.
- Batch sizes and job counts produce identical arrays. Tests reject malformed
  folds, missing rows, missing explicit global masks, incompatible options,
  empty training after gaps, memory limits, unknown policies and trial columns.
  Existing fixed-fraction and repeated-condition tests remain in the suite.
- AddressSanitizer/UndefinedBehaviorSanitizer pass for the prediction/fold unit
  executable, including the numeric global accumulator. Strict C99 compilation
  passes with warnings treated as errors for the new workflow (unused static
  definitions in existing AFNI headers are suppressed).
- No `goto` statements were introduced. Existing crossval function bodies,
  `3dDeconvolve`, regression/SVD/projection implementations and unrelated source
  files remain untouched. Temporary build/test evidence is in
  `/tmp/picard_phase3a_reference`. `git diff --check` passes.

Example for ordinary AFNI shared-stimulus matrices with run-local baselines:

```sh
3dPicard -input epi+orig -matrix X.xmat.1D \
         -cv_predict -fractions '1,.8,.6,.4,.2,0' \
         -cv_nuisance project -cv_rule min -prefix stats.predict
```

Purpose: Picard now supports both distinct tuning questions: predicting new
observations from shared regressors (3a), and reproducing condition responses
across trial-specific regressors (3b). Stop and prompt before stage 4, which
adds prepared-double elastic-net fitting and warm-start paths. Maker milestones
remain separate; this phase does not complete or authorize them.

## Phase 4 implementation record — 2026-11-06

User authorization: "go for phase 4". Scope is prepared-double lasso/elastic
net and their integration into Picard's existing fitting/tuning workflows.
No stage 5 or Maker milestone is included. Existing AFNI library implementations
outside the user's modules were not changed; the previously approved scoped
build/test registration edits add the shared regression source to Picard targets.
No `goto` statements were added.

Completed numerical contract:

- `thd_regression.c` now has one shared coordinate-descent kernel used by the
  original MVPA wrapper and new `THD_elastic_*` prepared-design/worker/path APIs.
  The wrapper retains MVPA centering, cold starts, parameter units and state
  layout. No second independent elastic-net implementation was introduced.
- Double-precision objective: `SSE/(2*n) + lambda *
  (rho*sum(abs(theta)) + (1-rho)*sum(theta^2)/2)`. Lambda is positive and
  finite; rho is in [0,1]. rho=1 is lasso; rho=0 uses the same normalized ridge
  objective with diagonal penalty `n*lambda`. Zero lambda is deliberately not
  part of this iterative API: fractional ridge at fraction=1 supplies the
  unregularized minimum-norm reference.
- Project X and Y against explicit unpenalized Z first. Preserve supplied
  predictor units by default. Optional standardization divides projected X
  columns by `sqrt(sum(X^2)/n)` (zero columns keep scale 1), with no additional
  centering or implicit intercept. Never scale Y. Export beta in original
  column units, then recover nuisance gamma from original Y-X*beta.
- Each immutable training plan prepares its own projection, scales and column
  norms. Workers own response/coefficient/residual scratch. Positive decreasing
  lambda paths warm-start within the same fold and response batch only; final
  selected fits start cold. No fitting-time allocation or all-voxel coefficient
  path tensor. Full SVDs remain available for rank and OLS-reference work;
  sparse fits additionally retain their effective design.
- Every coordinate sweep recomputes residuals and checks all KKT conditions.
  There is no active-set screening to bypass these checks. Nonfinite data,
  nonconvergence and callback rejection invalidate results/warm state and abort.
  Diagnostics include nonzero count, sweeps, maximum KKT residual, objective,
  and response-specific lambda_max in effective units. Pure ridge has no finite
  lambda_max (NaN in the API, -1 in the output map).

Executable/workflow:

- Added `-penalty ridge|lasso|enet` (default ridge), `-lambda`, `-lambdas`,
  `-l1_ratio` (enet default .5; lasso fixes 1), `-standardize`, `-solver_tol`
  (default 1e-8 absolute max KKT), and `-max_iterations` (default 10000).
  Exact option combinations and grid order are validated. Mixing is one
  explicit fixed choice per analysis; phase 4 does not jointly tune rho.
- Both shared-predictor CV and repeated-condition beta tuning accept lambda
  grids through the same existing scoring/fold machinery. Sparse minimum-loss
  ties choose the smallest lambda; prediction one-SE chooses the largest
  eligible lambda. All-undefined global selection also chooses the smallest
  lambda. OLS targets and reference loss remain explicitly unregularized.
  Existing fraction behavior and ridge bucket layout are preserved.
- Native HEAD/BRIK buckets remain the default. Sparse fits use ten labeled
  diagnostic sub-bricks: Lambda, L1Ratio, Nonzero, Iterations, KKT, Objective,
  LambdaMax, Rank, Standardized, Status. CV maps follow; `-save_cv` appends
  grid losses in the same bucket. Optional nuisance/fitted/residual products
  retain their existing roles. Nonconvergence writes no result datasets.
- `PICARD_ESTIMATOR`, `PICARD_LAMBDAS`, numerical-policy text, original row/run
  identities and existing fold metadata distinguish sparse fits from fractions.
  Help/source remain within 80 columns, with comments explaining preprocessing,
  ownership, optimality checks, warm starts and parameter direction.

Validation:

- Independent scikit-learn 1.9.1 `ElasticNet` fits with matched objective and
  `fit_intercept=False`, plus a direct linear solve for rho=0. Checked betas,
  nuisance recovery, fitted/residual series, objectives, lambda_max and KKT.
- Both tuning modes, fold-local scaling, unequal/censored runs, global and
  voxelwise selection, both prediction selection rules/nuisance policies,
  zero-response/global fallback, temporal gaps, batch/thread equality, default
  buckets/provenance and failure rollback. Wide/duplicate-column lasso compares
  fitted values/objectives rather than assuming unique coefficients. An
  unchanged 3dDeconvolve-produced matrix is also exercised in sparse mode.
- Hand-computable constant-predictor and nuisance-projection fits; warm/cold
  path agreement at rho=0,.5,1; zero columns/responses; invalid warm starts,
  changed batch size, failed-result invalidation and callback aborts.
- All 29 selected regression/MVPA/Picard tests passed. CMake and classic builds
  passed; address/undefined-behavior sanitizers and strict C99 warning checks
  passed on the solver and numerical workflow. Python integration tests now
  require scikit-learn in addition to NumPy/SciPy/nibabel. Dependencies were
  installed only in the build virtual environment for validation.

Examples:

```sh
3dPicard -input epi+orig -matrix X.xmat.1D -penalty lasso \
         -lambda .1 -standardize -prefix stats.lasso

3dPicard -input epi+orig -matrix X.xmat.1D -penalty enet \
         -l1_ratio .5 -lambdas '1,.3,.1,.03' -standardize \
         -cv_predict -cv_nuisance project -prefix stats.enet
```

Purpose: Picard can now select exact-zero coefficients or combine sparsity
with continuous shrinkage, while retaining explicit nuisance treatment,
original-unit outputs and the two distinct tuning questions. This does not
claim better results than ridge for every dataset, add penalized t/p maps,
or implement a full GLMsingle pipeline. Stop and prompt before phase 5:
candidate-design identity validation and HRF-selection policy. Maker remains
an independent, unstarted implementation roadmap.

## Phase 5 implementation record — 2026-11-06

User authorization: "go forth on phase 5", followed by approval of both
policies: HRF choice by unregularized full-data fit and general design choice
by held-out prediction. All scores remain explicitly selection scores, with
no independent outer-test performance claim. No phase 6 or Maker work is
included. Existing AFNI library implementations outside the user's modules
remain unchanged; scoped build/test registrations use the prior approval.
No `goto` statements were added.

Completed:

- New `thd_picard_candidates.[ch]` owns a numeric candidate-selection workflow
  over borrowed immutable fixed plans. It constructs its own ordinary Picard
  tuning plans; it does not implement another solver, fold engine or HRF
  generator. Plans validate compatible dimensions, roles, original row/run
  IDs, nuisance values and estimator/scaling settings. Worker results keep
  one winning coefficient/series set per response plus bounded scalar score
  grids. HRF final refits compact only responses that selected that design.
- New `picard_manifest.[ch]` owns versioned TSV/path and semantic identity
  validation at the file-adapter boundary. It checks exact column coverage,
  labels/roles, trial/condition/run/basis/item identities, declared run support,
  row/censor/run/TR metadata, stimulus labels/ranges and common duration and
  normalization declarations. Candidate-specific HRF formula/source metadata
  and penalized values may differ. Different basis dimensions are rejected.
- `-candidate_manifest FILE` replaces `-matrix` only for the collection mode;
  `-candidate_policy hrf|predict` is required. Ordinary single-matrix fitting
  and all phase 1-4 behaviors remain supported without a manifest or Maker.
  `-validate_only` reports each candidate's index, stable ID, rank and path.
- Stable candidate IDs are unique integers in [0,16777215], exactly representable
  in AFNI float sub-bricks. Exact ties select the smallest stable ID, independent
  of manifest order. Index maps remain explicitly zero-based manifest positions.

Scientific policies:

1. **HRF:** choose the candidate with minimum full-data unregularized residual
   MSE. With identical observations/nuisance terms and a common denominator,
   this agrees with choosing maximum fit R2. Require equal projected ranks
   rather than comparing different fitting complexities. Then apply the fixed
   fraction/lambda or repeated-condition tuning to the selected design. HRF
   choice is independent of the final penalty and its strength. Subsequent
   beta tuning is explicitly **conditional on full-data HRF selection**; it
   is not independent validation of that selection. This follows the HRF
   selection principle in [GLMsingle's implementation](https://github.com/cvnlab/GLMsingle/blob/main/glmsingle/glmsingle.py),
   not an end-to-end GLMsingle-equivalence claim.
2. **Prediction:** each design uses the same supplied folds and nuisance
   scoring policy but independently prepared training factors/scales. Apply
   its existing min/one-SE parameter rule, then choose the design with the
   smallest selected prediction MSE. One-SE is within a design; there is no
   assumed strong/weak ordering of HRFs. Validation data select both parameter
   and design, so MSE/R2 are selection scores, not outer-test performance.

Collection selection is voxelwise. `-cv_global` is explicitly rejected for
collections; a joint global design/parameter rule has not been defined.
HRF policy rejects `-cv_predict`, avoiding an implicit full-data HRF choice
masquerading as fold-local prediction. Prediction policy requires its existing
parameter grid and `-cv_predict`. Independent nested evaluation remains future
work and is neither implemented nor advertised by this phase.

### Candidate manifest v1: implemented consumer contract

Files are literal, unquoted TSV: spaces in paths/labels are preserved; tabs
separate fields; there is no shell expansion or escape processing. Relative
matrix/column-file paths resolve against the manifest directory and are
recorded as absolute paths. The first line is exactly
`# 3dPicard candidates 1`; the next noncomment line is the eight-field header.
Two through 4096 candidates are accepted. Blank lines and subsequent `#`
comments are allowed; empty fields, duplicate IDs, unsupported versions and
malformed rows are errors. Use `-` for absent provenance in prediction mode;
HRF mode requires explicit source/version/sampling/duration/normalization.

```tsv
# 3dPicard candidates 1
id	matrix	columns	source	version	sampling	duration	normalization
7	fast.xmat.1D	fast.columns.tsv	custom-HRFs	v1	0.1	boxcar:1s	HRF-peak:1
12	slow.xmat.1D	slow.columns.tsv	custom-HRFs	v1	0.1	boxcar:1s	HRF-peak:1
```

`sampling` records the source HRF sample interval in seconds, not the matrix
RowTR. Duration/normalization strings declare how already supplied matrices
were constructed and must agree across candidates. The fitter does not
convolve, resample or normalize an HRF based on these fields. Validation cannot
prove that an external generator truthfully described its waveforms or units.

Each columns file has eight tab-separated fields, one row per original matrix
column, in any row order. For example (expand to cover the entire matrix):

```tsv
column	label	role	trial	condition	run	basis	item
0	intercept	nuisance	-	-	-	0	-
1	face_trial0	penalized	0	0	0	0	face01
2	face_trial1	penalized	1	0	1	0	face01
```

Column/basis and supplied trial/condition/run IDs are zero-based integers.
Role is `penalized` or `nuisance`. Labels match the xmat exactly. Trial columns
require condition and run IDs; `(trial,basis)` pairs are unique, and multiple
basis components of a trial must agree on condition/run/item. A declared run
requires support only in that run. `-trial_info` must agree and retains its
one-column-per-trial/basis=0 contract. Conditions may be declared without a
trial for shared-condition predictors. Item is an optional supplied identity.

This is the concrete consumer contract for later Maker M6/M7 integration;
Maker itself, HRF-library generation and library redistribution remain unstarted.

### AFNI outputs and downstream identities

- Default remains one HEAD/BRIK bucket. `Picard_DesignID`, `Picard_DesignIndex`
  and `Picard_DesignMSE` follow the ordinary fit/tuning diagnostics.
  `-save_candidates` appends one `Picard_DesignMSE_ID...` per candidate in
  manifest order. These maps use -1 outside the mask. `-save_cv` contains the
  selected design's parameter-grid losses, preserving existing map meanings.
- Headers record the complete candidate catalog, column identity mapping,
  source manifest and selection policy. `PICARD_DESIGN_DIAGNOSTIC_START`
  locates the ID/index/MSE maps. `PICARD_MATRIX` explicitly indicates multiple
  designs. `PICARD_DIMENSIONS` describes the first candidate; the rank map and
  candidate catalog give actual selected/per-candidate ranks.
- Collection beta sidecars retain common coefficient/trial/basis/item identities
  and provide `CandidateManifest`, `DesignIDFile`, and `DesignIDVolume`.
  They use `Column`, not a misleading single `Xmat/XColumn` pair: different
  voxels can derive the same labeled beta from different candidate matrices.
  Missing numeric identity values are -1. Existing single-matrix sidecars are
  unchanged. Fitted/residual optional outputs reconstruct each voxel using its
  own selected design.
- Response and candidate score scratch is bounded by worker batch size.
  Immutable full/fold factors and batch-bounded solver scratch for each design
  remain resident, covered by the conservative memory estimate. No
  all-design/all-voxel coefficient tensor
  or claimed automatic speedup. Invalid candidate dimensions are rejected
  before factorization using the reference-sized collection budget.

Validation:

- All 30 selected Picard/regression/MVPA tests passed; classic-build Picard also
  passed the expanded integration suite. Existing regression/core source files
  were not edited in phase 5.
- Independent NumPy/SciPy/scikit-learn oracles cover both policies, fixed ridge
  endpoints and lasso/elastic net, min/one-SE, raw/conditional prediction,
  censored/unequal runs, per-voxel winners and exact ties, manifest reordering,
  mixed selected-response batches, masks, labels, metadata, sidecars and NIfTI.
- Synthetic convolved HRFs and run-local trial designs verify HRF selection
  followed by repeated-condition tuning, aligned event identities, selected
  parameter scores and final refits. Actual GLMsingle HRF-library waveform
  generation/parity remains a later Maker/orchestrator acceptance gate.
- Rejection cases cover malformed/versioned TSV, missing/duplicate columns and
  IDs, labels/roles/basis/trial mismatches, run/censor/TR differences, nuisance
  changes, HRF rank differences, incompatible options and memory limits.
- Candidate numeric and adapter fits, invalid metadata and partial setup cleanup
  passed address/undefined-behavior sanitizers. Strict C99 warning checks
  passed. Source/help width remains <=80 columns. Test evidence is in
  `/tmp/picard_phase5_reference` for this working session.

Examples:

```sh
3dPicard -input epi+orig -candidate_manifest hrfs.tsv \
         -candidate_policy hrf -fractions '1,.8,.5' \
         -trial_info trials.1D -prefix stats.hrf -save_candidates

3dPicard -input epi+orig -candidate_manifest designs.tsv \
         -candidate_policy predict -cv_predict -fractions '1,.8,.5' \
         -cv_nuisance project -prefix stats.designs
```

Purpose: voxels can now choose among explicitly aligned HRF/design candidates
while retaining trustworthy beta identities and AFNI-native output organization.
Stop and prompt before phase 6: banded-ridge adaptation, with prewhitening and
support-only refits treated as separately validated subchanges. Maker remains
separate; phase 5 does not complete M0/M1 or implement automatic noise-pool PCA.

## Phase 6 implementation record — 2026-11-06

User authorization: "go for phase 6", followed by the explicit choice of fixed
AR(1) coefficients with run/censor/fold boundaries respected and support-only
OLS refits for lasso/enet. Phase 7 and Maker milestones remain separate.

### 6a. Prepared banded ridge

- Added prepared-double banded designs/workspaces to `thd_regression.[ch]`.
  Legacy MVPA and Picard call one primal/dual banded kernel. Legacy centering,
  unnormalized penalty units and error text are preserved. Picard adds no
  implicit centering or scaling; its objective is
  `SSE/(2*n) + sum_j lambda[group[j]]*beta[j]^2/2` with unpenalized nuisance.
  Each training fold converts normalized lambdas to `n_train*lambda`.
- `-penalty banded -band_groups FILE -band_penalties FILE`: groups are an
  explicit integer row/column with one value per original matrix column,
  nuisance=-1, predictors=0..G-1. Every group occurs in the full design.
  Penalty-table rows are candidate vectors; columns are group IDs. Values
  must be finite and positive; duplicate vectors are rejected. One row is
  fixed fitting; multiple rows require prediction or repeat tuning.
- Training folds retain original group IDs when columns disappear, including
  groups with no remaining predictors. Candidate designs retain common group
  definitions. Fits support per-voxel and global prediction tuning, repeated
  conditions, and HRF/design selection through the existing workflows.
- Vector ordering is explicit: minimum-loss ties prefer smaller mean group
  lambda, then the earlier table row; the prediction one-SE rule prefers
  larger mean group lambda among eligible candidates. This is a convention
  for comparing vectors, not a claim of componentwise dominance.
- Native bucket diagnostics: BandIndex, Objective, KKT, Rank, one LambdaGroup
  map per group, Status, then existing tuning/design diagnostics. Headers
  contain the complete group/penalty tables and selection policies.

### 6b. Fixed AR(1) adapter

- `-prewhite_ar1 R`, finite `-1<R<1`, defaults to identity (zero). The user
  supplies R; Picard does not estimate noise or import REML/ARMA results.
- Transform X, nuisance Z and Y together before nuisance projection. Within
  a contiguous segment use `v[t]-R*v[t-1]`; the first sample is multiplied
  by `sqrt(1-R*R)`. Restart at run changes and nonconsecutive original row IDs,
  independently inside every training fold. No previous sample is borrowed
  from a held-out or censored row.
- Recover nuisance coefficients from `W(Y-XB)` against `WZ`. Check nuisance
  estimability using the same whitened training design/rank cutoff. Export
  fitted/residual data in original units, not whitened units.
- Prediction losses remain original-unit raw/conditional held-out errors;
  repeat losses remain original coefficient-unit errors. HRF choice compares
  original-unit residual MSE from unpenalized fits with the fixed whitening
  policy. These remain selection scores, not independent outer-test results.

### 6c. Support-only refits

- `-support_refit` requires lasso/enet with positive L1 mixing. Fit minimum-norm
  OLS on the exact nonzero sparse support in the same effective design space.
  Each training fold selects/refits its own support, and the final full-data
  fit selects/refits anew. Empty support yields nuisance-only predictions.
- Betas, nuisance coefficients, fitts and errts describe the refitted model.
  Sparse Objective/KKT/Nonzero retain their pre-refit meanings. SupportSize
  and SupportRank maps are inserted before Status. No t/F or ordinary
  post-selection inference is introduced.
- Factor one active design at a time using bounded temporary SVD allocations.
  Other fitting scratch is preallocated. The memory budget now includes
  whitening buffers, retained effective designs, banded factors and support
  scratch. Source comments explain transformations, ownership and dimensions.

### 6d. Reuse, measurements and validation

- Prepared primal designs reuse immutable X'X. Responses selecting an identical
  banded vector share one private factorization per batch; wide designs use
  the dual system. No all-voxels-by-candidates coefficient tensor is kept.
- Small local benchmark, 180 rows / 100 predictors / 256 voxels / two groups,
  three repetitions: median wall time was 0.151 s (batch 1, one worker),
  0.143 s (batch 32, one worker), 0.141 s (batch 32, two workers). Peak RSS was
  9.19, 9.56 and 10.09 MiB respectively. Output coefficients were identical.
  Startup dominates these small measurements; they do not establish a
  whole-brain speedup or general scaling claim.
- Fixed-AR(1) elastic-net benchmark, 120 rows / 20 predictors / 64 voxels,
  batch 16, three repetitions: median 0.117/0.114 s for one/two workers,
  peak RSS 7.78/7.94 MiB. With support refits: 0.120/0.114 s and 7.91/8.06 MiB.
  Coefficients were identical across worker counts for each fitting policy.
- All 30 selected MVPA/Picard/cross-validation tests pass. Added direct legacy
  versus prepared banded comparisons for primal/dual solves and copied-state/
  failed-fit contracts. Independent Python augmented-least-squares and
  scikit-learn oracles cover fixed/group-specific fits, wide/duplicate columns,
  positive/negative AR(1), empty supports, original-unit reconstruction,
  block/gap CV, one-SE/global selection, repeated-condition group remapping,
  candidate HRFs, invalid settings, and batch/thread invariance.
- CMake and classic Make builds pass; both binaries pass integration. Numeric
  and adapter address/undefined-behavior sanitizers pass, including new fits
  and failed setup cleanup. Strict C99 warning checks pass for Picard workflow
  and adapter modules. Source/help width checks and no-`goto` checks pass.
  Evidence for this working session is in `/tmp/picard_phase6_reference`.
- Only Picard/new regression source, headers, tests and this plan changed in
  this phase. No existing AFNI reader, 3dDeconvolve, REML/LSS implementation or
  additional build/test registration edit was needed.

Examples:

```sh
# groups.1D: 0 0 1 1 -1  (four predictors plus one nuisance column)
# penalties.1D: multiple rows, each holding lambda_group0 lambda_group1
3dPicard -input epi+orig -matrix X.xmat.1D -penalty banded \
         -band_groups groups.1D -band_penalties penalties.1D \
         -cv_predict -prewhite_ar1 .3 -prefix stats.banded

3dPicard -input epi+orig -matrix X.xmat.1D -penalty lasso \
         -lambdas '1,.3,.1,.03' -cv_predict -support_refit \
         -prewhite_ar1 .3 -prefix stats.support
```

Purpose: estimate different shrinkage for explicit predictor families, account
for a supplied temporal correlation model without crossing validation
boundaries, and evaluate sparse selection followed by unpenalized support
refitting. Prompt before phase 7 (noise-pool PCs, PC-count selection and the
single-trial orchestrator). Maker milestones remain pending.

## Phase 7 implementation record — 2026-11-06

User authorization: "go for phase 7", followed by explicit approval of the
staged workflow: supplied HRF selection, run-specific noise PCs, common PC-count
selection within an explicit mask, and final penalty tuning. All upstream
choices make the reported scores conditional selection scores; this is not an
independent outer-test assessment. Maker milestones remain separate.

Implemented:

- Added `-single_trial` to `3dPicard`. A supplied HRF candidate collection is
  selected by no-PC unregularized full-data residual MSE and then frozen. A
  single supplied matrix is also supported. Trial and nuisance columns must
  be run-local; explicit trial identities are required even for fixed fits.
- Noise voxels come from `-noise_pool MASK`, or from the explicit combination
  `-noise_model COARSE.xmat.1D -noise_search MASK -noise_r2 R`. The latter uses
  task-model improvement relative to the coarse model's nuisance-only SSE;
  undefined denominators are excluded. The coarse model must share sampling,
  censoring, run boundaries and TR, with no more columns than the trial design.
- New `thd_picard_single.[ch]` owns bounded-batch run-specific temporal PCA,
  augmented design construction and PC-count scoring. It reuses `thd_linproj`
  and the existing SVD wrapper; no MVPA PCA extraction or new eigensolver.
  Each run removes a constant plus supplied nuisance span, normalizes residual
  pool voxels to unit L2 norm and accumulates a temporal Gram matrix in fixed
  voxel order. PCs have unit sample SD and deterministic largest-entry sign.
  Gram rank uses `max(rank_tol^2,64*DBL_EPSILON*run_rows)` times the largest
  singular value; insufficient requested rank is an error, never clipping.
- `-noise_pcs K` fixes the count. Alternatively `-pc_counts 0,1,... -pc_mask MASK`
  selects a common count. The selection mask must lie inside analysis and
  outside the pool. Ordered cross-run condition-pair errors compare each
  candidate's OLS beta to a fixed no-PC OLS reference, divide by reference
  sample beta variance and pair count, then take the selection-mask median.
  Zero reference variance is excluded consistently at all counts; all undefined
  is an error. Choose the first count reaching best improvement / `-pc_stop`
  (default 1.05); no improvement selects zero. Fixed count and fixed penalty
  work without repeats. PC selection itself requires cross-run repeats.
- PCs are appended as run-local unpenalized columns, preserving original trial
  and penalty-group identities. Final fixed/repeated-condition fitting retains
  fractional ridge, lasso/enet, banded ridge, fixed AR(1), and sparse support
  refit options. The candidate API accepts locked HRF indices and original
  scores, invalidates failed results, and avoids exposing a stale score grid.
- The default remains one native AFNI HEAD/BRIK bucket. Four added maps record
  PC count, pool membership, supplied PC-selection membership and coarse-model
  R2. Membership maps cover the whole grid; PC count is -1 outside analysis;
  undefined/uncomputed R2 is -1. `-save_stages` adds labeled no-PC OLS (Stage B)
  and selected-PC OLS (Stage C) beta sub-bricks to that bucket; ordinary betas
  are the final Stage D. The optional fitted/residual buckets retain their
  existing meaning and original-time row maps.
- Headers record the retained-row/own-run PC array, run ranks, tested count/loss
  trend, masks, chosen count, score policy and augmented nuisance identities.
  Fixed counts have trend score -1 (not computed). The header also records the
  number of nonzero-variance selection voxels actually contributing to medians.
  Standard dimensions/roles/labels describe the augmented final design.
- All inputs and masks use AFNI dataset loading. New preparation and fitting
  failures preserve the output transaction. Memory preflight includes temporal
  covariance/SVD, augmented designs and tuning plans, workers, scalar selection
  tables, provenance and requested maps. It remains an estimate, not a hard RSS
  limit. `-validate_only` checks supplied designs/options; mask contents, the
  coarse model and learned PC rank require the actual data pass.
- Expanded help and source comments describe indexing, ownership, fixed HRF
  choices, PC normalization, scoring and stage lifetimes. Picard source/help
  stays within 80 columns. No `goto` statements were added.

Validation:

- CMake and classic Make builds passed. All **35 selected regression tests**
  passed, including existing MVPA, shared regression, xmat and Picard coverage.
  The existing full Picard integration suite also passed with the classic
  binary. Final focused Picard checks cover the latest source and fixtures.
- New numerical unit coverage checks analytic rank-one PCs, rank rejection,
  batch invariance, augmented identities, undefined variance, stopping policy
  and failed-state handling. Candidate tests exercise locked choices and
  original score retention, including invalid lock rejection.
- New `test_picard_single.py` independently checks HRF choices, PCA, the whole
  PC-count curve and chosen count, final repeat penalty tuning and coefficients,
  Stage B/C maps, native labels/provenance, batching/thread consistency, fixed
  counts, coarse-model pool R2, sparse/support/AR1 integration, censored-row PCA
  plus independently calculated AR1 fits, and banded-ridge group preservation.
  Failure cases include empty/overlapping masks, insufficient noise rank,
  malformed count policy, budget rejection, all-undefined selection, empty
  model-derived pools and rollback after a companion write failure.
- ASan/UBSan passed for numerical noise/selection and candidate units and the
  staged CLI integration, including partial failures. Strict C99 warning checks
  passed for the Picard modules. Existing external/AFNI warnings were not edited.
  Session evidence is under `/tmp/picard_phase7_reference`.

Pinned GLMsingle comparison:

- `make_picard_single_fixture.py` verifies SHA-256 digests for the unmodified
  sources at commit `082851f9db87ff30c49370f27365ac1112f18ff7`, then executes
  the actual [PCA block](https://github.com/cvnlab/GLMsingle/blob/082851f9db87ff30c49370f27365ac1112f18ff7/glmsingle/glmsingle.py),
  `calcbadness` and `select_noise_regressors`. Normal tests consume only the
  saved `picard_single_reference.json`; no network or GLMsingle installation.
- The controlled end-to-end adapter case uses one supplied design, explicit
  pool/selection masks, run intercepts, a single session and final fraction 1.
  NumPy supplies OLS; pinned upstream code supplies PCA, pair scores and PC
  stopping. It checks PC values, the score curve (upstream SSE / pair count),
  selected count, Stage B/final coefficients and fitted data. Float32 upstream
  projector/output differences have explicit comparison tolerances. A separate
  independent fixture exercises multiple supplied HRFs and final penalty CV.
- This is **not** a call to the complete `GLM_single.fit` pipeline or a claim
  of full GLMsingle equivalence. HRF/timing generation remains external; the
  pool/search/PC masks and coarse task model are explicit rather than its
  automatic intensity/R2/percentile policies. There are no implicit polynomial
  drifts, session-specific normalization, percent-signal conversion or beta
  autoscaling. Picard uses double PCA, a numerical rank gate and deterministic
  sign; degenerate PC eigenspaces need not have identical individual vectors.
  Zero reference variance is excluded rather than following upstream zerodiv
  behavior. Optional whitening/support/sparse/banded fits are Picard extensions.
  Final fractional tuning remains Picard's fold-local reference-norm policy.
  As in the candidate collection, HRF designs must retain equal projected
  ranks after augmentation; incompatible rank changes fail explicitly.

Purpose: one reproducible AFNI command now follows the chosen HRF through
noise removal and penalty selection, with intermediate betas available in the
same bucket and enough provenance to reconstruct the augmented design.
This completes Picard stage 7. Prompt before beginning a new milestone;
**1dXmatMaker M0/M1 are still pending**, not completed by this work.

## GLMsingle gap analysis and proposed noise-model (ARMA/REML) study — 2026-11-07

Status: **assessment and proposal only; nothing here is implemented or
authorized.** Written after comparing Picard (stages 1-7) and Maker (M0-M8)
against GLMsingle. GLMsingle facts below come from its eLife paper and wiki
and from memory of its options; they were not checked against its source.

### Remaining differences from GLMsingle

Picard covers HRF-library selection, run-specific noise PCs with PC-count
selection, and fractional ridge tuned by repeated-condition CV. Maker supplies
the pinned 20-HRF library. Not implemented (all are policy/convenience
differences, not deep methodology):

- Automatic noise-pool and PC-cutoff policies (Picard needs an explicit mask or
  a coarse model plus `-noise_r2`).
- Post-ridge per-voxel scale/offset ("autoscale"); deliberately omitted.
- Percent-BOLD conversion and session normalization; multi-session handling.
- Automatic drift polynomials (supplied in the xmat instead).
- Per-run R2 and per-run HRF-index maps.
- Fraction selection differs: GLMsingle scores a full-data fraction path;
  Picard uses fold-local fractions with raw coefficient-unit MSE. State this
  in any comparison.

Picard additionally offers lasso/elastic net, banded ridge, prediction CV with
one-SE and global selection, support refits, censoring and native provenance.
Only fixed, user-supplied AR(1) whitening exists; there is no noise estimation.

### Pilot simulation (effect-size and bias probe, not validation)

`sim_whiten_singletrial.py` (this directory, numpy only; run with argument
`a`, `b` or `c`). Design: 4 runs x 240 TRs, TR 1.5 s, jittered ISI 3-6 s, 301
single-trial regressors, per-run order-3 Legendre drift unpenalized, true trial
betas N(1, 0.5^2), noise SD 1, 80 realizations. Ridge lambda was the single best
value per method (averaged over realizations), so tuning is not a confound.
Whitening was per run.

| True noise | RMSE gain, oracle whitening | AR(1) from full single-trial residuals | AR(1) from ONOFF+drift residuals |
|---|---|---|---|
| AR(1) a=0.3 | 0.5% | 0.4% (a_hat 0.29) | 0.5% (a_hat 0.40) |
| AR(1) a=0.6 | 1.8% | 1.5% (a_hat 0.60) | 1.9% (a_hat 0.65) |
| ARMA(1,1) a=0.8, b=-0.5 | 2.9% | 0.2% (a_hat 0.11) | 2.2% (a_hat 0.46) |

Reading: whitening changes point estimates only slightly in this low-SNR toy.
The estimation stage is the hard part. A near-saturated trial design absorbs
slow noise, so REML on the full design can miss it (ARMA case), while an ONOFF
model overestimates when signal varies across trials (AR(1) cases). A third-party
report (fmrilss PR 21; relationship to 3dLSS unverified) describes the same
downward bias on a denser design. Limits: one SNR, one design, AR(1) whitening
only, oracle lambda, RMSE rather than inference or beta covariance. Results do
not transfer to real data.

### Recommendation

Pursue, framed around noise-model estimation and uncertainty rather than better
point estimates:

1. **Hyperparameters without repeats.** Ridge is a random-effects model; joint
   REML of the ridge variance and ARMA parameters could tune shrinkage on
   designs with no repeated conditions, which repeated-condition CV cannot.
   An unpenalized ONOFF column (already possible via roles) shrinks toward the
   mean response instead of zero.
2. **Trial-beta covariance** from whitened ridge, for downstream RSA/decoding
   (cf. ITEM, which handles trial covariance analytically).
3. **When noise estimation fails** under near-saturated designs, with a
   recommended two-stage estimator (noise model from a compact model, then
   single-trial fit).

### Preconditions and staging (proposed)

- Systematic prior-art search first (fmrilss, ITEM, fmriAR, Soch et al.); the
  searches behind this note were quick and not systematic.
- Voxelwise whitening breaks one-SVD-per-design reuse. Rough, unmeasured
  estimate: a 16-point AR(1) grid x 20 HRFs x 12 folds is about 4,000 SVDs
  (feasible); a full ARMA grid across HRFs, folds and PC counts is not. Freeze
  whitening per voxel after a pilot stage, use grid cells as 3dREMLfit does, and
  consider pooling across parcels.
- Cheapest first step: import voxelwise (a,b) from `3dREMLfit -Rvar` (a design
  mirroring its `-ABfile`), validated against fixed-AR(1) results.
- Estimate whitening within training folds, or disclose the leakage; never
  whiten across CV boundaries implicitly.
- Evaluate against GLMsingle on real test-retest data, not simulation alone.
  Scores remain selection scores, not outer-test performance.

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

The single-trial orchestrator is now `3dPicard -single_trial`: it consumes
supplied HRFs, estimates noise-pool PCs, selects PC counts and emits staged
outputs. Automatic HRF-library generation belongs to Maker. Fractional ridge
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

At the 2026-10-30 review, the completed numerical implementations were still in
`thd_mvpa.c`, with public interfaces in `thd_mvpa.h`. Phase 1 has now moved them
to `thd_regression.[ch]`; later phases added Picard projection and workflows.
The current consolidated layout places SVD/fractional fitting in
`thd_regression.c` and projection plus all Picard workflows in `thd_picard.c`.
The graph below describes the intended complete architecture. Use one regression
library for the estimator family, not a library per penalty.

```text
3dMVPA.c -> thd_mvpa -> thd_regression
                     -> thd_crossval / existing pattern and inference modules

3dPicard.c -> thd_picard -> thd_regression
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

| File/module | Owns | Does not own | Current boundary / remaining work |
|---|---|---|---|
| `thd_regression.[ch]` | Ridge and continuous-estimator fits, thin SVD/fractional paths, coefficients, fitted numeric state and diagnostics. | Dataset I/O, labels, censor/run interpretation, fold construction, candidate selection, inference. | Consolidated shared implementation; retain estimator units and existing MVPA behavior. |
| `thd_picard.[ch]` | Nuisance projection/recovery and all numeric Picard workflows: fixed/repeated/prediction fits, HRF selection, noise PCs, PC scoring and final refits. | AFNI dataset handles, command parsing, file writing, HRF generation, RSA statistics. | One implementation/header pair above shared solvers and folds; five numbered sections. |
| `thd_crossval.[ch]` | Index-based fold plans, nesting, validation, generic search traversal. | HRF/fraction interpretation, nuisance fits, repeated-condition beta score definitions. | Reuse; add metadata-only fold construction and a separate path/batch traversal when needed. Preserve current scalar-search API. |
| `thd_xmat.[ch]` | Matrix parsing/writing, labels/attributes, retained/full row mappings, run/column transformations, metadata adapters. | Penalized fitting, response-specific preprocessing, HRF candidate selection. | Reuse existing API; extend only general metadata/mapping needs. |
| `thd_mvpa.[ch]` | MVPA preprocessing, classifiers, selection/PCA workflows, CV policies, permutations, audits, activation patterns and frozen-model formats. | Ownership of the extracted regression implementations. | Include `thd_regression.h`, retain workflow APIs and serialized behavior, replace private solver-field access. |
| `3dPicard.c` | CLI, xmat/manifest ingestion, domain checks, mask and response batching, worker scheduling, AFNI outputs, provenance. | Numerical algorithms and candidate scoring logic. | New executable adapter; start with existing 3dDeconvolve matrices. |
| `thd_patterns`, `thd_sampleset`, `thd_datatable` | Existing spatial extraction and table/sample identities. | A mandatory ROI/searchlight representation for voxelwise fitting. | Reuse only where needed; begin with a local bounded response reader in the executable. |
| `thd_rsa`, `thd_simmatrix`, `thd_mapinfer`, `thd_permute` | Their established representational/inference responsibilities. | Picard's initial solver or model-selection engine. | No Picard-driven changes or new inference dependencies. |

### Stage 1: exact estimator extraction and compatibility surface (completed)

The following is the original extraction specification; the implementation
record above documents its completed API and validation.

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
and original-input-space intercepts. These accessors now fill the documented
`THD_regression_state` view. New public headers contain no AFNI datasets,
THD_sampleset, THD_foldplan or permutation types. The private
`thd_regression_internal.h` declares the SPD helper for LDA; it must not
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

The projection section of `thd_picard.h` exposes an opaque `THD_linproj` with prepare, project,
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

The executable loads full-timeline volume data and masks through AFNI's
`THD_open_dataset` / `DSET_load` interfaces, including NIfTI. It inherits the
loader's datatype conversion, scaling and nonfinite-value handling without
forcing `AFNI_FLOATSCAN`. Bounded retained-row Y batches are gathered from the
loaded AFNI bricks, promoting samples to double and applying their brick factors.
This cannot restore precision or nonfinite values already changed by the loader.
Any nonfinite values still present in selected retained samples or the loaded
mask are rejected by Picard. The xmat
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
`thd_regression`, `thd_crossval` and the required xmat/AFNI adapters; nuisance
projection now lives inside the consolidated `thd_picard` implementation.
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
| 1 (complete) | Created `thd_regression.[ch]`, moved ridge plus all four continuous estimators, replaced MVPA private-field access and updated build lists. | Passed 2026-11-06: ridge/EN/banded/PLS/SVR, LDA helper, CV, patterns, coefficients and saved-model coverage. No new numerical algorithm in this change. |
| 2 (complete) | Added SVD/fractional paths, `thd_linproj`, fixed-design `thd_picard`, and the `3dPicard.c` adapter with native AFNI buckets. | Passed 2026-11-06: existing Deconvolve matrices, fracridge/direct-model oracles, rank/endpoints, nuisance recovery, labels/maps, batch/thread invariance and memory-budget checks. |
| 3a (complete) | Added numeric time-series folds/path traversal in `thd_crossval` and shared-feature prediction workflow in `thd_picard_predict`. | Passed 2026-11-06: fold-local fitting, held-out perturbation, explicit raw/conditional nuisance policies, defined baseline, one-SE and deterministic voxel/global selection, undefined R2, full refits and original-time gap/fold metadata. |
| 3b (complete) | Added separate repeated-condition beta-tuning mode to `thd_picard` using explicit event identities and shared generic traversal where applicable. | Passed 2026-11-06: hand-computable fixture, independent fold/refit oracles, insufficient-repeat rejection, pinned GLMsingle scoring comparison and native bucket integration. |
| 4 (complete) | Shared EN kernel with prepared-double fits and warm paths; lasso/enet options, scaling, diagnostics and both tuning modes. | Passed 2026-11-06: scikit-learn/direct-solve objective/scaling/fold oracles, full KKT, warm/cold agreement, all 29 selected regression tests and sanitizer checks. |
| 5 (complete) | Numeric candidate collections, versioned identity manifests, explicit HRF and prediction selection policies, native ID/score maps. | Passed 2026-11-06: per-voxel IDs/beta identities, independent numerical selection/refit oracles, all 30 selected tests and sanitizer checks. Selection-only scores; no outer-test claim; ordinary matrices remain independent of Maker. |
| 6 (complete) | Shared prepared banded kernel, fixed AR(1) whitening, support-only sparse refits, bounded factor reuse and measurements. | Passed 2026-11-06: legacy parity, independent fixed/fold/refit oracles, run/censor/fold boundaries, group remapping, native diagnostics, all 30 selected tests, sanitizers and runtime/RSS/thread checks. |
| 7 (complete) | Staged single-trial mode with locked supplied HRFs, run-specific noise PCs, explicit PC-count selection and optional stage beta sub-bricks. | Passed 2026-11-06: pinned GLMsingle component/controlled end-to-end fixture, independent multi-HRF/final-tuning oracles, native provenance, all 35 selected tests, both builds and sanitizers. Full-package/automatic-policy differences documented above. |

Stages 3a and 3b are separate deliverables, not interchangeable implementations.
For a single-trial-first release, prioritize 3b and stage 5 over sparse penalties.
MVPA's numerical prerequisites and Picard stages 1, 2, 3a, 3b, 4, 5, 6 and 7 are complete.
Prompt before beginning a new milestone. Maker M0/M1 remain separate planned
milestones. Do not copy EN or
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
  See the 2026-11-07 gap analysis and ARMA/REML proposal above.
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

**2026-11-06 decision:** implement Maker in Rust, with development confined to
`src/pmolfese/1dXmatMaker/`. The dedicated
[1dXmatMaker plan](1dXmatMaker_PLAN.md) is now the working implementation
roadmap and preserves M0-M8. Extensive explanatory comments for non-Rust
developers are required at every milestone. Maker M0-M2 are complete as of
2026-11-06; see the [M0 report](1dXmatMaker/docs/validation/M0.md),
[M1 report](1dXmatMaker/docs/validation/M1.md), and
[M2 report](1dXmatMaker/docs/validation/M2.md). The standalone Rust partial
prototype builds baseline/nuisance and selected ordinary/IM event matrices.
M3-M8 remain planned. The 1dXmatTool sparse-ID follow-up is resolved; see the
[fix report](1dXmatMaker/docs/validation/SPARSE_STIM.md).

### Scope and compatibility boundary

The name is settled: **1dXmatMaker**. It is a standalone matrix producer
feeding **3dLSS, 3dREMLfit and 3dPicard**. The implementation language is
**Rust**. `1dXmatTool` retains inspection, validation, selection and splitting
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

Independent Rust basis and design modules may
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

The file interface and untouched 3dDeconvolve support the chosen standalone
Rust architecture. Reproduce basis/expression-parser behavior and xmat output
through explicit compatibility tests. Build/distribution support and extensive
comments for non-Rust maintainers are acceptance requirements; no speedup is
assumed from language alone. A narrow C compatibility component would require
a separate documented design decision, not an implicit AFNI dependency.

### Milestones and acceptance gates

M0-M8 are stable IDs for Maker work, separate from Picard's stage numbers.
M0-M2 are complete; M3-M8 are planned. Each milestone should end with reviewable outputs, documented
coverage and a validation report; partial coverage must not be advertised as
complete 3dDeconvolve compatibility.

| ID | Deliverable | Dependencies | Acceptance gate |
|---|---|---|---|
| M0 | Compatibility inventory, pinned reference corpus, Rust build policy, and input/output specification. | Picard contracts; fixtures may be collected before its initial release. | Every current basis family, variant/alias, timing mode and relevant option has a named test or an explicit unresolved item. Record AFNI reference revision, environment settings, precision and numerical tolerances. |
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

- Rust is chosen. Follow the dedicated Maker plan for independent basis/design
  modules, direct xmat serialization, build policy, and source organization.
  Complete basis and EXPR compatibility remain required; use actual AFNI
  readers and downstream fits to validate the file boundary.
- Extensive comments must explain algorithms, units, indexing, ownership,
  borrowing, and error handling to non-Rust developers. Review readability
  alongside scientific correctness at each milestone.
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
