# 3dMVPA Design and Implementation Plan

## Active Roadmap

Keep this ordered list at the top of the plan and update it whenever priorities
or completion state change.

1. **Prediction Metrics (Complete)** -- classification and continuous-prediction
   scoring now live in `thd_mvpa.[ch]`, with known-answer tests.
2. **Baseline Classifiers (Complete)** -- nearest-centroid fit/predict now
   supports Euclidean, correlation, and cosine decisions in `thd_mvpa.[ch]`.
3. **Fold Runner and Workspace (Complete)** -- fold plans, preprocessing,
   centroid fitting, held-out prediction, and scoring now share reusable buffers.
4. **Selected-Sample Pattern Extraction (Complete)** -- `THD_roi_samples()` now
   maps all or selected `THD_sampleset` rows into sample-by-voxel arrays.
5. **First 3dMVPA ROI Classification Slice (Complete)** -- the native program
   now runs direct or X-matrix samples through validated ROI classification.
6. **Volumetric Searchlights (Complete)** -- moving in-mask neighborhoods now
   use the same extraction, fold-local modeling, and map outputs as ROI analyses.
7. **Prediction Audit Tables (Complete)** -- opt-in row-level predictions and
   per-fold scores retain spatial-unit IDs, sample metadata, and class labels.
8. **Surface Searchlights (Complete)** -- geodesic neighborhoods now share the
   volume fold runner, maps, audit tables, and strict sparse-node validation.
9. **Blocked Permutation Inference (Complete)** -- exchangeability-aware label
   nulls now produce synchronized PVAL/Z/Q/FWE maps for ROI and searchlights.
10. **Shared Surface Domain Helpers (Complete)** -- `thd_patterns.[ch]` now owns
   exact dense/sparse node-domain comparison and sparse output-domain copying
   used by both `3dMVPA` and `3dRSA`.
11. **`thd_rsa.[ch]` RDM Construction Extraction (Complete)** -- the coherent
   RSA science library now owns condition centering, subject RDM construction,
   triangle conversion, mean RDMs, and reusable second-order RDM workspaces.
12. **Ridge Decoding and Encoding (Complete)** -- primal/dual multi-output
   ridge, training-only preprocessing/intercepts, continuous prediction audits,
   and ROI/volume/surface decoding plus voxel/pattern encoding are implemented.
13. **Nested Grouped Lambda Selection (Complete)** -- grouped inner folds are
   restricted to each outer training set; preprocessing is refitted inside
   every split, with candidate losses and selected penalties in audit tables.
14. **Continuous Permutation Inference (Complete)** -- explicit blocked target
   or design-row nulls refit preprocessing, tuning, and models, with empirical
   P/Z/Q/FWE outputs for decoding and voxel/pattern encoding.
15. **Ridge Classification (Complete)** -- 0/1 class-coded ridge responses,
   deterministic multiclass argmax decisions, grouped lambda selection,
   prediction/fold/tuning audits, and fully refitted label permutations share
   the existing classification runner and solver.
16. **Diagonal LDA (Complete)** -- fold-local class means and pooled diagonal
   variance, explicit class priors, and relative variance floors share the
   classification audits and fully refitted blocked nulls.
17. **Shrinkage LDA (Complete)** -- configurable spherical covariance shrinkage
   uses primal/sample-space solves for correlated high-dimensional patterns,
   with training-only fit diagnostics and known-answer stability tests.
18. **Descriptive Coefficient Maps (Complete)** -- opt-in full-data ridge
   refits export coefficients, intercepts, preprocessing, feature/response IDs,
   and volume/sparse-surface maps separately from held-out scores.
19. **Native Linear SVM (Complete)** -- native squared-hinge one-versus-rest
   classification, grouped training-only C selection, verified convergence,
   prediction/fold/tuning audits, and fully refitted blocked label inference.


20. **Classification and Subject-Summary Audits (Proposed; Next, awaiting user direction)**
   -- export confusion matrices, per-class recall, and explicit subject/stratum
   summaries without conflating pooled prediction with population inference.
21. **Explicit Train/Test Partition Options (Proposed)** -- expose the existing
   fixed-partition library machinery for prespecified independent test sets.
22. **Fold-Local Nuisance Adjustment (Proposed)** -- define and implement
   training-only covariate adjustment, with an explicit scientific target and
   compatible inference contract.
23. **Fold-Local Feature Selection (Proposed)** -- select informative features
   inside training splits, audit their original IDs, and tune feature counts
   within grouped inner CV.
24. **Fold-Local PCA (Proposed)** -- fit and apply training-only projections,
   with nested component-count selection and fully refitted permutation paths;
   expose principal-component regression by composing PCA with regression.
25. **Haufe/Activation-Pattern Outputs (Proposed)** -- add descriptive signal
   patterns with explicit scaling and interpretation contracts alongside
   existing decoder coefficient outputs.
26. **Paired Model Comparisons (Proposed)** -- compare models on identical
   held-out observations; specify and validate a separate null before adding
   inferential model contrasts.
27. **Saved-Model Application (Proposed)** -- persist fitted preprocessing,
   feature definitions, model parameters and labels for validated application
   to genuinely new datasets.


28. **Elastic-Net Regression (Proposed)** -- add sparse, regularized continuous
   decoding/encoding with nested penalty-strength and mixing-ratio selection.
29. **Regularized Logistic Regression (Proposed)** -- add L2-regularized binary
   and multinomial classification with probability-aware scoring and explicit
   calibration assessment; elastic-net classification penalties are a later
   extension.
30. **Partial Least Squares Regression (Proposed)** -- add supervised latent
   components for single/multiple-target decoding and multi-response encoding,
   with training-only fitting and nested component-count selection.
31. **Linear Support-Vector Regression (Proposed)** -- add epsilon-insensitive
   continuous decoding, with encoding support assessed for computational cost,
   nested C/epsilon selection and explicit solver convergence.
32. **Robust Huber Regression (Proposed; Later Priority)** -- add an optional
   robust loss for continuous decode/encode applications where large response
   residuals dominate squared-error fits, with nested regularization/loss tuning.
33. **Generalized Inner-CV Tuning (Proposed)** -- share estimator-independent
   grouped candidate evaluation, explicit selection objectives and audit schemas,
   including ACC/BACC classifier tuning rather than ridge-only squared loss.
34. **Response-Specific Ridge Penalties (Proposed)** -- optionally select one
   lambda per response within each outer training set, retaining shared-lambda
   behavior as the default and making the extra cost explicit.
35. **Banded Ridge Encoding (Proposed)** -- assign separate penalties to named
   predictor feature groups, with grouped nested selection and feature-space
   contribution summaries under an explicit interpretation contract.
36. **Reduced-Rank Regression (Proposed; Later Priority)** -- add a distinct
   low-rank multi-response regression option with nested rank/penalty selection
   and an explicit shared-response domain.
37. **Kernel Ridge and Multiple-Kernel Ridge (Proposed; Later Priority)** --
   stage explicit-kernel regression after linear baselines; validate training-only
   kernel preprocessing, nested kernel/penalty selection and memory limits.
38. **Cross-Condition Generalization (Proposed)** -- train on one specified
   task, condition or session domain and test another, with explicit sample,
   class/response and feature-domain matching.
39. **Time-Series Encoding (Proposed)** -- add run-aware lagged predictors,
   buffered temporal partitions and separately validated temporal nulls, without
   treating serially dependent rows as freely exchangeable samples.
40. **Elastic-Net Logistic Classification (Proposed; Later Priority)** --
   extend item 29's L2 baseline with sparse L1/L2 classification penalties,
   nested strength/mixing selection and audited convergence.

Items 20–40 are proposed scope, not implemented capabilities. Adding them to
this roadmap does not authorize implementation; obtain user direction before
starting item 20 and before advancing to each subsequent item.
Existing item numbers and item 20's next-step priority are preserved. Numbers
identify proposals, not permission to skip dependencies: item 33 should provide
shared tuning infrastructure for later estimators, and item 34 should precede
the multi-penalty encoding work in item 35.

Development cadence: complete one authorized roadmap item, explain why each
part matters and its success metric, report what fMRI users can now do, then
ask whether to proceed to the next item or expand the plan/scope. Completing
an item does not authorize starting the next one.

## Implementation Status

This checked-in copy is the living implementation record for the `somemvpa`
branch. Update it with each completed slice so that design decisions and code do
not drift apart.

Status as of 2026-10-15:

| Slice | Status | Implementation |
|---|---|---|
| X-matrix input adapter | Complete | `thd_xmat.[ch]`, quoting-aware `thd_datatable`, unit tests |
| Canonical sample organization | Complete | `thd_sampleset.[ch]`, direct/Xmat-expanded table tests |
| General fold plans | Complete | `thd_crossval.[ch]`, grouped/within-subject/fixed fold tests |
| Fold-local centering and scaling | Complete | `thd_mvpa.[ch]`, train-only fit/apply tests |
| Prediction metrics | Complete | `thd_mvpa.[ch]`, allocation-free classification/continuous scorers and known-answer tests |
| Baseline classifiers | Complete | `thd_mvpa.[ch]`, Euclidean/correlation/cosine centroid models and known-answer tests |
| Fold runner and workspace | Complete | `THD_mvpa_workspace`, allocation-free repeated centroid CV runs, row/fold diagnostics |
| Selected-sample pattern extraction | Complete | `THD_roi_samples()` preserves requested row/voxel order across datasets and scaled brick types, with bounded diagnostics |
| Native program and ROI classification | Complete | `3dMVPA.c`, CMake/Make targets, direct/Xmat input, grouped CV, fold-local centroid analysis, ROI table and AFNI map |
| Volumetric searchlights | Complete | `-searchlight` with SPHERE/RECT/RHDD/TOHD neighborhoods, one workspace per OpenMP worker, center-painted SCORE/NFEAT/NFOLD maps |
| Prediction audit tables | Complete | `-save_predictions` and `-save_fold_scores` write deterministic ROI/searchlight-keyed TSV files with metadata, decisions, folds, counts, ACC, and BACC |
| Surface searchlights | Complete | `-surf` exposes masked or whole-represented-domain geodesic neighborhoods, preserves sparse mesh node IDs, and rejects mismatched/malformed domains |
| Blocked permutation inference | Complete | `-nperm`, `-seed`, and explicit `-permute_within` blocks drive identity-inclusive label nulls, empirical above-chance/P/Z/Q maps, and synchronized max-statistic FWE |
| Shared surface domains | Complete | `THD_surface_domain_match()` compares dense and sparse node rows exactly; `THD_copy_surface_node_domain()` preserves sparse output node IDs without replacing result labels/statistics; both programs use the shared helpers |
| RSA RDM construction extraction | Complete | `thd_rsa.[ch]` owns dataset-independent condition centering/RDMs, exact triangle conversion, descriptive mean RDMs, and reusable second-order workspaces; direct known-answer, full oracle, and 22-subject real-data replay tests preserve numerical behavior exactly |
| Ridge decoding and encoding | Complete | `THD_ridge_workspace` and `THD_mvpa_ridge_workspace` provide primal/dual multi-output ridge, unpenalized intercepts, fold-local preprocessing, held-out training-mean baselines, and per-response/fold metrics; CLI decode/encode works for ROI, volume, and sparse-surface domains |
| Nested grouped lambda selection | Complete | `THD_nested_foldplan` in `thd_crossval.[ch]` constructs compact inner grouped plans; `thd_mvpa.[ch]` reuses tuning buffers, normalizes inner errors per response, refits selected models, and reports all candidate losses plus selected penalties |
| Continuous permutation inference | Complete | Joint blocked target/design-row nulls rerun preprocessing, nested tuning, and ridge fits; identity-inclusive empirical P/Z, spatial BH Q and synchronized max-statistic FWE work for ROI, volume, and sparse surfaces, with undefined-score safeguards and inference audits |
| Ridge classification | Complete | The shared classification runner uses 0/1 class indicators, multi-output ridge with unpenalized intercepts, and deterministic argmax decisions; fixed/grouped nested penalties, class coverage validation, observed audits and fully refitted blocked label inference work for ROI, volume and sparse surfaces |
| Diagonal and shrinkage LDA | Complete | `thd_mvpa.[ch]` owns allocation-free repeated LDA fits with pooled training covariance, equal/empirical priors, explicit variance floors, and primal/Woodbury shrinkage solves; shared ROI/volume/surface classification audits and synchronized label inference preserve original sample/node IDs |
| Descriptive coefficient maps | Complete | `-save_coefficients CPREFIX` adds an independent full-data ridge refit for decode, encode, and ridge classification; coefficient/intercept and preprocessing TSVs identify every predictor/response and map brick, with volume/sparse-surface output and explicit fixed refit penalties for tuned CV |
| Native linear SVM | Complete | `THD_svm_model` and `-method svm` provide deterministic squared-hinge classification with a penalized intercept, fixed/grouped nested C selection, explicit convergence diagnostics, shared ROI/volume/surface audits and fully refitted blocked label nulls |

## Purpose

`3dMVPA` will be a native AFNI C program for multivariate predictive analysis.
It will not be a wrapper around Python, R, Matlab, SVM-light, or a collection of
external command-line tools. It will own sample organization, grouped
cross-validation, leakage-safe preprocessing, prediction, permutation
inference, spatial map construction, and audit outputs.

The program will expose three related analysis directions:

```text
3dMVPA -mode classify ...
3dMVPA -mode decode   ...
3dMVPA -mode encode   ...
```

- `classify`: neural patterns predict a categorical target.
- `decode`: neural patterns predict one or more continuous targets.
- `encode`: an external design predicts neural responses.

Classification is technically categorical decoding, but it deserves a separate
mode because categorical targets require different estimators, scores, chance
definitions, diagnostics, and output tables.

`3dMVPA` should be more than a modern spelling of `3dsvm`. The latter trains or
applies an SVM to a 3D+t dataset using separately supplied labels and a legacy
model representation. `3dMVPA` should make samples, targets, subjects, runs,
folds, exchangeability blocks, and covariates explicit in one design and then
own the entire validation and inference procedure.

## Relationship to 3dRSA

`3dRSA` and `3dMVPA` should remain separate user-facing programs.

- `3dRSA` asks whether neural representational geometry agrees with one or more
  model geometries.
- `3dMVPA` asks whether held-out outcomes or neural responses can be predicted.

They nevertheless share a large substrate:

```text
datasets and sample metadata
        |
        v
ROI / searchlight feature extraction
        |
        v
sample x voxel matrix
        |
        +-----------------------+
        |                       |
        v                       v
   thd_rsa backend         thd_mvpa backend
   RDM construction        fit / predict / score
   RDM inference           CV / prediction inference
        |                       |
        +-----------+-----------+
                    v
          AFNI maps and text audits
```

The shared code should live in `thd_` libraries. Neither program should contain
private copies of sample indexing, spatial extraction, folds, numerical
estimators, permutation bookkeeping, or map correction.

The MVPA source boundary is intentionally compact:

- `thd_xmat.[ch]` adapts X-matrices and is also shared with `3dSRM`;
- `thd_sampleset.[ch]` owns canonical sample and dataset organization;
- `thd_crossval.[ch]` constructs and validates reusable fold plans;
- `thd_mvpa.[ch]` owns preprocessing, metrics, workspaces, and estimators;
- `3dMVPA.c` will own command-line policy and output orchestration.

Preprocessing and prediction metrics do not get separate `thd_` source modules.
Existing spatial and inference modules such as `thd_patterns`, `thd_permute`,
and `thd_mapinfer` remain shared rather than being duplicated for MVPA.

## What Already Exists

The AFNI tree already contains most of the neuroimaging and inference
infrastructure needed by `3dMVPA`.

### `thd_datatable.[ch]`

Already provides:

- file and inline table parsing;
- numeric and categorical columns;
- `Subj` and `InputFile` conventions;
- case-insensitive column lookup;
- row selection;
- Cartesian indexing over dimensions such as subject, run, and condition.

This remains the public metadata parser. A competing MVPA-only table parser is
unnecessary.

### `thd_xmat.[ch]`

The `somemvpa` branch now provides:

- parsing of the `X.xmat.1D` files written by `3dDeconvolve`;
- stimulus, item, run, censoring, and design-column metadata;
- lookup of the corresponding `_Coef` sub-bricks in a stats bucket;
- expansion of a subject-level table into an ordinary one-row-per-item
  `THD_datatable` through `THD_datatable_expand_xmat_items()`.

This is an optional input adapter for users whose samples are already described
by an AFNI design matrix. It must not become a requirement for `3dMVPA`: users
must also be able to specify every sample and its metadata directly in the
ordinary `-dataTable`.

### `thd_patterns.[ch]`

Already provides:

- atlas ROI construction;
- volumetric and surface searchlights;
- analysis-unit painting;
- ROI mean-series extraction;
- `[brick][voxel]` pattern extraction into caller-owned storage;
- selected canonical sample extraction into caller-owned `[sample][voxel]` storage;
- the RSA-specific run-aware `THD_runset` container.

The spatial abstraction is already correct for classification, decoding, and
encoding. It should remain focused on spatial geometry and extraction rather
than accumulating estimators or RSA inference.

### `thd_permute.[ch]`

Already provides:

- label permutations and sign flips;
- exchangeability blocks;
- exact and Monte Carlo relabeling sets;
- ordinary and stratified subject bootstraps;
- time-shift and phase-randomization families;
- generic and Freedman--Lane permutation drivers;
- empirical p-values, signed z scores, and max-statistic result structures.

For MVPA, a permutation must rerun the complete cross-validation pipeline, but
the permutation definitions and result accumulation should reuse this library.

### `thd_mapinfer.[ch]`

Already provides BH FDR, masked FDR, running max-null accumulation, and generic
memory planning.

### `thd_simmatrix.[ch]`

This remains the low-level similarity/RDM library. Some numerical components
are useful to MVPA—ranking, correlations, covariance shrinkage, and whitening—
but the file should not become a general machine-learning dumping ground.

### Code not suitable as the new core

`thd_lasso` and `thd_fitter` contain useful algorithms, but their global
configuration and persistent state do not provide a clean reentrant API for
threaded searchlights and nested CV. Algorithms may be adapted, but the new
model API must receive every option explicitly and use caller-owned or
workspace-owned state.

The legacy `3dsvm` implementation is useful as prior art, not as the foundation
of `3dMVPA`. The new program should not inherit SVM-light model files, kernel
machinery, real-time paths, or licensing constraints.

## Core Data Model

At every spatial analysis unit, define:

```text
N: neural matrix       [sample][voxel]
A: attribute matrix    [sample][variable]
```

The modes differ only in direction and outcome type:

| Mode | Predictors | Responses | Natural output |
|---|---|---|---|
| classify | `N` | categorical column of `A` | accuracy per ROI/searchlight |
| decode | `N` | numeric column(s) of `A` | predictive score per ROI/searchlight |
| encode | `A` | `N` | voxelwise or local-pattern predictive score |

The internal matrix layout should be row-major `[sample][feature]`, matching
the layout already produced by `THD_roi_pattern()`.

## Public Sample Input

The initial release will support two equal public input paths. Both produce the
same canonical `THD_datatable` with one row per independent sample before a
`THD_sampleset` is constructed. All later validation, extraction, CV, modeling,
inference, and audit output are therefore shared.

### Direct sample table

Users may specify samples and metadata by hand in an ordinary `-dataTable`. For
beta or trial volumes, `InputFile` may include an AFNI sub-brick selector:

```text
Subj  Run  Trial  Condition  Rating  Fold  InputFile
s01   r1   t001   face       0.72    r1    s01_stats+tlrc'[3]'
s01   r1   t002   house      0.31    r1    s01_stats+tlrc'[7]'
s01   r2   t003   face       0.68    r2    s01_stats+tlrc'[11]'
s01   r2   t004   house      0.28    r2    s01_stats+tlrc'[15]'
```

Advantages of keeping this a first-class public format:

- it is transparent and easy to inspect;
- it reuses `THD_datatable` without another nested parser;
- every sample can carry arbitrary metadata;
- trials are preserved rather than prematurely averaged into conditions;
- the same table works for all three modes.

This path does not require an X-matrix. Manually supplied `Run`, `Trial`,
`Condition`, `Fold`, target, design, and exchangeability columns remain under
the user's control.

### X-matrix expansion

As a convenience for AFNI regression output, users may instead provide one row
per subject with an `Xmat` column and the corresponding whole stats bucket in
`InputFile`:

```text
Subj  Group  Xmat                  InputFile
s01   ctl    s01.X.xmat.1D         s01_stats+tlrc
s02   ctl    s02.X.xmat.1D         s02_stats+tlrc
```

`THD_datatable_expand_xmat_items()` expands each selected stimulus column into
one sample row, supplies item and condition labels, supplies a run label when
the regressor belongs to exactly one run, and resolves the corresponding
`_Coef` sub-brick by label. Options such as `-xmat_stims` may restrict which
stimuli are expanded.

The expanded table is then processed exactly like the direct table. A regressor
that spans runs has no valid run label; requesting run-wise CV in that case must
fail clearly rather than inventing an assignment. `-show_table` and
`-show_samples` should display the expanded representation so the implicit
samples remain inspectable.

The user identifies roles on the command line:

```text
-target_column Condition
-target_columns Rating,Arousal
-design_columns Visual,Semantic,Motion
-cv_column Run
-subject_column Subj
-block_column ExchangeBlock
```

A third, compact container form can be added later for large beta-series
datasets:

```text
Subj Run InputFile            SampleTable
s01  r1  s01_r1_betas+tlrc   s01_r1_samples.tsv
s01  r2  s01_r2_betas+tlrc   s01_r2_samples.tsv
```

That form should also expand into the same canonical sample table; it must not
create a second analysis path.

## `thd_sampleset.[ch]`

The canonical observation map now preserves every sample and connects metadata
rows to datasets and bricks.

Conceptually:

```c
typedef struct {
    int nsample;
    int ndset;
    int nsubject;
    int nrun;

    THD_3dim_dataset **dset;
    int *sample_dset;
    int *sample_brick;
    int *subject;
    int *run;
    int *item;
    int *fold;

    char **sample_label;
    char **subject_label;
    char **run_label;

    THD_datatable *metadata;
} THD_sampleset;
```

Responsibilities:

- consume the canonical one-row-per-sample table regardless of whether it was
  written directly or produced by an input adapter such as `thd_xmat`;
- resolve dataset selectors once rather than repeatedly opening the same file;
- map each observation to dataset, brick, subject, run, item, and fold;
- encode categorical levels while retaining original labels;
- validate grids and sample identifiers;
- represent unbalanced runs and missing cells;
- retain original table row numbers for diagnostics and prediction output;
- provide compact views for selected subjects or samples;
- remain independent of classifiers, RDMs, and inference policy.

`THD_runset` should eventually become an RSA-specific view or adapter over
`THD_sampleset`. It currently averages repetitions into condition estimates,
which is useful for crossnobis RSA but is not the correct general sample model.

Repeated-run indexing, run factors, and run contrasts currently private to
`3dRSA` should move into `thd_sampleset`, `thd_datatable`, or a small general
repeated-measures helper where appropriate.

## `thd_patterns.[ch]` Additions

Keep `thd_patterns` responsible for spatial operations:

- validate that all datasets share the analysis domain;
- extract arbitrary selected samples from one or more datasets;
- construct ROI and searchlight feature lists;
- filter or exclude overlapping analysis units;
- paint one scalar per spatial unit.

A general selected-sample extraction primitive is now implemented as:

```c
int THD_roi_samples(
    const THD_sampleset *ss,
    const int *sample_index, int nsample,
    const intvec *vox,
    float *out,
    char *err, size_t errlen
);
```

`out` is `[sample][voxel]` and caller allocated, allowing one workspace per
OpenMP worker. A NULL `sample_index` extracts all canonical rows; an explicit
list preserves its order and may repeat rows for future resampling workflows.
The extractor applies AFNI brick scale factors, supports mixed scalar brick
types, validates sample/brick/voxel indices, and reports non-finite values with
the responsible canonical sample and voxel. It reuses the same typed-brick
copying path as the existing all-brick RSA extractor.

Feature filtering belongs in the fold-aware MVPA layer when it depends on the
data. Simple spatial validity checks may occur during extraction, but variance
thresholding, standardization, nuisance regression, and learned selection must
be fitted on training data only.

## `thd_mvpa.[ch]`

`thd_mvpa` is the dataset-agnostic numerical library. It now contains fold-local
preprocessing, prediction metrics, centroid estimators, and the first reusable
classification workspace/fold runner. It consumes arrays and fold plans, not
AFNI datasets or command-line options.

Conceptual API objects:

```c
typedef struct {
    int nsample;
    int npredictor;
    int nresponse;
    const float *x;
    const float *y;
} THD_mvpa_problem;

typedef struct {
    int estimator;
    float lambda;
    int fit_intercept;
    int standardize_x;
    int standardize_y;
    int class_weight;
} THD_mvpa_params;

typedef struct {
    float *prediction;
    float *decision;
    float *coefficient;
    float *fold_score;
    float score;
    int nvalid;
} THD_mvpa_result;
```

The library should provide:

- explicit reusable workspaces;
- fit, predict, and score operations;
- train-fold preprocessing and test-fold application;
- single- and multi-output responses;
- fold-level predictions and diagnostics;
- no static mutable configuration;
- deterministic results independent of OpenMP thread count.

`THD_mvpa_workspace_alloc()` reserves compact train/test matrices,
preprocessing statistics, centroid storage, predictions, decisions, and score
buffers once for a fold plan and maximum feature count. Repeated
`THD_mvpa_run_centroid_cv()` calls perform no allocation. This is the intended
one-workspace-per-OpenMP-worker model for ROI and searchlight execution.

The runner fits preprocessing and class centroids separately inside every
training fold, predicts only held-out rows, and restores predictions and
per-class decisions to original sample order. Fixed train/test plans mark
training-only rows with prediction `-1` and `NAN` decisions. It reports ACC and
BACC per fold plus pooled held-out confusion, recall, ACC, and BACC. The pooled
score is descriptive; within-subject or other stratified group inference must
use the preserved fold/stratum structure rather than silently treating all rows
as independent. Fold-context errors identify missing training classes,
non-finite inputs, undefined classifier geometry, and invalid run-centering
applications.

## Cross-Validation

Cross-validation must be explicit and hard to misuse. Fold construction may
live in `thd_crossval.[ch]` or as a distinct section of `thd_mvpa`.

Required initial forms:

```text
leave-one-run-out
leave-one-subject-out
leave-one-explicit-fold-out
grouped K-fold
leave-one-item-group-out
fixed train/test partition
```

Candidate public syntax:

```text
-cv run
-cv subject
-cv item
-cv_column Fold
```

Checks must include:

- nonempty training and test sets;
- all required classes in every training fold;
- optional requirement for all classes in every test fold;
- no subject, run, or item group split across folds when that group defines
  the generalization target;
- warnings for severe fold-wise class imbalance;
- a clear failure when a model cannot be estimated.

For multi-subject within-subject analysis, the natural default is leave-one-run-
out separately within each subject, followed by group-level aggregation.
Across-subject prediction should be explicitly requested, normally with
leave-one-subject-out CV.

Hyperparameter selection must use nested grouped CV. Inner folds are generated
only from the current outer training set. No preprocessing or tuning decision
may inspect the outer test data.

Implemented by `THD_nested_foldplan_grouped()`: each outer fold owns an inner
leave-group-out or balanced grouped K-fold plan with compact indexes into that
outer fold's training list. Other outer strata and held-out samples are never
included. These plans are shared by worker workspaces and reused across spatial
units. `-show_folds` displays both nesting levels with original sample indexes.

## Leakage-Safe Preprocessing

Every learned transformation follows the same contract:

```text
fit transformation on outer-training samples
apply fitted transformation to outer-training and outer-test samples
fit model on transformed training samples
score transformed held-out samples
```

This applies to:

- feature centering and scaling;
- target centering and scaling;
- variance filtering;
- univariate feature selection;
- PCA or other projections;
- nuisance regression;
- missing-value imputation;
- class weighting;
- selection of lambda, SVM C, shrinkage, or feature count.

Useful initial options:

```text
-scale none
-scale feature_zscore
-scale sample_center
-scale run_center
```

`thd_mvpa.[ch]` implements these as a fitted, dataset-independent matrix
transformation. `feature_zscore` learns feature means and population standard
deviations from the outer-training rows only; numerically constant training
features are mapped to zero. `sample_center` subtracts each sample's own feature
mean and therefore does not pool information across samples. `none` is an exact
copy.

`run_center` is deliberately strict: each feature mean for a run is learned
from training rows, and applying the transform to a run absent from training is
an error. In particular, the program must not silently estimate the mean of a
held-out run from outer-test samples. This makes `run_center` incompatible with
ordinary leave-one-run-out unless run means were estimated independently of the
analysis samples. The eventual command help must make that scientific change
and restriction explicit. In multi-subject analyses, the caller must supply
subject-by-run group IDs rather than pooling identically named runs across
subjects.

The fold convenience API accepts a `THD_foldplan`, fits only
`plan->train[fold]`, and writes compact transformed training and test matrices
in fold-index order. Lower-level fit/apply calls support nested CV and reuse of
one fitted transform across multiple applications.

The help must explicitly state whether each operation is fold-local and whether
it changes the scientific estimand.

## Estimators

### First end-to-end classifiers

Nearest-centroid and correlation classifiers remain the safest first executable
milestone:

- fast and low-memory;
- easy to validate analytically;
- no optimizer;
- natural for searchlights;
- a clear conceptual bridge from RSA to predictive analysis.

Centroid distance choices:

```text
euclidean
correlation
cosine
```

`THD_centroid_fit()` computes each class centroid only from the supplied compact
training matrix and requires every integer class ID to be represented.
`THD_centroid_predict()` returns class IDs and optional per-class decisions.
All decisions use a common higher-is-better convention: negative squared
Euclidean distance, Pearson correlation, or cosine similarity. Exact ties are
resolved deterministically in favor of the lowest class ID. Correlation and
cosine fits reject degenerate class centroids, and prediction rejects a
degenerate test pattern, rather than silently turning an undefined similarity
into a class assignment.

The standalone model owns its fitted arrays. The fold runner instead refits the
centroid storage held by `THD_mvpa_workspace`, so repeated spatial analyses do
not allocate model or prediction buffers.

### First general solver: ridge

Ridge regression is implemented in `thd_mvpa.[ch]` for continuous decoding,
encoding, and class-coded classification:

- classification through one-hot targets and argmax prediction;
- scalar or multivariate continuous decoding;
- multi-output encoding.

Both numerical forms are implemented:

```text
primal: (X'X + lambda I)^-1 X'Y       when predictors <= samples
dual:   X'(XX' + lambda I)^-1 Y       when predictors > samples
```

The dual formulation is essential for neuroimaging decoding, where voxels
often greatly outnumber samples. Encoding commonly has a small external design
and many voxel responses, allowing one factorization to predict every voxel.

`THD_ridge_alloc()` reserves the smaller Gram-matrix factor and multi-output
coefficient/scratch storage. `THD_ridge_fit()` chooses primal or dual form for
the active dimensions, centers predictors and responses on training rows,
and shares one double-precision Cholesky factor across all response columns.
The objective is `sum((Y - intercept - X B)^2) + lambda * sum(B^2)`; the
intercept is not penalized. Lambda must be finite and strictly positive.
The CLI defaults to fixed `-lambda 1`. `-lambda auto` enables nested grouped
selection over `0.01,0.1,1,10,100`; `-lambda_grid L1,L2,...` supplies a distinct
positive candidate set and also enables tuning. A numeric lambda and a grid
cannot be combined. `-inner_cv_column` defaults to the outer `-cv_column`;
`-inner_cv_k K` optionally balances whole groups instead of leaving each group
out. Each outer training set must contain at least two inner groups, or at
least K for grouped K-fold. Thus ordinary two-run leave-run-out CV cannot tune
on Run without another explicitly chosen inner grouping unit.

Selection minimizes the mean of response-wise `inner SSE / inner baseline SSE`,
where both sums pool inner held-out rows and the baseline uses each inner
training response mean. Responses with zero baseline SSE are excluded. This
gives each informative response equal weight and keeps response units from
dominating selection. A single penalty is chosen per outer fold/spatial unit
and shared across outputs. Exact loss ties choose the largest penalty. If all
responses are uninformative, the same deterministic choice is recorded with
zero valid responses and NaN losses. The selected model and preprocessing are
then refitted on all outer training rows.

`PREFIX.tuning.tsv` records every outer-fold candidate loss, informative
response count, inner-fold count, and selection indicator. Fold-score audits
append `Lambda`, `InnerNFold`, `InnerNValid`, and selected `InnerLoss`. Fixed
lambda runs report their penalty with zero inner folds and NaN inner loss.
Analytic loss/selection tests, outer held-out perturbations, unequal inner-group
means/scales, sparse surface response IDs, and worker-count comparisons guard
these contracts. All tuning lives in the existing `thd_crossval` and `thd_mvpa`
modules; no new library module is introduced.

`THD_mvpa_ridge_workspace` reuses compact train/test matrices, a preprocessor,
the ridge solver, predictions, baselines, and scorer scratch across spatial
units. Untested fixed-split rows retain NaN predictions. Each response receives
pooled held-out and per-fold Pearson, Spearman, cross-validated R2, MAE, and
RMSE scores. No target statistics are estimated from held-out responses.

### Descriptive full-data ridge coefficients (completed 2026-10-15)

Why: predictive scores establish generalization, while descriptive refits let
users inspect the fitted predictor/response relationships. They answer
complementary questions and must remain distinct scientific outputs.

| Part | Why are we doing this? | Success metric / validation |
|---|---|---|
| Numerical coefficient export | Expose the solver's fitted slopes and unpenalized intercept without duplicating numerical algorithms | Analytic multi-response slopes/intercepts agree; exported primal/dual slopes agree and reconstruct model predictions; failed fits cannot export stale models |
| Independent full-data refits | Let users request a descriptive model without contaminating held-out estimation or inference | Byte-identical prediction, continuous-result and inference tables, numerically identical score maps with/without exports; unchanged tuned predictions and candidate-loss audits |
| Spatial maps and identifiers | Locate encoding effects and decoding/class weights while retaining distinct overlapping models | Known-answer volume map values, sparse mesh IDs 0/2/5/9, explicit feature/response-to-map/brick mapping, and identical output values across one/two workers |
| Scaling and interpretation | Make weights interpretable in their actual predictor/response spaces | Known population-SD scaling and constant-feature flags, sample-centered analytic slopes, full-data run-group centers and sample/group mapping; explicit descriptive labels and no statistic codes |
| Penalty/output validation | Make the refit reproducible and protect held-out outputs | Tuned CV requires explicit positive `-coef_lambda`; fixed penalties default to `-lambda`; invalid modes/penalties, file collisions, and existing coefficient outputs are rejected |

`-save_coefficients CPREFIX` runs a separate final descriptive pass after the
observed CV/permutation outputs have been finalized. It fits preprocessing and
ridge on **all canonical samples**, pooled across CV strata, without producing
full-data performance scores. The numerical library adds
`THD_ridge_coefficients()` to export `[feature][response]` slopes and response
intercepts in the exact predictor space passed to the solver.

`-coef_lambda L` sets a finite positive refit penalty, with objective
`sum_squared_error + L * sum_squared_slopes`; intercepts are unpenalized.
For fixed CV penalties it defaults to `-lambda`. With `-lambda auto` or
`-lambda_grid`, an explicit `-coef_lambda` is mandatory. No outer-fold penalty
is silently reused or averaged, and this option never affects CV tuning.
Full-data penalty selection is outside this slice.

Outputs:

- `CPREFIX.coefficients.tsv` records spatial unit, feature and response IDs,
  slope, response intercept, refit penalty, total sample count, map prefix and
  zero-based map brick. Neural IDs are `voxel_INDEX` or `node_MESHID`; named
  numeric columns and class labels are preserved. Classification responses
  are the existing 0/1 class indicators, with scores that are not probabilities.
- `CPREFIX.preprocessing.tsv` records feature centers/scales and constant
  flags, or each full-data subject/run group mean and original sample/group
  membership. Weights satisfy `y = intercept + preprocessed_X * coefficient`.
  `feature_zscore` uses population SD and sets constant features to zero;
  responses stay in their original units. `sample_center` subtracts each
  sample's mean across that unit's predictors; `run_center` subtracts the
  recorded group feature means. These are transformed-space weights, not
  automatically converted raw-space weights.
- Encoding without searchlights writes `CPREFIX` with one coefficient brick
  per design column and `REFIT_INTERCEPT`, placing values at response
  voxels/nodes. Nonoverlapping ROI fits share the same map; voxel encoding
  avoids creating a dataset for every singleton response.
- Decoding and ridge classification write `CPREFIX.unit_ID` per ROI or
  searchlight, with one coefficient brick per target/class at neural feature
  locations. Searchlight encoding also keeps one dataset per unit, with design
  coefficient bricks and an intercept. Overlapping neighborhoods are never
  merged or center-painted. Surface map names end in `.niml.dset` and retain
  sparse mesh IDs. Values outside a unit are zero; use the TSV feature/response
  list to distinguish excluded locations from fitted zero weights.

Map labels begin with `REFIT_`, carry no inference statistic codes, and retain
command history. Searchlight exports can produce many full-domain datasets
and large tables, so they remain opt-in. Coefficient and score dataset paths
must differ, audit paths cannot collide, and existing outputs respect AFNI's
overwrite policy.

What fMRI users can now do: inspect how each external design variable weights
neural responses, or how neural features contribute to predicting continuous
targets and categorical class indicators, while retaining their original
held-out validation and inference. These predictive weights depend on
regularization, scaling and correlations among predictors. They are **not
activation/Haufe maps, independent voxel effects, or coefficient significance
tests**. Saved-model application and Haufe transforms remain future scope.

Example additions to an existing ridge analysis:

```text
-save_coefficients descriptive_weights
# With tuned CV, also specify a descriptive refit penalty:
-coef_lambda 1
```

Validation: `thd_mvpa:unit`, `3dMVPA:coefficients`, and the extended
`3dMVPA:surface-searchlight` test cover the contracts above; the existing
classification and ridge integration suites preserve prior behavior.

### Additional classifiers

Ridge classification is implemented with `-mode classify -method ridge`.
Centroid remains the classification default. `-distance` is centroid-only and
is rejected for ridge classification. Ridge uses the same positive fixed
`-lambda` or grouped nested `-lambda auto`/`-lambda_grid` options as continuous
ridge, without introducing another library module.

Each class is represented by its own 0/1 response column. Slopes are penalized,
the training-fitted intercept is not, and the largest held-out response score
determines the class. Exact ties choose the lowest encoded class ID (class
order is preserved in the audit header). DecisionClass columns are raw
regression scores, not probabilities; values are not clipped or calibrated.
Fits are unweighted, so constant patterns predict training class frequencies.
Balanced accuracy remains the default reported statistic; reporting BACC does
not imply class-weighted fitting.

Nested selection minimizes mean class-wise SSE/training-mean-baseline SSE on
inner held-out rows, not accuracy or BACC. This reuses the established
multi-response squared-loss objective and largest-penalty tie rule. Every outer
and inner training set must contain every class; invalid splits or null label
arrangements fail explicitly rather than dropping classes or permutations.
All preprocessing, class-coded responses, tuning and models are refitted under
each synchronized label null, with observed predictions and tuning restored
afterward. Existing empirical above-chance/P/Z/Q/FWE correction is unchanged.

`THD_mvpa_classification_set_ridge()` configures an existing
`THD_mvpa_workspace`; `THD_mvpa_run_classification_cv()` and its permutation
counterpart share confusion matrices, recall, pooled/fold metrics, and original
sample ordering across centroid and ridge. Existing centroid entry points stay
available. Ridge fold audits append Lambda/InnerNFold/InnerNValid/InnerLoss;
automatic `PREFIX.tuning.tsv` records all candidates for nested runs. Centroid
table layouts and default behavior remain unchanged.

Known-answer tests cover three-class decisions and losses, outer-test label
isolation, wide rank-deficient solves, training class frequencies, constant
pattern ties, missing inner classes, and malformed permutations. CLI tests
cover fixed/nested ROI fits, volume and sparse-surface nulls, restored observed
audits, exact surface probabilities, and thread-count reproducibility.

### Diagonal and shrinkage LDA

Both are implemented in the existing `thd_mvpa.[ch]` module and exposed as
`-mode classify -method lda_diag|lda_shrinkage`. Centroid remains the default.
The shared runner refits preprocessing and all LDA statistics inside each
outer training set, including every label permutation. At least one residual
degree of freedom is required: every outer training set must have more rows
than classes, and all classes must be present. Invalid fits fail explicitly.

Let R contain training rows minus their class means. The pooled covariance is
`S = R'R / (ntrain - nclass)`. Equal priors are the default; `-class_prior
empirical` uses training class proportions. Decisions are centered linear
discriminants `((x-global_mean)' inverse(Sreg) (class_mean-global_mean)) -
0.5*((class_mean-global_mean)' inverse(Sreg) (class_mean-global_mean)) +
log(prior)`, not probabilities. Centering subtracts only a common class score
relative to the usual uncentered discriminant. Exact ties choose the lowest
encoded class ID. Balanced accuracy is still the default reported metric.

`-variance_floor F` defaults to `1e-6` and accepts `[1e-12,1]`. The variance
reference T is the mean diagonal of the pooled within-class covariance. When
that is exactly zero, use mean overall training variance (denominator ntrain-1),
then 1 for wholly constant data. This preserves noiseless class separation
without an undefined inverse. Diagonal LDA uses `max(S[j,j], F*T)` per feature.
It does not form a full covariance or allocate quadratic feature storage.

Shrinkage LDA uses `(1-A)*S + max(A,F)*T*I`, where `-shrinkage A` defaults to
0.1 and accepts `(0,1]`. This is fixed, configurable shrinkage toward a spherical
target, not automatic shrinkage estimation or CV-selected shrinkage. A=1 is
purely spherical. The F lower bound protects the diagonal regularizer even
when A is extremely small. Feature scaling affects the spherical target, so
the chosen preprocessing is part of the model specification. Equal priors
and globally constant features produce deterministic ties; empirical priors
give the training-majority class for wholly constant patterns.

Fits normalize by sqrt(T) for numerical stability and use double precision.
The shrinkage model factors a feature-space matrix when features <= training
rows, or a sample-space residual Gram matrix and Woodbury solve otherwise.
Factor storage is quadratic in min(ntrain,nfeature), not always nfeature.
Repeated spatial/null fits allocate no new model or prediction buffers.

`THD_lda_alloc/fit/predict/free()` provide the dataset-independent estimator;
`THD_mvpa_classification_set_lda()` configures the shared classification
workspace. `-distance` is centroid-only, `-lambda`/inner-CV options are ridge-only,
and `-shrinkage` is restricted to lda_shrinkage. Prediction audit headers name
the method, prior, shrinkage and relative floor. LDA fold audits append
`Shrinkage`, `CovarianceScale` (T), `VarianceFloor` (absolute F*T), and
`CovarianceForm` (diagonal/primal/dual). Observed diagnostics are restored after
the synchronized label nulls. Existing P/Z/Q/FWE spatial families are unchanged.

Known-answer tests cover correlated 2x2 inverses, wide rank-deficient
sample-space solves, spherical limits, three-class discriminants, prior effects,
variance floors, zero residual variance, rejected zero residual degrees of
freedom, diagonal feature-scaling invariance away from active floors, large/small
measurement units, and held-out label/feature isolation.
ROI/volume/sparse-surface CLI tests cover both estimators, explicit priors,
restored prediction/fold audits, exact blocked surface p-values, and worker-count
reproducibility. No new library module or external runtime dependency is added.

### Native linear SVM (completed 2026-10-15)

Why: a margin-based classifier gives fMRI users another explicit predictive
model for categorical condition information, while retaining the same sample,
CV, spatial and inference definitions as the other classification methods.
The implementation is native C in `thd_mvpa.[ch]`; it adds no external runtime,
model format, library module, kernel or saved-model application path.

| Part | Why are we doing this? | Success metric / validation |
|---|---|---|
| Native squared-hinge solver | Provide a reproducible margin-based estimator with an explicit objective | Analytic symmetric margins, shifted-input slope/intercept solution, unequal-class constant-pattern solution, wide/rank-deficient matrices, and three-class argmax tests pass |
| Convergence and constant-pattern reduction | Prevent unfinished optimization from being reported as a valid scientific result, and handle flat brain regions efficiently | Every successful class fit passes its final absolute projected-gradient tolerance; insufficient sweeps fail and invalidate the model; identical nonzero patterns at large C match a closed-form solution in one step |
| Training-only fitting and grouped C selection | Choose regularization without using outer test patterns or labels | Perturbing outer held-out features/labels leaves the corresponding training transformation and selected model unchanged; tied inner scores select the smaller C; missing inner training classes fail clearly |
| Spatial output, decisions and audits | Make SVM results inspectable in the existing AFNI workflow | ROI scores, analytic volume margins (40/41), multiclass decision columns, selected C/candidate losses and solver diagnostics pass; sparse surface maps retain mesh IDs 0/2/5/9 |
| Fully refitted blocked nulls | Test the actual prediction procedure, including tuning, and preserve the observed scientific result | Null scores match explicit complete refits; exact blocked fixtures give PVAL/PFWE=.25; prediction/fold/tuning tables restore byte exactly after nulls; maps and audits agree across one/two workers |

Use `-mode classify -method svm`. Each class has a one-versus-rest binary
problem with signs +1 for that class and -1 for every other class:

```text
minimize  0.5 * (sum(weight^2) + intercept^2)
          + C * sum(max(0, 1 - sign * (X*weight + intercept))^2)
```

The intercept is **penalized**, equivalent to an appended constant feature of
one. The loss is a sum over training samples, not a sample-count-normalized
mean; no class weighting is applied. The default `-scale feature_zscore`
fits predictor means/population SD inside the training split. All existing
scaling modes remain available with their original training-only contracts.
Raw class margins use deterministic argmax; exact ties select the first
encoded class. Margins are not probabilities or calibrated confidence values.

`THD_svm_alloc/fit/predict/free()` own reusable per-worker model storage.
The solver uses deterministic cyclic dual coordinate descent with
nonnegative dual variables, diagonal term `1/(2*C)`, and maintained primal
weights/intercepts. No kernel/Gram matrix or per-fit allocation is needed.
After each sweep it recomputes every projected dual gradient using the final
weights. It accepts the model only when the largest absolute residual across
samples/classes is at most `-svm_tol` (default 1e-6, allowed [1e-12,1e-2]).
`-svm_max_iter` limits sweeps per class (default 10000); a failure invalidates
the model and reports the fold, class and remaining residual. Increasing the
limit, rescaling inputs or reducing C can address difficult fits.

Identical training patterns have an exact scalar reduction. With
`L = 1 + sum(pattern^2)`, the class margin is
`(npositive - nnegative) / (n + 1/(2*C*L))`, and its minimum-norm augmented
weight is `[pattern,1] * margin/L`. Both signs are present, so the solution
lies within the quadratic hinge region. This handles flat regions without
slow cancellation of large opposing dual variables, while checking the same
optimality residual.

`-svm_c C` supplies a fixed finite positive cost (default 1).
`-svm_c auto` uses the candidate grid 0.01,0.1,1,10,100;
`-svm_c_grid C1,C2,...` supplies distinct positive candidates and also enables
tuning. A numeric `-svm_c` cannot accompany a grid. Larger C gives margin
violations more weight relative to regularization; exact validation ties
choose the **smallest C**. Numerical precision and convergence checks apply
to every candidate rather than silently dropping failed fits.

Grouped tuning reuses `THD_nested_foldplan`, `-inner_cv_column` (default outer
CV grouping column) and `-inner_cv_k`. Each inner fold is restricted to the
outer training rows. It refits preprocessing and SVM from scratch and pools
held-out inner predictions. Selection minimizes `1-accuracy` or
`1-balanced_accuracy`, matching the chosen primary `-metric`. Every outer
and inner training set must contain all classes. The training-only model is
then refitted using the selected C before predicting the outer test rows.

The existing prediction audit retains one decision column per class, original
sample metadata and spatial IDs. Fold audits append:

```text
C  InnerNFold  InnerLoss  SolverIterations  SolverResidual
```

Fixed fits have zero inner folds and NaN inner loss. Solver diagnostics are
the largest sweep count and final residual over that outer fold's binary
class models. Tuned runs always write `PREFIX.tuning.tsv` with unit, outer
fold/stratum, inner fold count, candidate C, pooled loss and selection flag.
All preprocessing, C selection and fits repeat for each synchronized blocked
label permutation; the observed run is restored afterward. The existing
PVAL/Z/Q/FWE maps and correction-family definitions are retained.

Example additions to a classification command:

```text
-method svm -svm_c 1
# Or grouped training-only cost selection:
-method svm -svm_c auto
```

What fMRI users can now do: test whether held-out task conditions can be
predicted from ROI patterns or volume/geodesic surface neighborhoods using
a native linear SVM; inspect the class margins and fitted penalty for each
fold; and obtain blocked permutation inference for the complete analysis.
The same metadata and spatial definitions can be reused when running the
other estimators. Calibrated probabilities, kernels, class weighting, SVM
coefficient exports and saved-model application are outside this slice.
`-save_coefficients` remains ridge-only, preserving item 18's contract.

Validation: all six MVPA gates pass (`thd_mvpa:unit`, classification, ridge,
coefficient, SVM and surface integration). Both CMake and Make continue to
compile the existing `thd_mvpa` source; no new build dependency is introduced.
Roadmap items 1–19 are complete. Items 20–40 are proposed expansions; item 20
is next pending user direction. Further implementation follows the development
cadence recorded at the top of this plan.

Estimator proposals 28–32 cover elastic-net regression, regularized logistic
classification, PLS, linear SVR and Huber regression; see Phase 6 for their
mode-specific contracts and acceptance metrics. Phase 7 adds proposals 33–40
for generalized tuning, encoding extensions and generalization workflows.
These remain proposals, not implemented capabilities. Kernel models are
explicitly later-priority proposals, not current support. Real-time support,
legacy SVM model compatibility and transductive modes remain outside scope.

## Mode-Specific Behavior

### `-mode classify`

Required target:

```text
-target_column Condition
```

Initial metrics:

```text
accuracy
balanced_accuracy
per_class_recall
confusion_matrix
```

Later metrics:

```text
AUC
log_loss
dprime for well-defined binary designs
```

Balanced accuracy should be the default because ordinary accuracy can hide
class imbalance. Prediction tables should include the observed class,
predicted class, decision scores, fold, subject, run, and original table row.

`THD_mvpa_score_classification()` returns accuracy, balanced accuracy,
per-class recall, and a truth-row by prediction-column confusion matrix. A
class absent from the scored observations has `NAN` recall and is excluded from
the balanced-accuracy average; fold validation separately controls whether a
missing test class is permitted. Recall and confusion matrices are diagnostic
tables, not scalar spatial-map metrics.

### `-mode decode`

Required target:

```text
-target_column Rating
-target_columns Rating,Arousal
```

Initial metrics:

```text
pearson
spearman
cross_validated_r2
MAE
RMSE
```

For multivariate targets, retain per-target scores and offer an explicitly
defined aggregate rather than silently averaging incomparable scales.

Implemented: `-target_column` or comma-separated `-target_columns` selects
finite numeric columns, `-method ridge` is the continuous default, and
`-metric cross_validated_r2` is the default score. Each target has separate
SCORE/PEARSON/SPEARMAN/CVR2/MAE/RMSE/NVALID map bricks; NFEAT and NFOLD follow.
`PREFIX.continuous.tsv` preserves all metrics. Optional `-save_predictions`
and `-save_fold_scores` retain response identity and sample/fold provenance.
There is no implicit average across decoding targets.

`THD_mvpa_score_continuous()` returns Pearson correlation, tie-aware Spearman
correlation, cross-validated R-squared, MAE, and RMSE. Its R-squared denominator
uses caller-supplied predictions made without held-out outcomes, normally the
training-fold mean repeated over that fold's test rows. It must never use the
held-out target mean. Undefined correlations or a zero R-squared denominator
produce `NAN`; non-finite inputs are errors. The scorer accepts caller-owned
rank/index scratch space so repeated ROI and searchlight scoring does not
allocate memory.

### `-mode encode`

Required design:

```text
-design_columns Visual,Semantic,Motion
```

Encoding has two scientifically distinct spatial outputs:

```text
-score voxel
-score pattern
```

- `voxel` fits the external design to each masked voxel and writes a normal
  voxelwise predictive map. A searchlight is unnecessary.
- `pattern` predicts the multivoxel response in each ROI/searchlight and writes
  one score per spatial analysis unit.

Both paths are implemented with `-design_columns` and `-score voxel|pattern`
(default voxel). Voxel mode uses positive mask voxels/nodes as singleton units;
pattern mode uses ROIs or volume/geodesic surface searchlights and averages
each metric equally over responses where it is defined. NVALID reports the
number of responses with a defined primary metric, and NResponse is retained
in the continuous table. Undefined metrics stay NaN in text and become zero
in maps for AFNI compatibility; consult NVALID for the primary score. Encoding
prediction/fold audits preserve each individual voxel or mesh-node response.
Searchlight workers own independent reusable ridge workspaces; ordered audits
are deterministic across thread counts.

Continuous permutation inference is implemented as described below.
Descriptive full-data coefficient maps are implemented with
`-save_coefficients`, as described above.

This distinction prevents encoding from being forced into a decoding-shaped
output. Optional full-data refit coefficients may be written as one map per
design column, clearly labeled as descriptive refits rather than held-out
performance.

## Permutation Inference

Classification, decoding, and encoding implement:

```text
-nperm N
-seed S
-permute_within COLUMN
```

`-block_column` is an alias for `-permute_within`. The block column is explicit
and may be a subject, run, fold, or user-supplied composite exchangeability
unit. `none` requests an unrestricted null deliberately; inference never
silently assumes that all rows are exchangeable. The identity is included,
small groups are enumerated exactly, and classification relabelings are
checked for complete class coverage in every outer training fold before any
spatial analysis begins. Continuous exact enumeration uses row-index
permutations; equal-valued rows are not collapsed.

Continuous decoding jointly shuffles all numeric target columns as whole
rows; encoding jointly shuffles all design columns as whole rows. Neural rows,
folds, group IDs, and other sample metadata stay fixed. Attribute covariance is
preserved. Rows must be exchangeable within the selected blocks: this is not a
residual/nuisance-adjusted or autocorrelation-preserving time-series null.

Folds are built once. For every permutation:

```text
relabel the target within allowed exchangeability units
repeat fold-local preprocessing
repeat inner-CV tuning, when enabled
fit every outer-training model
predict every held-out sample
aggregate the requested statistic
```

The same permutation slot must be used at every spatial location. This makes a
synchronized mapwise max-statistic null possible and preserves thread-count
reproducibility.

The implemented classification bucket is:

```text
SCORE ABOVE_CHANCE PVAL ZSCR QVAL PFWE ZFWE NFEAT NFOLD
```

`ABOVE_CHANCE` is the observed selected metric minus its empirical null mean.
PVAL is its upper-tail identity-inclusive permutation probability, QVAL is BH
over the complete ROI/searchlight family, and PFWE/ZFWE use the synchronized
maximum score at every relabeling. ZSCR and ZFWE receive AFNI FIZT stat codes.
Without `-nperm`, output remains the original compact SCORE/NFEAT/NFOLD bucket.

Continuous output retains its seven observed metric/count bricks per target,
then adds these eight inference bricks per target, followed by NFEAT/NFOLD:

```text
NULL_MEAN ABOVE_NULL PVAL ZSCR QVAL PFWE ZFWE NNULL_VALID
```

`PREFIX.inference.tsv` records scores, null means, inference, finite-null counts,
and total identity-inclusive NPERM. One spatial correction family is used per
decoding target, or one for the selected encoding voxel/pattern statistic.
Targets are not jointly corrected. Encoding pattern scores use the existing
equal mean across responses with defined metrics in each relabeling.
Observed prediction/fold/tuning audits are restored after the null runs.

Pearson, Spearman, and CVR2 use upper tails in original metric units; MAE/RMSE
use lower tails. ABOVE_NULL is observed-minus-null-mean for higher-is-better
metrics and null-mean-minus-observed for errors. Positive always means better
than the empirical null. Correlations are not Fisher-transformed or forced to
have a zero-centered null. Max-statistic correction uses the same direction.

Undefined observed metrics receive P/Q/PFWE=1 and Z/ZFWE=0. For a defined
observation, undefined null metrics count as exceedances and conservatively
enter the synchronized family maximum; probabilities always use all NPERM
slots. NULL_MEAN uses finite scores, counted in NNULL_VALID. BH includes the
complete spatial family, including undefined observations at p=1. Undefined
scores remain NaN in text and zero in maps; consult NVALID/NNULL_VALID.

The user must be able to define exchangeability through subject, run, fold, or
an explicit block column. Invalid permutations should fail rather than silently
crossing nonexchangeable boundaries.

Score-specific null handling:

- classification: report accuracy and accuracy-minus-chance; infer against
  label permutations rather than a symmetric theoretical distribution;
- correlation: retain the metric's original units and use its empirical null;
- R2: use its permutation distribution rather than assuming symmetry around
  zero;
- paired model comparisons: drive both models with identical folds and
  permutation slots.

Reuse `THD_perm_result`, `THD_perm_emp_pvalue()`, `THD_perm_signed_z()`,
`THD_max_accum()`, and `THD_bh_fdr_masked()` for final map inference.

## Spatial Analysis and Output

### Searchlight

Volumetric searchlights use `SPHERE`, `RECT`, `RHDD`, or `TOHD` neighborhoods.
Surface searchlights use a plain positive geodesic radius and a SUMA-readable
mesh. An explicit surface mask and every sample dataset must have the same
dense or sparse node-row domain; without a mask, every node represented in the
sample datasets becomes a center and eligible feature. Sparse outputs retain
their mesh node indexes, and audit tables report mesh node IDs rather than
internal dataset rows.

For each center:

```text
collect in-mask neighboring voxels
extract sample x voxel matrix
run the selected CV backend
assign scalar results to the center voxel
```

### ROI

For each atlas parcel:

```text
extract sample x voxel matrix
run the selected CV backend
write a row to an ROI table
paint requested scalar results over the parcel
```

Default AFNI bucket contents should be small and predictable. Depending on
mode and requested inference, useful labels include:

```text
SCORE
PVAL
ZSCR
QVAL
PFWE
ZFWE
NFEAT
NFOLD
```

Mode-specific multi-value results, confusion matrices, fold diagnostics, and
row-level predictions belong in text/1D tables rather than hundreds of map
bricks.

Essential audit options:

```text
-show_table
-show_samples
-show_folds
-save_predictions FILE
-save_fold_scores FILE
```

The prediction table retains the spatial-unit identity and feature count,
original table row, sample label, subject, run, item, input-fold metadata, CV
test fold and stratum, observed and predicted class labels, correctness, and
one decision value per class. Rows
that a fixed split never tests remain present with `NA` predictions. The fold
table retains each spatial unit, fold and stratum labels, train/test counts,
ACC, and BACC. Both files are opt-in because row-level searchlight output can
be very large. These tables are part of the scientific result, not merely
debugging output.

## Aggressive 3dRSA Library Extraction

`3dRSA.c` is now large enough that waiting for duplicate code to appear in
`3dMVPA` would be counterproductive. Refactor behavior-preserving scientific
code before or alongside the new program.

### Current extraction assessment (2026-10-13)

Keep the library count controlled. `thd_rsa.[ch]` is the one coherent RSA
science library, not the start of separate libraries for every RSA feature.
Its first completed slice owns dataset-independent condition centering and RDM
construction, exact similarity/distance triangle conversion, descriptive mean
RDM construction, and reusable subject-RDM storage for second-order geometry.
`3dRSA` now supplies extracted condition patterns or crossnobis matrices to
that core and reuses one workspace per OpenMP worker. The memory preflight
accounts for this storage.

Later coherent groups may move Mantel/regression statistics, fitted-model CV,
commonality, and their reusable workspaces into the same library. Leave
`3dRSA` option policy, dataset/ROI extraction, residual whitening, scientific
null definitions, and output layout in the program. Resume MVPA development
before the next RSA extraction group.

Completed: `thd_patterns.[ch]` now contains the two small AFNI spatial helpers
that were duplicated between the programs: exact dense/sparse surface
node-domain comparison and copying a sparse surface node list into an output
dataset. They live beside ROI/searchlight construction and painting and did not
create another library.

The present `thd_permute` plus `thd_mapinfer` boundary is sufficient. Relabeling
sets, empirical result finalization, BH, and max accumulation are shared; each
program's scientific null, family definition, and output-brick layout should
remain program-owned. Likewise, do not extract CLI parsing, help text, audit
schemas, or a generic pipeline object merely because both programs have them.
Reconsider a shared progress/memory-preflight helper only after `3dMVPA` gains a
second concrete implementation with the same contract.

### Move to `thd_sampleset` or general table helpers

- repeated subject/run indexing;
- run-level factors and contrasts;
- label-level encoding;
- generic run normalization;
- selected dataset/brick mapping.

### Keep or move into `thd_patterns`

- ROI/searchlight construction and filtering;
- dataset-domain validation;
- selected-sample extraction;
- pattern centering helpers;
- analysis-unit painting.

### Create `thd_rsa.[ch]`

Move higher-level RSA science out of both `3dRSA.c` and `thd_simmatrix`:

- RSA model, contrast, commonality, and fitted-model descriptors;
- run-resolved Mantel tests and RDM regression;
- subject and condition CV for fitted RDM mixtures;
- subject RDM and second-order RDM construction;
- seed-model construction;
- RSA-specific subject and condition bootstrap procedures;
- temporal-null RSA routines;
- reusable RSA workspaces.

`thd_simmatrix` should remain the lower-level matrix/RDM arithmetic layer.
`thd_rsa` should own workflows built from those matrices.

### Move generic utilities carefully

- generic percentiles and sample summaries can live in an appropriate stats
  helper if AFNI does not already provide them;
- phase preparation/drawing belongs with `thd_phasefft` or `thd_permute`;
- progress rendering is worth a shared helper only if both programs use it;
- option parsing, usage text, program-specific error policy, output labels, and
  provenance remain in the front-end programs.

The goal is not an artificially tiny `3dRSA.c`. Its help and output definitions
will remain substantial. The goal is that it contain no private estimator, CV,
RDM-construction, bootstrap, or permutation algorithm that another program
could reasonably reuse or independently test.

## Testing Strategy

### Development-only reference checks

Independent MVPA implementations may be used from a developer environment to
compare numerical results, searchlight semantics, failure cases, and potential
feature directions. These are private development oracles, not dependencies:
distributed source, tests, help, provenance, and output metadata must not name,
import, or require an external MVPA toolbox. Checked-in tests should use
analytic answers or self-contained fixtures so every AFNI build remains
reproducible without a Python environment.

### Library unit tests

- direct sample tables and X-matrix-expanded sample tables producing equivalent
  `THD_sampleset` organization;
- selected-sample extraction preserving row and voxel order across datasets and
  scaled brick types;
- flat and compact sample-table expansion;
- dataset selector and brick indexing;
- categorical level encoding;
- grouped fold construction;
- rejection of leaking or incomplete folds;
- training-only scaling with known held-out values;
- primal and dual ridge agreement;
- single- and multi-output ridge solutions;
- centroid and correlation predictions;
- balanced accuracy and regression metrics;
- within-block permutation behavior;
- deterministic results across thread counts.

### Integration tests

- direct-table and X-matrix-expanded ROI classifications producing identical
  known ACC/BACC results and three-brick AFNI outputs;
- volumetric searchlights producing known center-wise ACC/BACC, neighborhood
  size, and fold-count maps without an impractically large location table,
  with numerically identical maps and byte-identical opt-in audit tables across
  one and two OpenMP workers;
- sparse surface searchlights producing known ACC/BACC/NFEAT/NFOLD values and
  mesh-node IDs, with identical explicit-all-node-mask and mask-free results;
- clear failures for mismatched surface node-row domains, node IDs outside the
  mesh, and shaped rather than geodesic surface-radius specifications;
- synthetic separable classes with a known spatial signal location;
- continuous target encoded in a known searchlight;
- encoding design with known responsive voxels;
- ROI result matching a direct C-level calculation;
- [x] chance data producing sane identity-inclusive permutation p-values;
- [x] max-FWE behavior under a synchronized null, identical across thread counts;
- prediction-table rows matching original metadata;
- clear failures for missing classes, impossible folds, grid mismatch, and
  non-finite input.

### 3dRSA regression protection

Library extraction must be checked against existing 3dRSA oracle and synthetic
tests. Extraction commits should not change numerical behavior except where a
separately documented bug is intentionally corrected.

On 2026-10-13, the first `thd_rsa` extraction was also replayed on the archived
22-subject Paranoia validation analysis (268 atlas ROIs, two behavioral models,
5,000 synchronized label permutations, fixed seed). The refactored run produced
a byte-identical `.rsa.1D` table and byte-identical decompressed six-brick AFNI
payload; geometry, labels, datum types, and statistic codes also matched. Only
expected provenance fields such as creation time, dataset ID, version, command
history, and output prefix differed in the `.HEAD` file.

## Implementation Sequence

### Phase 0: behavior-preserving shared-library work

1. [x] Bring `thd_xmat.[ch]`, its `thd_datatable` parsing improvements, and its
   tests from `SRMthings` into the shared-library baseline.
2. [x] Add `thd_sampleset.[ch]` and test both direct and X-matrix-expanded tables.
3. [x] Generalize selected-sample extraction in `thd_patterns` after the
   numerical fold runner is ready.
4. [x] Extend `thd_patterns.[ch]` with exact dense/sparse surface-domain
   comparison and sparse output node-domain copying; use both from `3dMVPA`
   and `3dRSA`.
5. [x] Add `thd_rsa.[ch]` and first move dataset-independent RDM construction
   and its reusable workspaces. Move further 3dRSA scientific helpers only in
   coherent groups.
6. [x] Keep 3dRSA output numerically identical for the first extraction group,
   verified by direct known-answer tests, the existing fast/full oracles, and a
   byte-identical 22-subject real-data validation replay.

### Phase 1: MVPA numerical core

1. [x] Add the general fold-plan object.
2. [x] Add leakage-safe centering/scaling.
3. [x] Add prediction metrics.
4. [x] Add nearest-centroid and correlation classifiers.
5. [x] Add reusable workspaces and the end-to-end fold runner.
6. [x] Add primal and dual multi-output ridge.
7. [x] Test current algorithms with small matrices having known solutions,
   primal/dual parity, rank-deficient predictors, training-only intercepts and
   baselines, and held-out response perturbations.

### Phase 2: first `3dMVPA` vertical slice

Deliver a scientifically useful classification program:

```text
3dMVPA
    -mode classify
    -dataTable @samples.tsv
    -target_column LABEL
    -cv_column RUN
    -mask MASK
    -searchlight SPHERE(R)
    -method centroid
    -metric balanced_accuracy
    -prefix PREFIX
```

Include `-show_samples`, `-show_folds`, row-level predictions, ACC/BACC/NFEAT/
NFOLD outputs, and strong design validation. The vertical slice must accept
either a direct one-row-per-sample table or a subject-level table with `Xmat`;
both must enter the same analysis path after table expansion.

Current Phase 2 progress:

1. [x] Add the native `3dMVPA` command shell with explicit classification,
   target, CV, sample-role, preprocessing, centroid-distance, metric, mask, and
   prefix options.
2. [x] Converge direct and X-matrix tables on one `THD_sampleset` analysis path.
3. [x] Run atlas/binary-mask ROIs through fold-local centroid classification
   and write ACC/BACC ROI tables plus SCORE/NFEAT/NFOLD AFNI maps.
4. [x] Expose canonical table, resolved sample, and fold audits on stdout.
5. [x] Add volumetric searchlight execution using the same extraction and
   workspace objects.
6. [x] Add row-level prediction and per-fold score files, retaining original
   sample metadata and class labels.
7. [x] Expose geodesic surface searchlights through the same analysis-unit
   interface, with explicit surface/domain validation and sparse-node tests.

### Phase 3: inference and fitted estimators

1. [x] Add blocked label permutations and PVAL/Z/Q/FWE outputs.
2. [x] Add continuous ridge decoding, including multiple numeric targets,
   per-target maps, and row/fold audit tables.
3. [x] Add nested grouped lambda selection and per-candidate/per-fold audits.
4. [x] Add continuous permutation inference with model/tuning refits, exact
   blocked known-answer checks, lower-tail errors, undefined-score safeguards,
   observed-audit restoration, and volume/surface thread reproducibility.
5. [x] Add ridge classification with fixed/nested penalties and full label
   inference through the shared classification runner.
6. [x] Add diagonal and shrinkage LDA, training-only covariance/prior fits,
   explicit floors, high-dimensional solves, fold diagnostics, and full label
   nulls through the existing classification runner.
7. [x] Add native linear squared-hinge SVM, deterministic one-versus-rest
   decisions, verified convergence, grouped nested C selection, fold/tuning
   diagnostics, and fully refitted blocked label inference.

### Phase 4: encoding

1. [x] Add voxelwise ridge encoding using the multi-output solver.
2. [x] Write requested descriptive full-data ridge coefficient maps, with
   feature/response IDs, intercepts, explicit preprocessing and refit penalties,
   and independent CV/inference outputs.
3. [x] Add ROI/volume/surface searchlight pattern prediction.
4. [x] Add encoding prediction and fold tables with response voxel/node IDs.
5. [x] Add encoding-appropriate joint design-row inference for exchangeable
   samples; temporal/residual nulls require separate future contracts.

### Phase 5: proposed roadmap expansion (items 20–27)

The user approved adding the following proposals on 2026-10-15. Implementation
has not started. Each slice must preserve existing CV/inference behavior,
record its validation evidence, and conclude with what fMRI users can now do
and how it helps their analyses. Prompt the user before the next slice.

| Item | Why are we doing this? | Completion metric | Intended benefit for fMRI users |
|---|---|---|---|
| 20. Classification and subject-summary audits | Reveal condition confusions and subject differences hidden by pooled performance | Exported confusion matrices reproduce ACC/BACC; per-class recall matches known answers; subject summaries match independent per-subject runs; unequal sample counts have an explicit weighting contract | Inspect which conditions are confused and how prediction varies across subjects/strata |
| 21. Explicit train/test partitions | Expose existing fixed-partition support for a prespecified independent test set | Only designated training rows affect preprocessing, tuning and fitting; training-only rows retain no held-out predictions; incompatible partition/tuning designs fail clearly | Evaluate a model on an explicitly reserved test set without inventing CV groups |
| 22. Fold-local nuisance adjustment | Assess prediction after accounting for specified motion, age, site or other covariates | Known confound-only fixtures lose predictable nuisance signal while independent genuine signal remains; held-out perturbations cannot change fitted adjustments; each supported inference path matches its stated null | Ask a clearly defined predictive question after specified covariate adjustment |
| 23. Fold-local feature selection | Retain informative features within larger patterns using a validated selection procedure | Selection uses training rows exclusively; selected voxel/node IDs are audited; feature-count tuning remains nested; permutations rerun selection | Reduce irrelevant features while retaining an inspectable mapping to the original neural domain |
| 24. Fold-local PCA / principal-component regression | Provide controlled dimensionality reduction and compose it with regression for correlated, high-dimensional patterns | Known projection/reconstruction and composed-regression predictions pass; PCA and regression preprocessing use training data only; component-count/penalty selection is nested; every permutation refits the pipeline | Predict through an audited lower-dimensional representation without needing a separate PCA implementation per estimator |
| 25. Haufe/activation-pattern outputs | Describe signal represented by a decoder using an explicit transformation beyond raw predictive weights | Synthetic correlated-feature examples recover expected patterns; scaling, covariance estimation and response conventions are specified; descriptive outputs stay separate from coefficient significance and held-out scores | Inspect descriptive neural signal patterns under a documented interpretation contract |
| 26. Paired model comparisons | Compare estimators or encoding designs using identical held-out observations | Paired score differences match independent calculations; sample/fold alignment is checked; any inferential extension has a separately specified null and null-calibration tests | Compare predictive performance on the same evaluation data |
| 27. Saved-model application | Apply a frozen fitted model to genuinely new data | Saved/reloaded predictions agree; preprocessing and feature order are preserved; incompatible feature domains, labels, dimensions and model versions fail clearly; applying the model never refits it | Reuse a validated model on new samples with a reproducible feature/preprocessing contract |

#### Starting point and boundaries

- Item 20 builds on existing internal confusion matrices/per-class recall and
  fold stratum metadata. These library capabilities are not yet dedicated CLI
  exports. Pooled sample scores, subject-level summaries and population
  inference must remain distinct; item 20 adds reporting, not a group test.
- Item 21 exposes `THD_foldplan_fixed()`, already implemented and unit tested
  in `thd_crossval.[ch]`. It must preserve the existing distinction between
  training-only rows and scored held-out rows.
- Item 22 must first define which predictors and/or responses are adjusted,
  which covariates are available at application time, and what scientific
  quantity is being predicted. Whole-dataset residualization is not a
  substitute for fitting the adjustment inside every applicable training
  split. Existing label/row permutations must not be assumed valid for a new
  nuisance-adjusted null without justification and validation.
- Items 23 and 24 are separate slices. Define their order relative to scaling
  and nuisance adjustment, training-only selection/projection, nested tuning,
  permutation refits and feature/component audit schemas before coding.
- Item 25 extends interpretation outputs; it must not relabel raw weights as
  activation patterns or imply coefficient/activation significance. Specify
  the applicable estimators and descriptive full-data versus fold-local
  covariance contract explicitly.
- Item 26 starts with aligned descriptive comparisons. Sharing folds and
  permutation slots alone does not establish a valid equal-performance null;
  inference requires a separately defined scientific contrast.
- Item 27 requires a versioned model/feature/preprocessing contract. Define
  supported estimators and spatial analyses before implementing persistence;
  a coefficient map alone is not a complete saved model.

#### Deferred scope requiring a separate proposal

Multi-target summary/family definitions, population inference across subjects,
class weighting, compact beta-series input and residual-null methods not
explicitly covered by item 22 remain deferred. Item 39 now proposes temporal
encoding/null contracts; it does not make current row-shuffle inference valid
for serially dependent samples. Estimators beyond proposals 28–40 remain
deferred. Adding those proposals
does not authorize implementation or imply support for general probability
calibration across all classifiers. Existing within-subject folds and pooled
prediction maps do not provide a dedicated population-inference workflow.

### Phase 6: proposed estimator expansion (items 28–32)

The user approved adding these model proposals on 2026-10-15. Item 20 remains
next; this phase does not reorder or start implementation. Huber regression
has its own item so its scientific scope and validation can be approved
separately from linear SVR. All implementations remain native and use the
existing AFNI sample, spatial, grouped CV and inference infrastructure.

| Item / mode | Why are we doing this? | Completion metric | Intended benefit for fMRI users |
|---|---|---|---|
| 28. Elastic-net regression / decode, encode | Combine shrinkage with sparse predictor selection for many correlated neural or stimulus features | Known sparse/correlated fixtures and solver optimality checks pass; penalty strength and L1/L2 mixing ratio are selected inside grouped inner CV; per-response CV-R2/MAE and selected original feature IDs are audited | Identify a sparse predictive subset of voxels or design variables while retaining held-out validation |
| 29. Regularized logistic regression / classify | Model binary/multiclass outcomes directly and provide class probabilities alongside decisions | Known-answer binary/multinomial fits, stable probability normalization and nested penalty selection pass; held-out log loss, Brier score, ACC/BACC and calibration assessment are available; any fitted calibration uses training data only | Evaluate class probabilities and probability-sensitive performance for categorical fMRI prediction |
| 30. PLS regression / decode, encode | Model predictor/response relationships through a small supervised latent space, including correlated multiple responses | Known latent-factor/projection and prediction fixtures pass; component fitting and count selection use training rows only; component-count tuning is nested; single/multiple-response scores and component diagnostics are audited | Predict behavioral targets or neural response patterns through lower-dimensional shared structure |
| 31. Linear SVR / decode first; encode subject to cost assessment | Compare epsilon-insensitive loss with squared-error regression, allowing errors within an explicit tolerance | Analytic solver/convergence checks and nested C/epsilon selection pass; any additional loss/penalty choices have explicit contracts; paired held-out MAE/RMSE matches independent calculations | Test a different error-tolerance model for continuous behavioral prediction, with possible voxelwise encoding after resource validation |
| 32. Huber regression / decode, encode; later priority | Reduce domination of the fit by unusually large response residuals in appropriate applications | Known robust-fit fixtures and convergence checks pass; injected-outlier tests demonstrate expected resistance while clean-data prediction is checked; regularization and Huber threshold tuning remain nested | Assess continuous prediction with a documented robust residual loss |

#### Model-specific scope and interpretation

- Keep ridge as the reference baseline. New models must support the same
  grouped generalization definitions, training-only preprocessing and audit
  identities. Evaluate performance on identical held-out observations; no
  model is assumed to improve prediction merely because it is more flexible.
- Item 28 starts with independent per-response elastic-net fits, including
  per-voxel encoding. Sparse coefficient selection describes predictive
  selection, not independent voxel activation or coefficient significance.
  A multitask penalty encouraging shared predictor selection across responses
  is a separate extension requiring an explicit scientific assumption and
  acceptance criteria; it is not silently substituted for independent fits.
- Item 29 starts with L2-regularized binary and multinomial logistic models.
  Define class coding, intercept treatment, penalty normalization and scoring
  before coding. Probability output alone does not establish calibration.
  Define multiclass Brier/log-loss conventions and any nested calibration
  protocol explicitly. Elastic-net logistic penalties are separately proposed
  in item 40 and do not expand the initial L2 implementation.
- Item 30 uses supervised PLS regression, distinct from item 24's unsupervised
  PCA. Specify the exact PLS variant, predictor/response centering and scaling,
  deflation, component/rank limits and multi-response objective. Joint encoding
  models introduce a shared response structure; their fit units and covariance
  domain must be explicit. Test across voxelwise versus pattern use cases
  rather than assuming these are interchangeable analyses.
- Item 31 is a regression estimator, distinct from the existing classification
  SVM. Specify its epsilon-insensitive objective, intercept treatment, penalty
  normalization, target units, multi-response behavior and convergence rule.
  First support decoding; assess the runtime/memory cost of separate fits over
  many encoding voxels before committing to encoding support.
- Item 32 requires a stated Huber loss/scale convention, intercept and penalty
  treatment, convergence diagnostics and response-wise behavior. It addresses
  residual outliers; it is not a replacement for nuisance adjustment or input
  quality control. Limit initial scope to documented continuous-response cases.

#### Shared acceptance gates

1. Known-answer and independent numerical checks, high-dimensional/collinear
   and constant-feature cases, and explicit failure on nonconvergence.
2. Outer held-out perturbation tests and grouped nested hyperparameter tests
   showing that preprocessing, selection and model fitting use training rows
   only. Report candidate losses, selected settings and solver diagnostics.
3. ROI, volume and sparse-surface prediction/audit validation for each supported
   mode, with deterministic results across worker counts.
4. Define the scientific null and score direction for each new model/metric.
   Supported permutations must refit preprocessing, hyperparameter selection
   and models, then restore observed audits. Probability-sensitive metrics
   require their own validated inference path before being advertised there.
5. Keep predictive scores, descriptive coefficients and activation-pattern
   interpretations distinct. New estimators do not automatically inherit
   item 18's coefficient-export contract or item 27's model-file support;
   explicitly implement and validate compatibility before exposing it.
6. Preserve existing estimator outputs and complete the relevant regression
   tests. Report measured benefits, limitations and what fMRI users can now do,
   then prompt before advancing to the next item or expanding scope.
7. Extend `thd_mvpa.[ch]` and `thd_crossval.[ch]` for coherent numerical and CV
   support rather than adding a library per estimator. Reuse extraction,
   workspaces, maps and audits. Keep independent-toolbox comparisons in the
   validation workflow, not the distributed program; nice Python comparisons
   and cap their numerical-library/core usage at five cores.

### Phase 7: proposed tuning, encoding and generalization extensions (items 33–40)

Added on 2026-10-15 after reconciling the model/capability discussion against
the existing roadmap. Existing proposals 20–32 and completed items 18–19 are
retained, not duplicated or restarted. Item 20 remains next pending user
direction. These additions authorize planning only, not implementation.

| Item / mode | Why are we doing this? | Completion metric | Intended benefit for fMRI users |
|---|---|---|---|
| 33. Generalized inner-CV tuning / all applicable modes | Avoid estimator-specific tuning machinery and make the selection objective explicit | Existing ridge/SVM defaults and audits remain unchanged; known candidate grids select the expected parameters for explicit ACC/BACC, continuous-error or probability-aware objectives; all selection uses grouped inner training/validation rows only | Choose a model using the performance criterion appropriate to the scientific question, with inspectable candidate results |
| 34. Response-specific ridge penalties / decode, encode | Allow targets or voxel responses with different signal/noise characteristics to receive different regularization | Analytic multi-response fixtures select different known penalties; held-out response perturbations cannot change selection; selected penalties and per-response losses are audited; shared-lambda results remain unchanged; runtime/memory costs are measured | Avoid forcing every response in a joint fit to use the same regularization while keeping the existing simpler option |
| 35. Banded ridge / encode first | Regularize named visual, semantic, acoustic or other feature spaces separately | Named group membership/order is validated; known block-penalty solutions and equivalent linear-kernel fits agree; grouped nested selection excludes outer tests; selected penalties and any descriptive contribution summaries are reproducible | Combine distinct stimulus feature spaces without requiring one penalty for all predictors |
| 36. Reduced-rank regression / multi-response decode, encode | Constrain the fitted coefficient matrix to exploit a specified low-dimensional response relationship | Known rank-constrained solutions/predictions agree; rank limits, centering/scaling and any ridge penalty are explicit; rank/penalty selection is nested; joint response domains and held-out per-response scores are audited | Evaluate a low-rank alternative to independent response fits without conflating it with PCA or PLS |
| 37. Kernel ridge / decode, encode; multiple-kernel extension later | Evaluate explicit nonlinear predictor relationships after establishing linear baselines | A linear kernel reproduces the matching ridge model; known nonlinear fixtures agree with an independent solver; train/test kernel centering and nested kernel parameters use training rows only; batch/memory limits and permutation refits are tested | Test a specified nonlinear relationship, or combine kernel-defined feature spaces, with the same held-out validation discipline |
| 38. Cross-condition generalization / classify, decode, encode as defined | Separate prediction within a domain from transfer across tasks, conditions or sessions | Explicit domain partitions prohibit training access to test samples; feature IDs and class/response mappings are checked; known transfer fixtures and disjointness checks pass; incompatible target domains fail clearly | Ask whether information learned in one experimental context predicts another |
| 39. Time-series encoding / encode first | Model delayed neural responses without crossing run boundaries or leaking neighboring time points | Lag alignment is verified on known impulse fixtures; delays never wrap across runs; buffered splits exclude the specified temporal neighborhood; lag/penalty selection remains nested; each supported temporal null has independent null-calibration tests | Fit continuous-stimulus encoding analyses with explicit hemodynamic timing and temporal-dependence assumptions |
| 40. Elastic-net logistic classification / classify; later priority | Add sparse class prediction after establishing an L2 probability model | The L2 limit matches item 29; sparse/correlated fixtures, convergence/optimality and nested strength/mixing selection pass; selected original feature IDs and probability-sensitive held-out metrics are audited | Explore sparse categorical predictors without interpreting selected voxels as independent activation effects |

#### Dependencies and scientific boundaries

- Item 33 generalizes candidate evaluation and objective selection, not the
  outer folds or inference family. Preserve current normalized squared-loss
  ridge tuning and SVM defaults unless users explicitly select another
  objective. Define aggregation, undefined-score handling and tie rules per
  objective; no selection or fitted calibration may use outer-test outcomes.
- Item 34 must distinguish per-response penalties within a multi-output fit
  from existing voxel encoding, which already creates singleton response
  units. Its main added value is independent penalties in shared multi-response
  fits, including pattern encoding. Do not silently replace a pattern-level
  tuning objective or alter target correction families.
- Item 35 defines predictor feature groups, not sample/CV groups. Start with
  explicit nonoverlapping named design-column groups; specify scaling, penalty
  normalization and whether penalties are shared or response-specific.
  Algebraically equivalent linear multiple-kernel computation may support
  banded ridge without implying support for arbitrary nonlinear kernels.
  Descriptive feature-space contribution summaries are not unique causal
  effects when groups are correlated; inferential comparisons belong to item 26.
- Item 36 is reduced-rank regression, not item 30's PLS algorithm or item 24's
  PCA-plus-regression pipeline. Specify the fitted rank constraint and shared
  response domain before implementation. Kernel models likewise do not inherit
  voxel-space coefficient exports or saved-model support automatically.
- Item 37 starts with a small explicit kernel catalog and measured CPU/memory
  cost. Define the training kernel and new-sample cross-kernel contract, and
  nest any bandwidth/degree selection. Multiple-kernel weights need their own
  validated selection and interpretation contract. No nonlinear model is
  presumed to improve generalization solely by being more flexible.
- Item 38 builds on item 21's explicit partitions and item 27's frozen-model
  application where appropriate; it does not require persistence for every
  transfer analysis. Define permitted overlap of subjects, items and sessions,
  and distinguish within-subject transfer from held-out-subject prediction.
- Item 39 requires its own scientific sampling and null contracts. Specify
  sample timing, censoring/missing rows, lag construction and temporal buffers;
  fit learned transforms inside each training split. Neither unrestricted nor
  within-run row shuffling automatically preserves temporal dependence.
- Item 40 follows item 29 and generalized grouped tuning. Sparse feature
  selection is predictive, not activation inference; probabilities still need
  calibration assessment. It is distinct from item 28's continuous-response
  elastic-net regression and its independent per-response penalty contract.

#### Reconciled proposals and priorities

- Existing items 22–27 already cover nuisance adjustment, feature selection,
  PCA, activation-pattern maps, paired comparisons and saved-model application.
  Principal-component regression extends item 24 rather than adding another
  PCA proposal. Item 29 already includes probability scoring/calibration, and
  item 30 already covers PLS; neither is duplicated here.
- Coefficient maps and linear classification SVM are complete (items 18–19),
  so recommendations to implement them are superseded by the current record.
- After the authorized audit/partition work, generalized tuning and logistic
  regression are high-priority classification extensions; response-specific
  and banded ridge are high-priority encoding extensions. Elastic net and
  fold-local transforms broaden continuous prediction. Reduced-rank and
  nonlinear kernel models remain later priorities.
- Random forests, boosting and neural networks are not proposed for the
  initial native estimator catalog. Revisit only for a concrete use case with
  justified validation, runtime and maintenance costs.

## Initial Command Surface

Classification:

```bash
3dMVPA                                      \
    -mode classify                         \
    -dataTable @samples.tsv                 \
    -target_column Condition                \
    -cv_column Run                          \
    -mask brain_mask+tlrc                   \
    -searchlight 'SPHERE(5)'                \
    -method centroid                        \
    -metric balanced_accuracy               \
    -prefix face_scene
```

The same analysis from AFNI regression output may use a subject-level table
with `Xmat` and whole stats buckets:

```bash
3dMVPA                                      \
    -mode classify                         \
    -dataTable @subjects_with_xmat.tsv      \
    -xmat_stims face,house                  \
    -target_column Condition                \
    -cv_column Run                          \
    -mask brain_mask+tlrc                   \
    -searchlight 'SPHERE(5)'                \
    -method centroid                        \
    -metric balanced_accuracy               \
    -prefix face_scene
```

These command forms differ only in how the canonical sample table is created.
After X-matrix expansion they use identical validation, folds, features,
estimators, and outputs.

Continuous decoding:

```bash
3dMVPA                                      \
    -mode decode                           \
    -dataTable @samples.tsv                 \
    -target_columns Rating,Arousal          \
    -cv_column Run                          \
    -mask atlas+tlrc                        \
    -method ridge                           \
    -lambda 1                               \
    -metric pearson                         \
    -prefix decoded_behavior
```

Voxelwise encoding:

```bash
3dMVPA                                      \
    -mode encode                           \
    -dataTable @samples.tsv                 \
    -design_columns Visual,Semantic,Motion  \
    -cv_column Run                          \
    -mask brain_mask+tlrc                   \
    -score voxel                            \
    -method ridge                           \
    -lambda 1                               \
    -metric cross_validated_r2              \
    -prefix encoding_model
```

The precise option names can evolve during implementation, but the concepts
should remain stable: direction, target/design roles, generalization unit,
estimator, score, exchangeability, and output.

## Design Principles

1. **Prediction is evaluated only on held-out data.** Descriptive full-data
   refits are separate outputs.
2. **Every learned operation is fold-local.** There are no convenient leakage
   exceptions.
3. **The public design is inspectable.** Users can print resolved samples,
   classes, folds, blocks, and exclusions before expensive computation.
4. **Inference repeats the real analysis.** Permutations do not bypass
   preprocessing, tuning, or model fitting.
5. **Spatial code and numerical code remain separate.** `thd_patterns` knows
   voxels; `thd_mvpa` knows arrays and folds.
6. **Manual and X-matrix input are both first-class.** A direct sample table,
   an X-matrix-expanded table, and later compact forms all converge on one
   canonical table and one `THD_sampleset`; none creates a separate analysis
   path.
7. **3dRSA and 3dMVPA are siblings, not modes of one oversized program.**
8. **Shared scientific algorithms live in testable `thd_` libraries.**
9. **Simple validated estimators come before a large method catalog.**
10. **Thread-count reproducibility is part of the API contract.**

The intended outcome is a coherent AFNI pattern-analysis stack: `3dRSA` for
representational geometry, `3dMVPA` for held-out prediction, and shared
libraries for sample organization, spatial extraction, cross-validation,
inference, and numerical modeling.

## Toolbox Comparison and Remaining Gaps (2026-10-15)

This is a capability and workflow assessment, not a numerical or speed
benchmark. It compares the implemented CLI/library surface in this checkout
against the primary documentation linked below. Documentation reviewed was
labeled PyMVPA 2.6.5.dev1 and BrainIAK 0.12; these labels are not assertions
about their newest releases or the versions installed locally. Proposed
roadmap items are not counted as implemented capabilities. Names and links
here are comparison references, not runtime dependencies or program branding.

### Overall position

`3dMVPA` is now a useful, focused native AFNI prediction program: centroid,
ridge, diagonal/shrinkage LDA and linear SVM classification; ridge continuous
decoding and voxel/pattern encoding; ROI, volume and geodesic surface analyses;
grouped nested ridge/SVM tuning; blocked empirical spatial inference; held-out
prediction/fold audits; and separate descriptive ridge coefficient exports.
It is not yet a replacement for a general-purpose MVPA research framework.

PyMVPA is the closer comparison for flexible classification/regression
pipelines. Its examples include surface searchlights, nested CV, estimator
adapters, kernels and sensitivity analyses. Our present estimator and pipeline
catalog is narrower. [PyMVPA analysis examples](https://www.pymvpa.org/examples.html).

BrainIAK is a broader computational-neuroimaging toolkit, not simply a rival
classifier catalog: it includes functional alignment, connectivity analysis,
event segmentation, reconstruction and Bayesian representational analyses.
Missing those methods is a scope difference, not evidence that our implemented
ridge or classifier is numerically inferior. [BrainIAK module index](https://brainiak.org/docs/py-modindex.html).

### Where we currently fall short

| Area | Current limitation in 3dMVPA | Comparison and planned response |
|---|---|---|
| Estimator breadth | No logistic/probability model, elastic net, PLS, SVR, Huber or nonlinear kernels yet | PyMVPA exposes additional learners and external estimator adapters, including sparse multinomial logistic regression. Existing proposals 28–32, 36–37 and 40 cover selected expansions, not blanket estimator compatibility. [Examples](https://www.pymvpa.org/examples.html), [sparse logistic model](https://www.pymvpa.org/generated/mvpa2.clfs.smlr.SMLR.html). |
| Pipeline composition | No CLI fold-local feature selection, PCA/PCR or nuisance adjustment; arbitrary user-defined measures cannot simply be plugged into the command | PyMVPA composes mappers and meta-classifiers within CV, and documents training-only feature selection. Items 22–24 and 33 are higher-value gaps than adding many nearly equivalent classifiers. [Meta-classifiers](https://www.pymvpa.org/tutorial_meta_classifiers.html), [training-only selection](https://www.pymvpa.org/tutorial_sensitivity.html). |
| Training/application workflow | CLI fixed independent partitions, frozen model reload/application and cross-condition transfer are still proposed | PyMVPA has reusable fitted learners/mappers and HDF5 object storage. Our coefficient exports are not complete reusable model files. Items 21, 27 and 38 need versioned feature/preprocessing contracts. [Object storage](https://www.pymvpa.org/generated/mvpa2.base.hdf5.html), [mapper composition](https://www.pymvpa.org/tutorial_meta_classifiers.html). |
| Functional alignment | No fold-local SRM or hyperalignment stage inside 3dMVPA | PyMVPA documents ROI/searchlight hyperalignment; BrainIAK offers SRM variants. The separate 3dSRM work does not establish leakage-safe integration with this classifier pipeline. Any integration needs its own training/application and held-out-subject contract; do not silently add it to current scope. [Hyperalignment](https://www.pymvpa.org/examples/hyperalignment.html), [functional alignment](https://brainiak.org/docs/brainiak.funcalign.html). |
| Specialized encoding/reconstruction | Ridge prediction of neural responses is not an inverted encoding model or circular stimulus reconstruction | BrainIAK documents 1D/2D basis-channel reconstruction. Items 34–37 improve predictive encoding but do not supply that distinct analysis; circular reconstruction would need a separate scientific proposal. [Inverted encoding models](https://brainiak.org/docs/brainiak.reconstruct.html). |
| Distributed/custom searchlights | OpenMP spatial workers on one machine; no MPI engine or user-supplied searchlight callback | BrainIAK distributes data across MPI ranks and runs user-defined voxel/block functions. This is a real flexibility/scaling gap; adding it should depend on measured workload requirements, not an assumption that native C is faster. [Distributed searchlights](https://brainiak.org/docs/brainiak.searchlight.html). |
| Reporting and interpretation | Confusion/recall are internal rather than dedicated exports; subject summaries, calibrated probabilities, Haufe patterns and inferential model comparisons remain incomplete | Items 20, 25–26 and 29 address these gaps. PyMVPA already demonstrates confusion/sensitivity outputs; raw coefficients are not activation maps. This is a usability/scientific-interpretation priority, not just a model-count issue. [Sensitivity and confusion outputs](https://www.pymvpa.org/tutorial_sensitivity.html). |
| Temporal designs and advanced nulls | Current continuous nulls require exchangeable rows; no general temporal, nuisance-residual or population-inference workflow | Existing row blocking must not be advertised as an autocorrelation model. Items 22 and 39 require explicit null validation; subject summaries alone do not provide population inference. This is our limitation, not a claim that either comparison toolbox automatically solves every design. |
| Validation and user ecosystem | Focused local known-answer/regression tests, but no comprehensive matched cross-toolbox benchmark or comparable body of application examples established here | Treat numerical agreement, null calibration, resource use and ease of use as separate validation goals. Documentation and reproducible real-data examples need to grow alongside features. Do not infer another project's maintenance status from documentation age. |

BrainIAK's FCMA, event segmentation, Bayesian RSA and other specialized models
are additional out-of-scope differences. Likewise, the existence of `3dRSA`
does not imply equivalence with Bayesian RSA. We should not reproduce an entire
toolkit inside one command. [BrainIAK module catalog](https://brainiak.org/docs/py-modindex.html).

### Where we have a concrete advantage for AFNI users

These are current integration strengths and workflow judgments, not claims of
unique algorithms or demonstrated improvements in scientific accuracy:

- **Direct AFNI regression-to-prediction workflow.** Manual sample tables and
  X-matrix expansion converge on the same sample/fold validation path, including
  selected coefficient sub-bricks. This can remove custom conversion/glue code
  for users already working with AFNI regression outputs.
- **A packaged end-to-end command for the supported analyses.** Extraction,
  fold-local fitting/tuning, synchronized blocked permutations, BH/max-statistic
  correction and AFNI statistical-map metadata are connected in one interface.
  This is an advantage in setup consistency, not proof that other toolboxes
  cannot implement the same analysis.
- **Native volume and sparse-surface output contracts.** Preserved mesh node
  IDs, domain validation and standard AFNI outputs make the existing surface
  path directly usable in an AFNI/SUMA workflow. PyMVPA also supports surfaces;
  surface analysis itself is not an exclusive capability. BrainIAK's documented
  searchlight interface is volumetric; custom extensions are a different claim.
  [PyMVPA surface examples](https://www.pymvpa.org/examples.html),
  [BrainIAK searchlight interface](https://brainiak.org/docs/brainiak.searchlight.html).
- **Consistent auditability across supported modes.** Sample identities,
  original feature/node IDs, held-out predictions, selected penalties, fold
  diagnostics and descriptive-refit labeling have explicit conventions and
  regression tests, including worker-count agreement and restoration after nulls.
  This is a strength of our contract, not evidence that other tools are
  irreproducible or less statistically careful.
- **No runtime Python scientific stack for these native analyses.** Users with
  the AFNI build do not need to assemble Python estimator pipelines. AFNI still
  has build/runtime dependencies, and external toolboxes remain valuable as
  independent validation references.

In particular, nested CV, training-only transformations, permutations and
surface searchlights are not new advantages invented by this program. PyMVPA
documents these capabilities. Our advantage is their constrained AFNI-native
packaging for the analyses we support. [Nested/surface examples](https://www.pymvpa.org/examples.html),
[permutation testing](https://www.pymvpa.org/examples/permutation_test.html),
[training-only feature selection](https://www.pymvpa.org/tutorial_sensitivity.html).

### What is not established, and how to test it

Do not claim better prediction accuracy, stronger statistical validity, lower
memory use or faster execution than either toolbox without matched benchmarks.
Compiled C is not sufficient evidence: other toolboxes also use optimized
numerical code and may offer different parallel/distributed strategies.

For a reproducible comparison, record exact installed versions and:

1. Match samples, voxel/node order, train/test folds, class encoding, scaling,
   intercept/prior conventions, covariance denominators and penalty objectives.
   Match squared-hinge versus hinge SVM loss and intercept regularization;
   similarly named estimators need not fit identical models.
2. Compare fixed-parameter decisions/predictions before comparing tuned models.
   Then match inner folds, candidate grids, objective aggregation and tie rules.
   Align CV-R2 baseline definitions rather than comparing score labels alone.
3. Supply the same permutation index arrays and correction families. Compare
   empirical counts and corrected probabilities, plus null-calibration fixtures;
   a common seed does not ensure common permutations across implementations.
4. Benchmark ROI, volume and sparse-surface workloads separately. Record setup
   and I/O time, fitting/null time, peak resident memory, thread/process counts,
   sample/feature sizes and null counts. A single-node comparison is not an
   MPI scalability comparison.
5. Include real AFNI-derived sample tables and reproducible synthetic nulls,
   not only easy separable examples. Keep reference-toolbox names/dependencies
   out of program code and runtime requirements. Nice Python runs and cap their
   numerical-library/core use at five cores, as requested.

The most valuable near-term gap closures remain item 20's reporting and item
21's explicit independent partitions, followed by training-only transforms,
generalized tuning/logistic classification, model application and the encoding
extensions already proposed. This comparison does not reorder the active
roadmap or authorize additional implementation.
