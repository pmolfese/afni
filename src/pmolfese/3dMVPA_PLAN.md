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
18. **Descriptive Coefficient Maps (Next)** -- expose requested full-data ridge
   coefficients separately from held-out scores, with feature/response IDs and
   explicit scaling and interpretation contracts.
19. **Native Linear SVM (Queued)** -- assess and add a validated native linear
   margin-based classifier after the current estimators and coefficient outputs.

## Implementation Status

This checked-in copy is the living implementation record for the `somemvpa`
branch. Update it with each completed slice so that design decisions and code do
not drift apart.

Status as of 2026-10-14:

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

After centroid, ridge, and LDA:

1. native linear L2 SVM with one-vs-rest multiclass composition;
2. multinomial ridge logistic regression if calibrated scores are needed.

Keep the initial SVM linear. Kernels, real-time support, SVM-light model
compatibility, and transductive modes are out of scope.

Later additions may include elastic net and sparse encoding/decoding, but only
after nested-CV behavior and high-dimensional numerical stability are covered
by tests.

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
Descriptive full-data coefficient maps remain queued.

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
7. [ ] Add native linear SVM if still scientifically useful.

### Phase 4: encoding

1. [x] Add voxelwise ridge encoding using the multi-output solver.
2. [ ] Write per-regressor descriptive refit coefficient maps when requested.
3. [x] Add ROI/volume/surface searchlight pattern prediction.
4. [x] Add encoding prediction and fold tables with response voxel/node IDs.
5. [x] Add encoding-appropriate joint design-row inference for exchangeable
   samples; temporal/residual nulls require separate future contracts.

### Phase 5: advanced features

- multi-target summaries;
- nuisance residualization;
- fold-local feature selection and PCA;
- elastic net;
- model comparison;
- interpretable activation/Haufe maps;
- saved-model application to genuinely new datasets.

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
