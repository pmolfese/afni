# 3dMVPA Design and Implementation Plan

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

The `SRMthings` branch already provides:

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

## New `thd_sampleset.[ch]`

The main missing input abstraction is a general observation map. It should
preserve every sample and connect metadata rows to datasets and bricks.

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

A useful general extraction primitive would resemble:

```c
int THD_roi_samples(
    THD_3dim_dataset **dset,
    const int *dset_index,
    const int *brick_index,
    int nsample,
    const intvec *vox,
    float *out
);
```

`out` is `[sample][voxel]` and caller allocated, allowing one workspace per
OpenMP worker.

Feature filtering belongs in the fold-aware MVPA layer when it depends on the
data. Simple spatial validity checks may occur during extraction, but variance
thresholding, standardization, nuisance regression, and learned selection must
be fitted on training data only.

## New `thd_mvpa.[ch]`

`thd_mvpa` will be a dataset-agnostic numerical library. It will consume arrays
and fold plans, not AFNI datasets or command-line options.

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

### First general solver: ridge

Ridge regression should be the first shared fitted estimator because one
implementation supports all three modes:

- classification through one-hot targets and argmax prediction;
- scalar or multivariate continuous decoding;
- multi-output encoding.

Implement both numerical forms:

```text
primal: (X'X + lambda I)^-1 X'Y       when predictors <= samples
dual:   X'(XX' + lambda I)^-1 Y       when predictors > samples
```

The dual formulation is essential for neuroimaging decoding, where voxels
often greatly outnumber samples. Encoding commonly has a small external design
and many voxel responses, allowing one factorization to predict every voxel.

### Additional classifiers

After centroid and ridge:

1. diagonal LDA;
2. shrinkage LDA;
3. native linear L2 SVM with one-vs-rest multiclass composition;
4. multinomial ridge logistic regression if calibrated scores are needed.

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

This distinction prevents encoding from being forced into a decoding-shaped
output. Optional full-data refit coefficients may be written as one map per
design column, clearly labeled as descriptive refits rather than held-out
performance.

## Permutation Inference

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

The user must be able to define exchangeability through subject, run, fold, or
an explicit block column. Invalid permutations should fail rather than silently
crossing nonexchangeable boundaries.

Score-specific null handling:

- classification: report accuracy and accuracy-minus-chance; infer against
  label permutations rather than a symmetric theoretical distribution;
- correlation: aggregate Fisher-z values and center the null at zero;
- R2: use its permutation distribution rather than assuming symmetry around
  zero;
- paired model comparisons: drive both models with identical folds and
  permutation slots.

Reuse `THD_perm_result`, `THD_perm_emp_pvalue()`, `THD_perm_signed_z()`,
`THD_max_accum()`, and `THD_bh_fdr_masked()` for final map inference.

## Spatial Analysis and Output

### Searchlight

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

The prediction table should retain original table row, sample label, subject,
run, fold, observed target, prediction, and decision score. This is part of the
scientific result, not merely debugging output.

## Aggressive 3dRSA Library Extraction

`3dRSA.c` is now large enough that waiting for duplicate code to appear in
`3dMVPA` would be counterproductive. Refactor behavior-preserving scientific
code before or alongside the new program.

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

### Library unit tests

- direct sample tables and X-matrix-expanded sample tables producing equivalent
  `THD_sampleset` organization;
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

- synthetic separable classes with a known spatial signal location;
- continuous target encoded in a known searchlight;
- encoding design with known responsive voxels;
- ROI result matching a direct C-level calculation;
- chance data producing sane permutation p-values;
- max-FWE behavior under a synchronized null;
- prediction-table rows matching original metadata;
- clear failures for missing classes, impossible folds, grid mismatch, and
  non-finite input.

### 3dRSA regression protection

Library extraction must be checked against existing 3dRSA oracle and synthetic
tests. Extraction commits should not change numerical behavior except where a
separately documented bug is intentionally corrected.

## Implementation Sequence

### Phase 0: behavior-preserving shared-library work

1. Bring `thd_xmat.[ch]`, its `thd_datatable` parsing improvements, and its
   tests from `SRMthings` into the shared-library baseline.
2. Add `thd_sampleset.[ch]` and test both direct and X-matrix-expanded tables.
3. Generalize selected-sample extraction in `thd_patterns`.
4. Add `thd_rsa.[ch]` and move existing 3dRSA scientific helpers in coherent
   groups.
5. Keep 3dRSA output numerically identical and validate against its oracle.

### Phase 1: MVPA numerical core

1. Add fold-plan and workspace objects.
2. Add leakage-safe centering/scaling.
3. Add prediction metrics.
4. Add nearest-centroid and correlation classifiers.
5. Add primal and dual multi-output ridge.
6. Test all algorithms with small matrices having known solutions.

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

### Phase 3: inference and fitted estimators

1. Add blocked label permutations and PVAL/Z/Q/FWE outputs.
2. Add ridge classification and continuous decoding.
3. Add nested grouped lambda selection.
4. Add diagonal and shrinkage LDA.
5. Add native linear SVM if still scientifically useful.

### Phase 4: encoding

1. Add voxelwise multi-output ridge encoding.
2. Write per-regressor descriptive refit coefficient maps when requested.
3. Add ROI/searchlight pattern prediction.
4. Add encoding-appropriate prediction tables and inference.

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
    -atlas atlas+tlrc                       \
    -method ridge                           \
    -lambda auto                            \
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
    -lambda auto                            \
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
