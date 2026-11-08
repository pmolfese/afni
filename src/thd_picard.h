#ifndef THD_PICARD_H
#define THD_PICARD_H
#include "thd_regression.h"

/* One public entry point for projection and Picard's numerical workflows.
   The numbered sections match thd_picard.c. Callers normally use opaque
   plans and the checked info/result views; they do not inspect plan fields.
   All matrix dimensions/ownership contracts are documented at each API. */

/* ========================================================================
   1. NUISANCE PROJECTION
   ======================================================================== */

/* Double-precision nuisance projection of exactly the supplied Z[n][k].
   No implicit intercept, centering, scaling, or normal equations. Relative
   SVD rank policy: s > rank_tol*smax in supplied units. Recovery returns the
   minimum-norm nuisance coefficients; redundant columns retain their mapping.
   k=0 is supported (identity projection, no nuisance coefficients).
   Immutable plans own their factors; a private workspace borrows the plan.
   Buffers are compact row-major with active q strides. Project permits exact
   input/output aliasing; solve requires separate output. No fit allocations.
   Projection makes a second pass to suppress roundoff, and treats residual
   norms <=64*DBL_EPSILON*max(n,k)*input_norm as numerical zero. */
typedef struct THD_linproj THD_linproj;
typedef struct THD_linproj_workspace THD_linproj_workspace;
THD_linproj *THD_linproj_prepare(const double *z, int n, int k, double rank_tol,
                                 char *err, int nerr);
int THD_linproj_rank(const THD_linproj *plan);
void THD_linproj_free(THD_linproj *plan);
THD_linproj_workspace *THD_linproj_workspace_alloc(const THD_linproj *plan,
                                                   int max_response, char *err,
                                                   int nerr);
void THD_linproj_workspace_free(THD_linproj_workspace *w);
int THD_linproj_project(THD_linproj_workspace *w, const double *y, int q,
                        double *projected, char *err, int nerr);
int THD_linproj_solve(THD_linproj_workspace *w, const double *y, int q,
                      double *gamma, char *err, int nerr);

/* ========================================================================
   2. FIXED FITS AND REPEATED-CONDITION TUNING
   ======================================================================== */

/* Fixed unit-marginal covariance. model: 0=white, 1=AR(1), 2=ARMA(1,1).
   epsilon[t]=a*epsilon[t-1]+u[t]+b*u[t-1], finite |a|,|b|<1.
   reset_gaps=0 preserves correlation at original TR lags; 1 restarts at gaps.
   Runs are independent. No estimation or interpolation is performed. */
typedef struct {
   int model, reset_gaps;
   double a, b;
} THD_picard_covariance;

/* Native grid REML for a supplied pilot mean model, independent responses.
   Shared parameters across independent run blocks, one pair per response.
   matrix[n][d], y[n][q] are row-major; row/run IDs describe estimation rows
   only. Redundant pilot columns are reduced to an orthonormal fixed span.
   Grid is [-0.9,0.9] by 0.1 (white: one point; AR1: b=0). No fallback.
   status: 0=interior, 1=grid boundary, 3=flat grid, 4=tied optimum.
   Zero variance fails. Ties within 1e-10*(1+|score|) favor smaller |a|+|b|.
   Return 0 on invalid/saturated design or numerical failure. Results are
   caller-owned q-element arrays; all estimates are conditional on the pilot.
   score omits candidate-independent likelihood constants; sigma2 is marginal.
   Memory O(n*n+n*(d+q)); factors are discarded after each grid cell. */
typedef struct {
   double a, b, sigma2, score, separation;
   int rank, status;
} THD_picard_reml_result;
int THD_picard_reml(const double *matrix, int n, int d,
                    const int *rows, const int *runs, const double *y, int q,
                    int model, int reset_gaps, double rank_tol,
                    THD_picard_reml_result *result, char *err, int nerr);

#define THD_PICARD_NUISANCE 0
#define THD_PICARD_PENALIZED 1
/* Fixed-design workflow. matrix[n][ncol] and roles[ncol] are copied.
   Every column has one explicit role, in original column order. No columns,
   intercepts or scales are added. Z projection precedes fractional ridge;
   original-unit nuisance coefficients are recovered from Y-X*B afterward.
   row_ids[n] (or identity) and run_ids[n] (or zero) are copied; rows must be
   strictly increasing within [0,nfull), run IDs nonnegative/nondecreasing.
   No file, AFNI dataset, cross-validation or inference dependencies. */
typedef struct THD_picard_plan THD_picard_plan;
typedef struct THD_picard_workspace THD_picard_workspace;
typedef struct {
   int n, ncol, p, k, nfull, xrank, zrank;
   const int *penalized_columns, *nuisance_columns, *row_ids, *run_ids;
   double rank_tol;
   int banded, ngroup, ncandidate, support_refit;
   double ar1;
   int covariance_enabled;
   THD_picard_covariance covariance;
   int elastic; /* 1: elastic net. Both flags zero: fractional ridge. */
} THD_picard_design_info;
typedef struct {
   int n, p, k, q;
   const double *beta, *gamma, *fitted, *residual;
   const THD_fraction_diagnostic *diagnostic; /* Only fractional ridge/OLS. */
   const THD_elastic_diagnostic
       *elastic_diagnostic; /* Selection-fit diagnostics. */
   const THD_banded_diagnostic *banded_diagnostic; /* Only banded fits. */
   const int *support_size, *support_rank; /* NULL unless support refit ran. */
} THD_picard_result;
/* Additional fitting policies. At most one of elastic/banded may be set.
   ar1 is fixed externally, finite in (-1,1); zero is identity. Whitening
   restarts at run changes and nonconsecutive original row IDs. Support-only
   OLS refits are permitted only for elastic/lasso with rho>0. Configuration
   arrays are copied; no data-driven noise estimation is performed. */
typedef struct {
   const THD_elastic_config *elastic;
   const THD_banded_config *banded;
   double ar1;
   int support_refit;
   /* Optional copied covariance; mutually exclusive with nonzero legacy ar1.
      NULL preserves legacy behavior. Plans own dense Cholesky factors. */
   const THD_picard_covariance *covariance;
} THD_picard_config;
THD_picard_plan *THD_picard_plan_create_config(
    const double *matrix, int n, int ncol, const int *roles, const int *row_ids,
    const int *run_ids, int nfull, double rank_tol,
    const THD_picard_config *config, char *err, int nerr);
THD_picard_plan *THD_picard_plan_create(const double *matrix, int n, int ncol,
                                        const int *roles, const int *row_ids,
                                        const int *run_ids, int nfull,
                                        double rank_tol, char *err, int nerr);
/* Same role/identity contract, with explicit positive lambda fits. Project X
   and Y against Z first, then optionally RMS-scale projected X; do NOT center
   implicitly. Each CV training plan learns its own scales. Export beta in
   original column units and recover gamma from the original Y-X*beta.
   The config is copied; no response-dependent preprocessing occurs here. */
THD_picard_plan *THD_picard_plan_create_elastic(
    const double *matrix, int n, int ncol, const int *roles, const int *row_ids,
    const int *run_ids, int nfull, double rank_tol,
    const THD_elastic_config *config, char *err, int nerr);
const THD_picard_design_info *THD_picard_plan_info(const THD_picard_plan *p);
void THD_picard_plan_free(THD_picard_plan *p);
/* Immutable plan must outlive its worker workspaces. All returned results
   borrow worker buffers and expire on any fit attempt or workspace free.
   Work is bounded by max_response. Ordinary fits allocate no scratch during
   fitting. Support refits factor one response-specific active design at a
   time using bounded temporary SVD allocations. */
THD_picard_workspace *THD_picard_workspace_alloc(const THD_picard_plan *p,
                                                 int max_response, char *err,
                                                 int nerr);
void THD_picard_workspace_free(THD_picard_workspace *w);
/* parameters[q] are fractions for the original constructor, positive lambdas
   for the elastic constructor, exact zero-based candidate indices for banded
   fits. This entry point always starts cold. */
int THD_picard_fit_batch(THD_picard_workspace *w, const double *y, int q,
                         const double *parameters, char *err, int nerr);
const THD_picard_result *THD_picard_get_result(const THD_picard_workspace *w);
/* HRF comparison after a successful unpenalized fit. Write q losses into
   caller-owned scores: ||L^-1 residual||^2/n for explicit covariance, ordinary
   residual MSE otherwise (including the historical -prewhite_ar1 policy).
   Reuses private scratch without altering the borrowed fit/result arrays.
   Compare only candidates sharing covariance and projected rank. This is a
   conditional fit criterion, not a REML score or held-out performance. */
int THD_picard_selection_mse(THD_picard_workspace *w, double *scores,
                             char *err, int nerr);

/* Repeated-condition tuning is a beta-replication score, not prediction of a
   held-out time series. One trial must describe each penalized column, in the
   plan's penalized-column order. IDs are nonnegative; trial IDs are unique.
   Every column (including nuisance) must have support in just one run.
   Trial columns must have nonzero support in their declared run.

   Each fold omits one run, fits all remaining trials jointly at each fraction,
   and compares each training beta with EVERY matching-condition held-out OLS
   beta. Pool squared errors over ordered pairs, then divide by pair count.
   Conditions without cross-run repeats do not contribute. Each run must have
   at least one usable pair. Ties choose the largest fraction (least shrinkage).
   Scores use supplied coefficient units, with no session normalization. They
   are selection losses, NOT independent performance estimates.

   Finally refit the full supplied design at each response's winning fraction.
   Fraction norms (and therefore lambdas) are local to each training design;
   the full refit has its own norm/lambda. No post-fit beta rescaling is done.
   Plans copy metadata and grids; the original fixed plan must outlive them.
   For an elastic full plan, the grid holds strictly decreasing positive lambdas
   at the plan's fixed rho. Paths warm-start within each fold only. OLS targets
   and baseline remain unregularized; exact ties choose the smallest lambda.
   Final refits start cold. The historical `fractions` info field then contains
   lambdas, not coefficient-norm fractions. For banded plans it contains
   candidate-table indices; each fold retains original group IDs even if a
   group has no remaining columns. Ties choose smaller mean group lambda,
   then earlier grid entry. Optional fixed AR(1) restarts within each training
   fold; support refits select only from that fold. Scores use original beta
   units, while the OLS targets use the same fixed whitening policy. */
typedef struct {
   int column, trial, condition, run;
} THD_picard_trial;
typedef struct THD_picard_repeat_plan THD_picard_repeat_plan;
typedef struct THD_picard_repeat_workspace THD_picard_repeat_workspace;
typedef struct {
   int nfold, ncandidate, nscored;
   double npair;
   const double *fractions, *fold_pairs;
   const int *fold_runs;
   const THD_picard_trial *trials;
} THD_picard_repeat_info;
typedef struct {
   int q;
   const THD_picard_result *fit;
   const int *candidate;
   const double *loss, *reference_loss;
   /* Candidate-major pooled MSE [ncandidate][q]. Fold MSE is
      [nfold][ncandidate+1][q], with OLS in the last candidate slot. */
   const double *candidate_loss, *fold_loss;
} THD_picard_repeat_result;
THD_picard_repeat_plan *THD_picard_repeat_plan_create(
    const THD_picard_plan *full, const THD_picard_trial *trials,
    const double *fractions, int count, char *err, int nerr);
const THD_picard_repeat_info *
THD_picard_repeat_plan_info(const THD_picard_repeat_plan *p);
void THD_picard_repeat_plan_free(THD_picard_repeat_plan *p);
THD_picard_repeat_workspace *
THD_picard_repeat_workspace_alloc(const THD_picard_repeat_plan *p,
                                  int max_response, char *err, int nerr);
void THD_picard_repeat_workspace_free(THD_picard_repeat_workspace *w);
/* Results expire on every attempt or free; support refits allocate as above. */
int THD_picard_repeat_fit(THD_picard_repeat_workspace *w, const double *y,
                          int q, char *err, int nerr);
const THD_picard_repeat_result *
THD_picard_repeat_get_result(const THD_picard_repeat_workspace *w);

/* ========================================================================
   3. TIME-SERIES PREDICTION AND GLOBAL SELECTION
   ======================================================================== */

struct THD_foldplan;

#define THD_PICARD_PREDICT_RAW 0
#define THD_PICARD_PREDICT_PROJECT 1

/* Shared-feature time-series CV, distinct from repeated-trial beta tuning.
   full must outlive the returned plan/workspaces. Folds and grid are copied.
   Every retained row must be tested exactly once. Training factors, nuisance
   fits and fractional reference norms use training rows only. All penalized
   columns must have nonzero training support in every fold.

   RAW predicts X_test beta + Z_test gamma using training coefficients only;
   test nuisance rows must lie in the row space of Z_train. PROJECT instead
   scores P(Z_test)Y_test against P(Z_test)X_test beta. This is conditional
   signal prediction, explicitly using held-out Y to remove a SUPPLIED nuisance
   span, never to learn beta or lambda. Neither mode adds/rescales columns.

   Baseline: training response mean, passed through the same test projection.
   Pool SSE over test rows; MSE=SSE/n and R2=1-SSE/baseline_SSE. Negative R2 is
   retained. Exactly zero baseline SSE gives NaN R2/status=1; MSE stays valid.
   These are selection scores, not independent outer-test performance.

   Minimum-loss ties prefer the largest fraction. one_se selects the smallest
   fraction with pooled MSE <= best MSE + SE(best). For weights w_f=ntest_f/n,
   SE^2=sum(w_f*(MSE_f-mean)^2)*sum(w_f^2)/(1-sum(w_f^2)). This weighted-fold
   heuristic is a selection rule, not a confidence interval (folds overlap).
   Final fits use the separately prepared full design.
   For an elastic full plan, the grid contains positive decreasing lambdas at
   fixed rho (including in the historical `fractions` info field). Each fold
   prepares its own RMS scales and warm-starts only within its response batch.
   Minimum-loss ties/all-undefined global data choose the smallest lambda;
   one_se chooses the largest eligible lambda. Final refits start cold. */
typedef struct THD_picard_predict_plan THD_picard_predict_plan;
typedef struct THD_picard_predict_workspace THD_picard_predict_workspace;
typedef struct {
   int nfold, ncandidate, nuisance_mode, one_se;
   const double *fractions;
} THD_picard_predict_info;
typedef struct {
   int q;
   const THD_picard_result *fit; /* NULL after score-only evaluation. */
   const int *candidate, *status;
   const double *loss, *baseline, *r2;
   const double *candidate_loss; /* [candidate][response] pooled MSE */
   const double *fold_loss;      /* [fold][candidate][response] MSE */
   const double *fold_baseline;  /* [fold][response] baseline MSE */
} THD_picard_predict_result;
THD_picard_predict_plan *THD_picard_predict_plan_create(
    const THD_picard_plan *full, const struct THD_foldplan *folds,
    const double *fractions, int count, int nuisance_mode, int one_se,
    char *err, int nerr);
void THD_picard_predict_plan_free(THD_picard_predict_plan *p);
const THD_picard_predict_info *
THD_picard_predict_plan_info(const THD_picard_predict_plan *p);
const struct THD_foldplan *
THD_picard_predict_folds(const THD_picard_predict_plan *p);
THD_picard_predict_workspace *
THD_picard_predict_workspace_alloc(const THD_picard_predict_plan *p,
                                   int max_response, char *err, int nerr);
void THD_picard_predict_workspace_free(THD_picard_predict_workspace *w);
/* No fit-time allocations. Results are borrowed until next attempt or free.
   score_only avoids full refitting (useful for streamed global selection).
   forced_index=-1 selects per response; otherwise use that candidate for all
   responses. Forced choice still reports each response's own CV diagnostics. */
int THD_picard_predict_fit(THD_picard_predict_workspace *w, const double *y,
                           int q, int score_only, int forced_index, char *err,
                           int nerr);
const THD_picard_predict_result *
THD_picard_predict_get_result(const THD_picard_predict_workspace *w);
/* Select from finite, nonnegative [fold][candidate] losses. Used by both
   response-specific and global selection. A global caller must explicitly
   define the voxel set and normalization before averaging these losses. */
int THD_picard_predict_select(const THD_picard_predict_plan *p,
                              const double *fold_losses, int *winner,
                              double *best_se, char *err, int nerr);
/* Bounded global aggregation belongs to the numeric workflow. Call add in
   ascending voxel order for batch/thread-invariant sums. The caller supplies
   the explicit voxel mask by choosing which responses enter score batches.
   Normalize each voxel's fold MSE by its pooled baseline MSE; exclude exactly
   zero baselines. If none are defined, select the largest fraction. */
typedef struct THD_picard_predict_global THD_picard_predict_global;
typedef struct {
   int candidate, nvalid;
   double best_se;
   const double *fold_loss; /* [fold][candidate] mean normalized loss */
} THD_picard_predict_global_result;
THD_picard_predict_global *
THD_picard_predict_global_alloc(const THD_picard_predict_plan *p, char *err,
                                int nerr);
void THD_picard_predict_global_free(THD_picard_predict_global *g);
int THD_picard_predict_global_add(THD_picard_predict_global *g,
                                  const THD_picard_predict_workspace *w,
                                  char *err, int nerr);
const THD_picard_predict_global_result *
THD_picard_predict_global_finish(THD_picard_predict_global *g, char *err,
                                 int nerr);

/* ========================================================================
   4. CANDIDATE DESIGNS AND LOCKED HRF CHOICES
   ======================================================================== */

#define THD_PICARD_SELECT_HRF 0
#define THD_PICARD_SELECT_PREDICT 1

typedef struct THD_picard_candidates THD_picard_candidates;
typedef struct THD_picard_candidates_workspace THD_picard_candidates_workspace;
typedef struct {
   int policy, nparameter, nuisance_mode, one_se;
   double fixed_parameter;
   const double *parameters;
   const THD_picard_trial *trials;
   const struct THD_foldplan *folds;
} THD_picard_candidates_config;

/* Numeric collection (one or more designs), with borrowed fixed plans that must
   outlive it. Copies IDs/configuration and prepares its own tuning plans. All
   designs must have identical dimensions, column roles/order, retained row/run
   IDs, nuisance values, and estimator settings. The file adapter additionally
   validates labels, TR, run boundaries and explicit coefficient identities.
   IDs are unique nonnegative integers <=16777215, exactly representable in
   float output maps. Exact score ties choose the smallest stable ID.

   HRF: select by minimum full-data unpenalized residual MSE (weighted by the
   common R^-1 for explicit covariance), requiring equal projected
   rank. This is equivalent to maximizing fit R2 with a common denominator;
   it is NOT held-out performance. Then fit the selected design at a fixed
   parameter or tune its repeated-condition betas. Full-data HRF choice means
   subsequent beta tuning is conditional on that choice, not independent CV.

   PREDICT: each design independently selects its parameter using the supplied
   folds and existing min/one-SE policy. Choose the design with the smallest
   selected prediction MSE. The one-SE rule applies WITHIN designs, never to
   an undefined ordering of HRFs. This reuses validation data for selection;
   no independent outer-test performance is estimated. Global selection is
   not part of this per-response collection interface. */
THD_picard_candidates *THD_picard_candidates_create(
    const THD_picard_plan *const *designs, const int *ids, int count,
    const THD_picard_candidates_config *config, char *err, int nerr);
void THD_picard_candidates_free(THD_picard_candidates *p);
THD_picard_candidates_workspace *
THD_picard_candidates_workspace_alloc(const THD_picard_candidates *p,
                                      int max_response, char *err, int nerr);
void THD_picard_candidates_workspace_free(THD_picard_candidates_workspace *w);

typedef struct {
   int q;
   const THD_picard_result *fit;
   const int *design_index, *design_id, *rank;
   const double *selection_loss;
   const double *design_loss; /* [design][response], policy-specific MSE */
   /* Selected design's tuning diagnostics; NULL for fixed-parameter HRF. */
   const int *parameter_index, *prediction_status;
   const double *loss, *baseline, *r2;
   const double *parameter_loss; /* [parameter][response] */
} THD_picard_candidates_result;
/* No fit-time allocation or design x response x coefficient tensor. Returned
   views borrow worker storage, invalidated on ANY fit attempt or destruction.
 */
int THD_picard_candidates_fit(THD_picard_candidates_workspace *w,
                              const double *y, int q, char *err, int nerr);
/* Refit previously selected HRFs without running selection again. indices[q]
   and scores[q] come from the upstream HRF-only stage. HRF policy only;
   design_loss is NULL because no new design comparison was performed. */
int THD_picard_candidates_fit_selected(THD_picard_candidates_workspace *w,
                                       const double *y, int q,
                                       const int *indices, const double *scores,
                                       char *err, int nerr);
const THD_picard_candidates_result *
THD_picard_candidates_get_result(const THD_picard_candidates_workspace *w);

/* ========================================================================
   5. NOISE PCS AND SINGLE-TRIAL PC SELECTION
   ======================================================================== */

/* Run-specific temporal PCA. The accumulator borrows an immutable base plan
   and receives pool responses in bounded batches. It removes a constant and
   the supplied run-local nuisance span, normalizes each residual time series
   to unit Euclidean norm, and accumulates temporal covariance per run.
   No sample or covariance is shared across runs. No task residualization,
   intensity scaling, implicit drift polynomials or whitening is added here. */
typedef struct THD_picard_noise THD_picard_noise;
THD_picard_noise *THD_picard_noise_create(const THD_picard_plan *base,
                                          int max_pc, int batch, char *err,
                                          int nerr);
int THD_picard_noise_add(THD_picard_noise *noise, const double *y, int q,
                         char *err, int nerr);
int THD_picard_noise_finish(THD_picard_noise *noise, char *err, int nerr);
void THD_picard_noise_free(THD_picard_noise *noise);
/* Borrowed views, available only after successful finish. pcs[n][max_pc]
   stores the components of each row's run (not a dense zero-padded design).
   Components have unit sample SD, with largest absolute entry positive.
   Rank cutoff is max(rank_tol^2,64*epsilon*run_rows) times largest eigenvalue.
   Insufficient nonzero rank is an error, never silent PC-count clipping. */
const double *THD_picard_noise_pcs(const THD_picard_noise *noise);
int THD_picard_noise_runs(const THD_picard_noise *noise);
const int *THD_picard_noise_ranks(const THD_picard_noise *noise);
/* Own a new design with count PCs per run appended as unpenalized columns.
   Original task column indices/roles/identities are unchanged. ols=1 uses
   fractional ridge's unpenalized reference and disables sparse support refit;
   otherwise preserve all original estimator/whitening settings. */
THD_picard_plan *THD_picard_noise_plan(const THD_picard_noise *noise,
                                       const THD_picard_plan *base, int count,
                                       int ols, char *err, int nerr);
/* Validate one run-local regressor per trial and run-local nuisance columns.
   Fixed staged fitting need not have repeats; tuning checks repeats itself. */
int THD_picard_single_trials(const THD_picard_plan *p,
                             const THD_picard_trial *trials, char *err,
                             int nerr);
/* PC-selection error: all ordered cross-run same-condition pairs against a
   fixed NO-PC reference, normalized by its sample beta variance. One session,
   no output beta rescaling. Undefined variance gives -1 (exclude from median).
   Returns mean squared pair error; multiply by pair count for GLMsingle SSE.
   These are conditional selection errors, not independent test performance. */
int THD_picard_single_loss(const THD_picard_trial *trials, int p,
                           const double *reference, const double *beta, int q,
                           double *loss, char *err, int nerr);
/* loss[] is an increasing-count sequence beginning with zero PCs. Choose the
   first count achieving at least 1/stop of the best improvement over zero.
   stop>=1; no improvement chooses zero. Identical policy to the pinned
   select_noise_regressors helper when its performance is negative loss. */
int THD_picard_single_select(const double *loss, int count, double stop);

/* ========================================================================
   PRIVATE IMPLEMENTATION CONTRACT
   ======================================================================== */

/* Only thd_picard.c and its AFNI adapter need this layout. Define the macro
   before first including this header. Keeping it here avoids a separate
   private header without exposing the representation to ordinary clients. */
#ifdef THD_PICARD_PRIVATE
/* Private immutable fixed-design storage shared by the numerical stages.
   Applications use the checked public info/result views instead. */
struct THD_picard_plan {
   THD_picard_design_info info;
   int *pen, *nuis, *rows, *runs;
   double *x, *z, *effective, *weighted_z;
   double *covariance_lower;
   THD_banded_design *banded;
   THD_banded_config band_config;
   int *groups;
   double *penalties;
   THD_linproj *projection;
   THD_regression_design *design;
   THD_elastic_design *elastic;
   THD_elastic_config config;
};
/* Internal path step: warm only within a single fold/response batch. */
int THD_picard_fit_step(THD_picard_workspace *w, const double *y, int q,
                        const double *parameters, int warm, char *err,
                        int nerr);
int THD_picard_fit_ols(THD_picard_workspace *w, const double *y, int q,
                       char *err, int nerr);
/* Constructor dispatch preserves preprocessing independently in each fold. */
THD_picard_plan *THD_picard_plan_like(const THD_picard_plan *full,
                                      const double *matrix, int n, int d,
                                      const int *roles, const int *rows,
                                      const int *penmap, const int *runs,
                                      int nfull, double tol, char *err,
                                      int nerr);
int THD_picard_grid_valid(const THD_picard_plan *p, const double *grid,
                          int count, char *err, int nerr);
/* Larger strength always means more regularization, regardless of units. */
double THD_picard_strength(const THD_picard_plan *p, double parameter);
#endif

#endif /* THD_PICARD_H */
