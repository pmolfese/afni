#ifndef _THD_MVPA_HEADER_
#define _THD_MVPA_HEADER_

#include "thd_crossval.h"
#include "thd_permute.h"

/*----------------------------------------------------------------------------
  Dataset-independent numerical support for multivariate prediction.  Fitted
  objects contain only statistics learned from supplied training rows;
  applying them never recomputes held-out statistics.
                                                     -- P Molfese, Oct 2026
------------------------------------------------------------------------------*/

#define THD_PREP_NONE             0
#define THD_PREP_FEATURE_ZSCORE   1
#define THD_PREP_SAMPLE_CENTER    2
#define THD_PREP_RUN_CENTER       3

#define THD_METRIC_ACCURACY              0
#define THD_METRIC_BALANCED_ACCURACY     1
#define THD_METRIC_PEARSON               2
#define THD_METRIC_SPEARMAN              3
#define THD_METRIC_CROSS_VALIDATED_R2    4
#define THD_METRIC_MAE                   5
#define THD_METRIC_RMSE                  6

#define THD_CENTROID_EUCLIDEAN           0
#define THD_CENTROID_CORRELATION         1
#define THD_CENTROID_COSINE              2

typedef struct {
   int method , nfeature , nfit , ngroup , nconstant ;
   int nfeature_alloc , ngroup_alloc ;

   int *group_value ;       /* [ngroup], training group IDs for run_center  */
   int *group_count ;       /* [ngroup], reusable run_center counts         */
   double *center ;         /* [nfeature] or [ngroup][nfeature]             */
   double *scale ;          /* [nfeature], feature_zscore only              */
   unsigned char *constant ;/* [nfeature], zero-variance training features  */
} THD_preprocessor ;

typedef struct {
   int nsample , nclass , ncorrect , nclass_present ;
   double accuracy , balanced_accuracy ;
} THD_mvpa_class_score ;

typedef struct {
   int nsample ;
   double pearson , spearman , cv_r2 , mae , rmse ;
} THD_mvpa_continuous_score ;

typedef struct {
   int distance , nclass , nfeature , nfit , nfeature_alloc ;
   int *class_count ;       /* [nclass]                                  */
   double *centroid ;       /* [nclass][nfeature]                         */
   double *centroid_mean ;  /* [nclass], correlation distance only       */
   double *centroid_norm ;  /* [nclass], correlation or cosine distance  */
} THD_centroid_model ;

typedef struct {
   const THD_foldplan *plan ; /* borrowed; must outlive the workspace       */
   const int *group_id ;      /* borrowed run_center groups, otherwise NULL */
   int nsample , nfold , nclass , nfeature_alloc , max_train , max_test ;
   int preprocess , ntest ;

   float *train_x , *test_x ;
   int *train_class , *test_class , *test_prediction ;
   double *test_decision ;

   THD_preprocessor *prep ;
   THD_centroid_model *centroid ;
   int *train_count ;
   int *fold_confusion ;
   double *fold_recall ;

   int *prediction ;        /* [nsample], -1 when a row was not tested      */
   double *decision ;       /* [nsample][nclass], NAN when not tested       */
   int *fold_nsample ;      /* [nfold]                                      */
   double *fold_accuracy ;  /* [nfold]                                      */
   double *fold_balanced_accuracy ; /* [nfold]                               */
   int *confusion ;         /* [nclass][nclass], truth rows                 */
   double *recall ;         /* [nclass]                                     */
   THD_mvpa_class_score score ;
} THD_mvpa_workspace ;

/*! Allocate buffers for repeated analyses with one fold plan.  max_feature is
    the largest pattern width that will be passed to the runner.  For
    run_center, group_id is borrowed and must also outlive the workspace. */
extern THD_mvpa_workspace * THD_mvpa_workspace_alloc(
                                      const THD_foldplan *plan,
                                      int max_feature, int nclass,
                                      int preprocess, const int *group_id,
                                      char *err, int nerr ) ;

/*! Run fold-local preprocessing, centroid fitting, held-out prediction, and
    classification scoring.  Results are workspace-owned and remain valid
    until the next call.  score pools all tested rows; fold scores preserve
    the inputs needed for later subject/stratum-aware aggregation.  group_id
    is required only for run_center and must match workspace allocation. */
extern int THD_mvpa_run_centroid_cv( const THD_foldplan *plan,
                                     int distance, const float *x,
                                     int nfeature, const int *class_id,
                                     const int *group_id,
                                     THD_mvpa_workspace *workspace,
                                     char *err, int nerr ) ;

/*! Repeat the complete centroid CV analysis under every synchronized label
    relabeling in pset.  null_score[pset->nperm] receives accuracy or balanced
    accuracy for each slot, including the identity in slot 0.  permuted_class
    is caller-owned scratch [plan->nsample].  On success the identity analysis
    is restored in workspace, so its predictions and fold diagnostics remain
    the observed result. */
extern int THD_mvpa_run_centroid_permutations(
                                     const THD_foldplan *plan,
                                     int distance, const float *x,
                                     int nfeature, const int *class_id,
                                     const int *group_id, int metric,
                                     const PERM_set *pset,
                                     int *permuted_class, float *null_score,
                                     THD_mvpa_workspace *workspace,
                                     char *err, int nerr ) ;

extern void THD_mvpa_workspace_free( THD_mvpa_workspace *workspace ) ;

/*! Centroid distances: euclidean, correlation, cosine. */
extern const char * THD_centroid_distance_name( int distance ) ;
extern int THD_centroid_distance_parse( const char *name, int *distance,
                                        char *err, int nerr ) ;

/*! Fit class centroids from a compact [nsample][nfeature] training matrix.
    class_id values must cover every integer in [0,nclass-1]. */
extern THD_centroid_model * THD_centroid_fit(
                                        int distance, const float *x,
                                        const int *class_id,
                                        int nsample, int nfeature, int nclass,
                                        char *err, int nerr ) ;

/*! Predict compact rows.  decision may be NULL; otherwise it receives
    [nsample][nclass] scores where larger is always better: negative squared
    Euclidean distance, Pearson correlation, or cosine similarity.  Exact
    ties resolve to the lowest class ID. */
extern int THD_centroid_predict( const THD_centroid_model *model,
                                 const float *x, int nsample, int nfeature,
                                 int *prediction, double *decision,
                                 char *err, int nerr ) ;

extern void THD_centroid_free( THD_centroid_model *model ) ;

/*! Scalar metric option names.  Per-class recall and the confusion matrix
    are diagnostics returned with classification scores, not scalar metrics. */
extern const char * THD_mvpa_metric_name( int metric ) ;
extern int THD_mvpa_metric_parse( const char *name, int *metric,
                                  char *err, int nerr ) ;

/*! Score integer class predictions.  confusion is caller-allocated
    [nclass][nclass], indexed [true class][predicted class], and recall is
    [nclass].  Recall for a class absent from truth is NAN and that class is
    excluded from balanced_accuracy. */
extern int THD_mvpa_score_classification(
                                  const int *truth, const int *prediction,
                                  int nsample, int nclass,
                                  int *confusion, double *recall,
                                  THD_mvpa_class_score *score,
                                  char *err, int nerr ) ;

/*! Score continuous predictions.  baseline[nsample] must contain predictions
    made without the held-out outcomes (normally each fold's training mean),
    so cv_r2 is 1-SSE(prediction)/SSE(baseline).  rank_scratch needs 2*nsample
    doubles and index_scratch needs 2*nsample ints.  Undefined correlations
    or cv_r2 are returned as NAN. */
extern int THD_mvpa_score_continuous(
                                  const float *truth, const float *prediction,
                                  const float *baseline, int nsample,
                                  double *rank_scratch, int *index_scratch,
                                  THD_mvpa_continuous_score *score,
                                  char *err, int nerr ) ;

/*! Exact public option names: none, feature_zscore, sample_center, run_center. */
extern const char * THD_preprocess_name( int method ) ;
extern int THD_preprocess_parse( const char *name, int *method,
                                 char *err, int nerr ) ;

/*! Fit using only fit_index[0..nfit-1].  x is [nsample][nfeature].
    group_id[nsample] is required only for run_center.  feature_zscore uses
    the population standard deviation of the training rows.  A feature that
    is constant to numerical precision in training is mapped to zero. */
extern THD_preprocessor * THD_preprocess_fit(
                                      int method, const float *x,
                                      int nsample, int nfeature,
                                      const int *fit_index, int nfit,
                                      const int *group_id,
                                      char *err, int nerr ) ;

/*! Apply a fitted transformation to selected rows.  out is caller-allocated
    [napply][nfeature], in apply_index order, and must not alias x.  run_center
    rejects a group not observed while fitting rather than estimating a mean
    from held-out samples. */
extern int THD_preprocess_apply( const THD_preprocessor *prep,
                                 const float *x, int nsample, int nfeature,
                                 const int *apply_index, int napply,
                                 const int *group_id, float *out,
                                 char *err, int nerr ) ;

/*! Fit on plan->train[fold], then transform train and test into compact caller-
    allocated matrices.  fitted may be NULL; otherwise the caller owns it. */
extern int THD_preprocess_fold( const THD_foldplan *plan, int fold, int method,
                                const float *x, int nfeature,
                                const int *group_id,
                                float *train_out, float *test_out,
                                THD_preprocessor **fitted,
                                char *err, int nerr ) ;

extern void THD_preprocess_free( THD_preprocessor *prep ) ;

#endif /* _THD_MVPA_HEADER_ */
