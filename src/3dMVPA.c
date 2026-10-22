#include "mrilib.h"
#include "thd_crossval.h"
#include "thd_mapinfer.h"
#include "thd_mvpa.h"
#include "thd_patterns.h"
#include "thd_permute.h"
#include "thd_sampleset.h"
#include "thd_xmat.h"

#include <errno.h>
#include <float.h>
#include <limits.h>

#ifdef USE_OMP
#include <omp.h>
#endif

#ifdef USE_SUMA
#include "SUMA_suma.h"
#endif

#define PROGRAM_NAME "3dMVPA"

/*----------------------------------------------------------------------------
  Classification and continuous prediction in ROIs and spatial neighborhoods.
  Compact input adapters end at THD_sampleset; everything after that is shared.
                                                     -- P Molfese, Oct 2026
------------------------------------------------------------------------------*/

static char *mvpa_options[]={
   "-mode","-dataTable","-dataTableFile","-trial_manifest","-xmat_stims",
   "-target_column","-target_columns","-design_columns","-lambda","-score",
   "-lambda_grid","-inner_cv_column","-inner_cv_k",
   "-shrinkage","-class_prior","-variance_floor",
   "-svm_c","-svm_c_grid","-svm_tol","-svm_max_iter",
   "-cv_column","-cv","-cv_k","-within_column",
   "-subject_column","-run_column","-item_column","-fold_column","-sample_column",
   "-mask","-roi_sel","-searchlight","-surf","-method","-distance","-scale","-metric","-prefix",
   "-nperm","-seed","-permute_within","-block_column",
   "-save_predictions","-save_fold_scores","-save_coefficients","-coef_lambda",
   "-show_table","-show_samples","-show_folds","-quiet","-help","-h",NULL
} ;

static void usage_3dMVPA( void )
{
   printf(
"\n"
"Usage: 3dMVPA [options]                                      ~1~\n"
"\n"
"Cross-validated classification, continuous decoding, and ridge encoding in\n"
"ROIs or volume/surface neighborhoods. Preprocessing, intercepts, and model\n"
"coefficients are fitted separately inside every training fold.\n"
"\n"
"Required options:                                           ~1~\n"
"\n"
"  -mode classify|decode|encode\n"
"\n"
"  -dataTable @TABLE\n"
"  -dataTableFile TABLE\n"
"  -dataTable Subj Run Condition InputFile ...\n"
"                         One row per sample. InputFile must select exactly\n"
"                         one sub-brick when its dataset has more than one.\n"
"                         A subject-level table with Xmat and whole-bucket\n"
"                         InputFile columns is expanded to the same sample\n"
"                         representation before analysis.\n"
"  -trial_manifest FILE   Alternative to -dataTable: one row per trial beta.\n"
"                         InputFile names a whole LSS or LSA beta bucket;\n"
"                         Volume is its zero-based sub-brick index. Supply\n"
"                         trial, run and target metadata explicitly. Optional\n"
"                         Xmat and XColumn must occur together; XColumn is\n"
"                         the zero-based original design-matrix column. These\n"
"                         indices are checked; metadata is retained unchanged.\n"
"                         No beta-label matching or Xmat expansion is applied.\n"
"                         Cannot accompany -dataTable, -dataTableFile, or\n"
"                         -xmat_stims.\n"
"\n"
"  -target_column NAME    Categorical classification or numeric decoding target.\n"
"  -target_columns A,B    Multiple numeric decoding targets.\n"
"  -design_columns A,B    Numeric predictors for encoding neural responses.\n"
"  -score voxel|pattern  Encoding output: voxel [default] predicts each masked\n"
"                         voxel; pattern averages per-response metrics within\n"
"                         each ROI/searchlight, excluding undefined values.\n"
"  -cv_column NAME        Leave one level of NAME out at a time.\n"
"  -mask DSET             Without -searchlight, positive values define ROIs and\n"
"                         a binary mask is one ROI. With -searchlight, nonzero\n"
"                         voxels define eligible centers and neighborhood voxels.\n"
"                         The mask and every sample must share a grid or exact\n"
"                         surface node-row domain. A surface searchlight may omit\n"
"                         this option to search every represented sample node.\n"
"\n"
"Cross-validation:                                          ~1~\n"
"\n"
"  -cv run|subject|item   Shorthand for -cv_column Run, Subj, or Item.\n"
"  -cv_k K               Assign whole -cv_column groups to K balanced folds\n"
"                         instead of leave-one-group-out.\n"
"  -within_column NAME    Build independent folds inside each NAME level.\n"
"                         For within-subject run CV, use Subj here.\n"
"\n"
"Sample roles:                                               ~1~\n"
"\n"
"  -subject_column NAME   Subject role [auto-detect Subj].\n"
"  -run_column NAME       Run role [auto-detect Run].\n"
"  -item_column NAME      Item/trial role [auto-detect Trial or Item].\n"
"  -fold_column NAME      Fold role [auto-detect Fold].\n"
"  -sample_column NAME    Unique sample labels [auto-detect Sample; otherwise\n"
"                         row000001, ...].\n"
"  -xmat_stims LIST       Comma-separated X-matrix stimuli to expand [all].\n"
"\n"
"Model options:                                              ~1~\n"
"\n"
"  -method METHOD         Classify: centroid [default], ridge, svm, lda_diag, or\n"
"                         lda_shrinkage. Ridge is required for decode/encode.\n"
"                         Ridge classification fits one 0/1 response per class\n"
"                         with argmax decisions; exact ties choose the first\n"
"                         encoded class. Decisions are not probabilities.\n"
"  -svm_c C|auto          Native linear SVM only: positive cost [1]. Each\n"
"                         one-vs-rest model minimizes .5*(||w||^2+b^2) +\n"
"                         C*sum(max(0,1-y*(Xw+b))^2). The intercept b IS\n"
"                         penalized (constant feature = 1); no class weights.\n"
"                         Larger C penalizes margin violations more strongly.\n"
"                         Raw class margins are not probabilities; argmax wins,\n"
"                         exact ties choose the first encoded class.\n"
"  -svm_c_grid C1,C2,... Positive distinct candidates; enables grouped nested\n"
"                         selection [auto default: 0.01,0.1,1,10,100]. Cannot\n"
"                         accompany numeric -svm_c. Uses -inner_cv_column and\n"
"                         -inner_cv_k; minimize 1-ACC/BACC (selected -metric)\n"
"                         over pooled inner held-out rows. Ties choose SMALLER\n"
"                         C. Every inner training set needs every class.\n"
"                         Scaling, selection and fitting repeat per permutation.\n"
"                         PREFIX.tuning.tsv records all candidate losses.\n"
"  -svm_tol T            Absolute projected-dual-gradient tolerance [1e-6],\n"
"                         range [1e-12,1e-2]. Deterministic cyclic coordinate\n"
"                         descent verifies final KKT residual for every class.\n"
"  -svm_max_iter N       Maximum sweeps per class [10000]. Nonconvergence is\n"
"                         an error; increase this limit, rescale, or reduce C.\n"
"                         Fold audits append C, InnerNFold, InnerLoss, maximum\n"
"                         class SolverIterations and SolverResidual. SVM is\n"
"                         classification-only; descriptive coefficient export\n"
"                         currently remains ridge-only.\n"
"  -shrinkage A           lda_shrinkage only: fixed covariance shrinkage in\n"
"                         (0,1] [0.1]. Use (1-A)*S + max(A,F)*T*I, where S is\n"
"                         pooled within-class covariance, T its mean variance,\n"
"                         and F the relative variance floor. No CV tuning.\n"
"  -variance_floor F      LDA only: relative variance floor [1e-6], allowed\n"
"                         range [1e-12,1]. lda_diag floors each variance at\n"
"                         F*T. If all within-class variances are zero, T uses\n"
"                         overall training variance, then 1 for constant data.\n"
"  -class_prior P         LDA only: equal [default] or empirical (training\n"
"                         class proportions). LDA needs more training rows\n"
"                         than classes in every outer fold. Slopes, means,\n"
"                         covariances and priors use training rows only.\n"
"                         DecisionClass scores are centered linear discriminants,\n"
"                         not probabilities; exact ties choose the first class.\n"
"                         Fold audits append Shrinkage, CovarianceScale,\n"
"                         VarianceFloor (absolute), and CovarianceForm.\n"
"  -lambda L              Positive ridge penalty [1]. Slopes are penalized;\n"
"                         the training-fitted intercept is unpenalized.\n"
"                         Use -lambda auto for nested grouped selection.\n"
"  -lambda_grid L1,L2,... Positive distinct candidates; enables nested tuning\n"
"                         [auto default: 0.01,0.1,1,10,100]. Cannot accompany\n"
"                         a numeric -lambda. One penalty per outer fold/unit\n"
"                         is shared across responses. Minimize mean response\n"
"                         SSE / training-mean-baseline SSE on inner held-out\n"
"                         rows; exclude constant responses. Exact ties or no\n"
"                         informative responses select the largest penalty.\n"
"                         For classification, responses are class indicators:\n"
"                         tuning uses normalized class-wise squared error,\n"
"                         not accuracy/BACC. All inner training sets must\n"
"                         contain every class, including under permutations.\n"
"  -inner_cv_column NAME Inner grouping column [outer -cv_column].\n"
"  -inner_cv_k K          Balance whole inner groups into K folds [default:\n"
"                         leave one group out]. Every outer training set must\n"
"                         contain at least two inner groups (or K if given).\n"
"                         Shared by ridge and SVM tuning. Preprocessing is\n"
"                         refitted within every inner split.\n"
"                         PREFIX.tuning.tsv records all candidate losses;\n"
"                         fold audits append selected penalties and losses.\n"
"  -distance D            Centroid only: euclidean [default], correlation, cosine.\n"
"  -scale S               feature_zscore [default], sample_center, run_center,\n"
"                         or none. All learned values are training-fold only.\n"
"                         run_center requires a run role and cannot be applied\n"
"                         to a run absent from that fold's training data.\n"
"  -metric M              balanced_accuracy [default] or accuracy. This selects\n"
"                         the SCORE map; the ROI table reports both.\n"
"                         Continuous modes default to cross_validated_r2 and\n"
"                         also accept pearson, spearman, MAE, RMSE. R2 uses\n"
"                         each fold's training response mean as its baseline.\n"
"                         All five metrics are saved per target; encoding\n"
"                         pattern summaries weight responses equally. Undefined\n"
"                         scores are NaN in tables and zero in maps; NVALID\n"
"                         reports the primary metric's defined response count.\n"
"\n"
"Permutation inference:                                    ~1~\n"
"\n"
"  -nperm N               Run N synchronized label/row relabelings [0: off]. The\n"
"                         identity is included. If the allowed group has at\n"
"                         most N distinct relabelings, it is enumerated exactly.\n"
"  -permute_within COL    Relabel only among rows sharing COL. Required\n"
"                         with -nperm; use 'none' only for an explicitly\n"
"                         unrestricted null. -block_column is an alias. Choose a\n"
"                         column that fully identifies the exchangeability unit,\n"
"                         such as Run, Subj, Fold, or a composite block column.\n"
"  -seed S                Monte-Carlo relabeling seed [1234567]. Exact tests do\n"
"                         not depend on this value. Classification must leave\n"
"                         every CV training fold with all target classes.\n"
"                         Continuous modes shuffle whole target rows (decode)\n"
"                         or design rows (encode), holding neural rows fixed.\n"
"                         All attribute columns move together; blocks and CV\n"
"                         metadata stay fixed. Rows must be exchangeable under\n"
"                         this null; serial time-series dependencies require\n"
"                         an appropriate design, not unrestricted row shuffling.\n"
"                         Preprocessing, tuning, and fitting repeat per null.\n"
"                         Each decoding target has a separate spatial family;\n"
"                         encoding tests its voxel or mean-pattern score.\n"
"                         MAE/RMSE use the lower tail; other metrics the upper.\n"
"                         PREFIX.inference.tsv reports all continuous results.\n"
"                         Buckets add NULL_MEAN, ABOVE_NULL, PVAL/ZSCR/QVAL,\n"
"                         PFWE/ZFWE, and NNULL_VALID per target. ABOVE_NULL is\n"
"                         positive for improvement, including lower errors.\n"
"                         Undefined observations have neutral p=1 and z=0;\n"
"                         undefined nulls count conservatively as exceedances.\n"
"\n"
"Classification uses the selected metric's upper-tail label-permutation null.\n"
"The same relabeling slot drives every ROI/searchlight center. Output adds the\n"
"empirical score above null-chance, PVAL/ZSCR, BH QVAL, and max-statistic\n"
"PFWE/ZFWE bricks. ZSCR and ZFWE carry AFNI FIZT statistical codes.\n"
"\n"
"Spatial and output options:                                ~1~\n"
"\n"
"  -roi_sel LIST          Restrict positive mask values with an AFNI int-list.\n"
"                         This option cannot be combined with -searchlight.\n"
"  -searchlight NBHD      Classify a moving volumetric neighborhood at every\n"
"                         nonzero mask voxel. NBHD is a radius in mm or one of\n"
"                         SPHERE(r), RECT(a,b,c), RHDD(r), or TOHD(r).\n"
"                         OpenMP builds use one workspace per worker; set\n"
"                         OMP_NUM_THREADS to control the worker count.\n"
"  -surf SURFACE          Make -searchlight a geodesic surface searchlight on\n"
"                         this mesh (.gii/.asc/.ply/...). With -surf, NBHD must\n"
"                         be one positive radius in mm. Inputs and an optional\n"
"                         mask must be NIML/GIFTI datasets with identical dense\n"
"                         or sparse node-row domains. The output preserves sparse\n"
"                         node indexes. Requires a build with SUMA support.\n"
"  -prefix PREFIX         Output dataset prefix [MVPA].\n"
"                         Without inference the dataset contains SCORE, NFEAT,\n"
"                         and NFOLD. -nperm adds ABOVE_CHANCE, PVAL, ZSCR,\n"
"                         QVAL, PFWE, and ZFWE before NFEAT and NFOLD.\n"
"                         ROI analysis also writes PREFIX.roi.tsv with SCORE,\n"
"                         ACC, BACC, NFEAT, and NFOLD for every ROI.\n"
"  -save_coefficients CPREFIX  Ridge only: separate descriptive full-data refit.\n"
"                         Writes CPREFIX.coefficients.tsv and .preprocessing.tsv.\n"
"                         Encoding without searchlights writes CPREFIX maps: one\n"
"                         brick per design column plus REFIT_INTERCEPT. Other\n"
"                         analyses write CPREFIX.unit_ID per ROI/searchlight: one\n"
"                         brick per decoding target/class (encoding: per design\n"
"                         column plus intercept). Weights occupy neural feature/\n"
"                         response voxels/nodes, preserving sparse node IDs.\n"
"                         Overlapping neighborhoods retain separate models;\n"
"                         this can create many datasets and large audit files.\n"
"                         TSV identifies unit, feature, response, map and brick;\n"
"                         intercepts and all fitted preprocessing are retained.\n"
"                         y = intercept + sum(preprocessed_X * coefficient).\n"
"                         Weights use the requested -scale predictor units;\n"
"                         responses are unscaled (class responses: 0/1).\n"
"                         feature_zscore uses full-data population SD; sample_\n"
"                         center subtracts each row mean; run_center subtracts\n"
"                         full-data subject/run means recorded in preprocessing.\n"
"                         These are descriptive predictive weights, NOT held-out\n"
"                         scores, activation/Haufe maps, or significance tests.\n"
"                         Correlated features and scaling affect weight values.\n"
"  -coef_lambda L         Positive full-data refit penalty [fixed -lambda].\n"
"                         Required with tuned CV; never changes CV penalties,\n"
"                         scores, predictions or permutation inference.\n"
"  -save_predictions FILE Write one held-out prediction row per sample and\n"
"                         spatial unit, including sample metadata, test fold,\n"
"                         observed/predicted class labels, and every class\n"
"                         decision value. Rows never tested are retained as NA.\n"
"  -save_fold_scores FILE Write ACC and BACC for every fold and spatial unit,\n"
"                         with stratum and train/test sample counts. These two\n"
"                         audit tables are opt-in because searchlight tables\n"
"                         can be very large.\n"
"  -show_table            Print the canonical table after input conversion.\n"
"  -show_samples          Print resolved datasets, bricks, and sample roles.\n"
"  -show_folds            Print every training/test row in every fold.\n"
"  -quiet                 Suppress informational messages.\n"
"\n"
"Example:                                                    ~1~\n"
"\n"
"  3dMVPA -mode classify -dataTable @samples.tsv             \\\n"
"         -target_column Condition -cv_column Run             \\\n"
"         -mask brain+tlrc -searchlight 'SPHERE(5)'           \\\n"
"         -nperm 5000 -permute_within Run -seed 314159        \\\n"
"         -prefix face_scene\n"
"\n"
"Without -nperm, the ROI table and maps remain descriptive held-out results.\n"
"\n"
"Trial betas from 3dLSS or LSA:                               ~1~\n"
"\n"
"3dMVPA analyzes existing beta maps; it does not fit the BOLD time series.\n"
"Choose either LSS estimates or jointly fitted LSA estimates upstream.\n"
"For example, create one LSS bucket per run using a design containing exactly\n"
"one -stim_times_IM stimulus with all trials in that run:\n"
"\n"
"  3dLSS -matrix X.r1.xmat.1D -input run1+tlrc -prefix LSS.r1\n"
"  3dLSS -matrix X.r2.xmat.1D -input run2+tlrc -prefix LSS.r2\n"
"\n"
"The input must be the full-length preprocessed time series corresponding to\n"
"that design; 3dLSS applies its GoodList censoring. It writes only the selected\n"
"IM stimulus's estimates, in design-column order. Volume 0 therefore maps to\n"
"that stimulus's first column, not necessarily design column 0.\n"
"\n"
"Create lss_trials.tsv (whitespace-separated; paths relative to the current\n"
"working directory). This small example assumes two trials per run:\n"
"\n"
"  Subj Run Item Condition Volume InputFile\n"
"  s01  1   t01  face      0      LSS.r1+tlrc\n"
"  s01  1   t02  house     1      LSS.r1+tlrc\n"
"  s01  2   t03  face      0      LSS.r2+tlrc\n"
"  s01  2   t04  house     1      LSS.r2+tlrc\n"
"\n"
"  3dMVPA -mode classify -trial_manifest lss_trials.tsv     \\\n"
"         -target_column Condition -cv_column Run           \\\n"
"         -mask rois+tlrc -show_table -show_samples -show_folds \\\n"
"         -prefix LSS_decode\n"
"\n"
"Condition comes from the original trial log: a single IM stimulus called\n"
"'trials' does not encode face/house labels. Supply continuous targets or\n"
"encoding predictors as additional manifest columns when needed.\n"
"\n"
"For design provenance, add BOTH Xmat and XColumn columns, for example:\n"
"\n"
"  Subj Run Item Condition Volume Xmat          XColumn InputFile\n"
"  s01  1   t01  face      0      X.r1.xmat.1D  2       LSS.r1+tlrc\n"
"  s01  1   t02  house     1      X.r1.xmat.1D  3       LSS.r1+tlrc\n"
"\n"
"Here columns 2 and 3 must actually be the corresponding trial regressors.\n"
"Range checks cannot establish that the chosen beta and regressor correspond.\n"
"Run and Condition remain explicit; neither is inferred from XColumn.\n"
"\n"
"For an LSA comparison, make lsa_trials.tsv with the same trial metadata but\n"
"InputFile naming the Deconvolve/REMLfit bucket and Volume selecting each\n"
"trial's actual _Coef brick (not a t-statistic, and not the design-column\n"
"index). Run the same command with that manifest and -prefix LSA_decode.\n"
"The existing -dataTable path with Xmat remains available for automatic\n"
"_Coef-label expansion of conventional regression buckets. LSS buckets use\n"
"different labels and should use the explicit manifest or sample table.\n"
"\n"
"Keep the trials, targets, folds and analysis settings matched when comparing\n"
"LSS and LSA. LSS can mix neighboring trials' responses; fold-local MVPA\n"
"preprocessing cannot remove upstream dependence. The example holds out whole\n"
"runs estimated separately. For multiple subjects, add -within_column Subj\n"
"for within-subject run CV. Estimator choice is an upstream modeling choice;\n"
"do not select it using the same held-out scores used for final evaluation.\n"
"\n") ;
}

static long mvpa_parse_long( const char *opt, const char *arg, long lo, long hi )
{
   char *ep ; long vv ; errno=0 ; vv=strtol(arg,&ep,10) ;
   if( errno==ERANGE || ep==arg || *ep!='\0' || vv<lo || vv>hi )
     ERROR_exit(PROGRAM_NAME ": %s needs an integer in [%ld,%ld], not '%s'",opt,lo,hi,arg) ;
   return vv ;
}

#ifdef USE_SUMA
static float mvpa_parse_positive_float( const char *opt, const char *arg )
{
   char *ep ; double vv ; errno=0 ; vv=strtod(arg,&ep) ;
   if( errno==ERANGE || ep==arg || *ep!='\0' || !isfinite(vv) ||
       vv<=0.0 || vv>FLT_MAX )
     ERROR_exit(PROGRAM_NAME ": %s needs one finite number > 0, not '%s'",opt,arg) ;
   return (float)vv ;
}
#endif

static void mvpa_write_roi_table( const char *fname, const THD_roilist *rl,
                                  const float *score, const float *acc,
                                  const float *bacc, const float *nfeat,
                                  const float *nfold, const char *metric,
                                  const char *target, char **class_label,
                                  int nclass, const PERM_result *pr,
                                  const float *qval,
                                  const char *permute_within )
{
   FILE *fp ; int rr,cc ;
   if( THD_deathcon() && THD_is_file((char *)fname) )
     ERROR_exit(PROGRAM_NAME ": ROI table '%s' already exists",fname) ;
   fp=fopen(fname,"w") ;
   if( fp==NULL ) ERROR_exit(PROGRAM_NAME ": cannot write ROI table '%s'",fname) ;
   fprintf(fp,"# 3dMVPA ROI classification; target=%s; primary_metric=%s\n",target,metric) ;
   fprintf(fp,"# Classes:") ;
   for( cc=0 ; cc<nclass ; cc++ ) fprintf(fp," %d=%s",cc,class_label[cc]) ;
   if( pr!=NULL )
     fprintf(fp,"\n# Label permutations: %d; within=%s; tail=upper; FWE_family=%d_ROIs",
             pr->nperm,permute_within,rl->nroi) ;
   fprintf(fp,"\nROI\tLabel\tSCORE\tACC\tBACC") ;
   if( pr!=NULL ) fputs("\tABOVE_CHANCE\tPVAL\tZSCR\tQVAL\tPFWE\tZFWE",fp) ;
   fputs("\tNFEAT\tNFOLD\n",fp) ;
   for( rr=0 ; rr<rl->nroi ; rr++ ){
     fprintf(fp,"%d\t%s\t%.9g\t%.9g\t%.9g",
             rl->val[rr],rl->lab[rr]?rl->lab[rr]:"-",score[rr],acc[rr],bacc[rr]) ;
     if( pr!=NULL )
       fprintf(fp,"\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g",
               pr->effect[rr],pr->p_unc[rr],pr->z_unc[rr],qval[rr],
               pr->p_fwe[rr],pr->z_fwe[rr]) ;
     fprintf(fp,"\t%.0f\t%.0f\n",nfeat[rr],nfold[rr]) ;
   }
   if( fclose(fp)!=0 ) ERROR_exit(PROGRAM_NAME ": error closing ROI table '%s'",fname) ;
}

static THD_3dim_dataset * mvpa_make_roi_map( THD_3dim_dataset *mset,
                                             const THD_roilist *rl,
                                             const char *prefix,
                                             const float *score,
                                             const float *nfeat,
                                             const float *nfold,
                                             const PERM_result *pr,
                                             const float *qval,
                                             const char *score_label,
                                             int argc, char **argv )
{
   const float *value[9] ;
   const char *label[9] ;
   char domain_err[512] ;
   int is_z[9]={0,0,0,0,0,0,0,0,0},nval=0,bb ;
   THD_3dim_dataset *oset=EDIT_empty_copy(mset) ;
   if( oset==NULL ) ERROR_exit(PROGRAM_NAME ": cannot allocate the output dataset") ;
   value[nval]=score ; label[nval++]=score_label ;
   if( pr!=NULL ){
     value[nval]=pr->effect ; label[nval++]="ABOVE_CHANCE" ;
     value[nval]=pr->p_unc  ; label[nval++]="PVAL" ;
     value[nval]=pr->z_unc  ; label[nval]="ZSCR" ; is_z[nval++]=1 ;
     value[nval]=qval       ; label[nval++]="QVAL" ;
     value[nval]=pr->p_fwe  ; label[nval++]="PFWE" ;
     value[nval]=pr->z_fwe  ; label[nval]="ZFWE" ; is_z[nval++]=1 ;
   }
   value[nval]=nfeat ; label[nval++]="NFEAT" ;
   value[nval]=nfold ; label[nval++]="NFOLD" ;
   EDIT_dset_items(oset,ADN_prefix,prefix,ADN_nvals,nval,ADN_ntt,0,
                   ADN_type,HEAD_FUNC_TYPE,ADN_func_type,FUNC_BUCK_TYPE,ADN_none) ;
   if( THD_deathcon() && THD_is_file(DSET_HEADNAME(oset)) )
     ERROR_exit(PROGRAM_NAME ": output dataset %s already exists",DSET_HEADNAME(oset)) ;
   for( bb=0 ; bb<nval ; bb++ ){
     EDIT_substitute_brick(oset,bb,MRI_float,NULL) ;
     THD_roilist_paint(DSET_ARRAY(oset,bb),rl,value[bb]) ;
     EDIT_BRICK_LABEL(oset,bb,(char *)label[bb]) ;
     if( is_z[bb] ) EDIT_BRICK_TO_FIZT(oset,bb) ;
     else EDIT_BRICK_TO_NOSTAT(oset,bb) ;
     EDIT_BRICK_FACTOR(oset,bb,0.0f) ;
   }
   THD_copy_labeltable_atr(oset->dblk,mset->dblk) ;
   if( !THD_copy_surface_node_domain(oset,mset,domain_err,sizeof(domain_err)) )
     ERROR_exit(PROGRAM_NAME ": %s",domain_err) ;
   tross_Copy_History(mset,oset) ; tross_Make_History(PROGRAM_NAME,argc,argv,oset) ;
   return oset ;
}

static void mvpa_tsv_text( FILE *fp, const char *text )
{
   const char *cc ; int quote=0 ;
   if( text==NULL ) text="NA" ;
   for( cc=text ; *cc!='\0' ; cc++ )
     if( *cc=='\t' || *cc=='\n' || *cc=='\r' || *cc=='\"' ){ quote=1 ; break ; }
   if( !quote ){ fputs(text,fp) ; return ; }
   fputc('\"',fp) ;
   for( cc=text ; *cc!='\0' ; cc++ ){
     if( *cc=='\"' ) fputc('\"',fp) ;
     fputc(*cc,fp) ;
   }
   fputc('\"',fp) ;
}

static FILE * mvpa_open_audit( const char *fname, const char *what )
{
   FILE *fp ;
   if( fname==NULL || fname[0]=='\0' )
     ERROR_exit(PROGRAM_NAME ": %s needs a nonempty filename",what) ;
   if( THD_deathcon() && THD_is_file((char *)fname) )
     ERROR_exit(PROGRAM_NAME ": %s file '%s' already exists",what,fname) ;
   fp=fopen(fname,"w") ;
   if( fp==NULL ) ERROR_exit(PROGRAM_NAME ": cannot write %s file '%s'",what,fname) ;
   return fp ;
}

static void mvpa_write_audit_header( FILE *fp, const char *kind,
                                     const char *target, char **class_label,
                                     int nclass, int searchlight, int surface )
{
   int cc ;
   fprintf(fp,"# 3dMVPA %s\n# Target\t",kind) ; mvpa_tsv_text(fp,target) ; fputc('\n',fp) ;
   for( cc=0 ; cc<nclass ; cc++ ){
     fprintf(fp,"# Class\t%d\t",cc) ; mvpa_tsv_text(fp,class_label[cc]) ; fputc('\n',fp) ;
   }
   if( surface )
     fputs("# Surface searchlight UnitID is the mesh node index; I/J/K are NA.\n",fp) ;
   else if( searchlight )
     fputs("# Searchlight UnitID and I/J/K are zero-based center indexes.\n",fp) ;
}

static void mvpa_write_spatial_prefix( FILE *fp, THD_3dim_dataset *mset,
                                       const THD_roilist *rl, int rr,
                                       int searchlight, int surface )
{
   int idx=rl->center[rr] ;
   fputs(surface?"Node\t":searchlight?"Voxel\t":"ROI\t",fp) ;
   fprintf(fp,"%d\t",rl->val[rr]) ;
   mvpa_tsv_text(fp,rl->lab[rr]?rl->lab[rr]:"NA") ; fputc('\t',fp) ;
   if( surface ) fputs("NA\tNA\tNA\t",fp) ;
   else if( searchlight ){
     int nx=DSET_NX(mset),nxy=nx*DSET_NY(mset) ;
     fprintf(fp,"%d\t%d\t%d\t",idx%nx,(idx%nxy)/nx,idx/nxy) ;
   } else fputs("NA\tNA\tNA\t",fp) ;
   fprintf(fp,"%d",rl->vox[rr].nar) ;
}

static const char * mvpa_level( char **label, int index )
{
   return (index>=0 && label!=NULL && label[index]!=NULL) ? label[index] : "NA" ;
}

static void mvpa_write_prediction_rows( FILE *fp, THD_3dim_dataset *mset,
                                        const THD_roilist *rl, int rr,
                                        int searchlight, int surface,
                                        const THD_sampleset *ss,
                                        const THD_foldplan *plan,
                                        const int *class_id, char **class_label,
                                        const THD_mvpa_workspace *workspace )
{
   int ii,cc ;
   for( ii=0 ; ii<ss->nsample ; ii++ ){
     int ff=plan->sample_test_fold[ii],pred=workspace->prediction[ii] ;
     mvpa_write_spatial_prefix(fp,mset,rl,rr,searchlight,surface) ;
     fprintf(fp,"\t%d\t",ss->table_row[ii]+1) ; mvpa_tsv_text(fp,ss->sample_label[ii]) ;
     fputc('\t',fp) ; mvpa_tsv_text(fp,mvpa_level(ss->subject_label,ss->subject[ii])) ;
     fputc('\t',fp) ; mvpa_tsv_text(fp,mvpa_level(ss->run_label,ss->run[ii])) ;
     fputc('\t',fp) ; mvpa_tsv_text(fp,mvpa_level(ss->item_label,ss->item[ii])) ;
     fputc('\t',fp) ; mvpa_tsv_text(fp,mvpa_level(ss->fold_label,ss->fold[ii])) ;
     fputc('\t',fp) ; mvpa_tsv_text(fp,(ff>=0)?plan->name[ff]:"NA") ;
     fputc('\t',fp) ;
     mvpa_tsv_text(fp,(ff>=0)?mvpa_level(plan->stratum_label,plan->fold_stratum[ff]):"NA") ;
     fputc('\t',fp) ; mvpa_tsv_text(fp,class_label[class_id[ii]]) ; fputc('\t',fp) ;
     if( pred>=0 ){
       mvpa_tsv_text(fp,class_label[pred]) ; fprintf(fp,"\t%d",pred==class_id[ii]) ;
     } else fputs("NA\tNA",fp) ;
     for( cc=0 ; cc<workspace->nclass ; cc++ ){
       double dd=workspace->decision[(size_t)ii*workspace->nclass+cc] ;
       if( isfinite(dd) ) fprintf(fp,"\t%.17g",dd) ; else fputs("\tNA",fp) ;
     }
     fputc('\n',fp) ;
   }
}

static void mvpa_write_fold_rows( FILE *fp, THD_3dim_dataset *mset,
                                  const THD_roilist *rl, int rr,
                                  int searchlight, int surface,
                                  const THD_foldplan *plan,
                                  const THD_mvpa_workspace *workspace )
{
   int ff ;
   for( ff=0 ; ff<plan->nfold ; ff++ ){
     mvpa_write_spatial_prefix(fp,mset,rl,rr,searchlight,surface) ; fputc('\t',fp) ;
     mvpa_tsv_text(fp,plan->name[ff]) ; fputc('\t',fp) ;
     mvpa_tsv_text(fp,mvpa_level(plan->stratum_label,plan->fold_stratum[ff])) ;
     fprintf(fp,"\t%d\t%d\t%.17g\t%.17g",plan->ntrain[ff],
             workspace->fold_nsample[ff],workspace->fold_accuracy[ff],
             workspace->fold_balanced_accuracy[ff]) ;
     if(workspace->ridge_cv){
       const THD_mvpa_ridge_workspace *rw=workspace->ridge_cv;double loss=NAN;int l;
       if(rw->nested)for(l=0;l<rw->nlambda;l++)if(rw->lambda_grid[l]==rw->selected_lambda[ff])
         loss=rw->inner_loss[(size_t)ff*rw->nlambda+l];
       fprintf(fp,"\t%.17g\t%d\t%d\t%.17g",rw->selected_lambda[ff],
               rw->nested?rw->nested->inner[ff]->nfold:0,rw->inner_nvalid[ff],loss);
     }
     if(workspace->svm){
       double loss=NAN;int l;
       if(workspace->svm_nested)for(l=0;l<workspace->svm_ngrid;l++)
         if(workspace->svm_grid[l]==workspace->svm_selected[ff])loss=workspace->svm_loss[(size_t)ff*workspace->svm_ngrid+l];
       fprintf(fp,"\t%.17g\t%d\t%.17g\t%d\t%.17g",workspace->svm_selected[ff],
         workspace->svm_nested?workspace->svm_nested->inner[ff]->nfold:0,loss,
         workspace->fold_svm_iterations[ff],workspace->fold_svm_residual[ff]);
     }
     if(workspace->lda){
       int form=workspace->fold_covariance_form[ff];
       fprintf(fp,"\t%.17g\t%.17g\t%.17g\t%s",workspace->lda_shrinkage,
               workspace->fold_covariance_scale[ff],workspace->fold_variance_floor[ff],
               form==THD_LDA_FORM_DIAGONAL?"diagonal":form==THD_LDA_FORM_PRIMAL?"primal":"dual");
     }
     fputc('\n',fp);
   }
}

static double mvpa_cont_value(const THD_mvpa_continuous_score *s,int metric)
{
   switch(metric){
      case THD_METRIC_PEARSON:return s->pearson;
      case THD_METRIC_SPEARMAN:return s->spearman;
      case THD_METRIC_CROSS_VALIDATED_R2:return s->cv_r2;
      case THD_METRIC_MAE:return s->mae;
      default:return s->rmse;
   }
}

static float *mvpa_numeric_columns(THD_sampleset *ss,const char *list,
                                    char ***names,int *count)
{
   char *copy,*cur,*next;int n=1,j,i,k;float *v;
   if(!list || !*list)ERROR_exit(PROGRAM_NAME ": numeric column list is empty");
   for(cur=(char *)list;*cur;cur++)if(*cur==',')n++;
   *names=calloc((size_t)n,sizeof(char *));
   v=malloc(sizeof(float)*(size_t)ss->nsample*n);copy=strdup(list);
   if(!*names || !v || !copy)ERROR_exit(PROGRAM_NAME ": cannot allocate numeric columns");
   cur=copy;
   for(j=0;j<n;j++){
      next=strchr(cur,',');if(next)*next++='\0';
      k=THD_datatable_column(ss->metadata,cur);
      if(k<0 || !ss->metadata->isnum[k])
         ERROR_exit(PROGRAM_NAME ": column '%s' must contain finite numeric values",cur);
      for(i=0;i<j;i++)if(!strcasecmp((*names)[i],cur))
         ERROR_exit(PROGRAM_NAME ": duplicate numeric column '%s'",cur);
      (*names)[j]=strdup(cur);
      for(i=0;i<ss->nsample;i++){
         float a=ss->metadata->val[k][ss->table_row[i]];
         if(!isfinite(a))ERROR_exit(PROGRAM_NAME ": column '%s', sample %d is not finite",cur,i+1);
         v[(size_t)i*n+j]=a;
      }cur=next;
   }
   free(copy);*count=n;return v;
}

static THD_roilist *mvpa_voxel_units(THD_3dim_dataset *mask)
{
   THD_roilist *atlas=THD_roilist_from_dset(mask,NULL),*units;
   int r,j,k=0,n=0;
   if(!atlas)ERROR_exit(PROGRAM_NAME ": encoding mask has no positive voxels/nodes");
   for(r=0;r<atlas->nroi;r++)n+=atlas->vox[r].nar;
   units=calloc(1,sizeof(*units));if(!units)ERROR_exit(PROGRAM_NAME ": cannot allocate voxel units");
   units->nroi=n;units->val=malloc(sizeof(int)*(size_t)n);
   units->center=malloc(sizeof(int)*(size_t)n);
   units->vox=calloc((size_t)n,sizeof(intvec));units->lab=calloc((size_t)n,sizeof(char *));
   if(!units->val || !units->center || !units->vox || !units->lab)
      ERROR_exit(PROGRAM_NAME ": cannot allocate voxel unit buffers");
   for(r=0;r<atlas->nroi;r++)for(j=0;j<atlas->vox[r].nar;j++,k++){
      int v=atlas->vox[r].ar[j];
      units->center[k]=v;units->val[k]=mask->dblk->node_list?mask->dblk->node_list[v]:v;
      units->vox[k].nar=1;units->vox[k].ar=malloc(sizeof(int));
      if(!units->vox[k].ar)ERROR_exit(PROGRAM_NAME ": cannot allocate voxel index");
      units->vox[k].ar[0]=v;
   }
   THD_roilist_free(atlas);return units;
}

static double *mvpa_penalty_grid(const char *text,int *count,const char *what)
{
   int n=1,i,j;char *copy=strdup(text),*cur,*next,*end;double *grid;
   if(!copy)ERROR_exit(PROGRAM_NAME ": cannot allocate penalty grid");
   for(cur=copy;*cur;cur++)if(*cur==',')n++;
   if(n<2)ERROR_exit(PROGRAM_NAME ": lambda grid needs at least two candidates");
   grid=malloc(sizeof(double)*(size_t)n);if(!grid)ERROR_exit(PROGRAM_NAME ": cannot allocate penalty candidates");
   cur=copy;
   for(i=0;i<n;i++){
      next=strchr(cur,',');if(next)*next++='\0';errno=0;grid[i]=strtod(cur,&end);
      if(errno || cur==end || *end || !isfinite(grid[i]) || grid[i]<=0)
         ERROR_exit(PROGRAM_NAME ": lambda grid values must be finite and positive");
      for(j=0;j<i;j++)if(grid[i]==grid[j])ERROR_exit(PROGRAM_NAME ": lambda grid contains duplicate penalties");
      cur=next;
   }
   free(copy);*count=n;return grid;
}

static void mvpa_continuous_analysis(THD_sampleset *ss,const THD_foldplan *plan,
 THD_3dim_dataset *mset,THD_roilist *rl,int encoding,const char *columns,
 int preprocess,const int *groups,double lambda,int metric,const char *prefix,
 const char *predfile,const char *foldfile,int nworker,int quiet,
 const THD_nested_foldplan *nested,const double *grid,int nl,const char *inner_column,
 const PERM_set *pset,const char *permute_within,
 int argc,char **argv)
{
   char **names=NULL,err[1024],label[256],*tablefile,runerr[1024]="";
   float *attr,**patterns,*values,*paint;int na,maxvox,r,o,k,b,nout,nbrick,i,failed=0;
   THD_mvpa_ridge_workspace **workers;THD_3dim_dataset *out;
   FILE *table,*pred=NULL,*fold=NULL,*tuning=NULL;char *tunefile=NULL;
   PERM_result **inference=NULL;
   float **permuted=NULL,**null_score=NULL,**max_null=NULL,*qval=NULL,*null_mean=NULL,*null_valid=NULL;
   unsigned char *valid=NULL;
   int np=pset?pset->nperm:0;double direction=metric==THD_METRIC_MAE || metric==THD_METRIC_RMSE?-1.0:1.0;
   const char *mlabel[]={"PEARSON","SPEARMAN","CVR2","MAE","RMSE"};
   const char *ilabel[]={"NULL_MEAN","ABOVE_NULL","PVAL","ZSCR","QVAL","PFWE","ZFWE","NNULL_VALID"};
   attr=mvpa_numeric_columns(ss,columns,&names,&na);maxvox=THD_roilist_maxvox(rl);
   nout=encoding?1:na;
   if(nout>(32760-2)/(pset?15:7))ERROR_exit(PROGRAM_NAME ": too many output targets");
   nbrick=nout*(pset?15:7)+2;
   workers=calloc((size_t)nworker,sizeof(*workers));patterns=calloc((size_t)nworker,sizeof(*patterns));
   if(!workers || !patterns)ERROR_exit(PROGRAM_NAME ": cannot allocate ridge workers");
   if(pset){
      if((size_t)np>SIZE_MAX/sizeof(float)/(size_t)nout)
         ERROR_exit(PROGRAM_NAME ": continuous null dimensions too large");
      inference=calloc((size_t)nout,sizeof(*inference));
      permuted=calloc((size_t)nworker,sizeof(*permuted));
      null_score=calloc((size_t)nworker,sizeof(*null_score));max_null=calloc((size_t)nworker,sizeof(*max_null));
      valid=calloc((size_t)nout*rl->nroi,1);
      qval=malloc(sizeof(float)*(size_t)nout*rl->nroi);
      null_mean=malloc(sizeof(float)*(size_t)nout*rl->nroi);
      null_valid=calloc((size_t)nout*rl->nroi,sizeof(float));
      if(!inference || !permuted || !null_score || !max_null || !valid || !qval || !null_mean || !null_valid)
         ERROR_exit(PROGRAM_NAME ": cannot allocate continuous inference");
      for(o=0;o<nout;o++){
         inference[o]=THD_perm_result_new(rl->nroi,np);
         if(!inference[o])ERROR_exit(PROGRAM_NAME ": cannot allocate continuous inference family");
         inference[o]->tail=PERM_TAIL_ONE;
         for(i=0;i<np;i++)inference[o]->max_null[i]=-FLT_MAX;
      }
   }
   for(i=0;i<nworker;i++){
      workers[i]=THD_mvpa_ridge_workspace_alloc(plan,encoding?na:maxvox,encoding?maxvox:na,
                                               preprocess,err,sizeof(err));
      patterns[i]=malloc(sizeof(float)*(size_t)ss->nsample*maxvox);
      if(!workers[i] || !patterns[i])ERROR_exit(PROGRAM_NAME ": cannot allocate ridge worker: %s",err);
      if(nested && !THD_mvpa_ridge_set_tuning(workers[i],nested,grid,nl,err,sizeof(err)))
         ERROR_exit(PROGRAM_NAME ": %s",err);
      if(pset){
         size_t j;
         permuted[i]=malloc(sizeof(float)*(size_t)ss->nsample*na);
         null_score[i]=malloc(sizeof(float)*(size_t)np*nout);
         max_null[i]=malloc(sizeof(float)*(size_t)np*nout);
         if(!permuted[i] || !null_score[i] || !max_null[i])ERROR_exit(PROGRAM_NAME ": cannot allocate continuous null worker");
         for(j=0;j<(size_t)np*nout;j++)max_null[i][j]=-FLT_MAX;
      }
   }
   if((size_t)rl->nroi>SIZE_MAX/sizeof(float)/(size_t)nbrick)
      ERROR_exit(PROGRAM_NAME ": continuous output dimensions too large");
   values=calloc((size_t)nbrick*rl->nroi,sizeof(float));
   paint=malloc(sizeof(float)*(size_t)rl->nroi);
   tablefile=malloc(strlen(prefix)+20);
   if(!values || !paint || !tablefile)ERROR_exit(PROGRAM_NAME ": out of memory in continuous analysis");
   sprintf(tablefile,"%s.continuous.tsv",prefix);
   if((predfile && !strcmp(predfile,tablefile)) || (foldfile && !strcmp(foldfile,tablefile)))
      ERROR_exit(PROGRAM_NAME ": continuous result and audit filenames must differ");
   table=mvpa_open_audit(tablefile,"continuous results");
   fprintf(table,"# 3dMVPA %s; lambda=%s; scale=%s; primary_metric=%s\n",
           encoding?"encode":"decode",nested?"nested":"fixed",THD_preprocess_name(preprocess),THD_mvpa_metric_name(metric));
   if(!nested)fprintf(table,"# fixed lambda=%.17g\n",lambda);
   fprintf(table,"# Pattern summaries: arithmetic mean over defined per-response metrics; NaN when none are defined.\n");
   fprintf(table,"# Undefined metrics are NaN here and zero in maps; NVALID counts defined primary-metric responses.\n");
   if(pset)fprintf(table,"# Inference in PREFIX.inference.tsv: %d identity-inclusive row permutations within %s; one spatial family per decoding target (one for encoding).\n",np,permute_within);
   fprintf(table,"Unit\tResponse\tNResponse\tNValid\tSCORE\tPEARSON\tSPEARMAN\tCVR2\tMAE\tRMSE\n");
   if(predfile){
      pred=mvpa_open_audit(predfile,"continuous predictions");
      fprintf(pred,"Unit\tResponse\tTableRow\tSample\tSubject\tRun\tCVFold\tObserved\tPredicted\tTrainingMean\tItem\tInputFold\tCVStratum\n");
   }
   if(foldfile){
      fold=mvpa_open_audit(foldfile,"continuous fold scores");
      fprintf(fold,"Unit\tResponse\tCVFold\tCVStratum\tNTrain\tNTest\tPEARSON\tSPEARMAN\tCVR2\tMAE\tRMSE\tLambda\tInnerNFold\tInnerNValid\tInnerLoss\n");
   }
   if(nested){
      tunefile=malloc(strlen(prefix)+20);if(!tunefile)ERROR_exit(PROGRAM_NAME ": cannot name tuning audit");
      sprintf(tunefile,"%s.tuning.tsv",prefix);
      if((predfile && !strcmp(predfile,tunefile)) || (foldfile && !strcmp(foldfile,tunefile)))
         ERROR_exit(PROGRAM_NAME ": tuning result and prediction/fold filenames must differ");
      tuning=mvpa_open_audit(tunefile,"lambda tuning audit");
      fprintf(tuning,"# inner_group=%s; objective=mean_response_SSE_over_training_mean_baseline_SSE; ties_or_no_valid_responses=largest_lambda\n",inner_column);
      fprintf(tuning,"Unit\tCVFold\tCVStratum\tInnerNFold\tInnerNValid\tLambda\tInnerLoss\tSelected\n");
   }
   AFNI_OMP_START;
#ifdef USE_OMP
#pragma omp parallel for num_threads(nworker) schedule(dynamic,1) ordered
#endif
   for(r=0;r<rl->nroi;r++){
      int q=encoding?rl->vox[r].nar:na,tid=0,ok,o,k,b,i,f;
      char localerr[1024];float *neural;THD_mvpa_ridge_workspace *w;
      const float *xx,*yy;
#ifdef USE_OMP
      tid=omp_get_thread_num();
#endif
      neural=patterns[tid];w=workers[tid];
      ok=THD_roi_samples(ss,NULL,ss->nsample,rl->vox+r,neural,localerr,sizeof(localerr));
      xx=encoding?attr:neural;yy=encoding?neural:attr;
      if(ok && pset)ok=THD_mvpa_run_ridge_permutations(w,xx,encoding?na:rl->vox[r].nar,yy,q,
         groups,lambda,pset,encoding,metric,encoding,permuted[tid],null_score[tid],localerr,sizeof(localerr));
      else if(ok)ok=THD_mvpa_run_ridge_cv(w,xx,encoding?na:rl->vox[r].nar,yy,q,groups,lambda,localerr,sizeof(localerr));
#ifdef USE_OMP
#pragma omp ordered
#endif
      {
      if(!ok){
         if(!failed){failed=1;snprintf(runerr,sizeof(runerr),"unit %d: %.900s",rl->val[r],localerr);}
      }else{
      if(pset)for(o=0;o<nout;o++){
         PERM_result *pr=inference[o];double mean=0;int ip,count=0;
         float obs=null_score[tid][o];
         valid[(size_t)o*rl->nroi+r]=isfinite(obs);
         pr->stat[r]=isfinite(obs)?(float)(direction*obs):0;
         for(ip=0;ip<np;ip++){
            float v=null_score[tid][(size_t)ip*nout+o];
            float directed=isfinite(v)?(float)(direction*v):isfinite(obs)?FLT_MAX:-FLT_MAX;
            if(isfinite(v)){mean+=v;count++;}
            if(isfinite(obs) && (!isfinite(v) || directed>=pr->stat[r]))pr->cnt_unc[r]++;
            if(directed>max_null[tid][(size_t)o*np+ip])max_null[tid][(size_t)o*np+ip]=directed;
         }
         null_mean[(size_t)o*rl->nroi+r]=count?(float)(mean/count):NAN;
         null_valid[(size_t)o*rl->nroi+r]=(float)count;
         pr->effect[r]=isfinite(obs) && count?(float)(direction*(obs-mean/count)):0;
      }
      if(tuning)for(f=0;f<plan->nfold;f++)for(i=0;i<nl;i++){
         fprintf(tuning,"%d\t",rl->val[r]);mvpa_tsv_text(tuning,plan->name[f]);
         fputc('\t',tuning);mvpa_tsv_text(tuning,plan->stratum_label[plan->fold_stratum[f]]);
         fprintf(tuning,"\t%d\t%d\t%.17g\t%.17g\t%d\n",nested->inner[f]->nfold,
                 w->inner_nvalid[f],grid[i],w->inner_loss[(size_t)f*nl+i],grid[i]==w->selected_lambda[f]);
      }
      for(o=0;o<nout;o++){
         double sum[6]={0,0,0,0,0,0};int count[6]={0,0,0,0,0,0},j;
         for(j=encoding?0:o;j<(encoding?q:o+1);j++)for(k=0;k<6;k++){
            double v=mvpa_cont_value(w->score+j,k==0?metric:THD_METRIC_PEARSON+k-1);
            if(isfinite(v)){sum[k]+=v;count[k]++;}
         }
         b=o*7;
         for(k=0;k<6;k++)values[(size_t)(b+k)*rl->nroi+r]=count[k]?(float)(sum[k]/count[k]):NAN;
         values[(size_t)(b+6)*rl->nroi+r]=(float)count[0];
         fprintf(table,"%d\t",rl->val[r]);mvpa_tsv_text(table,encoding?"pattern":names[o]);
         fprintf(table,"\t%d\t%d",encoding?q:1,count[0]);
         for(k=0;k<6;k++)fprintf(table,"\t%.9g",values[(size_t)(b+k)*rl->nroi+r]);
         fputc('\n',table);
      }
      values[(size_t)(nbrick-2)*rl->nroi+r]=(float)rl->vox[r].nar;
      values[(size_t)(nbrick-1)*rl->nroi+r]=(float)plan->nfold;
      for(o=0;o<q;o++){
         char response[256];
         if(encoding){
            int v=rl->vox[r].ar[o];
            int surface=DBLK_IS_NI_SURF_DSET(mset->dblk) || DBLK_IS_GIFTI(mset->dblk);
            int id=surface && mset->dblk->node_list?mset->dblk->node_list[v]:v;
            snprintf(response,sizeof(response),"%s_%d",surface?"node":"voxel",id);
         }
         else snprintf(response,sizeof(response),"%s",names[o]);
         if(pred)for(i=0;i<ss->nsample;i++){
            f=plan->sample_test_fold[i];
            fprintf(pred,"%d\t",rl->val[r]);mvpa_tsv_text(pred,response);
            fprintf(pred,"\t%d\t",ss->table_row[i]+1);mvpa_tsv_text(pred,ss->sample_label[i]);
            fputc('\t',pred);mvpa_tsv_text(pred,ss->subject_label[ss->subject[i]]);
            fputc('\t',pred);mvpa_tsv_text(pred,ss->run[i]>=0?ss->run_label[ss->run[i]]:NULL);
            fputc('\t',pred);mvpa_tsv_text(pred,f>=0?plan->name[f]:NULL);
            fprintf(pred,"\t%.9g\t%.9g\t%.9g\t",yy[(size_t)i*q+o],w->prediction[(size_t)i*q+o],w->baseline[(size_t)i*q+o]);
            mvpa_tsv_text(pred,ss->item[i]>=0?ss->item_label[ss->item[i]]:NULL);
            fputc('\t',pred);mvpa_tsv_text(pred,ss->fold[i]>=0?ss->fold_label[ss->fold[i]]:NULL);
            fputc('\t',pred);mvpa_tsv_text(pred,f>=0?plan->stratum_label[plan->fold_stratum[f]]:NULL);
            fputc('\n',pred);
         }
         if(fold)for(f=0;f<plan->nfold;f++){
            fprintf(fold,"%d\t",rl->val[r]);mvpa_tsv_text(fold,response);
            fputc('\t',fold);mvpa_tsv_text(fold,plan->name[f]);
            fputc('\t',fold);mvpa_tsv_text(fold,plan->stratum_label[plan->fold_stratum[f]]);
            fprintf(fold,"\t%d\t%d",plan->ntrain[f],plan->ntest[f]);
            for(k=0;k<5;k++)fprintf(fold,"\t%.9g",mvpa_cont_value(w->fold_score+(size_t)f*q+o,THD_METRIC_PEARSON+k));
            {
               double loss=NAN;int l;
               if(nested)for(l=0;l<nl;l++)if(grid[l]==w->selected_lambda[f])loss=w->inner_loss[(size_t)f*nl+l];
               fprintf(fold,"\t%.17g\t%d\t%d\t%.17g",w->selected_lambda[f],
                       nested?nested->inner[f]->nfold:0,w->inner_nvalid[f],loss);
            }
            fputc('\n',fold);
         }
      }
      }
      }
   }
   AFNI_OMP_END;
   if(ferror(table) || fclose(table))ERROR_exit(PROGRAM_NAME ": error writing continuous results");
   if(pred && (ferror(pred) || fclose(pred)))ERROR_exit(PROGRAM_NAME ": error writing continuous predictions");
   if(fold && (ferror(fold) || fclose(fold)))ERROR_exit(PROGRAM_NAME ": error writing continuous folds");
   if(tuning && (ferror(tuning) || fclose(tuning)))ERROR_exit(PROGRAM_NAME ": error writing tuning audit");
   if(failed)ERROR_exit(PROGRAM_NAME ": %s",runerr);
   if(pset){
      FILE *fp;char *fname=malloc(strlen(prefix)+20);int ww,ip;
      if(!fname)ERROR_exit(PROGRAM_NAME ": cannot name inference audit");
      sprintf(fname,"%s.inference.tsv",prefix);
      if((predfile && !strcmp(predfile,fname)) || (foldfile && !strcmp(foldfile,fname)))
         ERROR_exit(PROGRAM_NAME ": inference audit and prediction/fold filenames must differ");
      fp=mvpa_open_audit(fname,"continuous inference audit");
      fprintf(fp,"# null=%s_rows_within_%s; nperm=%d; exact=%d; tail=%s; correction_family=spatial_per_response\n",
              encoding?"design":"target",permute_within,np,pset->is_exact,direction>0?"upper":"lower");
      fprintf(fp,"# All preprocessing/tuning/models refitted per relabeling. Undefined observed metrics have p=q=pfwe=1 and z=0.\n");
      fprintf(fp,"# Undefined null metrics count as exceedances for defined observations and conservatively enter that slot's family maximum.\n");
      fprintf(fp,"Unit\tResponse\tSCORE\tNULL_MEAN\tABOVE_NULL\tPVAL\tZSCR\tQVAL\tPFWE\tZFWE\tNNULL_VALID\tNPERM\n");
      for(o=0;o<nout;o++){
         PERM_result *pr=inference[o];
         for(ww=0;ww<nworker;ww++)THD_max_accum(np,pr->max_null,max_null[ww]+(size_t)o*np);
         for(ip=0;ip<np;ip++)if(pr->max_null[ip]==-FLT_MAX)pr->max_null[ip]=0;
         THD_perm_result_finish(pr,valid+(size_t)o*rl->nroi);
         THD_bh_fdr_masked(rl->nroi,pr->p_unc,NULL,qval+(size_t)o*rl->nroi);
         for(r=0;r<rl->nroi;r++){
            float vv[8]={null_mean[(size_t)o*rl->nroi+r],pr->effect[r],pr->p_unc[r],pr->z_unc[r],
                         qval[(size_t)o*rl->nroi+r],pr->p_fwe[r],pr->z_fwe[r],null_valid[(size_t)o*rl->nroi+r]};
            for(k=0;k<8;k++)values[(size_t)(nout*7+o*8+k)*rl->nroi+r]=vv[k];
            fprintf(fp,"%d\t",rl->val[r]);mvpa_tsv_text(fp,encoding?"pattern":names[o]);
            fprintf(fp,"\t%.9g",values[(size_t)(o*7)*rl->nroi+r]);
            for(k=0;k<8;k++)fprintf(fp,"\t%.9g",vv[k]);
            fprintf(fp,"\t%d\n",np);
         }
      }
      if(ferror(fp) || fclose(fp))ERROR_exit(PROGRAM_NAME ": error writing inference audit");
      free(fname);
      if(!quiet)INFO_message(PROGRAM_NAME ": %d continuous inference famil%s, %d spatial units each, %d synchronized row relabelings",
                            nout,nout==1?"y":"ies",rl->nroi,np);
   }
   out=EDIT_empty_copy(mset);
   if(!out)ERROR_exit(PROGRAM_NAME ": cannot allocate continuous map");
   EDIT_dset_items(out,ADN_prefix,prefix,ADN_nvals,nbrick,ADN_ntt,0,
                  ADN_type,HEAD_FUNC_TYPE,ADN_func_type,FUNC_BUCK_TYPE,ADN_none);
   if(THD_deathcon() && THD_is_file(DSET_HEADNAME(out)))ERROR_exit(PROGRAM_NAME ": output dataset already exists");
   for(b=0;b<nbrick;b++){
      if(b>=nbrick-2)snprintf(label,sizeof(label),"%s",b==nbrick-2?"NFEAT":"NFOLD");
      else if(b<nout*7){
         o=b/7;k=b%7;
         snprintf(label,sizeof(label),"%s_%s",k==0?"SCORE":k==6?"NVALID":mlabel[k-1],encoding?"neural":names[o]);
      }
      else{
         o=(b-nout*7)/8;k=(b-nout*7)%8;
         snprintf(label,sizeof(label),"%s_%s",ilabel[k],encoding?"neural":names[o]);
      }
      EDIT_substitute_brick(out,b,MRI_float,NULL);
      for(r=0;r<rl->nroi;r++)paint[r]=isfinite(values[(size_t)b*rl->nroi+r])?values[(size_t)b*rl->nroi+r]:0.0f;
      THD_roilist_paint(DSET_ARRAY(out,b),rl,paint);
      EDIT_BRICK_LABEL(out,b,label);EDIT_BRICK_TO_NOSTAT(out,b);EDIT_BRICK_FACTOR(out,b,0.0f);
      if(pset && b>=nout*7 && b<nbrick-2 && ((b-nout*7)%8==3 || (b-nout*7)%8==6))EDIT_BRICK_TO_FIZT(out,b);
   }
   THD_copy_labeltable_atr(out->dblk,mset->dblk);
   if(!THD_copy_surface_node_domain(out,mset,err,sizeof(err)))ERROR_exit(PROGRAM_NAME ": %s",err);
   tross_Copy_History(mset,out);tross_Make_History(PROGRAM_NAME,argc,argv,out);
   if(DSET_write(out)!=True)ERROR_exit(PROGRAM_NAME ": failed writing continuous map");
   if(!quiet)INFO_message(PROGRAM_NAME ": ridge %s: %d units, %d folds; wrote %s",encoding?"encoding":"decoding",rl->nroi,plan->nfold,tablefile);
   DSET_delete(out);
   for(i=0;i<nworker;i++){THD_mvpa_ridge_workspace_free(workers[i]);free(patterns[i]);}
   if(pset){
      for(i=0;i<nworker;i++){free(permuted[i]);free(null_score[i]);free(max_null[i]);}
      for(o=0;o<nout;o++)THD_perm_result_free(inference[o]);
      free(inference);free(permuted);free(null_score);free(max_null);
      free(valid);free(qval);free(null_mean);free(null_valid);
   }
   free(workers);free(patterns);
   free(tunefile);
   THD_sampleset_free_levels(names,na);free(attr);free(values);free(paint);free(tablefile);
}

/* Descriptive refits are an independent final pass. They cannot replace CV
   predictions, tuning diagnostics or permutation statistics. Neural axes are
   mapped at their actual feature/response rows, never at searchlight centers. */
static THD_3dim_dataset *mvpa_coefficient_map(THD_3dim_dataset *mset,
 const char *prefix,int nb,int argc,char **argv)
{
   int b;char err[512];THD_3dim_dataset *out=EDIT_empty_copy(mset);
   if(!out || nb<1 || nb>32760)ERROR_exit(PROGRAM_NAME ": invalid coefficient map dimensions");
   EDIT_dset_items(out,ADN_prefix,prefix,ADN_nvals,nb,ADN_ntt,0,
                  ADN_type,HEAD_FUNC_TYPE,ADN_func_type,FUNC_BUCK_TYPE,ADN_none);
   if(THD_deathcon() && THD_is_file(DSET_HEADNAME(out)))
      ERROR_exit(PROGRAM_NAME ": coefficient dataset '%s' already exists",prefix);
   for(b=0;b<nb;b++){
      float *v=calloc((size_t)DSET_NVOX(mset),sizeof(float));
      if(!v)ERROR_exit(PROGRAM_NAME ": cannot allocate coefficient map");
      EDIT_substitute_brick(out,b,MRI_float,v);
      EDIT_BRICK_TO_NOSTAT(out,b);EDIT_BRICK_FACTOR(out,b,0.0f);
   }
   if(!THD_copy_surface_node_domain(out,mset,err,sizeof(err)))ERROR_exit(PROGRAM_NAME ": %s",err);
   tross_Copy_History(mset,out);tross_Make_History(PROGRAM_NAME,argc,argv,out);
   return out;
}

static void mvpa_neural_id(char *text,int size,THD_3dim_dataset *mset,int row)
{
   int surface=DBLK_IS_NI_SURF_DSET(mset->dblk) || DBLK_IS_GIFTI(mset->dblk);
   int id=surface && mset->dblk->node_list?mset->dblk->node_list[row]:row;
   snprintf(text,size,"%s_%d",surface?"node":"voxel",id);
}

static float mvpa_coefficient_float(double v)
{
   if(!isfinite(v) || fabs(v)>FLT_MAX)
      ERROR_exit(PROGRAM_NAME ": descriptive coefficient exceeds map float range");
   return (float)v;
}

static void mvpa_coefficient_prefix(char *name,size_t capacity,
 THD_3dim_dataset *mset,const char *prefix,int combined,int id)
{
   int surface=DBLK_IS_NI_SURF_DSET(mset->dblk) || DBLK_IS_GIFTI(mset->dblk);
   size_t stemlen=strlen(prefix);
   if(surface && stemlen>=10 && !strcmp(prefix+stemlen-10,".niml.dset"))stemlen-=10;
   if(combined)snprintf(name,capacity,"%.*s%s",(int)stemlen,prefix,surface?".niml.dset":"");
   else snprintf(name,capacity,"%.*s.unit_%d%s",(int)stemlen,prefix,id,surface?".niml.dset":"");
}

static void mvpa_descriptive_coefficients(THD_sampleset *ss,THD_3dim_dataset *mset,
 THD_roilist *rl,int encoding,int searchlight,const char *columns,
 const int *classes,char **class_names,int nc,int preprocess,const int *groups,
 double lambda,const char *prefix,int quiet,int argc,char **argv)
{
   char **names=NULL,err[1024],feature[256],response[256],label[256];
   char *fname=malloc(strlen(prefix)+32),*mapname=malloc(strlen(prefix)+64);
   float *attr,*neural,*transformed;double *coef,*intercept;
   int na,n=ss->nsample,maxvox=THD_roilist_maxvox(rl),maxp,maxq,*rows;
   int r,j,o,g,i,combined=encoding && !searchlight;
   FILE *table,*prepfile;THD_ridge_workspace *ridge;THD_3dim_dataset *out=NULL;
   if(!fname || !mapname)ERROR_exit(PROGRAM_NAME ": cannot name descriptive outputs");
   if(classes){
      na=nc;names=class_names;attr=calloc((size_t)n*nc,sizeof(float));
      if(!attr)ERROR_exit(PROGRAM_NAME ": cannot allocate descriptive class indicators");
      for(i=0;i<n;i++)attr[(size_t)i*nc+classes[i]]=1;
   }else attr=mvpa_numeric_columns(ss,columns,&names,&na);
   maxp=encoding?na:maxvox;maxq=encoding?maxvox:na;
   if((size_t)maxp>SIZE_MAX/sizeof(double)/(size_t)maxq ||
      (size_t)n>SIZE_MAX/sizeof(float)/(size_t)maxp)
      ERROR_exit(PROGRAM_NAME ": descriptive refit dimensions too large");
   rows=malloc(sizeof(int)*(size_t)n);neural=malloc(sizeof(float)*(size_t)n*maxvox);
   transformed=malloc(sizeof(float)*(size_t)n*maxp);
   coef=malloc(sizeof(double)*(size_t)maxp*maxq);intercept=malloc(sizeof(double)*(size_t)maxq);
   ridge=THD_ridge_alloc(n,maxp,maxq,err,sizeof(err));
   if(!rows || !neural || !transformed || !coef || !intercept || !ridge)
      ERROR_exit(PROGRAM_NAME ": cannot allocate descriptive refit");
   for(i=0;i<n;i++)rows[i]=i;
   sprintf(fname,"%s.coefficients.tsv",prefix);table=mvpa_open_audit(fname,"descriptive coefficients");
   fprintf(table,"# Descriptive full-data ridge refit; NOT held-out scores, activation/Haufe maps, or coefficient significance.\n");
   fprintf(table,"# scale=%s; lambda=%.17g; nfit=%d; objective=sum_squared_error+lambda*sum_squared_slopes; intercept_unpenalized\n",
           THD_preprocess_name(preprocess),lambda,n);
   fprintf(table,"# y = Intercept + sum(preprocessed_feature * Coefficient); responses retain original units%s.\n",
           classes?" (classification: 0/1 indicators; decisions are not probabilities)":"");
   fprintf(table,"# Preprocessing fitted on ALL samples; see PREFIX.preprocessing.tsv. Weights use transformed predictor units.\n");
   fprintf(table,"# Map bricks are zero-based; outside-unit values are zero. Overlapping fits have separate unit maps.\n");
   fprintf(table,"Unit\tFeature\tResponse\tCoefficient\tIntercept\tLambda\tNFit\tMapPrefix\tMapBrick\n");
   sprintf(fname,"%s.preprocessing.tsv",prefix);prepfile=mvpa_open_audit(fname,"descriptive preprocessing");
   fprintf(prepfile,"# scale=%s; feature_zscore: (x-center)/scale (population SD); constant features become zero.\n",THD_preprocess_name(preprocess));
   fprintf(prepfile,"# none: x unchanged; sample_center: subtract each sample's own mean across this unit's predictors.\n");
   fprintf(prepfile,"# run_center: subtract full-data subject/run feature means; Group is SubjectRunID below.\n");
   if(groups)for(i=0;i<n;i++){
      fprintf(prepfile,"# SampleGroup\t%d\t%d\t",ss->table_row[i]+1,groups[i]);
      mvpa_tsv_text(prepfile,ss->subject_label[ss->subject[i]]);fputc('\t',prepfile);
      mvpa_tsv_text(prepfile,ss->run_label[ss->run[i]]);fputc('\n',prepfile);
   }
   fprintf(prepfile,"Unit\tFeature\tGroup\tCenter\tScale\tConstant\n");
   if(combined){
      mvpa_coefficient_prefix(mapname,strlen(prefix)+64,mset,prefix,1,0);
      out=mvpa_coefficient_map(mset,mapname,na+1,argc,argv);
   }
   for(r=0;r<rl->nroi;r++){
      int p=encoding?na:rl->vox[r].nar,q=encoding?rl->vox[r].nar:na;
      const float *xx,*yy;THD_preprocessor *prep;
      if(!THD_roi_samples(ss,NULL,n,rl->vox+r,neural,err,sizeof(err)))ERROR_exit(PROGRAM_NAME ": %s",err);
      xx=encoding?attr:neural;yy=encoding?neural:attr;
      prep=THD_preprocess_fit(preprocess,xx,n,p,rows,n,groups,err,sizeof(err));
      if(!prep || !THD_preprocess_apply(prep,xx,n,p,rows,n,groups,transformed,err,sizeof(err)) ||
         !THD_ridge_fit(ridge,transformed,yy,n,p,q,lambda,THD_RIDGE_AUTO,err,sizeof(err)) ||
         !THD_ridge_coefficients(ridge,coef,intercept,err,sizeof(err)))
         ERROR_exit(PROGRAM_NAME ": descriptive unit %d: %s",rl->val[r],err);
      if(!combined){
         mvpa_coefficient_prefix(mapname,strlen(prefix)+64,mset,prefix,0,rl->val[r]);
         out=mvpa_coefficient_map(mset,mapname,encoding?na+1:na,argc,argv);
      }
      for(j=0;j<p;j++){
         const char *feature_id;
         if(!encoding)mvpa_neural_id(feature,sizeof(feature),mset,rl->vox[r].ar[j]);
         feature_id=encoding?names[j]:feature;
         for(g=0;g<(preprocess==THD_PREP_RUN_CENTER?prep->ngroup:1);g++){
            double center=0,scale=1;int constant=0;
            if(preprocess==THD_PREP_FEATURE_ZSCORE){center=prep->center[j];scale=prep->scale[j];constant=prep->constant[j];}
            if(preprocess==THD_PREP_RUN_CENTER)center=prep->center[(size_t)g*p+j];
            fprintf(prepfile,"%d\t",rl->val[r]);mvpa_tsv_text(prepfile,feature_id);
            if(preprocess==THD_PREP_RUN_CENTER)fprintf(prepfile,"\t%d",prep->group_value[g]);
            else fputs("\tNA",prepfile);
            fprintf(prepfile,"\t%.17g\t%.17g\t%d\n",center,scale,constant);
         }
         for(o=0;o<q;o++){
            int brick=encoding?j:o,row=rl->vox[r].ar[encoding?o:j];
            const char *response_id;
            if(encoding)mvpa_neural_id(response,sizeof(response),mset,row);
            response_id=encoding?response:names[o];
            fprintf(table,"%d\t",rl->val[r]);mvpa_tsv_text(table,feature_id);fputc('\t',table);mvpa_tsv_text(table,response_id);
            fprintf(table,"\t%.17g\t%.17g\t%.17g\t%d\t",coef[(size_t)j*q+o],intercept[o],lambda,n);
            mvpa_tsv_text(table,mapname);fprintf(table,"\t%d\n",brick);
            ((float *)DSET_ARRAY(out,brick))[row]=mvpa_coefficient_float(coef[(size_t)j*q+o]);
         }
      }
      for(i=0;i<na;i++){
         snprintf(label,sizeof(label),"REFIT_B%d_%s",i,names[i]);EDIT_BRICK_LABEL(out,i,label);
      }
      if(encoding){
         EDIT_BRICK_LABEL(out,na,"REFIT_INTERCEPT");
         for(o=0;o<q;o++)((float *)DSET_ARRAY(out,na))[rl->vox[r].ar[o]]=mvpa_coefficient_float(intercept[o]);
      }
      THD_preprocess_free(prep);
      if(!combined){
         if(DSET_write(out)!=True)ERROR_exit(PROGRAM_NAME ": failed writing coefficient map");
         DSET_delete(out);out=NULL;
      }
   }
   if(combined){
      if(DSET_write(out)!=True)ERROR_exit(PROGRAM_NAME ": failed writing coefficient map");
      DSET_delete(out);
   }
   if(ferror(table) || fclose(table))ERROR_exit(PROGRAM_NAME ": error writing coefficients");
   if(ferror(prepfile) || fclose(prepfile))ERROR_exit(PROGRAM_NAME ": error writing descriptive preprocessing");
   if(!quiet)INFO_message(PROGRAM_NAME ": descriptive full-data coefficients: %d units, lambda %.9g; %s.coefficients.tsv",rl->nroi,lambda,prefix);
   THD_ridge_free(ridge);free(rows);free(neural);free(transformed);free(coef);free(intercept);
   if(!classes)THD_sampleset_free_levels(names,na);
   free(attr);free(fname);free(mapname);
}

int main( int argc, char **argv )
{
   THD_datatable *tab=NULL,*expanded=NULL ; THD_sampleset *ss=NULL ;
   THD_sampleset_columns columns ; THD_foldplan *plan=NULL ;
   THD_nested_foldplan *nested=NULL;
   THD_3dim_dataset *mset=NULL,*oset=NULL ; THD_roilist *rl=NULL ;
   THD_mvpa_workspace **workspace=NULL ;
   PERM_scheme *perm_scheme=NULL ; PERM_set *perm_set=NULL ;
   PERM_result *perm_result=NULL ;
   char *mode=NULL,*target_column=NULL,*cv_column=NULL,*within_column=NULL ;
   char *subject_column=NULL,*run_column=NULL,*item_column=NULL ;
   char *fold_column=NULL,*sample_column=NULL,*xmat_stims=NULL ;
   int trial_manifest=0 ;
   char *maskname=NULL,*roi_sel=NULL,*searchlight=NULL,*surf_file=NULL ;
   char *prefix="MVPA",*method=NULL,*design_columns=NULL,*score_mode=NULL ;
   char *coef_prefix=NULL;double coef_lambda=0;int coef_lambda_given=0;
   int svm_class=0,svm_tune=0,svm_numeric=0,svm_options=0,svm_ngrid=0,svm_max_iter=10000;
   double svm_cost=1,svm_tolerance=1e-6,*svm_grid=NULL;char *svm_grid_text=NULL;
   double lambda=1.0;int continuous=0,encoding=0,voxel_encoding=0,ridge_class=0,distance_given=0;
   int lda_mode=-1,shrinkage_given=0,floor_given=0,empirical_prior=0;
   double shrinkage=.1,variance_floor=1e-6;char *class_prior=NULL;
   char *grid_text=NULL,*inner_column=NULL;int inner_k=0,tune=0,numeric_lambda=0,nlambda=0;
   double *lambda_grid=NULL;
   char *save_predictions=NULL,*save_fold_scores=NULL ;
   char *permute_within=NULL ; long seed=1234567L ;
   char *distance_name="euclidean",*scale_name="feature_zscore" ;
   char *metric_name=NULL,err[1024],*table_name=NULL ;
   char **class_label=NULL,**perm_block_label=NULL ;
   int *class_id=NULL,*group_id=NULL,*perm_block=NULL,**perm_class=NULL ;
   float **x=NULL,**perm_score=NULL,**perm_max=NULL ;
   float *score=NULL,*acc=NULL,*bacc=NULL,*nfeat=NULL,*nfold=NULL,*qval=NULL ;
   int nclass=0,distance=THD_CENTROID_EUCLIDEAN,preprocess=THD_PREP_FEATURE_ZSCORE ;
   int metric=THD_METRIC_BALANCED_ACCURACY,cv_k=0,nperm=0,nperm_block=0 ;
   int show_table=0,show_samples=0,show_folds=0,quiet=0,nopt=1,ii,rr,maxvox ;
   int nworker=1,failed=0,io_failed=0 ; char runerr[1024]="" ;
   FILE *prediction_fp=NULL,*fold_fp=NULL,*class_tuning=NULL ;char *class_tunefile=NULL;

   if( argc<2 ){ usage_3dMVPA() ; return 0 ; }
   mainENTRY(PROGRAM_NAME " main") ; machdep() ; AFNI_logger(PROGRAM_NAME,argc,argv) ;
   PRINT_VERSION(PROGRAM_NAME) ; AUTHOR("P Molfese") ; THD_check_AFNI_version(PROGRAM_NAME) ;

   while( nopt<argc && argv[nopt][0]=='-' ){
#define IS(s) (strcasecmp(argv[nopt],(s))==0)
#define NEEDARG(what) do{ if(++nopt>=argc)ERROR_exit(PROGRAM_NAME ": need %s after '%s'",(what),argv[nopt-1]); }while(0)
     if( IS("-help") || IS("-h") ){ usage_3dMVPA() ; return 0 ; }
     else if( IS("-mode") ){ NEEDARG("classify, decode, or encode") ; mode=argv[nopt++] ; }
     else if( IS("-target_column") ){ NEEDARG("a column name") ; target_column=argv[nopt++] ; }
     else if( IS("-target_columns") ){ NEEDARG("numeric column names") ; target_column=argv[nopt++] ; }
     else if( IS("-design_columns") ){ NEEDARG("numeric column names") ; design_columns=argv[nopt++] ; }
     else if( IS("-score") ){ NEEDARG("voxel or pattern") ; score_mode=argv[nopt++] ; }
     else if( IS("-lambda") ){
       char *end;NEEDARG("a positive penalty or auto");
       if(tune || numeric_lambda)ERROR_exit(PROGRAM_NAME ": -lambda was specified twice");
       if(!strcasecmp(argv[nopt],"auto"))tune=1;
       else{
         errno=0;lambda=strtod(argv[nopt],&end);numeric_lambda=1;
         if(errno || end==argv[nopt] || *end || !isfinite(lambda) || lambda<=0)
           ERROR_exit(PROGRAM_NAME ": -lambda must be finite and positive, or auto");
       }
       nopt++;
     }
     else if(IS("-svm_c")){
       char *end;NEEDARG("a positive C or auto");svm_options=1;
       if(svm_tune || svm_numeric)ERROR_exit(PROGRAM_NAME ": -svm_c was specified twice");
       if(!strcasecmp(argv[nopt],"auto"))svm_tune=1;
       else{
         errno=0;svm_cost=strtod(argv[nopt],&end);svm_numeric=1;
         if(errno || end==argv[nopt] || *end || !isfinite(svm_cost) || svm_cost<=0 || !isfinite(.5/svm_cost) || .5/svm_cost<=0)
           ERROR_exit(PROGRAM_NAME ": -svm_c must be finite and positive, or auto");
       }
       nopt++;
     }
     else if(IS("-svm_c_grid")){
       NEEDARG("comma-separated positive C values");svm_options=1;
       if(svm_grid_text)ERROR_exit(PROGRAM_NAME ": -svm_c_grid was specified twice");
       svm_grid_text=argv[nopt++];
     }
     else if(IS("-svm_tol")){
       char *end;NEEDARG("a tolerance in [1e-12,1e-2]");svm_options=1;
       errno=0;svm_tolerance=strtod(argv[nopt],&end);
       if(errno || end==argv[nopt] || *end || !isfinite(svm_tolerance) || svm_tolerance<1e-12 || svm_tolerance>1e-2)
         ERROR_exit(PROGRAM_NAME ": -svm_tol must be in [1e-12,1e-2]");
       nopt++;
     }
     else if(IS("-svm_max_iter")){
       NEEDARG("a positive sweep limit");svm_options=1;
       svm_max_iter=(int)mvpa_parse_long("-svm_max_iter",argv[nopt++],1,INT_MAX);
     }
     else if( IS("-lambda_grid") ){
       NEEDARG("comma-separated positive penalties");
       if(grid_text)ERROR_exit(PROGRAM_NAME ": -lambda_grid was specified twice");
       grid_text=argv[nopt++];
     }
     else if( IS("-inner_cv_column") ){NEEDARG("an inner grouping column");inner_column=argv[nopt++];}
     else if( IS("-inner_cv_k") ){NEEDARG("K >= 2");inner_k=(int)mvpa_parse_long("-inner_cv_k",argv[nopt++],2,INT_MAX);}
     else if( IS("-cv_column") ){ NEEDARG("a column name") ; cv_column=argv[nopt++] ; }
     else if( IS("-cv") ){
       NEEDARG("run, subject, or item") ;
            if( strcasecmp(argv[nopt],"run")==0 ) cv_column="Run" ;
       else if( strcasecmp(argv[nopt],"subject")==0 ) cv_column="Subj" ;
       else if( strcasecmp(argv[nopt],"item")==0 ) cv_column="Item" ;
       else ERROR_exit(PROGRAM_NAME ": -cv must be run, subject, or item") ;
       nopt++ ;
     }
     else if( IS("-cv_k") ){ NEEDARG("K >= 2") ; cv_k=(int)mvpa_parse_long("-cv_k",argv[nopt++],2,INT_MAX) ; }
     else if( IS("-within_column") ){ NEEDARG("a column name") ; within_column=argv[nopt++] ; }
     else if( IS("-subject_column") ){ NEEDARG("a column name") ; subject_column=argv[nopt++] ; }
     else if( IS("-run_column") ){ NEEDARG("a column name") ; run_column=argv[nopt++] ; }
     else if( IS("-item_column") ){ NEEDARG("a column name") ; item_column=argv[nopt++] ; }
     else if( IS("-fold_column") ){ NEEDARG("a column name") ; fold_column=argv[nopt++] ; }
     else if( IS("-sample_column") ){ NEEDARG("a column name") ; sample_column=argv[nopt++] ; }
     else if( IS("-xmat_stims") ){ NEEDARG("a comma-separated stimulus list") ; xmat_stims=argv[nopt++] ; }
     else if( IS("-mask") ){ NEEDARG("a dataset") ; maskname=argv[nopt++] ; }
     else if( IS("-roi_sel") ){ NEEDARG("an AFNI integer list") ; roi_sel=argv[nopt++] ; }
     else if( IS("-searchlight") ){ NEEDARG("a radius or neighborhood") ; searchlight=argv[nopt++] ; }
     else if( IS("-surf") ){ NEEDARG("a surface mesh") ; surf_file=argv[nopt++] ; }
     else if(IS("-shrinkage") || IS("-variance_floor")){
       int isshrink=IS("-shrinkage");char *end;double value;
       if(isshrink?shrinkage_given:floor_given)ERROR_exit(PROGRAM_NAME ": repeated LDA regularization option");
       NEEDARG("a finite regularization value");errno=0;value=strtod(argv[nopt],&end);
       if(errno || end==argv[nopt] || *end || !isfinite(value) ||
          (isshrink?(value<=0 || value>1):(value<1e-12 || value>1)))
         ERROR_exit(PROGRAM_NAME ": %s must be in %s",isshrink?"-shrinkage":"-variance_floor",isshrink?"(0,1]":"[1e-12,1]");
       nopt++;
       if(isshrink){shrinkage=value;shrinkage_given=1;}else{variance_floor=value;floor_given=1;}
     }
     else if(IS("-class_prior")){
       if(class_prior)ERROR_exit(PROGRAM_NAME ": -class_prior was specified twice");
       NEEDARG("equal or empirical");class_prior=argv[nopt++];
       if(strcasecmp(class_prior,"equal") && strcasecmp(class_prior,"empirical"))
         ERROR_exit(PROGRAM_NAME ": -class_prior must be equal or empirical");
       empirical_prior=!strcasecmp(class_prior,"empirical");
     }
     else if( IS("-method") ){ NEEDARG("centroid, ridge, svm, lda_diag, or lda_shrinkage") ; method=argv[nopt++] ; }
     else if( IS("-distance") ){ NEEDARG("euclidean, correlation, or cosine") ; distance_name=argv[nopt++] ;distance_given=1; }
     else if( IS("-scale") ){ NEEDARG("none, feature_zscore, sample_center, or run_center") ; scale_name=argv[nopt++] ; }
     else if( IS("-metric") ){ NEEDARG("a prediction metric") ; metric_name=argv[nopt++] ; }
     else if( IS("-nperm") ){ NEEDARG("a nonnegative count") ; nperm=(int)mvpa_parse_long("-nperm",argv[nopt++],0,INT_MAX) ; }
     else if( IS("-seed") ){ NEEDARG("an integer seed") ; seed=mvpa_parse_long("-seed",argv[nopt++],LONG_MIN,LONG_MAX) ; }
     else if( IS("-permute_within") || IS("-block_column") ){
       NEEDARG("an exchangeability-block column or none") ;
       if( permute_within!=NULL ) ERROR_exit(PROGRAM_NAME ": the permutation block was given twice") ;
       permute_within=argv[nopt++] ;
     }
     else if( IS("-prefix") ){ NEEDARG("an output prefix") ; prefix=argv[nopt++] ; }
     else if(IS("-save_coefficients")){
       NEEDARG("a separate descriptive output prefix");
       if(coef_prefix)ERROR_exit(PROGRAM_NAME ": -save_coefficients was specified twice");
       coef_prefix=argv[nopt++];
     }
     else if(IS("-coef_lambda")){
       char *end;NEEDARG("a positive descriptive refit penalty");
       if(coef_lambda_given)ERROR_exit(PROGRAM_NAME ": -coef_lambda was specified twice");
       errno=0;coef_lambda=strtod(argv[nopt],&end);coef_lambda_given=1;
       if(errno || end==argv[nopt] || *end || !isfinite(coef_lambda) || coef_lambda<=0)
         ERROR_exit(PROGRAM_NAME ": -coef_lambda must be finite and positive");
       nopt++;
     }
     else if( IS("-save_predictions") ){ NEEDARG("a filename") ; save_predictions=argv[nopt++] ; }
     else if( IS("-save_fold_scores") ){ NEEDARG("a filename") ; save_fold_scores=argv[nopt++] ; }
     else if( IS("-show_table") ){ show_table=1 ; nopt++ ; }
     else if( IS("-show_samples") ){ show_samples=1 ; nopt++ ; }
     else if( IS("-show_folds") ){ show_folds=1 ; nopt++ ; }
     else if( IS("-quiet") ){ quiet=1 ; nopt++ ; }
     else if( IS("-trial_manifest") ){
       NEEDARG("a trial manifest filename") ;
       if( tab!=NULL ) ERROR_exit(PROGRAM_NAME ": choose one of -trial_manifest and -dataTable/-dataTableFile") ;
       tab=THD_read_datatable_file(argv[nopt++]) ; trial_manifest=1 ;
     }
     else if( IS("-dataTableFile") ){
       NEEDARG("a table file") ;
       if( tab!=NULL ) ERROR_exit(PROGRAM_NAME ": the data table was given twice") ;
       tab=THD_read_datatable_file(argv[nopt++]) ;
     }
     else if( IS("-dataTable") ){
       int nused=0 ; NEEDARG("a table or @file") ;
       if( tab!=NULL ) ERROR_exit(PROGRAM_NAME ": the data table was given twice") ;
       tab=THD_read_datatable_args(argc,argv,nopt,mvpa_options,&nused) ; nopt+=nused ;
     }
     else ERROR_exit(PROGRAM_NAME ": unknown option '%s'",argv[nopt]) ;
#undef NEEDARG
#undef IS
   }
   if( nopt<argc ) ERROR_exit(PROGRAM_NAME ": unexpected argument '%s'",argv[nopt]) ;
   if( mode==NULL || (strcasecmp(mode,"classify") && strcasecmp(mode,"decode") && strcasecmp(mode,"encode")) )
     ERROR_exit(PROGRAM_NAME ": -mode must be classify, decode, or encode") ;
   encoding=!strcasecmp(mode,"encode");continuous=encoding || !strcasecmp(mode,"decode");
   if(!method)method=continuous?"ridge":"centroid";
   ridge_class=!continuous && !strcasecmp(method,"ridge");
   svm_class=!continuous && !strcasecmp(method,"svm");
   if(svm_options && !svm_class)ERROR_exit(PROGRAM_NAME ": SVM options require -mode classify -method svm");
   if(svm_grid_text && svm_numeric)ERROR_exit(PROGRAM_NAME ": -svm_c_grid cannot accompany numeric -svm_c");
   if(svm_grid_text)svm_tune=1;
   if(!continuous && !strcasecmp(method,"lda_diag"))lda_mode=THD_LDA_DIAGONAL;
   if(!continuous && !strcasecmp(method,"lda_shrinkage"))lda_mode=THD_LDA_SHRINKAGE;
   if(lda_mode<0 && (shrinkage_given || floor_given || class_prior))
     ERROR_exit(PROGRAM_NAME ": shrinkage, variance floor and class-prior options require LDA classification");
   if(shrinkage_given && lda_mode!=THD_LDA_SHRINKAGE)
     ERROR_exit(PROGRAM_NAME ": -shrinkage requires -method lda_shrinkage");
   if(grid_text && numeric_lambda)ERROR_exit(PROGRAM_NAME ": -lambda_grid cannot accompany a numeric -lambda");
   if(grid_text)tune=1;
   if((tune || numeric_lambda) && !continuous && !ridge_class)
     ERROR_exit(PROGRAM_NAME ": ridge lambda and inner-CV options require -method ridge");
   if((ridge_class || svm_class || lda_mode>=0) && distance_given)ERROR_exit(PROGRAM_NAME ": -distance is only used by centroid classification");
   if(!tune && !svm_tune && (inner_column || inner_k))ERROR_exit(PROGRAM_NAME ": inner-CV options require tuned ridge or SVM penalties");
   if(tune)lambda_grid=mvpa_penalty_grid(grid_text?grid_text:"0.01,0.1,1,10,100",&nlambda,"lambda");
   if(svm_tune)svm_grid=mvpa_penalty_grid(svm_grid_text?svm_grid_text:"0.01,0.1,1,10,100",&svm_ngrid,"SVM C");
   if(coef_lambda_given && !coef_prefix)ERROR_exit(PROGRAM_NAME ": -coef_lambda requires -save_coefficients");
   if(coef_prefix){
     if(strcasecmp(method,"ridge"))ERROR_exit(PROGRAM_NAME ": -save_coefficients requires -method ridge");
     if(!THD_filename_ok(coef_prefix) || !strcmp(coef_prefix,prefix))
       ERROR_exit(PROGRAM_NAME ": descriptive coefficient prefix must be legal and different from -prefix");
     if(tune && !coef_lambda_given)
       ERROR_exit(PROGRAM_NAME ": tuned CV requires an explicit -coef_lambda for descriptive full-data refits");
     if(!coef_lambda_given)coef_lambda=lambda;
   }
   if(!metric_name)metric_name=continuous?"cross_validated_r2":"balanced_accuracy";
   if(encoding){
     if(!design_columns)ERROR_exit(PROGRAM_NAME ": encoding requires -design_columns");
     if(target_column)ERROR_exit(PROGRAM_NAME ": encoding responses come from InputFile; use -design_columns");
     if(!score_mode)score_mode="voxel";
     if(strcasecmp(score_mode,"voxel") && strcasecmp(score_mode,"pattern"))
       ERROR_exit(PROGRAM_NAME ": encoding -score must be voxel or pattern");
     if(!strcasecmp(score_mode,"voxel")){
       if(searchlight || surf_file || roi_sel)ERROR_exit(PROGRAM_NAME ": voxel encoding uses -mask without -searchlight, -surf, or -roi_sel");
       voxel_encoding=1;
     }
   }else if(design_columns || score_mode)ERROR_exit(PROGRAM_NAME ": -design_columns and -score require -mode encode");
   if( tab==NULL ) ERROR_exit(PROGRAM_NAME ": -dataTable or -dataTableFile is required") ;
   if( !encoding && target_column==NULL ) ERROR_exit(PROGRAM_NAME ": -target_column or -target_columns is required") ;
   if( cv_column==NULL ) ERROR_exit(PROGRAM_NAME ": -cv_column or -cv is required") ;
   if( nperm==1 ) ERROR_exit(PROGRAM_NAME ": -nperm must be 0 or at least 2") ;
   if( nperm>0 && permute_within==NULL )
     ERROR_exit(PROGRAM_NAME ": -nperm requires -permute_within COLUMN (or explicit 'none')") ;
   if( nperm==0 && permute_within!=NULL )
     ERROR_exit(PROGRAM_NAME ": -permute_within is only used with -nperm") ;
   if( surf_file!=NULL && searchlight==NULL )
     ERROR_exit(PROGRAM_NAME ": -surf is only used with -searchlight") ;
   if( maskname==NULL && !(surf_file!=NULL && searchlight!=NULL) )
     ERROR_exit(PROGRAM_NAME ": -mask is required, except for a surface searchlight over every represented node") ;
   if( searchlight!=NULL && roi_sel!=NULL )
     ERROR_exit(PROGRAM_NAME ": -roi_sel cannot be combined with -searchlight") ;
   if( save_predictions!=NULL && save_fold_scores!=NULL &&
       strcmp(save_predictions,save_fold_scores)==0 )
     ERROR_exit(PROGRAM_NAME ": -save_predictions and -save_fold_scores need different files") ;
   if(continuous?strcasecmp(method,"ridge")!=0:(!ridge_class && !svm_class && lda_mode<0 && strcasecmp(method,"centroid")!=0))
     ERROR_exit(PROGRAM_NAME ": mode %s requires -method %s",mode,continuous?"ridge":"centroid, ridge, svm, lda_diag, or lda_shrinkage") ;
   if( !THD_filename_ok(prefix) ) ERROR_exit(PROGRAM_NAME ": illegal -prefix '%s'",prefix) ;
   if( !THD_centroid_distance_parse(distance_name,&distance,err,sizeof(err)) )
     ERROR_exit(PROGRAM_NAME ": %s",err) ;
   if( !THD_preprocess_parse(scale_name,&preprocess,err,sizeof(err)) )
     ERROR_exit(PROGRAM_NAME ": %s",err) ;
   if( !THD_mvpa_metric_parse(metric_name,&metric,err,sizeof(err)) )
     ERROR_exit(PROGRAM_NAME ": %s",err) ;
   if( !continuous && metric!=THD_METRIC_ACCURACY && metric!=THD_METRIC_BALANCED_ACCURACY )
     ERROR_exit(PROGRAM_NAME ": classification -metric must be accuracy or balanced_accuracy") ;
   if(continuous && metric<THD_METRIC_PEARSON)
     ERROR_exit(PROGRAM_NAME ": continuous metric must be pearson, spearman, cross_validated_r2, MAE, or RMSE");

   if( trial_manifest ){
     if( xmat_stims!=NULL ) ERROR_exit(PROGRAM_NAME ": -xmat_stims cannot accompany -trial_manifest") ;
     expanded=THD_datatable_trial_manifest(tab) ;
     THD_free_datatable(tab) ; tab=expanded ; expanded=NULL ;
     if( tab==NULL ) ERROR_exit(PROGRAM_NAME ": trial manifest conversion failed") ;
   } else if( THD_datatable_column(tab,"Xmat")>=0 ){
     const char *icol=item_column?item_column:"Item" ;
     const char *rcol=run_column?run_column:"Run" ; int have_run=0 ;
     expanded=THD_datatable_expand_xmat_items(tab,"Xmat",icol,rcol,xmat_stims,&have_run) ;
     THD_free_datatable(tab) ; tab=expanded ; expanded=NULL ;
     if( tab==NULL ) ERROR_exit(PROGRAM_NAME ": X-matrix table expansion failed") ;
   } else if( xmat_stims!=NULL )
     ERROR_exit(PROGRAM_NAME ": -xmat_stims requires an Xmat column") ;
   if( show_table ) THD_datatable_print(tab,stdout) ;

   THD_sampleset_detect_columns(tab,&columns) ;
   if( subject_column!=NULL ) columns.subject=subject_column ;
   if( run_column!=NULL ) columns.run=run_column ;
   if( item_column!=NULL ) columns.item=item_column ;
   if( fold_column!=NULL ) columns.fold=fold_column ;
   if( sample_column!=NULL ) columns.sample=sample_column ;
   ss=THD_sampleset_from_table(tab,&columns,err,sizeof(err)) ;
   THD_free_datatable(tab) ; tab=NULL ;
   if( ss==NULL ) ERROR_exit(PROGRAM_NAME ": %s",err) ;
   if( show_samples ) THD_sampleset_print(ss,stdout) ;

   if( !continuous && !THD_sampleset_encode_column(ss,target_column,&class_id,&class_label,&nclass,err,sizeof(err)) )
     ERROR_exit(PROGRAM_NAME ": %s",err) ;
   if( !continuous && nclass<2 ) ERROR_exit(PROGRAM_NAME ": target column '%s' has only %d class",target_column,nclass) ;
   plan=(cv_k>0) ? THD_foldplan_grouped_kfold(ss,cv_column,within_column,cv_k,err,sizeof(err))
                 : THD_foldplan_leave_group_out(ss,cv_column,within_column,err,sizeof(err)) ;
   if( plan==NULL ) ERROR_exit(PROGRAM_NAME ": %s",err) ;
   if( !continuous && !THD_foldplan_check_classes(plan,class_id,nclass,0,err,sizeof(err)) )
     ERROR_exit(PROGRAM_NAME ": %s",err) ;
   if( show_folds ) THD_foldplan_print(plan,ss,stdout) ;
   if(tune || svm_tune){
     int *inner_groups=NULL,ng=0;char **labels=NULL;
     if(!inner_column)inner_column=cv_column;
     if(!THD_sampleset_encode_column(ss,inner_column,&inner_groups,&labels,&ng,err,sizeof(err)))
       ERROR_exit(PROGRAM_NAME ": %s",err);
     nested=THD_nested_foldplan_grouped(plan,inner_groups,inner_k,err,sizeof(err));
     free(inner_groups);THD_sampleset_free_levels(labels,ng);
     if(!nested)ERROR_exit(PROGRAM_NAME ": %s",err);
     if(show_folds)THD_nested_foldplan_print(nested,ss,stdout);
     if(!quiet)INFO_message(PROGRAM_NAME ": nested %s selection over %d penalties; inner groups from %s",svm_class?"SVM C":"ridge",svm_class?svm_ngrid:nlambda,inner_column);
   }

   if( nperm>0 ){
     int ip,*check_class ; long long group_size ;
     perm_scheme=THD_perm_scheme_new(ss->nsample) ;
     if( perm_scheme==NULL ) ERROR_exit(PROGRAM_NAME ": cannot allocate the label-permutation scheme") ;
     perm_scheme->exchange=PERM_EE ; perm_scheme->blockmode=PERM_WITHIN_BLOCK ;
     perm_scheme->exact=0 ; perm_scheme->nperm=nperm ; perm_scheme->seed=seed ;
     if( strcasecmp(permute_within,"none")!=0 ){
       if( !THD_sampleset_encode_column(ss,permute_within,&perm_block,
                                        &perm_block_label,&nperm_block,
                                        err,sizeof(err)) )
         ERROR_exit(PROGRAM_NAME ": %s",err) ;
       if( THD_perm_scheme_set_blocks(perm_scheme,perm_block)!=nperm_block )
         ERROR_exit(PROGRAM_NAME ": cannot encode exchangeability blocks from '%s'",permute_within) ;
     } else nperm_block=1 ;
     if( !continuous && !THD_perm_scheme_set_eqclass(perm_scheme,class_id) )
       ERROR_exit(PROGRAM_NAME ": cannot encode equivalent target labels for permutation") ;
     group_size=THD_perm_group_size(perm_scheme) ;
     if( group_size!=PERM_COUNT_TOOBIG && group_size>0 && group_size<=nperm )
       perm_scheme->exact=1 ;
     perm_set=THD_perm_set_build(perm_scheme) ;
     THD_perm_scheme_free(perm_scheme) ; perm_scheme=NULL ;
     THD_sampleset_free_levels(perm_block_label,nperm_block) ; perm_block_label=NULL ;
     free(perm_block) ; perm_block=NULL ;
     if( perm_set==NULL ) ERROR_exit(PROGRAM_NAME ": cannot build the label-permutation set") ;
     if( perm_set->nperm<2 )
       ERROR_exit(PROGRAM_NAME ": -permute_within %s permits only the identity label arrangement",
                  permute_within) ;
     if(!continuous){
     check_class=(int *)malloc(sizeof(int)*(size_t)ss->nsample) ;
     if( check_class==NULL ) ERROR_exit(PROGRAM_NAME ": out of memory validating label permutations") ;
     for( ip=0 ; ip<perm_set->nperm ; ip++ ){
       for( ii=0 ; ii<ss->nsample ; ii++ ){
         int src=perm_set->perm[(size_t)ip*ss->nsample+ii] ;
         if( src<0 || src>=ss->nsample ||
             perm_set->sign[(size_t)ip*ss->nsample+ii]!=1 )
           ERROR_exit(PROGRAM_NAME ": invalid label relabeling %d at sample %d",ip+1,ii+1) ;
         check_class[ii]=class_id[src] ;
       }
       if( !THD_foldplan_check_classes(plan,check_class,nclass,0,err,sizeof(err)) )
         ERROR_exit(PROGRAM_NAME ": label relabeling %d violates CV class coverage: %s;\n"
                    "       choose a finer -permute_within column whose blocks preserve every training fold",
                    ip+1,err) ;
     }
     free(check_class) ;
     }
     if( !quiet )
       INFO_message(PROGRAM_NAME ": %d synchronized %s relabelings (%s; %d block%s from %s)",
                    perm_set->nperm,continuous?"row":"label",perm_set->is_exact?"exact":"Monte Carlo",
                    nperm_block,nperm_block==1?"":"s",permute_within) ;
   }

   if( preprocess==THD_PREP_RUN_CENTER ){
     if( ss->run==NULL || ss->nrun<1 )
       ERROR_exit(PROGRAM_NAME ": -scale run_center requires a Run role; use -run_column if needed") ;
     if( ss->nsubject>INT_MAX/ss->nrun )
       ERROR_exit(PROGRAM_NAME ": too many subject-by-run groups for run_center") ;
     group_id=(int *)malloc(sizeof(int)*(size_t)ss->nsample) ;
     if( group_id==NULL ) ERROR_exit(PROGRAM_NAME ": out of memory encoding subject-by-run groups") ;
     for( ii=0 ; ii<ss->nsample ; ii++ ){
       if( ss->subject[ii]<0 || ss->run[ii]<0 )
         ERROR_exit(PROGRAM_NAME ": sample %d has no subject/run role for run_center",ii+1) ;
       group_id[ii]=ss->subject[ii]*ss->nrun+ss->run[ii] ;
     }
   }

   mset=THD_open_dataset(maskname!=NULL?maskname:ss->dset_name[0]) ;
   if( mset==NULL )
     ERROR_exit(PROGRAM_NAME ": cannot open %s dataset '%s'",
                maskname!=NULL?"mask":"surface-domain",maskname!=NULL?maskname:ss->dset_name[0]) ;
   if( surf_file!=NULL || DBLK_IS_NI_SURF_DSET(mset->dblk) || DBLK_IS_GIFTI(mset->dblk) ){
     for( ii=0 ; ii<ss->ndset ; ii++ )
       if( !THD_surface_domain_match(mset,ss->dset[ii],err,sizeof(err)) )
         ERROR_exit(PROGRAM_NAME ": %s",err) ;
   } else if( !EQUIV_GRIDS(mset,ss->dset[0]) )
     ERROR_exit(PROGRAM_NAME ": mask '%s' is not on the sample dataset grid",maskname) ;
   DSET_load(mset) ;
   if( !DSET_LOADED(mset) )
     ERROR_exit(PROGRAM_NAME ": cannot load %s dataset '%s'",
                maskname!=NULL?"mask":"surface-domain",maskname!=NULL?maskname:ss->dset_name[0]) ;
   if(voxel_encoding){
     rl=mvpa_voxel_units(mset);
   }else if( searchlight!=NULL ){
     if( surf_file!=NULL ){
#ifdef USE_SUMA
       float radius=mvpa_parse_positive_float("surface -searchlight radius",searchlight) ;
       char FuncName[]=PROGRAM_NAME ;
       SUMA_STANDALONE_INIT ;
       SUMAg_DOv=SUMA_Alloc_DisplayObject_Struct(SUMA_MAX_DISPLAYABLE_OBJECTS) ;
       if( !quiet )
         INFO_message(PROGRAM_NAME ": geodesic surface searchlight, radius %.3g mm on '%s'",
                      radius,surf_file) ;
       rl=THD_roilist_searchlight_surf(surf_file,mset,radius,maskname==NULL,
                                      err,sizeof(err)) ;
       if( rl==NULL && err[0]!='\0' ) ERROR_exit(PROGRAM_NAME ": %s",err) ;
       if( rl==NULL || rl->nroi<1 )
         ERROR_exit(PROGRAM_NAME ": surface searchlight has no represented in-mask nodes") ;
#else
       ERROR_exit(PROGRAM_NAME ": -surf needs a build with surface support (USE_SUMA/libSUMA)") ;
#endif
     } else {
       MCW_cluster *nbhd ;
       if( DBLK_IS_NI_SURF_DSET(mset->dblk) || DBLK_IS_GIFTI(mset->dblk) )
         ERROR_exit(PROGRAM_NAME ": a surface -searchlight also requires -surf MESH") ;
       nbhd=THD_searchlight_parse(searchlight,fabsf(DSET_DX(mset)),
                                  fabsf(DSET_DY(mset)),fabsf(DSET_DZ(mset)),
                                  err,sizeof(err)) ;
       if( nbhd==NULL ) ERROR_exit(PROGRAM_NAME ": %s",err) ;
       if( !quiet ) INFO_message(PROGRAM_NAME ": searchlight neighborhood '%s' = %d voxel offsets",
                                 searchlight,nbhd->num_pt) ;
       rl=THD_roilist_searchlight(mset,nbhd) ;
       KILL_CLUSTER(nbhd) ;
       if( rl==NULL || rl->nroi<1 )
         ERROR_exit(PROGRAM_NAME ": mask has no nonzero searchlight centers") ;
     }
   } else {
     rl=THD_roilist_from_dset(mset,roi_sel) ;
     if( rl==NULL || rl->nroi<1 )
       ERROR_exit(PROGRAM_NAME ": mask has no selected positive-valued ROIs") ;
   }
   if(coef_prefix){
     const char *suffix[]={".coefficients.tsv",".preprocessing.tsv"};int c;
     char *cf=malloc(strlen(coef_prefix)+32);
     if(!cf)ERROR_exit(PROGRAM_NAME ": cannot name coefficient audits");
     for(c=0;c<2;c++){
       sprintf(cf,"%s%s",coef_prefix,suffix[c]);
       if((save_predictions && !strcmp(cf,save_predictions)) ||
          (save_fold_scores && !strcmp(cf,save_fold_scores)))
         ERROR_exit(PROGRAM_NAME ": coefficient and prediction/fold audit filenames must differ");
       if(THD_deathcon() && THD_is_file(cf))ERROR_exit(PROGRAM_NAME ": coefficient audit '%s' already exists",cf);
     }
     free(cf);
     {
       int c,combined=encoding && searchlight==NULL;
       char *name=malloc(strlen(coef_prefix)+64);
       THD_3dim_dataset *score_check=EDIT_empty_copy(mset),*coef_check=EDIT_empty_copy(mset);
       if(!name || !score_check || !coef_check)ERROR_exit(PROGRAM_NAME ": cannot check coefficient output paths");
       EDIT_dset_items(score_check,ADN_prefix,prefix,ADN_none);
       for(c=0;c<(combined?1:rl->nroi);c++){
         mvpa_coefficient_prefix(name,strlen(coef_prefix)+64,mset,coef_prefix,combined,rl->val[c]);
         EDIT_dset_items(coef_check,ADN_prefix,name,ADN_none);
         if(!strcmp(DSET_HEADNAME(score_check),DSET_HEADNAME(coef_check)))
           ERROR_exit(PROGRAM_NAME ": coefficient maps and held-out score dataset must differ");
         if(THD_deathcon() && THD_is_file(DSET_HEADNAME(coef_check)))
           ERROR_exit(PROGRAM_NAME ": coefficient dataset '%s' already exists",name);
       }
       DSET_delete(score_check);DSET_delete(coef_check);free(name);
     }
   }
   maxvox=THD_roilist_maxvox(rl) ;
   if( maxvox<1 || (size_t)ss->nsample>SIZE_MAX/sizeof(float)/(size_t)maxvox )
     ERROR_exit(PROGRAM_NAME ": sample-by-voxel workspace is too large") ;
#ifdef USE_OMP
   if( searchlight!=NULL || voxel_encoding ) nworker=omp_get_max_threads() ;
#endif
   if( nworker>rl->nroi ) nworker=rl->nroi ;
   for( ii=0 ; ii<ss->ndset ; ii++ ){
     DSET_load(ss->dset[ii]) ;
     if( !DSET_LOADED(ss->dset[ii]) )
       ERROR_exit(PROGRAM_NAME ": cannot load sample dataset '%s'",ss->dset_name[ii]) ;
   }

   if(continuous){
     mvpa_continuous_analysis(ss,plan,mset,rl,encoding,encoding?design_columns:target_column,
       preprocess,group_id,lambda,metric,prefix,save_predictions,save_fold_scores,
       nworker,quiet,nested,lambda_grid,nlambda,inner_column,perm_set,permute_within,argc,argv);
     if(coef_prefix)mvpa_descriptive_coefficients(ss,mset,rl,encoding,searchlight!=NULL,
       encoding?design_columns:target_column,NULL,NULL,0,preprocess,group_id,coef_lambda,
       coef_prefix,quiet,argc,argv);
     DSET_delete(mset);THD_roilist_free(rl);THD_nested_foldplan_free(nested);THD_foldplan_free(plan);
     free(lambda_grid);
     THD_perm_set_free(perm_set);
     THD_sampleset_free(ss);free(group_id);return 0;
   }

   x=(float **)calloc((size_t)nworker,sizeof(float *)) ;
   workspace=(THD_mvpa_workspace **)calloc((size_t)nworker,sizeof(*workspace)) ;
   if( perm_set!=NULL ){
     perm_class=(int **)calloc((size_t)nworker,sizeof(*perm_class)) ;
     perm_score=(float **)calloc((size_t)nworker,sizeof(*perm_score)) ;
     perm_max=(float **)calloc((size_t)nworker,sizeof(*perm_max)) ;
     perm_result=THD_perm_result_new(rl->nroi,perm_set->nperm) ;
     qval=(float *)malloc(sizeof(float)*(size_t)rl->nroi) ;
   }
   score=(float *)malloc(sizeof(float)*(size_t)rl->nroi) ;
   acc=(float *)malloc(sizeof(float)*(size_t)rl->nroi) ;
   bacc=(float *)malloc(sizeof(float)*(size_t)rl->nroi) ;
   nfeat=(float *)malloc(sizeof(float)*(size_t)rl->nroi) ;
   nfold=(float *)malloc(sizeof(float)*(size_t)rl->nroi) ;
   if( x==NULL || workspace==NULL || score==NULL || acc==NULL || bacc==NULL ||
       nfeat==NULL || nfold==NULL ||
       (perm_set!=NULL && (perm_class==NULL || perm_score==NULL || perm_max==NULL ||
                           perm_result==NULL || qval==NULL)) )
     ERROR_exit(PROGRAM_NAME ": out of memory allocating spatial analysis buffers") ;
   if( perm_result!=NULL ){
     perm_result->tail=PERM_TAIL_ONE ;
     for( ii=0 ; ii<perm_result->nperm ; ii++ ) perm_result->max_null[ii]=-FLT_MAX ;
   }
   for( ii=0 ; ii<nworker ; ii++ ){
     int ip ;
     x[ii]=(float *)malloc(sizeof(float)*(size_t)ss->nsample*maxvox) ;
     if( x[ii]==NULL ) ERROR_exit(PROGRAM_NAME ": out of memory allocating worker patterns") ;
     workspace[ii]=THD_mvpa_workspace_alloc(plan,maxvox,nclass,preprocess,
                                             group_id,err,sizeof(err)) ;
     if( workspace[ii]==NULL ) ERROR_exit(PROGRAM_NAME ": %s",err) ;
     if(ridge_class && !THD_mvpa_classification_set_ridge(workspace[ii],lambda,nested,
                                     lambda_grid,nlambda,err,sizeof(err)))
       ERROR_exit(PROGRAM_NAME ": %s",err);
     if(svm_class && !THD_mvpa_classification_set_svm(workspace[ii],svm_cost,svm_tolerance,
       svm_max_iter,metric,nested,svm_grid,svm_ngrid,err,sizeof(err)))ERROR_exit(PROGRAM_NAME ": %s",err);
     if(lda_mode>=0 && !THD_mvpa_classification_set_lda(workspace[ii],lda_mode,shrinkage,
                          variance_floor,empirical_prior,err,sizeof(err)))ERROR_exit(PROGRAM_NAME ": %s",err);
     if( perm_set!=NULL ){
       perm_class[ii]=(int *)malloc(sizeof(int)*(size_t)ss->nsample) ;
       perm_score[ii]=(float *)malloc(sizeof(float)*(size_t)perm_set->nperm) ;
       perm_max[ii]=(float *)malloc(sizeof(float)*(size_t)perm_set->nperm) ;
       if( perm_class[ii]==NULL || perm_score[ii]==NULL || perm_max[ii]==NULL )
         ERROR_exit(PROGRAM_NAME ": out of memory allocating permutation worker buffers") ;
       for( ip=0 ; ip<perm_set->nperm ; ip++ ) perm_max[ii][ip]=-FLT_MAX ;
     }
   }
   if( searchlight==NULL ){
     table_name=(char *)malloc(strlen(prefix)+9) ;
     if( table_name==NULL ) ERROR_exit(PROGRAM_NAME ": out of memory naming the ROI table") ;
     sprintf(table_name,"%s.roi.tsv",prefix) ;
     if( save_predictions!=NULL && strcmp(save_predictions,table_name)==0 )
       ERROR_exit(PROGRAM_NAME ": -save_predictions cannot replace ROI table '%s'",table_name) ;
     if( save_fold_scores!=NULL && strcmp(save_fold_scores,table_name)==0 )
       ERROR_exit(PROGRAM_NAME ": -save_fold_scores cannot replace ROI table '%s'",table_name) ;
   }
   if( save_predictions!=NULL ){
     prediction_fp=mvpa_open_audit(save_predictions,"prediction audit") ;
     mvpa_write_audit_header(prediction_fp,"row-level held-out predictions",
                             target_column,class_label,nclass,searchlight!=NULL,
                             surf_file!=NULL) ;
     if(svm_class)fprintf(prediction_fp,"# Method: svm; loss=squared_hinge; intercept_penalized=1; C=%.17g; tuning=%s; tolerance=%.17g; max_iter=%d; decisions=raw_one_vs_rest_margins_not_probabilities\n",svm_cost,nested?"nested":"fixed",svm_tolerance,svm_max_iter);
     if(ridge_class)fputs("# Method: ridge; DecisionClass columns are 0/1-regression scores, not probabilities.\n",prediction_fp);
     if(lda_mode>=0)fprintf(prediction_fp,"# Method: %s; class_prior=%s; shrinkage=%.17g; relative_variance_floor=%.17g; decisions=centered_linear_discriminants_not_probabilities\n",
                           method,empirical_prior?"empirical":"equal",lda_mode==THD_LDA_DIAGONAL?0:shrinkage,variance_floor);
     fputs("UnitType\tUnitID\tUnitLabel\tI\tJ\tK\tNFeat\tTableRow\tSample\tSubj\tRun\tItem\tInputFold\tCVFold\tCVStratum\tObserved\tPredicted\tCorrect",prediction_fp) ;
     for( ii=0 ; ii<nclass ; ii++ ) fprintf(prediction_fp,"\tDecisionClass%d",ii) ;
     fputc('\n',prediction_fp) ;
   }
   if( save_fold_scores!=NULL ){
     fold_fp=mvpa_open_audit(save_fold_scores,"fold-score audit") ;
     mvpa_write_audit_header(fold_fp,"per-fold classification scores",
                             target_column,class_label,nclass,searchlight!=NULL,
                             surf_file!=NULL) ;
     if(svm_class)fprintf(fold_fp,"# Method: svm; loss=squared_hinge; intercept_penalized=1; tolerance=%.17g; max_iter=%d; inner_objective=1-%s; ties=smallest_C\n",svm_tolerance,svm_max_iter,THD_mvpa_metric_name(metric));
     if(lda_mode>=0)fprintf(fold_fp,"# Method: %s; class_prior=%s; covariance_denominator=ntrain-nclass; relative_variance_floor=%.17g\n",
                           method,empirical_prior?"empirical":"equal",variance_floor);
     fputs("UnitType\tUnitID\tUnitLabel\tI\tJ\tK\tNFeat\tCVFold\tStratum\tNTrain\tNTest\tACC\tBACC",fold_fp) ;
     if(svm_class)fputs("\tC\tInnerNFold\tInnerLoss\tSolverIterations\tSolverResidual",fold_fp);
     if(ridge_class)fputs("\tLambda\tInnerNFold\tInnerNValid\tInnerLoss",fold_fp);
     if(lda_mode>=0)fputs("\tShrinkage\tCovarianceScale\tVarianceFloor\tCovarianceForm",fold_fp);
     fputc('\n',fold_fp);
   }
   if((ridge_class || svm_class) && nested){
     class_tunefile=malloc(strlen(prefix)+20);
     if(!class_tunefile)ERROR_exit(PROGRAM_NAME ": cannot name classification tuning audit");
     sprintf(class_tunefile,"%s.tuning.tsv",prefix);
     if((save_predictions && !strcmp(save_predictions,class_tunefile)) ||
        (save_fold_scores && !strcmp(save_fold_scores,class_tunefile)))
       ERROR_exit(PROGRAM_NAME ": tuning audit and prediction/fold filenames must differ");
     class_tuning=mvpa_open_audit(class_tunefile,"classification tuning audit");
     if(svm_class){
       fprintf(class_tuning,"# method=svm; inner_group=%s; objective=1-%s; ties=smallest_C; intercept_penalized=1\n",inner_column,THD_mvpa_metric_name(metric));
       fputs("Unit\tCVFold\tCVStratum\tInnerNFold\tC\tInnerLoss\tSelected\n",class_tuning);
     }else{
     fprintf(class_tuning,"# method=ridge; coding=one_hot_0_1; inner_group=%s; objective=mean_class_SSE_over_training_mean_baseline_SSE; ties=largest_lambda\n",inner_column);
     fputs("Unit\tCVFold\tCVStratum\tInnerNFold\tInnerNValid\tLambda\tInnerLoss\tSelected\n",class_tuning);
     }
   }
   if( searchlight!=NULL && !quiet )
     INFO_message(PROGRAM_NAME ": %d centers, at most %d features, %d worker%s",
                  rl->nroi,maxvox,nworker,nworker==1?"":"s") ;

   AFNI_OMP_START ;
#ifdef USE_OMP
#pragma omp parallel for if(searchlight != NULL && rl->nroi > 1) \
                         num_threads(nworker) schedule(dynamic,1) ordered
#endif
   for( rr=0 ; rr<rl->nroi ; rr++ ){
     int tid=0,ok ; char localerr[1024] ;
#ifdef USE_OMP
     tid=omp_get_thread_num() ;
#endif
     ok=THD_roi_samples(ss,NULL,ss->nsample,rl->vox+rr,x[tid],
                        localerr,sizeof(localerr)) ;
     if( ok && perm_set!=NULL )
       ok=THD_mvpa_run_classification_permutations(plan,distance,x[tid],rl->vox[rr].nar,
                                             class_id,group_id,metric,perm_set,
                                             perm_class[tid],perm_score[tid],
                                             workspace[tid],localerr,sizeof(localerr)) ;
     else if( ok )
       ok=THD_mvpa_run_classification_cv(plan,distance,x[tid],rl->vox[rr].nar,
                                   class_id,group_id,workspace[tid],
                                   localerr,sizeof(localerr)) ;
     if( !ok ){
#ifdef USE_OMP
#pragma omp critical(mvpa_searchlight_error)
#endif
       {
         if( !failed ){
           failed=1 ;
           snprintf(runerr,sizeof(runerr),"%s %d: %.900s",
                    surf_file?"surface searchlight mesh node":
                    searchlight?"searchlight center voxel":"ROI",
                    rl->val[rr],localerr) ;
         }
       }
     } else {
       acc[rr]=(float)workspace[tid]->score.accuracy ;
       bacc[rr]=(float)workspace[tid]->score.balanced_accuracy ;
       score[rr]=(metric==THD_METRIC_ACCURACY)?acc[rr]:bacc[rr] ;
       nfeat[rr]=(float)rl->vox[rr].nar ; nfold[rr]=(float)plan->nfold ;
       if( perm_result!=NULL ){
         double mean=0.0 ; int ip ;
         perm_result->stat[rr]=score[rr] ;
         for( ip=0 ; ip<perm_set->nperm ; ip++ ){
           float vv=perm_score[tid][ip] ; mean+=vv ;
           if( vv>=score[rr] ) perm_result->cnt_unc[rr]++ ;
           if( vv>perm_max[tid][ip] ) perm_max[tid][ip]=vv ;
         }
         perm_result->effect[rr]=score[rr]-(float)(mean/(double)perm_set->nperm) ;
       }
     }
#ifdef USE_OMP
#pragma omp ordered
#endif
     {
       if( ok && prediction_fp!=NULL )
         mvpa_write_prediction_rows(prediction_fp,mset,rl,rr,searchlight!=NULL,
                                    surf_file!=NULL,
                                    ss,plan,class_id,class_label,workspace[tid]) ;
       if( ok && fold_fp!=NULL )
         mvpa_write_fold_rows(fold_fp,mset,rl,rr,searchlight!=NULL,
                              surf_file!=NULL,plan,workspace[tid]) ;
       if(ok && class_tuning && svm_class){
         int f,l;const THD_mvpa_workspace *sw=workspace[tid];
         for(f=0;f<plan->nfold;f++)for(l=0;l<svm_ngrid;l++){
           fprintf(class_tuning,"%d\t",rl->val[rr]);mvpa_tsv_text(class_tuning,plan->name[f]);
           fputc('\t',class_tuning);mvpa_tsv_text(class_tuning,plan->stratum_label[plan->fold_stratum[f]]);
           fprintf(class_tuning,"\t%d\t%.17g\t%.17g\t%d\n",nested->inner[f]->nfold,
             svm_grid[l],sw->svm_loss[(size_t)f*svm_ngrid+l],svm_grid[l]==sw->svm_selected[f]);
         }
       }
       if(ok && class_tuning && ridge_class){
         int f,l;const THD_mvpa_ridge_workspace *rw=workspace[tid]->ridge_cv;
         for(f=0;f<plan->nfold;f++)for(l=0;l<nlambda;l++){
           fprintf(class_tuning,"%d\t",rl->val[rr]);mvpa_tsv_text(class_tuning,plan->name[f]);
           fputc('\t',class_tuning);mvpa_tsv_text(class_tuning,plan->stratum_label[plan->fold_stratum[f]]);
           fprintf(class_tuning,"\t%d\t%d\t%.17g\t%.17g\t%d\n",nested->inner[f]->nfold,
             rw->inner_nvalid[f],lambda_grid[l],rw->inner_loss[(size_t)f*nlambda+l],lambda_grid[l]==rw->selected_lambda[f]);
         }
       }
     }
   }
   AFNI_OMP_END ;
   if( prediction_fp!=NULL ){
     if( ferror(prediction_fp) ) io_failed=1 ;
     if( fclose(prediction_fp)!=0 ) io_failed=1 ;
     prediction_fp=NULL ;
   }
   if( fold_fp!=NULL ){
     if( ferror(fold_fp) ) io_failed=1 ;
     if( fclose(fold_fp)!=0 ) io_failed=1 ;
     fold_fp=NULL ;
   }
   if(class_tuning && (ferror(class_tuning) || fclose(class_tuning)))io_failed=1;
   if( failed ) ERROR_exit(PROGRAM_NAME ": %s",runerr) ;
   if( io_failed ) ERROR_exit(PROGRAM_NAME ": error writing prediction/fold audit table") ;

   if( perm_result!=NULL ){
     int ip,ww ;
     for( ww=0 ; ww<nworker ; ww++ )
       THD_max_accum(perm_set->nperm,perm_result->max_null,perm_max[ww]) ;
     for( ip=0 ; ip<perm_set->nperm ; ip++ )
       if( perm_result->max_null[ip]==-FLT_MAX ) perm_result->max_null[ip]=0.0f ;
     THD_perm_result_finish(perm_result,NULL) ;
     THD_bh_fdr_masked(rl->nroi,perm_result->p_unc,NULL,qval) ;
   }

   oset=mvpa_make_roi_map(mset,rl,prefix,score,nfeat,nfold,perm_result,qval,
                          metric==THD_METRIC_ACCURACY?"ACC":"BACC",argc,argv) ;
   if( searchlight==NULL ){
     mvpa_write_roi_table(table_name,rl,score,acc,bacc,nfeat,nfold,
                          THD_mvpa_metric_name(metric),target_column,class_label,nclass,
                          perm_result,qval,permute_within) ;
   }
   if( DSET_write(oset)!=True ) ERROR_exit(PROGRAM_NAME ": failed writing output dataset") ;
   if( !quiet ){
     INFO_message(PROGRAM_NAME ": classified %d samples in %d %s%s with %d folds",
                  ss->nsample,rl->nroi,
                  surf_file?"surface searchlight mesh node":
                  searchlight?"searchlight center":"ROI",
                  rl->nroi==1?"":"s",plan->nfold) ;
     if( table_name!=NULL ) INFO_message(PROGRAM_NAME ": wrote ROI table %s",table_name) ;
     if( save_predictions!=NULL ) INFO_message(PROGRAM_NAME ": wrote predictions %s",save_predictions) ;
     if( save_fold_scores!=NULL ) INFO_message(PROGRAM_NAME ": wrote fold scores %s",save_fold_scores) ;
     if( perm_result!=NULL )
       INFO_message(PROGRAM_NAME ": inference family: %d spatial units x %d synchronized relabelings",
                    rl->nroi,perm_result->nperm) ;
     INFO_message(PROGRAM_NAME ": wrote dataset %s",DSET_BRIKNAME(oset)) ;
   }

   if(coef_prefix)mvpa_descriptive_coefficients(ss,mset,rl,0,searchlight!=NULL,
     target_column,class_id,class_label,nclass,preprocess,group_id,coef_lambda,
     coef_prefix,quiet,argc,argv);
   DSET_delete(oset) ; DSET_delete(mset) ; THD_roilist_free(rl) ;
   for( ii=0 ; ii<nworker ; ii++ ){
     THD_mvpa_workspace_free(workspace[ii]) ; free(x[ii]) ;
     if( perm_set!=NULL ){
       free(perm_class[ii]) ; free(perm_score[ii]) ; free(perm_max[ii]) ;
     }
   }
   free(workspace) ; free(x) ; THD_nested_foldplan_free(nested);THD_foldplan_free(plan) ;
   free(perm_class) ; free(perm_score) ; free(perm_max) ;
   THD_perm_result_free(perm_result) ; THD_perm_set_free(perm_set) ;
   THD_sampleset_free_levels(class_label,nclass) ; free(class_id) ; free(group_id) ;
   THD_sampleset_free(ss) ; free(score) ; free(acc) ; free(bacc) ;
   free(nfeat) ; free(nfold) ; free(qval) ; free(table_name) ;
   free(lambda_grid);free(svm_grid);free(class_tunefile);
   return 0 ;
}
