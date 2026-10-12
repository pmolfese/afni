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
  First vertical slice: categorical prediction in ROIs and volume searchlights.
  Compact input adapters end at THD_sampleset; everything after that is shared.
                                                     -- P Molfese, Oct 2026
------------------------------------------------------------------------------*/

static char *mvpa_options[]={
   "-mode","-dataTable","-dataTableFile","-xmat_stims",
   "-target_column","-cv_column","-cv","-cv_k","-within_column",
   "-subject_column","-run_column","-item_column","-fold_column","-sample_column",
   "-mask","-roi_sel","-searchlight","-surf","-method","-distance","-scale","-metric","-prefix",
   "-nperm","-seed","-permute_within","-block_column",
   "-save_predictions","-save_fold_scores",
   "-show_table","-show_samples","-show_folds","-quiet","-help","-h",NULL
} ;

static void usage_3dMVPA( void )
{
   printf(
"\n"
"Usage: 3dMVPA [options]                                      ~1~\n"
"\n"
"Cross-validated multivariate prediction from AFNI datasets.  This first\n"
"release slice performs categorical nearest-centroid classification in atlas\n"
"ROIs or moving volume/surface neighborhoods.  Learned preprocessing and\n"
"class centroids are fitted separately inside every training fold.\n"
"\n"
"Required options:                                           ~1~\n"
"\n"
"  -mode classify         Classification is currently the implemented mode.\n"
"\n"
"  -dataTable @TABLE\n"
"  -dataTableFile TABLE\n"
"  -dataTable Subj Run Condition InputFile ...\n"
"                         One row per sample. InputFile must select exactly\n"
"                         one sub-brick when its dataset has more than one.\n"
"                         A subject-level table with Xmat and whole-bucket\n"
"                         InputFile columns is expanded to the same sample\n"
"                         representation before analysis.\n"
"\n"
"  -target_column NAME    Categorical outcome column.\n"
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
"  -method centroid       Nearest class centroid [currently the only method].\n"
"  -distance D            euclidean [default], correlation, or cosine.\n"
"  -scale S               feature_zscore [default], sample_center, run_center,\n"
"                         or none. All learned values are training-fold only.\n"
"                         run_center requires a run role and cannot be applied\n"
"                         to a run absent from that fold's training data.\n"
"  -metric M              balanced_accuracy [default] or accuracy. This selects\n"
"                         the SCORE map; the ROI table reports both.\n"
"\n"
"Permutation inference:                                    ~1~\n"
"\n"
"  -nperm N               Run N synchronized label relabelings [0: off]. The\n"
"                         identity is included. If the allowed group has at\n"
"                         most N distinct relabelings, it is enumerated exactly.\n"
"  -permute_within COL    Shuffle labels only among rows sharing COL. Required\n"
"                         with -nperm; use 'none' only for an explicitly\n"
"                         unrestricted null. -block_column is an alias. Choose a\n"
"                         column that fully identifies the exchangeability unit,\n"
"                         such as Run, Subj, Fold, or a composite block column.\n"
"  -seed S                Monte-Carlo relabeling seed [1234567]. Exact tests do\n"
"                         not depend on this value. Every relabeling must leave\n"
"                         every CV training fold with all target classes.\n"
"\n"
"Inference uses the selected metric and an upper-tail label-permutation null.\n"
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
"  -save_predictions FILE Write one held-out prediction row per sample and\n"
"                         spatial unit, including sample metadata, test fold,\n"
"                         observed/predicted class labels, and every class\n"
"                         decision value. Rows never tested are retained as NA.\n"
"  -save_fold_scores FILE Write ACC and BACC for every fold and spatial unit,\n"
"                         with stratum and train/test sample counts. These two\n"
"                         audit tables are opt-in because searchlight tables\n"
"                         can be very large.\n"
"  -show_table            Print the canonical table after optional Xmat expansion.\n"
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

static int mvpa_is_surface_dset( const THD_3dim_dataset *dset )
{
   return dset!=NULL &&
          (DBLK_IS_NI_SURF_DSET(dset->dblk) || DBLK_IS_GIFTI(dset->dblk)) ;
}

static int mvpa_surface_domains_match( THD_3dim_dataset *aa,
                                       THD_3dim_dataset *bb,
                                       char *err, size_t errlen )
{
   int ii,na,nb,sa,sb ;
   if( !mvpa_is_surface_dset(aa) || !mvpa_is_surface_dset(bb) ){
     snprintf(err,errlen,"surface searchlight inputs must all be NIML/GIFTI surface datasets") ;
     return 0 ;
   }
   na=DSET_NVOX(aa) ; nb=DSET_NVOX(bb) ;
   if( na!=nb ){
     snprintf(err,errlen,"surface node-row domains differ: '%.350s' has %d rows but '%.350s' has %d",
              DSET_HEADNAME(aa),na,DSET_HEADNAME(bb),nb) ;
     return 0 ;
   }
   sa=(aa->dblk->node_list!=NULL) ; sb=(bb->dblk->node_list!=NULL) ;
   if( (sa && aa->dblk->nnodes!=na) || (sb && bb->dblk->nnodes!=nb) ){
     snprintf(err,errlen,"surface node-row metadata is malformed in '%.430s' or '%.430s'",
              DSET_HEADNAME(aa),DSET_HEADNAME(bb)) ;
     return 0 ;
   }
   for( ii=0 ; ii<na ; ii++ ){
     int an=sa?aa->dblk->node_list[ii]:ii , bn=sb?bb->dblk->node_list[ii]:ii ;
     if( an!=bn ){
       snprintf(err,errlen,
                "surface node-row domains differ at row %d: '%.350s' has node %d, '%.350s' has node %d",
                ii+1,DSET_HEADNAME(aa),an,DSET_HEADNAME(bb),bn) ;
       return 0 ;
     }
   }
   return 1 ;
}

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
   if( (DBLK_IS_NI_SURF_DSET(mset->dblk)||DBLK_IS_GIFTI(mset->dblk)) &&
       mset->dblk->nnodes>0 && mset->dblk->node_list!=NULL ){
     int nn=mset->dblk->nnodes ;
     oset->dblk->node_list=(int *)RwcMalloc(sizeof(int)*(size_t)nn) ;
     oset->dblk->nnodes=nn ;
     memcpy(oset->dblk->node_list,mset->dblk->node_list,sizeof(int)*(size_t)nn) ;
   }
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
     fprintf(fp,"\t%d\t%d\t%.17g\t%.17g\n",plan->ntrain[ff],
             workspace->fold_nsample[ff],workspace->fold_accuracy[ff],
             workspace->fold_balanced_accuracy[ff]) ;
   }
}

int main( int argc, char **argv )
{
   THD_datatable *tab=NULL,*expanded=NULL ; THD_sampleset *ss=NULL ;
   THD_sampleset_columns columns ; THD_foldplan *plan=NULL ;
   THD_3dim_dataset *mset=NULL,*oset=NULL ; THD_roilist *rl=NULL ;
   THD_mvpa_workspace **workspace=NULL ;
   PERM_scheme *perm_scheme=NULL ; PERM_set *perm_set=NULL ;
   PERM_result *perm_result=NULL ;
   char *mode=NULL,*target_column=NULL,*cv_column=NULL,*within_column=NULL ;
   char *subject_column=NULL,*run_column=NULL,*item_column=NULL ;
   char *fold_column=NULL,*sample_column=NULL,*xmat_stims=NULL ;
   char *maskname=NULL,*roi_sel=NULL,*searchlight=NULL,*surf_file=NULL ;
   char *prefix="MVPA",*method="centroid" ;
   char *save_predictions=NULL,*save_fold_scores=NULL ;
   char *permute_within=NULL ; long seed=1234567L ;
   char *distance_name="euclidean",*scale_name="feature_zscore" ;
   char *metric_name="balanced_accuracy",err[1024],*table_name=NULL ;
   char **class_label=NULL,**perm_block_label=NULL ;
   int *class_id=NULL,*group_id=NULL,*perm_block=NULL,**perm_class=NULL ;
   float **x=NULL,**perm_score=NULL,**perm_max=NULL ;
   float *score=NULL,*acc=NULL,*bacc=NULL,*nfeat=NULL,*nfold=NULL,*qval=NULL ;
   int nclass=0,distance=THD_CENTROID_EUCLIDEAN,preprocess=THD_PREP_FEATURE_ZSCORE ;
   int metric=THD_METRIC_BALANCED_ACCURACY,cv_k=0,nperm=0,nperm_block=0 ;
   int show_table=0,show_samples=0,show_folds=0,quiet=0,nopt=1,ii,rr,maxvox ;
   int nworker=1,failed=0,io_failed=0 ; char runerr[1024]="" ;
   FILE *prediction_fp=NULL,*fold_fp=NULL ;

   if( argc<2 ){ usage_3dMVPA() ; return 0 ; }
   mainENTRY(PROGRAM_NAME " main") ; machdep() ; AFNI_logger(PROGRAM_NAME,argc,argv) ;
   PRINT_VERSION(PROGRAM_NAME) ; AUTHOR("P Molfese") ; THD_check_AFNI_version(PROGRAM_NAME) ;

   while( nopt<argc && argv[nopt][0]=='-' ){
#define IS(s) (strcasecmp(argv[nopt],(s))==0)
#define NEEDARG(what) do{ if(++nopt>=argc)ERROR_exit(PROGRAM_NAME ": need %s after '%s'",(what),argv[nopt-1]); }while(0)
     if( IS("-help") || IS("-h") ){ usage_3dMVPA() ; return 0 ; }
     else if( IS("-mode") ){ NEEDARG("classify") ; mode=argv[nopt++] ; }
     else if( IS("-target_column") ){ NEEDARG("a column name") ; target_column=argv[nopt++] ; }
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
     else if( IS("-method") ){ NEEDARG("centroid") ; method=argv[nopt++] ; }
     else if( IS("-distance") ){ NEEDARG("euclidean, correlation, or cosine") ; distance_name=argv[nopt++] ; }
     else if( IS("-scale") ){ NEEDARG("none, feature_zscore, sample_center, or run_center") ; scale_name=argv[nopt++] ; }
     else if( IS("-metric") ){ NEEDARG("accuracy or balanced_accuracy") ; metric_name=argv[nopt++] ; }
     else if( IS("-nperm") ){ NEEDARG("a nonnegative count") ; nperm=(int)mvpa_parse_long("-nperm",argv[nopt++],0,INT_MAX) ; }
     else if( IS("-seed") ){ NEEDARG("an integer seed") ; seed=mvpa_parse_long("-seed",argv[nopt++],LONG_MIN,LONG_MAX) ; }
     else if( IS("-permute_within") || IS("-block_column") ){
       NEEDARG("an exchangeability-block column or none") ;
       if( permute_within!=NULL ) ERROR_exit(PROGRAM_NAME ": the permutation block was given twice") ;
       permute_within=argv[nopt++] ;
     }
     else if( IS("-prefix") ){ NEEDARG("an output prefix") ; prefix=argv[nopt++] ; }
     else if( IS("-save_predictions") ){ NEEDARG("a filename") ; save_predictions=argv[nopt++] ; }
     else if( IS("-save_fold_scores") ){ NEEDARG("a filename") ; save_fold_scores=argv[nopt++] ; }
     else if( IS("-show_table") ){ show_table=1 ; nopt++ ; }
     else if( IS("-show_samples") ){ show_samples=1 ; nopt++ ; }
     else if( IS("-show_folds") ){ show_folds=1 ; nopt++ ; }
     else if( IS("-quiet") ){ quiet=1 ; nopt++ ; }
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
   if( mode==NULL || strcasecmp(mode,"classify")!=0 )
     ERROR_exit(PROGRAM_NAME ": this release slice requires '-mode classify'") ;
   if( tab==NULL ) ERROR_exit(PROGRAM_NAME ": -dataTable or -dataTableFile is required") ;
   if( target_column==NULL ) ERROR_exit(PROGRAM_NAME ": -target_column is required") ;
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
   if( strcasecmp(method,"centroid")!=0 )
     ERROR_exit(PROGRAM_NAME ": only '-method centroid' is implemented") ;
   if( !THD_filename_ok(prefix) ) ERROR_exit(PROGRAM_NAME ": illegal -prefix '%s'",prefix) ;
   if( !THD_centroid_distance_parse(distance_name,&distance,err,sizeof(err)) )
     ERROR_exit(PROGRAM_NAME ": %s",err) ;
   if( !THD_preprocess_parse(scale_name,&preprocess,err,sizeof(err)) )
     ERROR_exit(PROGRAM_NAME ": %s",err) ;
   if( !THD_mvpa_metric_parse(metric_name,&metric,err,sizeof(err)) )
     ERROR_exit(PROGRAM_NAME ": %s",err) ;
   if( metric!=THD_METRIC_ACCURACY && metric!=THD_METRIC_BALANCED_ACCURACY )
     ERROR_exit(PROGRAM_NAME ": classification -metric must be accuracy or balanced_accuracy") ;

   if( THD_datatable_column(tab,"Xmat")>=0 ){
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

   if( !THD_sampleset_encode_column(ss,target_column,&class_id,&class_label,&nclass,err,sizeof(err)) )
     ERROR_exit(PROGRAM_NAME ": %s",err) ;
   if( nclass<2 ) ERROR_exit(PROGRAM_NAME ": target column '%s' has only %d class",target_column,nclass) ;
   plan=(cv_k>0) ? THD_foldplan_grouped_kfold(ss,cv_column,within_column,cv_k,err,sizeof(err))
                 : THD_foldplan_leave_group_out(ss,cv_column,within_column,err,sizeof(err)) ;
   if( plan==NULL ) ERROR_exit(PROGRAM_NAME ": %s",err) ;
   if( !THD_foldplan_check_classes(plan,class_id,nclass,0,err,sizeof(err)) )
     ERROR_exit(PROGRAM_NAME ": %s",err) ;
   if( show_folds ) THD_foldplan_print(plan,ss,stdout) ;

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
     if( !THD_perm_scheme_set_eqclass(perm_scheme,class_id) )
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
     if( !quiet )
       INFO_message(PROGRAM_NAME ": %d synchronized label relabelings (%s; %d block%s from %s)",
                    perm_set->nperm,perm_set->is_exact?"exact":"Monte Carlo",
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
   if( surf_file!=NULL ){
     for( ii=0 ; ii<ss->ndset ; ii++ )
       if( !mvpa_surface_domains_match(mset,ss->dset[ii],err,sizeof(err)) )
         ERROR_exit(PROGRAM_NAME ": %s",err) ;
   } else if( !EQUIV_GRIDS(mset,ss->dset[0]) )
     ERROR_exit(PROGRAM_NAME ": mask '%s' is not on the sample dataset grid",maskname) ;
   DSET_load(mset) ;
   if( !DSET_LOADED(mset) )
     ERROR_exit(PROGRAM_NAME ": cannot load %s dataset '%s'",
                maskname!=NULL?"mask":"surface-domain",maskname!=NULL?maskname:ss->dset_name[0]) ;
   if( searchlight!=NULL ){
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
       if( mvpa_is_surface_dset(mset) )
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
   maxvox=THD_roilist_maxvox(rl) ;
   if( maxvox<1 || (size_t)ss->nsample>SIZE_MAX/sizeof(float)/(size_t)maxvox )
     ERROR_exit(PROGRAM_NAME ": sample-by-voxel workspace is too large") ;
#ifdef USE_OMP
   if( searchlight!=NULL ) nworker=omp_get_max_threads() ;
#endif
   if( nworker>rl->nroi ) nworker=rl->nroi ;
   for( ii=0 ; ii<ss->ndset ; ii++ ){
     DSET_load(ss->dset[ii]) ;
     if( !DSET_LOADED(ss->dset[ii]) )
       ERROR_exit(PROGRAM_NAME ": cannot load sample dataset '%s'",ss->dset_name[ii]) ;
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
     fputs("UnitType\tUnitID\tUnitLabel\tI\tJ\tK\tNFeat\tTableRow\tSample\tSubj\tRun\tItem\tInputFold\tCVFold\tCVStratum\tObserved\tPredicted\tCorrect",prediction_fp) ;
     for( ii=0 ; ii<nclass ; ii++ ) fprintf(prediction_fp,"\tDecisionClass%d",ii) ;
     fputc('\n',prediction_fp) ;
   }
   if( save_fold_scores!=NULL ){
     fold_fp=mvpa_open_audit(save_fold_scores,"fold-score audit") ;
     mvpa_write_audit_header(fold_fp,"per-fold classification scores",
                             target_column,class_label,nclass,searchlight!=NULL,
                             surf_file!=NULL) ;
     fputs("UnitType\tUnitID\tUnitLabel\tI\tJ\tK\tNFeat\tCVFold\tStratum\tNTrain\tNTest\tACC\tBACC\n",fold_fp) ;
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
       ok=THD_mvpa_run_centroid_permutations(plan,distance,x[tid],rl->vox[rr].nar,
                                             class_id,group_id,metric,perm_set,
                                             perm_class[tid],perm_score[tid],
                                             workspace[tid],localerr,sizeof(localerr)) ;
     else if( ok )
       ok=THD_mvpa_run_centroid_cv(plan,distance,x[tid],rl->vox[rr].nar,
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

   DSET_delete(oset) ; DSET_delete(mset) ; THD_roilist_free(rl) ;
   for( ii=0 ; ii<nworker ; ii++ ){
     THD_mvpa_workspace_free(workspace[ii]) ; free(x[ii]) ;
     if( perm_set!=NULL ){
       free(perm_class[ii]) ; free(perm_score[ii]) ; free(perm_max[ii]) ;
     }
   }
   free(workspace) ; free(x) ; THD_foldplan_free(plan) ;
   free(perm_class) ; free(perm_score) ; free(perm_max) ;
   THD_perm_result_free(perm_result) ; THD_perm_set_free(perm_set) ;
   THD_sampleset_free_levels(class_label,nclass) ; free(class_id) ; free(group_id) ;
   THD_sampleset_free(ss) ; free(score) ; free(acc) ; free(bacc) ;
   free(nfeat) ; free(nfold) ; free(qval) ; free(table_name) ;
   return 0 ;
}
