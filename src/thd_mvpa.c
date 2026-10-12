#include "thd_mvpa.h"

#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <strings.h>

static int pp_error( char *err, int nerr, const char *fmt, ... )
{
   va_list ap ;
   if( err!=NULL && nerr>0 ){
     va_start(ap,fmt) ; vsnprintf(err,(size_t)nerr,fmt,ap) ; va_end(ap) ;
   }
   return 0 ;
}

const char * THD_mvpa_metric_name( int metric )
{
   switch( metric ){
     case THD_METRIC_ACCURACY:           return "accuracy" ;
     case THD_METRIC_BALANCED_ACCURACY:  return "balanced_accuracy" ;
     case THD_METRIC_PEARSON:            return "pearson" ;
     case THD_METRIC_SPEARMAN:           return "spearman" ;
     case THD_METRIC_CROSS_VALIDATED_R2: return "cross_validated_r2" ;
     case THD_METRIC_MAE:                return "mae" ;
     case THD_METRIC_RMSE:               return "rmse" ;
   }
   return NULL ;
}

int THD_mvpa_metric_parse( const char *name, int *metric, char *err, int nerr )
{
   int mm ;
   if( err!=NULL && nerr>0 ) err[0]='\0' ;
   if( name==NULL || metric==NULL )
     return pp_error(err,nerr,"metric parsing needs a name and output") ;
   for( mm=THD_METRIC_ACCURACY ; mm<=THD_METRIC_RMSE ; mm++ )
     if( strcasecmp(name,THD_mvpa_metric_name(mm))==0 ){ *metric=mm ; return 1 ; }
   return pp_error(err,nerr,"unknown MVPA metric '%s'",name) ;
}

int THD_mvpa_score_classification( const int *truth, const int *prediction,
                                    int nsample, int nclass,
                                    int *confusion, double *recall,
                                    THD_mvpa_class_score *score,
                                    char *err, int nerr )
{
   int ii,cc,ncorrect=0,npresent=0 ; double recall_sum=0.0 ;
   if( err!=NULL && nerr>0 ) err[0]='\0' ;
   if( truth==NULL || prediction==NULL || confusion==NULL || recall==NULL ||
       score==NULL || nsample<1 || nclass<2 )
     return pp_error(err,nerr,"classification scoring needs samples, at least two classes, and output storage") ;
   if( (size_t)nclass > SIZE_MAX/sizeof(int)/(size_t)nclass )
     return pp_error(err,nerr,"classification confusion matrix is too large") ;
   memset(confusion,0,sizeof(int)*(size_t)nclass*nclass) ;
   for( ii=0 ; ii<nsample ; ii++ ){
     int tt=truth[ii],pp=prediction[ii] ;
     if( tt<0 || tt>=nclass )
       return pp_error(err,nerr,"true class %d at sample %d is outside [0,%d]",tt,ii+1,nclass-1) ;
     if( pp<0 || pp>=nclass )
       return pp_error(err,nerr,"predicted class %d at sample %d is outside [0,%d]",pp,ii+1,nclass-1) ;
     confusion[(size_t)tt*nclass+pp]++ ;
     if( tt==pp ) ncorrect++ ;
   }
   for( cc=0 ; cc<nclass ; cc++ ){
     int pp,total=0 ;
     for( pp=0 ; pp<nclass ; pp++ ) total+=confusion[(size_t)cc*nclass+pp] ;
     if( total==0 ) recall[cc]=NAN ;
     else{
       recall[cc]=(double)confusion[(size_t)cc*nclass+cc]/(double)total ;
       recall_sum+=recall[cc] ; npresent++ ;
     }
   }
   score->nsample=nsample ; score->nclass=nclass ; score->ncorrect=ncorrect ;
   score->nclass_present=npresent ;
   score->accuracy=(double)ncorrect/(double)nsample ;
   score->balanced_accuracy=npresent>0 ? recall_sum/(double)npresent : NAN ;
   return 1 ;
}

static void mvpa_rank_values( const float *value, int nvalue, double *rank,
                              int *order, int *temp )
{
   int ii,width ;
   for( ii=0 ; ii<nvalue ; ii++ ) order[ii]=ii ;
   for( width=1 ; width<nvalue ; width*=2 ){
     int start ;
     for( start=0 ; start<nvalue ; start+=2*width ){
       int aa=start,bb=start+width,end=start+2*width,out=start ;
       if( bb>nvalue ) bb=nvalue ;
       if( end>nvalue ) end=nvalue ;
       {
         int amid=bb ;
         while( aa<amid && bb<end ){
           int ia=order[aa],ib=order[bb] ;
           if( value[ia]<value[ib] || (value[ia]==value[ib] && ia<ib) ) temp[out++]=order[aa++] ;
           else temp[out++]=order[bb++] ;
         }
         while( aa<amid ) temp[out++]=order[aa++] ;
         while( bb<end ) temp[out++]=order[bb++] ;
       }
     }
     for( ii=0 ; ii<nvalue ; ii++ ) order[ii]=temp[ii] ;
     if( width>nvalue/2 ) break ;
   }
   for( ii=0 ; ii<nvalue ; ){
     int jj=ii+1 ; double average ;
     while( jj<nvalue && value[order[jj]]==value[order[ii]] ) jj++ ;
     average=0.5*((double)ii+1.0+(double)jj) ;
     while( ii<jj ){ rank[order[ii]]=average ; ii++ ; }
   }
}

static double mvpa_corr_float( const float *aa, const float *bb, int nn )
{
   int ii ; double ma=0.0,mb=0.0,ab=0.0,a2=0.0,b2=0.0 ;
   for( ii=0 ; ii<nn ; ii++ ){ ma+=aa[ii] ; mb+=bb[ii] ; }
   ma/=(double)nn ; mb/=(double)nn ;
   for( ii=0 ; ii<nn ; ii++ ){
     double da=(double)aa[ii]-ma,db=(double)bb[ii]-mb ;
     ab+=da*db ; a2+=da*da ; b2+=db*db ;
   }
   return (a2>0.0 && b2>0.0) ? ab/sqrt(a2*b2) : NAN ;
}

static double mvpa_corr_double( const double *aa, const double *bb, int nn )
{
   int ii ; double ma=0.0,mb=0.0,ab=0.0,a2=0.0,b2=0.0 ;
   for( ii=0 ; ii<nn ; ii++ ){ ma+=aa[ii] ; mb+=bb[ii] ; }
   ma/=(double)nn ; mb/=(double)nn ;
   for( ii=0 ; ii<nn ; ii++ ){
     double da=aa[ii]-ma,db=bb[ii]-mb ;
     ab+=da*db ; a2+=da*da ; b2+=db*db ;
   }
   return (a2>0.0 && b2>0.0) ? ab/sqrt(a2*b2) : NAN ;
}

int THD_mvpa_score_continuous( const float *truth, const float *prediction,
                                const float *baseline, int nsample,
                                double *rank_scratch, int *index_scratch,
                                THD_mvpa_continuous_score *score,
                                char *err, int nerr )
{
   int ii ; double sae=0.0,sse=0.0,baseline_sse=0.0 ;
   if( err!=NULL && nerr>0 ) err[0]='\0' ;
   if( truth==NULL || prediction==NULL || baseline==NULL || rank_scratch==NULL ||
       index_scratch==NULL || score==NULL || nsample<1 )
     return pp_error(err,nerr,"continuous scoring needs nonempty inputs and scratch/output storage") ;
   for( ii=0 ; ii<nsample ; ii++ ){
     double ee,bb ;
     if( !isfinite(truth[ii]) || !isfinite(prediction[ii]) || !isfinite(baseline[ii]) )
       return pp_error(err,nerr,"continuous score input at sample %d is not finite",ii+1) ;
     ee=(double)truth[ii]-prediction[ii] ;
     bb=(double)truth[ii]-baseline[ii] ;
     sae+=fabs(ee) ; sse+=ee*ee ; baseline_sse+=bb*bb ;
   }
   mvpa_rank_values(truth,nsample,rank_scratch,index_scratch,index_scratch+nsample) ;
   mvpa_rank_values(prediction,nsample,rank_scratch+nsample,index_scratch,index_scratch+nsample) ;
   score->nsample=nsample ;
   score->pearson=mvpa_corr_float(truth,prediction,nsample) ;
   score->spearman=mvpa_corr_double(rank_scratch,rank_scratch+nsample,nsample) ;
   score->cv_r2=baseline_sse>0.0 ? 1.0-sse/baseline_sse : NAN ;
   score->mae=sae/(double)nsample ; score->rmse=sqrt(sse/(double)nsample) ;
   return 1 ;
}

void THD_preprocess_free( THD_preprocessor *p )
{
   if( p==NULL ) return ;
   free(p->group_value) ; free(p->center) ; free(p->scale) ;
   free(p->constant) ; free(p) ;
}

static THD_preprocessor * pp_fail( THD_preprocessor *p, char *err, int nerr,
                                   const char *fmt, ... )
{
   va_list ap ;
   if( err!=NULL && nerr>0 ){
     va_start(ap,fmt) ; vsnprintf(err,(size_t)nerr,fmt,ap) ; va_end(ap) ;
   }
   THD_preprocess_free(p) ; return NULL ;
}

const char * THD_preprocess_name( int method )
{
   switch( method ){
     case THD_PREP_NONE:           return "none" ;
     case THD_PREP_FEATURE_ZSCORE: return "feature_zscore" ;
     case THD_PREP_SAMPLE_CENTER:  return "sample_center" ;
     case THD_PREP_RUN_CENTER:     return "run_center" ;
   }
   return NULL ;
}

int THD_preprocess_parse( const char *name, int *method, char *err, int nerr )
{
   int mm ;
   if( err!=NULL && nerr>0 ) err[0]='\0' ;
   if( name==NULL || method==NULL )
     return pp_error(err,nerr,"preprocessing needs a method name and output") ;
   for( mm=THD_PREP_NONE ; mm<=THD_PREP_RUN_CENTER ; mm++ )
     if( strcasecmp(name,THD_preprocess_name(mm))==0 ){ *method=mm ; return 1 ; }
   return pp_error(err,nerr,"unknown preprocessing method '%s'",name) ;
}

static int pp_check_rows( const float *x, int nsample, int nfeature,
                          const int *index, int nindex, int **seen_out,
                          char *err, int nerr )
{
   int ii,jj,*seen ;
   if( x==NULL || nsample<1 || nfeature<1 || index==NULL || nindex<1 )
     return pp_error(err,nerr,"preprocessing needs a nonempty matrix and sample list") ;
   seen=(int *)calloc((size_t)nsample,sizeof(int)) ;
   if( seen==NULL ) return pp_error(err,nerr,"out of memory validating preprocessing rows") ;
   for( ii=0 ; ii<nindex ; ii++ ){
     int ss=index[ii] ;
     if( ss<0 || ss>=nsample || seen[ss] ){
       free(seen) ; return pp_error(err,nerr,"invalid or duplicate preprocessing sample index %d",ss) ;
     }
     seen[ss]=1 ;
     for( jj=0 ; jj<nfeature ; jj++ ) if( !isfinite(x[(size_t)ss*nfeature+jj]) ){
       free(seen) ; return pp_error(err,nerr,"sample %d, feature %d is not finite",ss+1,jj+1) ;
     }
   }
   if( seen_out!=NULL ) *seen_out=seen ; else free(seen) ;
   return 1 ;
}

static int pp_group_index( const THD_preprocessor *p, int value )
{
   int gg ;
   for( gg=0 ; gg<p->ngroup ; gg++ ) if( p->group_value[gg]==value ) return gg ;
   return -1 ;
}

THD_preprocessor * THD_preprocess_fit( int method, const float *x,
                                        int nsample, int nfeature,
                                        const int *fit_index, int nfit,
                                        const int *group_id,
                                        char *err, int nerr )
{
   THD_preprocessor *p ; int ii,jj,*seen=NULL ;
   if( err!=NULL && nerr>0 ) err[0]='\0' ;
   if( THD_preprocess_name(method)==NULL ){
     pp_error(err,nerr,"invalid preprocessing method %d",method) ; return NULL ;
   }
   if( !pp_check_rows(x,nsample,nfeature,fit_index,nfit,&seen,err,nerr) ) return NULL ;
   free(seen) ;
   if( method==THD_PREP_RUN_CENTER && group_id==NULL ){
     pp_error(err,nerr,"run_center needs one group ID per sample") ; return NULL ;
   }
   p=(THD_preprocessor *)calloc(1,sizeof(*p)) ;
   if( p==NULL ){ pp_error(err,nerr,"out of memory fitting preprocessing") ; return NULL ; }
   p->method=method ; p->nfeature=nfeature ; p->nfit=nfit ;

   if( method==THD_PREP_NONE || method==THD_PREP_SAMPLE_CENTER ) return p ;

   if( method==THD_PREP_FEATURE_ZSCORE ){
     p->ngroup=1 ;
     p->center=(double *)calloc((size_t)nfeature,sizeof(double)) ;
     p->scale=(double *)malloc(sizeof(double)*(size_t)nfeature) ;
     p->constant=(unsigned char *)calloc((size_t)nfeature,sizeof(unsigned char)) ;
     if( p->center==NULL || p->scale==NULL || p->constant==NULL )
       return pp_fail(p,err,nerr,"out of memory fitting feature_zscore") ;
     for( ii=0 ; ii<nfit ; ii++ ){
       const float *row=x+(size_t)fit_index[ii]*nfeature ;
       for( jj=0 ; jj<nfeature ; jj++ ) p->center[jj]+=row[jj] ;
     }
     for( jj=0 ; jj<nfeature ; jj++ ) p->center[jj]/=(double)nfit ;
     for( jj=0 ; jj<nfeature ; jj++ ){
       double ss=0.0,tol ;
       for( ii=0 ; ii<nfit ; ii++ ){
         double d=(double)x[(size_t)fit_index[ii]*nfeature+jj]-p->center[jj] ; ss+=d*d ;
       }
       p->scale[jj]=sqrt(ss/(double)nfit) ;
       tol=DBL_EPSILON*fmax(1.0,fabs(p->center[jj])) ;
       if( p->scale[jj]<=tol ){
         p->constant[jj]=1 ; p->nconstant++ ; p->scale[jj]=1.0 ;
       }
     }
     return p ;
   }

   p->group_value=(int *)malloc(sizeof(int)*(size_t)nfit) ;
   if( p->group_value==NULL ) return pp_fail(p,err,nerr,"out of memory fitting run_center") ;
   for( ii=0 ; ii<nfit ; ii++ ){
     int value=group_id[fit_index[ii]],gg ;
     if( value<0 ) return pp_fail(p,err,nerr,"sample %d has no valid run/group ID",fit_index[ii]+1) ;
     for( gg=0 ; gg<p->ngroup ; gg++ ) if( p->group_value[gg]==value ) break ;
     if( gg==p->ngroup ) p->group_value[p->ngroup++]=value ;
   }
   if( (size_t)p->ngroup > SIZE_MAX/sizeof(double)/(size_t)nfeature )
     return pp_fail(p,err,nerr,"run_center dimensions are too large") ;
   p->center=(double *)calloc((size_t)p->ngroup*nfeature,sizeof(double)) ;
   {
     int *count=(int *)calloc((size_t)p->ngroup,sizeof(int)) ;
     if( p->center==NULL || count==NULL ){
       free(count) ; return pp_fail(p,err,nerr,"out of memory fitting run_center") ;
     }
     for( ii=0 ; ii<nfit ; ii++ ){
       int ss=fit_index[ii],gg=pp_group_index(p,group_id[ss]) ;
       count[gg]++ ;
       for( jj=0 ; jj<nfeature ; jj++ ) p->center[(size_t)gg*nfeature+jj]+=x[(size_t)ss*nfeature+jj] ;
     }
     for( ii=0 ; ii<p->ngroup ; ii++ ) for( jj=0 ; jj<nfeature ; jj++ )
       p->center[(size_t)ii*nfeature+jj]/=(double)count[ii] ;
     free(count) ;
   }
   return p ;
}

int THD_preprocess_apply( const THD_preprocessor *p, const float *x,
                           int nsample, int nfeature,
                           const int *apply_index, int napply,
                           const int *group_id, float *out,
                           char *err, int nerr )
{
   int ii,jj ;
   if( err!=NULL && nerr>0 ) err[0]='\0' ;
   if( p==NULL || out==NULL || p->nfeature!=nfeature )
     return pp_error(err,nerr,"invalid fitted preprocessing object or output shape") ;
   if( out==x ) return pp_error(err,nerr,"preprocessing output must not alias its input matrix") ;
   if( p->method==THD_PREP_RUN_CENTER && group_id==NULL )
     return pp_error(err,nerr,"run_center needs one group ID per sample") ;
   if( !pp_check_rows(x,nsample,nfeature,apply_index,napply,NULL,err,nerr) ) return 0 ;
   for( ii=0 ; ii<napply ; ii++ ){
     int ss=apply_index[ii],gg=-1 ; const float *row=x+(size_t)ss*nfeature ;
     float *dest=out+(size_t)ii*nfeature ; double rowmean=0.0 ;
     if( p->method==THD_PREP_RUN_CENTER ){
       gg=pp_group_index(p,group_id[ss]) ;
       if( gg<0 ) return pp_error(err,nerr,"sample %d has run/group %d, which was not present in training",ss+1,group_id[ss]) ;
     } else if( p->method==THD_PREP_SAMPLE_CENTER ){
       for( jj=0 ; jj<nfeature ; jj++ ) rowmean+=row[jj] ;
       rowmean/=(double)nfeature ;
     }
     for( jj=0 ; jj<nfeature ; jj++ ){
       double value ;
       if( p->method==THD_PREP_NONE ) value=row[jj] ;
       else if( p->method==THD_PREP_FEATURE_ZSCORE )
         value=p->constant[jj]?0.0:((double)row[jj]-p->center[jj])/p->scale[jj] ;
       else if( p->method==THD_PREP_SAMPLE_CENTER ) value=(double)row[jj]-rowmean ;
       else value=(double)row[jj]-p->center[(size_t)gg*nfeature+jj] ;
       if( !isfinite(value) || fabs(value)>FLT_MAX )
         return pp_error(err,nerr,"sample %d, feature %d transforms outside the float range",ss+1,jj+1) ;
       dest[jj]=(float)value ;
     }
   }
   return 1 ;
}

int THD_preprocess_fold( const THD_foldplan *plan, int fold, int method,
                          const float *x, int nfeature, const int *group_id,
                          float *train_out, float *test_out,
                          THD_preprocessor **fitted, char *err, int nerr )
{
   THD_preprocessor *p ;
   if( err!=NULL && nerr>0 ) err[0]='\0' ;
   if( fitted!=NULL ) *fitted=NULL ;
   if( !THD_foldplan_validate(plan,err,nerr) ) return 0 ;
   if( fold<0 || fold>=plan->nfold ) return pp_error(err,nerr,"invalid fold index %d",fold) ;
   if( train_out==NULL || test_out==NULL ) return pp_error(err,nerr,"fold preprocessing needs train and test outputs") ;
   p=THD_preprocess_fit(method,x,plan->nsample,nfeature,plan->train[fold],
                        plan->ntrain[fold],group_id,err,nerr) ;
   if( p==NULL ) return 0 ;
   if( !THD_preprocess_apply(p,x,plan->nsample,nfeature,plan->train[fold],
                             plan->ntrain[fold],group_id,train_out,err,nerr) ||
       !THD_preprocess_apply(p,x,plan->nsample,nfeature,plan->test[fold],
                             plan->ntest[fold],group_id,test_out,err,nerr) ){
     THD_preprocess_free(p) ; return 0 ;
   }
   if( fitted!=NULL ) *fitted=p ; else THD_preprocess_free(p) ;
   return 1 ;
}
