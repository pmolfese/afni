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

void THD_centroid_free( THD_centroid_model *m )
{
   if( m==NULL ) return ;
   free(m->class_count) ; free(m->centroid) ;
   free(m->centroid_mean) ; free(m->centroid_norm) ; free(m) ;
}

static THD_centroid_model * centroid_alloc( int nclass, int max_feature,
                                            char *err, int nerr )
{
   THD_centroid_model *m ;
   if( nclass<2 || max_feature<1 ||
       (size_t)nclass > SIZE_MAX/sizeof(double)/(size_t)max_feature ){
     pp_error(err,nerr,"invalid or excessive centroid dimensions") ; return NULL ;
   }
   m=(THD_centroid_model *)calloc(1,sizeof(*m)) ;
   if( m==NULL ){ pp_error(err,nerr,"out of memory allocating centroid model") ; return NULL ; }
   m->nclass=nclass ; m->nfeature_alloc=max_feature ;
   m->class_count=(int *)calloc((size_t)nclass,sizeof(int)) ;
   m->centroid=(double *)calloc((size_t)nclass*max_feature,sizeof(double)) ;
   m->centroid_mean=(double *)calloc((size_t)nclass,sizeof(double)) ;
   m->centroid_norm=(double *)calloc((size_t)nclass,sizeof(double)) ;
   if( m->class_count==NULL || m->centroid==NULL ||
       m->centroid_mean==NULL || m->centroid_norm==NULL ){
     THD_centroid_free(m) ;
     pp_error(err,nerr,"out of memory allocating centroid model") ; return NULL ;
   }
   return m ;
}

const char * THD_centroid_distance_name( int distance )
{
   switch( distance ){
     case THD_CENTROID_EUCLIDEAN:   return "euclidean" ;
     case THD_CENTROID_CORRELATION: return "correlation" ;
     case THD_CENTROID_COSINE:      return "cosine" ;
   }
   return NULL ;
}

int THD_centroid_distance_parse( const char *name, int *distance,
                                 char *err, int nerr )
{
   int dd ;
   if( err!=NULL && nerr>0 ) err[0]='\0' ;
   if( name==NULL || distance==NULL )
     return pp_error(err,nerr,"centroid distance parsing needs a name and output") ;
   for( dd=THD_CENTROID_EUCLIDEAN ; dd<=THD_CENTROID_COSINE ; dd++ )
     if( strcasecmp(name,THD_centroid_distance_name(dd))==0 ){
       *distance=dd ; return 1 ;
     }
   return pp_error(err,nerr,"unknown centroid distance '%s'",name) ;
}

static int centroid_refit( THD_centroid_model *m, int distance,
                           const float *x, const int *class_id,
                           int nsample, int nfeature, char *err, int nerr )
{
   int ii,jj,cc ;
   if( err!=NULL && nerr>0 ) err[0]='\0' ;
   if( m==NULL || m->class_count==NULL || m->centroid==NULL ||
       m->centroid_mean==NULL || m->centroid_norm==NULL ||
       THD_centroid_distance_name(distance)==NULL )
     return pp_error(err,nerr,"invalid centroid model or distance") ;
   if( x==NULL || class_id==NULL || nsample<1 || nfeature<1 ||
       nfeature>m->nfeature_alloc )
     return pp_error(err,nerr,"centroid fitting matrix exceeds the model workspace") ;
   if( m->nclass>nsample )
     return pp_error(err,nerr,"centroid fitting needs at least one sample in every class") ;
   if( distance==THD_CENTROID_CORRELATION && nfeature<2 )
     return pp_error(err,nerr,"correlation distance needs at least two features") ;

   m->distance=distance ; m->nfeature=nfeature ; m->nfit=nsample ;
   memset(m->class_count,0,sizeof(int)*(size_t)m->nclass) ;
   memset(m->centroid,0,sizeof(double)*(size_t)m->nclass*nfeature) ;
   memset(m->centroid_mean,0,sizeof(double)*(size_t)m->nclass) ;
   memset(m->centroid_norm,0,sizeof(double)*(size_t)m->nclass) ;
   for( ii=0 ; ii<nsample ; ii++ ){
     cc=class_id[ii] ;
     if( cc<0 || cc>=m->nclass )
       return pp_error(err,nerr,"training class %d at sample %d is outside [0,%d]",cc,ii+1,m->nclass-1) ;
     m->class_count[cc]++ ;
     for( jj=0 ; jj<nfeature ; jj++ ){
       float vv=x[(size_t)ii*nfeature+jj] ;
       if( !isfinite(vv) )
         return pp_error(err,nerr,"training sample %d, feature %d is not finite",ii+1,jj+1) ;
       m->centroid[(size_t)cc*nfeature+jj]+=(double)vv ;
     }
   }
   for( cc=0 ; cc<m->nclass ; cc++ ){
     if( m->class_count[cc]==0 )
       return pp_error(err,nerr,"training data have no samples for class %d",cc) ;
     for( jj=0 ; jj<nfeature ; jj++ )
       m->centroid[(size_t)cc*nfeature+jj]/=(double)m->class_count[cc] ;
   }
   if( distance==THD_CENTROID_EUCLIDEAN ) return 1 ;
   for( cc=0 ; cc<m->nclass ; cc++ ){
     const double *cent=m->centroid+(size_t)cc*nfeature ;
     double mean=0.0,ss=0.0,rawss=0.0 ;
     if( distance==THD_CENTROID_CORRELATION ){
       for( jj=0 ; jj<nfeature ; jj++ ) mean+=cent[jj] ;
       mean/=(double)nfeature ; m->centroid_mean[cc]=mean ;
     }
     for( jj=0 ; jj<nfeature ; jj++ ){
       double vv=cent[jj],dd=vv-mean ; ss+=dd*dd ; rawss+=vv*vv ;
     }
     if( ss<=DBL_EPSILON*rawss ){
       if( distance==THD_CENTROID_CORRELATION )
         return pp_error(err,nerr,"class %d centroid is constant for correlation distance",cc) ;
       return pp_error(err,nerr,"class %d centroid has zero norm for cosine distance",cc) ;
     }
     m->centroid_norm[cc]=sqrt(ss) ;
   }
   return 1 ;
}

THD_centroid_model * THD_centroid_fit( int distance, const float *x,
                                        const int *class_id,
                                        int nsample, int nfeature, int nclass,
                                        char *err, int nerr )
{
   THD_centroid_model *m=centroid_alloc(nclass,nfeature,err,nerr) ;
   if( m==NULL ) return NULL ;
   if( !centroid_refit(m,distance,x,class_id,nsample,nfeature,err,nerr) ){
     THD_centroid_free(m) ; return NULL ;
   }
   return m ;
}

int THD_centroid_predict( const THD_centroid_model *m,
                          const float *x, int nsample, int nfeature,
                          int *prediction, double *decision,
                          char *err, int nerr )
{
   int ii,jj,cc ;
   if( err!=NULL && nerr>0 ) err[0]='\0' ;
   if( m==NULL || m->centroid==NULL || m->class_count==NULL || x==NULL ||
       prediction==NULL || nsample<1 || nfeature!=m->nfeature || m->nclass<2 )
     return pp_error(err,nerr,"invalid centroid model, prediction matrix, or output shape") ;
   if( m->distance!=THD_CENTROID_EUCLIDEAN &&
       (m->centroid_mean==NULL || m->centroid_norm==NULL) )
     return pp_error(err,nerr,"centroid model is missing normalization statistics") ;
   for( ii=0 ; ii<nsample ; ii++ ){
     const float *row=x+(size_t)ii*nfeature ;
     double mean=0.0,norm=0.0,best=0.0 ; int best_class=0 ;
     for( jj=0 ; jj<nfeature ; jj++ ){
       if( !isfinite(row[jj]) )
         return pp_error(err,nerr,"prediction sample %d, feature %d is not finite",ii+1,jj+1) ;
       if( m->distance==THD_CENTROID_CORRELATION ) mean+=row[jj] ;
     }
     if( m->distance==THD_CENTROID_CORRELATION ) mean/=(double)nfeature ;
     if( m->distance!=THD_CENTROID_EUCLIDEAN ){
       double rawss=0.0 ;
       for( jj=0 ; jj<nfeature ; jj++ ){
         double vv=row[jj],dd=vv-mean ; norm+=dd*dd ; rawss+=vv*vv ;
       }
       if( norm<=DBL_EPSILON*rawss ){
         if( m->distance==THD_CENTROID_CORRELATION )
           return pp_error(err,nerr,"prediction sample %d is constant for correlation distance",ii+1) ;
         return pp_error(err,nerr,"prediction sample %d has zero norm for cosine distance",ii+1) ;
       }
       norm=sqrt(norm) ;
     }
     for( cc=0 ; cc<m->nclass ; cc++ ){
       const double *cent=m->centroid+(size_t)cc*nfeature ; double value=0.0 ;
       if( m->distance==THD_CENTROID_EUCLIDEAN ){
         for( jj=0 ; jj<nfeature ; jj++ ){
           double dd=(double)row[jj]-cent[jj] ; value-=dd*dd ;
         }
       } else {
         double cmean=m->centroid_mean[cc] ;
         for( jj=0 ; jj<nfeature ; jj++ )
           value+=((double)row[jj]-mean)*(cent[jj]-cmean) ;
         value/=norm*m->centroid_norm[cc] ;
         if( value>1.0 ) value=1.0 ; else if( value < -1.0 ) value=-1.0 ;
       }
       if( !isfinite(value) )
         return pp_error(err,nerr,"decision score is not finite for sample %d, class %d",ii+1,cc) ;
       if( decision!=NULL ) decision[(size_t)ii*m->nclass+cc]=value ;
       if( cc==0 || value>best ){ best=value ; best_class=cc ; }
     }
     prediction[ii]=best_class ;
   }
   return 1 ;
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

static void mvpa_class_score_from_confusion( const int *confusion,
                                              int nsample, int nclass,
                                              double *recall,
                                              THD_mvpa_class_score *score )
{
   int cc,ncorrect=0,npresent=0 ; double recall_sum=0.0 ;
   for( cc=0 ; cc<nclass ; cc++ ){
     int pp,total=0 ;
     for( pp=0 ; pp<nclass ; pp++ ) total+=confusion[(size_t)cc*nclass+pp] ;
     ncorrect+=confusion[(size_t)cc*nclass+cc] ;
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
}

int THD_mvpa_score_classification( const int *truth, const int *prediction,
                                    int nsample, int nclass,
                                    int *confusion, double *recall,
                                    THD_mvpa_class_score *score,
                                    char *err, int nerr )
{
   int ii ;
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
   }
   mvpa_class_score_from_confusion(confusion,nsample,nclass,recall,score) ;
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
   free(p->group_value) ; free(p->group_count) ;
   free(p->center) ; free(p->scale) ;
   free(p->constant) ; free(p) ;
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

static THD_preprocessor * pp_alloc( int max_feature, int max_group,
                                    char *err, int nerr )
{
   THD_preprocessor *p ;
   if( max_feature<1 || max_group<1 ||
       (size_t)max_group > SIZE_MAX/sizeof(double)/(size_t)max_feature ){
     pp_error(err,nerr,"invalid or excessive preprocessing workspace dimensions") ;
     return NULL ;
   }
   p=(THD_preprocessor *)calloc(1,sizeof(*p)) ;
   if( p==NULL ){ pp_error(err,nerr,"out of memory allocating preprocessing") ; return NULL ; }
   p->nfeature_alloc=max_feature ; p->ngroup_alloc=max_group ;
   p->group_value=(int *)calloc((size_t)max_group,sizeof(int)) ;
   p->group_count=(int *)calloc((size_t)max_group,sizeof(int)) ;
   p->center=(double *)calloc((size_t)max_group*max_feature,sizeof(double)) ;
   p->scale=(double *)calloc((size_t)max_feature,sizeof(double)) ;
   p->constant=(unsigned char *)calloc((size_t)max_feature,sizeof(unsigned char)) ;
   if( p->group_value==NULL || p->group_count==NULL || p->center==NULL ||
       p->scale==NULL || p->constant==NULL ){
     THD_preprocess_free(p) ;
     pp_error(err,nerr,"out of memory allocating preprocessing") ; return NULL ;
   }
   return p ;
}

static int pp_refit( THD_preprocessor *p, int method, const float *x,
                     int nsample, int nfeature, const int *fit_index, int nfit,
                     const int *group_id, int check_rows, char *err, int nerr )
{
   int ii,jj ;
   if( p==NULL || x==NULL || fit_index==NULL || nfit<1 || nfeature<1 ||
       nfeature>p->nfeature_alloc || THD_preprocess_name(method)==NULL )
     return pp_error(err,nerr,"preprocessing input exceeds its fitted workspace") ;
   if( check_rows && !pp_check_rows(x,nsample,nfeature,fit_index,nfit,NULL,err,nerr) )
     return 0 ;
   if( method==THD_PREP_RUN_CENTER && group_id==NULL )
     return pp_error(err,nerr,"run_center needs one group ID per sample") ;
   p->method=method ; p->nfeature=nfeature ; p->nfit=nfit ;
   p->ngroup=0 ; p->nconstant=0 ;
   if( method==THD_PREP_NONE || method==THD_PREP_SAMPLE_CENTER ) return 1 ;

   if( method==THD_PREP_FEATURE_ZSCORE ){
     p->ngroup=1 ;
     memset(p->center,0,sizeof(double)*(size_t)nfeature) ;
     memset(p->constant,0,sizeof(unsigned char)*(size_t)nfeature) ;
     for( ii=0 ; ii<nfit ; ii++ ){
       const float *row=x+(size_t)fit_index[ii]*nfeature ;
       for( jj=0 ; jj<nfeature ; jj++ ) p->center[jj]+=row[jj] ;
     }
     for( jj=0 ; jj<nfeature ; jj++ ) p->center[jj]/=(double)nfit ;
     for( jj=0 ; jj<nfeature ; jj++ ){
       double ss=0.0,tol ;
       for( ii=0 ; ii<nfit ; ii++ ){
         double dd=(double)x[(size_t)fit_index[ii]*nfeature+jj]-p->center[jj] ;
         ss+=dd*dd ;
       }
       p->scale[jj]=sqrt(ss/(double)nfit) ;
       tol=DBL_EPSILON*fmax(1.0,fabs(p->center[jj])) ;
       if( p->scale[jj]<=tol ){
         p->constant[jj]=1 ; p->nconstant++ ; p->scale[jj]=1.0 ;
       }
     }
     return 1 ;
   }

   for( ii=0 ; ii<nfit ; ii++ ){
     int value=group_id[fit_index[ii]],gg ;
     if( value<0 )
       return pp_error(err,nerr,"sample %d has no valid run/group ID",fit_index[ii]+1) ;
     for( gg=0 ; gg<p->ngroup ; gg++ ) if( p->group_value[gg]==value ) break ;
     if( gg==p->ngroup ){
       if( p->ngroup>=p->ngroup_alloc )
         return pp_error(err,nerr,"run_center has more groups than its workspace") ;
       p->group_value[p->ngroup++]=value ;
     }
   }
   memset(p->center,0,sizeof(double)*(size_t)p->ngroup*nfeature) ;
   memset(p->group_count,0,sizeof(int)*(size_t)p->ngroup) ;
   for( ii=0 ; ii<nfit ; ii++ ){
     int ss=fit_index[ii],gg=pp_group_index(p,group_id[ss]) ;
     p->group_count[gg]++ ;
     for( jj=0 ; jj<nfeature ; jj++ )
       p->center[(size_t)gg*nfeature+jj]+=x[(size_t)ss*nfeature+jj] ;
   }
   for( ii=0 ; ii<p->ngroup ; ii++ ) for( jj=0 ; jj<nfeature ; jj++ )
     p->center[(size_t)ii*nfeature+jj]/=(double)p->group_count[ii] ;
   return 1 ;
}

THD_preprocessor * THD_preprocess_fit( int method, const float *x,
                                        int nsample, int nfeature,
                                        const int *fit_index, int nfit,
                                        const int *group_id,
                                        char *err, int nerr )
{
   THD_preprocessor *p ; int ii,jj,max_group=1 ;
   if( err!=NULL && nerr>0 ) err[0]='\0' ;
   if( THD_preprocess_name(method)==NULL || x==NULL || fit_index==NULL ||
       nsample<1 || nfeature<1 || nfit<1 ){
     pp_error(err,nerr,"preprocessing needs a valid method, matrix, and sample list") ;
     return NULL ;
   }
   if( !pp_check_rows(x,nsample,nfeature,fit_index,nfit,NULL,err,nerr) ) return NULL ;
   if( method==THD_PREP_RUN_CENTER ){
     if( group_id==NULL ){
       pp_error(err,nerr,"run_center needs one group ID per sample") ; return NULL ;
     }
     max_group=0 ;
     for( ii=0 ; ii<nfit ; ii++ ){
       int value=group_id[fit_index[ii]] ;
       if( value<0 ){
         pp_error(err,nerr,"sample %d has no valid run/group ID",fit_index[ii]+1) ;
         return NULL ;
       }
       for( jj=0 ; jj<ii ; jj++ ) if( group_id[fit_index[jj]]==value ) break ;
       if( jj==ii ) max_group++ ;
     }
   }
   p=pp_alloc(nfeature,max_group,err,nerr) ;
   if( p==NULL ) return NULL ;
   if( !pp_refit(p,method,x,nsample,nfeature,fit_index,nfit,group_id,0,err,nerr) ){
     THD_preprocess_free(p) ; return NULL ;
   }
   return p ;
}

static int pp_apply( const THD_preprocessor *p, const float *x,
                     int nsample, int nfeature, const int *apply_index,
                     int napply, const int *group_id, float *out,
                     int check_rows, char *err, int nerr )
{
   int ii,jj ;
   if( p==NULL || out==NULL || p->nfeature!=nfeature )
     return pp_error(err,nerr,"invalid fitted preprocessing object or output shape") ;
   if( out==x ) return pp_error(err,nerr,"preprocessing output must not alias its input matrix") ;
   if( p->method==THD_PREP_RUN_CENTER && group_id==NULL )
     return pp_error(err,nerr,"run_center needs one group ID per sample") ;
   if( check_rows && !pp_check_rows(x,nsample,nfeature,apply_index,napply,NULL,err,nerr) )
     return 0 ;
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

int THD_preprocess_apply( const THD_preprocessor *p, const float *x,
                           int nsample, int nfeature,
                           const int *apply_index, int napply,
                           const int *group_id, float *out,
                           char *err, int nerr )
{
   if( err!=NULL && nerr>0 ) err[0]='\0' ;
   return pp_apply(p,x,nsample,nfeature,apply_index,napply,group_id,out,1,err,nerr) ;
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

void THD_mvpa_workspace_free( THD_mvpa_workspace *w )
{
   if( w==NULL ) return ;
   free(w->train_x) ; free(w->test_x) ;
   free(w->train_class) ; free(w->test_class) ; free(w->test_prediction) ;
   free(w->test_decision) ; free(w->train_count) ;
   free(w->fold_confusion) ; free(w->fold_recall) ;
   free(w->prediction) ; free(w->decision) ; free(w->fold_nsample) ;
   free(w->fold_accuracy) ; free(w->fold_balanced_accuracy) ;
   free(w->confusion) ; free(w->recall) ;
   THD_preprocess_free(w->prep) ; THD_centroid_free(w->centroid) ; free(w) ;
}

static THD_mvpa_workspace * workspace_fail( THD_mvpa_workspace *w,
                                             char *err, int nerr,
                                             const char *message )
{
   THD_mvpa_workspace_free(w) ; pp_error(err,nerr,"%s",message) ; return NULL ;
}

THD_mvpa_workspace * THD_mvpa_workspace_alloc( const THD_foldplan *plan,
                                                int max_feature, int nclass,
                                                int preprocess,
                                                const int *group_id,
                                                char *err, int nerr )
{
   THD_mvpa_workspace *w ; int ff,ii,jj,max_group=1 ;
   if( err!=NULL && nerr>0 ) err[0]='\0' ;
   if( !THD_foldplan_validate(plan,err,nerr) ) return NULL ;
   if( max_feature<1 || nclass<2 || THD_preprocess_name(preprocess)==NULL ){
     pp_error(err,nerr,"workspace needs features, at least two classes, and valid preprocessing") ;
     return NULL ;
   }
   if( preprocess==THD_PREP_RUN_CENTER && group_id==NULL ){
     pp_error(err,nerr,"run_center workspace needs one group ID per sample") ;
     return NULL ;
   }
   w=(THD_mvpa_workspace *)calloc(1,sizeof(*w)) ;
   if( w==NULL ){ pp_error(err,nerr,"out of memory allocating MVPA workspace") ; return NULL ; }
   w->plan=plan ; w->nsample=plan->nsample ; w->nfold=plan->nfold ;
   w->group_id=(preprocess==THD_PREP_RUN_CENTER)?group_id:NULL ;
   w->nclass=nclass ; w->nfeature_alloc=max_feature ; w->preprocess=preprocess ;
   for( ff=0 ; ff<plan->nfold ; ff++ ){
     if( plan->ntrain[ff]>w->max_train ) w->max_train=plan->ntrain[ff] ;
     if( plan->ntest[ff]>w->max_test ) w->max_test=plan->ntest[ff] ;
   }
   if( (size_t)w->max_train > SIZE_MAX/sizeof(float)/(size_t)max_feature ||
       (size_t)w->max_test > SIZE_MAX/sizeof(float)/(size_t)max_feature ||
       (size_t)w->max_test > SIZE_MAX/sizeof(double)/(size_t)nclass ||
       (size_t)w->nsample > SIZE_MAX/sizeof(double)/(size_t)nclass ||
       (size_t)nclass > SIZE_MAX/sizeof(int)/(size_t)nclass )
     return workspace_fail(w,err,nerr,"MVPA workspace dimensions are too large") ;

   w->train_x=(float *)malloc(sizeof(float)*(size_t)w->max_train*max_feature) ;
   w->test_x=(float *)malloc(sizeof(float)*(size_t)w->max_test*max_feature) ;
   w->train_class=(int *)malloc(sizeof(int)*(size_t)w->max_train) ;
   w->test_class=(int *)malloc(sizeof(int)*(size_t)w->max_test) ;
   w->test_prediction=(int *)malloc(sizeof(int)*(size_t)w->max_test) ;
   w->test_decision=(double *)malloc(sizeof(double)*(size_t)w->max_test*nclass) ;
   w->train_count=(int *)malloc(sizeof(int)*(size_t)nclass) ;
   w->fold_confusion=(int *)malloc(sizeof(int)*(size_t)nclass*nclass) ;
   w->fold_recall=(double *)malloc(sizeof(double)*(size_t)nclass) ;
   w->prediction=(int *)malloc(sizeof(int)*(size_t)w->nsample) ;
   w->decision=(double *)malloc(sizeof(double)*(size_t)w->nsample*nclass) ;
   w->fold_nsample=(int *)malloc(sizeof(int)*(size_t)w->nfold) ;
   w->fold_accuracy=(double *)malloc(sizeof(double)*(size_t)w->nfold) ;
   w->fold_balanced_accuracy=(double *)malloc(sizeof(double)*(size_t)w->nfold) ;
   w->confusion=(int *)malloc(sizeof(int)*(size_t)nclass*nclass) ;
   w->recall=(double *)malloc(sizeof(double)*(size_t)nclass) ;
   if( w->train_x==NULL || w->test_x==NULL || w->train_class==NULL ||
       w->test_class==NULL || w->test_prediction==NULL || w->test_decision==NULL ||
       w->train_count==NULL || w->fold_confusion==NULL ||
       w->fold_recall==NULL || w->prediction==NULL || w->decision==NULL ||
       w->fold_nsample==NULL || w->fold_accuracy==NULL ||
       w->fold_balanced_accuracy==NULL || w->confusion==NULL || w->recall==NULL )
     return workspace_fail(w,err,nerr,"out of memory allocating MVPA workspace buffers") ;
   if( preprocess==THD_PREP_RUN_CENTER ){
     max_group=0 ;
     for( ii=0 ; ii<w->nsample ; ii++ ){
       if( group_id[ii]<0 ){
         THD_mvpa_workspace_free(w) ;
         pp_error(err,nerr,"sample %d has no valid run/group ID",ii+1) ; return NULL ;
       }
       for( jj=0 ; jj<ii ; jj++ ) if( group_id[jj]==group_id[ii] ) break ;
       if( jj==ii ) max_group++ ;
     }
   }
   w->prep=pp_alloc(max_feature,max_group,err,nerr) ;
   if( w->prep==NULL ){ THD_mvpa_workspace_free(w) ; return NULL ; }
   w->centroid=centroid_alloc(nclass,max_feature,err,nerr) ;
   if( w->centroid==NULL ){ THD_mvpa_workspace_free(w) ; return NULL ; }
   return w ;
}

int THD_mvpa_run_centroid_cv( const THD_foldplan *plan,
                               int distance, const float *x,
                               int nfeature, const int *class_id,
                               const int *group_id, THD_mvpa_workspace *w,
                               char *err, int nerr )
{
   int ff,ii,jj,total=0 ; char why[512] ; THD_mvpa_class_score fold_score ;
   if( err!=NULL && nerr>0 ) err[0]='\0' ;
   if( w==NULL || plan==NULL || plan!=w->plan || x==NULL || class_id==NULL ||
       nfeature<1 || nfeature>w->nfeature_alloc )
     return pp_error(err,nerr,"CV inputs do not match the allocated MVPA workspace") ;
   if( THD_centroid_distance_name(distance)==NULL )
     return pp_error(err,nerr,"invalid centroid distance %d",distance) ;
   if( w->preprocess==THD_PREP_RUN_CENTER && group_id!=w->group_id )
     return pp_error(err,nerr,"run_center group IDs do not match the allocated workspace") ;
   for( ii=0 ; ii<w->nsample ; ii++ ){
     if( class_id[ii]<0 || class_id[ii]>=w->nclass )
       return pp_error(err,nerr,"sample %d has invalid class index %d",ii+1,class_id[ii]) ;
     for( jj=0 ; jj<nfeature ; jj++ )
       if( !isfinite(x[(size_t)ii*nfeature+jj]) )
         return pp_error(err,nerr,"sample %d, feature %d is not finite",ii+1,jj+1) ;
   }
   for( ff=0 ; ff<w->nfold ; ff++ ){
     memset(w->train_count,0,sizeof(int)*(size_t)w->nclass) ;
     for( ii=0 ; ii<plan->ntrain[ff] ; ii++ )
       w->train_count[class_id[plan->train[ff][ii]]]++ ;
     for( ii=0 ; ii<w->nclass ; ii++ ) if( w->train_count[ii]==0 )
       return pp_error(err,nerr,"fold %s has no class %d in its training set",plan->name[ff],ii+1) ;
   }

   for( ii=0 ; ii<w->nsample ; ii++ ){
     w->prediction[ii]=-1 ;
     for( jj=0 ; jj<w->nclass ; jj++ )
       w->decision[(size_t)ii*w->nclass+jj]=NAN ;
   }
   memset(w->confusion,0,sizeof(int)*(size_t)w->nclass*w->nclass) ;
   w->ntest=0 ;
   for( ff=0 ; ff<w->nfold ; ff++ ){
     int ntrain=plan->ntrain[ff],ntest=plan->ntest[ff] ;
     if( !pp_refit(w->prep,w->preprocess,x,w->nsample,nfeature,
                   plan->train[ff],ntrain,group_id,0,why,sizeof(why)) )
       return pp_error(err,nerr,"fold %s preprocessing fit: %s",plan->name[ff],why) ;
     if( !pp_apply(w->prep,x,w->nsample,nfeature,plan->train[ff],ntrain,
                   group_id,w->train_x,0,why,sizeof(why)) ||
         !pp_apply(w->prep,x,w->nsample,nfeature,plan->test[ff],ntest,
                   group_id,w->test_x,0,why,sizeof(why)) )
       return pp_error(err,nerr,"fold %s preprocessing apply: %s",plan->name[ff],why) ;
     for( ii=0 ; ii<ntrain ; ii++ ) w->train_class[ii]=class_id[plan->train[ff][ii]] ;
     for( ii=0 ; ii<ntest ; ii++ ) w->test_class[ii]=class_id[plan->test[ff][ii]] ;
     if( !centroid_refit(w->centroid,distance,w->train_x,w->train_class,
                         ntrain,nfeature,why,sizeof(why)) )
       return pp_error(err,nerr,"fold %s classifier fit: %s",plan->name[ff],why) ;
     if( !THD_centroid_predict(w->centroid,w->test_x,ntest,nfeature,
                               w->test_prediction,w->test_decision,why,sizeof(why)) )
       return pp_error(err,nerr,"fold %s prediction: %s",plan->name[ff],why) ;
     if( !THD_mvpa_score_classification(w->test_class,w->test_prediction,
                                        ntest,w->nclass,w->fold_confusion,
                                        w->fold_recall,&fold_score,why,sizeof(why)) )
       return pp_error(err,nerr,"fold %s scoring: %s",plan->name[ff],why) ;
     w->fold_nsample[ff]=ntest ; w->fold_accuracy[ff]=fold_score.accuracy ;
     w->fold_balanced_accuracy[ff]=fold_score.balanced_accuracy ;
     for( ii=0 ; ii<ntest ; ii++ ){
       int ss=plan->test[ff][ii],truth=w->test_class[ii],pred=w->test_prediction[ii] ;
       w->prediction[ss]=pred ;
       for( jj=0 ; jj<w->nclass ; jj++ )
         w->decision[(size_t)ss*w->nclass+jj]=w->test_decision[(size_t)ii*w->nclass+jj] ;
       w->confusion[(size_t)truth*w->nclass+pred]++ ; total++ ;
     }
   }
   w->ntest=total ;
   mvpa_class_score_from_confusion(w->confusion,total,w->nclass,w->recall,&w->score) ;
   return 1 ;
}

int THD_mvpa_run_centroid_permutations(
                               const THD_foldplan *plan,
                               int distance, const float *x,
                               int nfeature, const int *class_id,
                               const int *group_id, int metric,
                               const PERM_set *pset,
                               int *permuted_class, float *null_score,
                               THD_mvpa_workspace *w,
                               char *err, int nerr )
{
   int ip,ii ; char why[512] ;
   if( err!=NULL && nerr>0 ) err[0]='\0' ;
   if( plan==NULL || class_id==NULL || pset==NULL || permuted_class==NULL ||
       null_score==NULL || w==NULL || pset->nobs!=plan->nsample ||
       pset->nperm<1 || pset->perm==NULL || pset->sign==NULL )
     return pp_error(err,nerr,"permutation CV needs matching folds, relabelings, scratch, and workspace") ;
   if( metric!=THD_METRIC_ACCURACY && metric!=THD_METRIC_BALANCED_ACCURACY )
     return pp_error(err,nerr,"classification permutations need accuracy or balanced_accuracy") ;

   for( ip=0 ; ip<pset->nperm ; ip++ ){
     const int *perm=pset->perm+(size_t)ip*pset->nobs ;
     const signed char *sign=pset->sign+(size_t)ip*pset->nobs ;
     for( ii=0 ; ii<pset->nobs ; ii++ ){
       int src=perm[ii] ;
       if( src<0 || src>=pset->nobs || sign[ii]!=1 )
         return pp_error(err,nerr,"invalid label relabeling %d at sample %d",ip+1,ii+1) ;
       permuted_class[ii]=class_id[src] ;
     }
     if( !THD_mvpa_run_centroid_cv(plan,distance,x,nfeature,permuted_class,
                                    group_id,w,why,sizeof(why)) )
       return pp_error(err,nerr,"label relabeling %d: %s",ip+1,why) ;
     null_score[ip]=(float)(metric==THD_METRIC_ACCURACY ? w->score.accuracy
                                                        : w->score.balanced_accuracy) ;
     if( !isfinite(null_score[ip]) )
       return pp_error(err,nerr,"label relabeling %d produced a non-finite score",ip+1) ;
   }

   if( !THD_mvpa_run_centroid_cv(plan,distance,x,nfeature,class_id,group_id,
                                  w,why,sizeof(why)) )
     return pp_error(err,nerr,"restoring observed labels after permutation: %s",why) ;
   return 1 ;
}
