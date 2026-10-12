#include "thd_mvpa.h"

#include <float.h>
#include <math.h>

#define CHECK(COND,MSG) do{ if(!(COND)){ fprintf(stderr,"FAIL: %s\n",MSG); fail=1; goto done; } }while(0)
#define CLOSE(A,B) (fabs((double)(A)-(double)(B)) < 1.0e-6)

static void fixed_plan( THD_foldplan *p, int nsample, int *train, int ntrain,
                        int *test, int ntest, int *sample_fold )
{
   static char *name[1]={"test"},*stratum_label[1]={"All"} ;
   static int fold_stratum[1]={0},ntr[1],nte[1] ;
   static int *tr[1],*te[1] ;
   memset(p,0,sizeof(*p)) ; ntr[0]=ntrain ; nte[0]=ntest ; tr[0]=train ; te[0]=test ;
   p->kind=THD_CV_FIXED ; p->nsample=nsample ; p->nfold=1 ; p->nstratum=1 ;
   p->name=name ; p->fold_stratum=fold_stratum ; p->stratum_label=stratum_label ;
   p->ntrain=ntr ; p->ntest=nte ; p->train=tr ; p->test=te ; p->sample_test_fold=sample_fold ;
}

int main( void )
{
   float zmat[]={ 1,10,5, 3,14,5, 102,212,500, -98,-188,-500 } ;
   float smat[]={ 9,9,9, 2,4,6, 1,2,3, 8,6,4 } ;
   float rmat[]={ 1,3, 3,7, 10,20, 14,24, 4,9, 11,25 } ;
   float badmat[]={ 1,2,3, 4,NAN,6, 7,8,9, 10,11,12 } ;
   float overflow_mat[]={ 1.0f,1.0000001192092896f,FLT_MAX,0.0f } ;
   float train_out[12],test_out[12],apply_out[12] ;
   int tr2[]={0,1},te2[]={2,3},sf2[]={-1,-1,0,0} ;
   int tr4[]={0,1,2,3},te4[]={4,5},sf4[]={-1,-1,-1,-1,0,0} ;
   int group[]={10,10,20,20,10,20},bad_group[]={10,10,20,20,10,30} ;
   int class_truth[]={0,0,0,0,1,1},class_pred[]={0,0,1,0,1,0} ;
   int confusion[9],metric=-1,rank_index[8] ; double recall[3],rank_scratch[8] ;
   float cont_truth[]={1,2,3,4},cont_pred[]={1,3,2,5},cont_base[]={2.5,2.5,2.5,2.5} ;
   float tie_truth[]={1,1,2,3},tie_pred[]={1,2,2,3} ;
   float flat[]={2,2,2,2},zero_base[]={2,2,2,2} ;
   THD_mvpa_class_score cscore ; THD_mvpa_continuous_score dscore ;
   THD_foldplan plan ; THD_preprocessor *p=NULL ; char err[512] ; int fail=0,method=-1 ;

   fixed_plan(&plan,4,tr2,2,te2,2,sf2) ;
   CHECK(THD_preprocess_fold(&plan,0,THD_PREP_FEATURE_ZSCORE,zmat,3,NULL,
                             train_out,test_out,&p,err,sizeof(err)),err) ;
   CHECK(CLOSE(p->center[0],2) && CLOSE(p->center[1],12) && CLOSE(p->center[2],5),
         "feature_zscore learns training means only") ;
   CHECK(CLOSE(p->scale[0],1) && CLOSE(p->scale[1],2) && p->nconstant==1,
         "feature_zscore learns training scales and constant features") ;
   CHECK(CLOSE(train_out[0],-1) && CLOSE(train_out[1],-1) && CLOSE(train_out[2],0) &&
         CLOSE(train_out[3],1) && CLOSE(train_out[4],1) && CLOSE(train_out[5],0),
         "training rows use their fitted z-score") ;
   CHECK(CLOSE(test_out[0],100) && CLOSE(test_out[1],100) && CLOSE(test_out[2],0) &&
         CLOSE(test_out[3],-100) && CLOSE(test_out[4],-100) && CLOSE(test_out[5],0),
         "held-out values do not influence z-score parameters") ;
   THD_preprocess_free(p) ; p=NULL ;
   CHECK(THD_preprocess_fit(THD_PREP_FEATURE_ZSCORE,badmat,4,3,tr2,2,NULL,err,sizeof(err))==NULL &&
         strstr(err,"not finite")!=NULL,"non-finite training data are rejected") ;
   p=THD_preprocess_fit(THD_PREP_FEATURE_ZSCORE,overflow_mat,4,1,tr2,2,NULL,err,sizeof(err)) ; CHECK(p!=NULL,err) ;
   CHECK(!THD_preprocess_apply(p,overflow_mat,4,1,te2,2,NULL,apply_out,err,sizeof(err)) &&
         strstr(err,"outside the float range")!=NULL,"transformation overflow is rejected") ;
   THD_preprocess_free(p) ; p=NULL ;

   p=THD_preprocess_fit(THD_PREP_SAMPLE_CENTER,smat,4,3,tr2,2,NULL,err,sizeof(err)) ; CHECK(p!=NULL,err) ;
   CHECK(THD_preprocess_apply(p,smat,4,3,te2,2,NULL,apply_out,err,sizeof(err)),err) ;
   CHECK(CLOSE(apply_out[0],-1) && CLOSE(apply_out[1],0) && CLOSE(apply_out[2],1) &&
         CLOSE(apply_out[3],2) && CLOSE(apply_out[4],0) && CLOSE(apply_out[5],-2),
         "sample_center uses each sample independently") ;
   THD_preprocess_free(p) ; p=NULL ;

   p=THD_preprocess_fit(THD_PREP_NONE,smat,4,3,tr2,2,NULL,err,sizeof(err)) ; CHECK(p!=NULL,err) ;
   CHECK(THD_preprocess_apply(p,smat,4,3,te2,2,NULL,apply_out,err,sizeof(err)),err) ;
   CHECK(!memcmp(apply_out,smat+6,6*sizeof(float)),"none copies selected rows exactly") ;
   THD_preprocess_free(p) ; p=NULL ;

   fixed_plan(&plan,6,tr4,4,te4,2,sf4) ;
   CHECK(THD_preprocess_fold(&plan,0,THD_PREP_RUN_CENTER,rmat,2,group,
                             train_out,test_out,&p,err,sizeof(err)),err) ;
   CHECK(p->ngroup==2 && CLOSE(p->center[0],2) && CLOSE(p->center[1],5) &&
         CLOSE(p->center[2],12) && CLOSE(p->center[3],22),
         "run_center learns one training mean per run") ;
   CHECK(CLOSE(test_out[0],2) && CLOSE(test_out[1],4) &&
         CLOSE(test_out[2],-1) && CLOSE(test_out[3],3),
         "run_center applies training run means to held-out rows") ;
   CHECK(!THD_preprocess_apply(p,rmat,6,2,te4,2,bad_group,apply_out,err,sizeof(err)) &&
         strstr(err,"not present in training")!=NULL,
         "run_center rejects a held-out-only run") ;
   THD_preprocess_free(p) ; p=NULL ;

   CHECK(THD_preprocess_parse("FEATURE_ZSCORE",&method,err,sizeof(err)) && method==THD_PREP_FEATURE_ZSCORE,
         "preprocessing option parsing is case-insensitive") ;
   CHECK(!THD_preprocess_parse("mystery",&method,err,sizeof(err)) && strstr(err,"mystery")!=NULL,
         "unknown preprocessing option is rejected") ;

   CHECK(THD_mvpa_score_classification(class_truth,class_pred,6,2,confusion,recall,
                                        &cscore,err,sizeof(err)),err) ;
   CHECK(confusion[0]==3 && confusion[1]==1 && confusion[2]==1 && confusion[3]==1,
         "classification confusion matrix uses truth rows and prediction columns") ;
   CHECK(cscore.ncorrect==4 && cscore.nclass_present==2 && CLOSE(cscore.accuracy,2.0/3.0) &&
         CLOSE(recall[0],0.75) && CLOSE(recall[1],0.5) && CLOSE(cscore.balanced_accuracy,0.625),
         "classification accuracy and class-balanced recall are correct") ;
   CHECK(THD_mvpa_score_classification(class_truth,class_pred,6,3,confusion,recall,
                                        &cscore,err,sizeof(err)),err) ;
   CHECK(cscore.nclass_present==2 && isnan(recall[2]) && CLOSE(cscore.balanced_accuracy,0.625),
         "classes absent from truth are excluded from balanced accuracy") ;
   class_pred[5]=3 ;
   CHECK(!THD_mvpa_score_classification(class_truth,class_pred,6,3,confusion,recall,
                                         &cscore,err,sizeof(err)) && strstr(err,"outside")!=NULL,
         "invalid predicted classes are rejected") ;
   class_pred[5]=0 ;

   CHECK(THD_mvpa_score_continuous(cont_truth,cont_pred,cont_base,4,rank_scratch,rank_index,
                                    &dscore,err,sizeof(err)),err) ;
   CHECK(CLOSE(dscore.pearson,5.5/sqrt(43.75)) && CLOSE(dscore.spearman,0.8),
         "Pearson and Spearman correlations match known answers") ;
   CHECK(CLOSE(dscore.cv_r2,0.4) && CLOSE(dscore.mae,0.75) && CLOSE(dscore.rmse,sqrt(0.75)),
         "cross-validated R2 and absolute/squared errors match known answers") ;
   CHECK(THD_mvpa_score_continuous(tie_truth,tie_pred,cont_base,4,rank_scratch,rank_index,
                                    &dscore,err,sizeof(err)) && CLOSE(dscore.spearman,5.0/6.0),
         "Spearman scoring assigns average ranks to ties") ;
   CHECK(THD_mvpa_score_continuous(flat,flat,zero_base,4,rank_scratch,rank_index,
                                    &dscore,err,sizeof(err)),err) ;
   CHECK(isnan(dscore.pearson) && isnan(dscore.spearman) && isnan(dscore.cv_r2) &&
         CLOSE(dscore.mae,0) && CLOSE(dscore.rmse,0),
         "undefined continuous scores are NAN while error metrics remain defined") ;
   cont_pred[2]=NAN ;
   CHECK(!THD_mvpa_score_continuous(cont_truth,cont_pred,cont_base,4,rank_scratch,rank_index,
                                     &dscore,err,sizeof(err)) && strstr(err,"not finite")!=NULL,
         "non-finite continuous predictions are rejected") ;
   cont_pred[2]=2 ;

   CHECK(THD_mvpa_metric_parse("BALANCED_ACCURACY",&metric,err,sizeof(err)) &&
         metric==THD_METRIC_BALANCED_ACCURACY &&
         strcmp(THD_mvpa_metric_name(metric),"balanced_accuracy")==0,
         "metric option parsing is case-insensitive") ;
   CHECK(!THD_mvpa_metric_parse("mystery",&metric,err,sizeof(err)) && strstr(err,"mystery")!=NULL,
         "unknown metric option is rejected") ;

done:
   THD_preprocess_free(p) ;
   if( fail ) return 1 ;
   puts("PASS thd_mvpa preprocessing and prediction metrics") ; return 0 ;
}
