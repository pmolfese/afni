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

static void paired_three_fold_plan( THD_foldplan *p )
{
   static char *name[3]={"pair1","pair2","pair3"},*stratum_label[1]={"All"} ;
   static char group_column[]="Pair" ;
   static int fold_stratum[3]={0,0,0},ntr[3]={4,4,4},nte[3]={2,2,2} ;
   static int tr0[]={2,3,4,5},tr1[]={0,1,4,5},tr2[]={0,1,2,3} ;
   static int te0[]={0,1},te1[]={2,3},te2[]={4,5},sample_fold[]={0,0,1,1,2,2} ;
   static int *tr[3]={tr0,tr1,tr2},*te[3]={te0,te1,te2} ;
   memset(p,0,sizeof(*p)) ;
   p->kind=THD_CV_LEAVE_GROUP_OUT ; p->nsample=6 ; p->nfold=3 ; p->nstratum=1 ;
   p->group_column=group_column ; p->name=name ; p->fold_stratum=fold_stratum ;
   p->stratum_label=stratum_label ; p->ntrain=ntr ; p->ntest=nte ;
   p->train=tr ; p->test=te ; p->sample_test_fold=sample_fold ;
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
   float euc_train[]={0,0, 0,2, 4,0, 4,2},euc_test[]={1,1, 3,2, 2,1} ;
   float corr_train[]={1,2,3, 2,4,6, 3,2,1, 6,4,2} ;
   float corr_test[]={10,20,30, 9,6,3},constant_test[]={5,5,5} ;
   float constant_train[]={1,1,1, 1,1,1, 1,2,3, 2,4,6} ;
   float cos_train[]={1,0, 2,0, 0,1, 0,2},cos_test[]={0.8f,0.2f, 0.1f,0.9f} ;
   float zero_pattern[]={0,0} ;
   float cv_x[]={0,0, 10,0, 1,1, 9,1, -1,2, 11,2} ;
   int cv_class[]={0,1,0,1,0,1},cv_group[]={0,0,1,1,2,2} ;
   int train_class[]={0,0,1,1},predicted[3],distance=-1 ; double decision[6] ;
   THD_mvpa_class_score cscore ; THD_mvpa_continuous_score dscore ;
   THD_centroid_model *cm=NULL ;
   THD_mvpa_workspace *ws=NULL ;
   THD_foldplan plan,cvplan ; THD_preprocessor *p=NULL ; char err[512] ; int fail=0,method=-1 ;

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

   cm=THD_centroid_fit(THD_CENTROID_EUCLIDEAN,euc_train,train_class,4,2,2,err,sizeof(err)) ;
   CHECK(cm!=NULL,err) ;
   CHECK(cm->class_count[0]==2 && cm->class_count[1]==2 &&
         CLOSE(cm->centroid[0],0) && CLOSE(cm->centroid[1],1) &&
         CLOSE(cm->centroid[2],4) && CLOSE(cm->centroid[3],1),
         "nearest-centroid fitting computes class means") ;
   CHECK(THD_centroid_predict(cm,euc_test,3,2,predicted,decision,err,sizeof(err)),err) ;
   CHECK(predicted[0]==0 && predicted[1]==1 && predicted[2]==0 &&
         CLOSE(decision[0],-1) && CLOSE(decision[1],-9) &&
         CLOSE(decision[2],-10) && CLOSE(decision[3],-2) &&
         CLOSE(decision[4],-4) && CLOSE(decision[5],-4),
         "Euclidean centroid scores and deterministic ties match known answers") ;
   THD_centroid_free(cm) ; cm=NULL ;

   cm=THD_centroid_fit(THD_CENTROID_CORRELATION,corr_train,train_class,4,3,2,err,sizeof(err)) ;
   CHECK(cm!=NULL,err) ;
   CHECK(THD_centroid_predict(cm,corr_test,2,3,predicted,decision,err,sizeof(err)),err) ;
   CHECK(predicted[0]==0 && predicted[1]==1 && CLOSE(decision[0],1) &&
         CLOSE(decision[1],-1) && CLOSE(decision[2],-1) && CLOSE(decision[3],1),
         "correlation-to-centroid predictions match known answers") ;
   CHECK(!THD_centroid_predict(cm,constant_test,1,3,predicted,decision,err,sizeof(err)) &&
         strstr(err,"constant")!=NULL,"correlation rejects a constant prediction pattern") ;
   THD_centroid_free(cm) ; cm=NULL ;
   CHECK(THD_centroid_fit(THD_CENTROID_CORRELATION,constant_train,train_class,4,3,2,
                          err,sizeof(err))==NULL && strstr(err,"constant")!=NULL,
         "correlation rejects a constant training centroid") ;

   cm=THD_centroid_fit(THD_CENTROID_COSINE,cos_train,train_class,4,2,2,err,sizeof(err)) ;
   CHECK(cm!=NULL,err) ;
   CHECK(THD_centroid_predict(cm,cos_test,2,2,predicted,NULL,err,sizeof(err)) &&
         predicted[0]==0 && predicted[1]==1,
         "cosine centroid prediction supports label-only output") ;
   CHECK(!THD_centroid_predict(cm,zero_pattern,1,2,predicted,decision,err,sizeof(err)) &&
         strstr(err,"zero norm")!=NULL,"cosine rejects a zero-norm prediction pattern") ;
   THD_centroid_free(cm) ; cm=NULL ;

   CHECK(THD_centroid_fit(THD_CENTROID_EUCLIDEAN,euc_train,train_class,4,2,3,
                          err,sizeof(err))==NULL && strstr(err,"class 2")!=NULL,
         "centroid fitting rejects a missing training class") ;
   CHECK(THD_centroid_distance_parse("CORRELATION",&distance,err,sizeof(err)) &&
         distance==THD_CENTROID_CORRELATION,
         "centroid distance parsing is case-insensitive") ;
   CHECK(!THD_centroid_distance_parse("mystery",&distance,err,sizeof(err)) &&
         strstr(err,"mystery")!=NULL,"unknown centroid distance is rejected") ;

   paired_three_fold_plan(&cvplan) ;
   ws=THD_mvpa_workspace_alloc(&cvplan,2,2,THD_PREP_FEATURE_ZSCORE,NULL,err,sizeof(err)) ;
   CHECK(ws!=NULL,err) ;
   {
     float *train_buffer=ws->train_x ; double *decision_buffer=ws->decision ;
     CHECK(THD_mvpa_run_centroid_cv(&cvplan,THD_CENTROID_EUCLIDEAN,cv_x,2,
                                    cv_class,NULL,ws,err,sizeof(err)),err) ;
     CHECK(ws->ntest==6 && ws->score.nsample==6 && ws->score.ncorrect==6 &&
           CLOSE(ws->score.accuracy,1) && CLOSE(ws->score.balanced_accuracy,1),
           "fold runner computes perfect out-of-fold scores") ;
     CHECK(!memcmp(ws->prediction,cv_class,6*sizeof(int)) &&
           ws->confusion[0]==3 && ws->confusion[1]==0 &&
           ws->confusion[2]==0 && ws->confusion[3]==3 &&
           CLOSE(ws->recall[0],1) && CLOSE(ws->recall[1],1),
           "fold runner restores predictions and diagnostics to original row order") ;
     for( method=0 ; method<3 ; method++ )
       CHECK(ws->fold_nsample[method]==2 && CLOSE(ws->fold_accuracy[method],1) &&
             CLOSE(ws->fold_balanced_accuracy[method],1),
             "fold runner records per-fold sample counts and scores") ;
     for( method=0 ; method<12 ; method++ )
       CHECK(isfinite(ws->decision[method]),"tested rows receive finite per-class decisions") ;
     CHECK(THD_mvpa_run_centroid_cv(&cvplan,THD_CENTROID_EUCLIDEAN,cv_x,2,
                                    cv_class,NULL,ws,err,sizeof(err)),err) ;
     CHECK(ws->train_x==train_buffer && ws->decision==decision_buffer,
           "repeated CV runs reuse workspace buffers") ;
     {
       int perm_index[]={0,1,2,3,4,5, 1,0,2,3,4,5} ;
       signed char perm_sign[]={1,1,1,1,1,1, 1,1,1,1,1,1} ;
       int perm_class[6] ; float null_score[2] ; PERM_set pset ;
       memset(&pset,0,sizeof(pset)) ; pset.nperm=2 ; pset.nobs=6 ;
       pset.perm=perm_index ; pset.sign=perm_sign ;
       CHECK(THD_mvpa_run_centroid_permutations(
                 &cvplan,THD_CENTROID_EUCLIDEAN,cv_x,2,cv_class,NULL,
                 THD_METRIC_BALANCED_ACCURACY,&pset,perm_class,null_score,
                 ws,err,sizeof(err)),err) ;
       CHECK(CLOSE(null_score[0],1) && isfinite(null_score[1]),
             "permutation runner returns the identity and relabeled CV scores") ;
       CHECK(CLOSE(ws->score.balanced_accuracy,1) &&
             !memcmp(ws->prediction,cv_class,6*sizeof(int)),
             "permutation runner restores observed predictions and scores") ;
       CHECK(!THD_mvpa_run_centroid_permutations(
                 &cvplan,THD_CENTROID_EUCLIDEAN,cv_x,2,cv_class,NULL,
                 THD_METRIC_PEARSON,&pset,perm_class,null_score,
                 ws,err,sizeof(err)) && strstr(err,"accuracy")!=NULL,
             "classification permutation runner rejects continuous metrics") ;
     }
   }
   THD_mvpa_workspace_free(ws) ; ws=NULL ;

   fixed_plan(&plan,6,tr4,4,te4,2,sf4) ;
   ws=THD_mvpa_workspace_alloc(&plan,2,2,THD_PREP_NONE,NULL,err,sizeof(err)) ;
   CHECK(ws!=NULL,err) ;
   CHECK(THD_mvpa_run_centroid_cv(&plan,THD_CENTROID_EUCLIDEAN,cv_x,2,
                                  cv_class,NULL,ws,err,sizeof(err)),err) ;
   CHECK(ws->ntest==2 && ws->score.nsample==2 && ws->prediction[4]==0 &&
         ws->prediction[5]==1 && ws->prediction[0]==-1 && isnan(ws->decision[0]),
         "fixed CV marks training-only rows and scores only held-out rows") ;
   THD_mvpa_workspace_free(ws) ; ws=NULL ;

   ws=THD_mvpa_workspace_alloc(&cvplan,2,2,THD_PREP_RUN_CENTER,cv_group,err,sizeof(err)) ;
   CHECK(ws!=NULL,err) ;
   CHECK(!THD_mvpa_run_centroid_cv(&cvplan,THD_CENTROID_EUCLIDEAN,cv_x,2,
                                   cv_class,cv_group,ws,err,sizeof(err)) &&
         strstr(err,"fold pair1")!=NULL && strstr(err,"not present in training")!=NULL,
         "fold runner reports run_center leakage restrictions with fold context") ;
   THD_mvpa_workspace_free(ws) ; ws=NULL ;

done:
   THD_mvpa_workspace_free(ws) ;
   THD_centroid_free(cm) ;
   THD_preprocess_free(p) ;
   if( fail ) return 1 ;
   puts("PASS thd_mvpa preprocessing, metrics, classifiers, CV, and permutation runner") ; return 0 ;
}
