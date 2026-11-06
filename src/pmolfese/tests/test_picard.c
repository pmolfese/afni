#include "thd_picard.h"
#include <math.h>
#include <stdio.h>
#define CHECK(c,m) do{if(!(c)){fprintf(stderr,"line %d: %s\n",__LINE__,m);return 1;}}while(0)
static int exercise(THD_picard_workspace *w)
{
   double y[]={1,3,3,6,5,9},f[]={.5,1};char err[256];const THD_picard_result *r;int i;
   CHECK(!THD_picard_get_result(w),"unfitted result invalid");
   CHECK(THD_picard_fit_batch(w,y,2,f,err,sizeof(err)),err);r=THD_picard_get_result(w);CHECK(r,"fitted result");
   CHECK(fabs(r->beta[0]-1)<1e-9 && fabs(r->beta[1]-3)<1e-9,"independent per-response fractions");
   CHECK(fabs(r->gamma[0]-2)<1e-9 && fabs(r->gamma[1]-3)<1e-9,"nuisance recovered from Y-XB");
   CHECK(fabs(r->diagnostic[0].penalty-2)<1e-8 && fabs(r->diagnostic[0].lambda-2.0/3)<1e-8,"normalized penalty units");
   for(i=0;i<6;i++)CHECK(fabs(r->fitted[i]+r->residual[i]-y[i])<1e-12,"original-unit reconstruction");
   y[0]=NAN;CHECK(!THD_picard_fit_batch(w,y,2,f,err,sizeof(err)) && !THD_picard_get_result(w),"failed fit invalidates result");
   y[0]=1;f[0]=0;f[1]=1;
   CHECK(THD_picard_fit_batch(w,y,2,f,err,sizeof(err)),err);r=THD_picard_get_result(w);
   CHECK(r->beta[0]==0 && isinf(r->diagnostic[0].penalty),"zero fraction endpoint");
   CHECK(fabs(r->gamma[0]-3)<1e-12,"zero fraction leaves nuisance fitted");
   for(i=0;i<6;i++)y[i]=7;
   CHECK(THD_picard_fit_batch(w,y,2,f,err,sizeof(err)),err);r=THD_picard_get_result(w);
   CHECK(r->diagnostic[0].status==THD_FRACTION_ZERO_REFERENCE && isnan(r->diagnostic[0].achieved_fraction),"zero reference policy");
   CHECK(r->beta[0]==0 && fabs(r->gamma[0]-7)<1e-12,"nuisance-only response");
   return 0;
}
static int sink(void *context,int c,const THD_fraction_result *r,char *err,int nerr)
{
   int *count=context;(void)err;(void)nerr;
   if(c!=*count || r->q!=1 || fabs(r->coef[0]-r->coef[1])>1e-10)return 0;
   (*count)++;return 1;
}
static int reject(void *context,int c,const THD_fraction_result *r,char *err,int nerr)
{(void)context;(void)c;(void)r;(void)err;(void)nerr;return 0;}
static int fractional(THD_regression_workspace *w)
{
   double y[]={2,4},f[]={0,.3,1};char err[256];int count=0;const THD_fraction_result *r;
   CHECK(THD_regression_fraction_path(w,y,1,f,3,sink,&count,err,sizeof(err)) && count==3,"streamed wide rank-deficient path");
   r=THD_regression_fraction_result(w);
   CHECK(r->diagnostic[0].rank==1 && fabs(r->coef[0]-1)<1e-10 && fabs(r->coef[1]-1)<1e-10,"minimum-norm endpoint");
   CHECK(!THD_regression_fraction_path(w,y,1,f,3,reject,NULL,err,sizeof(err)) && !THD_regression_fraction_result(w),"callback rejection invalidates path");
   CHECK(THD_regression_fit_fraction(w,y,1,.3,err,sizeof(err)),err);r=THD_regression_fraction_result(w);
   CHECK(fabs(r->coef[0]-.3)<1e-9 && fabs(r->diagnostic[0].achieved_fraction-.3)<1e-10,"exact norm fraction");
   y[0]=NAN;CHECK(!THD_regression_fit_fraction(w,y,1,.3,err,sizeof(err)) && !THD_regression_fraction_result(w),"nonfinite fit invalidation");
   return 0;
}
static int rank_policy(void)
{
   double x[]={1,0,0,1e-12},y[]={4,7e-12};int i,fail=0;char err[256];
   for(i=0;i<2;i++){
      THD_regression_design *d=THD_regression_design_prepare(x,2,2,i?1e-14:1e-10,err,sizeof(err));
      THD_regression_workspace *w;const THD_fraction_result *r;
      CHECK(d,err);w=THD_regression_workspace_alloc(d,1,err,sizeof(err));
      if(!w){THD_regression_design_free(d);return 1;}
      if(!THD_regression_fit_fraction(w,y,1,1,err,sizeof(err)))fail=1;
      r=THD_regression_fraction_result(w);
      if(!r || r->diagnostic[0].rank!=i+1 || fabs(r->coef[0]-4)>1e-12 ||
         fabs(r->coef[1]-(i?7:0))>1e-12)fail=1;
      THD_regression_workspace_free(w);THD_regression_design_free(d);
   }
   CHECK(!fail,"relative SVD cutoff controls estimable reference span");return 0;
}
/* Orthogonal identity design: beta(f)=f*y. This gives a hand-computable
   ordered-pair loss without using an independent implementation of ridge. */
static int repeated_conditions(void)
{
   double x[] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
   double y[] = {4,1, 0,2, 0,1, 2,2}, fractions[] = {.5,0,1};
   int roles[] = {1,1,1,1}, runs[] = {0,0,1,1};
   THD_picard_trial trials[] = {{0,10,0,0},{1,11,1,0},
                               {2,12,0,1},{3,13,1,1}};
   char err[256];
   THD_picard_plan *p = THD_picard_plan_create(x,4,4,roles,NULL,runs,4,
                                              1e-10,err,sizeof(err));
   THD_picard_repeat_plan *rp;
   THD_picard_repeat_workspace *w;
   const THD_picard_repeat_result *r;
   int i;
   CHECK(p,err);
   rp=THD_picard_repeat_plan_create(p,trials,fractions,3,err,sizeof(err));
   CHECK(rp,err);
   CHECK(THD_picard_repeat_plan_info(rp)->npair==4,"ordered pair counts");
   /* Metadata and grid belong to the plan, not to the caller. */
   trials[0].condition=999; fractions[0]=.9;
   w=THD_picard_repeat_workspace_alloc(rp,2,err,sizeof(err));
   CHECK(w,err);
   CHECK(THD_picard_repeat_fit(w,y,2,err,sizeof(err)),err);
   r=THD_picard_repeat_get_result(w);
   CHECK(r->candidate[0]==1 && r->candidate[1]==2,"per-response winners");
   CHECK(fabs(r->loss[0]-5)<1e-10 && r->loss[1]<1e-20,"known pair MSE");
   CHECK(fabs(r->reference_loss[0]-10)<1e-10,"OLS reference MSE");
   CHECK(fabs(r->fold_loss[2]-8)<1e-10 &&
         fabs(r->fold_loss[10]-2)<1e-10,"fold-major loss layout");
   CHECK(fabs(r->candidate_loss[0]-6.25)<1e-9,"half-fraction MSE");
   for(i=0;i<4;i++) {
      CHECK(r->fit->beta[2*i]==0,"winning zero fraction refit");
      CHECK(fabs(r->fit->beta[2*i+1]-y[2*i+1])<1e-10,"OLS winning refit");
   }
   /* A zero response has equal losses at every candidate: choose fraction 1.
      A smaller second batch also checks packed response strides. */
   for(i=0;i<4;i++)y[i]=0;
   CHECK(THD_picard_repeat_fit(w,y,1,err,sizeof(err)),err);
   CHECK(THD_picard_repeat_get_result(w)->candidate[0]==2,"tie prefers 1");
   y[0]=NAN;
   CHECK(!THD_picard_repeat_fit(w,y,1,err,sizeof(err)) &&
         !THD_picard_repeat_get_result(w),"failed tuning invalidates result");
   THD_picard_repeat_workspace_free(w);
   THD_picard_repeat_plan_free(rp);
   for(i=0;i<4;i++)trials[i].condition=i;
   CHECK(!THD_picard_repeat_plan_create(p,trials,fractions,3,err,sizeof(err)),
         "reject no cross-run repeats");
   trials[2].trial=trials[0].trial;
   CHECK(!THD_picard_repeat_plan_create(p,trials,fractions,3,err,sizeof(err)),
         "reject duplicate trial IDs");
   THD_picard_plan_free(p);
   return 0;
}

/* Centering is supplied by the explicit nuisance intercept; scaling is then
   applied only to projected X. The known one-feature solution makes nuisance
   recovery and coefficient back-transformation independently checkable. */
static int sparse_workflow(void)
{
   double x[] = {1,0, 1,1, 1,2}, y[] = {1,3,5}, lambda[] = {.5};
   int roles[] = {0,1}, scaled;
   char err[256];
   for (scaled = 0; scaled < 2; scaled++) {
      THD_elastic_config cfg = {1,1e-12,1000,scaled};
      THD_picard_plan *p = THD_picard_plan_create_elastic(
         x,3,2,roles,NULL,NULL,3,1e-10,&cfg,err,256);
      THD_picard_workspace *w;
      const THD_picard_result *r;
      double expected = scaled ? 2-.5/sqrt(2.0/3) : 1.25;
      CHECK(p,err);
      w = THD_picard_workspace_alloc(p,1,err,256);
      CHECK(w,err);
      CHECK(THD_picard_fit_batch(w,y,1,lambda,err,256),err);
      r = THD_picard_get_result(w);
      CHECK(r && !r->diagnostic && r->elastic_diagnostic,
            "sparse diagnostic type");
      CHECK(fabs(r->beta[0]-expected) < 1e-10,"original-unit sparse beta");
      CHECK(fabs(r->gamma[0]-(3-expected)) < 1e-10,"unpenalized intercept");
      y[1] = NAN;
      CHECK(!THD_picard_fit_batch(w,y,1,lambda,err,256) &&
            !THD_picard_get_result(w),"nonfinite response invalidates result");
      y[1] = 3;
      CHECK(THD_picard_fit_batch(w,y,1,lambda,err,256),err);
      THD_picard_workspace_free(w);
      THD_picard_plan_free(p);
   }
   return 0;
}

int main(void)
{
   double matrix[]={1,0,1,1,1,2},wide[]={1,1,0,2,2,0};
   int roles[]={0,1},rows[]={0,2,4},runs[]={0,0,1},failed;char err[256];
   THD_picard_plan *p=THD_picard_plan_create(matrix,3,2,roles,rows,runs,5,1e-10,err,sizeof(err));
   THD_picard_workspace *w;THD_regression_design *d;THD_regression_workspace *rw;
   CHECK(p,err);matrix[1]=999;rows[0]=99;roles[1]=0;
   CHECK(THD_picard_plan_info(p)->row_ids[0]==0 && THD_picard_plan_info(p)->penalized_columns[0]==1,"plan owns maps");
   w=THD_picard_workspace_alloc(p,2,err,sizeof(err));
   if(!w){THD_picard_plan_free(p);return 1;}
   failed=exercise(w);THD_picard_workspace_free(w);THD_picard_plan_free(p);
   CHECK(!THD_picard_plan_create(matrix,3,2,roles,rows,runs,5,1e-10,err,sizeof(err)),"reject invalid design roles/maps");
   d=THD_regression_design_prepare(wide,2,3,1e-10,err,sizeof(err));CHECK(d,err);
   rw=THD_regression_workspace_alloc(d,1,err,sizeof(err));
   if(!rw){THD_regression_design_free(d);return 1;}
   failed+=fractional(rw);THD_regression_workspace_free(rw);THD_regression_design_free(d);
   failed+=rank_policy();
   failed+=repeated_conditions();
   failed+=sparse_workflow();
   if(!failed)puts("PASS Picard recovery, ownership, rank, endpoints, fraction paths, repeat tuning and invalidation");
   return failed;
}
