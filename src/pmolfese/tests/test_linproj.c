#include "thd_picard.h"
#include <math.h>
#include <stdio.h>
#define CHECK(c,m) do{if(!(c)){fprintf(stderr,"line %d: %s\n",__LINE__,m);return 1;}}while(0)
static int exercise(THD_linproj_workspace *w)
{
   double y[]={1,4,3,4,5,4},out[6],gamma[4];char err[256];int i;
   CHECK(THD_linproj_project(w,y,2,out,err,sizeof(err)),err);
   CHECK(fabs(out[0]+2)<1e-12 && fabs(out[2])<1e-12 && fabs(out[4]-2)<1e-12,"intercept span projection");
   for(i=0;i<3;i++)CHECK(out[2*i+1]==0,"constant response projected to numerical zero");
   CHECK(THD_linproj_solve(w,y,2,gamma,err,sizeof(err)),err);
   CHECK(fabs(gamma[0]-.6)<1e-12 && fabs(gamma[2]-1.2)<1e-12,"redundant Z minimum-norm recovery");
   CHECK(THD_linproj_project(w,out,2,out,err,sizeof(err)),err);
   CHECK(fabs(out[0]+2)<1e-12 && fabs(out[4]-2)<1e-12,"idempotent in-place projection");
   y[0]=NAN;CHECK(!THD_linproj_project(w,y,2,out,err,sizeof(err)),"nonfinite response rejected");
   CHECK(!THD_linproj_project(w,out,3,out,err,sizeof(err)),"capacity enforced");
   return 0;
}
int main(void)
{
   double z[]={1,2,1,2,1,2},y[]={1,2,4},out[3];char err[256];int fail;
   THD_linproj *p=THD_linproj_prepare(z,3,2,1e-10,err,sizeof(err));THD_linproj_workspace *w;
   CHECK(p,err);CHECK(THD_linproj_rank(p)==1,"redundant nuisance rank");
   w=THD_linproj_workspace_alloc(p,2,err,sizeof(err));
   if(!w){THD_linproj_free(p);fprintf(stderr,"%s\n",err);return 1;}
   fail=exercise(w);THD_linproj_workspace_free(w);THD_linproj_free(p);
   p=THD_linproj_prepare(NULL,3,0,1e-10,err,sizeof(err));CHECK(p,err);
   w=THD_linproj_workspace_alloc(p,1,err,sizeof(err));
   if(!w){THD_linproj_free(p);return 1;}
   if(!THD_linproj_project(w,y,1,out,err,sizeof(err)) || out[0]!=1 || out[1]!=2 || out[2]!=4 ||
      !THD_linproj_solve(w,y,1,NULL,err,sizeof(err)))fail=1;
   THD_linproj_workspace_free(w);THD_linproj_free(p);
   CHECK(!THD_linproj_prepare(z,3,2,0,err,sizeof(err)),"invalid rank tolerance rejected");
   if(!fail)puts("PASS nuisance rank, minimum-norm recovery, projection, and no implicit intercept");
   return fail;
}
