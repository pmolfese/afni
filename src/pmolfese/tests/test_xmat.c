/* Unit tests for thd_xmat.c: a design matrix written here, in the form 3dDeconvolve writes, with every value
   known, read back and asked the questions analysis programs ask of it.  The item expansion against a real
   3dDeconvolve bucket is checked end to end in tests/scripts/check_3dSRM.py. */

#include "mrilib.h"
#include "thd_xmat.h"

#define CHECK(COND,MSG) do{ if(!(COND)){ fprintf(stderr,"FAIL: %s\n",MSG); return 1; } }while(0)

/* 3 runs (starts 0, 6, 12; 18 time points), time points 4, 5 and 13 censored, so 15 rows.  Columns:
   2 baseline, face#0 (run 1 only), face#1 (run 2 only), face#2 (runs 2 and 3), house#0 (everywhere). */
static const int good[15] = { 0,1,2,3,6,7,8,9,10,11,12,14,15,16,17 } ;

static double reg( int c , int t )                /* the regressors on the full timeline */
{
   switch( c ){
     case 0: return 1.0 ;
     case 1: return (double)t ;
     case 2: return (t==1) ? 1.0 : (t==2) ? 0.6 : 0.0 ;            /* face#0: run 1, peak at t=1 */
     case 3: return (t==8) ? 2.0 : (t==9) ? 0.8 : 0.0 ;            /* face#1: run 2, peak at t=8 */
     case 4: return (t==11) ? 1.0 : (t==12) ? 0.9 : 0.0 ;          /* face#2: t=11 (run 2), t=12 (run 3) */
     case 5: return (t==2 || t==15) ? 4.0 : (t==9) ? 1.0 : 0.0 ;   /* house: in runs 1, 2, 3 */
   }
   return 0.0 ;
}

static int write_xmat( const char *fn )
{
   FILE *fp=fopen(fn,"w") ; int r , c ;
   if( fp==NULL ) return 1 ;
   fprintf(fp,"# <matrix\n#  ni_type = \"6*double\"\n#  ni_dimen = \"15\"\n"
              "#  ColumnLabels = \"Run#1Pol#0 ; Run#1Pol#1 ; face#0 ; face#1 ; face#2 ; house#0\"\n"
              "#  ColumnGroups = \"2@-1,3@1,2\"\n#  RowTR = \"2.5\"\n"
              "#  GoodList = \"0..3,6..12,14..17\"\n#  NRowFull = \"18\"\n#  RunStart = \"0,6,12\"\n"
              "#  Nstim = \"2\"\n#  StimBots = \"2,5\"\n#  StimTops = \"4,5\"\n#  StimLabels = \"face ; house\"\n"
              "#  BasisOption_000001 = \"-stim_times_IM\"\n#  BasisFormula_000001 = \"GAM\"\n"
              "#  BasisOption_000002 = \"-stim_times\"\n#  BasisFormula_000002 = \"GAM\"\n"
              "#  CommandLine = \"3dDeconvolve -input x+orig -polort 1\"\n# >\n") ;
   for( r=0 ; r<15 ; r++ ){
     for( c=0 ; c<6 ; c++ ) fprintf(fp,"%s%g",c?" ":"",reg(c,good[r])) ;
     fprintf(fp,"\n") ;
   }
   fprintf(fp,"# </matrix>\n") ;
   fclose(fp) ;
   return 0 ;
}

int main( int argc , char **argv )
{
   char fn[512] , err[512] , *dir=getenv("TMPDIR") ; THD_xmat *xm ; int ii , r , *lab , n ;
   unsigned char *cens , *want ; THD_xmat_item *it ;
   snprintf(fn,sizeof(fn),"%s/thd_xmat_unit_%d.xmat.1D",dir?dir:".",(int)getpid()) ;
   CHECK(write_xmat(fn)==0,"write the test matrix") ;

   xm=THD_xmat_read(fn,err,sizeof(err)) ;
   CHECK(xm!=NULL,err) ;
   CHECK(xm->nrow==15 && xm->ncol==6 && xm->nfull==18,"rows, columns, NRowFull") ;
   for( r=0 ; r<15 ; r++ ) CHECK(xm->good[r]==good[r],"GoodList ranges (0..3,6..12,...) decoded") ;
   CHECK(xm->nrun==3 && xm->runstart[1]==6 && xm->runstart[2]==12,"RunStart") ;
   CHECK(fabs(xm->tr-2.5f)<1e-6,"RowTR") ;
   CHECK(!strcmp(xm->collab[0],"Run#1Pol#0") && !strcmp(xm->collab[4],"face#2"),"ColumnLabels, trimmed") ;
   CHECK(xm->nstim==2 && !strcmp(xm->stim[0].label,"face") && xm->stim[0].bot==2 && xm->stim[0].top==4,"stimulus face") ;
   CHECK(!strcmp(xm->stim[1].label,"house") && xm->stim[1].bot==5 && xm->stim[1].top==5,"stimulus house") ;
   CHECK(!strcmp(xm->stim[0].option,"-stim_times_IM") && !strcmp(xm->stim[1].formula,"GAM"),"basis option and formula") ;
   CHECK(fabs(xm->mat[(size_t)3*15+6]-2.0f)<1e-6,"the matrix: face#1 at row 6 (time point 8)") ;
   /* every attribute as key/value */
   CHECK(THD_xmat_attr(xm,"commandline")!=NULL && !strcmp(THD_xmat_attr(xm,"CommandLine"),"3dDeconvolve -input x+orig -polort 1"),
         "an attribute by name, case-insensitive") ;
   CHECK(THD_xmat_attr(xm,"ColumnGroups")!=NULL && THD_xmat_attr(xm,"NoSuchThing")==NULL,"any attribute; NULL if absent") ;
   /* runs, censoring */
   CHECK(THD_xmat_run_of(xm,0)==0 && THD_xmat_run_of(xm,5)==0 && THD_xmat_run_of(xm,6)==1 && THD_xmat_run_of(xm,17)==2
         && THD_xmat_run_of(xm,18)==-1,"the run of a time point") ;
   CHECK(THD_xmat_run_length(xm,0)==6 && THD_xmat_run_length(xm,2)==6,"run lengths") ;
   cens=THD_xmat_censored(xm) ;
   CHECK(THD_xmat_ncensored(xm)==3 && cens[4] && cens[5] && cens[13] && !cens[3] && !cens[14],"the censored time points") ;
   free(cens) ;
   /* which stimuli */
   want=THD_xmat_stims_wanted(xm,"house,nothing",&n) ;
   CHECK(n==1 && !want[0] && want[1],"a list of stimuli") ;
   free(want) ;
   /* each time point's condition: strongest regressor, scaled to its peak, at half of it */
   lab=THD_xmat_condition_labels(xm,0.5,NULL) ;
   CHECK(lab[1]==0 && lab[2]==1,"t=1 face (its peak); t=2 house (1.0 of its peak beats face's 0.3)") ;
   CHECK(lab[8]==0 && lab[9]==-1,"t=8 face; t=9: face 0.4 and house 0.25 of their peaks, below half") ;
   CHECK(lab[4]==-1 && lab[13]==-1,"censored time points are unlabelled") ;
   CHECK(lab[15]==1 && lab[0]==-1,"t=15 house; nothing at t=0") ;
   free(lab) ;
   want=THD_xmat_stims_wanted(xm,"face",NULL) ; lab=THD_xmat_condition_labels(xm,0.5,want) ;
   CHECK(lab[1]==0 && lab[2]==-1 && lab[15]==-1,"without house: t=2 has face at 0.3 of its peak (2.0), below half") ;
   free(lab) ; free(want) ;
   /* items: one per stimulus column, with the run its regressor lies in */
   it=THD_xmat_items(xm,NULL,&n) ;
   CHECK(n==4,"four items: face#0, face#1, face#2, house") ;
   CHECK(!strcmp(it[0].key,"face#0") && !strcmp(it[0].coef_label,"face#0_Coef") && !strcmp(it[0].condition,"face"),"a trial item") ;
   CHECK(!strcmp(it[3].key,"house") && !strcmp(it[3].coef_label,"house#0_Coef"),"a one-column stimulus is named by its label") ;
   CHECK(it[0].run==0 && it[1].run==1 && it[2].run==-1 && it[3].run==-1,"each item's run (-1 when it spans runs)") ;
   THD_xmat_items_free(it,n) ;
   /* Direct library callers cannot bypass the range checks in the reader. */
   xm->stim[1].bot=2 ; xm->stim[1].top=5 ;
   it=THD_xmat_items(xm,NULL,&n) ;
   CHECK(it==NULL && n==0,"overlapping manually constructed stimulus ranges are rejected") ;
   THD_xmat_free(xm) ;

   /* summaries as table columns */
   { char *tok[]={ "Subj","Group","Xmat","InputFile" , "s1","ctl",fn,"a+orig" , "s2","pat",fn,"b+orig" } ;
     THD_datatable *dt=THD_parse_datatable(tok,12,"a test") , *dx ;
     float *v ;
     CHECK(dt!=NULL,"a small table") ;
     dx=THD_datatable_add_xmat_columns(dt,"Xmat") ;
     CHECK(dx!=NULL && dx->ncol==10 && dx->icol_input==9,"six columns added, InputFile still last") ;
     v=THD_datatable_values(dx,"Xmat.ncensored") ; CHECK(v!=NULL && v[0]==3.0f && v[1]==3.0f,"Xmat.ncensored, numeric") ;
     v=THD_datatable_values(dx,"Xmat.frac_censored") ; CHECK(v!=NULL && fabs(v[0]-3.0/18.0)<1e-5,"Xmat.frac_censored") ;
     v=THD_datatable_values(dx,"Xmat.TR") ; CHECK(v!=NULL && fabs(v[0]-2.5f)<1e-6,"Xmat.TR") ;
     v=THD_datatable_values(dx,"Xmat.nrun") ; CHECK(v!=NULL && v[1]==3.0f,"Xmat.nrun") ;
     v=THD_datatable_values(dx,"Xmat.ntime") ; CHECK(v!=NULL && v[0]==18.0f,"Xmat.ntime") ;
     v=THD_datatable_values(dx,"Xmat.nstim") ; CHECK(v!=NULL && v[0]==2.0f,"Xmat.nstim") ;
     CHECK(!strcmp(DT_CELL(dx,1,0),"s2") && !strcmp(dx->fname[1],"b+orig"),"the original columns kept") ;
     CHECK(THD_datatable_add_xmat_columns(dx,"Xmat")==dx,"adding them twice changes nothing") ;
   }

   /* failures say why */
   { FILE *fp ; char bad[600] ;
     CHECK(THD_xmat_read("no_such_file.xmat.1D",err,sizeof(err))==NULL && strstr(err,"cannot read")!=NULL,"a missing file") ;
     snprintf(bad,sizeof(bad),"%s.bad",fn) ;
     fp=fopen(bad,"w") ;
     fprintf(fp,"# <matrix\n#  ni_type = \"1*double\"\n#  ni_dimen = \"3\"\n#  GoodList = \"0..4\"\n#  NRowFull = \"6\"\n# >\n1\n1\n1\n# </matrix>\n") ;
     fclose(fp) ;
     CHECK(THD_xmat_read(bad,err,sizeof(err))==NULL && strstr(err,"GoodList")!=NULL,"a GoodList that does not match the rows") ;
     fp=fopen(bad,"w") ;
     fprintf(fp,"# <matrix\n# ni_type = \"2*double\"\n# ni_dimen = \"2\"\n"
                "# StimLabels = \"A ; B\"\n# StimBots = \"0,0\"\n# StimTops = \"1,1\"\n"
                "# >\n1 0\n0 1\n# </matrix>\n") ;
     fclose(fp) ;
     CHECK(THD_xmat_read(bad,err,sizeof(err))==NULL && strstr(err,"ranges overlap")!=NULL,
           "overlapping file ranges fail before item expansion can overrun ncol storage") ;
     remove(bad) ; }
   remove(fn) ;
   (void)argc ; (void)argv ; (void)ii ;
   printf("thd_xmat tests passed\n") ;
   return 0 ;
}
