#include "mrilib.h"
#include "thd_patterns.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int parser_test( void )
{
   const char *spec[5]={"2","SPHERE(2)","rect(1,2,3)","RHDD(2)","TOHD(2)"} ;
   const char *bad[]={"SPHERE(0)","SPHERE(-1)","SPHERE(nan)","SPHERE(2)junk",
                      "RECT(1,2)","RECT(1,0,2)","RECT(1,2,inf)","2junk","nan"} ;
   MCW_cluster *cc[5]={NULL,NULL,NULL,NULL,NULL} ;
   char err[256] ;
   int ii,fail=0 ;

   for( ii=0 ; ii<5 ; ii++ ){
     cc[ii]=THD_searchlight_parse(spec[ii],1.0f,1.0f,1.0f,err,sizeof(err)) ;
     if( cc[ii]==NULL || cc[ii]->num_pt<1 ){
       fprintf(stderr,"FAIL parser %s: %s\n",spec[ii],err) ; fail=1 ;
     }
   }
   if( cc[0]!=NULL && cc[1]!=NULL && cc[0]->num_pt!=cc[1]->num_pt ){
     fprintf(stderr,"FAIL bare radius differs from SPHERE radius\n") ; fail=1 ;
   }
   for( ii=0 ; ii<5 ; ii++ ) if( cc[ii]!=NULL ) KILL_CLUSTER(cc[ii]) ;

   cc[0]=THD_searchlight_parse("SPHERE(0)",1.0f,1.0f,1.0f,err,sizeof(err)) ;
   if( cc[0]!=NULL || strstr(err,"finite")==NULL || strstr(err,"> 0")==NULL ){
     fprintf(stderr,"FAIL invalid-radius diagnostic: %s\n",err) ; fail=1 ;
   }
   for( ii=0 ; ii<(int)(sizeof(bad)/sizeof(bad[0])) ; ii++ ){
     cc[0]=THD_searchlight_parse(bad[ii],1.0f,1.0f,1.0f,err,sizeof(err)) ;
     if( cc[0]!=NULL ){
       fprintf(stderr,"FAIL malformed parser accepted %s\n",bad[ii]) ;
       KILL_CLUSTER(cc[0]) ; fail=1 ;
     }
   }
   cc[0]=THD_searchlight_parse("BANANA(2)",1.0f,1.0f,1.0f,err,sizeof(err)) ;
   if( cc[0]!=NULL || strstr(err,"unknown -searchlight neighborhood")==NULL ){
     fprintf(stderr,"FAIL unknown-shape diagnostic: %s\n",err) ; fail=1 ;
   }
   return fail ;
}

static int paint_test( void )
{
   THD_roilist rl ;
   int v0[2]={0,2},center[2]={-1,3} ;
   intvec vox[2] ;
   float brick[4]={-9.0f,-9.0f,-9.0f,-9.0f} ;
   float values[2]={1.5f,2.5f} ;
   float want[4]={1.5f,-9.0f,1.5f,2.5f} ;
   int ii ;

   memset(&rl,0,sizeof(rl)) ; memset(vox,0,sizeof(vox)) ;
   vox[0].nar=2 ; vox[0].ar=v0 ;
   rl.nroi=2 ; rl.vox=vox ; rl.center=center ;
   THD_roilist_paint(brick,&rl,values) ;
   for( ii=0 ; ii<4 ; ii++ ) if( brick[ii]!=want[ii] ){
     fprintf(stderr,"FAIL paint[%d]: got %.9g, wanted %.9g\n",
             ii,brick[ii],want[ii]) ;
     return 1 ;
   }
   return 0 ;
}

static THD_3dim_dataset * sample_pattern_dset( int second, double **double_brick )
{
   THD_3dim_dataset *ds=EDIT_empty_copy(NULL) ;
   THD_ivec3 nxyz ; THD_fvec3 org,del ;

   LOAD_IVEC3(nxyz,4,1,1) ; LOAD_FVEC3(org,0.0f,0.0f,0.0f) ;
   LOAD_FVEC3(del,2.0f,2.0f,2.0f) ;
   EDIT_dset_items(ds,ADN_nxyz,nxyz,ADN_xyzorg,org,ADN_xyzdel,del,
                   ADN_nvals,2,ADN_type,HEAD_FUNC_TYPE,
                   ADN_func_type,FUNC_BUCK_TYPE,ADN_none) ;
   ds->daxes->xxorient=ORI_R2L_TYPE ; ds->daxes->yyorient=ORI_A2P_TYPE ;
   ds->daxes->zzorient=ORI_I2S_TYPE ;
   if( !second ){
     float *a=(float *)malloc(4*sizeof(float)) ;
     short *b=(short *)malloc(4*sizeof(short)) ;
     int ii ;
     for( ii=0 ; ii<4 ; ii++ ){ a[ii]=(float)(ii+1) ; b[ii]=(short)(10*(ii+1)) ; }
     EDIT_substitute_brick(ds,0,MRI_float,a) ;
     EDIT_substitute_brick(ds,1,MRI_short,b) ; EDIT_BRICK_FACTOR(ds,1,0.5f) ;
   } else {
     int *a=(int *)malloc(4*sizeof(int)) ;
     double *b=(double *)malloc(4*sizeof(double)) ;
     int ii ;
     for( ii=0 ; ii<4 ; ii++ ){ a[ii]=100*(ii+1) ; b[ii]=(double)ii+1.5 ; }
     EDIT_substitute_brick(ds,0,MRI_int,a) ; EDIT_substitute_brick(ds,1,MRI_double,b) ;
     if( double_brick!=NULL ) *double_brick=b ;
   }
   return ds ;
}

static int sample_extraction_test( void )
{
   THD_3dim_dataset *dset[2]={NULL,NULL} ; THD_sampleset ss ;
   char *dname[2]={"first","second"},*label[4]={"s1","s2","s3","s4"} ;
   int sample_dset[4]={0,1,0,1},sample_brick[4]={1,0,0,1} ;
   int varray[2]={3,1},selected[3]={3,0,3},bad_sample[1]={4},ii,fail=0 ;
   intvec vox ; double *double_brick=NULL ; char err[256] ;
   float out[8],legacy[4] ;
   float want_all[8]={20,10, 400,200, 4,2, 4.5f,2.5f} ;
   float want_selected[6]={4.5f,2.5f, 20,10, 4.5f,2.5f} ;

   memset(&ss,0,sizeof(ss)) ; memset(&vox,0,sizeof(vox)) ;
   dset[0]=sample_pattern_dset(0,NULL) ; dset[1]=sample_pattern_dset(1,&double_brick) ;
   if( dset[0]==NULL || dset[1]==NULL || double_brick==NULL ){
     fprintf(stderr,"FAIL construct sample extraction datasets\n") ; fail=1 ; goto done ;
   }
   ss.nsample=4 ; ss.ndset=2 ; ss.nvox=4 ; ss.dset=dset ; ss.dset_name=dname ;
   ss.sample_dset=sample_dset ; ss.sample_brick=sample_brick ; ss.sample_label=label ;
   vox.nar=2 ; vox.ar=varray ;

   if( !THD_roi_samples(&ss,NULL,4,&vox,out,err,sizeof(err)) ){
     fprintf(stderr,"FAIL all-sample extraction: %s\n",err) ; fail=1 ; goto done ;
   }
   for( ii=0 ; ii<8 ; ii++ ) if( out[ii]!=want_all[ii] ){
     fprintf(stderr,"FAIL all-sample extraction[%d]: got %.9g, wanted %.9g\n",
             ii,out[ii],want_all[ii]) ; fail=1 ; goto done ;
   }
   if( !THD_roi_samples(&ss,selected,3,&vox,out,err,sizeof(err)) ){
     fprintf(stderr,"FAIL selected-sample extraction: %s\n",err) ; fail=1 ; goto done ;
   }
   for( ii=0 ; ii<6 ; ii++ ) if( out[ii]!=want_selected[ii] ){
     fprintf(stderr,"FAIL selected-sample extraction[%d]: got %.9g, wanted %.9g\n",
             ii,out[ii],want_selected[ii]) ; fail=1 ; goto done ;
   }

   THD_roi_pattern(dset[0],&vox,legacy) ;
   if( legacy[0]!=4.0f || legacy[1]!=2.0f || legacy[2]!=20.0f || legacy[3]!=10.0f ){
     fprintf(stderr,"FAIL legacy all-brick pattern extraction changed\n") ; fail=1 ; goto done ;
   }
   if( THD_roi_samples(&ss,NULL,3,&vox,out,err,sizeof(err)) || strstr(err,"all 4")==NULL ){
     fprintf(stderr,"FAIL all-sample count validation: %s\n",err) ; fail=1 ; goto done ;
   }
   if( THD_roi_samples(&ss,bad_sample,1,&vox,out,err,sizeof(err)) ||
       strstr(err,"canonical index 4")==NULL ){
     fprintf(stderr,"FAIL sample bounds validation: %s\n",err) ; fail=1 ; goto done ;
   }
   varray[1]=4 ;
   if( THD_roi_samples(&ss,selected,3,&vox,out,err,sizeof(err)) ||
       strstr(err,"voxel 2")==NULL || strstr(err,"index 4")==NULL ){
     fprintf(stderr,"FAIL voxel bounds validation: %s\n",err) ; fail=1 ; goto done ;
   }
   varray[1]=1 ; double_brick[1]=NAN ;
   if( THD_roi_samples(&ss,selected,1,&vox,out,err,sizeof(err)) ||
       strstr(err,"s4")==NULL || strstr(err,"not finite")==NULL ){
     fprintf(stderr,"FAIL non-finite sample validation: %s\n",err) ; fail=1 ; goto done ;
   }

done:
   if( dset[0]!=NULL ) DSET_delete(dset[0]) ;
   if( dset[1]!=NULL ) DSET_delete(dset[1]) ;
   return fail ;
}

static THD_3dim_dataset * surface_domain_dset( int nrow, const int *nodes,
                                                int nnodes )
{
   THD_3dim_dataset *ds=EDIT_empty_copy(NULL) ;
   THD_ivec3 nxyz ; THD_fvec3 org,del ; int ii ;
   LOAD_IVEC3(nxyz,nrow,1,1) ; LOAD_FVEC3(org,0.0f,0.0f,0.0f) ;
   LOAD_FVEC3(del,1.0f,1.0f,1.0f) ;
   EDIT_dset_items(ds,ADN_nxyz,nxyz,ADN_xyzorg,org,ADN_xyzdel,del,
                   ADN_nvals,1,ADN_none) ;
   ds->dblk->diskptr->storage_mode=STORAGE_BY_NI_SURF_DSET ;
   if( nodes!=NULL ){
     ds->dblk->node_list=(int *)RwcMalloc(sizeof(int)*(size_t)nrow) ;
     for( ii=0 ; ii<nrow ; ii++ ) ds->dblk->node_list[ii]=nodes[ii] ;
     ds->dblk->nnodes=nnodes ;
   }
   return ds ;
}

static int surface_domain_test( void )
{
   int identity[4]={0,1,2,3},sparse[4]={0,2,5,9},wrong[4]={0,2,5,8} ;
   THD_3dim_dataset *dense=NULL,*dense_sparse=NULL,*aa=NULL,*bb=NULL,*bad=NULL,*dst=NULL ;
   char err[256] ; int fail=0 ;

   dense=surface_domain_dset(4,NULL,0) ;
   dense_sparse=surface_domain_dset(4,identity,4) ;
   aa=surface_domain_dset(4,sparse,4) ; bb=surface_domain_dset(4,sparse,4) ;
   bad=surface_domain_dset(4,wrong,4) ; dst=surface_domain_dset(4,NULL,0) ;
   if( dense==NULL || dense_sparse==NULL || aa==NULL || bb==NULL || bad==NULL || dst==NULL ){
     fprintf(stderr,"FAIL construct surface-domain datasets\n") ; fail=1 ; goto done ;
   }
   if( !THD_surface_domain_match(dense,dense_sparse,err,sizeof(err)) ){
     fprintf(stderr,"FAIL dense/sparse identity domain comparison: %s\n",err) ; fail=1 ;
   }
   if( !THD_surface_domain_match(aa,bb,err,sizeof(err)) ){
     fprintf(stderr,"FAIL equal sparse domain comparison: %s\n",err) ; fail=1 ;
   }
   if( THD_surface_domain_match(aa,bad,err,sizeof(err)) ||
       strstr(err,"differ at row 4")==NULL ){
     fprintf(stderr,"FAIL unequal sparse domain diagnostic: %s\n",err) ; fail=1 ;
   }
   bad->dblk->nnodes=3 ;
   if( THD_surface_domain_match(aa,bad,err,sizeof(err)) ||
       strstr(err,"malformed")==NULL ){
     fprintf(stderr,"FAIL malformed sparse domain diagnostic: %s\n",err) ; fail=1 ;
   }
   if( !THD_copy_surface_node_domain(dst,aa,err,sizeof(err)) ){
     fprintf(stderr,"FAIL sparse surface domain copy: %s\n",err) ; fail=1 ;
   } else if( dst->dblk->nnodes!=4 || dst->dblk->node_list==NULL ||
              memcmp(dst->dblk->node_list,sparse,sizeof(sparse))!=0 ){
     fprintf(stderr,"FAIL copied sparse surface domain contents\n") ; fail=1 ;
   } else {
     aa->dblk->node_list[1]=7 ;
     if( dst->dblk->node_list[1]!=2 ){
       fprintf(stderr,"FAIL surface domain copy aliases source storage\n") ; fail=1 ;
     }
     aa->dblk->node_list[1]=2 ;
   }
   if( !THD_copy_surface_node_domain(dst,dense,err,sizeof(err)) ||
       dst->dblk->node_list!=NULL || dst->dblk->nnodes!=0 ){
     fprintf(stderr,"FAIL dense surface domain clears sparse destination: %s\n",err) ; fail=1 ;
   }
   {
     int unordered[4]={9,2,5,0} ;
     memcpy(aa->dblk->node_list,unordered,sizeof(unordered)) ;
     memcpy(bb->dblk->node_list,unordered,sizeof(unordered)) ;
     if( !THD_surface_domain_match(aa,bb,err,sizeof(err)) ||
         memcmp(aa->dblk->node_list,unordered,sizeof(unordered))!=0 ||
         !THD_copy_surface_node_domain(dst,aa,err,sizeof(err)) ||
         memcmp(dst->dblk->node_list,unordered,sizeof(unordered))!=0 ){
       fprintf(stderr,"FAIL valid unsorted node rows must be preserved: %s\n",err) ; fail=1 ;
     }
     aa->dblk->node_list[3]=bb->dblk->node_list[3]=2 ;
     if( THD_surface_domain_match(aa,bb,err,sizeof(err)) || strstr(err,"repeats node 2")==NULL ){
       fprintf(stderr,"FAIL matching duplicate domains must be rejected: %s\n",err) ; fail=1 ;
     }
     if( THD_copy_surface_node_domain(dst,aa,err,sizeof(err)) ||
         strstr(err,"repeats node 2")==NULL ||
         memcmp(dst->dblk->node_list,unordered,sizeof(unordered))!=0 ){
       fprintf(stderr,"FAIL duplicate copy must fail without modifying destination: %s\n",err) ; fail=1 ;
     }
   }

done:
   if( dense!=NULL ) DSET_delete(dense) ;
   if( dense_sparse!=NULL ) DSET_delete(dense_sparse) ;
   if( aa!=NULL ) DSET_delete(aa) ; if( bb!=NULL ) DSET_delete(bb) ;
   if( bad!=NULL ) DSET_delete(bad) ; if( dst!=NULL ) DSET_delete(dst) ;
   return fail ;
}

int main( void )
{
   int fail=parser_test()+paint_test()+sample_extraction_test()+surface_domain_test() ;
   if( fail ) return 1 ;
   puts("PASS thd_patterns geometry, domains, painting, and selected-sample extraction") ;
   return 0 ;
}
