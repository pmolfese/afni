#include "thd_rsa.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

static int closef( float aa, float bb )
{
   return fabsf(aa-bb)<1.0e-5f ;
}

static int condition_rdm_test( void )
{
   float raw[6]={1,10, 3,14, 5,16} ;
   float centered[6],want[6]={-2,-10.0f/3.0f, 0,2.0f/3.0f, 2,8.0f/3.0f} ;
   float euclid[6] ; THD_simmat *sm=NULL ; char err[256] ; int ii,fail=0 ;

   memcpy(centered,raw,sizeof(raw)) ;
   if( !THD_rsa_center_condition_patterns(3,2,centered,err,sizeof(err)) ){
     fprintf(stderr,"FAIL condition centering: %s\n",err) ; return 1 ;
   }
   for( ii=0 ; ii<6 ; ii++ ) if( !closef(centered[ii],want[ii]) ){
     fprintf(stderr,"FAIL centered pattern %d: got %.9g, wanted %.9g\n",
             ii,centered[ii],want[ii]) ; fail=1 ;
   }

   memcpy(euclid,raw,sizeof(raw)) ;
   sm=THD_rsa_condition_rdm(3,2,euclid,SIM_EUCLID,1,err,sizeof(err)) ;
   if( sm==NULL ){
     fprintf(stderr,"FAIL Euclidean condition RDM: %s\n",err) ; fail=1 ;
   } else {
     if( memcmp(euclid,raw,sizeof(raw))!=0 ){
       fprintf(stderr,"FAIL Euclidean condition centering was not an exact no-op\n") ; fail=1 ;
     }
     if( !sm->is_dist || !closef(sm->mat[1],sqrtf(20.0f)) ||
         !closef(sm->mat[2],sqrtf(52.0f)) ){
       fprintf(stderr,"FAIL Euclidean condition RDM values\n") ; fail=1 ;
     }
     THD_simmat_free(sm) ;
   }
   return fail ;
}

static THD_simmat * triangle_matrix( const float tri[3], int is_dist )
{
   THD_simmat *sm=THD_simmat_new(3) ;
   if( sm==NULL ) return NULL ;
   THD_tri_to_simmat(3,(float *)tri,sm) ; sm->is_dist=is_dist ; return sm ;
}

static int second_order_test( void )
{
   const float s0[3]={1.0f,0.5f,0.0f} ;
   const float s1[3]={1.0f,0.0f,-1.0f} ;
   const float s2[3]={0.5f,0.5f,0.5f} ;
   const float *src[3]={s0,s1,s2} ;
   float all[9],tri[3],meanwant[3] ;
   THD_simmat *inner[3]={NULL,NULL,NULL},*outer=NULL,*outer_sim=NULL,*mean=NULL ;
   THD_rsa_rdm_workspace *ws=NULL ; char err[256]={0} ; int ss,aa,fail=0 ;

   ws=THD_rsa_rdm_workspace_new(4,5) ;
   if( ws==NULL || !THD_rsa_rdm_workspace_reset(ws,3,3,err,sizeof(err)) ){
     fprintf(stderr,"FAIL RSA RDM workspace: %s\n",err) ; fail=1 ; goto done ;
   }
   for( ss=0 ; ss<3 ; ss++ ){
     inner[ss]=triangle_matrix(src[ss],0) ;
     if( inner[ss]==NULL ||
         !THD_rsa_simmat_to_triangle(inner[ss],tri,1,err,sizeof(err)) ||
         !THD_rsa_rdm_workspace_store(ws,ss,inner[ss],1,err,sizeof(err)) ){
       fprintf(stderr,"FAIL storing subject RDM %d: %s\n",ss+1,err) ; fail=1 ; goto done ;
     }
     for( aa=0 ; aa<3 ; aa++ ){
       float want=1.0f-src[ss][aa] ;
       if( !closef(tri[aa],want) ){
         fprintf(stderr,"FAIL canonical distance subject %d entry %d\n",ss+1,aa+1) ;
         fail=1 ;
       }
       all[ss*3+aa]=src[ss][aa] ;
     }
   }
   outer=THD_rsa_second_order_rdm(ws,SIM_EUCLID,err,sizeof(err)) ;
   if( outer==NULL ){
     fprintf(stderr,"FAIL second-order RDM: %s\n",err) ; fail=1 ;
   } else if( !closef(outer->mat[1],sqrtf(1.25f)) ||
              !closef(outer->mat[2],sqrtf(0.5f)) ||
              !closef(outer->mat[5],sqrtf(2.75f)) ){
     fprintf(stderr,"FAIL second-order Euclidean geometry: %.9g %.9g %.9g\n",
             outer->mat[1],outer->mat[2],outer->mat[5]) ; fail=1 ;
   }

   if( !THD_rsa_rdm_workspace_reset(ws,3,3,err,sizeof(err)) ){
     fprintf(stderr,"FAIL resetting RSA RDM workspace: %s\n",err) ; fail=1 ;
   } else {
     for( ss=0 ; ss<3 ; ss++ )
       if( !THD_rsa_rdm_workspace_store(ws,ss,inner[ss],0,err,sizeof(err)) ){
         fprintf(stderr,"FAIL storing similarity RDM %d: %s\n",ss+1,err) ; fail=1 ;
       }
     outer_sim=THD_rsa_second_order_rdm(ws,SIM_EUCLID,err,sizeof(err)) ;
     if( outer_sim==NULL || !closef(outer_sim->mat[1],sqrtf(1.25f)) ||
         !closef(outer_sim->mat[2],sqrtf(0.5f)) ||
         !closef(outer_sim->mat[5],sqrtf(2.75f)) ){
       fprintf(stderr,"FAIL reusable similarity-to-distance scratch conversion\n") ; fail=1 ;
     }
   }

   for( aa=0 ; aa<3 ; aa++ )
     meanwant[aa]=(s0[aa]+s1[aa]+s2[aa])/3.0f ;
   mean=THD_rsa_mean_subject_rdms(3,3,all,0,err,sizeof(err)) ;
   if( mean==NULL ){
     fprintf(stderr,"FAIL mean subject RDM: %s\n",err) ; fail=1 ;
   } else {
     THD_simmat_to_tri(mean,tri) ;
     if( mean->is_dist ){
       fprintf(stderr,"FAIL mean subject RDM changed similarity sense\n") ; fail=1 ;
     }
     for( aa=0 ; aa<3 ; aa++ ) if( !closef(tri[aa],meanwant[aa]) ){
       fprintf(stderr,"FAIL mean subject RDM entry %d\n",aa+1) ; fail=1 ;
     }
   }

done:
   for( ss=0 ; ss<3 ; ss++ ) THD_simmat_free(inner[ss]) ;
   THD_simmat_free(outer) ; THD_simmat_free(outer_sim) ; THD_simmat_free(mean) ;
   THD_rsa_rdm_workspace_free(ws) ; return fail ;
}

static int masked_rdm_test( void )
{
   char in[]="/tmp/afni_masked_rdm_in_XXXXXX",out[]="/tmp/afni_masked_rdm_out_XXXXXX" ;
   int fd1=mkstemp(in),fd2=mkstemp(out),fail=0 ; FILE *fp ;
   THD_simmat *sm=NULL,*back=NULL ; THD_rdm_ws *ws=NULL ;
   float neural[6]={0.0f,1.0f,9.0f,1.0f,9.0f,9.0f},r ; int perm[4]={1,2,3,0} ;
   if( fd1<0 || fd2<0 ){ fprintf(stderr,"FAIL masked RDM temp files\n") ; return 1 ; }
   close(fd1) ; close(fd2) ;
   fp=fopen(in,"w") ;
   fprintf(fp,"0 0 1 .\n0 0 1 .\n1 1 0 .\n. . . 0\n") ; fclose(fp) ;
   sm=THD_simmat_read_1D(in,4) ;
   if( sm==NULL || !THD_simmat_has_missing(sm) || THD_simmat_nvalid(sm)!=3 ){
     fprintf(stderr,"FAIL reading dot-masked RDM\n") ; fail=1 ; goto done ;
   }
   ws=THD_rdm_ws_new(4,1) ; r=THD_simmat_corr_triangle(neural,sm,NULL,CMP_PEARSON,ws) ;
   if( !closef(r,1.0f) ){
     fprintf(stderr,"FAIL masked RDM finite-pair correlation: %.9g\n",r) ; fail=1 ;
   }
   r=THD_simmat_corr_triangle(neural,sm,perm,CMP_PEARSON,ws) ;
   if( !closef(r,-0.5f) ){
     fprintf(stderr,"FAIL masked RDM row/column relabeling: %.9g\n",r) ; fail=1 ;
   }
   if( !THD_simmat_write_1D(out,sm) ){
     fprintf(stderr,"FAIL writing dot-masked RDM\n") ; fail=1 ; goto done ;
   }
   back=THD_simmat_read_1D(out,4) ;
   if( back==NULL || THD_simmat_nvalid(back)!=3 || !isnan(back->mat[3]) ||
       !isnan(back->mat[12]) ){
     fprintf(stderr,"FAIL dot-masked RDM round trip\n") ; fail=1 ;
   }
done:
   THD_simmat_free(sm) ; THD_simmat_free(back) ; THD_rdm_ws_free(ws) ;
   unlink(in) ; unlink(out) ; return fail ;
}

int main( void )
{
   int fail=condition_rdm_test()+second_order_test()+masked_rdm_test() ;
   if( fail ) return 1 ;
   puts("PASS thd_rsa condition, subject, mean, and second-order RDM construction") ;
   return 0 ;
}
