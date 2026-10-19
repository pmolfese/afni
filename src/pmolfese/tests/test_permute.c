#include "mrilib.h"
#include "thd_permute.h"

/* Nonascending, repeated labels in BOTH blocks. The literal identity is not
   the sorted-key representative used internally to enumerate the group. */
static int check_set( int exchange, int blockmode, int equivalent,
                      int exact, int expected )
{
   int blocks[6]={0,0,0,1,1,1},labels[6]={2,0,2,1,2,1};
   PERM_scheme *scheme=THD_perm_scheme_new(6);
   PERM_set *set=NULL;
   int ip,jp,i,fail=1;
   if(!scheme) return 1;
   scheme->exchange=exchange;scheme->blockmode=blockmode;
   scheme->exact=exact;scheme->nperm=exact?1000:2;scheme->seed=42;
   if(THD_perm_scheme_set_blocks(scheme,blocks)!=2) goto done;
   if(equivalent && !THD_perm_scheme_set_eqclass(scheme,labels)) goto done;
   set=THD_perm_set_build(scheme);
   if(!set || set->nperm!=expected || set->is_exact!=exact) goto done;
   for(i=0;i<6;i++) if(set->perm[i]!=i || set->sign[i]!=1) goto done;
   for(ip=0;ip<set->nperm;ip++){
      int seen[6]={0};
      for(i=0;i<6;i++){
         int src=set->perm[ip*6+i],sign=set->sign[ip*6+i];
         if(src<0 || src>=6 || seen[src]++ || (sign!=1 && sign!=-1)) goto done;
         if(blockmode==PERM_WITHIN_BLOCK && blocks[src]!=blocks[i]) goto done;
         if(exchange==PERM_EE && sign!=1) goto done;
         if(exchange==PERM_ISE && src!=i) goto done;
      }
      if(!exact) continue; /* random draws may repeat */
      for(jp=0;jp<ip;jp++){
         for(i=0;i<6;i++){
            int a=set->perm[ip*6+i],b=set->perm[jp*6+i];
            if((equivalent?labels[a]!=labels[b]:a!=b) ||
               set->sign[ip*6+i]!=set->sign[jp*6+i]) break;
         }
         if(i==6) goto done;
      }
   }
   fail=0;
done:
   if(fail) fprintf(stderr,"FAIL exact identity/coverage: exchange=%d blocks=%d eq=%d exact=%d\n",
                    exchange,blockmode,equivalent,exact);
   THD_perm_set_free(set);THD_perm_scheme_free(scheme);
   return fail;
}

int main(void)
{
   int fail=0;
   fail+=check_set(PERM_EE,PERM_WITHIN_BLOCK,1,1,9);
   fail+=check_set(PERM_ISE,PERM_WITHIN_BLOCK,1,1,64);
   fail+=check_set(PERM_BOTH,PERM_WITHIN_BLOCK,1,1,576);
   fail+=check_set(PERM_EE,PERM_WITHIN_BLOCK,0,1,36);
   fail+=check_set(PERM_EE,PERM_WHOLE_BLOCK,1,1,2);
   fail+=check_set(PERM_BOTH,PERM_WHOLE_BLOCK,1,1,8);
   fail+=check_set(PERM_EE,PERM_BOTH_BLOCK,1,1,18);
   fail+=check_set(PERM_EE,PERM_WITHIN_BLOCK,1,0,2);
   if(fail) return 1;
   puts("PASS permutation identity, bijections, and exact distinct-label coverage");
   return 0;
}
