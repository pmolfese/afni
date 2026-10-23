/* Unit tests for the file syntax handled by thd_datatable.c. */

#include "mrilib.h"
#include "thd_datatable.h"

#define CHECK(COND,MSG) do{ if(!(COND)){ fprintf(stderr,"FAIL: %s\n",MSG); return 1; } }while(0)

static int write_table( const char *fn )
{
   FILE *fp=fopen(fn,"w") ;
   if( fp==NULL ) return 1 ;
   fprintf(fp,
     "# a full-line comment\n"
     "Subj Run InputFile Note\n"
     "s01 1 \"stats+tlrc[face#0_Coef]\" \"quoted # value\" # trailing comment\n"
     "s02 2 stats+tlrc[house#1_Coef] plain#value\n"
     "s03 3 stats+tlrc[cat#2_Coef] \"#not-a-comment\"\n") ;
   fclose(fp) ;
   return 0 ;
}

static int check_table( THD_datatable *dt )
{
   int inote ;
   CHECK(dt!=NULL,"read the table") ;
   inote=THD_datatable_column(dt,"Note") ;
   CHECK(dt->nrow==3 && dt->ncol==4,"comments do not change the table shape") ;
   CHECK(inote==3,"find the Note column") ;
   CHECK(!strcmp(dt->fname[0],"stats+tlrc[face#0_Coef]"),
         "a hash inside a quoted AFNI selector is data") ;
   CHECK(!strcmp(dt->fname[1],"stats+tlrc[house#1_Coef]"),
         "a hash inside an unquoted AFNI selector is data") ;
   CHECK(!strcmp(DT_CELL(dt,0,inote),"quoted # value"),
         "quoted whitespace and hash are data") ;
   CHECK(!strcmp(DT_CELL(dt,1,inote),"plain#value"),
         "an embedded unquoted hash is data") ;
   CHECK(!strcmp(DT_CELL(dt,2,inote),"#not-a-comment"),
         "a quoted value may begin with a hash") ;
   return 0 ;
}

int main( int argc , char **argv )
{
   char fn[512],printed[512],*dir=getenv("TMPDIR") ;
   THD_datatable *dt,*roundtrip ; FILE *fp ;

   snprintf(fn,sizeof(fn),"%s/thd_datatable_unit_%d.txt",
            dir?dir:".",(int)getpid()) ;
   snprintf(printed,sizeof(printed),"%s/thd_datatable_roundtrip_%d.txt",
            dir?dir:".",(int)getpid()) ;
   CHECK(write_table(fn)==0,"write the test table") ;

   dt=THD_read_datatable_file(fn) ;
   CHECK(check_table(dt)==0,"parse hashes and comments") ;

   fp=fopen(printed,"w") ;
   CHECK(fp!=NULL,"open the round-trip table") ;
   THD_datatable_print(dt,fp) ; fclose(fp) ;
   roundtrip=THD_read_datatable_file(printed) ;
   CHECK(check_table(roundtrip)==0,"printed tables preserve hashes") ;

   THD_free_datatable(roundtrip) ; THD_free_datatable(dt) ;
   remove(printed) ; remove(fn) ;
   (void)argc ; (void)argv ;
   printf("thd_datatable tests passed\n") ;
   return 0 ;
}
