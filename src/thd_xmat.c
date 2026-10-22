#include "mrilib.h"
#include "thd_xmat.h"
#include <errno.h>

/*----------------------------------------------------------------------------
  Design matrices (X.xmat.1D).  See thd_xmat.h for the format and the API.
                                                     -- P Molfese, Oct 2026
------------------------------------------------------------------------------*/

static void xm_err( char *err , int nerr , const char *fmt , const char *a )
{
   if( err != NULL && nerr > 0 ) snprintf(err,(size_t)nerr,fmt,a) ;
}

/*! A ';'-separated string attribute, each entry trimmed. */
static NI_str_array * xm_strings( NI_element *nel , const char *name )
{
   char *a = NI_get_attribute_nocase(nel,(char *)name) ; NI_str_array *sar ; int ii ;
   if( a == NULL ) return NULL ;
   sar = NI_decode_string_list(a,";") ;
   if( sar == NULL ) return NULL ;
   for( ii=0 ; ii < sar->num ; ii++ ){
     char *p = sar->str[ii] , *e ;
     while( isspace((unsigned char)*p) ) p++ ;
     memmove(sar->str[ii],p,strlen(p)+1) ;
     e = sar->str[ii] + strlen(sar->str[ii]) ;
     while( e > sar->str[ii] && isspace((unsigned char)e[-1]) ) *(--e) = '\0' ;
   }
   return sar ;
}

static int * xm_ints( NI_element *nel , const char *name , int *n )
{
   char *a = NI_get_attribute_nocase(nel,(char *)name) ; NI_int_array *ia ; int *v ;
   *n = 0 ;
   if( a == NULL ) return NULL ;
   ia = NI_decode_int_list(a,";,") ;
   if( ia == NULL ) return NULL ;
   *n = ia->num ; v = (int *)malloc(sizeof(int)*(size_t)(ia->num>0?ia->num:1)) ;
   memcpy(v,ia->ar,sizeof(int)*(size_t)ia->num) ;
   NI_delete_int_array(ia) ;
   return v ;
}

THD_xmat * THD_xmat_read( const char *fname , char *err , int nerr )
{
   NI_element *nel ; THD_xmat *xm ; NI_str_array *sar ; int ii , c , r , nb , nt ;
   int *bot , *top ; char *a ;

   if( fname == NULL ){ xm_err(err,nerr,"%s","no file name") ; return NULL ; }
   nel = (NI_element *)NI_read_element_fromfile((char *)fname) ;
   if( nel == NULL || nel->type != NI_ELEMENT_TYPE ){
     xm_err(err,nerr,"cannot read '%s' as an X.xmat.1D (a NIML matrix from 3dDeconvolve)",fname) ;
     if( nel != NULL ) NI_free_element(nel) ;
     return NULL ;
   }
   xm = (THD_xmat *)calloc(1,sizeof(THD_xmat)) ;
   xm->fname = strdup(fname) ;
   xm->ncol = nel->vec_num ; xm->nrow = nel->vec_len ;
   if( xm->ncol < 1 || xm->nrow < 1 ){
     xm_err(err,nerr,"'%s' holds no matrix",fname) ; NI_free_element(nel) ; THD_xmat_free(xm) ; return NULL ;
   }

   /* the matrix, as float, column by column */
   xm->mat = (float *)malloc(sizeof(float)*(size_t)xm->ncol*(size_t)xm->nrow) ;
   for( c=0 ; c < xm->ncol ; c++ ){
     float *col = xm->mat + (size_t)c*xm->nrow ;
     for( r=0 ; r < xm->nrow ; r++ ){
       switch( nel->vec_typ[c] ){
         case NI_FLOAT:  col[r] = ((float  *)nel->vec[c])[r] ; break ;
         case NI_DOUBLE: col[r] = (float)((double *)nel->vec[c])[r] ; break ;
         case NI_INT:    col[r] = (float)((int    *)nel->vec[c])[r] ; break ;
         default:
           xm_err(err,nerr,"'%s' has a matrix column that is not numeric",fname) ;
           NI_free_element(nel) ; THD_xmat_free(xm) ; return NULL ;
       }
     }
   }

   /* every header attribute, as written */
   xm->nattr = nel->attr_num ;
   xm->attr_name  = (char **)calloc((size_t)(xm->nattr>0?xm->nattr:1),sizeof(char *)) ;
   xm->attr_value = (char **)calloc((size_t)(xm->nattr>0?xm->nattr:1),sizeof(char *)) ;
   for( ii=0 ; ii < xm->nattr ; ii++ ){
     xm->attr_name[ii]  = strdup(nel->attr_lhs[ii] ? nel->attr_lhs[ii] : "") ;
     xm->attr_value[ii] = strdup(nel->attr_rhs[ii] ? nel->attr_rhs[ii] : "") ;
   }

   /* column labels */
   xm->collab = (char **)calloc((size_t)xm->ncol,sizeof(char *)) ;
   sar = xm_strings(nel,"ColumnLabels") ;
   for( c=0 ; c < xm->ncol ; c++ ){
     char b[32] ;
     if( sar != NULL && c < sar->num ) xm->collab[c] = strdup(sar->str[c]) ;
     else { snprintf(b,sizeof(b),"col#%d",c) ; xm->collab[c] = strdup(b) ; }
   }
   if( sar != NULL ) NI_delete_str_array(sar) ;

   /* the timeline: full length, the rows kept, the runs */
   a = NI_get_attribute_nocase(nel,"NRowFull") ;
   xm->nfull = (a != NULL) ? (int)strtod(a,NULL) : xm->nrow ;
   xm->good = xm_ints(nel,"GoodList",&nt) ;
   if( xm->good == NULL ){
     xm->good = (int *)malloc(sizeof(int)*(size_t)xm->nrow) ;
     for( r=0 ; r < xm->nrow ; r++ ) xm->good[r] = r ;
     nt = xm->nrow ;
   }
   if( nt != xm->nrow || xm->nfull < xm->nrow ){
     xm_err(err,nerr,"'%s': GoodList and NRowFull do not match the matrix",fname) ;
     NI_free_element(nel) ; THD_xmat_free(xm) ; return NULL ;
   }
   for( r=0 ; r < xm->nrow ; r++ )
     if( xm->good[r] < 0 || xm->good[r] >= xm->nfull || (r > 0 && xm->good[r] <= xm->good[r-1]) ){
       xm_err(err,nerr,"'%s': GoodList has bad values",fname) ;
       NI_free_element(nel) ; THD_xmat_free(xm) ; return NULL ;
     }
   xm->runstart = xm_ints(nel,"RunStart",&xm->nrun) ;
   if( xm->runstart == NULL || xm->nrun < 1 ){
     free(xm->runstart) ; xm->runstart = (int *)calloc(1,sizeof(int)) ; xm->nrun = 1 ;
   }
   for( ii=0 ; ii < xm->nrun ; ii++ )
     if( xm->runstart[ii] < 0 || xm->runstart[ii] >= xm->nfull || (ii == 0 && xm->runstart[ii] != 0)
         || (ii > 0 && xm->runstart[ii] <= xm->runstart[ii-1]) ){
       xm_err(err,nerr,"'%s': RunStart has bad values",fname) ;
       NI_free_element(nel) ; THD_xmat_free(xm) ; return NULL ;
     }
   a = NI_get_attribute_nocase(nel,"RowTR") ;
   xm->tr = (a != NULL) ? (float)strtod(a,NULL) : 0.0f ;

   /* the stimuli */
   sar = xm_strings(nel,"StimLabels") ;
   bot = xm_ints(nel,"StimBots",&nb) ; top = xm_ints(nel,"StimTops",&nt) ;
   if( sar != NULL && bot != NULL && top != NULL && nb == sar->num && nt == sar->num ){
     xm->nstim = sar->num ;
     xm->stim = (THD_xmat_stim *)calloc((size_t)(xm->nstim>0?xm->nstim:1),sizeof(THD_xmat_stim)) ;
     for( ii=0 ; ii < xm->nstim ; ii++ ){
       char an[64] ;
       xm->stim[ii].label = strdup(sar->str[ii]) ;
       xm->stim[ii].bot = bot[ii] ; xm->stim[ii].top = top[ii] ;
       if( bot[ii] < 0 || top[ii] < bot[ii] || top[ii] >= xm->ncol ){
         xm_err(err,nerr,"'%s': a stimulus has columns outside the matrix",fname) ;
         NI_delete_str_array(sar) ; free(bot) ; free(top) ; NI_free_element(nel) ; THD_xmat_free(xm) ; return NULL ;
       }
       snprintf(an,sizeof(an),"BasisOption_%06d",ii+1) ; a = NI_get_attribute_nocase(nel,an) ;
       xm->stim[ii].option = (a != NULL) ? strdup(a) : NULL ;
       snprintf(an,sizeof(an),"BasisFormula_%06d",ii+1) ; a = NI_get_attribute_nocase(nel,an) ;
       xm->stim[ii].formula = (a != NULL) ? strdup(a) : NULL ;
     }
   }
   if( sar != NULL ) NI_delete_str_array(sar) ;
   free(bot) ; free(top) ;
   NI_free_element(nel) ;
   return xm ;
}

void THD_xmat_free( THD_xmat *xm )
{
   int ii ;
   if( xm == NULL ) return ;
   free(xm->fname) ; free(xm->mat) ; free(xm->good) ; free(xm->runstart) ;
   if( xm->collab != NULL ){ for( ii=0 ; ii < xm->ncol ; ii++ ) free(xm->collab[ii]) ; free(xm->collab) ; }
   if( xm->stim != NULL ){
     for( ii=0 ; ii < xm->nstim ; ii++ ){ free(xm->stim[ii].label) ; free(xm->stim[ii].option) ; free(xm->stim[ii].formula) ; }
     free(xm->stim) ;
   }
   if( xm->attr_name != NULL ){
     for( ii=0 ; ii < xm->nattr ; ii++ ){ free(xm->attr_name[ii]) ; free(xm->attr_value[ii]) ; }
     free(xm->attr_name) ; free(xm->attr_value) ;
   }
   free(xm) ;
}

const char * THD_xmat_attr( const THD_xmat *xm , const char *name )
{
   int ii ;
   if( xm == NULL || name == NULL ) return NULL ;
   for( ii=0 ; ii < xm->nattr ; ii++ ) if( strcasecmp(xm->attr_name[ii],name) == 0 ) return xm->attr_value[ii] ;
   return NULL ;
}

int THD_xmat_run_of( const THD_xmat *xm , int t )
{
   int r ;
   if( xm == NULL || t < 0 || t >= xm->nfull ) return -1 ;
   for( r=0 ; r+1 < xm->nrun && t >= xm->runstart[r+1] ; r++ ) ;
   return r ;
}

int THD_xmat_run_length( const THD_xmat *xm , int r )
{
   if( xm == NULL || r < 0 || r >= xm->nrun ) return 0 ;
   return ((r+1 < xm->nrun) ? xm->runstart[r+1] : xm->nfull) - xm->runstart[r] ;
}

unsigned char * THD_xmat_censored( const THD_xmat *xm )
{
   unsigned char *c ; int r ;
   if( xm == NULL ) return NULL ;
   c = (unsigned char *)malloc((size_t)xm->nfull) ; memset(c,1,(size_t)xm->nfull) ;
   for( r=0 ; r < xm->nrow ; r++ ) c[xm->good[r]] = 0 ;
   return c ;
}

int THD_xmat_ncensored( const THD_xmat *xm )
{
   return (xm == NULL) ? 0 : xm->nfull - xm->nrow ;
}

unsigned char * THD_xmat_stims_wanted( const THD_xmat *xm , const char *list , int *nfound )
{
   unsigned char *w ; int s , nf=0 ;
   if( xm == NULL ) return NULL ;
   w = (unsigned char *)calloc((size_t)(xm->nstim>0?xm->nstim:1),1) ;
   for( s=0 ; s < xm->nstim ; s++ ){
     if( list == NULL ) w[s] = 1 ;
     else {
       char *copy=strdup(list) , *tok ;
       for( tok=strtok(copy,",") ; tok != NULL && !w[s] ; tok=strtok(NULL,",") ) w[s] = (strcmp(tok,xm->stim[s].label) == 0) ;
       free(copy) ;
     }
     nf += w[s] ;
   }
   if( nfound != NULL ) *nfound = nf ;
   return w ;
}

int * THD_xmat_condition_labels( const THD_xmat *xm , double frac , const unsigned char *wanted )
{
   int *which , s , c , r , t ; double *best ;
   if( xm == NULL ) return NULL ;
   which = (int *)malloc(sizeof(int)*(size_t)xm->nfull) ; best = (double *)malloc(sizeof(double)*(size_t)xm->nfull) ;
   for( t=0 ; t < xm->nfull ; t++ ){ which[t] = -1 ; best[t] = -1.0 ; }
   for( s=0 ; s < xm->nstim ; s++ ){
     double mx=0.0 ;
     if( wanted != NULL && !wanted[s] ) continue ;
     for( c=xm->stim[s].bot ; c <= xm->stim[s].top ; c++ )
       for( r=0 ; r < xm->nrow ; r++ ){ double v=fabs(xm->mat[(size_t)c*xm->nrow+r]) ; if( v > mx ) mx = v ; }
     if( mx <= 0.0 ) continue ;
     for( r=0 ; r < xm->nrow ; r++ ){
       double v=0.0 ; t = xm->good[r] ;
       for( c=xm->stim[s].bot ; c <= xm->stim[s].top ; c++ ){ double u=fabs(xm->mat[(size_t)c*xm->nrow+r])/mx ; if( u > v ) v = u ; }
       if( v >= frac && v > best[t] ){ best[t] = v ; which[t] = s ; }
     }
   }
   free(best) ;
   return which ;
}

THD_xmat_item * THD_xmat_items( const THD_xmat *xm , const unsigned char *wanted , int *nitem )
{
   THD_xmat_item *it ; int s , c , r , n=0 ;
   *nitem = 0 ;
   if( xm == NULL ) return NULL ;
   it = (THD_xmat_item *)calloc((size_t)(xm->ncol>0?xm->ncol:1),sizeof(THD_xmat_item)) ;
   for( s=0 ; s < xm->nstim ; s++ ){
     if( wanted != NULL && !wanted[s] ) continue ;
     for( c=xm->stim[s].bot ; c <= xm->stim[s].top ; c++ ){
       THD_xmat_item *q = it + n ; char b[THD_MAX_NAME] ; double mx=0.0 ; int run=-1 , one=1 ;
       q->key = strdup( (xm->stim[s].top > xm->stim[s].bot) ? xm->collab[c] : xm->stim[s].label ) ;
       q->condition = strdup(xm->stim[s].label) ;
       snprintf(b,sizeof(b),"%s_Coef",xm->collab[c]) ; q->coef_label = strdup(b) ;
       q->stim = s ; q->column = c ;
       for( r=0 ; r < xm->nrow ; r++ ){ double v=fabs(xm->mat[(size_t)c*xm->nrow+r]) ; if( v > mx ) mx = v ; }
       for( r=0 ; r < xm->nrow && one ; r++ ) if( fabs(xm->mat[(size_t)c*xm->nrow+r]) > 1.0e-6*mx ){
         int rr = THD_xmat_run_of(xm,xm->good[r]) ;
         if( run < 0 ) run = rr ; else if( rr != run ) one = 0 ;
       }
       q->run = (one && mx > 0.0 && xm->nrun > 1) ? run : -1 ;
       n++ ;
     }
   }
   *nitem = n ;
   return it ;
}

void THD_xmat_items_free( THD_xmat_item *it , int nitem )
{
   int ii ;
   if( it == NULL ) return ;
   for( ii=0 ; ii < nitem ; ii++ ){ free(it[ii].key) ; free(it[ii].condition) ; free(it[ii].coef_label) ; }
   free(it) ;
}

/*============================================================================*/
/*  Tables with an Xmat column                                                */
/*============================================================================*/

typedef struct { int ntok , ncap ; char **tok ; } xm_tokens ;

static void xm_add( xm_tokens *t , const char *s )
{
   if( t->ntok >= t->ncap ){ t->ncap = 2*t->ncap + 64 ; t->tok = (char **)realloc(t->tok,sizeof(char *)*(size_t)t->ncap) ; }
   t->tok[t->ntok++] = strdup(s) ;
}

static THD_datatable * xm_table( xm_tokens *t , char *source )
{
   THD_datatable *out = THD_parse_datatable(t->tok,t->ntok,source) ; int ii ;
   for( ii=0 ; ii < t->ntok ; ii++ ) free(t->tok[ii]) ;
   free(t->tok) ;
   return out ;
}

static THD_xmat * xm_row_xmat( THD_datatable *dt , int row , int icol )
{
   const char *fn = DT_CELL(dt,row,icol) ; char err[512] ; THD_xmat *xm ;
   if( fn[0] == '\0' || strcmp(fn,"-") == 0 || strcasecmp(fn,"NA") == 0 )
     ERROR_exit("Subj %s: no Xmat on row %d of the table",dt->subj[row],row+1) ;
   xm = THD_xmat_read(fn,err,sizeof(err)) ;
   if( xm == NULL ) ERROR_exit("Subj %s: %s",dt->subj[row],err) ;
   return xm ;
}

static int xm_manifest_index( const char *s , const char *name , int row , int limit )
{
   char *end ; long v ;
   errno = 0 ; v = strtol(s,&end,10) ;
   if( errno || end == s || *end || v < 0 || v >= limit )
     ERROR_exit("trial manifest row %d: %s '%s' must be a zero-based integer in 0..%d",
                row+1,name,s,limit-1) ;
   return (int)v ;
}

THD_datatable * THD_datatable_trial_manifest( THD_datatable *dt )
{
   xm_tokens t = {0,0,NULL} ; int iv,ix,ic,r,c ;
   THD_3dim_dataset *ds=NULL ; THD_xmat *xm=NULL ;
   char *last_input=NULL,*last_xmat=NULL ;
   if( dt == NULL ) return NULL ;
   iv = THD_datatable_column(dt,"Volume") ;
   ix = THD_datatable_column(dt,"Xmat") ;
   ic = THD_datatable_column(dt,"XColumn") ;
   if( dt->icol_input < 0 || iv < 0 )
     ERROR_exit("trial manifest requires InputFile and Volume columns") ;
   if( (ix < 0) != (ic < 0) )
     ERROR_exit("trial manifest: optional Xmat and XColumn columns must occur together") ;
   /* The canonical parser requires InputFile last. Keep provenance as metadata. */
   for( c=0 ; c < dt->ncol ; c++ ) if( c != dt->icol_input ) xm_add(&t,dt->cname[c]) ;
   xm_add(&t,"InputFile") ;
   for( r=0 ; r < dt->nrow ; r++ ){
     const char *fn = DT_CELL(dt,r,dt->icol_input) ; int b ; char *sel ;
     if( !*fn || strchr(fn,'[') || !strcmp(fn,"-") || !strcasecmp(fn,"NA") )
       ERROR_exit("trial manifest row %d: InputFile must name a whole beta bucket",r+1) ;
     if( last_input == NULL || strcmp(last_input,fn) ){
       if( ds != NULL ) DSET_delete(ds) ;
       ds = THD_open_dataset((char *)fn) ; CHECK_OPEN_ERROR(ds,fn) ;
       last_input = (char *)fn ;
     }
     b = xm_manifest_index(DT_CELL(dt,r,iv),"Volume",r,DSET_NVALS(ds)) ;
     if( ix >= 0 ){
       const char *xf = DT_CELL(dt,r,ix) ; char err[512] ;
       if( last_xmat == NULL || strcmp(last_xmat,xf) ){
         THD_xmat_free(xm) ; xm = THD_xmat_read(xf,err,sizeof(err)) ;
         if( xm == NULL ) ERROR_exit("trial manifest row %d: %s",r+1,err) ;
         last_xmat = (char *)xf ;
       }
       (void)xm_manifest_index(DT_CELL(dt,r,ic),"XColumn",r,xm->ncol) ;
     }
     for( c=0 ; c < dt->ncol ; c++ ) if( c != dt->icol_input ) xm_add(&t,DT_CELL(dt,r,c)) ;
     sel = (char *)malloc(strlen(fn)+32) ; sprintf(sel,"%s[%d]",fn,b) ;
     xm_add(&t,sel) ; free(sel) ;
   }
   if( ds != NULL ) DSET_delete(ds) ;
   THD_xmat_free(xm) ;
   return xm_table(&t,"trial manifest") ;
}

THD_datatable * THD_datatable_expand_xmat_items( THD_datatable *dt , const char *xmat_col ,
                                                 const char *item_col , const char *run_col ,
                                                 const char *stims , int *have_run_out )
{
   int icol , r , c , j , ncarry=0 , *carry , have_run=1 , *nit ; THD_xmat_item **it ;
   xm_tokens t = {0,0,NULL} ; THD_datatable *out ;
   if( dt == NULL || xmat_col == NULL ) return NULL ;
   icol = THD_datatable_column(dt,(char *)xmat_col) ;
   if( icol < 0 ) ERROR_exit("the table has no '%s' column",xmat_col) ;
   if( dt->icol_input < 0 ) ERROR_exit("the data table needs an InputFile column") ;
   if( dt->icol_subj < 0 ) ERROR_exit("a table with an %s column needs a Subj column",xmat_col) ;
   for( c=0 ; c < dt->ncol ; c++ )
     if( !strcasecmp(dt->cname[c],item_col) || !strcasecmp(dt->cname[c],"Condition") || !strcasecmp(dt->cname[c],run_col) )
       ERROR_exit("with an %s column the items, their conditions, and runs come from the design matrix: remove the '%s' column",
                  xmat_col,dt->cname[c]) ;
   carry = (int *)malloc(sizeof(int)*(size_t)dt->ncol) ;
   for( c=0 ; c < dt->ncol ; c++ ) if( c != icol && c != dt->icol_input ) carry[ncarry++] = c ;
   nit = (int *)calloc((size_t)dt->nrow,sizeof(int)) ;
   it  = (THD_xmat_item **)calloc((size_t)dt->nrow,sizeof(THD_xmat_item *)) ;

   /* each row: its items, and the bucket volume of each (found by label) */
   for( r=0 ; r < dt->nrow ; r++ ){
     const char *bn = DT_CELL(dt,r,dt->icol_input) ; THD_xmat *xm ; THD_3dim_dataset *ds ;
     unsigned char *want ; int nf ;
     if( strchr(bn,'[') != NULL )
       ERROR_exit("Subj %s: with an %s column InputFile is the whole stats bucket, without a selector ('%s')",dt->subj[r],xmat_col,bn) ;
     xm = xm_row_xmat(dt,r,icol) ;
     if( xm->nstim < 1 ) ERROR_exit("Subj %s: '%s' has no stimulus columns (StimLabels, StimBots, StimTops)",dt->subj[r],xm->fname) ;
     want = THD_xmat_stims_wanted(xm,stims,&nf) ;
     if( nf < 1 ){
       if( stims != NULL ) ERROR_exit("Subj %s: no stimulus of '%s' is in the list '%s'",dt->subj[r],xm->fname,stims) ;
       ERROR_exit("Subj %s: '%s' has no stimulus",dt->subj[r],xm->fname) ;
     }
     it[r] = THD_xmat_items(xm,want,&nit[r]) ;
     ds = THD_open_dataset((char *)bn) ; CHECK_OPEN_ERROR(ds,bn) ;
     for( j=0 ; j < nit[r] ; j++ ){
       int b ;
       for( b=0 ; b < DSET_NVALS(ds) ; b++ )
         if( DSET_BRICK_LABEL(ds,b) != NULL && strcmp(DSET_BRICK_LABEL(ds,b),it[r][j].coef_label) == 0 ) break ;
       if( b == DSET_NVALS(ds) )
         ERROR_exit("Subj %s: '%s' has no volume labelled %s; is it the stats bucket of the regression in '%s'?",
                    dt->subj[r],bn,it[r][j].coef_label,xm->fname) ;
       it[r][j].column = b ;                          /* from here: the bucket volume */
       if( it[r][j].run < 0 ) have_run = 0 ;
     }
     DSET_delete(ds) ; free(want) ; THD_xmat_free(xm) ;
   }

   /* the item table */
   for( j=0 ; j < ncarry ; j++ ) xm_add(&t,dt->cname[carry[j]]) ;
   xm_add(&t,item_col) ; xm_add(&t,"Condition") ;
   if( have_run ) xm_add(&t,run_col) ;
   xm_add(&t,"InputFile") ;
   for( r=0 ; r < dt->nrow ; r++ ) for( j=0 ; j < nit[r] ; j++ ){
     char sel[THD_MAX_NAME+32] , rn[24] ; int q ;
     for( q=0 ; q < ncarry ; q++ ) xm_add(&t,DT_CELL(dt,r,carry[q])) ;
     xm_add(&t,it[r][j].key) ; xm_add(&t,it[r][j].condition) ;
     if( have_run ){ snprintf(rn,sizeof(rn),"%d",it[r][j].run+1) ; xm_add(&t,rn) ; }
     snprintf(sel,sizeof(sel),"%s[%d]",DT_CELL(dt,r,dt->icol_input),it[r][j].column) ; xm_add(&t,sel) ;
   }
   for( r=0 ; r < dt->nrow ; r++ ) THD_xmat_items_free(it[r],nit[r]) ;
   free(it) ; free(nit) ; free(carry) ;
   out = xm_table(&t,"the items of the Xmat column") ;
   if( out == NULL ) ERROR_exit("could not build the item table from the %s column",xmat_col) ;
   if( have_run_out != NULL ) *have_run_out = have_run ;
   return out ;
}

THD_datatable * THD_datatable_add_xmat_columns( THD_datatable *dt , const char *xmat_col )
{
   static const char *nm[6] = { "Xmat.TR" , "Xmat.nrun" , "Xmat.ntime" , "Xmat.ncensored" , "Xmat.frac_censored" , "Xmat.nstim" } ;
   int icol , r , c , q ; xm_tokens t = {0,0,NULL} ; THD_datatable *out ;
   char **lastf=NULL ; double (*val)[6] ;
   if( dt == NULL || xmat_col == NULL ) return NULL ;
   icol = THD_datatable_column(dt,(char *)xmat_col) ;
   if( icol < 0 ) ERROR_exit("the table has no '%s' column",xmat_col) ;
   if( dt->icol_input < 0 ) ERROR_exit("the data table needs an InputFile column") ;
   for( q=0 ; q < 6 ; q++ ) if( THD_datatable_column(dt,(char *)nm[q]) >= 0 ) return dt ;   /* already there */
   val = (double (*)[6])malloc(sizeof(double)*6*(size_t)dt->nrow) ;
   lastf = (char **)calloc((size_t)dt->nrow,sizeof(char *)) ;
   for( r=0 ; r < dt->nrow ; r++ ){
     int prev=-1 ;
     for( q=0 ; q < r && prev < 0 ; q++ ) if( strcmp(DT_CELL(dt,q,icol),DT_CELL(dt,r,icol)) == 0 ) prev = q ;
     if( prev >= 0 ){ memcpy(val[r],val[prev],sizeof(double)*6) ; continue ; }  /* rows sharing a file */
     { THD_xmat *xm = xm_row_xmat(dt,r,icol) ;
       val[r][0] = xm->tr ; val[r][1] = xm->nrun ; val[r][2] = xm->nfull ; val[r][3] = THD_xmat_ncensored(xm) ;
       val[r][4] = (xm->nfull > 0) ? (double)THD_xmat_ncensored(xm)/xm->nfull : 0.0 ; val[r][5] = xm->nstim ;
       THD_xmat_free(xm) ; }
   }
   for( c=0 ; c < dt->ncol ; c++ ) if( c != dt->icol_input ) xm_add(&t,dt->cname[c]) ;
   for( q=0 ; q < 6 ; q++ ) xm_add(&t,nm[q]) ;
   xm_add(&t,"InputFile") ;
   for( r=0 ; r < dt->nrow ; r++ ){
     char b[64] ;
     for( c=0 ; c < dt->ncol ; c++ ) if( c != dt->icol_input ) xm_add(&t,DT_CELL(dt,r,c)) ;
     for( q=0 ; q < 6 ; q++ ){ snprintf(b,sizeof(b),"%.6g",val[r][q]) ; xm_add(&t,b) ; }
     xm_add(&t,DT_CELL(dt,r,dt->icol_input)) ;
   }
   free(val) ; free(lastf) ;
   out = xm_table(&t,dt->source != NULL ? dt->source : "the table") ;
   if( out == NULL ) ERROR_exit("could not add the %s summary columns",xmat_col) ;
   return out ;
}
