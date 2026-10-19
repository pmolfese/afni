#ifndef _THD_XMAT_HEADER_
#define _THD_XMAT_HEADER_
/*----------------------------------------------------------------------------
  Design matrices: the X.xmat.1D files that 3dDeconvolve writes (and 3dREMLfit
  reads), and that afni_proc.py keeps for every subject.  Such a file is a
  NIML element in 1D comments: the matrix (one row per UNCENSORED time point,
  one column per regressor) and header attributes that say how it was made:

      ColumnLabels  "Run#1Pol#0 ; ... ; face#0 ; face#1 ; house#0"
      RunStart      "0,466,932"       where each run starts, full timeline
      NRowFull      "1290"            time points before censoring
      GoodList      "3,5..157,..."    the full-timeline index of each row
      RowTR         "2"
      StimLabels    "face ; house"    and StimBots/StimTops: their columns
      BasisOption_000001 "-stim_times_IM" , BasisFormula_000001 "GAM" , ...
      CommandLine   "3dDeconvolve -input ..."

  THD_xmat_read() parses it once.  Every attribute stays available as a
  key/value pair (THD_xmat_attr), and the helpers below answer the questions
  analysis programs ask of it: where runs start, which time points were
  censored, which condition a time point belongs to, and which _Coef volume of
  the stats bucket holds each regressor's beta.

  The THD_datatable_*xmat* functions at the end connect it to '-dataTable'
  tables (thd_datatable.h) with an Xmat column naming each row's file.  They
  live here, not in thd_datatable.c, so that a program that does not use
  design matrices does not have to link this file.
                                                     -- P Molfese, Oct 2026
------------------------------------------------------------------------------*/

#include "thd_datatable.h"

typedef struct {
   char *label ;         /* StimLabels entry                                  */
   int   bot , top ;     /* its columns of the matrix, inclusive              */
   char *option ;        /* BasisOption (-stim_times, -stim_times_IM, ...), or NULL */
   char *formula ;       /* BasisFormula (GAM, TENT(0,12,7), ...), or NULL    */
} THD_xmat_stim ;

typedef struct {
   char   *fname ;
   int     nrow , ncol ; /* rows: the uncensored time points, as in the file  */
   float  *mat ;         /* column c is mat + c*nrow                          */
   char  **collab ;      /* [ncol] ColumnLabels                               */
   int     nfull ;       /* NRowFull (nrow if absent)                         */
   int    *good ;        /* [nrow] GoodList: full-timeline index of each row  */
   int     nrun ;
   int    *runstart ;    /* [nrun] RunStart on the full timeline ({0} if absent) */
   float   tr ;          /* RowTR, or 0                                       */
   int     nstim ;
   THD_xmat_stim *stim ; /* [nstim]                                           */
   int     nattr ;       /* every header attribute, as written                */
   char  **attr_name , **attr_value ;
} THD_xmat ;

/*! One item of a stats bucket: a stimulus column of the matrix and its beta. */
typedef struct {
   char *key ;           /* the stimulus label for a one-column stimulus
                            (house), else the column label (face#3)           */
   char *condition ;     /* the stimulus label                                */
   char *coef_label ;    /* the bucket volume that holds its beta (face#3_Coef) */
   int   stim , column ; /* indexes into stim[] and the matrix's columns      */
   int   run ;           /* the run its regressor lies in, or -1 if it spans
                            runs (or the matrix has none)                     */
} THD_xmat_item ;

/*! Read a design matrix.  NULL on failure, with the reason in err (nerr
    chars; may be NULL). */
extern THD_xmat * THD_xmat_read( const char *fname , char *err , int nerr ) ;
extern void       THD_xmat_free( THD_xmat *xm ) ;

/*! A header attribute by name (case-insensitive), as written; NULL if absent. */
extern const char * THD_xmat_attr( const THD_xmat *xm , const char *name ) ;

/*! The run of full-timeline time point t (0-based), and run r's length. */
extern int THD_xmat_run_of( const THD_xmat *xm , int t ) ;
extern int THD_xmat_run_length( const THD_xmat *xm , int r ) ;

/*! [nfull] 1 where censored (not in GoodList), else 0.  Caller frees. */
extern unsigned char * THD_xmat_censored( const THD_xmat *xm ) ;
extern int THD_xmat_ncensored( const THD_xmat *xm ) ;

/*! Which stimuli: a comma list of labels (NULL: all).  [nstim] 0/1, caller
    frees; *nfound receives how many of the list exist (may be NULL). */
extern unsigned char * THD_xmat_stims_wanted( const THD_xmat *xm , const char *list , int *nfound ) ;

/*! [nfull] the condition of each full-timeline time point: the index of the
    wanted stimulus whose regressor, scaled to its own peak (the largest of its
    columns), is largest there, if that reaches frac; -1 otherwise, and at
    censored time points.  wanted NULL is all.  Caller frees.  The regressors
    are already convolved, so the labels follow the BOLD response. */
extern int * THD_xmat_condition_labels( const THD_xmat *xm , double frac , const unsigned char *wanted ) ;

/*! The items of the wanted stimuli, in column order: one per stimulus
    column.  *nitem receives the count. Invalid/out-of-bounds or overlapping
    stimulus ranges return NULL with count zero, as does allocation failure.
    Free with THD_xmat_items_free. */
extern THD_xmat_item * THD_xmat_items( const THD_xmat *xm , const unsigned char *wanted , int *nitem ) ;
extern void THD_xmat_items_free( THD_xmat_item *it , int nitem ) ;

/*--- tables with an Xmat column (each row names its X.xmat.1D) -------------*/

/*! One row per subject (Subj, Xmat, InputFile = the stats bucket of that
    regression, other columns) becomes one row per item: the subject's other
    columns, then item_col (the key), Condition, run_col when every item of
    every subject lies inside one run, and InputFile as bucket[index] of the
    item's _Coef volume (found by label).  stims as THD_xmat_stims_wanted.
    *have_run (may be NULL) says whether run_col was made.  Errors are fatal
    and name the subject.  The input table is not changed. */
extern THD_datatable * THD_datatable_expand_xmat_items( THD_datatable *dt , const char *xmat_col ,
                                                        const char *item_col , const char *run_col ,
                                                        const char *stims , int *have_run ) ;

/*! Each row's design matrix summarized in new numeric columns, before
    InputFile: Xmat.TR, Xmat.nrun, Xmat.ntime (time points before censoring),
    Xmat.ncensored, Xmat.frac_censored, Xmat.nstim.  Rows sharing a file share
    the values; a row whose Xmat is empty, '-' or NA is an error (the columns
    must be numeric).  The input table is not changed. */
extern THD_datatable * THD_datatable_add_xmat_columns( THD_datatable *dt , const char *xmat_col ) ;

#endif
