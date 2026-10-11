#ifndef _THD_SAMPLESET_HEADER_
#define _THD_SAMPLESET_HEADER_

#include "mrilib.h"
#include "thd_datatable.h"

/*----------------------------------------------------------------------------
  A canonical one-row-per-sample view of a data table.  Input adapters such as
  thd_xmat expand their compact form before this object is built; direct and
  expanded tables therefore use exactly the same dataset, metadata, and row
  bookkeeping.

  Dataset selectors are resolved once.  dset[] holds distinct base datasets,
  while sample_dset[] and sample_brick[] locate every observation.  The table
  is deep-copied, so the caller retains ownership of the input table.
                                                     -- P Molfese, Oct 2026
------------------------------------------------------------------------------*/

typedef struct {
   const char *subject ;  /* NULL: all samples belong to one subject "All"  */
   const char *run ;      /* NULL: no run role                              */
   const char *item ;     /* NULL: no item/trial role                       */
   const char *fold ;     /* NULL: no explicit fold role                    */
   const char *sample ;   /* NULL: synthesize row000001, ...                 */
} THD_sampleset_columns ;

typedef struct {
   int nsample , ndset , nsubject , nrun , nitem , nfold , nvox ;

   THD_3dim_dataset **dset ;  /* [ndset], headers open, common grid          */
   char             **dset_name ;
   int *sample_dset ;         /* [nsample] index into dset                   */
   int *sample_brick ;        /* [nsample] brick of the base dataset        */

   int *subject , *run , *item , *fold ;  /* [nsample], -1 when absent      */
   char **subject_label , **run_label , **item_label , **fold_label ;

   char **sample_label ;      /* [nsample], always present and unique        */
   int   *table_row ;         /* [nsample], zero-based canonical table row   */

   THD_datatable *metadata ;  /* owned deep copy of the canonical table      */
   char *source ;             /* original table source for diagnostics       */
} THD_sampleset ;

/*! Detect conventional role columns that are present: Subj, Run, Trial (or
    Item), Fold, and Sample.  Missing roles are left NULL. */
extern void THD_sampleset_detect_columns( const THD_datatable *dt,
                                          THD_sampleset_columns *columns ) ;

/*! Build and validate a sample set.  Every non-NULL role name must exist.
    InputFile is required and must select exactly one brick per row.  On
    failure, returns NULL and writes a concise reason to err when supplied. */
extern THD_sampleset * THD_sampleset_from_table(
                                      const THD_datatable *dt,
                                      const THD_sampleset_columns *columns,
                                      char *err, int nerr ) ;

/*! Encode any metadata column in first-appearance order.  The returned code
    and labels belong to the caller.  Use THD_sampleset_free_levels(). */
extern int THD_sampleset_encode_column( const THD_sampleset *ss,
                                        const char *column,
                                        int **code, char ***level, int *nlevel,
                                        char *err, int nerr ) ;
extern void THD_sampleset_free_levels( char **level, int nlevel ) ;

/*! A deterministic audit table suitable for a future -show_samples option. */
extern void THD_sampleset_print( const THD_sampleset *ss, FILE *fp ) ;

extern void THD_sampleset_free( THD_sampleset *ss ) ;

#endif /* _THD_SAMPLESET_HEADER_ */
