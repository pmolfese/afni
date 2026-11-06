#ifndef PICARD_MANIFEST_H
#define PICARD_MANIFEST_H
#include "thd_xmat.h"
#include "thd_picard.h"

/* File-adapter metadata, deliberately separate from numerical plan state. */
typedef struct {
   char *label, *item;
   int role, trial, condition, run, basis;
} PicardColumnIdentity;
typedef struct {
   int id, ncolumn;
   char *matrix, *columns, *source, *version, *sampling, *duration,
       *normalization;
   PicardColumnIdentity *identity;
} PicardCandidateInput;
typedef struct {
   int count;
   char *path;
   PicardCandidateInput *entry;
   size_t metadata_bytes;
} PicardManifest;
/* Versioned, unquoted TSV. Relative matrix/column paths resolve against the
   manifest's directory; spaces are literal. No shell expansion or escaping.
   Header: id matrix columns source version sampling duration normalization.
   The first line is exactly '# 3dPicard candidates 1'. */
PicardManifest *picard_manifest_read(const char *path, int hrf_policy,
                                     char *err, int nerr);
/* Validate one complete identity table against a prepared matrix's roles and
   optional one-column-per-trial mapping. Header: column label role trial
   condition run basis item. Numeric IDs are zero-based; '-' means absent for
   trial/condition/run/item. Role is 'penalized' or 'nuisance'. */
int picard_manifest_read_columns(PicardManifest *m, int index,
                                 const THD_xmat *xm, const int *roles,
                                 const THD_picard_trial *trials, char *err,
                                 int nerr);
/* Semantic comparison includes every column, all row/run mappings, RowTR,
   stimulus labels/ranges and a common duration/normalization convention.
   Matrix values may differ in penalized columns only (numeric plan enforces
   the unchanged nuisance values). Candidate-specific basis formula/source
   metadata may differ, but dimension and basis-component IDs cannot. */
int picard_manifest_compare(const PicardManifest *m, int index,
                            const THD_xmat *reference,
                            const THD_xmat *candidate, char *err, int nerr);
void picard_manifest_free(PicardManifest *m);
#endif
