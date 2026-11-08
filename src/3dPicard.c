/* 3dPicard: the AFNI adapter for fixed or tuned penalized regression.

   This file owns command-line parsing, dataset I/O, identity metadata,
   scheduling and output transactions. thd_picard owns the scientific workflow;
   thd_picard also owns projection; thd_regression owns shared solvers.

   Ownership follows the App object: it owns inputs, plans, worker buffers and
   output buckets. Plans are immutable after construction and can be shared by
   threads. Each Worker owns mutable scratch and a private error message.
   Numerical result pointers are borrowed only until the next fit attempt.

   Numerical shorthand: n = retained time points, p = penalized columns,
   k = nuisance columns, q = voxel responses in the current batch. Numerical
   matrices are row-major; AFNI output bricks each hold one spatial map.

   Suggested reading path:
     1. main() connects preparation, loading, fitting and writing.
     2. prepare_design() translates AFNI matrix metadata into a numeric plan.
     3. load_inputs() and fit_block() map spatial/time indices into batches.
     4. thd_picard.c explains projection, regression and nuisance recovery.
     5. prepare_outputs() and write_outputs() package and save the results.
     6. single_prepare() adds the optional staged single-trial workflow.
     7. The manifest section above App owns candidate-table parsing.

   Staged fitting has two distinct lifetimes. The serial prepass estimates
   noise PCs, freezes each voxel's original HRF choice, and selects one PC
   count. It builds new immutable designs with appended nuisance columns.
   Workers then fit those designs using the stored HRF indices. Intermediate
   OLS maps are optional; scalar choices and provenance always survive until
   output. thd_picard.c section 5 owns the PCA/scoring math, while this adapter
   decides which spatial voxels enter each stage.

   A "plan" is reusable work determined by the design, such as a matrix
   factorization. A "workspace" holds changing values for one batch of voxels.
   Keeping these separate lets threads share expensive design calculations
   while keeping response data and intermediate results private.
*/
#define THD_PICARD_PRIVATE
#include "mrilib.h"
#include <stdarg.h>
#include "thd_xmat.h"
#include "thd_picard.h"
#include "thd_crossval.h"
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <stdint.h>
#include <sys/stat.h>
#ifdef USE_OMP
#include <omp.h>
#endif

#define PROGRAM "3dPicard"
#define NDIAG 7
#define NELASTIC_DIAG 10
#define NTUNE_DIAG 3
#define NPRED_DIAG 5
#define NDESIGN_DIAG 3
/* Dataset slots precede text slots so output allocation and rollback can use
   the same names array. These are internal slots, not sub-brick indices. */
enum {
   OUT_MAIN,
   OUT_FITTED,
   OUT_RESIDUAL,
   DATASET_COUNT,
   OUT_BETAS_TSV = DATASET_COUNT,
   OUT_ROWS_1D,
   OUTPUT_COUNT
};
static const char *diag_labels[NDIAG + NTUNE_DIAG] = {
    "Picard_Fraction",      "Picard_Lambda",           "Picard_PenaltyXtX",
    "Picard_ReferenceNorm", "Picard_AchievedFraction", "Picard_Rank",
    "Picard_Status",        "Picard_CV_MSE",           "Picard_CV_OLS_MSE",
    "Picard_CV_Index"};
static const char *elastic_labels[NELASTIC_DIAG] = {
    "Picard_Lambda",     "Picard_L1Ratio", "Picard_Nonzero",
    "Picard_Iterations", "Picard_KKT",     "Picard_Objective",
    "Picard_LambdaMax",  "Picard_Rank",    "Picard_Standardized",
    "Picard_Status"};
static const char *prediction_labels[NPRED_DIAG] = {
    "Picard_PredMSE", "Picard_PredBaseline", "Picard_PredR2",
    "Picard_PredIndex", "Picard_PredStatus"};
/* Parsed settings. String fields borrow argv; grid is heap-owned. */
typedef struct {
   char *matrix, *input, *mask, *prefix, *pen, *nuis, *fitts, *errts;
   char *trial_info, *fraction_list, *lambda_list, *penalty;
   int single, save_stages, fixed_pcs, pc_stop_set;
   char *noise_pool, *noise_model, *noise_search, *pc_mask, *pc_list;
   double noise_r2, pc_stop;
   int *pc_counts, npc;
   int elastic, rho_set, solver_set;
   int banded, ngroup, nband, support_refit;
   char *group_file, *band_file;
   double *band_values;
   double ar1;
   int legacy_ar1_set, covariance_set, covariance_parameters, covariance_gap;
   THD_picard_covariance covariance;
   char *covariance_map, *covariance_source;
   char *covariance_estimate, *covariance_rows;
   double lambda;
   THD_elastic_config solver;
   double *grid;
   int ngrid, save_cv;
   char *candidate_manifest, *candidate_policy;
   int save_candidates, hrf_policy;
   int predict, block_size, gap, gap_set, project, one_se, global;
   char *fold_file, *cv_nuisance, *cv_rule;
   double fraction, tol, memory;
   int batch, jobs, validate, save_nuis, save_fit, save_res, manifest;
} Options;

/* ========================================================================
   CANDIDATE MANIFEST TYPES AND PARSING
   ======================================================================== */

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
static PicardManifest *picard_manifest_read(const char *path, int hrf_policy,
                                            char *err, int nerr);
/* Validate one complete identity table against a prepared matrix's roles and
   optional one-column-per-trial mapping. Header: column label role trial
   condition run basis item. Numeric IDs are zero-based; '-' means absent for
   trial/condition/run/item. Role is 'penalized' or 'nuisance'. */
static int picard_manifest_read_columns(PicardManifest *m, int index,
                                        const THD_xmat *xm, const int *roles,
                                        const THD_picard_trial *trials,
                                        char *err, int nerr);
/* Semantic comparison includes every column, all row/run mappings, RowTR,
   stimulus labels/ranges and a common duration/normalization convention.
   Matrix values may differ in penalized columns only (numeric plan enforces
   the unchanged nuisance values). Candidate-specific basis formula/source
   metadata may differ, but dimension and basis-component IDs cannot. */
static int picard_manifest_compare(const PicardManifest *m, int index,
                                   const THD_xmat *reference,
                                   const THD_xmat *candidate, char *err,
                                   int nerr);
static void picard_manifest_free(PicardManifest *m);

static int manifest_error(char *err, int nerr, const char *format, ...)
{
   va_list ap;
   if (err && nerr > 0) {
      va_start(ap, format);
      vsnprintf(err, (size_t)nerr, format, ap);
      va_end(ap);
   }
   return 0;
}

/* Return 1 for a valid line, 0 at EOF, -1 for malformed input/I/O failure. */
static int manifest_read_line(FILE *f, char **line, size_t *capacity, char *err,
                              int nerr)
{
   ssize_t n = getline(line, capacity, f);
   if (n < 0)
      return ferror(f)
                 ? (manifest_error(err, nerr, "cannot read candidate metadata"),
                    -1)
                 : 0;
   if (n > 65536 || strlen(*line) != (size_t)n) {
      manifest_error(err, nerr, "metadata line exceeds 65536 bytes or has NUL");
      return -1;
   }
   while (n && ((*line)[n - 1] == '\n' || (*line)[n - 1] == '\r'))
      (*line)[--n] = 0;
   return 1;
}

/* Empty fields and extra tabs are errors, not implicit missing values. */
/* Split exactly the expected number of tab-separated fields in place.
   Spaces are data (including in paths), not delimiters; there is no quoting
   or shell expansion. Returned field pointers borrow the current line buffer.
 */
static int manifest_fields(char *line, char **values, int count)
{
   int i;
   for (i = 0; i < count; i++) {
      char *next = strchr(line, '\t');
      if (next)
         *next++ = 0;
      if (!*line || (i < count - 1 && !next) || (i == count - 1 && next))
         return 0;
      values[i] = line;
      line = next;
   }
   return 1;
}

/* Parse a complete nonnegative identifier, optionally accepting the literal
   missing-value marker. Reject trailing text and overflow so a malformed ID
   cannot silently name a different trial, condition, run or candidate. */
static int manifest_integer(const char *s, int missing, int *value)
{
   char *end;
   long v;
   if (missing && !strcmp(s, "-")) {
      *value = -1;
      return 1;
   }
   if (*s < '0' || *s > '9')
      return 0;
   errno = 0;
   v = strtol(s, &end, 10);
   if (errno || *end || v < 0 || v > INT_MAX)
      return 0;
   *value = (int)v;
   return 1;
}

/* Canonical absolute paths preserve meaning when the output is inspected
   from another working directory. realpath checks existence without writing. */
/* Resolve relative file paths against the manifest directory, not the
   process working directory. The returned path is owned by the manifest;
   absolute paths pass through without shell interpretation. */
static char *manifest_resolve(const char *manifest, const char *field)
{
   const char *slash = strrchr(manifest, '/');
   size_t prefix = slash ? (size_t)(slash - manifest + 1) : 0;
   char *joined, *absolute;
   if (*field == '/')
      return realpath(field, NULL);
   if (prefix > SIZE_MAX - strlen(field) - 1)
      return NULL;
   joined = malloc(prefix + strlen(field) + 1);
   if (!joined)
      return NULL;
   memcpy(joined, manifest, prefix);
   strcpy(joined + prefix, field);
   absolute = realpath(joined, NULL);
   free(joined);
   return absolute;
}

static void picard_manifest_free(PicardManifest *m)
{
   int c, j;
   if (!m)
      return;
   for (c = 0; c < m->count; c++) {
      PicardCandidateInput *e = m->entry + c;
      free(e->matrix);
      free(e->columns);
      free(e->source);
      free(e->version);
      free(e->sampling);
      free(e->duration);
      free(e->normalization);
      if (e->identity)
         for (j = 0; j < e->ncolumn; j++) {
            free(e->identity[j].label);
            free(e->identity[j].item);
         }
      free(e->identity);
   }
   free(m->entry);
   free(m->path);
   free(m);
}

/* Parse a bounded row without mutating any already prepared design. */
/* Validate one candidate record and transfer copies into the growing owned
   table. Count metadata bytes alongside allocations for the adapter budget.
   The common manifest destructor also handles a partially filled last entry. */
static int manifest_add_candidate(PicardManifest *m, char **v, int hrf,
                                  char *err, int nerr)
{
   PicardCandidateInput *array, *e;
   int id, c;
   char *end;
   double sampling;
   if (!manifest_integer(v[0], 0, &id) || id > 16777215)
      return manifest_error(err, nerr, "candidate ID must be in [0,16777215]");
   for (c = 0; c < m->count; c++)
      if (m->entry[c].id == id)
         return manifest_error(err, nerr, "duplicate candidate ID %d", id);
   if (m->count >= 4096)
      return manifest_error(err, nerr,
                            "candidate manifest exceeds 4096 designs");
   sampling = strtod(v[5], &end);
   if (strcmp(v[5], "-") && (*end || !isfinite(sampling) || sampling <= 0))
      return manifest_error(err, nerr,
                            "sampling must be positive seconds or '-'");
   if (hrf && (!strcmp(v[3], "-") || !strcmp(v[4], "-") || !strcmp(v[5], "-") ||
               !strcmp(v[6], "-") || !strcmp(v[7], "-")))
      return manifest_error(err, nerr,
                            "HRF source/version/sampling/duration/"
                            "normalization must be explicit");
   array = realloc(m->entry, (size_t)(m->count + 1) * sizeof(*array));
   if (!array)
      return manifest_error(err, nerr, "cannot allocate candidate row");
   m->entry = array;
   e = m->entry + m->count++;
   memset(e, 0, sizeof(*e));
   e->id = id;
   e->matrix = manifest_resolve(m->path, v[1]);
   e->columns = manifest_resolve(m->path, v[2]);
   e->source = strdup(v[3]);
   e->version = strdup(v[4]);
   e->sampling = strdup(v[5]);
   e->duration = strdup(v[6]);
   e->normalization = strdup(v[7]);
   if (!e->matrix || !e->columns || !e->source || !e->version || !e->sampling ||
       !e->duration || !e->normalization)
      return manifest_error(err, nerr,
                            "candidate %d: missing file or allocation "
                            "failure resolving matrix/columns",
                            id);
   m->metadata_bytes += sizeof(*e) + strlen(e->matrix) + strlen(e->columns) +
                        strlen(e->source) + strlen(e->version) +
                        strlen(e->sampling) + strlen(e->duration) +
                        strlen(e->normalization) + 7;
   return 1;
}

static PicardManifest *picard_manifest_read(const char *path, int hrf,
                                            char *err, int nerr)
{
   PicardManifest *m = calloc(1, sizeof(*m));
   FILE *f;
   char *line = NULL, *v[8];
   size_t capacity = 0;
   int status, ok, header = 0;
   if (!m) {
      manifest_error(err, nerr, "cannot allocate manifest");
      return NULL;
   }
   m->path = path ? realpath(path, NULL) : NULL;
   f = m->path ? fopen(m->path, "r") : NULL;
   if (!f) {
      manifest_error(err, nerr, "cannot open candidate manifest %s",
                     path ? path : "(null)");
      picard_manifest_free(m);
      return NULL;
   }
   ok = manifest_read_line(f, &line, &capacity, err, nerr) == 1 &&
        !strcmp(line, "# 3dPicard candidates 1");
   if (!ok)
      manifest_error(err, nerr, "expected '# 3dPicard candidates 1'");
   while (ok &&
          (status = manifest_read_line(f, &line, &capacity, err, nerr)) != 0) {
      if (status < 0) {
         ok = 0;
         break;
      }
      if (!*line || *line == '#')
         continue;
      if (!header) {
         ok = !strcmp(line, "id\tmatrix\tcolumns\tsource\tversion\tsampling\t"
                            "duration\tnormalization");
         if (!ok)
            manifest_error(err, nerr, "invalid candidate TSV header");
         header = 1;
      } else {
         ok = manifest_fields(line, v, 8);
         if (!ok)
            manifest_error(err, nerr, "candidate row needs 8 TSV fields");
         if (ok)
            ok = manifest_add_candidate(m, v, hrf, err, nerr);
      }
   }
   if (ok && m->count < 2)
      ok = manifest_error(err, nerr, "need >=2 candidates");
   free(line);
   fclose(f);
   if (!ok) {
      picard_manifest_free(m);
      return NULL;
   }
   return m;
}

/* Parse into column-indexed storage, so file row order carries no identity. */
/* Validate a single identity-table row against its matrix column and role.
   Trial metadata must agree with the supplied trial table when present; rows
   may arrive in any order, but duplicate column assignments are rejected. */
static int manifest_column_row(PicardCandidateInput *e, const THD_xmat *xm,
                               const int *roles, char **v, char *err, int nerr)
{
   PicardColumnIdentity *ci;
   int column, role;
   if (!manifest_integer(v[0], 0, &column) || column >= xm->ncol)
      return manifest_error(err, nerr, "column identity index out of range");
   role = !strcmp(v[2], "penalized") ? 1 : !strcmp(v[2], "nuisance") ? 0 : -1;
   ci = e->identity + column;
   if (ci->label || strcmp(v[1], xm->collab[column]) || role != roles[column])
      return manifest_error(err, nerr,
                            "duplicate column or label/role mismatch "
                            "in candidate %d column %d",
                            e->id, column);
   if (!manifest_integer(v[3], 1, &ci->trial) ||
       !manifest_integer(v[4], 1, &ci->condition) ||
       !manifest_integer(v[5], 1, &ci->run) ||
       !manifest_integer(v[6], 0, &ci->basis) ||
       (ci->trial >= 0 && (ci->condition < 0 || ci->run < 0 || !role)) ||
       (!role && ci->condition >= 0) || ci->run >= xm->nrun)
      return manifest_error(err, nerr,
                            "invalid trial/condition/run/basis identity");
   ci->role = role;
   ci->label = strdup(v[1]);
   ci->item = strdup(v[7]);
   return (ci->label && ci->item) ||
          manifest_error(err, nerr, "cannot allocate column identity");
}

static int picard_manifest_read_columns(PicardManifest *m, int index,
                                        const THD_xmat *xm, const int *roles,
                                        const THD_picard_trial *trials,
                                        char *err, int nerr)
{
   PicardCandidateInput *e;
   char *line = NULL, *v[8];
   size_t capacity = 0;
   FILE *f;
   int status, ok = 1, header = 0, i, t = 0;
   if (!m || !xm || !roles || index < 0 || index >= m->count)
      return manifest_error(err, nerr, "invalid candidate column inputs");
   e = m->entry + index;
   if (e->identity)
      return manifest_error(err, nerr, "candidate identities already loaded");
   f = fopen(e->columns, "r");
   if (!f)
      return manifest_error(err, nerr, "cannot open columns %s", e->columns);
   e->ncolumn = xm->ncol;
   e->identity = calloc((size_t)xm->ncol, sizeof(*e->identity));
   if (!e->identity)
      ok = manifest_error(err, nerr, "cannot allocate identities");
   while (ok &&
          (status = manifest_read_line(f, &line, &capacity, err, nerr)) != 0) {
      if (status < 0) {
         ok = 0;
         break;
      }
      if (!*line || *line == '#')
         continue;
      if (!header) {
         ok = !strcmp(
             line, "column\tlabel\trole\ttrial\tcondition\trun\tbasis\titem");
         if (!ok)
            manifest_error(err, nerr, "invalid columns TSV header");
         header = 1;
      } else {
         ok = manifest_fields(line, v, 8);
         if (!ok)
            manifest_error(err, nerr, "column row needs 8 TSV fields");
         if (ok)
            ok = manifest_column_row(e, xm, roles, v, err, nerr);
      }
   }
   for (i = 0; i < xm->ncol && ok; i++) {
      const PicardColumnIdentity *ci = e->identity + i;
      if (!ci->label)
         ok = manifest_error(err, nerr, "missing column identity %d", i);
      else if (trials && roles[i]) {
         if (ci->trial != trials[t].trial ||
             ci->condition != trials[t].condition || ci->run != trials[t].run ||
             ci->basis != 0)
            ok = manifest_error(err, nerr,
                                "column identity disagrees with trial_info");
         t++;
      }
      if (ok && ci->run >= 0) {
         int row;
         for (row = 0; row < xm->nrow && ok; row++)
            if (THD_xmat_run_of(xm, xm->good[row]) != ci->run &&
                xm->mat[(size_t)i * xm->nrow + row] != 0)
               ok = manifest_error(
                   err, nerr,
                   "declared column run disagrees with matrix support");
      }
      if (ok && ci->trial >= 0) {
         int previous;
         for (previous = 0; previous < i && ok; previous++) {
            const PicardColumnIdentity *other = e->identity + previous;
            if (other->trial == ci->trial &&
                (other->basis == ci->basis ||
                 other->condition != ci->condition || other->run != ci->run ||
                 strcmp(other->item, ci->item)))
               ok = manifest_error(
                   err, nerr, "duplicate or inconsistent trial/basis identity");
         }
      }
      if (ok)
         m->metadata_bytes +=
             sizeof(*ci) + strlen(ci->label) + strlen(ci->item) + 2;
   }
   free(line);
   fclose(f);
   return ok;
}

static int manifest_same_text(const char *a, const char *b)
{
   return (!a && !b) || (a && b && !strcmp(a, b));
}

static int picard_manifest_compare(const PicardManifest *m, int index,
                                   const THD_xmat *a, const THD_xmat *b,
                                   char *err, int nerr)
{
   const PicardCandidateInput *base = m->entry, *entry = m->entry + index;
   const char *ta = THD_xmat_attr(a, "RowTR"), *tb = THD_xmat_attr(b, "RowTR");
   int j;
   if (a->nrow != b->nrow || a->nfull != b->nfull || a->ncol != b->ncol ||
       a->nrun != b->nrun || a->tr != b->tr || a->nstim != b->nstim ||
       ((!ta != !tb) || (ta && tb && strtod(ta, NULL) != strtod(tb, NULL))) ||
       memcmp(a->good, b->good, (size_t)a->nrow * sizeof(int)) ||
       memcmp(a->runstart, b->runstart, (size_t)a->nrun * sizeof(int)))
      return manifest_error(err, nerr,
                            "candidate row/run/TR/dimension metadata differs");
   if (strcmp(base->duration, entry->duration) ||
       strcmp(base->normalization, entry->normalization))
      return manifest_error(err, nerr,
                            "candidate duration/normalization differs");
   for (j = 0; j < a->nstim; j++)
      if (!manifest_same_text(a->stim[j].label, b->stim[j].label) ||
          a->stim[j].bot != b->stim[j].bot || a->stim[j].top != b->stim[j].top)
         return manifest_error(err, nerr,
                               "candidate stimulus identities differ");
   for (j = 0; j < a->ncol; j++) {
      const PicardColumnIdentity *x = base->identity + j,
                                 *y = entry->identity + j;
      if (strcmp(x->label, y->label) || strcmp(x->item, y->item) ||
          x->role != y->role || x->trial != y->trial ||
          x->condition != y->condition || x->run != y->run ||
          x->basis != y->basis)
         return manifest_error(err, nerr,
                               "candidate coefficient identity differs "
                               "at column %d",
                               j);
   }
   return 1;
}

/* ========================================================================
   APPLICATION STATE AND COMMAND-LINE INTERFACE
   ======================================================================== */

/* One worker per OpenMP thread, with no mutable state shared by fits.
   One fitting-mode workspace plus optional stage-C OLS scratch is allocated. */
typedef struct {
   THD_picard_candidates_workspace *candidates;
   THD_picard_workspace *work;
   THD_picard_repeat_workspace *repeat;
   THD_picard_predict_workspace *predict;
   double *y, *fractions;
   int *locked_hrf;
   double *locked_loss;
   THD_picard_workspace **stage_c;
   int nstage;
   int failed;
   char err[512];
} Worker;
/* Data-derived stages are prepared before allocating the final buckets.
   Scalar maps are spatially indexed; beta_b is allocated only when requested.
   Final plans own appended PC columns while the original plans remain intact
   for trial labels, the pre-PC HRF comparison and provenance. */
typedef struct {
   int nh, chosen, nrun, added, nscore;
   unsigned char *pool, *select;
   int *hrf;
   double *hrf_loss, *design_loss, *noise_r2, *beta_b, *trend;
   THD_picard_noise *noise;
   THD_picard_plan **final;
   THD_picard_candidates *collection;
} Single;

/* Zero-initialized owner: free_app also handles partial setup failures.
   "Owns" means this object must eventually free the allocation. "Borrows"
   means it may use another object's allocation but must not free it. Most
   preparation functions attach allocations here as soon as they succeed;
   main can then use one cleanup call regardless of which stage failed. */
typedef struct {
   Options opt;
   Single *single;
   THD_xmat *xm;
   THD_picard_plan *plan;
   PicardManifest *manifest;
   THD_picard_plan **candidate_designs;
   THD_picard_candidates *candidate_plan;
   THD_picard_repeat_plan *repeat;
   THD_picard_predict_plan *predict;
   THD_foldplan *folds;
   int global_index;
   THD_picard_trial *trials;
   THD_3dim_dataset *input, *mask, *out[DATASET_COUNT];
   char *names[OUTPUT_COUNT];
   int *groups; /* Penalized-column group IDs for banded fits. */
   int *roles;  /* One penalized/nuisance flag per original matrix column. */
   int *runs;   /* One original run ID per retained matrix row. */
   int *vox;    /* Compact selected-voxel index -> full spatial index. */
   int nvox, nselected;
   double *matrix; /* Temporary row-major copy, released after plan creation. */
   Worker *workers;
} App;

static void free_app(App *a);
static int single_prepare(App *a, char *err, int nerr);
static int single_output(App *a, char *err, int nerr);
static int single_start(const App *a);

/* Write a caller-owned error buffer and return false for concise checks. */
static int failure(char *err, int nerr, const char *msg)
{
   if (err && nerr > 0)
      snprintf(err, (size_t)nerr, "%s", msg);
   return 0;
}

/* Keep both source literals and rendered help within 80 columns. */
// clang-format off
static void help(void)
{
   puts(
      "3dPicard -- voxelwise penalized regression\n"
      "=========================================\n"
      "Fit or tune coefficients with unpenalized nuisance columns.\n"
      "\n"

      "Usage:\n"
      "  3dPicard -input DATA -matrix X.xmat.1D -fraction F [options]\n"
      "\n"

      "STAGED SINGLE-TRIAL FITTING -- ENGAGE\n"
      "-----------------------------------\n"
      "  -single_trial         Enable supplied-HRF -> noise PCs -> penalty.\n"
      "  -noise_pool MASK      Explicit spatial pool for noise estimation.\n"
      "  -noise_model X        Alternative: supplied coarse task xmat.1D.\n"
      "  -noise_search MASK    Search mask required with noise_model.\n"
      "  -noise_r2 R           Keep search voxels with task R2 < R.\n"
      "                        R must be in [0,1].\n"
      "  -noise_pcs K          Fixed nonnegative PC count per run; OR:\n"
      "  -pc_counts LIST       Increasing integer counts starting at 0.\n"
      "                        Example: 0,1,2,3. Requires cross-run repeats.\n"
      "  -pc_mask MASK         Explicit PC-selection voxels, inside analysis\n"
      "                        mask and outside the selected noise pool.\n"
      "  -pc_stop S            First count reaching best improvement/S.\n"
      "                        S >= 1; default 1.05. No improvement: 0 PCs.\n"
      "  -save_stages          Add StageB_ and StageC_ beta sub-bricks to the\n"
      "                        main bucket. Ordinary betas are final Stage D.\n"
      "\n"
      "Use either noise_pool OR noise_model/noise_search/noise_r2, and either\n"
      "noise_pcs OR pc_counts/pc_mask. All masks use the input volume grid.\n"
      "The pool may extend outside the analysis mask. Generic prediction CV\n"
      "is a separate mode. trial_info is required even with a fixed penalty.\n"
      "Every trial and supplied nuisance column must be local to one run.\n"
      "A fixed count/fixed penalty does not require repeated conditions.\n"
      "\n"
      "A candidate_manifest with candidate_policy hrf selects each voxel's\n"
      "HRF by unregularized full-data fit BEFORE adding PCs. This HRF choice\n"
      "is then frozen. A single matrix uses its supplied HRF throughout.\n"
      "The coarse noise model must share rows, censoring, runs and TR with\n"
      "the trial design and have no more columns. Its task R2 is improvement\n"
      "relative to its supplied nuisance-only fit. Undefined R2 is excluded.\n"
      "\n"
      "PCA removes a constant and the supplied run-local nuisance span from\n"
      "pool time series, normalizes each voxel to unit L2 norm, and extracts\n"
      "run-specific temporal PCs. Each PC has unit sample SD. No task fit,\n"
      "whitening or automatic drift polynomials are used to estimate PCs.\n"
      "Requested PCs exceeding any run's numerical noise rank cause an error.\n"
      "Selected PCs enter as unpenalized, run-local nuisance columns; final\n"
      "fitting preserves the chosen estimator and any fixed AR(1) policy.\n"
      "\n"
      "PC-count scoring compares cross-run, same-condition OLS trial betas\n"
      "with a FIXED no-PC OLS reference. Pair MSE is divided by the reference\n"
      "sample beta variance, then the median is taken over pc_mask voxels.\n"
      "Zero-variance voxels are excluded; all undefined is an error. This is\n"
      "one-session scoring. Later penalty tuning uses the existing repeat\n"
      "policy with the selected PCs fixed. These are CONDITIONAL SELECTION\n"
      "scores, not independent test performance: HRF and pool choices used\n"
      "these data, and each run contributed its own noise-pool time series.\n"
      "\n"
      "Stage B = chosen HRF, no PCs, unregularized fit.\n"
      "Stage C = chosen HRF, selected PCs, unregularized fit.\n"
      "Stage D = chosen HRF, selected PCs, final fixed or tuned penalty.\n"
      "The main bucket adds Picard_PCCount, Picard_NoisePool,\n"
      "Picard_PCSelection and Picard_NoiseModelR2. Pool/selection flags cover\n"
      "the whole grid; count is -1 outside analysis. Undefined/uncomputed R2\n"
      "is -1. Header attributes record PC values, run ranks, the count-loss\n"
      "trend, masks, policies and augmented nuisance-column identities.\n"
      "validate_only checks supplied designs/options, not learned PC rank,\n"
      "mask contents or the coarse noise model; those require the data pass.\n"
      "\n"
      "This follows GLMsingle's staged idea, not its entire automatic policy.\n"
      "Supply preconvolved HRFs, trial identities, drift columns and masks.\n"
      "There is no automatic intensity pool, session normalization or beta\n"
      "autoscaling. PCA uses double precision and explicit rank checks.\n"
      "The final fractional penalty uses Picard's fold-local repeat policy.\n"
      "Consult the phase-7 reference record for precisely tested equivalence.\n"
      "\n"
      "Example (run-local designs and repeated conditions already supplied):\n"
      "  3dPicard -input runs+orig -matrix trials.xmat.1D                 \\\n"
      "    -single_trial -trial_info trials.1D -noise_pool noise+orig     \\\n"
      "    -pc_counts 0,1,2,3 -pc_mask reliable+orig -fractions 1,.8,.5   \\\n"
      "    -save_stages -prefix picard\n"
      "\n"
      "PENALTIES -- THERE ARE FOUR LIGHTS\n"
      "--------------------------------\n"
      "  -penalty ridge|lasso|enet|banded  Default: ridge (fractional\n"
      "                            parameterization).\n"
      "  -lambda L                 Fixed positive lambda for lasso/enet.\n"
      "  -lambdas LIST             Positive, strictly decreasing tuning grid.\n"
      "                            Example: '1,.3,.1,.03'. Use -cv_predict or\n"
      "                            -trial_info FILE, as with -fractions.\n"
      "  -l1_ratio R               Elastic-net mixing in [0,1]; default: 0.5.\n"
      "                            Lasso fixes R=1. Enet with R=0 is ordinary\n"
      "                            ridge parameterized by lambda, not\n"
      "                            fraction.\n"
      "  -standardize              Scale projected X columns to unit RMS.\n"
      "                            Lasso/enet only; default keeps input\n"
      "                            units.\n"
      "  -solver_tol T             Absolute maximum KKT error; default: 1e-8.\n"
      "  -max_iterations N         Maximum coordinate sweeps; default: 10000.\n"
      "\n"
      "For lasso/enet the objective is:\n"
      "  SSE/(2*n) + lambda * (R*sum(abs(theta)) + (1-R)*sum(theta^2)/2).\n"
      "R=1 can select exact zeros; 0<R<1 combines sparsity and shrinkage.\n"
      "Lambda is sklearn's alpha, not its l1_ratio. Lambda must be positive;\n"
      "use ridge -fraction 1 for the unregularized minimum-norm reference.\n"
      "\n"
      "Project X and Y against supplied nuisance Z first. Optional scaling\n"
      "uses sqrt(sum(projected_X^2)/n), independently within each training\n"
      "fold. Zero columns keep scale 1. There is no implicit centering,\n"
      "intercept, or response scaling. theta is in this optimization space;\n"
      "output beta is back-transformed to supplied column units. Recover\n"
      "unpenalized nuisance coefficients from original Y-X*beta afterward.\n"
      "\n"
      "Paths warm-start only within each fold/response batch. Every sweep\n"
      "checks ALL KKT conditions. Nonconvergence fails before output.\n"
      "Final selected fits start cold on the complete design. Changing R is\n"
      "an explicit model choice; this release tunes lambda at one fixed R.\n"
      "\n"
      "Both tuning modes below accept lambdas instead of fractions. Exact\n"
      "minimum-loss ties choose the smallest lambda. The prediction 1se rule\n"
      "chooses the LARGEST eligible lambda. All-undefined global scores\n"
      "choose the smallest lambda. OLS targets stay unregularized.\n"
      "No independence claim is made for tuning scores. More regularization\n"
      "does not guarantee a more scientifically meaningful sparse model.\n"
      "\n"
      "Sparse fits replace the seven ridge diagnostic maps with ten maps:\n"
      "  Picard_Lambda, Picard_L1Ratio, Picard_Nonzero, Picard_Iterations,\n"
      "  Picard_KKT, Picard_Objective, Picard_LambdaMax, Picard_Rank,\n"
      "  Picard_Standardized, Picard_Status.\n"
      "KKT/objective use the optimization units above. LambdaMax is the\n"
      "all-zero threshold at R>0, otherwise -1. Rank describes projected X,\n"
      "not selected support. Status=0 for converged voxels, -1 outside mask.\n"
      "CV maps follow these diagnostics in the SAME bucket; -save_cv appends\n"
      "one loss map per lambda. Headers record PICARD_LAMBDAS and policies.\n"
      "The fraction-specific explanations below apply to -penalty ridge.\n"
      "\n"

      "BANDED RIDGE -- DIFFERENT DEPARTMENTS, DIFFERENT ORDERS\n"
      "----------------------------------------------------\n"
      "  -penalty banded       Separate ridge penalties for predictor groups.\n"
      "  -band_groups FILE    One integer per ORIGINAL matrix column, in one\n"
      "                       row or column. Nuisance=-1; predictors=0..G-1.\n"
      "                       Every group must occur in the full design.\n"
      "  -band_penalties FILE Whitespace table: rows=candidate vectors,\n"
      "                       columns=groups 0..G-1. Positive finite lambdas.\n"
      "                       # comments allowed; duplicate vectors rejected.\n"
      "One row gives a fixed fit. Multiple rows require -cv_predict or\n"
      "-trial_info. Do not combine with fraction(s), lambda(s), standardize,\n"
      "l1_ratio or solver options. Group identities are explicit, independent\n"
      "of stimulus labels. Candidate indices follow table order, starting at\n"
      "0.\n"
      "\n"
      "Minimize SSE/(2*n) + sum_j lambda[group[j]]*beta[j]^2/2 after removing\n"
      "nuisance. Each training fold uses n_train*lambda in its normal\n"
      "equations.\n"
      "CV ties prefer smaller MEAN group lambda, then earlier table row.\n"
      "The 1se rule prefers larger mean group lambda among eligible vectors.\n"
      "This is an explicit ordering convention for otherwise unordered\n"
      "vectors.\n"
      "\n"
      "Banded diagnostic maps: Picard_BandIndex, Picard_Objective,\n"
      "Picard_KKT,\n"
      "Picard_Rank, Picard_LambdaGroup000... and Picard_Status (0=success).\n"
      "CV/design maps follow in the same bucket. Headers record the group map\n"
      "and complete penalty table. No extra penalty-map datasets are needed.\n"
      "\n"
      "WHITENING AND SUPPORT REFITS -- CHECK THE INERTIAL DAMPERS\n"
      "-------------------------------------------------------\n"
      "  -covariance white|ar1|arma11  Fixed unit-marginal covariance.\n"
      "  -covariance_ar A      AR parameter, default 0, strictly (-1,1).\n"
      "  -covariance_ma B      MA parameter (arma11 only), default 0.\n"
      "  -covariance_gaps preserve|reset  Default preserve: actual TR lags.\n"
      "  -covariance_estimate XMAT  Estimate native pilot REML per voxel.\n"
      "  -covariance_rows FILE  Estimation rows: original IDs, one per line.\n"
      "Default estimation uses all retained rows. Pilot must share the full\n"
      "timeline/censor/run map. AR/MA grid is [-.9,.9] by .1; AR1 fixes B=0.\n"
      "Parameters are shared across runs per voxel; no spatial pooling.\n"
      "Native mode requires fixed penalty, not internal CV or staged fits.\n"
      "For tuning, estimate separately on each training fold using rows,\n"
      "freeze parameters, then fit validation data with the resulting maps.\n"
      "Saturated pilots or zero residual variance fail explicitly. Outputs\n"
      "include marginal Sigma2, REML score, pilot rank, status and AR/MA.\n"
      "Status: 0 interior, 1 boundary, 3 flat, 4 tied; outside mask -1.\n"
      "  -covariance_map DSET  AR1: one a brick; ARMA: a,b first two bricks.\n"
      "                       Accepts REMLfit -Rvar; other bricks ignored.\n"
      "  -covariance_source TEXT  Optional external estimation provenance.\n"
      "Without covariance_estimate, supply parameters. MA uses +B*u[t-1].\n"
      "Runs are independent; censored samples stay excluded, never filled.\n"
      "Exact dense Cholesky factors are memory-budgeted and shared by equal\n"
      "parameter groups. No parameter rounding. Map mode excludes staged\n"
      "single_trial and global tuning (future integration); scalar mode\n"
      "supports existing workflows. Map estimates remain externally supplied\n"
      "in every fold: scores are conditional, not independent noise tuning.\n"
      "Invalid loaded parameters/geometry fail; no white fallback. AFNI\n"
      "loaders may repair NaN/Inf; repaired values cannot be detected.\n"
      "New covariance options cannot be combined with prewhite_ar1.\n"
      "New AR1 penalty units differ from the legacy filter: identical fits\n"
      "with matching gaps require lambda_new=lambda_old/(1-A*A).\n"
      "Covariance map validation needs input/mask data even in validate_only.\n"
      "\n"
      "  -prewhite_ar1 R       Externally fixed AR(1) coefficient in (-1,1).\n"
      "                       Default 0: identity. No noise estimation.\n"
      "Whiten X, nuisance Z and Y BEFORE nuisance projection. Within a\n"
      "contiguous segment use v[t]-R*v[t-1]; at its start use\n"
      "sqrt(1-R*R)*v[t].\n"
      "Restart at run changes and censor gaps, separately within every fold.\n"
      "Training never reads a held-out response as its previous sample.\n"
      "Fitts/errts and selection scores stay in original data/coefficient\n"
      "units;\n"
      "conditional prediction still projects the held-out nuisance span.\n"
      "HRF selection uses unpenalized fits with the chosen whitening policy,\n"
      "then compares original-unit residual MSE. It remains selection-only.\n"
      "This is fixed AR(1), not voxel-specific ARMA estimation or REML\n"
      "fitting.\n"
      "\n"
      "  -support_refit       Lasso/enet with l1_ratio>0 only. Fit minimum-\n"
      "                       norm\n"
      "                       OLS using the exact nonzero selected\n"
      "                       predictors,\n"
      "                       with the same whitening/nuisance\n"
      "                       transformation.\n"
      "Every training fold selects its own support; final selection/refitting\n"
      "uses all retained rows. Empty support gives nuisance-only predictions.\n"
      "Output betas/fitts/errts are REFITTED. Sparse Objective/KKT/Nonzero\n"
      "maps\n"
      "still describe the penalized selection fit. Picard_SupportSize and\n"
      "Picard_SupportRank appear before Status. This adds no t/F inference:\n"
      "choosing a support from these data does not justify ordinary OLS\n"
      "tests.\n"
      "Support SVD scratch is allocated for one voxel at a time; other\n"
      "fitting\n"
      "scratch is preallocated. Rank-deficient supports use minimum norm.\n"
      "\n"

      "DEDICATION -- THE INNER LIGHT\n"
      "----------------------------\n"
      "Dedicated to Captain Jean-Luc Picard's curiosity and humanity,\n"
      "and to the memory of Sir Patrick Stewart, who brought him to life\n"
      "and who passed away on October 28, 2026.\n"
      "\n"

      "In 'The Inner Light', a lifetime seen through another's eyes leaves\n"
      "Picard with a new perspective and the music of a Ressikan flute.\n"
      "That willingness to see differently inspires this program: explore\n"
      "what regularization can reveal about an uncertain fit, while keeping\n"
      "the assumptions, the evidence, and the limits in view.\n"
      "\n"

      "  'The first duty of every Starfleet officer is to the truth'\n"
      "                                      -- Picard, 'The First Duty'\n"
      "\n"

      "INPUT AND DESIGN -- CHART THE COURSE\n"
      "-----------------------------------\n"
      "  -input DATA          Volume time series. Required for fitting.\n"
      "  -matrix X.xmat.1D    AFNI design matrix. Required.\n"
      "  -fraction F          Fixed coefficient-norm fraction in [0,1].\n"
      "                       Use this OR a -fractions tuning mode.\n"
      "  -mask DATA           Optional matching one-volume mask.\n"
      "                       Nonzero mask voxels are kept.\n"
      "\n"

      "  -penalty_cols LIST   Penalized columns: zero-based AFNI integer list\n"
      "                       or 'all'.\n"
      "  -nuisance_cols LIST  Nuisance columns: zero-based AFNI integer list\n"
      "                       or 'none'.\n"
      "\n"

      "Supply one column list and its complement gets the other role.\n"
      "Supply both lists and they must be disjoint and cover all columns.\n"
      "Otherwise, ColumnGroups > 0 are penalized; groups <= 0 are nuisance.\n"
      "Missing or malformed ColumnGroups require explicit column roles.\n"
      "\n"

      "GoodList, NRowFull, and RunStart refer to original acquisition time.\n"
      "The input volume count must match NRowFull. Censored matrices require\n"
      "explicit GoodList and NRowFull; there is no extra censoring option.\n"
      "Input and mask must be unselected datasets (no bracket selectors).\n"
      "Volume data only. All formats use the standard AFNI dataset loader.\n"
      "This includes its datatype conversion, scaling and NaN/Inf handling.\n"
      "Picard honors AFNI reader settings (including AFNI_FLOATSCAN).\n"
      "Fitting uses doubles after loading; it cannot restore precision or\n"
      "nonfinite values already changed by the loader. Values still\n"
      "nonfinite in the mask or retained samples are rejected.\n"
      "\n"

      "OUTPUT -- THE CAPTAIN'S LOG\n"
      "--------------------------\n"
      "  -prefix P            Main bucket prefix; default: Picard (HEAD/BRIK)\n"
      "                       Use P.nii or P.nii.gz to request NIfTI.\n"
      "  -save_nuisance       Include nuisance coefficients in that bucket.\n"
      "  -fitts PREFIX        Write a fitted retained-row dataset.\n"
      "  -errts PREFIX        Write a residual retained-row dataset.\n"
      "  -save_fitted         Use P_fitts in the main output format.\n"
      "  -save_residuals      Use P_errts in the main output format.\n"
      "  -write_manifest      Write P.betas.tsv and P.rows.1D sidecars.\n"
      "\n"

      "By default, write one AFNI bucket: P+orig.HEAD/BRIK (or input view).\n"
      "Sub-bricks contain penalized coefficients, optional nuisance\n"
      "coefficients, then Picard_* diagnostics listed below.\n"
      "Coefficient labels use the original column labels plus _Coef.\n"
      "\n"

      "The log travels with the ship: PICARD_* header attributes and history\n"
      "record column roles, matrix source, solver settings, and row/run maps.\n"
      "Text sidecars are optional. Existing outputs are never overwritten.\n"
      "\n"

      "Fitted and residual sub-bricks contain retained rows only, labelled\n"
      "by original time. These are buckets, not a regular-TR full timeline.\n"
      "PICARD_ROW_IDS maps their rows; censored predictions are not invented.\n"
      "\n"

      "MODEL AND DIAGNOSTICS -- THERE ARE FOUR LIGHTS!\n"
      "---------------------------------------------\n"
      "The supplied design is the design: no implicit intercept, detrending,\n"
      "whitening, or scaling. Include baseline/intercept columns\n"
      "in the matrix with nuisance roles.\n"
      "\n"

      "The fraction compares the penalized coefficient norm to a minimum-\n"
      "norm least-squares reference after nuisance projection, in supplied\n"
      "coefficient units. Fraction 1 gives that unregularized reference;\n"
      "fraction 0 gives zero penalized coefficients.\n"
      "\n"

      "Rank retains singular values s > rank_tol*smax. Projection removes\n"
      "numerical cancellation at 64*DBL_EPSILON*max(n,k) relative to the\n"
      "supplied vector norm. Redundant nuisance columns receive minimum-norm\n"
      "coefficients; their fitted signal remains identifiable.\n"
      "\n"

      "Diagnostic sub-brick labels:\n"
      "  Picard_Fraction          Fixed or selected norm fraction.\n"
      "  Picard_Lambda            Normalized ridge penalty, PenaltyXtX/n.\n"
      "  Picard_PenaltyXtX        Penalty added to projected X'X.\n"
      "  Picard_ReferenceNorm     Minimum-norm reference coefficient norm.\n"
      "  Picard_AchievedFraction  Achieved fraction of the reference norm.\n"
      "  Picard_Rank              Rank of the projected penalized design.\n"
      "  Picard_Status            Fit status code (see below).\n"
      "\n"

      "Lambda uses the objective SSE/(2n) + lambda*||beta||^2/2, where n is\n"
      "the number of retained rows. A value of -1 denotes an infinite\n"
      "penalty (fraction 0) or an undefined achieved fraction (zero\n"
      "reference norm). A zero reference yields zero betas and lambda = 0.\n"
      "\n"

      "Status codes:  0 = normal fit       1 = zero reference\n"
      "               2 = zero endpoint  -1 = outside mask\n"
      "\n"

      "No conventional penalized t- or p-value maps are produced.\n"
      "\n"

      "REPEATED-CONDITION TUNING -- COMPARE THE CAPTAIN'S LOGS\n"
      "----------------------------------------------------\n"
      "  -fractions LIST      Comma-separated fractions (>=2 unique).\n"
      "                       Example: '1,.8,.6,.4,.2,0'. Replaces -fraction.\n"
      "  -trial_info FILE     Repeat-tuning identities. Four integer columns:\n"
      "                         column  trial  condition  run\n"
      "                       Zero-based column and RunStart run index.\n"
      "                       Nonnegative IDs. Optional # comments.\n"
      "  -save_cv             Append candidate MSE maps to the main bucket:\n"
      "                       Picard_CV_Grid000, Grid001, etc.\n"
      "\n"

      "Describe every penalized column exactly once. Trial IDs must be\n"
      "unique; matching condition IDs identify repeats across runs. One\n"
      "column per trial is supported; multi-basis trials require a later\n"
      "extension. File row order is immaterial. Each trial column must have\n"
      "nonzero support confined to its specified run. Nuisance columns must\n"
      "also be run-local: split a shared intercept or nuisance vector into\n"
      "separate run columns before tuning.\n"
      "\n"

      "Each retained run must contain a condition repeated in another run.\n"
      "Unmatched conditions contribute no score but remain in the fitted\n"
      "model. Without cross-run repeats, supply a fixed -fraction instead.\n"
      "\n"

      "Leave one run out. Fit the other runs jointly at each fraction,\n"
      "including only their trial/nuisance columns. Compare EACH\n"
      "training-trial beta with each held-out unregularized beta for the same\n"
      "condition. Pool squared errors across ordered pairs and divide by the\n"
      "total number of pairs. Unequal folds therefore receive pair-count\n"
      "weights, not equal run weights.\n"
      "\n"

      "Choose the smallest MSE separately for each voxel. Exact ties choose\n"
      "the largest fraction. Refit the full design at that fraction for final\n"
      "outputs. The fraction reference and lambda belong to each fit's own\n"
      "design; the final full-design lambda need not equal any training-fold\n"
      "lambda.\n"
      "\n"

      "Additional diagnostic sub-bricks (always included when tuning):\n"
      "  Picard_CV_MSE         Selected candidate's pooled pair MSE.\n"
      "  Picard_CV_OLS_MSE     OLS pair MSE (even if 1 is off-grid).\n"
      "  Picard_CV_Index       Zero-based index into the fraction grid.\n"
      "These maps and optional grid MSE maps use -1 outside the mask.\n"
      "\n"

      "Scores are in supplied coefficient units. Use comparable trial scaling\n"
      "and response units across runs. There is no session normalization or\n"
      "post-fit beta rescaling. Zero responses tie and select the largest\n"
      "fraction. These are SELECTION losses, not independent prediction\n"
      "scores or p-values. This implements GLMsingle's repeated-trial\n"
      "pair-scoring idea with fold-local fraction fits; it is not the full\n"
      "GLMsingle pipeline. Add -single_trial for the explicit staged HRF\n"
      "and noise-PC workflow. The first duty is still truth.\n"
      "\n"

      "PICARD_TRIAL_INFO records column/trial/condition/run identities;\n"
      "PICARD_FRACTIONS records grid order; PICARD_CV_FOLDS records each\n"
      "omitted run and its pair count. -write_manifest adds trial identities\n"
      "to the beta TSV. All identities otherwise travel in the main bucket\n"
      "header.\n"
      "\n"

      "SHARED-PREDICTOR CV -- CHART A COURSE THROUGH NEW DATA\n"
      "---------------------------------------------------\n"
      "  -cv_predict          Tune shared predictors on held-out time points.\n"
      "                       Use -fractions LIST; omit -trial_info.\n"
      "  -cv_block_size N     Hold out blocks of N original time steps.\n"
      "                       Default: leave one complete run out.\n"
      "  -cv_folds FILE       Two columns: full-row ID, fold ID. Include each\n"
      "                       retained row once; # comments allowed. Choose\n"
      "                       this OR -cv_block_size.\n"
      "  -cv_gap N            Exclude training rows within N original time\n"
      "                       steps of test rows in the same run. Default: 0.\n"
      "  -cv_nuisance MODE    train (default) or project; see below.\n"
      "  -cv_rule RULE        min (default) or 1se; see below.\n"
      "  -cv_global           One fraction over the explicit -mask. Default:\n"
      "                       select separately per voxel.\n"
      "\n"
      "This mode needs shared predictor columns with training support in\n"
      "every fold. Trial-specific columns belong to repeated-condition tuning\n"
      "above. Each retained row is tested once. Blocks start at the ORIGINAL\n"
      "RunStart; censoring never compresses block boundaries or gap\n"
      "distances. Fold IDs in -cv_folds are arbitrary nonnegative integers,\n"
      "ordered by first occurrence. There must be >=2 nonempty test folds\n"
      "with nonempty training sets.\n"
      "\n"
      "Nuisance mode train predicts X_test*beta + Z_test*gamma using training\n"
      "coefficients only. Test nuisance rows must be estimable from training.\n"
      "A baseline unique to the held-out run therefore causes a clear error.\n"
      "\n"
      "Nuisance mode project explicitly scores P(Z_test)*Y_test against\n"
      "P(Z_test)*X_test*beta. It removes the supplied held-out nuisance span\n"
      "from both target and prediction. This is CONDITIONAL signal\n"
      "prediction, using held-out Y to remove nuisance effects. It does not\n"
      "learn beta, lambda, or new nuisance vectors from held-out data. No\n"
      "mode adds or rescales columns.\n"
      "\n"
      "Score each candidate by total held-out SSE divided by retained row\n"
      "count. The baseline is the training response mean through the same\n"
      "test projection. R2 = 1 - SSE/baseline_SSE; negative R2 is retained.\n"
      "These are tuning scores, not independently evaluated outer-test\n"
      "performance. No nested evaluation or conventional inferential t/p maps\n"
      "are produced in this stage.\n"
      "\n"
      "Rule min selects the smallest pooled MSE; exact ties choose the\n"
      "largest fraction. Rule 1se chooses the SMALLEST fraction within best\n"
      "MSE + SE(best). For fold weights w = ntest/n, the weighted-fold\n"
      "heuristic is:\n"
      "\n"
      "  SE^2 = sum(w*(fold_MSE-mean)^2)*sum(w^2)/(1-sum(w^2)).\n"
      "\n"
      "Overlapping training folds make this a selection rule, not a\n"
      "confidence interval. The selected fraction is refit on the complete\n"
      "supplied design.\n"
      "\n"
      "Global selection requires an explicit mask. Normalize each voxel's\n"
      "fold losses by its pooled baseline MSE, then average over masked\n"
      "voxels with positive baseline MSE. Exclude zero baselines; if all are\n"
      "zero, choose the largest fraction. A serial prepass preserves voxel\n"
      "accumulation order; final fitting uses -jobs normally. Global\n"
      "selection need not run faster.\n"
      "\n"
      "Prediction diagnostics replace the repeat-tuning diagnostics:\n"
      "  Picard_PredMSE       Selected candidate's per-voxel held-out MSE.\n"
      "  Picard_PredBaseline  Per-voxel training-mean baseline MSE.\n"
      "  Picard_PredR2        Selected prediction R2 (raw or conditional).\n"
      "  Picard_PredIndex     Zero-based index into the fraction grid.\n"
      "  Picard_PredStatus    0 = defined R2; 1 = zero baseline.\n"
      "\n"
      "Status 1 stores R2 as 0. All these maps use -1 outside the mask.\n"
      "Status distinguishes missing R2 from a valid negative R2. -save_cv\n"
      "appends candidate MSEs to the main bucket. PICARD_PRED_* header\n"
      "attributes record folds, rule, nuisance policy, and any global\n"
      "selection/normalization. The captain keeps a reproducible log.\n"
      "\n"

      "CANDIDATE DESIGNS -- INFINITE DIVERSITY\n"
      "-------------------------------------\n"
      "  -candidate_manifest FILE   TSV input; replaces -matrix.\n"
      "  -candidate_policy hrf      Choose by full-data OLS fit MSE,\n"
      "                             then fit/tune the selected design.\n"
      "  -candidate_policy predict  Tune each design with -cv_predict;\n"
      "                             choose its smallest selected MSE.\n"
      "  -save_candidates           Add each design's MSE to the bucket.\n"
      "\n"
      "HRF policy uses minimum OLS residual MSE before applying any\n"
      "requested penalty. With the common nuisance span and response\n"
      "denominator, this selects the same winner as maximum fit R2.\n"
      "Candidates must have equal projected rank. Use a fixed\n"
      "fraction/lambda, or a parameter grid with -trial_info. HRF choice\n"
      "uses ALL retained data. Subsequent beta tuning is conditional on\n"
      "that choice; it is NOT independent validation of the selected HRF.\n"
      "This is not the full GLMsingle pipeline.\n"
      "\n"
      "Prediction policy requires -cv_predict and a fraction/lambda grid.\n"
      "Every candidate prepares its own training factors/scales. Each\n"
      "applies the requested min/1se rule to its parameter grid; the design\n"
      "with the smallest resulting prediction MSE wins. One-SE applies\n"
      "WITHIN designs. The same validation data select parameters and\n"
      "designs: reported MSE and R2 are SELECTION scores, not independent\n"
      "outer-test performance.\n"
      "\n"
      "Candidate selection is voxelwise; -cv_global is not supported for a\n"
      "collection. HRF policy cannot be combined with -cv_predict. Exact\n"
      "score ties choose the smallest stable ID, independent of manifest\n"
      "order. No HRFs are generated here: supply already constructed\n"
      "candidate matrices. Ordinary -matrix fitting remains independent of\n"
      "this collection interface.\n"
      "\n"
      "The first manifest line must be exactly:\n"
      "  # 3dPicard candidates 1\n"
      "The next noncomment line is the literal TAB-separated header:\n"
      "  id matrix columns source version sampling duration normalization\n"
      "Each row has those eight fields. Paths resolve relative to the\n"
      "manifest; spaces are literal. No quoting, shell expansion or escape\n"
      "processing. IDs are unique integers 0..16777215; 2..4096 designs are\n"
      "supported. source/version identify the HRF library or custom design\n"
      "revision. sampling is the source HRF sample interval in seconds;\n"
      "duration and normalization describe how the supplied matrices were\n"
      "constructed. These are provenance declarations, not instructions to\n"
      "convolve/rescale. HRF policy requires all five fields. Prediction\n"
      "may use '-' for unknown provenance. Duration and normalization must\n"
      "agree across candidates.\n"
      "\n"
      "Each columns file has this literal TAB-separated header:\n"
      "  column label role trial condition run basis item\n"
      "Provide every matrix column exactly once, in any row order. column\n"
      "and basis are nonnegative integer IDs. role is penalized or\n"
      "nuisance; label must exactly match the xmat label.\n"
      "trial/condition/run are integer IDs or '-' when absent; item is a\n"
      "supplied identifier or '-'. Trial columns need a condition and run;\n"
      "(trial,basis) pairs must be unique. Declared run-local columns must\n"
      "have no nonzero support in other runs. When using -trial_info, its\n"
      "trial/condition/run IDs must agree and basis must be 0.\n"
      "\n"
      "All designs must have matching dimensions, retained/full row\n"
      "mappings, run boundaries, TR, column roles/labels/identities, and\n"
      "stimulus ranges. Nuisance VALUES must also match. Different numbers\n"
      "of basis components are rejected. Candidate HRF formulas and\n"
      "penalized values may differ. Keep coefficient units comparable\n"
      "through the normalization declaration. These checks validate\n"
      "supplied identities, not the scientific accuracy of an external\n"
      "generator's HRF/provenance declarations.\n"
      "\n"
      "The main bucket gains Picard_DesignID, Picard_DesignIndex (manifest\n"
      "row, zero-based), and Picard_DesignMSE after the usual\n"
      "diagnostics/CV maps. -save_candidates appends Picard_DesignMSE_ID...\n"
      "maps in manifest order. All these maps use -1 outside the mask.\n"
      "Headers store the complete candidate catalog, column identities and\n"
      "selection policy. No extra score files are required. Betas keep\n"
      "common column labels for every voxel.\n"
      "\n"
      "For collections, -write_manifest writes Column and\n"
      "CandidateManifest, plus DesignIDFile/DesignIDVolume. It does not\n"
      "claim a single Xmat for voxelwise mixed designs.\n"
      "Trial/Condition/Run/Basis/Item remain attached to the common beta\n"
      "identity; missing numeric identities are written -1. The supplied\n"
      "matrix determines the meaning of each beta. There are no unannounced\n"
      "changes of identity aboard this ship.\n"
      "\n"

      "COMPUTATION -- ENGINEERING REPORT\n"
      "--------------------------------\n"
      "  -batch_size N        Responses per worker; default: 64.\n"
      "  -jobs N              Worker count; default: 1.\n"
      "                       More than one requires an OpenMP build.\n"
      "  -memory_mb M         Payload/work-buffer budget, MB; default: 2048.\n"
      "  -rank_tol T          Relative singular-value cutoff; default: 1e-10.\n"
      "\n"

      "Output arrays are float; factorization and fitting use doubles.\n"
      "The existing xmat reader parses floats before conversion to double.\n"
      "Input payload and requested output maps remain resident in memory.\n"
      "Response scratch is bounded by batch_size*jobs.\n"
      "\n"

      "memory_mb is a conservative payload/scratch estimate, not a hard\n"
      "process-memory limit: AFNI allocator bookkeeping and compression\n"
      "buffers are excluded. Even the Enterprise has finite resources.\n"
      "\n"

      "CHECKS AND EXAMPLES -- MAKE IT SO\n"
      "-------------------------------\n"
      "  -validate_only       Validate/show design; write nothing.\n"
      "                       Input is optional in this mode.\n"
      "  -show_design         Synonym for -validate_only.\n"
      "  -help                Show this help (also --help).\n"
      "\n"

      "Inspect the course before engaging:\n"
      "  3dPicard -matrix X.xmat.1D -validate_only\n"
      "\n"

      "Fit a fraction of 0.5; write the main bucket and residual bucket:\n"
      "  3dPicard -input epi+orig -matrix X.xmat.1D -fraction .5 \\\n"
      "           -prefix stats.picard -errts errts.picard\n"
      "\n"

      "Tune using explicit trial identities; keep scores in the main bucket:\n"
      "  3dPicard -input epi+orig -matrix X.xmat.1D \\\n"
      "           -fractions '1,.8,.6,.4,.2,0' -trial_info trials.1D \\\n"
      "           -prefix stats.picard -save_cv\n"
      "\n"

      "Predict shared stimulus effects with supplied run-local baselines:\n"
      "  3dPicard -input epi+orig -matrix X.xmat.1D \\\n"
      "           -cv_predict -fractions '1,.8,.6,.4,.2,0' \\\n"
      "           -cv_nuisance project -prefix stats.predict\n"
      "\n"
      "Fit a lasso model using explicit nuisance columns from the matrix:\n"
      "  3dPicard -input epi+orig -matrix X.xmat.1D -penalty lasso \\\n"
      "           -lambda .1 -standardize -prefix stats.lasso\n"
      "\n"
      "Tune elastic net for shared predictors, with conditional scoring:\n"
      "  3dPicard -input epi+orig -matrix X.xmat.1D -penalty enet \\\n"
      "           -l1_ratio .5 -lambdas '1,.3,.1,.03' -standardize \\\n"
      "           -cv_predict -cv_nuisance project -prefix stats.enet\n"
      "\n"
      "Select a supplied HRF, then tune repeated-condition betas:\n"
      "  3dPicard -input epi+orig -candidate_manifest hrfs.tsv \\\n"
      "           -candidate_policy hrf -fractions '1,.8,.5' \\\n"
      "           -trial_info trials.1D -prefix stats.hrf -save_candidates\n"
      "\n"
      "Compare shared-predictor designs using conditional prediction:\n"
      "  3dPicard -input epi+orig -candidate_manifest designs.tsv \\\n"
      "           -candidate_policy predict -cv_predict \\\n"
      "           -fractions '1,.8,.5' \\\n"
      "           -cv_nuisance project -prefix stats.designs\n"
      "\n"
      "Two conditions and two runs: 1dXmatMaker -> 3dPicard\n"
      "\n"
      "  Assume two preprocessed runs, run1.nii.gz and run2.nii.gz, each with\n"
      "  120 volumes and TR = 2 seconds. Each timing file has one line per run;\n"
      "  onsets are in seconds relative to that run's start. Use new output names.\n"
      "\n"
      "    # Create timing files: line 1 = run 1, line 2 = run 2.\n"
      "    printf '10 30 50 70 90 110 130 150 170 190\\n12 32 52 72 92 112 132 152 172 192\\n' \\\n"
      "      > frequent.1D\n"
      "    printf '60 120 180\\n62 122 182\\n' > target.1D\n"
      "\n"
      "    # Generate condition regressors and separate polynomial baselines per run.\n"
      "    1dXmatMaker -run_lengths 120,120 -TR 2 -polort 2 \\\n"
      "      -num_stimts 2 -local_times \\\n"
      "      -stim_times 1 frequent.1D GAM -stim_label 1 frequent \\\n"
      "      -stim_times 2 target.1D GAM -stim_label 2 target \\\n"
      "      -x1D frequent_target.xmat.1D\n"
      "\n"
      "    # Quote the input list so AFNI concatenates the runs in the stated order.\n"
      "    3dPicard -input 'run1.nii.gz run2.nii.gz' \\\n"
      "      -matrix frequent_target.xmat.1D \\\n"
      "      -fraction 0.6 \\\n"
      "      -prefix frequent_target.picard.nii.gz \\\n"
      "      -fitts frequent_target.fitts.nii.gz \\\n"
      "      -errts frequent_target.errts.nii.gz\n"
      "\n"
      "  The matrix has eight columns: six unpenalized baseline columns (constant,\n"
      "  linear and quadratic for each run), plus the frequent and target columns.\n"
      "  Picard estimates one coefficient per condition shared across both runs;\n"
      "  it uses the baseline columns already in the matrix rather than adding polort.\n"
      "  The input order must match the timing-file rows and -run_lengths. Run\n"
      "  boundaries come from the matrix. Quoted-list concatenation uses the standard\n"
      "  AFNI reader settings (AFNI_PATH_SPACES_OK must not be YES).\n"
      "  The fraction 0.6 is an illustrative fixed choice, not a tuned value.\n"
      "\n"
      "End of captain's log. The Ressikan flute is safely in its case.\n"
   );
}
// clang-format on

/* Accept a complete finite number; reject suffixes and range errors. */
static int parse_finite_double(const char *text, double *v)
{
   char *end;
   errno = 0;
   *v = strtod(text, &end);
   return !errno && end != text && !*end && isfinite(*v);
}

/* Parse counts without truncating a fractional value or overflowing int. */
static int parse_positive_integer(const char *text, int *v)
{
   double d;
   if (!parse_finite_double(text, &d) || d < 1 || d > INT_MAX || d != floor(d))
      return 0;
   *v = (int)d;
   return 1;
}

/* Estimator diagnostics precede the unchanged tuning diagnostics. */
static int base_diagnostic_count(const App *a)
{
   if (a->opt.banded)
      return 5 + a->opt.ngroup;
   return a->opt.elastic ? NELASTIC_DIAG + 2 * a->opt.support_refit : NDIAG;
}

/* Generated PC coefficients follow the original nuisance coefficients. */
static int saved_nuisance_count(const App *a)
{
   return a->opt.save_nuis ? THD_picard_plan_info(a->plan)->k +
                                 (a->single ? a->single->added : 0)
                           : 0;
}
static int single_extra_count(const App *a)
{
   return a->opt.single
              ? 4 + (a->opt.save_stages
                         ? 2 * (a->plan ? THD_picard_plan_info(a->plan)->p
                                        : a->xm->ncol)
                         : 0)
              : 0;
}

/* Count additional maps without changing the fixed-fraction bucket layout. */
static int diagnostic_count(const App *a)
{
   return single_extra_count(a) + base_diagnostic_count(a) +
          (a->opt.predict ? NPRED_DIAG
           : a->opt.ngrid ? NTUNE_DIAG
                          : 0) +
          (a->opt.save_cv ? a->opt.ngrid : 0) +
          (a->manifest ? NDESIGN_DIAG +
                             (a->opt.save_candidates ? a->manifest->count : 0)
                       : 0);
}

/* Parse an explicit grid; reject empty fields rather than silently changing
   the candidate set. The numerical plan also validates range and uniqueness. */
static int parse_fraction_grid(Options *o, char *err, int nerr)
{
   char *copy, *cursor;
   const char *list = o->elastic ? o->lambda_list : o->fraction_list;
   int count = 1, i;
   if (!list)
      return 1;
   for (cursor = (char *)list; *cursor; ++cursor)
      count += *cursor == ',';
   if (count < 2)
      return failure(err, nerr, "tuning grid needs >=2 values");
   o->grid = calloc((size_t)count, sizeof(double));
   copy = strdup(list);
   if (!o->grid || !copy) {
      free(copy);
      return failure(err, nerr, "cannot allocate fraction grid");
   }
   cursor = copy;
   for (i = 0; i < count; ++i) {
      char *next = strchr(cursor, ',');
      if (next)
         *next++ = '\0';
      if (!parse_finite_double(cursor, o->grid + i) ||
          (o->elastic ? o->grid[i] <= 0 : (o->grid[i] < 0 || o->grid[i] > 1))) {
         free(copy);
         return failure(
             err, nerr,
             "grid needs comma-separated [0,1] fractions or positive lambdas");
      }
      cursor = next;
   }
   free(copy);
   o->ngrid = count;
   return 1;
}

/* Small, explicit numeric tables use doubles, not the float xmat parser.
   Whitespace separates columns; # starts a comment. Reject ragged rows,
   nonfinite values and oversized lines before interpreting group identities.
   Successful output belongs to the caller; failures free partial storage. */
static int read_numeric_table(const char *path, double **values, int *nr,
                              int *nc, double memory, char *err, int nerr)
{
   FILE *f = fopen(path, "r");
   char line[65536];
   double *data = NULL;
   size_t used = 0, capacity = 0;
   int rows = 0, cols = 0, ok = 1;
   if (!f)
      return failure(err, nerr, "cannot open banded numeric table");
   while (ok && fgets(line, sizeof(line), f)) {
      char *cursor = line, *end;
      int count = 0;
      if (!strchr(line, '\n') && !feof(f)) {
         ok = failure(err, nerr, "banded table line too long");
         break;
      }
      while (ok) {
         double v;
         while (isspace((unsigned char)*cursor))
            cursor++;
         if (!*cursor || *cursor == '#')
            break;
         errno = 0;
         v = strtod(cursor, &end);
         if (end == cursor || errno || !isfinite(v) ||
             (*end && *end != '#' && !isspace((unsigned char)*end))) {
            ok = failure(err, nerr, "invalid/nonfinite banded table value");
            break;
         }
         if (used == capacity) {
            size_t next = capacity ? capacity * 2 : 64;
            double *grown;
            if (next > INT_MAX || next > SIZE_MAX / sizeof(double) ||
                (long double)next * sizeof(double) > memory * 1048576.0L) {
               ok = failure(err, nerr, "banded table exceeds memory_mb");
               break;
            }
            grown = realloc(data, next * sizeof(double));
            if (!grown) {
               ok = failure(err, nerr, "cannot allocate banded table");
               break;
            }
            data = grown;
            capacity = next;
         }
         data[used++] = v;
         count++;
         cursor = end;
      }
      if (count && ok) {
         if (rows && count != cols)
            ok = failure(err, nerr, "ragged banded table");
         else {
            cols = count;
            rows++;
         }
      }
   }
   if (ferror(f))
      ok = failure(err, nerr, "cannot read banded table");
   fclose(f);
   if (ok && !rows)
      ok = failure(err, nerr, "empty banded table");
   if (!ok) {
      free(data);
      return 0;
   }
   *values = data;
   *nr = rows;
   *nc = cols;
   return 1;
}

/* Parse arguments into borrowed argv strings; only the numeric grid is owned.
   Validate option combinations before touching input datasets. */
static int parse_options(int argc, char **argv, Options *o, char *err, int nerr)
{
   int i;
   memset(o, 0, sizeof(*o));
   o->fraction = NAN;
   o->noise_r2 = NAN;
   o->fixed_pcs = -1;
   o->pc_stop = 1.05;
   o->lambda = NAN;
   o->penalty = "ridge";
   o->solver.rho = 0.5;
   o->solver.tolerance = 1.e-8;
   o->solver.max_iterations = 10000;
   o->tol = THD_FRACTION_RANK_TOL;
   o->prefix = "Picard";
   o->batch = 64;
   o->jobs = 1;
   o->memory = 2048;
   for (i = 1; i < argc; i++) {
      char *arg = argv[i], *v;
      if (!strcmp(arg, "-validate_only") || !strcmp(arg, "-show_design")) {
         o->validate = 1;
         continue;
      }
      if (!strcmp(arg, "-save_candidates")) {
         o->save_candidates = 1;
         continue;
      }
      if (!strcmp(arg, "-single_trial")) {
         o->single = 1;
         continue;
      }
      if (!strcmp(arg, "-save_stages")) {
         o->save_stages = 1;
         continue;
      }
      if (!strcmp(arg, "-support_refit")) {
         o->support_refit = 1;
         continue;
      }
      if (!strcmp(arg, "-standardize")) {
         o->solver.standardize = 1;
         o->solver_set = 1;
         continue;
      }
      if (!strcmp(arg, "-cv_predict")) {
         o->predict = 1;
         continue;
      }
      if (!strcmp(arg, "-cv_global")) {
         o->global = 1;
         continue;
      }
      if (!strcmp(arg, "-save_cv")) {
         o->save_cv = 1;
         continue;
      }
      if (!strcmp(arg, "-write_manifest")) {
         o->manifest = 1;
         continue;
      }
      if (!strcmp(arg, "-save_nuisance")) {
         o->save_nuis = 1;
         continue;
      }
      if (!strcmp(arg, "-save_fitted")) {
         o->save_fit = 1;
         continue;
      }
      if (!strcmp(arg, "-save_residuals")) {
         o->save_res = 1;
         continue;
      }
      if (i + 1 >= argc)
         return failure(err, nerr, "option needs a value (see -help)");
      v = argv[++i];
      if (!strcmp(arg, "-noise_pool"))
         o->noise_pool = v;
      else if (!strcmp(arg, "-noise_model"))
         o->noise_model = v;
      else if (!strcmp(arg, "-noise_search"))
         o->noise_search = v;
      else if (!strcmp(arg, "-pc_mask"))
         o->pc_mask = v;
      else if (!strcmp(arg, "-pc_counts"))
         o->pc_list = v;
      else if (!strcmp(arg, "-noise_pcs")) {
         double value;
         if (!parse_finite_double(v, &value) || value < 0 || value > INT_MAX ||
             value != floor(value))
            return failure(err, nerr, "invalid noise_pcs");
         o->fixed_pcs = (int)value;
      } else if (!strcmp(arg, "-noise_r2")) {
         if (!parse_finite_double(v, &o->noise_r2) || o->noise_r2 < 0 ||
             o->noise_r2 > 1)
            return failure(err, nerr, "noise_r2 must be in [0,1]");
      } else if (!strcmp(arg, "-pc_stop")) {
         if (!parse_finite_double(v, &o->pc_stop) || o->pc_stop < 1)
            return failure(err, nerr, "pc_stop must be >=1");
         o->pc_stop_set = 1;
      } else if (!strcmp(arg, "-matrix"))
         o->matrix = v;
      else if (!strcmp(arg, "-candidate_manifest"))
         o->candidate_manifest = v;
      else if (!strcmp(arg, "-candidate_policy")) {
         if (strcmp(v, "hrf") && strcmp(v, "predict"))
            return failure(err, nerr,
                           "candidate_policy must be hrf or predict");
         o->candidate_policy = v;
         o->hrf_policy = !strcmp(v, "hrf");
      } else if (!strcmp(arg, "-penalty")) {
         if (strcmp(v, "ridge") && strcmp(v, "lasso") && strcmp(v, "enet") &&
             strcmp(v, "banded"))
            return failure(err, nerr,
                           "penalty must be ridge, lasso, enet or banded");
         o->penalty = v;
         o->elastic = !strcmp(v, "lasso") || !strcmp(v, "enet");
         o->banded = !strcmp(v, "banded");
      } else if (!strcmp(arg, "-band_groups"))
         o->group_file = v;
      else if (!strcmp(arg, "-band_penalties"))
         o->band_file = v;
      else if (!strcmp(arg, "-covariance")) {
         if (!strcmp(v, "white")) o->covariance.model = 0;
         else if (!strcmp(v, "ar1")) o->covariance.model = 1;
         else if (!strcmp(v, "arma11")) o->covariance.model = 2;
         else return failure(err, nerr, "covariance: white|ar1|arma11");
         o->covariance_set = 1;
      } else if (!strcmp(arg, "-covariance_ar") ||
                 !strcmp(arg, "-covariance_ma")) {
         double *value = !strcmp(arg, "-covariance_ar") ?
                            &o->covariance.a : &o->covariance.b;
         if (!parse_finite_double(v, value) || fabs(*value) >= 1)
            return failure(err, nerr, "covariance parameter outside (-1,1)");
         o->covariance_parameters = 1;
      } else if (!strcmp(arg, "-covariance_gaps")) {
         if (strcmp(v, "preserve") && strcmp(v, "reset"))
            return failure(err, nerr, "covariance_gaps: preserve|reset");
         o->covariance.reset_gaps = !strcmp(v, "reset");
         o->covariance_gap = 1;
      } else if (!strcmp(arg, "-covariance_estimate"))
         o->covariance_estimate = v;
      else if (!strcmp(arg, "-covariance_rows"))
         o->covariance_rows = v;
      else if (!strcmp(arg, "-covariance_map"))
         o->covariance_map = v;
      else if (!strcmp(arg, "-covariance_source"))
         o->covariance_source = v;
      else if (!strcmp(arg, "-prewhite_ar1")) {
         o->legacy_ar1_set = 1;
         if (!parse_finite_double(v, &o->ar1) || fabs(o->ar1) >= 1)
            return failure(err, nerr, "prewhite_ar1 must be in (-1,1)");
      } else if (!strcmp(arg, "-lambda")) {
         if (!parse_finite_double(v, &o->lambda) || o->lambda <= 0)
            return failure(err, nerr, "lambda must be finite and positive");
      } else if (!strcmp(arg, "-lambdas"))
         o->lambda_list = v;
      else if (!strcmp(arg, "-l1_ratio")) {
         if (!parse_finite_double(v, &o->solver.rho) || o->solver.rho < 0 ||
             o->solver.rho > 1)
            return failure(err, nerr, "l1_ratio must be in [0,1]");
         o->rho_set = 1;
      } else if (!strcmp(arg, "-solver_tol")) {
         if (!parse_finite_double(v, &o->solver.tolerance) ||
             o->solver.tolerance <= 0)
            return failure(err, nerr, "solver_tol must be positive");
         o->solver_set = 1;
      } else if (!strcmp(arg, "-max_iterations")) {
         if (!parse_positive_integer(v, &o->solver.max_iterations) ||
             o->solver.max_iterations == INT_MAX)
            return failure(err, nerr, "invalid max_iterations");
         o->solver_set = 1;
      } else if (!strcmp(arg, "-fractions"))
         o->fraction_list = v;
      else if (!strcmp(arg, "-cv_folds"))
         o->fold_file = v;
      else if (!strcmp(arg, "-cv_block_size")) {
         if (!parse_positive_integer(v, &o->block_size))
            return failure(err, nerr,
                           "cv_block_size must be a positive integer");
      } else if (!strcmp(arg, "-cv_gap")) {
         double value;
         if (!parse_finite_double(v, &value) || value < 0 || value > INT_MAX ||
             value != floor(value))
            return failure(err, nerr, "cv_gap must be a nonnegative integer");
         o->gap = (int)value;
         o->gap_set = 1;
      } else if (!strcmp(arg, "-cv_nuisance")) {
         o->cv_nuisance = v;
         if (strcmp(v, "train") && strcmp(v, "project"))
            return failure(err, nerr, "cv_nuisance must be train or project");
         o->project = !strcmp(v, "project");
      } else if (!strcmp(arg, "-cv_rule")) {
         o->cv_rule = v;
         if (strcmp(v, "min") && strcmp(v, "1se"))
            return failure(err, nerr, "cv_rule must be min or 1se");
         o->one_se = !strcmp(v, "1se");
      } else if (!strcmp(arg, "-trial_info"))
         o->trial_info = v;
      else if (!strcmp(arg, "-input"))
         o->input = v;
      else if (!strcmp(arg, "-mask"))
         o->mask = v;
      else if (!strcmp(arg, "-prefix"))
         o->prefix = v;
      else if (!strcmp(arg, "-fitts")) {
         o->fitts = v;
         o->save_fit = 1;
      } else if (!strcmp(arg, "-errts")) {
         o->errts = v;
         o->save_res = 1;
      } else if (!strcmp(arg, "-penalty_cols"))
         o->pen = v;
      else if (!strcmp(arg, "-nuisance_cols"))
         o->nuis = v;
      else if (!strcmp(arg, "-fraction")) {
         if (!parse_finite_double(v, &o->fraction) || o->fraction < 0 ||
             o->fraction > 1)
            return failure(err, nerr, "fraction must be in [0,1]");
      } else if (!strcmp(arg, "-rank_tol")) {
         if (!parse_finite_double(v, &o->tol) || o->tol <= 0 || o->tol >= 1)
            return failure(err, nerr, "rank_tol must be in (0,1)");
      } else if (!strcmp(arg, "-memory_mb")) {
         if (!parse_finite_double(v, &o->memory) || o->memory <= 0)
            return failure(err, nerr, "memory_mb must be positive");
      } else if (!strcmp(arg, "-batch_size")) {
         if (!parse_positive_integer(v, &o->batch))
            return failure(err, nerr, "batch_size must be a positive integer");
      } else if (!strcmp(arg, "-jobs")) {
         if (!parse_positive_integer(v, &o->jobs))
            return failure(err, nerr, "jobs must be a positive integer");
      } else {
         snprintf(err, (size_t)nerr, "unknown option %s (see -help)", arg);
         return 0;
      }
   }
   if ((!o->elastic && (isfinite(o->lambda) || o->lambda_list || o->rho_set ||
                        o->solver_set)) ||
       (o->elastic && (isfinite(o->fraction) || o->fraction_list)))
      return failure(err, nerr,
                     "ridge uses fraction(s); lasso/enet use lambda(s) and "
                     "solver options");
   if (!strcmp(o->penalty, "lasso")) {
      if (o->rho_set && o->solver.rho != 1)
         return failure(err, nerr, "lasso requires l1_ratio=1");
      o->solver.rho = 1;
   }
   if ((!o->covariance_set && (o->covariance_parameters ||
        o->covariance_gap || o->covariance_map || o->covariance_source ||
        o->covariance_estimate || o->covariance_rows)) ||
       (o->covariance_set && o->legacy_ar1_set))
      return failure(err, nerr,
          "covariance options require -covariance; not with prewhite_ar1");
   if (o->covariance_set &&
       ((o->covariance.model == 0 && o->covariance.a != 0) ||
        (o->covariance.model < 2 && o->covariance.b != 0)))
      return failure(err, nerr, "parameters incompatible with covariance");
   if (o->covariance_map && (o->covariance_parameters ||
       o->covariance.model == 0 || o->single || o->global))
      return failure(err, nerr,
          "covariance_map: no scalar parameters, white, single_trial or "
          "global tuning; staged voxel covariance is deferred to NG2");
   if (o->covariance_source && !o->covariance_map)
      return failure(err, nerr, "covariance_source requires covariance_map");
   if (o->support_refit && (!o->elastic || o->solver.rho <= 0))
      return failure(err, nerr,
                     "support_refit requires lasso/enet with l1_ratio>0");
   if (o->banded) {
      int c, g, j;
      if (!o->group_file || !o->band_file || isfinite(o->fraction) ||
          o->fraction_list)
         return failure(
             err, nerr,
             "banded needs band_groups and band_penalties, not fractions");
      if (!read_numeric_table(o->band_file, &o->band_values, &o->nband,
                              &o->ngroup, o->memory, err, nerr))
         return 0;
      if (o->nband > 16777215)
         return failure(err, nerr, "too many banded candidates");
      for (c = 0; c < o->nband; c++) {
         for (g = 0; g < o->ngroup; g++)
            if (o->band_values[(size_t)c * o->ngroup + g] <= 0)
               return failure(err, nerr, "banded penalties must be positive");
         for (j = 0; j < c; j++)
            if (!memcmp(o->band_values + (size_t)c * o->ngroup,
                        o->band_values + (size_t)j * o->ngroup,
                        sizeof(double) * (size_t)o->ngroup))
               return failure(err, nerr, "duplicate banded penalty vector");
      }
      if (o->nband == 1)
         o->fraction = 0;
      else {
         o->grid = calloc((size_t)o->nband, sizeof(double));
         if (!o->grid)
            return failure(err, nerr, "cannot allocate banded grid");
         o->ngrid = o->nband;
         for (c = 0; c < o->nband; c++)
            o->grid[c] = c;
      }
   } else if (o->group_file || o->band_file)
      return failure(err, nerr,
                     "band_groups/band_penalties require penalty banded");
   /* One scalar/grid travels through common scheduling; its units are
      determined by the explicitly recorded estimator, never inferred. */
   if (o->elastic)
      o->fraction = o->lambda;
   if (!o->banded && !parse_fraction_grid(o, err, nerr))
      return 0;
   if ((o->ngrid && isfinite(o->fraction)) ||
       (!o->single && (!!o->ngrid != (!!o->trial_info || o->predict))) ||
       (o->trial_info && o->predict) || (o->save_cv && !o->ngrid))
      return failure(
          err, nerr,
          "use either a fixed fraction/lambda or its tuning grid "
          "with -trial_info FILE or -cv_predict; -save_cv requires tuning");
   if ((!o->predict && (o->fold_file || o->block_size || o->gap_set ||
                        o->cv_nuisance || o->cv_rule || o->global)) ||
       (o->fold_file && o->block_size))
      return failure(
          err, nerr,
          "CV options require -cv_predict; choose blocks OR fold file");
   if (o->global && !o->mask)
      return failure(err, nerr, "cv_global requires an explicit -mask");
   if (o->candidate_manifest) {
      if (o->matrix || !o->candidate_policy || o->global)
         return failure(
             err, nerr,
             "candidate_manifest replaces matrix, requires candidate_policy, "
             "and does not support cv_global");
      if ((o->hrf_policy && o->predict) || (!o->hrf_policy && !o->predict))
         return failure(
             err, nerr,
             "candidate policy hrf excludes cv_predict; predict requires it");
   } else if (o->candidate_policy || o->save_candidates)
      return failure(err, nerr, "candidate options require candidate_manifest");
   if ((!o->matrix && !o->candidate_manifest) ||
       (!o->validate && (!o->input || (!isfinite(o->fraction) && !o->ngrid))))
      return failure(err, nerr,
                     "matrix, input and fraction/lambda or grid are required");
   if (o->single) {
      if (!o->trial_info || o->predict || o->global ||
          (o->candidate_manifest && !o->hrf_policy) ||
          (!!o->noise_pool == !!o->noise_model) ||
          (o->noise_model && (!o->noise_search || !isfinite(o->noise_r2))) ||
          (o->noise_pool && (o->noise_search || isfinite(o->noise_r2))) ||
          (!!o->pc_list == (o->fixed_pcs >= 0)))
         return failure(err, nerr,
                        "single_trial needs trial_info, one pool policy and "
                        "noise_pcs OR pc_counts; "
                        "prediction mode is separate");
      if (o->pc_list) {
         char *copy = strdup(o->pc_list), *cursor = copy, *next;
         int count = 1, j;
         if (!copy)
            return failure(err, nerr, "cannot parse PC counts");
         for (j = 0; copy[j]; j++)
            count += copy[j] == ',';
         o->pc_counts = calloc((size_t)count, sizeof(int));
         o->npc = count;
         if (!o->pc_counts) {
            free(copy);
            return failure(err, nerr, "cannot allocate PC counts");
         }
         for (j = 0; j < count; j++) {
            double value;
            next = strchr(cursor, ',');
            if (next)
               *next++ = 0;
            if (!parse_finite_double(cursor, &value) || value < 0 ||
                value > INT_MAX || value != floor(value) ||
                (j && value <= o->pc_counts[j - 1])) {
               free(copy);
               return failure(err, nerr, "PC counts must strictly increase");
            }
            o->pc_counts[j] = (int)value;
            cursor = next;
         }
         free(copy);
         if (count < 2 || o->pc_counts[0] != 0 || !o->pc_mask)
            return failure(
                err, nerr,
                "PC selection needs counts starting at zero and pc_mask");
      } else {
         if (o->pc_mask || o->pc_stop_set)
            return failure(err, nerr,
                           "pc_mask/pc_stop require PC-count selection");
         o->pc_counts = calloc(1, sizeof(int));
         o->npc = 1;
         if (!o->pc_counts)
            return failure(err, nerr, "cannot allocate fixed PC count");
         o->pc_counts[0] = o->fixed_pcs;
      }
   } else if (o->noise_pool || o->noise_model || o->noise_search ||
              o->pc_mask || o->pc_list || o->fixed_pcs >= 0 ||
              isfinite(o->noise_r2) || o->save_stages || o->pc_stop_set)
      return failure(err, nerr, "noise/PC/stage options require single_trial");
   if (o->covariance_rows && !o->covariance_estimate)
      return failure(err, nerr, "covariance_rows requires covariance_estimate");
   if (o->covariance_estimate &&
       (o->covariance_map || o->covariance_parameters || o->single ||
        o->candidate_manifest || o->ngrid || o->global || !o->input))
      return failure(err, nerr,
          "native covariance requires input and fixed penalty; no maps, "
          "scalar parameters, candidates, staged mode or internal CV. "
          "For fold-local tuning use covariance_rows and explicit folds.");
   if (o->mask && !o->input)
      return failure(err, nerr, "mask requires input");
   if ((o->input && strchr(o->input, '[')) || (o->mask && strchr(o->mask, '[')))
      return failure(err, nerr, "input and mask must be unselected datasets");
   if (!THD_filename_ok(o->prefix) ||
       (o->fitts && !THD_filename_ok(o->fitts)) ||
       (o->errts && !THD_filename_ok(o->errts)))
      return failure(err, nerr, "invalid output prefix");
#ifndef USE_OMP
   if (o->jobs > 1)
      return failure(err, nerr, "this build has no OpenMP; use -jobs 1");
#endif
   return 1;
}

/* Expand one AFNI integer selector into flags, rejecting duplicates and
   out-of-range columns so roles cannot silently overlap. */
static int mark_selected_columns(const char *text, int n, int *flags, char *err,
                                 int nerr)
{
   NI_int_array *a;
   int j;
   if (!strcmp(text, "none"))
      return 1;
   if (!strcmp(text, "all")) {
      for (j = 0; j < n; j++)
         flags[j] = 1;
      return 1;
   }
   a = NI_decode_int_list((char *)text, " ,;");
   if (!a || a->num < 1) {
      if (a)
         NI_delete_int_array(a);
      return failure(err, nerr, "invalid column index list");
   }
   for (j = 0; j < a->num; j++) {
      int c = a->ar[j];
      if (c < 0 || c >= n || flags[c]) {
         NI_delete_int_array(a);
         return failure(err, nerr, "duplicate/out-of-range column index");
      }
      flags[c] = 1;
   }
   NI_delete_int_array(a);
   return 1;
}

/* Partition every design column into penalized or nuisance roles. Explicit
   lists override ColumnGroups, with a missing list filled by the complement. */
static int assign_column_roles(App *a, char *err, int nerr)
{
   int c, n = a->xm->ncol, np = 0;
   int *pen = calloc((size_t)n, sizeof(int)),
       *nuis = calloc((size_t)n, sizeof(int));
   int ok = pen && nuis;
   if (!ok) {
      free(pen);
      free(nuis);
      return failure(err, nerr, "cannot allocate column roles");
   }
   if (a->opt.pen || a->opt.nuis) {
      if (a->opt.pen)
         ok = mark_selected_columns(a->opt.pen, n, pen, err, nerr);
      if (ok && a->opt.nuis)
         ok = mark_selected_columns(a->opt.nuis, n, nuis, err, nerr);
      for (c = 0; c < n && ok; c++) {
         if (!a->opt.pen)
            pen[c] = !nuis[c];
         if (!a->opt.nuis)
            nuis[c] = !pen[c];
         if (pen[c] + nuis[c] != 1)
            ok = failure(err, nerr,
                         "column roles must be disjoint and cover the matrix");
      }
   } else {
      const char *g = THD_xmat_attr(a->xm, "ColumnGroups");
      NI_int_array *groups = g ? NI_decode_int_list((char *)g, " ,;") : NULL;
      if (!groups || groups->num != n)
         ok = failure(
             err, nerr,
             "missing/malformed ColumnGroups: supply explicit column roles");
      if (ok)
         for (c = 0; c < n; c++)
            pen[c] = groups->ar[c] > 0;
      if (groups)
         NI_delete_int_array(groups);
   }
   for (c = 0; c < n && ok; c++) {
      a->roles[c] = pen[c];
      np += pen[c];
   }
   free(pen);
   free(nuis);
   if (ok && !np)
      ok = failure(err, nerr, "at least one column must be penalized");
   return ok;
}

/* Limit inputs to scalar real-valued AFNI bricks, excluding complex/RGB data.
 */
static int is_real_datum(int type)
{
   return type == MRI_byte || type == MRI_short || type == MRI_int ||
          type == MRI_float || type == MRI_double;
}

/* Read one sample from an already-loaded AFNI sub-brick. The loader has
   handled file layout, byte order and any NIfTI slope/intercept conversion.
   Only the AFNI brick factor remains: zero means "no scaling", not zero data.

   Promote the stored sample to double before multiplying the factor. This
   avoids further float rounding here, but cannot recover precision lost in
   the file or AFNI loader. We never reopen the source file or repair samples;
   the numerical layer rejects any nonfinite retained values still present. */
static double read_voxel(THD_3dim_dataset *d, int t, int v)
{
   double x = 0, fac = DSET_BRICK_FACTOR(d, t);
   void *p = DSET_ARRAY(d, t);
   switch (DSET_BRICK_TYPE(d, t)) {
   case MRI_byte:
      x = ((byte *)p)[v];
      break;
   case MRI_short:
      x = ((short *)p)[v];
      break;
   case MRI_int:
      x = ((int *)p)[v];
      break;
   case MRI_float:
      x = ((float *)p)[v];
      break;
   case MRI_double:
      x = ((double *)p)[v];
      break;
   default:
      return NAN;
   }
   return x * (fac == 0 ? 1 : fac);
}

/* Groups describe ORIGINAL design columns, including -1 for nuisance.
   Convert to the compact penalized-column order only after roles are known.
   Explicit groups avoid conflating stimulus IDs with penalty families. */
static int prepare_groups(App *a, char *err, int nerr)
{
   double *table = NULL;
   int nr, nc, j, g, np = 0, at = 0, ok = 1;
   if (!a->opt.banded)
      return 1;
   if (!read_numeric_table(a->opt.group_file, &table, &nr, &nc, a->opt.memory,
                           err, nerr))
      return 0;
   if ((nr != 1 && nc != 1) || (size_t)nr * nc != (size_t)a->xm->ncol)
      ok = failure(err, nerr,
                   "band_groups must be one value per original column");
   for (j = 0; j < a->xm->ncol; j++)
      np += a->roles[j] != 0;
   if (a->opt.ngroup > np)
      ok = failure(err, nerr, "more band groups than predictors");
   a->groups = calloc((size_t)np, sizeof(int));
   if (!a->groups)
      ok = failure(err, nerr, "cannot allocate banded groups");
   for (j = 0; j < a->xm->ncol && ok; j++) {
      double value = table[j];
      if (value != floor(value) || (!a->roles[j] && value != -1) ||
          (a->roles[j] && (value < 0 || value >= a->opt.ngroup)))
         ok = failure(err, nerr,
                      "band_groups: nuisance=-1; predictor ID in 0..G-1");
      else if (a->roles[j])
         a->groups[at++] = (int)value;
   }
   for (g = 0; g < a->opt.ngroup && ok; g++) {
      for (j = 0; j < np && a->groups[j] != g; j++)
         ;
      if (j == np)
         ok = failure(err, nerr, "empty band group in full design");
   }
   free(table);
   return ok;
}

/* Check geometry, type and full timeline length before loading voxel data.
   Opening headers first lets the budget check run before large allocations. */
static int open_input_headers(App *a, char *err, int nerr)
{
   int t;
   if (!a->opt.input)
      return 1;
   a->input = THD_open_dataset(a->opt.input);
   if (!a->input || (DSET_IS_NI_SURF_DSET(a->input) ||
                     DSET_IS_GIFTI(a->input) || DSET_IS_1D(a->input)))
      return failure(err, nerr, "cannot open volume input");
   if (DSET_NVALS(a->input) != a->xm->nfull)
      return failure(err, nerr, "input time axis must match matrix NRowFull");
   a->nvox = DSET_NVOX(a->input);
   for (t = 0; t < DSET_NVALS(a->input); t++)
      if (!is_real_datum(DSET_BRICK_TYPE(a->input, t)))
         return failure(err, nerr, "input must contain real numeric bricks");
   if (a->opt.mask) {
      a->mask = THD_open_dataset(a->opt.mask);
      if (!a->mask || DSET_IS_NI_SURF_DSET(a->mask) || DSET_IS_GIFTI(a->mask) ||
          DSET_IS_1D(a->mask) || DSET_NVALS(a->mask) != 1 ||
          !is_real_datum(DSET_BRICK_TYPE(a->mask, 0)) ||
          !EQUIV_GRIDS(a->input, a->mask) ||
          THD_dataset_mismatch(a->input, a->mask))
         return failure(err, nerr,
                        "mask must be one real volume on the input grid");
   }
   return 1;
}

/* Bound resident image/map payloads, immutable design factors and private
   worker scratch. Long double arithmetic avoids overflowing size estimates.
   This deliberately excludes AFNI allocator and compression overhead. */
static int check_memory_budget(App *a, char *err, int nerr)
{
   int j, np = 0, k;
   long double n = a->xm->nrow, d = a->xm->ncol, v = a->nvox, bytes;
   for (j = 0; j < a->xm->ncol; j++)
      np += a->roles[j];
   k = a->xm->ncol - np;
   bytes = 1048576 + 128 * n * d;
   if (!a->opt.validate)
      bytes += 8 * v * a->xm->nfull + 32 * v +
               8.0L * a->opt.batch * a->opt.jobs * (6 * n + 4 * d + 32) +
               4 * v *
                   (np + diagnostic_count(a) + (a->opt.save_nuis ? k : 0) +
                    (a->opt.save_fit + a->opt.save_res) * n);
   /* Fold designs are immutable and shared; each worker owns fold scratch.
      Use the full size for every fold as a conservative pre-allocation bound.
    */
   if (a->opt.ngrid) {
      long double folds = a->folds ? a->folds->nfold : a->xm->nrun;
      long double q = a->opt.batch;
      bytes += 160 * folds * n * d + 64 * d + 16 * a->opt.ngrid;
      bytes += 16 * folds * n + 32 * folds * a->opt.ngrid;
      bytes += 8.0L * a->opt.jobs * folds * a->opt.ngrid;
      if (a->opt.global) /* Serial scoring workspace for the global prepass. */
         bytes +=
             8 * q *
             (folds * (8 * n + 4 * d + 64 + a->opt.ngrid) + 8 * n + 4 * d + 64);
      if (a->opt.predict)
         bytes += 8 * q * a->opt.jobs * (folds * (3 * n + 8) + n);
      if (!a->opt.validate)
         bytes += 8 * q * a->opt.jobs *
                  (folds * (6 * n + 4 * d + 32) + n + d + 8 + a->opt.ngrid +
                   folds * (a->opt.ngrid + 1));
   }
   if (a->opt.elastic) {
      long double plans =
          1 + (a->opt.ngrid ? (a->folds ? a->folds->nfold : a->xm->nrun) : 0);
      bytes += plans * (8 * n * d + 16 * d);
      bytes += plans * (a->opt.jobs + (a->opt.global ? 1 : 0)) *
               (8 * n + a->opt.batch * (16 * d + 128));
   }
   {
      long double plans =
          1 + (a->opt.ngrid ? (a->folds ? a->folds->nfold : a->xm->nrun) : 0);
      long double workers = a->opt.jobs + (a->opt.global ? 1 : 0);
      long double size = n < np ? n : np, q = a->opt.batch;
      bytes += plans * (8 * n * d + 8 * workers * n * q);
      if (a->opt.banded)
         bytes +=
             plans * (8 * n * np + 8 * size * size + 16 * np +
                      16.0L * a->opt.nband * a->opt.ngroup +
                      workers * (8 * size * size + 8 * size * q + 8 * n * q +
                                 16 * np * q + 64 * q + 8 * n));
      if (a->opt.support_refit)
         bytes +=
             plans * workers *
             (128 * n * np + 128 * np * np + 16 * np * q + 16 * n + 16 * q);
   }
   if (a->manifest) {
      long double count = a->manifest->count;
      long double shared =
          a->opt.validate
              ? 0
              : 8 * v * a->xm->nfull + 32 * v +
                    4 * v *
                        (np + diagnostic_count(a) + (a->opt.save_nuis ? k : 0) +
                         (a->opt.save_fit + a->opt.save_res) * n);
      /* Include the reference tuning plan used for metadata as well as each
         candidate's factors and workers. Image/output payloads stay shared. */
      bytes = shared + (bytes - shared) * (count + 1);
      bytes += a->manifest->metadata_bytes + count * d * 256;
      bytes += 8.0L * a->opt.jobs * a->opt.batch *
               (count * (6 * n + 4 * d + 32) + 2 * n + d + a->opt.ngrid +
                count + 32);
   }
   if (a->opt.single) {
      long double nr = a->xm->nrun, maxpc = a->opt.pc_counts[a->opt.npc - 1];
      long double added = nr * maxpc, h = a->manifest ? a->manifest->count : 1;
      if (maxpc > n || added > INT_MAX - a->xm->ncol)
         return failure(err, nerr, "noise PC count exceeds design dimensions");
      bytes +=
          80 * n * n + 32 * n * d + 16 * n * a->opt.batch +
          v * (64 + 8 * a->opt.npc + 8 * h + 16 * np * a->opt.save_stages) +
          160 * h * (nr + 1) * n * (d + added) +
          64 * h * (nr + 1) * a->opt.jobs * a->opt.batch * (n + d + added) +
          64 * n * maxpc + 4 * v * added * a->opt.save_nuis;
   }
   if (a->opt.covariance_set) {
      long double folds = a->opt.ngrid ?
          (a->folds ? a->folds->nfold : a->xm->nrun) : 0;
      long double candidates = a->manifest ? a->manifest->count + 1 : 1;
      /* Full/training factors plus staged augmentation and temporary plans. */
      bytes += 16*n*n*(folds+2)*candidates*(a->opt.single ? 4 : 1);
      if (a->opt.covariance_map || a->opt.covariance_estimate)
         /* One active group plus reference plans and map-loading reserve. */
         bytes = 2*bytes + 128*v + a->opt.memory*1048576.0L/4;
   }
   if (bytes > a->opt.memory * 1048576.0L) {
      snprintf(err, (size_t)nerr,
               "estimated payload/work buffers %.1Lf MB exceed -memory_mb %.1f",
               bytes / 1048576, a->opt.memory);
      return 0;
   }
   return 1;
}

/* Check the original metadata spelling as well as the xmat reader result.
   A permissive legacy decoder must not turn malformed censor maps into a
   different scientific design without reporting an error. */
static int check_integer_metadata(const THD_xmat *xm, const char *key,
                                  const int *expected, int count, char *err,
                                  int nerr)
{
   const char *text = THD_xmat_attr(xm, key);
   NI_int_array *a;
   int i, ok = 1;
   if (!text)
      return 1;
   for (i = 0; text[i]; i++) {
      unsigned char c = (unsigned char)text[i];
      if (!isdigit(c) && !isspace(c) && !strchr(",;+-().@", c))
         ok = 0;
      if (c == '.' && (i == 0 || text[i - 1] != '.') && text[i + 1] != '.')
         ok = 0;
   }
   a = ok ? NI_decode_int_list((char *)text, " ,;") : NULL;
   if (!a || a->num != count)
      ok = 0;
   if (ok && expected)
      for (i = 0; i < count; i++)
         if (a->ar[i] != expected[i])
            ok = 0;
   if (a)
      NI_delete_int_array(a);
   if (!ok) {
      snprintf(err, (size_t)nerr, "malformed %s metadata", key);
      return 0;
   }
   return 1;
}

/* The file is deliberately small and explicit: four nonnegative integer
   columns, with optional # comments. Sort by the design's penalized columns,
   never by condition or trial ID, so beta identities survive output ordering.
 */
static int read_trial_info(App *a, char *err, int nerr)
{
   const THD_picard_design_info *d = THD_picard_plan_info(a->plan);
   FILE *fp = fopen(a->opt.trial_info, "r");
   char line[4096];
   int line_number = 0, ok = 1, i;
   if (!fp)
      return failure(err, nerr, "cannot open -trial_info");
   a->trials = calloc((size_t)d->p, sizeof(*a->trials));
   if (!a->trials) {
      fclose(fp);
      return failure(err, nerr, "cannot allocate trial identities");
   }
   for (i = 0; i < d->p; ++i)
      a->trials[i].column = -1;
   while (ok && fgets(line, sizeof(line), fp)) {
      char *cursor = line, *comment = strchr(line, '#');
      int values[4], count = 0, index = -1;
      ++line_number;
      if (!strchr(line, '\n') && !feof(fp)) {
         ok = 0;
         break;
      }
      if (comment)
         *comment = '\0';
      while (*cursor && ok) {
         char *end;
         long value;
         while (isspace((unsigned char)*cursor))
            ++cursor;
         if (!*cursor)
            break;
         errno = 0;
         value = strtol(cursor, &end, 10);
         if (count == 4 || errno || end == cursor ||
             (*end && !isspace((unsigned char)*end)) || value < 0 ||
             value > INT_MAX) {
            ok = 0;
            break;
         }
         values[count++] = (int)value;
         cursor = end;
      }
      if (!count && ok)
         continue;
      if (count != 4) {
         ok = 0;
         break;
      }
      for (i = 0; i < d->p; ++i)
         if (d->penalized_columns[i] == values[0])
            index = i;
      if (index < 0 || a->trials[index].column != -1) {
         ok = 0;
         break;
      }
      a->trials[index].column = values[0];
      a->trials[index].trial = values[1];
      a->trials[index].condition = values[2];
      a->trials[index].run = values[3];
   }
   if (ferror(fp))
      ok = 0;
   if (fclose(fp))
      ok = 0;
   for (i = 0; i < d->p; ++i)
      if (a->trials[i].column == -1)
         ok = 0;
   if (!ok) {
      snprintf(err, (size_t)nerr,
               "invalid/incomplete trial_info near line %d; "
               "expected: column trial condition run (zero-based integers)",
               line_number);
      return 0;
   }
   return 1;
}

/* Explicit fold tables use original full-timeline rows, not compressed
   retained-row indices. Every retained row occurs exactly once. */
static int *read_prediction_groups(const App *a, char *err, int nerr)
{
   FILE *fp = fopen(a->opt.fold_file, "r");
   int n = a->xm->nrow, *groups = malloc((size_t)n * sizeof(int)), i, ok = 1;
   char line[4096];
   if (!fp || !groups) {
      if (fp)
         fclose(fp);
      free(groups);
      failure(err, nerr, "cannot open/allocate cv_folds");
      return NULL;
   }
   for (i = 0; i < n; ++i)
      groups[i] = -1;
   while (ok && fgets(line, sizeof(line), fp)) {
      char *cursor = line, *comment = strchr(line, '#');
      int fields[2], count = 0, index = -1;
      if (!strchr(line, '\n') && !feof(fp)) {
         ok = 0;
         break;
      }
      if (comment)
         *comment = '\0';
      while (*cursor && ok) {
         char *end;
         long value;
         while (isspace((unsigned char)*cursor))
            ++cursor;
         if (!*cursor)
            break;
         errno = 0;
         value = strtol(cursor, &end, 10);
         if (count == 2 || errno || end == cursor || value < 0 ||
             value > INT_MAX || (*end && !isspace((unsigned char)*end))) {
            ok = 0;
            break;
         }
         fields[count++] = (int)value;
         cursor = end;
      }
      if (!count && ok)
         continue;
      if (count != 2) {
         ok = 0;
         break;
      }
      for (i = 0; i < n; ++i)
         if (a->xm->good[i] == fields[0])
            index = i;
      if (index < 0 || groups[index] != -1) {
         ok = 0;
         break;
      }
      groups[index] = fields[1];
   }
   if (ferror(fp))
      ok = 0;
   if (fclose(fp))
      ok = 0;
   for (i = 0; i < n; ++i)
      if (groups[i] < 0)
         ok = 0;
   if (!ok) {
      free(groups);
      failure(err, nerr,
              "cv_folds needs each retained full-row ID "
              "exactly once, followed by a nonnegative fold ID");
      return NULL;
   }
   return groups;
}

/* Count groups before allocating their explicit O(rows*folds) index lists.
   This early bound complements the later design/worker payload budget. */
static int prepare_prediction_folds(App *a, char *err, int nerr)
{
   int *groups = NULL, *first, n = a->xm->nrow, nf = 0, i, f;
   if (!a->opt.predict)
      return 1;
   if (a->opt.fold_file) {
      groups = read_prediction_groups(a, err, nerr);
      if (!groups)
         return 0;
   }
   first = calloc((size_t)n, sizeof(int));
   if (!first) {
      free(groups);
      return failure(err, nerr, "cannot count folds");
   }
   for (i = 0; i < n; ++i) {
      for (f = 0; f < nf; ++f) {
         int j = first[f];
         if (groups ? groups[i] == groups[j]
                    : (a->runs[i] == a->runs[j] &&
                       (!a->opt.block_size ||
                        (a->xm->good[i] - a->xm->runstart[a->runs[i]]) /
                                a->opt.block_size ==
                            (a->xm->good[j] - a->xm->runstart[a->runs[j]]) /
                                a->opt.block_size)))
            break;
      }
      if (f == nf)
         first[nf++] = i;
   }
   free(first);
   if (16 * (long double)n * nf > a->opt.memory * 1048576.0L) {
      free(groups);
      return failure(err, nerr, "fold indices exceed memory budget");
   }
   a->folds = THD_foldplan_time_series(
       n, a->xm->good, a->runs, a->xm->nrun, a->xm->runstart, groups,
       a->opt.block_size, a->opt.gap, err, nerr);
   free(groups);
   return a->folds != NULL;
}

/* Validate the xmat contract, assign roles, then build immutable numerical
   plans. AFNI stores matrices column-major; thd_picard expects row-major.
   The temporary double matrix is freed once the plan has copied it. */
static int prepare_design(App *a, const App *reference, char *err, int nerr)
{
   int i, j;
   const char *full, *good;
   if (!a->xm)
      a->xm = THD_xmat_read(a->opt.matrix, err, nerr);
   if (!a->xm)
      return 0;
   full = THD_xmat_attr(a->xm, "NRowFull");
   good = THD_xmat_attr(a->xm, "GoodList");
   if ((full && !good && a->xm->nfull != a->xm->nrow) || (good && !full))
      return failure(err, nerr,
                     "ambiguous censoring: GoodList and NRowFull must be "
                     "supplied together");
   if (full) {
      double value;
      if (!parse_finite_double(full, &value) || value != a->xm->nfull)
         return failure(err, nerr,
                        "NRowFull must be an exact positive integer");
   }
   if (!check_integer_metadata(a->xm, "GoodList", a->xm->good, a->xm->nrow, err,
                               nerr) ||
       !check_integer_metadata(a->xm, "RunStart", a->xm->runstart, a->xm->nrun,
                               err, nerr) ||
       (!(a->opt.pen || a->opt.nuis) &&
        !check_integer_metadata(a->xm, "ColumnGroups", NULL, a->xm->ncol, err,
                                nerr)))
      return 0;
   {
      const char *labels = THD_xmat_attr(a->xm, "ColumnLabels"),
                 *tr = THD_xmat_attr(a->xm, "RowTR");
      NI_str_array *list =
          labels ? NI_decode_string_list((char *)labels, ";") : NULL;
      int ok = !labels || (list && list->num == a->xm->ncol);
      double value;
      if (list)
         NI_delete_str_array(list);
      if (!ok)
         return failure(err, nerr, "ColumnLabels count does not match matrix");
      if (tr && (!parse_finite_double(tr, &value) || value <= 0))
         return failure(err, nerr, "RowTR must be finite and positive");
   }
   for (j = 0; j < a->xm->ncol; j++) {
      if (!a->xm->collab[j] || !*a->xm->collab[j])
         return failure(err, nerr, "empty column label");
      for (i = 0; i < j; i++)
         if (!strcmp(a->xm->collab[i], a->xm->collab[j]))
            return failure(err, nerr, "duplicate column label");
   }
   a->roles = calloc((size_t)a->xm->ncol, sizeof(int));
   a->runs = calloc((size_t)a->xm->nrow, sizeof(int));
   if (!a->roles || !a->runs)
      return failure(err, nerr, "cannot allocate design maps");
   for (i = 0; i < a->xm->nrow; ++i)
      a->runs[i] = THD_xmat_run_of(a->xm, a->xm->good[i]);
   if (!prepare_prediction_folds(a, err, nerr))
      return 0;
   if (!assign_column_roles(a, err, nerr) || !prepare_groups(a, err, nerr))
      return 0;
   /* Reject incompatible shapes before allocating any candidate factors.
      The collection budget was calculated using the reference dimensions. */
   if (reference &&
       (a->xm->nrow != reference->xm->nrow ||
        a->xm->ncol != reference->xm->ncol ||
        a->xm->nfull != reference->xm->nfull ||
        a->xm->nrun != reference->xm->nrun ||
        memcmp(a->roles, reference->roles, (size_t)a->xm->ncol * sizeof(int))))
      return failure(err, nerr, "candidate dimensions/column roles differ");
   if (!open_input_headers(a, err, nerr) || !check_memory_budget(a, err, nerr))
      return 0;
   if ((size_t)a->xm->nrow > SIZE_MAX / sizeof(double) / (size_t)a->xm->ncol)
      return failure(err, nerr, "matrix too large");
   a->matrix = malloc(sizeof(double) * (size_t)a->xm->nrow * a->xm->ncol);
   if (!a->matrix)
      return failure(err, nerr, "cannot allocate numeric matrix");
   /* The xmat reader stores one whole column at a time. Numerical plans use
      one whole row at a time: offset = row * number_of_columns + column.
      This is a layout conversion, not a transpose of the statistical model.
      Promotion to double cannot restore precision lost by the text reader. */
   for (i = 0; i < a->xm->nrow; i++) {
      a->runs[i] = THD_xmat_run_of(a->xm, a->xm->good[i]);
      for (j = 0; j < a->xm->ncol; j++)
         a->matrix[(size_t)i * a->xm->ncol + j] =
             a->xm->mat[(size_t)j * a->xm->nrow + i];
   }
   {
      THD_banded_config bc = {a->opt.ngroup, a->opt.nband, a->groups,
                              a->opt.band_values};
      THD_picard_config cfg = {a->opt.elastic ? &a->opt.solver : NULL,
                               a->opt.banded ? &bc : NULL, a->opt.ar1,
                               a->opt.support_refit,
                               a->opt.covariance_set ?
                                   &a->opt.covariance : NULL};
      a->plan = THD_picard_plan_create_config(
          a->matrix, a->xm->nrow, a->xm->ncol, a->roles, a->xm->good, a->runs,
          a->xm->nfull, a->opt.tol, &cfg, err, nerr);
   }
   /* The plan copied its inputs. Release the temporary conversion now,
      and clear its pointer so later App cleanup cannot free it twice. */
   free(a->matrix);
   a->matrix = NULL;
   if (!a->plan)
      return 0;
   if (a->opt.single &&
       (!read_trial_info(a, err, nerr) ||
        !THD_picard_single_trials(a->plan, a->trials, err, nerr)))
      return 0;
   if (a->opt.predict) {
      a->global_index = -1;
      a->predict = THD_picard_predict_plan_create(
          a->plan, a->folds, a->opt.grid, a->opt.ngrid, a->opt.project,
          a->opt.one_se, err, nerr);
      if (!a->predict)
         return 0;
   } else if (a->opt.ngrid) {
      if (!a->trials && !read_trial_info(a, err, nerr))
         return 0;
      a->repeat = THD_picard_repeat_plan_create(a->plan, a->trials, a->opt.grid,
                                                a->opt.ngrid, err, nerr);
      if (!a->repeat)
         return 0;
   }
   return 1;
}

/* Read the collection before opening its first matrix, so the aggregate
   memory check accounts for every candidate before factorization begins. */
static int prepare_candidate_manifest(App *a, char *err, int nerr)
{
   if (!a->opt.candidate_manifest)
      return 1;
   a->manifest = picard_manifest_read(a->opt.candidate_manifest,
                                      a->opt.hrf_policy, err, nerr);
   if (!a->manifest)
      return 0;
   a->candidate_designs =
       calloc((size_t)a->manifest->count, sizeof(*a->candidate_designs));
   if (!a->candidate_designs)
      return failure(err, nerr, "cannot allocate candidate design references");
   a->opt.matrix = a->manifest->entry[0].matrix;
   return 1;
}

/* Each additional matrix passes the same strict xmat/role/numeric checks as a
   standalone fit. The temporary App owns only that matrix's preparation.
   Column identities are validated before any input response is loaded. */
static int prepare_candidate_designs(App *a, char *err, int nerr)
{
   THD_picard_candidates_config config = {0};
   int c, ok = 1, *ids;
   if (!a->manifest)
      return 1;
   a->candidate_designs[0] =
       a->plan; /* Borrowed; free_app owns it separately. */
   if (!picard_manifest_read_columns(a->manifest, 0, a->xm, a->roles, a->trials,
                                     err, nerr))
      return 0;
   ids = calloc((size_t)a->manifest->count, sizeof(int));
   if (!ids)
      return failure(err, nerr, "cannot allocate candidate IDs");
   ids[0] = a->manifest->entry[0].id;
   for (c = 1; c < a->manifest->count && ok; c++) {
      App candidate = {0};
      /* Copy settings, then clear owned/irrelevant fields before cleanup
         can see them. This temporary App prepares a design, not another
         dataset load or a second copy of the tuning grid. */
      candidate.opt = a->opt;
      candidate.opt.matrix = a->manifest->entry[c].matrix;
      candidate.opt.input = candidate.opt.mask = NULL;
      candidate.opt.single = 0;
      candidate.opt.grid = NULL;
      candidate.opt.ngrid = candidate.opt.predict = 0;
      candidate.opt.trial_info = NULL;
      candidate.opt.validate = 1;
      ok = prepare_design(&candidate, a, err, nerr);
      ok =
          ok &&
          picard_manifest_read_columns(a->manifest, c, candidate.xm,
                                       candidate.roles, a->trials, err, nerr) &&
          picard_manifest_compare(a->manifest, c, a->xm, candidate.xm, err,
                                  nerr);
      if (ok) {
         /* Transfer ownership: clearing the source pointer prevents the
            temporary App's destructor from freeing our retained plan. */
         a->candidate_designs[c] = candidate.plan;
         candidate.plan = NULL;
         ids[c] = a->manifest->entry[c].id;
      }
      /* The settings copy borrowed the outer penalty table. */
      candidate.opt.band_values = NULL;
      candidate.opt.pc_counts = NULL;
      free_app(&candidate);
   }
   config.policy =
       a->opt.hrf_policy ? THD_PICARD_SELECT_HRF : THD_PICARD_SELECT_PREDICT;
   config.nparameter = a->opt.ngrid;
   config.fixed_parameter = a->opt.fraction;
   if (a->opt.validate && !isfinite(config.fixed_parameter))
      config.fixed_parameter =
          1; /* No response is fitted in validation mode. */
   config.parameters = a->opt.grid;
   config.trials = a->trials;
   config.folds = a->folds;
   config.nuisance_mode = a->opt.project;
   config.one_se = a->opt.one_se;
   if (ok) {
      a->candidate_plan = THD_picard_candidates_create(
          (const THD_picard_plan *const *)a->candidate_designs, ids,
          a->manifest->count, &config, err, nerr);
      ok = a->candidate_plan != NULL;
   }
   free(ids);
   return ok && check_memory_budget(a, err, nerr);
}

/* Three selection maps follow all ordinary fit/tuning maps, preserving their
   existing offsets. Optional per-design MSE maps follow these three maps. */
static int design_diagnostic_start(const App *a)
{
   return diagnostic_count(a) - single_extra_count(a) - NDESIGN_DIAG -
          (a->opt.save_candidates ? a->manifest->count : 0);
}

/* Treat any existing path, including a dangling symlink, as a collision. */
static int path_exists(const char *name)
{
   struct stat st;
   return lstat(name, &st) == 0;
}

/* Allocate a derived prefix while preserving an explicitly requested NIfTI
   suffix for datasets and removing it for optional text manifests. */
static char *companion_name(const char *prefix, const char *suffix,
                            int keep_format)
{
   size_t n = strlen(prefix), base = n;
   const char *ext = "";
   char *name;
   if (n >= 7 && !strcmp(prefix + n - 7, ".nii.gz")) {
      base = n - 7;
      ext = ".nii.gz";
   } else if (n >= 4 && !strcmp(prefix + n - 4, ".nii")) {
      base = n - 4;
      ext = ".nii";
   }
   name = malloc(base + strlen(suffix) + strlen(ext) + 1);
   if (name)
      snprintf(name, base + strlen(suffix) + strlen(ext) + 1, "%.*s%s%s",
               (int)base, prefix, suffix, keep_format ? ext : "");
   return name;
}

/* Check both the header and common compressed brick paths before writing. */
static int dataset_exists(THD_3dim_dataset *d)
{
   const char *extensions[] = {"", ".gz", ".bz2", ".Z"};
   char *path;
   int i, found;
   if (path_exists(DSET_HEADNAME(d)))
      return 1;
   path = malloc(strlen(DSET_BRIKNAME(d)) + 5);
   if (!path)
      return 1;
   found = 0;
   for (i = 0; i < 4; i++) {
      sprintf(path, "%s%s", DSET_BRIKNAME(d), extensions[i]);
      if (path_exists(path))
         found = 1;
   }
   free(path);
   return found;
}

/* Remove files belonging to a newly created dataset during rollback only.
   Existing outputs are rejected before this function can be reached. */
static void remove_dataset(THD_3dim_dataset *d)
{
   const char *extensions[] = {"", ".gz", ".bz2", ".Z"};
   char *path;
   int i;
   remove(DSET_HEADNAME(d));
   path = malloc(strlen(DSET_BRIKNAME(d)) + 5);
   if (!path)
      return;
   for (i = 0; i < 4; i++) {
      sprintf(path, "%s%s", DSET_BRIKNAME(d), extensions[i]);
      remove(path);
   }
   free(path);
}

/* Record exact train/test membership, including gap exclusions. This is
   selection-CV provenance, not a claim of independent outer evaluation. */
static int attach_prediction_provenance(App *a, THD_3dim_dataset *out,
                                        char *err, int nerr)
{
   char *grid = NULL, *folds = NULL;
   int f, c, j;
   if (!a->predict)
      return 1;
   for (c = 0; c < a->opt.ngrid; ++c)
      grid = THD_zzprintf(grid, "%s%.17g", c ? "," : "", a->opt.grid[c]);
   for (f = 0; f < a->folds->nfold; ++f) {
      folds =
          THD_zzprintf(folds, "%s%s train:", f ? ";" : "", a->folds->name[f]);
      for (j = 0; j < a->folds->ntrain[f]; ++j)
         folds = THD_zzprintf(folds, "%s%d", j ? "," : "",
                              a->xm->good[a->folds->train[f][j]]);
      folds = THD_zzprintf(folds, " test:");
      for (j = 0; j < a->folds->ntest[f]; ++j)
         folds = THD_zzprintf(folds, "%s%d", j ? "," : "",
                              a->xm->good[a->folds->test[f][j]]);
   }
   if (!grid || !folds) {
      free(grid);
      free(folds);
      return failure(err, nerr, "cannot allocate CV provenance");
   }
   THD_set_string_atr(out->dblk,
                      a->opt.banded    ? "PICARD_BAND_INDICES"
                      : a->opt.elastic ? "PICARD_LAMBDAS"
                                       : "PICARD_FRACTIONS",
                      grid);
   THD_set_string_atr(out->dblk, "PICARD_PRED_FOLDS", folds);
   THD_set_string_atr(out->dblk, "PICARD_PRED_NUISANCE",
                      a->opt.project ? "project" : "train");
   THD_set_string_atr(out->dblk, "PICARD_PRED_RULE",
                      a->opt.one_se ? "1se" : "min");
   THD_set_string_atr(out->dblk, "PICARD_PRED_SCOPE",
                      a->opt.global ? "global_mask" : "per_voxel");
   THD_set_string_atr(
       out->dblk, "PICARD_TUNING",
       "prediction-v1; shared columns; training-only fits/scales/penalty "
       "preparation; "
       "train nuisance predicts from training gamma; project nuisance scores "
       "P(Ztest)Ytest against P(Ztest)Xtest beta (conditional signal). "
       "Baseline is training Y mean through the same test projection. "
       "MSE pools squared errors by test-row counts; R2=1-SSE/baselineSSE. "
       "Negative R2 retained; zero baseline: R2=0, PredStatus=1; normal=0; "
       "outside mask=-1. Min ties choose weakest penalty. 1se chooses "
       "strongest "
       "penalty within best weighted mean + weighted-fold SE. "
       "SE^2=sum(w*(foldMSE-mean)^2)*sum(w^2)/(1-sum(w^2)), w=ntest/n. "
       "Global: explicit mask, exclude zero-baseline voxels, normalize each "
       "voxel's fold losses by its pooled baseline MSE, then average voxels. "
       "No valid global voxels: weakest penalty. These are selection scores, "
       "not independent outer-test performance. Final full-design refit. "
       "Folds record ORIGINAL full-row IDs; gaps never use compressed time.");
   free(grid);
   free(folds);
   return 1;
}

/* Store a self-contained ID/path/policy and column identity record in every
   output. A voxelwise design choice cannot be represented by one Xmat path. */
static int attach_candidate_provenance(App *a, THD_3dim_dataset *out, char *err,
                                       int nerr)
{
   char *catalog = NULL, *columns = NULL;
   const char *policy;
   int c, j, start;
   if (!a->manifest)
      return 1;
   for (c = 0; c < a->manifest->count; c++) {
      const PicardCandidateInput *e = a->manifest->entry + c;
      catalog = THD_zzprintf(
          catalog, "%s%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d",
          c ? "\n"
            : "index\tid\tmatrix\tcolumns\tsource\tversion\tsampling\t"
              "duration\tnormalization\trank\n",
          c, e->id, e->matrix, e->columns, e->source, e->version, e->sampling,
          e->duration, e->normalization,
          THD_picard_plan_info(a->candidate_designs[c])->xrank);
   }
   for (j = 0; j < a->xm->ncol; j++) {
      const PicardColumnIdentity *ci = a->manifest->entry[0].identity + j;
      columns = THD_zzprintf(
          columns, "%s%d\t%s\t%s\t%d\t%d\t%d\t%d\t%s",
          j ? "\n"
            : "column\tlabel\trole\ttrial\tcondition\trun\tbasis\titem\n",
          j, ci->label, ci->role ? "penalized" : "nuisance", ci->trial,
          ci->condition, ci->run, ci->basis, ci->item);
   }
   if (!catalog || !columns) {
      free(catalog);
      free(columns);
      return failure(err, nerr, "cannot allocate candidate provenance");
   }
   policy =
       a->opt.hrf_policy
           ? "candidates-v1 HRF: minimum full-data unregularized residual MSE; "
             "equal projected rank; common nuisance values; no CV for HRF "
             "choice. "
             "Then fixed-parameter fit or repeated-condition tuning on the "
             "chosen "
             "design. Subsequent beta CV is CONDITIONAL on full-data HRF "
             "selection. "
             "OLS HRF MSE is a fit-selection score, not held-out performance."
           : "candidates-v1 prediction: tune the parameter separately within "
             "each "
             "design using the same folds/nuisance/selection rule; then choose "
             "the "
             "design with smallest selected prediction MSE. One-SE applies "
             "only "
             "within designs. Validation data are reused to select design and "
             "parameter; neither MSE nor R2 is independent outer-test "
             "performance.";
   THD_set_string_atr(out->dblk, "PICARD_CANDIDATES", catalog);
   THD_set_string_atr(out->dblk, "PICARD_CANDIDATE_COLUMNS", columns);
   THD_set_string_atr(out->dblk, "PICARD_CANDIDATE_MANIFEST",
                      a->manifest->path);
   THD_set_string_atr(out->dblk, "PICARD_CANDIDATE_POLICY",
                      a->opt.candidate_policy);
   THD_set_string_atr(out->dblk, "PICARD_CANDIDATE_SELECTION", (char *)policy);
   THD_set_string_atr(
       out->dblk, "PICARD_CANDIDATE_MAPS",
       "DesignID=stable manifest ID; DesignIndex=zero-based manifest row; "
       "DesignMSE=policy-specific selection loss; exact ties choose smallest "
       "stable ID; outside mask=-1. Optional DesignMSE_ID maps follow manifest "
       "order. Coefficient columns share the recorded identities for ALL "
       "designs. Rank in PICARD_DIMENSIONS describes the first candidate; "
       "the Picard_Rank map and catalog describe actual selected/candidate "
       "ranks.");
   if (out == a->out[OUT_MAIN]) {
      start = THD_picard_plan_info(a->plan)->p +
              (a->opt.save_nuis ? THD_picard_plan_info(a->plan)->k : 0) +
              design_diagnostic_start(a);
      THD_set_int_atr(out->dblk, "PICARD_DESIGN_DIAGNOSTIC_START", 1, &start);
   }
   tross_Append_History(out, (char *)policy);
   free(catalog);
   free(columns);
   return 1;
}

/* Record numerical policy and identity mappings in each dataset header.
   Header metadata keeps the default workflow self-contained in AFNI buckets. */
static int attach_provenance(App *a, THD_3dim_dataset *out, char *err, int nerr)
{
   const THD_picard_design_info *d = THD_picard_plan_info(a->plan);
   char *text = NULL, *labels = NULL;
   int j, diag = d->p + saved_nuisance_count(a),
          dims[] = {d->n, d->nfull, d->p, d->k, d->xrank, d->zrank};
   text = THD_zzprintf(
       NULL,
       "3dPicard phase3a-v1: mode=%s; fixed_fraction=%.17g; rank_tol=%.17g; "
       "fraction_tol=%.17g; "
       "batch_size=%d; jobs=%d. "
       "Preserve supplied design units; no added intercept. Project X/Y "
       "against supplied Z, "
       "fit penalized coefficients, recover minimum-norm nuisance coefficients "
       "from Y-XB. "
       "Reference: minimum-norm LS in retained SVD span s>rank_tol*smax. "
       "Lambda=PenaltyXtX/n in SSE/(2n)+lambda*||beta||^2/2. "
       "Zero reference: beta=0, lambda=0, achieved fraction undefined (-1). "
       "Fraction zero: infinite penalty (-1). Status "
       "0=normal,1=zero_reference,2=zero_endpoint,-1=outside_mask. "
       "Projection cancellation cutoff: 64*DBL_EPSILON*max(n,k)*input_norm. "
       "SVD: existing EISPACK; bounded bisection <=256 iterations. "
       "xmat parsed as float; numerical work double; output float. "
       "Fitts/errts contain retained rows ONLY; original rows/runs in "
       "PICARD_ROW_IDS/PICARD_RUN_IDS. "
       "Dimensions attribute: "
       "retained_rows,full_rows,penalized_columns,nuisance_columns,projected_"
       "rank,nuisance_rank.",
       a->predict  ? "prediction"
       : a->repeat ? "repeated_condition"
                   : "fixed",
       a->opt.ngrid ? -1 : a->opt.fraction, a->opt.tol, THD_FRACTION_TOL,
       a->opt.batch, a->opt.jobs);
   if (a->opt.elastic) {
      free(text);
      text = THD_zzprintf(
          NULL,
          "3dPicard phase4-v1: penalty=%s; fixed_lambda=%.17g; rho=%.17g; "
          "standardize=%d; KKT_tol=%.17g; max_iterations=%d; rank_tol=%.17g. "
          "Objective SSE/(2n)+lambda*(rho*L1(theta)+(1-rho)*L2(theta)^2/2). "
          "Project X/Y against supplied unpenalized Z first. No implicit "
          "centering/intercept or Y scaling. Optional projected-X RMS scales "
          "use training rows only; zero columns keep scale 1. Export beta "
          "in supplied units and recover gamma from original Y-X*beta. "
          "LambdaMax=max(abs(Xeffective'Yprojected))/(n*rho), undefined at "
          "rho=0 (-1). Full KKT check each sweep; failure aborts output. "
          "Warm paths within fold/batch only; full refits cold. "
          "Status=0 converged,-1 outside mask. Objective/KKT use effective "
          "optimization units; Rank uses original projected X and rank_tol. "
          "xmat parsed as float; numerical work double; output float. "
          "Fitts/errts contain retained rows only; original IDs in headers. "
          "Dimensions: retained_rows,full_rows,penalized_columns,"
          "nuisance_columns,projected_rank,nuisance_rank.",
          a->opt.penalty, a->opt.ngrid ? -1 : a->opt.lambda, a->opt.solver.rho,
          a->opt.solver.standardize, a->opt.solver.tolerance,
          a->opt.solver.max_iterations, a->opt.tol);
   }
   if (a->opt.banded) {
      char *groups = NULL, *penalties = NULL;
      int c, g;
      free(text);
      text = THD_zzprintf(
          NULL,
          "3dPicard phase6-v1 banded: SSE/(2n)+sum_j "
          "lambda[group[j]]*beta[j]^2/2. "
          "No implicit centering/scaling. Nuisance is unpenalized. "
          "XtX penalties=n_train*lambda. Positive group penalties only. "
          "Candidate index is zero-based table row. Ties prefer smaller mean "
          "group lambda, then earlier table row. One-SE prefers larger mean "
          "group lambda; vector ordering is a convention, not a partial-order "
          "claim. Status=0 success,-1 outside mask. Rank=projected design "
          "rank. "
          "Double computation; float xmat input/output maps.");
      for (j = 0; j < d->p; j++)
         groups = THD_zzprintf(groups, "%s%d:%d", j ? "," : "",
                               d->penalized_columns[j], a->groups[j]);
      for (c = 0; c < a->opt.nband; c++)
         for (g = 0; g < a->opt.ngroup; g++)
            penalties =
                THD_zzprintf(penalties, "%s%.17g",
                             g   ? " "
                             : c ? ";"
                                 : "",
                             a->opt.band_values[(size_t)c * a->opt.ngroup + g]);
      if (!groups || !penalties) {
         free(groups);
         free(penalties);
         free(text);
         return failure(err, nerr, "cannot allocate banded provenance");
      }
      THD_set_string_atr(out->dblk, "PICARD_BAND_GROUPS", groups);
      THD_set_string_atr(out->dblk, "PICARD_BAND_PENALTIES", penalties);
      free(groups);
      free(penalties);
   }
   text = THD_zzprintf(
       text,
       " LegacyFixedAR1=%.17g (inactive with PICARD_COVARIANCE): "
       "whiten before nuisance projection; restart at run "
       "changes and nonconsecutive original rows, separately inside folds. "
       "Segment starts multiply by sqrt(1-rho^2); other rows use y-rho*y_prev. "
       "No noise estimation. Fitts/errts and selection losses remain in "
       "original units. SupportRefit=%d: exact nonzero sparse support, "
       "minimum-norm OLS on that support in the same effective fitting space. "
       "Support learned anew in each training fold and full refit. "
       "Sparse Objective/KKT/Nonzero describe PRE-refit selection; "
       "SupportSize/SupportRank describe the refit. No post-selection "
       "inference.",
       a->opt.ar1, a->opt.support_refit);
   if (a->opt.covariance_set) {
      char *cov = THD_zzprintf(NULL,
          "model=%s scalar_a=%.17g scalar_b=%.17g gaps=%s "
          "normalization=unit_marginal "
          "MA_sign=plus factor=exact_dense_cholesky "
          "map=%s source=%s selection=conditional_on_covariance_parameters "
          "voxel_parameters_override_scalars=1 "
          "parameter_bricks=CovarianceAR,CovarianceMA",
          a->opt.covariance.model == 0 ? "white" :
          a->opt.covariance.model == 1 ? "ar1" : "arma11",
          a->opt.covariance.a, a->opt.covariance.b,
          a->opt.covariance.reset_gaps ? "reset" : "preserve",
          a->opt.covariance_map ? a->opt.covariance_map : "none",
          a->opt.covariance_source ? a->opt.covariance_source :
          a->opt.covariance_map ? "unknown_external" :
          a->opt.covariance_estimate ? "native_pilot_REML" : "fixed_scalar");
      if (!cov) {
         free(text);
         return failure(err, nerr, "cannot allocate covariance metadata");
      }
      THD_set_string_atr(out->dblk, "PICARD_COVARIANCE", cov);
      free(cov);
   }
   for (j = 0; j < d->ncol; j++)
      labels = THD_zzprintf(labels, "%s%s", j ? " ; " : "", a->xm->collab[j]);
   if (!text || !labels) {
      free(text);
      free(labels);
      return failure(err, nerr, "cannot allocate header provenance");
   }
   if (a->repeat) {
      const THD_picard_repeat_info *info =
          THD_picard_repeat_plan_info(a->repeat);
      char *grid = NULL, *trials = NULL, *folds = NULL;
      for (j = 0; j < info->ncandidate; ++j)
         grid = THD_zzprintf(grid, "%s%.17g", j ? "," : "", info->fractions[j]);
      for (j = 0; j < d->p; ++j) {
         const THD_picard_trial *t = info->trials + j;
         trials = THD_zzprintf(trials, "%s%d %d %d %d", j ? ";" : "", t->column,
                               t->trial, t->condition, t->run);
      }
      for (j = 0; j < info->nfold; ++j)
         folds = THD_zzprintf(folds, "%s%d %.0f", j ? ";" : "",
                              info->fold_runs[j], info->fold_pairs[j]);
      if (!grid || !trials || !folds) {
         free(grid);
         free(trials);
         free(folds);
         free(text);
         free(labels);
         return failure(err, nerr, "cannot allocate tuning provenance");
      }
      THD_set_string_atr(out->dblk,
                         a->opt.banded    ? "PICARD_BAND_INDICES"
                         : a->opt.elastic ? "PICARD_LAMBDAS"
                                          : "PICARD_FRACTIONS",
                         grid);
      THD_set_string_atr(out->dblk, "PICARD_TRIAL_INFO", trials);
      THD_set_string_atr(out->dblk, "PICARD_TRIAL_SOURCE", a->opt.trial_info);
      THD_set_string_atr(out->dblk, "PICARD_CV_FOLDS", folds);
      THD_set_string_atr(
          out->dblk, "PICARD_TUNING",
          "repeat-v2: leave-one-run-out; fold-local nuisance, scaling and "
          "penalty preparation; "
          "all ordered cross-run same-condition trial pairs against held-out "
          "OLS; "
          "pooled MSE in coefficient units, pair-weighted; no session "
          "normalization; "
          "exact ties choose weakest penalty; final full-design refit; "
          "no autoscaling. CV losses are selection losses, not outer-test "
          "scores. "
          "TRIAL_INFO rows: column trial condition run (zero-based). "
          "CV_FOLDS rows: held-out run, pair count. CV_Index indexes "
          "FRACTIONS (ridge) or LAMBDAS (lasso/enet). "
          "CV maps outside mask=-1. This pair-scoring workflow is not "
          "full GLMsingle equivalence.");
      free(grid);
      free(trials);
      free(folds);
   }
   THD_set_string_atr(out->dblk, "PICARD_MODEL", text);
   THD_set_string_atr(out->dblk, "PICARD_ESTIMATOR", a->opt.penalty);
   THD_set_string_atr(out->dblk, "PICARD_MATRIX",
                      a->manifest ? "MULTIPLE: see PICARD_CANDIDATES"
                                  : a->opt.matrix);
   THD_set_string_atr(out->dblk, "PICARD_INPUT", a->opt.input);
   THD_set_string_atr(out->dblk, "PICARD_MASK",
                      a->opt.mask ? a->opt.mask : "NONE");
   if (a->xm->tr > 0)
      THD_set_float_atr(out->dblk, "PICARD_ROW_TR", 1, &a->xm->tr);
   {
      const char *keys[] = {"StimLabels", "StimBots", "StimTops"};
      char *names[] = {"PICARD_STIM_LABELS", "PICARD_STIM_BOTS",
                       "PICARD_STIM_TOPS"};
      for (j = 0; j < 3; j++) {
         const char *value = THD_xmat_attr(a->xm, keys[j]);
         if (value)
            THD_set_string_atr(out->dblk, names[j], (char *)value);
      }
   }
   THD_set_string_atr(out->dblk, "PICARD_COLUMN_LABELS", labels);
   THD_set_int_atr(out->dblk, "PICARD_DIMENSIONS", 6, dims);
   THD_set_int_atr(out->dblk, "PICARD_ROW_IDS", d->n, (int *)d->row_ids);
   THD_set_int_atr(out->dblk, "PICARD_RUN_IDS", d->n, (int *)d->run_ids);
   THD_set_int_atr(out->dblk, "PICARD_RUN_START", a->xm->nrun, a->xm->runstart);
   THD_set_int_atr(out->dblk, "PICARD_COLUMN_ROLES", d->ncol, a->roles);
   THD_set_int_atr(out->dblk, "PICARD_PENALIZED_COLUMNS", d->p,
                   (int *)d->penalized_columns);
   if (d->k)
      THD_set_int_atr(out->dblk, "PICARD_NUISANCE_COLUMNS", d->k,
                      (int *)d->nuisance_columns);
   if (out == a->out[OUT_MAIN])
      THD_set_int_atr(out->dblk, "PICARD_DIAGNOSTIC_START", 1, &diag);
   THD_set_string_atr(out->dblk, "PICARD_DATASET_ROLE",
                      out == a->out[OUT_MAIN] ? "coefficients_and_diagnostics"
                      : out == a->out[OUT_FITTED] ? "fitted_retained_rows"
                                                  : "residual_retained_rows");
   if (!attach_prediction_provenance(a, out, err, nerr) ||
       !attach_candidate_provenance(a, out, err, nerr)) {
      free(text);
      free(labels);
      return 0;
   }
   tross_Append_History(out, text);
   free(text);
   free(labels);
   return 1;
}

/* Allocate float buckets and assign labels before starting worker threads.
   Fitted/residual buckets represent retained rows, not a regular-TR series.
   No files are written until every voxel has fitted successfully. */
static int prepare_outputs(App *a, int argc, char **argv, char *err, int nerr)
{
   const THD_picard_design_info *d = THD_picard_plan_info(a->plan);
   int nk = saved_nuisance_count(a);
   int counts[] = {d->p + nk + diagnostic_count(a) +
                      (a->opt.covariance_estimate ? 6 : 2*a->opt.covariance_set),
                   a->opt.save_fit ? d->n : 0,
                   a->opt.save_res ? d->n : 0};
   int t, b, j;
   if (a->opt.save_nuis && !nk)
      return failure(
          err, nerr,
          "save_nuisance requested but design has no nuisance columns");
   a->names[OUT_MAIN] = strdup(a->opt.prefix);
   if (a->opt.save_fit)
      a->names[OUT_FITTED] = a->opt.fitts
                                 ? strdup(a->opt.fitts)
                                 : companion_name(a->opt.prefix, "_fitts", 1);
   if (a->opt.save_res)
      a->names[OUT_RESIDUAL] = a->opt.errts
                                   ? strdup(a->opt.errts)
                                   : companion_name(a->opt.prefix, "_errts", 1);
   if (a->opt.manifest) {
      a->names[OUT_BETAS_TSV] = companion_name(a->opt.prefix, ".betas.tsv", 0);
      a->names[OUT_ROWS_1D] = companion_name(a->opt.prefix, ".rows.1D", 0);
   }
   for (t = 0; t < OUTPUT_COUNT; t++)
      if ((t < DATASET_COUNT && counts[t]) ||
          (t >= DATASET_COUNT && a->opt.manifest)) {
         if (!a->names[t])
            return failure(err, nerr, "cannot allocate output name");
         if (t >= DATASET_COUNT && path_exists(a->names[t])) {
            snprintf(err, (size_t)nerr, "output already exists: %s",
                     a->names[t]);
            return 0;
         }
      }
   for (t = 0; t < DATASET_COUNT; t++)
      if (counts[t]) {
         a->out[t] = EDIT_empty_copy(a->input);
         if (!a->out[t])
            return failure(err, nerr, "cannot allocate output dataset");
         /* Explicitly choose native storage unless the caller asks for NIfTI.
          */
         a->out[t]->dblk->diskptr->storage_mode =
             STRING_HAS_SUFFIX(a->names[t], ".nii") ||
                     STRING_HAS_SUFFIX(a->names[t], ".nii.gz")
                 ? STORAGE_BY_NIFTI
                 : STORAGE_BY_BRICK;
         EDIT_dset_items(a->out[t], ADN_prefix, a->names[t], ADN_nvals,
                         counts[t], ADN_ntt, 0, ADN_type, HEAD_FUNC_TYPE,
                         ADN_func_type, FUNC_BUCK_TYPE, ADN_datum_all,
                         MRI_float, ADN_none);
         if (dataset_exists(a->out[t])) {
            snprintf(err, (size_t)nerr, "output already exists: %.*s",
                     nerr > 24 ? nerr - 24 : 0, DSET_HEADNAME(a->out[t]));
            return 0;
         }
         for (j = 0; j < t; j++)
            if (a->out[j] &&
                !strcmp(DSET_HEADNAME(a->out[j]), DSET_HEADNAME(a->out[t])))
               return failure(err, nerr, "output prefixes collide");
         for (b = 0; b < counts[t]; b++) {
            char row[64], *coef_label = NULL;
            const char *label;
            EDIT_substitute_brick(a->out[t], b, MRI_float, NULL);
            if (!DSET_ARRAY(a->out[t], b))
               return failure(err, nerr, "cannot allocate output brick");
            memset(DSET_ARRAY(a->out[t], b), 0,
                   sizeof(float) * (size_t)a->nvox);
            if (t == OUT_MAIN && b < d->p + nk) {
               if (b >= d->p + d->k) {
                  int pc = b - d->p - d->k;
                  coef_label = THD_zzprintf(NULL, "NoisePC_Run%03d_PC%03d_Coef",
                                            pc / a->single->chosen,
                                            pc % a->single->chosen);
               } else {
                  int c = b < d->p ? d->penalized_columns[b]
                                   : d->nuisance_columns[b - d->p];
                  coef_label = THD_zzprintf(NULL, "%s_Coef", a->xm->collab[c]);
               }
               if (!coef_label)
                  return failure(err, nerr,
                                 "cannot allocate coefficient label");
               label = coef_label;
            } else if (t == OUT_MAIN) {
               int diagnostic = b - d->p - nk;
               int extra = a->predict ? NPRED_DIAG : NTUNE_DIAG;
               if (diagnostic >= diagnostic_count(a)) {
                  static const char *native[] = {"Picard_CovarianceSigma2",
                      "Picard_CovarianceREML", "Picard_CovarianceRank",
                      "Picard_CovarianceStatus", "Picard_CovarianceAR",
                      "Picard_CovarianceMA"};
                  label = native[diagnostic-diagnostic_count(a) +
                                 (a->opt.covariance_estimate ? 0 : 4)];
               }
               else if (a->single && diagnostic >= single_start(a)) {
                  static const char *names[] = {
                      "Picard_PCCount", "Picard_NoisePool",
                      "Picard_PCSelection", "Picard_NoiseModelR2"};
                  int index = diagnostic - single_start(a);
                  if (index < 4)
                     label = names[index];
                  else {
                     index -= 4;
                     coef_label = THD_zzprintf(
                         NULL, "Stage%c_%s_Coef", index < d->p ? 'B' : 'C',
                         a->xm->collab[d->penalized_columns[index % d->p]]);
                     if (!coef_label)
                        return failure(err, nerr, "cannot label stage");
                     label = coef_label;
                  }
               } else if (a->manifest &&
                          diagnostic >= design_diagnostic_start(a)) {
                  static const char *names[] = {"Picard_DesignID",
                                                "Picard_DesignIndex",
                                                "Picard_DesignMSE"};
                  int index = diagnostic - design_diagnostic_start(a);
                  if (index < NDESIGN_DIAG)
                     label = names[index];
                  else {
                     snprintf(row, sizeof(row), "Picard_DesignMSE_ID%d",
                              a->manifest->entry[index - NDESIGN_DIAG].id);
                     label = row;
                  }
               } else if (diagnostic < base_diagnostic_count(a)) {
                  if (a->opt.banded) {
                     const char *band_labels[] = {"Picard_BandIndex",
                                                  "Picard_Objective",
                                                  "Picard_KKT", "Picard_Rank"};
                     if (diagnostic < 4)
                        label = band_labels[diagnostic];
                     else if (diagnostic == 4 + a->opt.ngroup)
                        label = "Picard_Status";
                     else {
                        snprintf(row, sizeof(row), "Picard_LambdaGroup%03d",
                                 diagnostic - 4);
                        label = row;
                     }
                  } else if (a->opt.support_refit && diagnostic >= 9) {
                     const char *refit_labels[] = {"Picard_SupportSize",
                                                   "Picard_SupportRank",
                                                   "Picard_Status"};
                     label = refit_labels[diagnostic - 9];
                  } else
                     label = a->opt.elastic ? elastic_labels[diagnostic]
                                            : diag_labels[diagnostic];
               } else if (a->predict &&
                          diagnostic < base_diagnostic_count(a) + NPRED_DIAG)
                  label =
                      prediction_labels[diagnostic - base_diagnostic_count(a)];
               else if (diagnostic < base_diagnostic_count(a) + extra)
                  label = diag_labels[NDIAG + diagnostic -
                                      base_diagnostic_count(a)];
               else {
                  snprintf(row, sizeof(row), "Picard_CV_Grid%03d",
                           diagnostic - base_diagnostic_count(a) - extra);
                  label = row;
               }
            } else {
               snprintf(row, sizeof(row), "t%06d", a->xm->good[b]);
               label = row;
            }
            EDIT_BRICK_LABEL(a->out[t], b, (char *)label);
            free(coef_label);
            EDIT_BRICK_FACTOR(a->out[t], b, 0);
            EDIT_BRICK_TO_NOSTAT(a->out[t], b);
         }
         tross_Copy_History(a->input, a->out[t]);
         tross_Make_History(PROGRAM, argc, argv, a->out[t]);
         if (!attach_provenance(a, a->out[t], err, nerr))
            return 0;
      }
   /* Coefficient maps start at zero. Diagnostic -1 marks voxels that were
      never fitted, so a masked-out voxel cannot look like a successful fit
      whose coefficients happen to be zero. fit_block replaces these marks. */
   for (b = 0; b < a->nvox; b++) {
      if (a->opt.covariance_set)
         for (j = 0; j < (a->opt.covariance_estimate ? 6 : 2); j++)
            ((float *)DSET_ARRAY(a->out[OUT_MAIN],
                 d->p + nk + diagnostic_count(a) + j))[b] = -1;
      ((float *)DSET_ARRAY(a->out[OUT_MAIN],
                           d->p + nk + base_diagnostic_count(a) - 1))[b] = -1;
      if (a->opt.ngrid || a->manifest || a->single)
         for (j = base_diagnostic_count(a);
              j < (a->single ? single_start(a) + 4 : diagnostic_count(a)); ++j)
            ((float *)DSET_ARRAY(a->out[OUT_MAIN], d->p + nk + j))[b] = -1;
   }
   return 1;
}

/* Load every format through AFNI's normal dataset interface, including
   NIfTI and native HEAD/BRIK. Picard deliberately inherits AFNI's datatype
   conversion and nonfinite-value handling, including environment settings.
   Keep loader fixes in AFNI itself so other programs receive the same fix.

   Opening a header does not allocate its voxel arrays. After DSET_load, check
   each array before workers can dereference it. App owns this loaded storage;
   DSET_delete in free_app releases it on both success and failure. */
static int load_source(THD_3dim_dataset *d, char *err, int nerr)
{
   int t;
   DSET_load(d);
   if (!DSET_LOADED(d))
      return failure(err, nerr, "cannot load input dataset");
   for (t = 0; t < DSET_NVALS(d); t++) {
      if (!DSET_ARRAY(d, t))
         return failure(err, nerr, "cannot load input dataset brick");
      if (!is_real_datum(DSET_BRICK_TYPE(d, t)))
         return failure(err, nerr, "loaded input must contain real bricks");
   }
   return 1;
}

/* Build a compact voxel-index list from the mask, then load the response
   payload. Test mask values after AFNI loading: a repaired zero excludes a
   voxel, any other finite nonzero value includes it. Empty masks are errors.
   The compact list saves fitting work without changing spatial output order. */
static int load_inputs(App *a, char *err, int nerr)
{
   int i;
   if (a->mask && !load_source(a->mask, err, nerr))
      return 0;
   a->vox = calloc((size_t)a->nvox, sizeof(int));
   if (!a->vox)
      return failure(err, nerr, "cannot allocate voxel map");
   for (i = 0; i < a->nvox; i++) {
      double m = a->mask ? read_voxel(a->mask, 0, i) : 1;
      if (!isfinite(m))
         return failure(err, nerr, "mask contains nonfinite values");
      if (m != 0)
         a->vox[a->nselected++] = i;
   }
   if (!a->nselected)
      return failure(err, nerr, "mask is empty");
   return load_source(a->input, err, nerr);
}

/* Open auxiliary masks through AFNI too. Mask ownership ends here; only the
   compact inclusion flags survive. Geometry and nonfinite checks match the
   main mask path, including the loader's normal NaN/Inf handling. */
static unsigned char *single_mask(App *a, const char *path, char *err, int nerr)
{
   THD_3dim_dataset *d = THD_open_dataset((char *)path);
   unsigned char *mask = NULL;
   int i, count = 0, ok = 1;
   if (!d || DSET_NVALS(d) != 1 || !EQUIV_GRIDS(a->input, d) ||
       THD_dataset_mismatch(a->input, d) || DSET_IS_1D(d) || DSET_IS_GIFTI(d) ||
       DSET_IS_NI_SURF_DSET(d) || strchr(path, '['))
      ok = failure(err, nerr,
                   "noise/PC mask must be one unselected volume on input grid");
   if (ok)
      ok = load_source(d, err, nerr);
   if (ok) {
      mask = calloc((size_t)a->nvox, 1);
      ok = mask != NULL;
      if (!ok)
         failure(err, nerr, "cannot allocate noise/PC mask");
   }
   for (i = 0; i < a->nvox && ok; i++) {
      double value = read_voxel(d, 0, i);
      if (!isfinite(value))
         ok = failure(err, nerr, "nonfinite noise/PC mask");
      else
         count += (mask[i] = value != 0);
   }
   if (ok && !count)
      ok = failure(err, nerr, "empty noise/PC mask");
   if (d)
      DSET_delete(d);
   if (!ok) {
      free(mask);
      return NULL;
   }
   return mask;
}
static THD_picard_plan *single_base(App *a, int c)
{
   return a->manifest ? a->candidate_designs[c] : a->plan;
}
/* Gather arbitrary spatial indices, preserving the same retained-row mapping
   as ordinary fit_block. Noise voxels need not be in the analysis mask. */
static void single_gather(App *a, const int *vox, int q, double *y)
{
   int i, o;
   for (i = 0; i < a->xm->nrow; i++)
      for (o = 0; o < q; o++)
         y[(size_t)i * q + o] = read_voxel(a->input, a->xm->good[i], vox[o]);
}
static int single_start(const App *a)
{
   return diagnostic_count(a) - single_extra_count(a);
}

/* A supplied coarse task design replaces GLMsingle's internally generated
   ON/OFF model. Its declared nuisance span defines the R2 denominator. An
   explicit search mask replaces automatic intensity/Gaussian-mixture choices.
   Require identical sampling and no more columns than the trial design. */
static int single_auto_pool(App *a, char *err, int nerr)
{
   App model = {0};
   THD_picard_workspace *work = NULL;
   THD_linproj_workspace *projection = NULL;
   double *y = NULL, *residual = NULL;
   int *vox = NULL;
   int i, o, start, q, count = 0, kept = 0, ok = 1;
   model.opt.matrix = a->opt.noise_model;
   model.opt.penalty = "ridge";
   model.opt.fraction = 1;
   model.opt.tol = a->opt.tol;
   model.opt.batch = a->opt.batch;
   model.opt.jobs = 1;
   model.opt.memory = a->opt.memory;
   model.opt.validate = 1;
   model.xm = THD_xmat_read(model.opt.matrix, err, nerr);
   if (!model.xm)
      ok = 0;
   if (ok && (model.xm->nrow != a->xm->nrow ||
              model.xm->nfull != a->xm->nfull || model.xm->ncol > a->xm->ncol ||
              model.xm->nrun != a->xm->nrun || model.xm->tr != a->xm->tr ||
              memcmp(model.xm->good, a->xm->good,
                     sizeof(int) * (size_t)a->xm->nrow) ||
              memcmp(model.xm->runstart, a->xm->runstart,
                     sizeof(int) * (size_t)a->xm->nrun)))
      ok = failure(err, nerr,
                   "noise_model sampling/shape must match the trial design");
   if (ok)
      ok = prepare_design(&model, NULL, err, nerr);
   if (ok) {
      work = THD_picard_workspace_alloc(model.plan, a->opt.batch, err, nerr);
      projection = THD_linproj_workspace_alloc(model.plan->projection,
                                               a->opt.batch, err, nerr);
      y = calloc((size_t)a->xm->nrow * a->opt.batch, sizeof(double));
      residual = calloc((size_t)a->xm->nrow * a->opt.batch, sizeof(double));
      vox = calloc((size_t)a->nvox, sizeof(int));
      ok = work && projection && y && residual && vox;
   }
   if (ok)
      for (i = 0; i < a->nvox; i++)
         if (a->single->pool[i])
            vox[count++] = i;
   for (start = 0; start < count && ok; start += a->opt.batch) {
      const THD_picard_result *fit;
      q = count - start;
      if (q > a->opt.batch)
         q = a->opt.batch;
      single_gather(a, vox + start, q, y);
      ok = THD_picard_fit_ols(work, y, q, err, nerr) &&
           THD_linproj_project(projection, y, q, residual, err, nerr);
      if (!ok)
         break;
      fit = THD_picard_get_result(work);
      for (o = 0; o < q; o++) {
         long double sse = 0, baseline = 0;
         double r2 = -1;
         int v = vox[start + o];
         for (i = 0; i < a->xm->nrow; i++) {
            long double value = fit->residual[(size_t)i * q + o];
            sse += value * value;
            value = residual[(size_t)i * q + o];
            baseline += value * value;
         }
         if (baseline > 0)
            r2 = (double)(1 - sse / baseline);
         if (!isfinite(r2)) {
            ok = failure(err, nerr, "noise-model R2 exceeds range");
            break;
         }
         a->single->noise_r2[v] = r2;
         a->single->pool[v] = baseline > 0 && r2 < a->opt.noise_r2;
         kept += a->single->pool[v];
      }
   }
   if (ok && !kept)
      ok = failure(err, nerr, "noise R2 threshold selected an empty pool");
   THD_picard_workspace_free(work);
   THD_linproj_workspace_free(projection);
   free(y);
   free(residual);
   free(vox);
   free_app(&model);
   return ok;
}
static int single_compare_double(const void *left, const void *right)
{
   double a = *(const double *)left, b = *(const double *)right;
   return (a > b) - (a < b);
}
/* Pick HRFs once using pre-PC unpenalized fits. Each count then uses that
   fixed choice, a fixed zero-PC beta reference and the same selection voxels.
   Only scalar selection losses are retained across count/voxel combinations;
   coefficient arrays stay bounded by batch size. */
static int single_select(App *a, char *err, int nerr)
{
   Single *s = a->single;
   const THD_picard_design_info *d = THD_picard_plan_info(a->plan);
   THD_picard_workspace **base = calloc((size_t)s->nh, sizeof(*base));
   THD_picard_plan **pcplan = calloc((size_t)s->nh, sizeof(*pcplan));
   THD_picard_workspace **pcwork = calloc((size_t)s->nh, sizeof(*pcwork));
   double *y = calloc((size_t)d->n * a->opt.batch, sizeof(double));
   double *reference = calloc((size_t)d->p * a->opt.batch, sizeof(double));
   double *beta = calloc((size_t)d->p * a->opt.batch, sizeof(double));
   double *loss = calloc((size_t)a->opt.batch, sizeof(double));
   double *scores = calloc((size_t)a->opt.npc * a->nselected, sizeof(double));
   double *median = calloc((size_t)a->nselected, sizeof(double));
   int c, h, i, j, o, start, q,
       ok = base && pcplan && pcwork && y && reference && beta && loss &&
            scores && median;
   if (!ok)
      failure(err, nerr, "cannot allocate staged selection buffers");
   for (h = 0; h < s->nh && ok; h++) {
      ok = THD_picard_single_trials(single_base(a, h), a->trials, err, nerr);
      if (ok)
         base[h] = THD_picard_workspace_alloc(single_base(a, h), a->opt.batch,
                                              err, nerr);
      ok = ok && base[h];
   }
   for (start = 0; start < a->nselected && ok; start += a->opt.batch) {
      q = a->nselected - start;
      if (q > a->opt.batch)
         q = a->opt.batch;
      single_gather(a, a->vox + start, q, y);
      for (h = 0; h < s->nh && ok; h++) {
         const THD_picard_result *fit;
         ok = THD_picard_fit_ols(base[h], y, q, err, nerr);
         if (!ok)
            break;
         fit = THD_picard_get_result(base[h]);
         for (o = 0; o < q; o++) {
            int v = a->vox[start + o];
            long double mse = 0;
            int id = a->manifest ? a->manifest->entry[h].id : 0;
            int prior = a->manifest ? a->manifest->entry[s->hrf[v]].id : 0;
            for (i = 0; i < d->n; i++) {
               long double e = fit->residual[(size_t)i * q + o];
               mse += e * e / d->n;
            }
            if (!isfinite(mse) || mse > DBL_MAX) {
               ok = failure(err, nerr, "HRF MSE exceeds range");
               break;
            }
            if (s->design_loss)
               s->design_loss[(size_t)h * a->nvox + v] = (double)mse;
            if (!h || (double)mse < s->hrf_loss[v] ||
                ((double)mse == s->hrf_loss[v] && id < prior)) {
               s->hrf[v] = h;
               s->hrf_loss[v] = (double)mse;
               if (s->beta_b)
                  for (j = 0; j < d->p; j++)
                     s->beta_b[(size_t)j * a->nvox + v] =
                         fit->beta[(size_t)j * q + o];
            }
         }
      }
   }
   /* With one explicitly fixed count there is no data-driven PC selection. */
   for (c = 0; c < a->opt.npc && a->opt.npc > 1 && ok; c++) {
      for (h = 0; h < s->nh && ok; h++) {
         pcplan[h] = THD_picard_noise_plan(s->noise, single_base(a, h),
                                           a->opt.pc_counts[c], 1, err, nerr);
         if (pcplan[h])
            pcwork[h] =
                THD_picard_workspace_alloc(pcplan[h], a->opt.batch, err, nerr);
         ok = pcwork[h] != NULL;
      }
      for (start = 0; start < a->nselected && ok; start += a->opt.batch) {
         q = a->nselected - start;
         if (q > a->opt.batch)
            q = a->opt.batch;
         single_gather(a, a->vox + start, q, y);
         for (h = 0; h < s->nh && ok; h++) {
            const THD_picard_result *zero, *current;
            ok = THD_picard_fit_ols(base[h], y, q, err, nerr) &&
                 THD_picard_fit_ols(pcwork[h], y, q, err, nerr);
            if (!ok)
               break;
            zero = THD_picard_get_result(base[h]);
            current = THD_picard_get_result(pcwork[h]);
            for (o = 0; o < q; o++)
               if (s->hrf[a->vox[start + o]] == h)
                  for (j = 0; j < d->p; j++) {
                     reference[(size_t)j * q + o] =
                         zero->beta[(size_t)j * q + o];
                     beta[(size_t)j * q + o] = current->beta[(size_t)j * q + o];
                  }
         }
         if (ok)
            ok = THD_picard_single_loss(a->trials, d->p, reference, beta, q,
                                        loss, err, nerr);
         if (ok)
            for (o = 0; o < q; o++)
               scores[(size_t)c * a->nselected + start + o] = loss[o];
      }
      for (h = 0; h < s->nh; h++) {
         THD_picard_workspace_free(pcwork[h]);
         pcwork[h] = NULL;
         THD_picard_plan_free(pcplan[h]);
         pcplan[h] = NULL;
      }
   }
   if (ok && a->opt.npc > 1) {
      /* Eligibility depends only on the fixed zero-PC reference variance.
         Reuse that same voxel set for every count: changing denominators or
         membership with the candidate would change what the trend measures. */
      for (c = 0; c < a->opt.npc; c++) {
         int count = 0;
         for (i = 0; i < a->nselected; i++)
            if (s->select[a->vox[i]] && scores[i] >= 0)
               median[count++] = scores[(size_t)c * a->nselected + i];
         if (!count) {
            ok = failure(err, nerr,
                         "PC mask has no defined reference beta variance");
            break;
         }
         s->nscore = count;
         qsort(median, (size_t)count, sizeof(double), single_compare_double);
         s->trend[c] = count % 2
                           ? median[count / 2]
                           : median[count / 2 - 1] / 2 + median[count / 2] / 2;
      }
      if (ok) {
         c = THD_picard_single_select(s->trend, a->opt.npc, a->opt.pc_stop);
         if (c < 0)
            ok = failure(err, nerr, "invalid PC selection trend");
         else
            s->chosen = a->opt.pc_counts[c];
      }
   } else if (ok)
      s->chosen = a->opt.pc_counts[0];
   for (h = 0; h < s->nh; h++) {
      if (base)
         THD_picard_workspace_free(base[h]);
      if (pcwork)
         THD_picard_workspace_free(pcwork[h]);
      if (pcplan)
         THD_picard_plan_free(pcplan[h]);
   }
   free(base);
   free(pcwork);
   free(pcplan);
   free(y);
   free(reference);
   free(beta);
   free(loss);
   free(scores);
   free(median);
   return ok;
}
/* All learned settings become immutable before the OpenMP voxel-fitting loop.
   The final candidate collection receives locked HRFs, so added PCs cannot
   change their identities. App owns partial setup on every failure path. */
static int single_prepare(App *a, char *err, int nerr)
{
   Single *s;
   THD_picard_candidates_config cfg = {0};
   double *y = NULL;
   int *poolvox = NULL, *ids = NULL;
   int i, c, q, start, npool = 0, ok = 1,
                       maxpc = a->opt.pc_counts[a->opt.npc - 1];
   a->single = s = calloc(1, sizeof(*s));
   if (!s)
      return failure(err, nerr, "cannot allocate staged workflow");
   s->nh = a->manifest ? a->manifest->count : 1;
   s->pool = single_mask(
       a, a->opt.noise_pool ? a->opt.noise_pool : a->opt.noise_search, err,
       nerr);
   if (!s->pool)
      return 0;
   s->noise_r2 = malloc(sizeof(double) * (size_t)a->nvox);
   s->select = calloc((size_t)a->nvox, 1);
   s->hrf = calloc((size_t)a->nvox, sizeof(int));
   s->hrf_loss = calloc((size_t)a->nvox, sizeof(double));
   s->trend = calloc((size_t)a->opt.npc, sizeof(double));
   s->final = calloc((size_t)s->nh, sizeof(*s->final));
   if (a->opt.save_candidates)
      s->design_loss = calloc((size_t)s->nh * a->nvox, sizeof(double));
   if (a->opt.save_stages)
      s->beta_b = calloc((size_t)a->plan->info.p * a->nvox, sizeof(double));
   if (!s->noise_r2 || !s->select || !s->hrf || !s->hrf_loss || !s->trend ||
       !s->final || (a->opt.save_candidates && !s->design_loss) ||
       (a->opt.save_stages && !s->beta_b))
      return failure(err, nerr, "cannot allocate staged maps");
   for (i = 0; i < a->opt.npc; i++)
      s->trend[i] = -1; /* Fixed counts have no measured selection score. */
   for (i = 0; i < a->nvox; i++)
      s->noise_r2[i] = -1;
   if (a->opt.noise_model && !single_auto_pool(a, err, nerr))
      return 0;
   if (a->opt.pc_mask) {
      unsigned char *analysis = calloc((size_t)a->nvox, 1);
      free(s->select);
      s->select = single_mask(a, a->opt.pc_mask, err, nerr);
      if (!analysis || !s->select) {
         free(analysis);
         return 0;
      }
      for (i = 0; i < a->nselected; i++)
         analysis[a->vox[i]] = 1;
      for (i = 0; i < a->nvox; i++)
         if (s->select[i] && (!analysis[i] || s->pool[i]))
            ok = 0;
      free(analysis);
      if (!ok)
         return failure(
             err, nerr,
             "pc_mask must lie inside analysis mask and outside noise pool");
   }
   s->noise = THD_picard_noise_create(a->plan, maxpc, a->opt.batch, err, nerr);
   if (!s->noise)
      return 0;
   s->nrun = THD_picard_noise_runs(s->noise);
   poolvox = calloc((size_t)a->nvox, sizeof(int));
   y = calloc((size_t)a->xm->nrow * a->opt.batch, sizeof(double));
   if (!poolvox || !y)
      ok = failure(err, nerr, "cannot allocate pool batch");
   if (ok)
      for (i = 0; i < a->nvox; i++)
         if (s->pool[i])
            poolvox[npool++] = i;
   for (start = 0; start < npool && ok; start += a->opt.batch) {
      q = npool - start;
      if (q > a->opt.batch)
         q = a->opt.batch;
      single_gather(a, poolvox + start, q, y);
      ok = THD_picard_noise_add(s->noise, y, q, err, nerr);
   }
   free(poolvox);
   free(y);
   if (ok)
      ok = THD_picard_noise_finish(s->noise, err, nerr) &&
           single_select(a, err, nerr);
   if (!ok)
      return 0;
   s->added = s->chosen * s->nrun;
   ids = calloc((size_t)s->nh, sizeof(int));
   if (!ids)
      return failure(err, nerr, "cannot allocate staged HRF IDs");
   for (c = 0; c < s->nh && ok; c++) {
      ids[c] = a->manifest ? a->manifest->entry[c].id : 0;
      s->final[c] = THD_picard_noise_plan(s->noise, single_base(a, c),
                                          s->chosen, 0, err, nerr);
      ok = s->final[c] != NULL;
   }
   cfg.policy = THD_PICARD_SELECT_HRF;
   cfg.nparameter = a->opt.ngrid;
   cfg.fixed_parameter = a->opt.fraction;
   cfg.parameters = a->opt.grid;
   cfg.trials = a->opt.ngrid ? a->trials : NULL;
   if (ok)
      s->collection =
          THD_picard_candidates_create((const THD_picard_plan *const *)s->final,
                                       ids, s->nh, &cfg, err, nerr);
   free(ids);
   if (!s->collection)
      return 0;
   INFO_message(
       PROGRAM
       ": staged fit: %d pool voxels; %d PCs per run; %d selection voxels",
       npool, s->chosen, s->nscore);
   return 1;
}

/* Reject values that cannot be represented in AFNI float output bricks.
   Each worker owns its error buffer and writes only its assigned voxels. */
static int store_value(THD_3dim_dataset *d, int brick, int v, double value,
                       Worker *w)
{
   if (!isfinite(value) || fabs(value) > FLT_MAX) {
      w->failed = 1;
      snprintf(w->err, sizeof(w->err), "output exceeds float range at voxel %d",
               v);
      return 0;
   }
   ((float *)DSET_ARRAY(d, brick))[v] = (float)value;
   return 1;
}

/* Record learned settings with the same bucket as the final coefficients.
   PC values are retained-row major, with each row carrying its own run's PCs.
   The augmented column identities below explain saved nuisance coefficients;
   the original task-column identities and trial manifest remain unchanged. */
static int single_output(App *a, char *err, int nerr)
{
   Single *s = a->single;
   const THD_picard_design_info *d = THD_picard_plan_info(s->final[0]);
   const double *pcs = THD_picard_noise_pcs(s->noise);
   int maxpc = a->opt.pc_counts[a->opt.npc - 1], i, j, t;
   int base = d->p + saved_nuisance_count(a) + single_start(a);
   int dims[] = {d->n, d->nfull, d->p, d->k, d->xrank, d->zrank};
   int *roles = calloc((size_t)d->ncol, sizeof(int));
   char *values = NULL, *trend = NULL, *labels = NULL, *policy = NULL;
   Worker check = {0};
   if (!roles)
      return failure(err, nerr, "cannot allocate staged provenance");
   memcpy(roles, a->roles, sizeof(int) * (size_t)a->xm->ncol);
   for (i = 0; i < a->xm->ncol; i++)
      labels = THD_zzprintf(labels, "%s%s", i ? ";" : "", a->xm->collab[i]);
   for (i = 0; i < s->added; i++)
      labels = THD_zzprintf(labels, ";NoisePC_Run%03d_PC%03d", i / s->chosen,
                            i % s->chosen);
   for (i = 0; i < d->n; i++)
      for (j = 0; j < maxpc; j++)
         values = THD_zzprintf(values, "%s%.17g", i || j ? " " : "",
                               pcs[(size_t)i * maxpc + j]);
   for (i = 0; i < a->opt.npc; i++)
      trend = THD_zzprintf(trend, "%s%d %.17g", i ? ";" : "",
                           a->opt.pc_counts[i], s->trend[i]);
   policy = THD_zzprintf(
       NULL,
       "single-trial-v1; HRF=full-data no-PC unregularized fit, frozen; "
       "pool=%s; search=%s; coarse_model=%s; "
       "task_incremental_R2_threshold=%.17g; "
       "selection_mask=%s; pc_policy=%s; pc_stop=%.17g; "
       "defined_selection_voxels=%d; "
       "PCA=per-run constant+supplied nuisance removal, unit-L2 pool voxels, "
       "temporal Gram SVD, unit sample SD, positive largest-entry sign; "
       "PC storage=retained-row-major own-run components, max_count=%d; "
       "PC loss=ordered cross-run condition pairs vs fixed no-PC OLS betas, "
       "mean squared error divided by no-PC sample beta variance; "
       "zero variance excluded; one session; median across selection mask; "
       "first count reaching best improvement/pc_stop; no improvement=0 PCs; "
       "final penalty=existing fixed/repeat policy conditional on learned "
       "settings; "
       "ALL scores are conditional selection scores, NOT independent test "
       "performance; "
       "StageB=no-PC OLS; StageC=selected-PC OLS; ordinary "
       "coefficients=StageD; "
       "external designs/HRFs; no automatic timing convolution, intensity "
       "pool, "
       "polynomials, session normalization or beta autoscaling; "
       "column labels append run-ordinal PCs after supplied design columns.",
       a->opt.noise_pool ? a->opt.noise_pool : "coarse_model_R2",
       a->opt.noise_search ? a->opt.noise_search : "NONE",
       a->opt.noise_model ? a->opt.noise_model : "NONE",
       a->opt.noise_model ? a->opt.noise_r2 : -1,
       a->opt.pc_mask ? a->opt.pc_mask : "NONE",
       a->opt.npc > 1 ? "repeated-condition selection" : "fixed",
       a->opt.pc_stop, s->nscore, maxpc);
   if (!labels || !trend || !policy || (maxpc && !values)) {
      free(roles);
      free(labels);
      free(values);
      free(trend);
      free(policy);
      return failure(err, nerr, "cannot allocate staged provenance text");
   }
   for (t = 0; t < DATASET_COUNT; t++)
      if (a->out[t]) {
         THD_datablock *b = a->out[t]->dblk;
         THD_set_string_atr(b, "PICARD_SINGLE_TRIAL", policy);
         THD_set_string_atr(b, "PICARD_PC_TREND", trend);
         if (maxpc)
            THD_set_string_atr(b, "PICARD_NOISE_PCS", values);
         THD_set_int_atr(b, "PICARD_PC_COUNT", 1, &s->chosen);
         THD_set_int_atr(b, "PICARD_PC_RANKS", s->nrun,
                         (int *)THD_picard_noise_ranks(s->noise));
         THD_set_int_atr(b, "PICARD_DIMENSIONS", 6, dims);
         THD_set_int_atr(b, "PICARD_COLUMN_ROLES", d->ncol, roles);
         THD_set_string_atr(b, "PICARD_COLUMN_LABELS", labels);
         if (d->k)
            THD_set_int_atr(b, "PICARD_NUISANCE_COLUMNS", d->k,
                            (int *)d->nuisance_columns);
         tross_Append_History(a->out[t], policy);
      }
   for (i = 0; i < a->nvox; i++) {
      store_value(a->out[OUT_MAIN], base + 1, i, s->pool[i], &check);
      store_value(a->out[OUT_MAIN], base + 2, i, s->select[i], &check);
      store_value(a->out[OUT_MAIN], base + 3, i, s->noise_r2[i], &check);
   }
   for (i = 0; i < a->nselected; i++) {
      int v = a->vox[i];
      store_value(a->out[OUT_MAIN], base, v, s->chosen, &check);
      if (a->opt.save_stages)
         for (j = 0; j < d->p; j++)
            store_value(a->out[OUT_MAIN], base + 4 + j, v,
                        s->beta_b[(size_t)j * a->nvox + v], &check);
   }
   free(roles);
   free(labels);
   free(values);
   free(trend);
   free(policy);
   return check.failed ? failure(err, nerr, check.err) : 1;
}

/* Gather retained time points into a row-major [time][response] batch,
   run fixed fitting or repeated-condition selection, then scatter the
   borrowed numerical results into spatial AFNI bricks before reusing scratch.
 */
static void fit_block(App *a, Worker *w, int start)
{
   int q = a->nselected - start, i, j, o;
   if (q > a->opt.batch)
      q = a->opt.batch;
   const THD_picard_result *r;
   const THD_picard_repeat_result *tuned = NULL;
   const THD_picard_predict_result *predicted = NULL;
   const THD_picard_candidates_result *chosen = NULL;
   if (w->failed)
      return;
   /* There are two independent index maps here. good[i] maps a retained
      matrix row to its original acquisition time. vox[start+o] maps a batch
      column to a spatial voxel. Censored samples are never passed to fits.
      Use q (this batch's width), not batch capacity, as the row stride. */
   for (i = 0; i < a->xm->nrow; i++)
      for (o = 0; o < q; o++)
         w->y[(size_t)i * q + o] =
             read_voxel(a->input, a->xm->good[i], a->vox[start + o]);
   /* Each workflow returns the same final-fit view. Tuning workflows do
      their selection first and then refit the chosen setting on all retained
      rows. The extra result pointers below only expose mode-specific scores. */
   if (w->candidates) {
      int ok;
      if (a->single) {
         for (o = 0; o < q; o++) {
            int v = a->vox[start + o];
            w->locked_hrf[o] = a->single->hrf[v];
            w->locked_loss[o] = a->single->hrf_loss[v];
         }
         ok = THD_picard_candidates_fit_selected(w->candidates, w->y, q,
                                                 w->locked_hrf, w->locked_loss,
                                                 w->err, sizeof(w->err));
      } else
         ok = THD_picard_candidates_fit(w->candidates, w->y, q, w->err,
                                        sizeof(w->err));
      if (!ok) {
         w->failed = 1;
         return;
      }
      chosen = THD_picard_candidates_get_result(w->candidates);
      r = chosen->fit;
   } else if (w->predict) {
      if (!THD_picard_predict_fit(w->predict, w->y, q, 0, a->global_index,
                                  w->err, sizeof(w->err))) {
         w->failed = 1;
         return;
      }
      predicted = THD_picard_predict_get_result(w->predict);
      r = predicted->fit;
   } else if (w->repeat) {
      if (!THD_picard_repeat_fit(w->repeat, w->y, q, w->err, sizeof(w->err))) {
         w->failed = 1;
         return;
      }
      tuned = THD_picard_repeat_get_result(w->repeat);
      r = tuned->fit;
   } else {
      if (!THD_picard_fit_batch(w->work, w->y, q, w->fractions, w->err,
                                sizeof(w->err))) {
         w->failed = 1;
         return;
      }
      r = THD_picard_get_result(w->work);
   }
   /* Results borrow workspace storage. Copy them into output maps before
      the next batch overwrites that storage. All workers share the maps but
      write disjoint spatial indices, so no per-voxel locking is needed. */
   for (o = 0; o < q; o++) {
      int v = a->vox[start + o];
      /* Some mathematically undefined diagnostics are NaN/Inf in the
         numerical API. Translate those documented cases to -1 for AFNI
         maps; store_value still rejects accidental nonfinite fit outputs. */
      double diagnostics[NELASTIC_DIAG + 2] = {0};
      if (a->opt.covariance_set) {
         int end = DSET_NVALS(a->out[OUT_MAIN]);
         if (!store_value(a->out[OUT_MAIN], end-2, v,
                          a->opt.covariance.a, w) ||
             !store_value(a->out[OUT_MAIN], end-1, v,
                          a->opt.covariance.b, w))
            return;
      }
      if (r->banded_diagnostic) {
         const THD_banded_diagnostic *dg = r->banded_diagnostic + o;
         diagnostics[0] = dg->candidate;
         diagnostics[1] = dg->objective;
         diagnostics[2] = dg->gradient;
         diagnostics[3] =
             chosen ? chosen->rank[o] : THD_picard_plan_info(a->plan)->xrank;
      } else if (r->elastic_diagnostic) {
         const THD_elastic_diagnostic *dg = r->elastic_diagnostic + o;
         diagnostics[0] = dg->lambda;
         diagnostics[1] = dg->rho;
         diagnostics[2] = dg->solver.nonzero;
         diagnostics[3] = dg->solver.iterations;
         diagnostics[4] = dg->solver.residual;
         diagnostics[5] = dg->solver.objective;
         diagnostics[6] = isnan(dg->lambda_max) ? -1 : dg->lambda_max;
         diagnostics[7] =
             chosen ? chosen->rank[o] : THD_picard_plan_info(a->plan)->xrank;
         diagnostics[8] = a->opt.solver.standardize;
         diagnostics[9] = 0; /* Failed convergence aborts before output. */
         if (r->support_size) {
            diagnostics[9] = r->support_size[o];
            diagnostics[10] = r->support_rank[o];
         }
      } else {
         const THD_fraction_diagnostic *dg = r->diagnostic + o;
         diagnostics[0] = dg->fraction;
         diagnostics[1] = isinf(dg->lambda) ? -1 : dg->lambda;
         diagnostics[2] = isinf(dg->penalty) ? -1 : dg->penalty;
         diagnostics[3] = dg->reference_norm;
         diagnostics[4] =
             isnan(dg->achieved_fraction) ? -1 : dg->achieved_fraction;
         diagnostics[5] = dg->rank;
         diagnostics[6] = dg->status;
      }
      for (j = 0; j < r->p; j++)
         if (!store_value(a->out[OUT_MAIN], j, v, r->beta[(size_t)j * q + o],
                          w))
            return;
      for (j = 0; j < base_diagnostic_count(a); j++) {
         double value;
         if (a->opt.banded && j >= 4)
            value = j == 4 + a->opt.ngroup
                        ? 0
                        : a->opt.band_values[(size_t)r->banded_diagnostic[o]
                                                     .candidate *
                                                 a->opt.ngroup +
                                             j - 4];
         else
            value = diagnostics[j];
         if (!store_value(a->out[OUT_MAIN],
                          r->p + (a->opt.save_nuis ? r->k : 0) + j, v, value,
                          w))
            return;
      }
      if (tuned || (chosen && a->repeat)) {
         int base =
             r->p + (a->opt.save_nuis ? r->k : 0) + base_diagnostic_count(a);
         double scores[] = {
             chosen ? chosen->loss[o] : tuned->loss[o],
             chosen ? chosen->baseline[o] : tuned->reference_loss[o],
             chosen ? chosen->parameter_index[o] : tuned->candidate[o]};
         for (j = 0; j < NTUNE_DIAG; ++j)
            if (!store_value(a->out[OUT_MAIN], base + j, v, scores[j], w))
               return;
         if (a->opt.save_cv)
            for (j = 0; j < a->opt.ngrid; ++j)
               if (!store_value(
                       a->out[OUT_MAIN], base + NTUNE_DIAG + j, v,
                       (chosen ? chosen->parameter_loss
                               : tuned->candidate_loss)[(size_t)j * q + o],
                       w))
                  return;
      }
      if (predicted || (chosen && a->predict)) {
         int base =
             r->p + (a->opt.save_nuis ? r->k : 0) + base_diagnostic_count(a);
         int status =
             chosen ? chosen->prediction_status[o] : predicted->status[o];
         double scores[] = {
             chosen ? chosen->loss[o] : predicted->loss[o],
             chosen ? chosen->baseline[o] : predicted->baseline[o],
             status ? 0 : (chosen ? chosen->r2[o] : predicted->r2[o]),
             chosen ? chosen->parameter_index[o] : predicted->candidate[o],
             status};
         for (j = 0; j < NPRED_DIAG; ++j)
            if (!store_value(a->out[OUT_MAIN], base + j, v, scores[j], w))
               return;
         if (a->opt.save_cv)
            for (j = 0; j < a->opt.ngrid; ++j)
               if (!store_value(
                       a->out[OUT_MAIN], base + NPRED_DIAG + j, v,
                       (chosen ? chosen->parameter_loss
                               : predicted->candidate_loss)[(size_t)j * q + o],
                       w))
                  return;
      }
      if (chosen && a->manifest) {
         int base =
             r->p + (a->opt.save_nuis ? r->k : 0) + design_diagnostic_start(a);
         double scores[] = {chosen->design_id[o], chosen->design_index[o],
                            chosen->selection_loss[o]};
         for (j = 0; j < NDESIGN_DIAG; j++)
            if (!store_value(a->out[OUT_MAIN], base + j, v, scores[j], w))
               return;
         if (a->opt.save_candidates)
            for (j = 0; j < a->manifest->count; j++)
               if (!store_value(
                       a->out[OUT_MAIN], base + NDESIGN_DIAG + j, v,
                       a->single
                           ? a->single->design_loss[(size_t)j * a->nvox + v]
                           : chosen->design_loss[(size_t)j * q + o],
                       w))
                  return;
      }
      if (a->opt.save_nuis)
         for (j = 0; j < r->k; j++)
            if (!store_value(a->out[OUT_MAIN], r->p + j, v,
                             r->gamma[(size_t)j * q + o], w))
               return;
      if (a->out[OUT_FITTED])
         for (i = 0; i < r->n; i++)
            if (!store_value(a->out[OUT_FITTED], i, v,
                             r->fitted[(size_t)i * q + o], w))
               return;
      if (a->out[OUT_RESIDUAL])
         for (i = 0; i < r->n; i++)
            if (!store_value(a->out[OUT_RESIDUAL], i, v,
                             r->residual[(size_t)i * q + o], w))
               return;
   }
   if (a->single && a->opt.save_stages) {
      int h, base = r->p + saved_nuisance_count(a) + single_start(a) + 4 + r->p;
      for (h = 0; h < a->single->nh; h++) {
         const THD_picard_result *ols;
         if (!THD_picard_fit_ols(w->stage_c[h], w->y, q, w->err,
                                 sizeof(w->err))) {
            w->failed = 1;
            return;
         }
         ols = THD_picard_get_result(w->stage_c[h]);
         for (o = 0; o < q; o++)
            if (a->single->hrf[a->vox[start + o]] == h)
               for (j = 0; j < ols->p; j++)
                  if (!store_value(a->out[OUT_MAIN], base + j,
                                   a->vox[start + o],
                                   ols->beta[(size_t)j * q + o], w))
                     return;
      }
   }
}

/* Global selection is a serial, bounded-memory prepass. Accumulate voxels in
   spatial-index order regardless of batch size or the final OpenMP schedule.
   Normalize by each voxel's pooled training-mean baseline MSE; exclude exactly
   zero baselines. The second pass refits every masked voxel at one fraction. */
static int choose_global_fraction(App *a, char *err, int nerr)
{
   THD_picard_predict_workspace *work;
   THD_picard_predict_global *global;
   const THD_picard_predict_global_result *result = NULL;
   double *y;
   int start, i, o, ok = 1;
   if (!a->opt.global)
      return 1;
   y = calloc((size_t)a->xm->nrow * a->opt.batch, sizeof(double));
   work =
       THD_picard_predict_workspace_alloc(a->predict, a->opt.batch, err, nerr);
   global = THD_picard_predict_global_alloc(a->predict, err, nerr);
   if (!y || !work || !global)
      ok = failure(err, nerr, "cannot allocate global prediction prepass");
   for (start = 0; start < a->nselected && ok; start += a->opt.batch) {
      int q = a->nselected - start;
      if (q > a->opt.batch)
         q = a->opt.batch;
      for (i = 0; i < a->xm->nrow; ++i)
         for (o = 0; o < q; ++o)
            y[(size_t)i * q + o] =
                read_voxel(a->input, a->xm->good[i], a->vox[start + o]);
      ok = THD_picard_predict_fit(work, y, q, 1, -1, err, nerr) &&
           THD_picard_predict_global_add(global, work, err, nerr);
   }
   if (ok) {
      result = THD_picard_predict_global_finish(global, err, nerr);
      ok = result != NULL;
   }
   if (ok) {
      char *description;
      size_t cell, cells = (size_t)a->folds->nfold * a->opt.ngrid;
      a->global_index = result->candidate;
      description = THD_zzprintf(
          NULL,
          "index=%d; parameter=%.17g; defined_voxels=%d; best_fold_SE=%.17g; "
          "normalization=each_voxel_pooled_baseline_MSE; "
          "zero_baseline=excluded; all_undefined=weakest_penalty; "
          "fold_major_normalized_MSE=",
          result->candidate, a->opt.grid[result->candidate], result->nvalid,
          result->best_se);
      for (cell = 0; cell < cells; ++cell)
         description = THD_zzprintf(description, "%s%.17g", cell ? "," : "",
                                    result->fold_loss[cell]);
      if (!description)
         ok = failure(err, nerr, "cannot record global CV provenance");
      else {
         for (i = 0; i < DATASET_COUNT; ++i)
            if (a->out[i]) {
               THD_set_string_atr(a->out[i]->dblk, "PICARD_PRED_GLOBAL",
                                  description);
               tross_Append_History(a->out[i], description);
            }
         free(description);
      }
   }
   THD_picard_predict_global_free(global);
   THD_picard_predict_workspace_free(work);
   free(y);
   return ok;
}

/* Allocate one mutable workspace per worker around shared immutable plans.
   Static scheduling gives disjoint voxel batches; workers never write files.
   Report any failure before the transaction reaches the output-writing stage.
 */
static int fit_voxels(App *a, char *err, int nerr)
{
   int t, b, blocks;
   if (!choose_global_fraction(a, err, nerr))
      return 0;
   a->workers = calloc((size_t)a->opt.jobs, sizeof(*a->workers));
   if (!a->workers)
      return failure(err, nerr, "cannot allocate workers");
   for (t = 0; t < a->opt.jobs; t++) {
      Worker *w = a->workers + t;
      int q = a->opt.batch;
      if (a->candidate_plan || a->single)
         w->candidates = THD_picard_candidates_workspace_alloc(
             a->single ? a->single->collection : a->candidate_plan, q, err,
             nerr);
      else if (a->predict)
         w->predict =
             THD_picard_predict_workspace_alloc(a->predict, q, err, nerr);
      else if (a->repeat)
         w->repeat = THD_picard_repeat_workspace_alloc(a->repeat, q, err, nerr);
      else
         w->work = THD_picard_workspace_alloc(a->plan, q, err, nerr);
      if (a->single) {
         int c;
         w->nstage = a->single->nh;
         w->locked_hrf = calloc((size_t)q, sizeof(int));
         w->locked_loss = calloc((size_t)q, sizeof(double));
         if (a->opt.save_stages) {
            w->stage_c = calloc((size_t)w->nstage, sizeof(*w->stage_c));
            if (!w->stage_c)
               return failure(err, nerr, "cannot allocate stage workers");
            for (c = 0; c < w->nstage; c++) {
               w->stage_c[c] = THD_picard_workspace_alloc(a->single->final[c],
                                                          q, err, nerr);
               if (!w->stage_c[c])
                  return 0;
            }
         }
         if (!w->locked_hrf || !w->locked_loss)
            return failure(err, nerr, "cannot allocate locked HRFs");
      }
      w->y = calloc((size_t)a->xm->nrow * q, sizeof(double));
      w->fractions = calloc((size_t)q, sizeof(double));
      if ((!w->work && !w->repeat && !w->predict && !w->candidates) || !w->y ||
          !w->fractions)
         return failure(err, nerr, "cannot allocate response batch");
      for (b = 0; b < q; b++)
         w->fractions[b] = a->opt.fraction;
   }
   /* Integer ceiling division includes a possibly shorter final batch.
      A worker can process several batches; its scratch is reused each time.
      A failed worker skips later batches and main suppresses all file writes.
    */
   blocks = a->nselected / a->opt.batch + (a->nselected % a->opt.batch != 0);
#ifdef USE_OMP
#pragma omp parallel for num_threads(a->opt.jobs) schedule(static)
#endif
   for (b = 0; b < blocks; b++) {
      int tid = 0;
#ifdef USE_OMP
      tid = omp_get_thread_num();
#endif
      fit_block(a, a->workers + tid, b * a->opt.batch);
   }
   for (t = 0; t < a->opt.jobs; t++)
      if (a->workers[t].failed)
         return failure(err, nerr, a->workers[t].err);
   return 1;
}

/* Escape quotes in a TSV field using doubled quotes, as CSV readers expect. */
static void write_quoted_field(FILE *f, const char *s)
{
   fputc('"', f);
   for (; *s; s++) {
      if (*s == '"')
         fputc('"', f);
      fputc(*s, f);
   }
   fputc('"', f);
}

/* Write an explicitly requested sidecar with exclusive creation. Trial
   identities are joined by penalized-column order, never guessed from labels.
 */
static int write_manifest(App *a, int which, char *err, int nerr)
{
   FILE *f = fopen(a->names[which], "wx");
   const THD_picard_design_info *d = THD_picard_plan_info(a->plan);
   int i, j, ok;
   if (!f)
      return failure(err, nerr, "cannot create metadata output");
   if (which == OUT_BETAS_TSV && a->manifest) {
      int design_volume =
          d->p + saved_nuisance_count(a) + design_diagnostic_start(a);
      fputs("InputFile\tVolume\tColumn\tLabel\tTrial\tCondition\tRun\tBasis"
            "\tItem\tCandidateManifest\tDesignIDFile\tDesignIDVolume\n",
            f);
      for (i = 0; i < d->p; i++) {
         int c = d->penalized_columns[i];
         const PicardColumnIdentity *ci = a->manifest->entry[0].identity + c;
         write_quoted_field(f, DSET_HEADNAME(a->out[OUT_MAIN]));
         fprintf(f, "\t%d\t%d\t", i, c);
         write_quoted_field(f, ci->label);
         fprintf(f, "\t%d\t%d\t%d\t%d\t", ci->trial, ci->condition, ci->run,
                 ci->basis);
         write_quoted_field(f, ci->item);
         fputc('\t', f);
         write_quoted_field(f, a->manifest->path);
         fputc('\t', f);
         write_quoted_field(f, DSET_HEADNAME(a->out[OUT_MAIN]));
         fprintf(f, "\t%d\n", design_volume);
      }
   } else if (which == OUT_BETAS_TSV) {
      fputs("InputFile\tVolume\tXmat\tXColumn\tLabel\tStimulus", f);
      if (a->repeat)
         fputs("\tTrial\tCondition\tRun", f);
      fputc('\n', f);
      for (i = 0; i < d->p; i++) {
         int c = d->penalized_columns[i];
         const char *stim = "-";
         for (j = 0; j < a->xm->nstim; j++)
            if (c >= a->xm->stim[j].bot && c <= a->xm->stim[j].top)
               stim = a->xm->stim[j].label;
         write_quoted_field(f, DSET_HEADNAME(a->out[OUT_MAIN]));
         fprintf(f, "\t%d\t", i);
         write_quoted_field(f, a->opt.matrix);
         fprintf(f, "\t%d\t", c);
         write_quoted_field(f, a->xm->collab[c]);
         fputc('\t', f);
         write_quoted_field(f, stim);
         if (a->repeat)
            fprintf(f, "\t%d\t%d\t%d", a->trials[i].trial,
                    a->trials[i].condition, a->trials[i].run);
         fputc('\n', f);
      }
   } else if (which == OUT_ROWS_1D) {
      fputs("# Columns: RetainedRow FullRow Run\n", f);
      for (i = 0; i < d->n; i++)
         fprintf(f, "%d %d %d\n", i, d->row_ids[i], d->run_ids[i]);
   }
   ok = !ferror(f);
   if (fclose(f))
      ok = 0;
   if (!ok) {
      remove(a->names[which]);
      return failure(err, nerr, "failed writing metadata");
   }
   return 1;
}

/* Commit outputs only after fitting succeeds. Track files created here and
   roll them back on a later write failure, preserving pre-existing paths. */
static int write_outputs(App *a, char *err, int nerr)
{
   int i, created[OUTPUT_COUNT] = {0}, ok = 1;
   for (i = 0; i < DATASET_COUNT; i++)
      if (a->out[i] && dataset_exists(a->out[i]))
         return failure(err, nerr,
                        "output appeared during fitting; refusing overwrite");
   for (i = DATASET_COUNT; i < OUTPUT_COUNT; i++)
      if (a->names[i] && path_exists(a->names[i]))
         return failure(err, nerr,
                        "manifest appeared during fitting; refusing overwrite");
   /* Sidecars use exclusive creation. Dataset writers can create more than
      one file (HEAD plus BRIK), so mark that dataset before trying the write:
      a failed write may already have created a partial file to clean up. */
   for (i = DATASET_COUNT; i < OUTPUT_COUNT && ok; i++)
      if (a->names[i]) {
         ok = write_manifest(a, i, err, nerr);
         created[i] = ok;
      }
   for (i = 0; i < DATASET_COUNT && ok; i++)
      if (a->out[i]) {
         created[i] = 1;
         if (DSET_write(a->out[i]) != True)
            ok = failure(err, nerr, "failed writing Picard dataset");
      }
   /* This is best-effort rollback of our own outputs, not an atomic commit
      across files. Never remove a path that this attempt did not create. */
   if (!ok)
      for (i = 0; i < OUTPUT_COUNT; i++)
         if (created[i]) {
            if (i < DATASET_COUNT)
               remove_dataset(a->out[i]);
            else
               remove(a->names[i]);
         }
   return ok;
}

/* Release partially or fully initialized state. Workers borrow their plans,
   so destroy workers first, then repeat plans, then the full numerical plan. */
static void free_app(App *a)
{
   int i;
   if (a->workers)
      for (i = 0; i < a->opt.jobs; i++) {
         THD_picard_candidates_workspace_free(a->workers[i].candidates);
         THD_picard_workspace_free(a->workers[i].work);
         THD_picard_repeat_workspace_free(a->workers[i].repeat);
         THD_picard_predict_workspace_free(a->workers[i].predict);
         free(a->workers[i].y);
         free(a->workers[i].fractions);
         free(a->workers[i].locked_hrf);
         free(a->workers[i].locked_loss);
         if (a->workers[i].stage_c) {
            int c;
            for (c = 0; c < a->workers[i].nstage; c++)
               THD_picard_workspace_free(a->workers[i].stage_c[c]);
            free(a->workers[i].stage_c);
         }
      }
   free(a->workers);
   for (i = 0; i < DATASET_COUNT; i++)
      if (a->out[i])
         DSET_delete(a->out[i]);
   if (a->input)
      DSET_delete(a->input);
   if (a->mask)
      DSET_delete(a->mask);
   for (i = 0; i < OUTPUT_COUNT; i++)
      free(a->names[i]);
   if (a->single) {
      Single *s = a->single;
      int c;
      THD_picard_candidates_free(s->collection);
      if (s->final)
         for (c = 0; c < s->nh; c++)
            THD_picard_plan_free(s->final[c]);
      THD_picard_noise_free(s->noise);
      free(s->final);
      free(s->pool);
      free(s->select);
      free(s->hrf);
      free(s->hrf_loss);
      free(s->design_loss);
      free(s->noise_r2);
      free(s->beta_b);
      free(s->trend);
      free(s);
   }
   THD_picard_candidates_free(a->candidate_plan);
   if (a->candidate_designs && a->manifest)
      for (i = 1; i < a->manifest->count; i++)
         THD_picard_plan_free(a->candidate_designs[i]);
   free(a->candidate_designs);
   picard_manifest_free(a->manifest);
   THD_picard_repeat_plan_free(a->repeat);
   THD_picard_predict_plan_free(a->predict);
   THD_foldplan_free(a->folds);
   free(a->trials);
   free(a->opt.pc_counts);
   free(a->opt.grid);
   free(a->opt.band_values);
   free(a->groups);
   THD_picard_plan_free(a->plan);
   THD_xmat_free(a->xm);
   free(a->roles);
   free(a->runs);
   free(a->vox);
   free(a->matrix);
}

/* Exact parameter grouping bounds resident factors to one covariance group
   (plus the metadata reference plans). Stable voxel IDs break ties, making
   batching/thread choices immaterial. Imported parameters are not rounded. */
typedef struct {
   double ar, ma;
   int voxel;
   THD_picard_reml_result estimate;
} CovarianceVoxel;
static int covariance_voxel_compare(const void *left, const void *right)
{
   const CovarianceVoxel *a = left, *b = right;
   if (a->ar != b->ar) return a->ar < b->ar ? -1 : 1;
   if (a->ma != b->ma) return a->ma < b->ma ? -1 : 1;
   return (a->voxel > b->voxel) - (a->voxel < b->voxel);
}

/* Estimate on explicitly selected original rows, never on an implicitly
   whitened/full-response series. Pilot columns are a fixed mean-model span;
   compact and full-trial models are both supplied as ordinary xmat files.
   Responses are batched; no candidate-by-voxel factor cache is retained. */
static int estimate_covariance(App *a, CovarianceVoxel **values,
                               int validate, char *err, int nerr)
{
   THD_xmat *pilot = THD_xmat_read(a->opt.covariance_estimate, err, nerr);
   int *rows = NULL, *runs = NULL;
   double *table = NULL, *matrix = NULL, *y = NULL;
   THD_picard_reml_result *result = NULL;
   int n = 0, nr = 0, nc = 0, i, j, at, q, ok = pilot != NULL;
   if (ok && (pilot->nrow != a->xm->nrow ||
       pilot->nfull != a->xm->nfull || pilot->nrun != a->xm->nrun ||
       memcmp(pilot->good, a->xm->good, sizeof(int)*(size_t)pilot->nrow)))
      ok = failure(err, nerr, "covariance pilot timeline mismatch");
   if (ok) {
      const char *pt = THD_xmat_attr(pilot, "RowTR");
      const char *ft = THD_xmat_attr(a->xm, "RowTR");
      double pv, fv;
      if ((pt == NULL) != (ft == NULL) ||
          (pt && (!parse_finite_double(pt, &pv) ||
                  !parse_finite_double(ft, &fv) || pv <= 0 || pv != fv)))
         ok = failure(err, nerr, "covariance pilot TR mismatch");
   }
   for (i = 0; ok && i < pilot->nrow; i++)
      if (THD_xmat_run_of(pilot, pilot->good[i]) != a->runs[i])
         ok = failure(err, nerr, "covariance pilot runs differ");
   if (ok && a->opt.covariance_rows)
      ok = read_numeric_table(a->opt.covariance_rows, &table, &nr, &nc,
                              a->opt.memory, err, nerr);
   if (ok && table && nc != 1)
      ok = failure(err, nerr,
                   "covariance_rows needs one original row per line");
   if (ok) {
      n = table ? nr : pilot->nrow;
      if (n < 2 || 32.0L*n*((long double)n+pilot->ncol+a->opt.batch) >
          a->opt.memory*1048576.0L/2)
         ok = failure(err, nerr, "native covariance rows/memory invalid");
   }
   if (ok) {
      rows = calloc((size_t)n, sizeof(int));
      runs = calloc((size_t)n, sizeof(int));
      matrix = calloc((size_t)n*pilot->ncol, sizeof(double));
      y = calloc((size_t)n*a->opt.batch, sizeof(double));
      result = calloc((size_t)a->opt.batch, sizeof(*result));
      *values = calloc((size_t)a->nselected, sizeof(**values));
      ok = rows && runs && matrix && y && result && *values;
      if (!ok) failure(err, nerr, "cannot allocate native covariance buffers");
   }
   for (i = at = 0; ok && i < n; i++) {
      double row = table ? table[i] : pilot->good[i];
      if (!isfinite(row) || row < 0 || row >= pilot->nfull ||
          row != floor(row) || (i && row <= rows[i-1])) {
         ok = failure(err, nerr, "invalid covariance estimation rows");
         break;
      }
      while (at < pilot->nrow && pilot->good[at] < row) at++;
      if (at == pilot->nrow || pilot->good[at] != row) {
         ok = failure(err, nerr, "covariance row is not retained");
         break;
      }
      rows[i] = (int)row;
      runs[i] = a->runs[at];
      for (j = 0; j < pilot->ncol; j++)
         matrix[(size_t)i*pilot->ncol+j] =
             pilot->mat[(size_t)j*pilot->nrow+at];
   }
   for (at = 0; ok && at < a->nselected; at += q) {
      q = a->nselected-at;
      if (q > a->opt.batch) q = a->opt.batch;
      for (i = 0; i < n; i++)
         for (j = 0; j < q; j++)
            y[(size_t)i*q+j] = read_voxel(a->input, rows[i], a->vox[at+j]);
      ok = THD_picard_reml(matrix, n, pilot->ncol, rows, runs, y, q,
                           a->opt.covariance.model,
                           a->opt.covariance.reset_gaps, a->opt.tol,
                           result, err, nerr);
      for (j = 0; ok && j < q; j++) {
         CovarianceVoxel *v = *values+at+j;
         v->ar = result[j].a; v->ma = result[j].b;
         v->voxel = a->vox[at+j]; v->estimate = result[j];
         if (!validate) {
            Worker check = {0};
            int base = DSET_NVALS(a->out[OUT_MAIN])-6;
            ok = store_value(a->out[OUT_MAIN], base, v->voxel,
                              result[j].sigma2, &check) &&
                 store_value(a->out[OUT_MAIN], base+1, v->voxel,
                              result[j].score, &check) &&
                 store_value(a->out[OUT_MAIN], base+2, v->voxel,
                              result[j].rank, &check) &&
                 store_value(a->out[OUT_MAIN], base+3, v->voxel,
                              result[j].status, &check);
            if (!ok) failure(err, nerr, check.err);
         }
      }
   }
   if (ok && !validate) {
      char *description = THD_zzprintf(NULL,
          "native_grid_REML pilot=%s grid=-0.9:0.1:0.9 "
          "pooling=per_voxel shared_across_runs=1 shared_variance=1 df=%d "
          "zero_variance=error status=0:interior,1:boundary,3:flat,4:tied "
          "score=logdetR+logdetBtRiB+df*logRSS_over_df rows=",
          a->opt.covariance_estimate, n-(*values)[0].estimate.rank);
      for (i = 0; i < n; i++)
         description = THD_zzprintf(description, "%s%d", i ? "," : "", rows[i]);
      if (!description)
         ok = failure(err, nerr, "cannot record REML provenance");
      for (j = 0; ok && j < DATASET_COUNT; j++)
         if (a->out[j]) THD_set_string_atr(a->out[j]->dblk,
                                          "PICARD_COVARIANCE_ESTIMATION",
                                          description);
      free(description);
   }
   free(rows);
   free(runs);
   free(table);
   free(matrix);
   free(y);
   free(result);
   THD_xmat_free(pilot);
   return ok;
}

/* Validate all selected map values before fitting any group. Dataset loading
   follows AFNI_FLOATSCAN; original values repaired by AFNI cannot be detected
   afterward. Extra Rvar bricks (lag-one correlation, sigma, likelihood) are
   not used as covariance parameters. No files are written on failure. */
static int fit_covariance_groups(App *a, int validate, char *err, int nerr)
{
   THD_3dim_dataset *map = NULL;
   CovarianceVoxel *values = NULL;
   int i, j, start, stop, groups = 0, ok = 1;
   int nb = a->opt.covariance.model == 1 ? 1 : 2;
   if (a->opt.covariance_estimate)
      ok = estimate_covariance(a, &values, validate, err, nerr);
   else {
      if (!a->input || strchr(a->opt.covariance_map, '['))
         return failure(err, nerr,
                        "covariance_map needs input and an unselected dataset");
      map = THD_open_dataset(a->opt.covariance_map);
      if (!map || DSET_IS_1D(map) || DSET_IS_GIFTI(map) ||
          DSET_IS_NI_SURF_DSET(map) || DSET_NVALS(map) < nb ||
          !EQUIV_GRIDS(a->input, map) || THD_dataset_mismatch(a->input, map))
         ok = failure(err, nerr, "covariance map grid/brick count mismatch");
      if (ok) {
         /* Validate every loaded brick type and account for the entire dataset
            payload before DSET_load, including ignored Rvar diagnostic bricks. */
         long double size = 0;
         for (i = 0; i < DSET_NVALS(map); i++) {
            if (!is_real_datum(DSET_BRICK_TYPE(map, i)))
               ok = failure(err, nerr, "covariance map must be real numeric");
            size += (long double)DSET_NVOX(map)*
                    mri_datum_size(DSET_BRICK_TYPE(map, i));
         }
         if (size > a->opt.memory*1048576.0L/4)
            ok = failure(err, nerr, "covariance map exceeds reserved memory");
      }
      if (ok) {
         DSET_load(map);
         if (!DSET_LOADED(map))
            ok = failure(err, nerr, "cannot load covariance map");
      }
      if (ok) {
         values = calloc((size_t)a->nselected, sizeof(*values));
         if (!values) ok = failure(err, nerr, "cannot allocate covariance groups");
      }
      for (i = 0; ok && i < a->nselected; i++) {
         values[i].voxel = a->vox[i];
         values[i].ar = read_voxel(map, 0, a->vox[i]);
         values[i].ma = nb == 2 ? read_voxel(map, 1, a->vox[i]) : 0;
         if (!isfinite(values[i].ar) || !isfinite(values[i].ma) ||
             fabs(values[i].ar) >= 1 || fabs(values[i].ma) >= 1)
            ok = failure(err, nerr, "invalid selected covariance map parameter");
      }
      if (ok && !validate) {
         ATR_string *history = THD_find_string_atr(map->dblk, "HISTORY_NOTE");
         for (j = 0; j < DATASET_COUNT; j++)
            if (a->out[j]) {
               THD_set_string_atr(a->out[j]->dblk, "PICARD_COVARIANCE_SOURCE_ID",
                                  map->idcode.str);
               if (history && history->ch)
                  THD_set_string_atr(a->out[j]->dblk,
                                     "PICARD_COVARIANCE_SOURCE_HISTORY",
                                     history->ch);
            }
      }
      if (map) DSET_delete(map);
   }
   if (ok) qsort(values, (size_t)a->nselected, sizeof(*values),
                 covariance_voxel_compare);
   for (start = 0; ok && start < a->nselected; start = stop) {
      App g = {0};
      stop = start+1;
      while (stop < a->nselected && values[stop].ar == values[start].ar &&
             values[stop].ma == values[start].ma) stop++;
      groups++;
      g.opt = a->opt;
      g.opt.covariance_map = NULL;
      g.opt.covariance_estimate = g.opt.covariance_rows = NULL;
      g.opt.covariance_source = NULL;
      g.opt.covariance.a = values[start].ar;
      g.opt.covariance.b = values[start].ma;
      ok = prepare_candidate_manifest(&g, err, nerr) &&
           prepare_design(&g, NULL, err, nerr) &&
           prepare_candidate_designs(&g, err, nerr);
      if (ok && !validate) {
         if (g.input) DSET_delete(g.input);
         g.input = a->input; /* Borrow loaded data and output buckets. */
         for (j = 0; j < DATASET_COUNT; j++) g.out[j] = a->out[j];
         g.nselected = stop-start;
         g.vox = calloc((size_t)g.nselected, sizeof(int));
         if (!g.vox) ok = failure(err, nerr, "cannot allocate group voxels");
         for (i = 0; ok && i < g.nselected; i++)
            g.vox[i] = values[start+i].voxel;
         if (ok) ok = fit_voxels(&g, err, nerr);
         g.input = NULL;
         for (j = 0; j < DATASET_COUNT; j++) g.out[j] = NULL;
      }
      /* Options borrow the parent's heap arrays. All other group resources
         are owned by g, including reopened header-only input/mask handles. */
      g.opt.grid = NULL;
      g.opt.band_values = NULL;
      g.opt.pc_counts = NULL;
      free_app(&g);
   }
   if (ok) {
      char summary[128];
      snprintf(summary, sizeof(summary),
               "groups=%d exact_parameters=1 resident_groups=1", groups);
      INFO_message("Covariance: %s", summary);
      if (!validate)
         for (j = 0; j < DATASET_COUNT; j++)
            if (a->out[j]) THD_set_string_atr(a->out[j]->dblk,
                                            "PICARD_COVARIANCE_GROUPS",
                                            summary);
   }
   free(values);
   return ok;
}

/* Adapter pipeline: parse -> manifest/designs -> allocate -> load -> fit ->
   write. Validation stops after plan creation. All paths converge on ordinary
   cleanup and return; no jump-based error handling is needed. */
int main(int argc, char **argv)
{
   App a;
   char err[1024] = "Picard failed";
   int ok;
   if (argc == 1 || !strcmp(argv[1], "-help") || !strcmp(argv[1], "--help")) {
      help();
      return 0;
   }
   memset(&a, 0, sizeof(a));
   mainENTRY(PROGRAM);
   machdep();
   AFNI_logger(PROGRAM, argc, argv);
   /* Leave AFNI_FLOATSCAN and other reader settings under AFNI/user control.
      Both input and mask use the same loader as other AFNI programs. */
   /* Small spatial grids must remain datasets, not become text files. */
   AFNI_setenv("AFNI_WRITE_1D_AS_PREFIX=YES");
   /* && evaluates left to right and stops at the first failure. Later
      stages therefore never consume an incomplete plan. Every failure still
      reaches free_app below; zero initialization makes that cleanup safe. */
   ok = parse_options(argc, argv, &a.opt, err, sizeof(err)) &&
        prepare_candidate_manifest(&a, err, sizeof(err)) &&
        prepare_design(&a, NULL, err, sizeof(err)) &&
        prepare_candidate_designs(&a, err, sizeof(err));
   if (ok && a.opt.validate &&
       (a.opt.covariance_map || a.opt.covariance_estimate))
      ok = load_inputs(&a, err, sizeof(err)) &&
           fit_covariance_groups(&a, 1, err, sizeof(err));
   if (ok && a.opt.validate) {
      const THD_picard_design_info *d = THD_picard_plan_info(a.plan);
      int j;
      printf("Rows=%d FullRows=%d Runs=%d Penalized=%d Nuisance=%d XRank=%d "
             "ZRank=%d\n",
             d->n, d->nfull, a.xm->nrun, d->p, d->k, d->xrank, d->zrank);
      if (a.manifest) {
         printf("CandidateDesigns=%d Policy=%s; scores are selection-only\n",
                a.manifest->count, a.opt.candidate_policy);
         for (j = 0; j < a.manifest->count; j++)
            printf("CandidateIndex=%d ID=%d XRank=%d Matrix=%s\n", j,
                   a.manifest->entry[j].id,
                   THD_picard_plan_info(a.candidate_designs[j])->xrank,
                   a.manifest->entry[j].matrix);
      }
      if (a.opt.single)
         printf("SingleTrial: supplied design/options only; pool masks, "
                "coarse model and learned PC rank require the data pass.\n");
      printf("Penalty=%s", a.opt.penalty);
      if (a.opt.elastic)
         printf(" L1Ratio=%.9g Standardize=%d KKT_tol=%.9g MaxIterations=%d",
                a.opt.solver.rho, a.opt.solver.standardize,
                a.opt.solver.tolerance, a.opt.solver.max_iterations);
      printf(" AR1=%.9g SupportRefit=%d", a.opt.ar1, a.opt.support_refit);
      if (a.opt.covariance_set)
         printf(" Covariance=%d a=%.9g b=%.9g ResetGaps=%d Map=%s",
                a.opt.covariance.model, a.opt.covariance.a,
                a.opt.covariance.b, a.opt.covariance.reset_gaps,
                a.opt.covariance_map ? a.opt.covariance_map : "none");
      if (a.opt.banded)
         printf(" BandGroups=%d PenaltyVectors=%d", a.opt.ngroup, a.opt.nband);
      printf("\n");
      if (a.predict) {
         printf("PredictionCV: Folds=%d Candidates=%d Nuisance=%s Rule=%s\n",
                a.folds->nfold, a.opt.ngrid,
                a.opt.project ? "project" : "train",
                a.opt.one_se ? "1se" : "min");
      }
      if (a.repeat) {
         const THD_picard_repeat_info *info =
             THD_picard_repeat_plan_info(a.repeat);
         printf("RepeatTuning: Folds=%d Candidates=%d ScoredTrials=%d "
                "Pairs=%.0f\n",
                info->nfold, info->ncandidate, info->nscored, info->npair);
      }
      for (j = 0; j < d->ncol; j++)
         printf("%d\t%s\t%s\n", j, a.roles[j] ? "penalized" : "nuisance",
                a.xm->collab[j]);
   } else if (ok) {
      /* Output preparation allocates in-memory buckets and metadata only.
         No files are written until every selected voxel has fitted safely. */
      ok = (a.opt.single ? load_inputs(&a, err, sizeof(err)) &&
                               single_prepare(&a, err, sizeof(err))
                         : 1) &&
           prepare_outputs(&a, argc, argv, err, sizeof(err)) &&
           (a.opt.single ? single_output(&a, err, sizeof(err))
                         : load_inputs(&a, err, sizeof(err))) &&
           ((a.opt.covariance_map || a.opt.covariance_estimate) ?
                fit_covariance_groups(&a, 0, err, sizeof(err)) :
                fit_voxels(&a, err, sizeof(err))) &&
           write_outputs(&a, err, sizeof(err));
      if (ok)
         INFO_message(PROGRAM ": fit %d voxels; wrote %s", a.nselected,
                      DSET_HEADNAME(a.out[OUT_MAIN]));
   }
   if (!ok)
      ERROR_message(PROGRAM ": %s", err);
   free_app(&a);
   return ok ? 0 : 1;
}
