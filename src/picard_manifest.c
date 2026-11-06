#include "picard_manifest.h"
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
static int read_line(FILE *f, char **line, size_t *capacity, char *err,
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
static int fields(char *line, char **values, int count)
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

static int integer(const char *s, int missing, int *value)
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
static char *resolve(const char *manifest, const char *field)
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

void picard_manifest_free(PicardManifest *m)
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
static int add_candidate(PicardManifest *m, char **v, int hrf, char *err,
                         int nerr)
{
   PicardCandidateInput *array, *e;
   int id, c;
   char *end;
   double sampling;
   if (!integer(v[0], 0, &id) || id > 16777215)
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
   e->matrix = resolve(m->path, v[1]);
   e->columns = resolve(m->path, v[2]);
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

PicardManifest *picard_manifest_read(const char *path, int hrf, char *err,
                                     int nerr)
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
   ok = read_line(f, &line, &capacity, err, nerr) == 1 &&
        !strcmp(line, "# 3dPicard candidates 1");
   if (!ok)
      manifest_error(err, nerr, "expected '# 3dPicard candidates 1'");
   while (ok && (status = read_line(f, &line, &capacity, err, nerr)) != 0) {
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
         ok = fields(line, v, 8);
         if (!ok)
            manifest_error(err, nerr, "candidate row needs 8 TSV fields");
         if (ok)
            ok = add_candidate(m, v, hrf, err, nerr);
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
static int column_row(PicardCandidateInput *e, const THD_xmat *xm,
                      const int *roles, char **v, char *err, int nerr)
{
   PicardColumnIdentity *ci;
   int column, role;
   if (!integer(v[0], 0, &column) || column >= xm->ncol)
      return manifest_error(err, nerr, "column identity index out of range");
   role = !strcmp(v[2], "penalized") ? 1 : !strcmp(v[2], "nuisance") ? 0 : -1;
   ci = e->identity + column;
   if (ci->label || strcmp(v[1], xm->collab[column]) || role != roles[column])
      return manifest_error(err, nerr,
                            "duplicate column or label/role mismatch "
                            "in candidate %d column %d",
                            e->id, column);
   if (!integer(v[3], 1, &ci->trial) || !integer(v[4], 1, &ci->condition) ||
       !integer(v[5], 1, &ci->run) || !integer(v[6], 0, &ci->basis) ||
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

int picard_manifest_read_columns(PicardManifest *m, int index,
                                 const THD_xmat *xm, const int *roles,
                                 const THD_picard_trial *trials, char *err,
                                 int nerr)
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
   while (ok && (status = read_line(f, &line, &capacity, err, nerr)) != 0) {
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
         ok = fields(line, v, 8);
         if (!ok)
            manifest_error(err, nerr, "column row needs 8 TSV fields");
         if (ok)
            ok = column_row(e, xm, roles, v, err, nerr);
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

static int same_text(const char *a, const char *b)
{
   return (!a && !b) || (a && b && !strcmp(a, b));
}

int picard_manifest_compare(const PicardManifest *m, int index,
                            const THD_xmat *a, const THD_xmat *b, char *err,
                            int nerr)
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
      if (!same_text(a->stim[j].label, b->stim[j].label) ||
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
