/* Picard's numerical workflows, collected in one implementation file.

   Reading map (search for the numbered section banners):
     1. Nuisance projection: remove/recover the supplied Z-column span.
     2. Fixed fitting and repeated-condition beta tuning.
     3. Held-out time-series prediction and global parameter selection.
     4. Candidate design/HRF selection and locked-choice refitting.
     5. Run-specific noise PCs and single-trial PC-count scoring.

   The common model is Y = X B + Z G + error. Only B is penalized. X is
   [time][task], Z is [time][nuisance], and Y is [time][response]. All arrays
   use compact row-major indexing; the active batch width q is the stride,
   even when the allocation can hold more responses.

   Plans own immutable design calculations. Workspaces borrow their plans
   and own mutable response scratch; a plan must outlive all its workspaces.
   Result views borrow scratch and expire on the next fit attempt. Only a
   successful fit publishes a valid result. Free functions tolerate NULL and
   partial construction, so failure paths unwind without goto statements.

   thd_regression.c supplies shared solvers and the existing EISPACK SVD
   adapter. thd_crossval supplies generic folds/path traversal. Dataset I/O,
   spatial masks, option policy and output transactions live in 3dPicard.c.
   These source-file boundaries do not change the statistical stage order.
*/
#define THD_PICARD_PRIVATE
#include "thd_picard.h"
#include "thd_regression_internal.h"
#include "thd_crossval.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ========================================================================
   1. NUISANCE PROJECTION
   ======================================================================== */

/* Remove or recover the part of a response explained by nuisance columns.

   For Z = U S V^T, retained columns of U form an orthonormal basis for Z's
   span. Projection removes U U^T Y; coefficient recovery uses V S^+ U^T Y.
   A singular-value cutoff ignores directions too small to estimate reliably.
   This avoids directly forming or inverting Z^T Z and supports redundant Z.

   The plan owns the decomposition; each worker owns its temporary U^T Y.
   n = rows, k = nuisance columns, q = response columns. Arrays are row-major.
   No file reading, intercept insertion, centering or sample repair occurs. */

struct THD_linproj {
   int n, k;
   THD_numeric_svd *svd;
};
struct THD_linproj_workspace {
   const THD_linproj *plan;
   int maxq;
   double *rot, *norm;
};
/* Use caller-owned error storage; return false to propagate the failure. */
static int linproj_error(char *err, int nerr, const char *msg)
{
   if (err && nerr > 0)
      snprintf(err, (size_t)nerr, "%s", msg);
   return 0;
}
/* Safe for NULL/partial construction; workers must be destroyed first. */
void THD_linproj_free(THD_linproj *p)
{
   if (p) {
      THD_numeric_svd_free(p->svd);
      free(p);
   }
}
/* Copy/factor Z once. k=0 is valid: projection then leaves Y unchanged and
   nuisance recovery has no coefficients to write. The caller keeps Z. */
/* Own an SVD of exactly the supplied nuisance design. With k=0 there is
   nothing to remove, so projection is identity and nuisance recovery is empty.
   Redundant nuisance columns are retained in the output coefficient mapping. */
THD_linproj *THD_linproj_prepare(const double *z, int n, int k, double tol,
                                 char *err, int nerr)
{
   THD_linproj *p;
   if (n < 1 || k < 0 || !isfinite(tol) || tol <= 0 || tol >= 1 || (k && !z)) {
      linproj_error(err, nerr, "invalid nuisance design");
      return NULL;
   }
   p = calloc(1, sizeof(*p));
   if (!p) {
      linproj_error(err, nerr, "cannot allocate nuisance projection");
      return NULL;
   }
   p->n = n;
   p->k = k;
   if (k) {
      p->svd = THD_numeric_svd_prepare(z, n, k, tol, err, nerr);
      if (!p->svd) {
         THD_linproj_free(p);
         return NULL;
      }
   }
   return p;
}
/* Count estimable nuisance directions, which may be fewer than columns. */
int THD_linproj_rank(const THD_linproj *p)
{
   return !p ? -1 : p->svd ? p->svd->rank : 0;
}
/* Release private scratch without freeing the borrowed decomposition. */
void THD_linproj_workspace_free(THD_linproj_workspace *w)
{
   if (w) {
      free(w->rot);
      free(w->norm);
      free(w);
   }
}
/* Allocate enough scratch for q responses. r is the number of singular
   directions stored; a no-nuisance plan needs only placeholder storage. */
THD_linproj_workspace *THD_linproj_workspace_alloc(const THD_linproj *p, int q,
                                                   char *err, int nerr)
{
   THD_linproj_workspace *w;
   size_t r = p && p->svd ? (size_t)p->svd->r : 1;
   if (!p || q < 1 || (size_t)q > SIZE_MAX / sizeof(double) / r) {
      linproj_error(err, nerr, "invalid projection workspace capacity");
      return NULL;
   }
   w = calloc(1, sizeof(*w));
   if (!w) {
      linproj_error(err, nerr, "cannot allocate projection workspace");
      return NULL;
   }
   w->plan = p;
   w->maxq = q;
   w->rot = calloc(r * q, sizeof(double));
   w->norm = calloc((size_t)q, sizeof(double));
   if (!w->rot || !w->norm) {
      THD_linproj_workspace_free(w);
      linproj_error(err, nerr, "cannot allocate projection buffers");
      return NULL;
   }
   return w;
}
/* Compute U^T Y: dot each retained basis direction with each response.
   Discarded directions get zero coordinates. Validate even when Z is absent
   so missing nuisance columns do not bypass the finite-data contract. */
/* Compute U^T Y once into workspace scratch, using only retained singular
   directions. Both projection and coefficient recovery start from this shared
   reduction, which is small when the nuisance span has low rank. */
static int linproj_rotate(THD_linproj_workspace *w, const double *y, int q,
                          char *err, int nerr)
{
   const THD_numeric_svd *s;
   int i, k, o;
   if (!w || !y || q < 1 || q > w->maxq)
      return linproj_error(err, nerr, "invalid projection inputs");
   for (i = 0; i < w->plan->n; i++)
      for (o = 0; o < q; o++)
         if (!isfinite(y[(size_t)i * q + o]))
            return linproj_error(err, nerr, "nonfinite nuisance response");
   s = w->plan->svd;
   if (!s)
      return 1;
   for (k = 0; k < s->r; k++)
      for (o = 0; o < q; o++) {
         double v = 0;
         if (s->s[k] > s->cutoff)
            for (i = 0; i < s->n; i++)
               v += s->u[(size_t)i * s->r + k] * y[(size_t)i * q + o];
         if (!isfinite(v))
            return linproj_error(err, nerr,
                                 "nuisance projection exceeds double range");
         w->rot[(size_t)k * q + o] = v;
      }
   return 1;
}
/* Write the component orthogonal to Z into out; in-place out==y is allowed.
   The first subtraction removes the nuisance component. The second pass
   removes numerical residue caused by floating-point rounding in the first. */
int THD_linproj_project(THD_linproj_workspace *w, const double *y, int q,
                        double *out, char *err, int nerr)
{
   const THD_numeric_svd *s;
   int i, k, o, pass;
   if (!out || !linproj_rotate(w, y, q, err, nerr))
      return linproj_error(err, nerr,
                           "invalid/nonfinite nuisance projection inputs");
   for (o = 0; o < q; o++) {
      double norm = 0;
      for (i = 0; i < w->plan->n; i++)
         norm = hypot(norm, y[(size_t)i * q + o]);
      w->norm[o] = norm;
      if (!isfinite(norm))
         return linproj_error(err, nerr,
                              "projection response norm exceeds double range");
   }
   if (out != y)
      memcpy(out, y, sizeof(double) * (size_t)w->plan->n * q);
   s = w->plan->svd;
   if (!s)
      return 1;
   for (pass = 0; pass < 2; pass++) {
      if (pass && !linproj_rotate(w, out, q, err, nerr))
         return 0;
      for (i = 0; i < s->n; i++)
         for (o = 0; o < q; o++) {
            double v = 0;
            for (k = 0; k < s->r; k++)
               v += s->u[(size_t)i * s->r + k] * w->rot[(size_t)k * q + o];
            out[(size_t)i * q + o] -= v;
            if (!isfinite(out[(size_t)i * q + o]))
               return linproj_error(err, nerr,
                                    "projection result exceeds double range");
         }
   }
   /* A response entirely in Z's span should project to zero. Roundoff can
      leave a tiny vector that a later solver could amplify. Only zero values
      below a dimension-scaled machine-precision bound, not a model threshold.
    */
   for (o = 0; o < q; o++) {
      double norm = 0;
      for (i = 0; i < s->n; i++)
         norm = hypot(norm, out[(size_t)i * q + o]);
      if (norm <= 64 * DBL_EPSILON * (s->n > s->p ? s->n : s->p) * w->norm[o])
         for (i = 0; i < s->n; i++)
            out[(size_t)i * q + o] = 0;
   }
   return 1;
}
/* Recover minimum-norm coefficients G = V S^+ U^T Y. Divide only by
   retained singular values; zeroing the others avoids unstable division.
   gamma is a separate k-by-q output array, not the n-by-q response buffer. */
int THD_linproj_solve(THD_linproj_workspace *w, const double *y, int q,
                      double *gamma, char *err, int nerr)
{
   const THD_numeric_svd *s;
   int j, k, o;
   if (!linproj_rotate(w, y, q, err, nerr))
      return 0;
   s = w->plan->svd;
   if (!s)
      return 1;
   if (!gamma || gamma == y)
      return linproj_error(
          err, nerr, "nuisance recovery requires separate coefficient output");
   for (j = 0; j < s->p; j++)
      for (o = 0; o < q; o++) {
         double v = 0;
         for (k = 0; k < s->r; k++)
            if (s->s[k] > s->cutoff)
               v += s->v[(size_t)j * s->r + k] *
                    (w->rot[(size_t)k * q + o] / s->s[k]);
         if (!isfinite(v))
            return linproj_error(err, nerr,
                                 "nuisance coefficients exceed double range");
         gamma[(size_t)j * q + o] = v;
      }
   return 1;
}

/* ========================================================================
   2. FIXED FITS AND REPEATED-CONDITION TUNING
   ======================================================================== */

/* Numerical workflow shared by the CLI's fixed and tuned fits.

   Let Y be n rows of time samples by q voxel responses, X be n by p penalized
   predictors, and Z be n by k nuisance predictors. Fit Y = X B + Z G + error.
   Only B is penalized. Removing Z's span from both X and Y lets the regression
   solver estimate B without penalizing nuisance effects; afterward we recover
   G from Y-XB. With fixed AR(1), apply W to X, Z and Y first; recover G
   using W(Y-XB), then reconstruct original-unit outputs. See picard_prepare()
   and picard_fit_model().

   Plans own design-dependent arrays and are shared read-only. Workspaces
   borrow a plan and own all changing batch storage. Every numerical array is
   row-major; element (row, column) has offset row * width + column. */

/* Mutable state for one active voxel batch. The plan is borrowed; coefficient,
   projected-response and fitted/residual arrays are owned here. Alternative
   solver workspaces implement the same result layout, keeping callers separate
   from estimator-specific storage and convergence diagnostics. */
struct THD_picard_workspace {
   const THD_picard_plan *plan;
   int maxq, valid;
   THD_linproj_workspace *projection;
   THD_regression_workspace *regression;
   THD_elastic_workspace *elastic;
   THD_banded_workspace *banded;
   double *whitened, *support_x, *support_y, *support_beta;
   int *support_size, *support_rank;
   double *projected, *gamma, *fitted, *residual;
   THD_picard_result result;
};
/* Returning false lets callers report an error without allocating a string.
   The caller owns err; a NULL buffer simply disables the explanatory text. */
static int picard_error(char *err, int nerr, const char *msg)
{
   if (err && nerr > 0)
      snprintf(err, (size_t)nerr, "%s", msg);
   return 0;
}
/* Check multiplication before allocation to prevent integer wraparound.
   A missing nuisance dimension still gets a harmless placeholder allocation,
   so callers can use one allocation-failure check for all their buffers. */
static double *picard_array(size_t n, size_t p)
{
   if (!p)
      p = 1;
   return n > SIZE_MAX / sizeof(double) / p ? NULL
                                            : calloc(n * p, sizeof(double));
}
/* Accept NULL and partially constructed plans, including failed preparation.
   Workspaces must be freed first because they borrow the plan's factors. */
void THD_picard_plan_free(THD_picard_plan *p)
{
   if (!p)
      return;
   free(p->pen);
   free(p->nuis);
   free(p->rows);
   free(p->runs);
   free(p->x);
   free(p->z);
   free(p->effective);
   free(p->weighted_z);
   free(p->covariance_lower);
   free(p->groups);
   free(p->penalties);
   THD_banded_design_free(p->banded);
   THD_linproj_free(p->projection);
   THD_regression_design_free(p->design);
   THD_elastic_design_free(p->elastic);
   free(p);
}
/* Build the retained-observation correlation, then its lower Cholesky
   factor. Unit marginal variance fixes penalty units. No missing responses
   enter this operation. Dense storage is deliberate for the initial exact
   implementation; adapter budgets include every full/fold/candidate factor. */
static int picard_covariance_prepare(THD_picard_plan *p, char *err, int nerr)
{
   int i, j, k, start = 0, n = p->info.n;
   const THD_picard_covariance *c = &p->info.covariance;
   double first = (c->a+c->b)*(1+c->a*c->b) /
                  (1+2*c->a*c->b+c->b*c->b);
   if (!p->info.covariance_enabled || c->model == 0)
      return 1;
   p->covariance_lower = picard_array(n, n);
   if (!p->covariance_lower)
      return picard_error(err, nerr, "cannot allocate covariance factor");
   for (i = 0; i < n; i++) {
      if (i && (p->runs[i] != p->runs[i-1] ||
                (c->reset_gaps && p->rows[i] != p->rows[i-1]+1)))
         start = i;
      for (j = start; j <= i; j++) {
         double value = i == j ? 1 :
             first * pow(c->a, p->rows[i]-p->rows[j]-1);
         for (k = start; k < j; k++)
            value -= p->covariance_lower[(size_t)i*n+k] *
                     p->covariance_lower[(size_t)j*n+k];
         if (!isfinite(value) || (i == j && value <= 0))
            return picard_error(err, nerr,
                               "covariance is not positive definite");
         p->covariance_lower[(size_t)i*n+j] = i == j ? sqrt(value) :
             value / p->covariance_lower[(size_t)j*n+j];
      }
   }
   return 1;
}

/* Apply an externally fixed stationary AR(1) transform. A fresh segment
   starts at each run/censor/fold gap; its first value gets sqrt(1-rho^2).
   Contiguous values use y[t]-rho*y[t-1]. Reverse traversal permits in-place
   use without overwriting the previous sample before it is read. Never read
   across missing original rows, even if retained rows are adjacent in RAM. */
static void picard_whiten(const THD_picard_plan *p, const double *input, int q,
                          double *output)
{
   int i, o;
   double rho = p->info.ar1, first = sqrt(1 - rho * rho);
   if (p->covariance_lower) {
      int j, n = p->info.n;
      for (i = 0; i < n; i++)
         for (o = 0; o < q; o++) {
            double value = input[(size_t)i*q+o];
            for (j = 0; j < i; j++)
               value -= p->covariance_lower[(size_t)i*n+j] *
                        output[(size_t)j*q+o];
            output[(size_t)i*q+o] =
                value / p->covariance_lower[(size_t)i*n+i];
         }
      return;
   }
   for (i = p->info.n - 1; i >= 0; i--) {
      int adjacent =
          i && p->runs[i] == p->runs[i - 1] && p->rows[i] == p->rows[i - 1] + 1;
      for (o = 0; o < q; o++)
         output[(size_t)i * q + o] =
             adjacent ? input[(size_t)i * q + o] -
                            rho * input[(size_t)(i - 1) * q + o]
                      : first * input[(size_t)i * q + o];
   }
}

/* Profile Gaussian restricted likelihood on a FIXED mean-model span.
   With orthonormal pilot basis B, -2 log L (up to common constants) is
   log|R| + log|B'R^-1 B| + (n-r)*log(RSS/(n-r)). Using B rather than raw
   redundant columns makes rank handling and candidate comparisons coherent.
   This is ordinary pilot REML, not joint estimation of ridge shrinkage. */
int THD_picard_reml(const double *matrix, int n, int d,
                    const int *rows, const int *runs, const double *y, int q,
                    int model, int reset_gaps, double tol,
                    THD_picard_reml_result *result, char *err, int nerr)
{
   THD_numeric_svd *pilot = NULL;
   double *basis = NULL, *weighted = NULL, *wy = NULL, *residual = NULL;
   double *second = NULL, *worst = NULL;
   int rank = 0, i, j, c, o, ai, bi, ok = 1;
   int bound = model == 0 ? 0 : 9;
   if (!matrix || !rows || !runs || !y || !result || n < 2 || d < 1 ||
       q < 1 || model < 0 || model > 2 ||
       (reset_gaps != 0 && reset_gaps != 1))
      return picard_error(err, nerr, "invalid REML inputs");
   for (i = 0; i < n; i++) {
      if (rows[i] < 0 || runs[i] < 0 ||
          (i && (rows[i] <= rows[i-1] || runs[i] < runs[i-1])))
         return picard_error(err, nerr, "invalid REML row/run identities");
      for (o = 0; o < q; o++)
         if (!isfinite(y[(size_t)i*q+o]))
            return picard_error(err, nerr, "nonfinite REML response");
   }
   pilot = THD_numeric_svd_prepare(matrix, n, d, tol, err, nerr);
   if (!pilot) return 0;
   rank = pilot->rank;
   if (rank >= n) {
      THD_numeric_svd_free(pilot);
      return picard_error(err, nerr, "saturated REML pilot: no residual df");
   }
   basis = picard_array(n, rank);
   weighted = picard_array(n, rank);
   wy = picard_array(n, q);
   residual = picard_array(n, q);
   second = picard_array(q, 1);
   worst = picard_array(q, 1);
   ok = basis && weighted && wy && residual && second && worst;
   if (!ok) picard_error(err, nerr, "cannot allocate REML scratch");
   for (j = c = 0; ok && j < pilot->r; j++)
      if (pilot->s[j] > pilot->cutoff) {
         for (i = 0; i < n; i++)
            basis[(size_t)i*rank+c] = pilot->u[(size_t)i*pilot->r+j];
         c++;
      }
   THD_numeric_svd_free(pilot);
   for (o = 0; ok && o < q; o++) {
      result[o] = (THD_picard_reml_result){0};
      result[o].score = second[o] = INFINITY;
      worst[o] = -INFINITY;
      result[o].rank = rank;
   }
   for (ai = -bound; ok && ai <= bound; ai++)
      for (bi = model == 2 ? -9 : 0;
           ok && bi <= (model == 2 ? 9 : 0); bi++) {
         THD_picard_plan transform = {0};
         THD_linproj *projection = NULL;
         THD_linproj_workspace *work = NULL;
         double determinant = 0;
         transform.info.n = n;
         transform.info.covariance_enabled = 1;
         transform.info.covariance = (THD_picard_covariance){
             model, reset_gaps, ai/10.0, bi/10.0};
         transform.rows = (int *)rows;
         transform.runs = (int *)runs;
         ok = picard_covariance_prepare(&transform, err, nerr);
         if (ok) {
            picard_whiten(&transform, basis, rank, weighted);
            picard_whiten(&transform, y, q, wy);
            projection = THD_linproj_prepare(weighted, n, rank, tol,
                                             err, nerr);
            work = THD_linproj_workspace_alloc(projection, q, err, nerr);
            ok = work && THD_linproj_rank(projection) == rank &&
                 THD_linproj_project(work, wy, q, residual, err, nerr);
            if (!ok) picard_error(err, nerr,
                                 "REML weighted pilot lost rank");
         }
         if (ok) {
            if (transform.covariance_lower)
               for (i = 0; i < n; i++)
                  determinant += 2*log(transform.covariance_lower[
                                           (size_t)i*n+i]);
            if (rank)
               for (j = 0; j < projection->svd->r; j++)
                  determinant += 2*log(projection->svd->s[j]);
            for (o = 0; o < q; o++) {
               double rss = 0, score, variance;
               for (i = 0; i < n; i++)
                  rss += residual[(size_t)i*q+o]*residual[(size_t)i*q+o];
               variance = rss/(n-rank);
               if (!(variance > 0) || !isfinite(variance)) {
                  ok = picard_error(err, nerr,
                                    "REML response has zero/invalid "
                                    "residual variance");
                  break;
               }
               score = determinant+(n-rank)*log(variance);
               if (score < result[o].score -
                       1e-10*(1+fabs(score)) ||
                   (fabs(score-result[o].score) <= 1e-10*(1+fabs(score)) &&
                    (abs(ai)+abs(bi))/10.0 <
                        fabs(result[o].a)+fabs(result[o].b))) {
                  second[o] = result[o].score;
                  result[o].score = score;
                  result[o].sigma2 = variance;
                  result[o].a = ai/10.0;
                  result[o].b = bi/10.0;
                  result[o].status = (bound && abs(ai) == 9) ||
                                      (model == 2 && abs(bi) == 9);
               } else if (score < second[o]) second[o] = score;
               if (score > worst[o]) worst[o] = score;
            }
         }
         THD_linproj_workspace_free(work);
         THD_linproj_free(projection);
         free(transform.covariance_lower);
      }
   for (o = 0; ok && o < q; o++) {
      result[o].separation = isfinite(second[o]) ?
          fmax(0, second[o]-result[o].score) : 0;
      if (bound && result[o].separation <=
          1e-10*(1+fabs(result[o].score))) result[o].status = 4;
      if (bound && worst[o]-result[o].score <=
          1e-10*(1+fabs(result[o].score))) result[o].status = 3;
   }
   free(basis);
   free(weighted);
   free(wy);
   free(residual);
   free(second);
   free(worst);
   return ok;
}

/* Whiten X and Z together, then remove the whitened nuisance span from X.
   Keep original X/Z for output reconstruction and fold creation. The effective
   X is shared by regression preparation and optional active-support refits. */
static int picard_prepare(THD_picard_plan *p, char *err, int nerr)
{
   int n = p->info.n, np = p->info.p;
   double *z = picard_array(n, p->info.k);
   THD_linproj_workspace *w;
   p->effective = picard_array(n, np);
   if (!z || !p->effective) {
      free(z);
      return picard_error(err, nerr, "cannot allocate effective design");
   }
   picard_whiten(p, p->x, np, p->effective);
   picard_whiten(p, p->z, p->info.k, z);
   p->projection =
       THD_linproj_prepare(z, n, p->info.k, p->info.rank_tol, err, nerr);
   p->weighted_z = z;
   if (!p->projection)
      return 0;
   w = THD_linproj_workspace_alloc(p->projection, np, err, nerr);
   if (w && THD_linproj_project(w, p->effective, np, p->effective, err, nerr))
      p->design = THD_regression_design_prepare(p->effective, n, np,
                                                p->info.rank_tol, err, nerr);
   THD_linproj_workspace_free(w);
   if (!p->design)
      return 0;
   p->info.xrank = THD_regression_design_rank(p->design);
   p->info.zrank = THD_linproj_rank(p->projection);
   return 1;
}
/* Validate the design and copy it into owned X/Z arrays. roles partitions
   columns; rows/runs preserve acquisition identities after censoring. The
   pen/nuis maps retain original column numbers for coefficient labels.
   Input arrays can be released after this function returns successfully. */
/* Validate and copy column roles plus original row/run identities before
   constructing factors. Split the supplied matrix into X and Z without
   reordering their internal column identities; those maps later restore labels
   and nuisance coefficients to their original design columns. */
static THD_picard_plan *picard_create_base(const double *matrix, int n, int d,
                                           const int *roles, const int *rows,
                                           const int *runs, int nfull,
                                           double tol, double ar1,
                                           const THD_picard_covariance *cov,
                                           char *err,
                                           int nerr)
{
   THD_picard_plan *p;
   int i, j, np = 0, nk = 0;
   if (!matrix || !roles || n < 1 || d < 1 || nfull < n || !isfinite(tol) ||
       tol <= 0 || tol >= 1 ||
       (size_t)n > SIZE_MAX / sizeof(double) / (size_t)d) {
      picard_error(err, nerr, "invalid Picard design inputs");
      return NULL;
   }
   for (j = 0; j < d; j++) {
      if (roles[j] == THD_PICARD_PENALIZED)
         np++;
      else if (roles[j] == THD_PICARD_NUISANCE)
         nk++;
      else {
         picard_error(err, nerr,
                      "every column must have exactly one explicit role");
         return NULL;
      }
   }
   if (!np) {
      picard_error(err, nerr, "Picard needs at least one penalized column");
      return NULL;
   }
   for (i = 0; i < n; i++) {
      if ((rows && (rows[i] < 0 || rows[i] >= nfull ||
                    (i && rows[i] <= rows[i - 1]))) ||
          (runs && (runs[i] < 0 || (i && runs[i] < runs[i - 1])))) {
         picard_error(err, nerr, "invalid retained row/run identities");
         return NULL;
      }
      for (j = 0; j < d; j++)
         if (!isfinite(matrix[(size_t)i * d + j])) {
            picard_error(err, nerr, "nonfinite Picard design");
            return NULL;
         }
   }
   p = calloc(1, sizeof(*p));
   if (!p) {
      picard_error(err, nerr, "cannot allocate Picard plan");
      return NULL;
   }
   p->info.n = n;
   p->info.ncol = d;
   p->info.p = np;
   p->info.k = nk;
   p->info.nfull = nfull;
   p->info.rank_tol = tol;
   p->info.ar1 = ar1;
   p->info.covariance_enabled = cov != NULL;
   if (cov)
      p->info.covariance = *cov;
   p->pen = calloc((size_t)np, sizeof(int));
   p->nuis = calloc((size_t)(nk ? nk : 1), sizeof(int));
   p->rows = calloc((size_t)n, sizeof(int));
   p->runs = calloc((size_t)n, sizeof(int));
   p->x = picard_array(n, np);
   p->z = picard_array(n, nk);
   if (!p->pen || !p->nuis || !p->rows || !p->runs || !p->x || !p->z) {
      THD_picard_plan_free(p);
      picard_error(err, nerr, "cannot allocate Picard design buffers");
      return NULL;
   }
   /* Split columns without changing their relative order within each role.
      col is a position in X or Z, while j is a position in the input matrix. */
   np = nk = 0;
   for (j = 0; j < d; j++) {
      int col = roles[j] ? np++ : nk++;
      if (roles[j])
         p->pen[col] = j;
      else
         p->nuis[col] = j;
      for (i = 0; i < n; i++) {
         if (roles[j])
            p->x[(size_t)i * p->info.p + col] = matrix[(size_t)i * d + j];
         else
            p->z[(size_t)i * p->info.k + col] = matrix[(size_t)i * d + j];
      }
   }
   for (i = 0; i < n; i++) {
      p->rows[i] = rows ? rows[i] : i;
      p->runs[i] = runs ? runs[i] : 0;
   }
   p->info.penalized_columns = p->pen;
   p->info.nuisance_columns = p->nuis;
   p->info.row_ids = p->rows;
   p->info.run_ids = p->runs;
   if (!picard_covariance_prepare(p, err, nerr) ||
       !picard_prepare(p, err, nerr)) {
      THD_picard_plan_free(p);
      return NULL;
   }
   return p;
}
/* One construction path guarantees whitening precedes nuisance projection
   and every penalty sees the same effective design. The legacy convenience
   constructors below delegate here with identity whitening. */
THD_picard_plan *THD_picard_plan_create_config(const double *matrix, int n,
                                               int d, const int *roles,
                                               const int *rows, const int *runs,
                                               int nfull, double tol,
                                               const THD_picard_config *cfg,
                                               char *err, int nerr)
{
   THD_picard_plan *p;
   THD_picard_config identity = {0};
   if (!cfg)
      cfg = &identity;
   if (!isfinite(cfg->ar1) || fabs(cfg->ar1) >= 1 ||
       (cfg->elastic && cfg->banded) ||
       (cfg->support_refit != 0 && cfg->support_refit != 1) ||
       (cfg->support_refit && (!cfg->elastic || cfg->elastic->rho <= 0))) {
      picard_error(err, nerr,
                   "invalid whitening/support/penalty configuration");
      return NULL;
   }
   if (cfg->covariance) {
      const THD_picard_covariance *c = cfg->covariance;
      if (cfg->ar1 != 0 || c->model < 0 || c->model > 2 ||
          (c->reset_gaps != 0 && c->reset_gaps != 1) ||
          !isfinite(c->a) || !isfinite(c->b) ||
          fabs(c->a) >= 1 || fabs(c->b) >= 1 ||
          (c->model == 0 && c->a != 0) || (c->model < 2 && c->b != 0)) {
         picard_error(err, nerr, "invalid stationary covariance");
         return NULL;
      }
   }
   p = picard_create_base(matrix, n, d, roles, rows, runs, nfull, tol,
                          cfg->ar1, cfg->covariance, err, nerr);
   if (!p)
      return NULL;
   p->info.support_refit = cfg->support_refit;
   if (cfg->elastic) {
      p->elastic = THD_elastic_prepare(p->effective, n, p->info.p, cfg->elastic,
                                       err, nerr);
      if (!p->elastic) {
         THD_picard_plan_free(p);
         return NULL;
      }
      p->config = *cfg->elastic;
      p->info.elastic = 1;
   }
   if (cfg->banded) {
      const THD_banded_config *bc = cfg->banded;
      p->banded = THD_banded_prepare(p->effective, n, p->info.p, bc, err, nerr);
      if (!p->banded) {
         THD_picard_plan_free(p);
         return NULL;
      }
      p->groups = calloc((size_t)p->info.p, sizeof(int));
      p->penalties = picard_array(bc->ncandidate, bc->ngroup);
      if (!p->groups || !p->penalties) {
         THD_picard_plan_free(p);
         picard_error(err, nerr, "cannot copy banded settings");
         return NULL;
      }
      memcpy(p->groups, bc->group, sizeof(int) * (size_t)p->info.p);
      memcpy(p->penalties, bc->penalty,
             sizeof(double) * (size_t)bc->ncandidate * bc->ngroup);
      p->band_config = *bc;
      p->band_config.group = p->groups;
      p->band_config.penalty = p->penalties;
      p->info.banded = 1;
      p->info.ngroup = bc->ngroup;
      p->info.ncandidate = bc->ncandidate;
   }
   return p;
}
THD_picard_plan *THD_picard_plan_create(const double *matrix, int n, int d,
                                        const int *roles, const int *rows,
                                        const int *runs, int nfull, double tol,
                                        char *err, int nerr)
{
   return THD_picard_plan_create_config(matrix, n, d, roles, rows, runs, nfull,
                                        tol, NULL, err, nerr);
}
THD_picard_plan *THD_picard_plan_create_elastic(
    const double *matrix, int n, int d, const int *roles, const int *rows,
    const int *runs, int nfull, double tol, const THD_elastic_config *config,
    char *err, int nerr)
{
   THD_picard_config cfg = {config, NULL, 0, 0, NULL};
   if (!config) {
      picard_error(err, nerr, "missing elastic-net configuration");
      return NULL;
   }
   return THD_picard_plan_create_config(matrix, n, d, roles, rows, runs, nfull,
                                        tol, &cfg, err, nerr);
}

/* Preserve estimator configuration while rebuilding all training factors. */
THD_picard_plan *THD_picard_plan_like(const THD_picard_plan *full,
                                      const double *matrix, int n, int d,
                                      const int *roles, const int *rows,
                                      const int *penmap, const int *runs,
                                      int nfull, double tol, char *err,
                                      int nerr)
{
   THD_picard_config cfg = {full->elastic ? &full->config : NULL, NULL,
                            full->info.ar1, full->info.support_refit,
                            full->info.covariance_enabled ?
                                &full->info.covariance : NULL};
   THD_banded_config bc = full->band_config;
   THD_picard_plan *result;
   int *groups = NULL, j, np = 0;
   if (full->banded) {
      for (j = 0; j < d; j++)
         np += roles[j] != 0;
      groups = calloc((size_t)np, sizeof(int));
      if (!groups) {
         picard_error(err, nerr, "cannot map training groups");
         return NULL;
      }
      for (j = 0; j < np; j++)
         groups[j] = full->groups[penmap ? penmap[j] : j];
      bc.group = groups;
      cfg.banded = &bc;
   }
   result = THD_picard_plan_create_config(matrix, n, d, roles, rows, runs,
                                          nfull, tol, &cfg, err, nerr);
   free(groups);
   return result;
}

/* Sparse grids run strongest to weakest to enable warm coordinate paths. */
int THD_picard_grid_valid(const THD_picard_plan *p, const double *grid,
                          int count, char *err, int nerr)
{
   int i, j;
   if (p->banded) {
      for (i = 0; i < count; i++) {
         if (!isfinite(grid[i]) || grid[i] < 0 ||
             grid[i] >= p->info.ncandidate || grid[i] != floor(grid[i]))
            return picard_error(err, nerr, "invalid banded candidate index");
         for (j = 0; j < i; j++)
            if (grid[i] == grid[j])
               return picard_error(err, nerr,
                                   "duplicate banded candidate index");
      }
      return 1;
   }
   for (i = 0; i < count; i++) {
      if (!isfinite(grid[i]) ||
          (p->elastic ? grid[i] <= 0 : (grid[i] < 0 || grid[i] > 1)))
         return picard_error(err, nerr, "invalid fraction/lambda grid value");
      if (p->elastic && i && grid[i] >= grid[i - 1])
         return picard_error(err, nerr, "lambda grid must strictly decrease");
      for (j = 0; j < i; j++)
         if (grid[i] == grid[j])
            return picard_error(err, nerr, "duplicate tuning parameter");
   }
   return 1;
}

double THD_picard_strength(const THD_picard_plan *p, double parameter)
{
   if (p->banded) {
      int g;
      double strength = 0;
      for (g = 0; g < p->info.ngroup; g++)
         strength += p->penalties[(size_t)(int)parameter * p->info.ngroup + g] /
                     p->info.ngroup;
      return strength;
   }
   return p->elastic ? parameter : -parameter;
}

const THD_picard_design_info *THD_picard_plan_info(const THD_picard_plan *p)
{
   return p ? &p->info : NULL;
}
/* Free only this worker's mutable storage; its borrowed plan survives. */
void THD_picard_workspace_free(THD_picard_workspace *w)
{
   if (!w)
      return;
   THD_linproj_workspace_free(w->projection);
   THD_regression_workspace_free(w->regression);
   THD_elastic_workspace_free(w->elastic);
   THD_banded_workspace_free(w->banded);
   free(w->whitened);
   free(w->support_x);
   free(w->support_y);
   free(w->support_beta);
   free(w->support_size);
   free(w->support_rank);
   free(w->projected);
   free(w->gamma);
   free(w->fitted);
   free(w->residual);
   free(w);
}
/* Reserve capacity for up to q simultaneous voxel responses. Fit calls may
   use fewer columns, including the final partial batch. Nothing allocated
   here is shared between workers, even though their plan may be identical. */
THD_picard_workspace *THD_picard_workspace_alloc(const THD_picard_plan *p,
                                                 int q, char *err, int nerr)
{
   THD_picard_workspace *w;
   if (!p || q < 1) {
      picard_error(err, nerr, "invalid Picard worker capacity");
      return NULL;
   }
   w = calloc(1, sizeof(*w));
   if (!w) {
      picard_error(err, nerr, "cannot allocate Picard worker");
      return NULL;
   }
   w->plan = p;
   w->maxq = q;
   w->projection = THD_linproj_workspace_alloc(p->projection, q, err, nerr);
   w->regression = THD_regression_workspace_alloc(p->design, q, err, nerr);
   if (p->elastic)
      w->elastic = THD_elastic_workspace_alloc(p->elastic, q, err, nerr);
   if (p->banded)
      w->banded = THD_banded_workspace_alloc(p->banded, q, err, nerr);
   w->whitened = picard_array(p->info.n, q);
   if (p->info.support_refit) {
      w->support_x = picard_array(p->info.n, p->info.p);
      w->support_y = picard_array(p->info.n, 1);
      w->support_beta = picard_array(p->info.p, q);
      w->support_size = calloc((size_t)q, sizeof(int));
      w->support_rank = calloc((size_t)q, sizeof(int));
   }
   w->projected = picard_array(p->info.n, q);
   w->fitted = picard_array(p->info.n, q);
   w->residual = picard_array(p->info.n, q);
   w->gamma = picard_array(p->info.k ? p->info.k : 1, q);
   if (!w->whitened || (p->banded && !w->banded) ||
       (p->info.support_refit &&
        (!w->support_x || !w->support_y || !w->support_beta ||
         !w->support_size || !w->support_rank)) ||
       (p->elastic && !w->elastic) || !w->projection || !w->regression ||
       !w->projected || !w->fitted || !w->residual || !w->gamma) {
      THD_picard_workspace_free(w);
      picard_error(err, nerr, "cannot allocate Picard worker buffers");
      return NULL;
   }
   return w;
}
/* Re-estimate only the sparse solution's exact nonzero support with OLS.
   Selection occurred in this same training plan; CV never borrows full-data
   support. Factor at most one n-by-support design at a time. Empty support
   yields zero B; nuisance recovery below still runs normally. Diagnostics
   from the sparse solver remain its pre-refit diagnostics, not OLS scores. */
static int picard_refit_support(THD_picard_workspace *w, const double *beta,
                                int q, char *err, int nerr)
{
   const THD_picard_plan *p = w->plan;
   int i, j, o, at, ns, ok = 1;
   memset(w->support_beta, 0, sizeof(double) * (size_t)p->info.p * q);
   for (o = 0; o < q && ok; o++) {
      THD_regression_design *design;
      THD_regression_workspace *work;
      const THD_fraction_result *r;
      ns = 0;
      for (j = 0; j < p->info.p; j++)
         ns += beta[(size_t)j * q + o] != 0;
      w->support_size[o] = ns;
      w->support_rank[o] = 0;
      if (!ns)
         continue;
      for (i = 0; i < p->info.n; i++) {
         at = 0;
         for (j = 0; j < p->info.p; j++)
            if (beta[(size_t)j * q + o] != 0)
               w->support_x[(size_t)i * ns + at++] =
                   p->effective[(size_t)i * p->info.p + j];
         w->support_y[i] = w->projected[(size_t)i * q + o];
      }
      design = THD_regression_design_prepare(w->support_x, p->info.n, ns,
                                             p->info.rank_tol, err, nerr);
      work =
          design ? THD_regression_workspace_alloc(design, 1, err, nerr) : NULL;
      ok = work &&
           THD_regression_fit_fraction(work, w->support_y, 1, 1, err, nerr);
      if (ok) {
         r = THD_regression_fraction_result(work);
         at = 0;
         for (j = 0; j < p->info.p; j++)
            if (beta[(size_t)j * q + o] != 0)
               w->support_beta[(size_t)j * q + o] = r->coef[at++];
         w->support_rank[o] = THD_regression_design_rank(design);
      }
      THD_regression_workspace_free(work);
      THD_regression_design_free(design);
   }
   return ok;
}

/* Shared reconstruction keeps nuisance recovery identical across estimators.
   OLS uses the original minimum-norm reference, even for a sparse plan. */
static int picard_fit_model(THD_picard_workspace *w, const double *y, int q,
                            const double *parameters, int warm, int ols,
                            char *err, int nerr)
{
   const THD_picard_plan *p;
   const THD_fraction_result *r = NULL;
   const THD_elastic_result *e = NULL;
   const THD_banded_result *b = NULL;
   const double *beta;
   int i, j, o, was_valid = w && w->valid;
   /* Invalidate old results before any operation can fail. Otherwise a
      caller could mistake a previous batch's coefficients for this batch. */
   if (w)
      w->valid = 0;
   if (!w || !y || (!parameters && !ols) || q < 1 || q > w->maxq)
      return picard_error(err, nerr, "invalid Picard response batch");
   p = w->plan;
   if (warm && p->elastic && (!was_valid || !w->result.elastic_diagnostic))
      return picard_error(err, nerr,
                          "Picard warm step needs a preceding sparse fit");
   /* Project Y using the same nuisance basis used for X in prepare().
      This is not automatic demeaning: an intercept is removed only when
      the supplied nuisance columns span a constant vector. */
   picard_whiten(p, y, q, w->whitened);
   if (!THD_linproj_project(w->projection, w->whitened, q, w->projected, err,
                            nerr))
      return 0;
   if (p->banded && !ols) {
      if (!THD_banded_fit(w->banded, w->projected, q, parameters, err, nerr))
         return 0;
      b = THD_banded_get_result(w->banded);
      beta = b->coef;
   } else if (p->elastic && !ols) {
      if (!THD_elastic_fit(w->elastic, w->projected, q, parameters, warm, err,
                           nerr))
         return 0;
      e = THD_elastic_get_result(w->elastic);
      beta = e->coef;
   } else {
      int ok = ols ? THD_regression_fit_fraction(w->regression, w->projected, q,
                                                 1, err, nerr)
                   : THD_regression_fit_fractions(w->regression, w->projected,
                                                  q, parameters, err, nerr);
      if (!ok)
         return 0;
      r = THD_regression_fraction_result(w->regression);
      beta = r->coef;
   }
   if (p->info.support_refit && !ols) {
      if (!picard_refit_support(w, beta, q, err, nerr))
         return 0;
      beta = w->support_beta;
   }
   /* Reconstruct X B using the original X, not the projected X. Temporarily
      store Y-XB in residual so the next step can solve for nuisance effects. */
   for (i = 0; i < p->info.n; i++)
      for (o = 0; o < q; o++) {
         double v = 0;
         for (j = 0; j < p->info.p; j++)
            v += p->x[(size_t)i * p->info.p + j] * beta[(size_t)j * q + o];
         w->fitted[(size_t)i * q + o] = v;
         w->residual[(size_t)i * q + o] = y[(size_t)i * q + o] - v;
      }
   /* Solve G = (WZ)^+ W(Y-XB); + is the minimum-norm pseudoinverse.
      This also works with duplicate nuisance columns. Add ZG back to get
      fitted = XB+ZG and residual = Y-fitted in original response units. */
   picard_whiten(p, w->residual, q, w->whitened);
   if (!THD_linproj_solve(w->projection, w->whitened, q, w->gamma, err, nerr))
      return 0;
   for (i = 0; i < p->info.n; i++)
      for (o = 0; o < q; o++) {
         double v = w->fitted[(size_t)i * q + o];
         for (j = 0; j < p->info.k; j++)
            v += p->z[(size_t)i * p->info.k + j] * w->gamma[(size_t)j * q + o];
         if (!isfinite(v) || !isfinite(y[(size_t)i * q + o] - v))
            return picard_error(err, nerr,
                                "Picard prediction exceeds double range");
         w->fitted[(size_t)i * q + o] = v;
         w->residual[(size_t)i * q + o] = y[(size_t)i * q + o] - v;
      }
   /* Publish borrowed views only after every calculation succeeded. No
      coefficients are copied here: the next fit may overwrite these arrays. */
   w->result.n = p->info.n;
   w->result.p = p->info.p;
   w->result.k = p->info.k;
   w->result.q = q;
   w->result.beta = beta;
   w->result.gamma = w->gamma;
   w->result.fitted = w->fitted;
   w->result.residual = w->residual;
   w->result.diagnostic = r ? r->diagnostic : NULL;
   w->result.elastic_diagnostic = e ? e->diagnostic : NULL;
   w->result.banded_diagnostic = b ? b->diagnostic : NULL;
   w->result.support_size =
       p->info.support_refit && !ols ? w->support_size : NULL;
   w->result.support_rank =
       p->info.support_refit && !ols ? w->support_rank : NULL;
   w->valid = 1;
   return 1;
}
/* Public independent fit; never warm-start across response batches. */
int THD_picard_fit_batch(THD_picard_workspace *w, const double *y, int q,
                         const double *parameters, char *err, int nerr)
{
   return picard_fit_model(w, y, q, parameters, 0, 0, err, nerr);
}

/* Only tuning code calls this with warm=1 after the preceding grid point. */
int THD_picard_fit_step(THD_picard_workspace *w, const double *y, int q,
                        const double *parameters, int warm, char *err, int nerr)
{
   return picard_fit_model(w, y, q, parameters, warm, 0, err, nerr);
}

/* Force the unpenalized minimum-norm reference, including for an elastic-net
   plan. HRF selection and repeat tuning use this common reference model. */
int THD_picard_fit_ols(THD_picard_workspace *w, const double *y, int q,
                       char *err, int nerr)
{
   return picard_fit_model(w, y, q, NULL, 0, 1, err, nerr);
}

/* Score the residual in the SAME geometry used to estimate the coefficients.
   A correlated-noise GLS fit minimizes e'R^-1 e, not e'e. With R fixed across
   HRFs, its determinant is common, so weighted SSE provides a comparable
   criterion. We deliberately do not re-estimate R per HRF: restricted
   likelihoods from different mean-model spans cannot simply be ranked.

   Residuals exported to users remain in original data units. Only a private
   scratch copy is whitened here; a caller can still read fitted/residual maps
   after scoring. Keeping this operation outside fit_model avoids an extra
   dense triangular solve in ordinary fits that do not compare HRFs. */
int THD_picard_selection_mse(THD_picard_workspace *w, double *scores,
                             char *err, int nerr)
{
   const double *residual;
   int row, response, n, q;
   if (!w || !w->valid || !scores)
      return picard_error(err, nerr, "HRF scoring needs a successful fit");
   n = w->result.n;
   q = w->result.q;
   residual = w->residual;
   if (w->plan->info.covariance_enabled) {
      picard_whiten(w->plan, residual, q, w->whitened);
      residual = w->whitened;
   }
   for (response = 0; response < q; response++) {
      long double mse = 0;
      for (row = 0; row < n; row++) {
         long double value = residual[(size_t)row*q+response];
         mse += value*value/n;
      }
      if (!isfinite(mse) || mse > DBL_MAX)
         return picard_error(err, nerr, "HRF selection loss exceeds range");
      scores[response] = (double)mse;
   }
   return 1;
}

/* NULL means no successful current fit. The returned view belongs to w and
   remains usable only until its next fit attempt or destruction. */
const THD_picard_result *THD_picard_get_result(const THD_picard_workspace *w)
{
   return w && w->valid ? &w->result : NULL;
}

/* A fold owns only a training design. Held-out targets come from full-design
   OLS: the strict run-local column contract makes that OLS block-separable.
   In particular, we never fit a training fraction using held-out responses. */
typedef struct {
   THD_picard_plan *train;
   int *rows, *columns;
} RepeatFold;
struct THD_picard_repeat_plan {
   const THD_picard_plan *full;
   THD_picard_repeat_info info;
   THD_picard_trial *trials;
   double *fractions, *pairs;
   int *runs;
   RepeatFold *folds;
};
struct THD_picard_repeat_workspace {
   const THD_picard_repeat_plan *plan;
   int maxq, valid;
   THD_picard_workspace *full, **train;
   double *y, *reference, *fractions, *loss, *baseline, *candidate_loss;
   double *fold_loss;
   int *candidate;
   THD_picard_repeat_result result;
};

/* Every destructor accepts partial construction, so failures unwind with
   ordinary returns rather than jumps to a shared cleanup label. */
void THD_picard_repeat_plan_free(THD_picard_repeat_plan *p)
{
   int f;
   if (!p)
      return;
   if (p->folds)
      for (f = 0; f < p->info.nfold; ++f) {
         THD_picard_plan_free(p->folds[f].train);
         free(p->folds[f].rows);
         free(p->folds[f].columns);
      }
   free(p->folds);
   free(p->runs);
   free(p->pairs);
   free(p->trials);
   free(p->fractions);
   free(p);
}

/* Find the sole run with a nonzero entry. -1 means an all-zero column;
   -2 means support spans runs, which would couple their nuisance estimates. */
static int picard_column_run(const THD_picard_plan *p, const double *x,
                             int stride, int column)
{
   int i, run = -1;
   for (i = 0; i < p->info.n; ++i)
      if (x[(size_t)i * stride + column] != 0) {
         if (run >= 0 && run != p->runs[i])
            return -2;
         run = p->runs[i];
      }
   return run;
}

/* Copy a training submatrix, dropping the omitted run's trial and nuisance
   columns. Keep original row identities and a local-to-full beta mapping. */
static int repeat_prepare_fold(THD_picard_repeat_plan *p, int f, char *err,
                               int nerr)
{
   const THD_picard_plan *full = p->full;
   RepeatFold *fold = p->folds + f;
   int i, j, n = 0, np = 0, nk = 0, at, run = p->runs[f];
   int *roles, *row_ids, *run_ids, *nuisance;
   double *matrix;
   for (i = 0; i < full->info.n; ++i)
      n += full->runs[i] != run;
   for (j = 0; j < full->info.p; ++j)
      np += p->trials[j].run != run;
   for (j = 0; j < full->info.k; ++j)
      nk += picard_column_run(full, full->z, full->info.k, j) != run;
   fold->rows = calloc((size_t)n, sizeof(int));
   fold->columns = calloc((size_t)np, sizeof(int));
   roles = calloc((size_t)(np + nk), sizeof(int));
   row_ids = calloc((size_t)n, sizeof(int));
   run_ids = calloc((size_t)n, sizeof(int));
   nuisance = calloc((size_t)(nk ? nk : 1), sizeof(int));
   matrix = picard_array(n, np + nk);
   if (fold->rows && fold->columns && roles && row_ids && run_ids && nuisance &&
       matrix) {
      at = 0;
      for (j = 0; j < full->info.p; ++j)
         if (p->trials[j].run != run) {
            roles[at] = THD_PICARD_PENALIZED;
            fold->columns[at++] = j;
         }
      at = 0;
      for (j = 0; j < full->info.k; ++j)
         if (picard_column_run(full, full->z, full->info.k, j) != run)
            nuisance[at++] = j;
      at = 0;
      for (i = 0; i < full->info.n; ++i)
         if (full->runs[i] != run) {
            fold->rows[at] = i;
            row_ids[at] = full->rows[i];
            run_ids[at] = full->runs[i];
            for (j = 0; j < np; ++j)
               matrix[(size_t)at * (np + nk) + j] =
                   full->x[(size_t)i * full->info.p + fold->columns[j]];
            for (j = 0; j < nk; ++j)
               matrix[(size_t)at * (np + nk) + np + j] =
                   full->z[(size_t)i * full->info.k + nuisance[j]];
            ++at;
         }
      fold->train = THD_picard_plan_like(
          full, matrix, n, np + nk, roles, row_ids, fold->columns, run_ids,
          full->info.nfull, full->info.rank_tol, err, nerr);
   } else {
      picard_error(err, nerr,
                   "cannot allocate repeated-condition training design");
   }
   free(roles);
   free(row_ids);
   free(run_ids);
   free(nuisance);
   free(matrix);
   return fold->train != NULL;
}

/* Validate identities and run isolation before allocating factorizations.
   We deliberately require explicit metadata instead of guessing conditions
   from stimulus labels: an IM label is not a condition identity contract. */
THD_picard_repeat_plan *THD_picard_repeat_plan_create(
    const THD_picard_plan *full, const THD_picard_trial *trials,
    const double *fractions, int count, char *err, int nerr)
{
   THD_picard_repeat_plan *p;
   int i, j, f, nf = 0, scored = 0;
   if (!full || !trials || !fractions || count < 2 || count == INT_MAX) {
      picard_error(err, nerr,
                   "repeat tuning needs a design, trials and >=2 fractions");
      return NULL;
   }
   if (!THD_picard_grid_valid(full, fractions, count, err, nerr))
      return NULL;
   for (i = 0; i < full->info.p; ++i) {
      if (trials[i].column != full->pen[i] || trials[i].trial < 0 ||
          trials[i].condition < 0 || trials[i].run < 0 ||
          picard_column_run(full, full->x, full->info.p, i) != trials[i].run) {
         picard_error(err, nerr,
                      "trial mapping must cover penalized columns with "
                      "nonzero support confined to the declared run");
         return NULL;
      }
      for (j = 0; j < i; ++j)
         if (trials[i].trial == trials[j].trial) {
            picard_error(err, nerr,
                         "trial IDs must be unique (one column per trial)");
            return NULL;
         }
   }
   for (j = 0; j < full->info.k; ++j)
      if (picard_column_run(full, full->z, full->info.k, j) == -2) {
         picard_error(err, nerr,
                      "repeat tuning requires run-local nuisance columns; "
                      "split shared nuisance columns by run");
         return NULL;
      }
   for (i = 0; i < full->info.n; ++i)
      nf += i == 0 || full->runs[i] != full->runs[i - 1];
   if (nf < 2) {
      picard_error(err, nerr,
                   "repeat tuning requires at least two retained runs");
      return NULL;
   }
   p = calloc(1, sizeof(*p));
   if (!p) {
      picard_error(err, nerr, "cannot allocate repeat plan");
      return NULL;
   }
   p->full = full;
   p->info.nfold = nf;
   p->info.ncandidate = count;
   p->trials = calloc((size_t)full->info.p, sizeof(*trials));
   p->fractions = picard_array(count, 1);
   p->pairs = picard_array(nf, 1);
   p->runs = calloc((size_t)nf, sizeof(int));
   p->folds = calloc((size_t)nf, sizeof(*p->folds));
   if (!p->trials || !p->fractions || !p->pairs || !p->runs || !p->folds) {
      THD_picard_repeat_plan_free(p);
      picard_error(err, nerr, "cannot allocate repeat metadata");
      return NULL;
   }
   memcpy(p->trials, trials, (size_t)full->info.p * sizeof(*trials));
   memcpy(p->fractions, fractions, (size_t)count * sizeof(double));
   f = 0;
   for (i = 0; i < full->info.n; ++i)
      if (i == 0 || full->runs[i] != full->runs[i - 1])
         p->runs[f++] = full->runs[i];
   for (i = 0; i < full->info.p; ++i) {
      int matched = 0;
      for (j = 0; j < full->info.p; ++j)
         if (trials[i].run != trials[j].run &&
             trials[i].condition == trials[j].condition) {
            for (f = 0; f < nf; ++f)
               if (p->runs[f] == trials[i].run)
                  p->pairs[f] += 1;
            matched = 1;
         }
      scored += matched;
   }
   for (f = 0; f < nf; ++f) {
      if (p->pairs[f] == 0) {
         THD_picard_repeat_plan_free(p);
         picard_error(
             err, nerr,
             "each run needs a condition repeated in another "
             "run; use a fixed fraction when repeats are insufficient");
         return NULL;
      }
      p->info.npair += p->pairs[f];
      if (!repeat_prepare_fold(p, f, err, nerr)) {
         THD_picard_repeat_plan_free(p);
         return NULL;
      }
   }
   p->info.nscored = scored;
   p->info.trials = p->trials;
   p->info.fractions = p->fractions;
   p->info.fold_pairs = p->pairs;
   p->info.fold_runs = p->runs;
   return p;
}

const THD_picard_repeat_info *
THD_picard_repeat_plan_info(const THD_picard_repeat_plan *p)
{
   return p ? &p->info : NULL;
}

void THD_picard_repeat_workspace_free(THD_picard_repeat_workspace *w)
{
   int f;
   if (!w)
      return;
   if (w->train)
      for (f = 0; f < w->plan->info.nfold; ++f)
         THD_picard_workspace_free(w->train[f]);
   THD_picard_workspace_free(w->full);
   free(w->train);
   free(w->y);
   free(w->reference);
   free(w->fractions);
   free(w->loss);
   free(w->baseline);
   free(w->candidate_loss);
   free(w->fold_loss);
   free(w->candidate);
   free(w);
}

/* Share immutable SVDs across workers, but keep every response, beta, and
   loss buffer private. Memory scales with batch size, never all voxels. */
THD_picard_repeat_workspace *
THD_picard_repeat_workspace_alloc(const THD_picard_repeat_plan *p, int q,
                                  char *err, int nerr)
{
   THD_picard_repeat_workspace *w;
   int f, ok;
   if (!p || q < 1) {
      picard_error(err, nerr, "invalid repeat worker");
      return NULL;
   }
   w = calloc(1, sizeof(*w));
   if (!w) {
      picard_error(err, nerr, "cannot allocate repeat worker");
      return NULL;
   }
   w->plan = p;
   w->maxq = q;
   w->full = THD_picard_workspace_alloc(p->full, q, err, nerr);
   w->train = calloc((size_t)p->info.nfold, sizeof(*w->train));
   w->y = picard_array(p->full->info.n, q);
   w->reference = picard_array(p->full->info.p, q);
   w->fractions = picard_array(q, 1);
   w->loss = picard_array(q, 1);
   w->baseline = picard_array(q, 1);
   w->candidate = calloc((size_t)q, sizeof(int));
   w->candidate_loss = picard_array(p->info.ncandidate, q);
   w->fold_loss =
       picard_array((size_t)p->info.nfold * (p->info.ncandidate + 1), q);
   ok = w->full && w->train && w->y && w->reference && w->fractions &&
        w->loss && w->baseline && w->candidate && w->candidate_loss &&
        w->fold_loss;
   for (f = 0; f < p->info.nfold && ok; ++f) {
      w->train[f] = THD_picard_workspace_alloc(p->folds[f].train, q, err, nerr);
      ok = w->train[f] != NULL;
   }
   if (!ok) {
      THD_picard_repeat_workspace_free(w);
      picard_error(err, nerr, "cannot allocate repeat worker buffers");
      return NULL;
   }
   return w;
}

/* Score ordered training/held-out pairs, not condition averages. Repeated
   presentations within a run therefore contribute their actual multiplicity.
   Long-double accumulation avoids intermediate overflow when a finite MSE
   can still be represented by a double. */
static int repeat_score(THD_picard_repeat_workspace *w, int f, int c,
                        const double *beta, int q, char *err, int nerr)
{
   const THD_picard_repeat_plan *p = w->plan;
   const RepeatFold *fold = p->folds + f;
   int i, j, o;
   for (o = 0; o < q; ++o) {
      long double sum = 0;
      for (i = 0; i < p->full->info.p; ++i)
         if (p->trials[i].run == p->runs[f])
            for (j = 0; j < fold->train->info.p; ++j)
               if (p->trials[fold->columns[j]].condition ==
                   p->trials[i].condition) {
                  long double delta = (long double)beta[(size_t)j * q + o] -
                                      w->reference[(size_t)i * q + o];
                  sum += delta * delta;
               }
      sum /= p->pairs[f];
      if (!isfinite(sum) || sum > DBL_MAX)
         return picard_error(err, nerr,
                             "repeat tuning loss exceeds double range");
      w->fold_loss[((size_t)f * (p->info.ncandidate + 1) + c) * q + o] =
          (double)sum;
   }
   return 1;
}

/* Evaluate OLS targets once, stream training paths, select per response,
   then refit the original full design. Selection never borrows final betas
   from a fold with missing trials or a different coefficient norm. */
int THD_picard_repeat_fit(THD_picard_repeat_workspace *w, const double *y,
                          int q, char *err, int nerr)
{
   const THD_picard_repeat_plan *p;
   const THD_picard_result *r;
   int f, c, i, o, nc;
   if (w)
      w->valid = 0;
   if (!w || !y || q < 1 || q > w->maxq)
      return picard_error(err, nerr,
                          "invalid repeated-condition response batch");
   p = w->plan;
   nc = p->info.ncandidate;
   for (o = 0; o < q; ++o)
      w->fractions[o] = 1;
   if (!THD_picard_fit_ols(w->full, y, q, err, nerr))
      return 0;
   r = THD_picard_get_result(w->full);
   memcpy(w->reference, r->beta, (size_t)r->p * q * sizeof(double));
   for (f = 0; f < p->info.nfold; ++f) {
      const RepeatFold *fold = p->folds + f;
      for (i = 0; i < fold->train->info.n; ++i)
         for (o = 0; o < q; ++o)
            w->y[(size_t)i * q + o] = y[(size_t)fold->rows[i] * q + o];
      /* The extra slot is an OLS baseline even if 1 is absent from the grid. */
      for (c = 0; c <= nc; ++c) {
         for (o = 0; o < q; ++o)
            w->fractions[o] = c == nc ? 1 : p->fractions[c];
         int ok = c == nc ? THD_picard_fit_ols(w->train[f], w->y, q, err, nerr)
                          : THD_picard_fit_step(w->train[f], w->y, q,
                                                w->fractions, c > 0, err, nerr);
         if (!ok)
            return 0;
         r = THD_picard_get_result(w->train[f]);
         if (!repeat_score(w, f, c, r->beta, q, err, nerr))
            return 0;
      }
   }
   for (o = 0; o < q; ++o) {
      int best = 0;
      for (c = 0; c <= nc; ++c) {
         long double sum = 0;
         for (f = 0; f < p->info.nfold; ++f)
            sum += (long double)p->pairs[f] *
                   w->fold_loss[((size_t)f * (nc + 1) + c) * q + o];
         sum /= p->info.npair;
         if (!isfinite(sum) || sum > DBL_MAX)
            return picard_error(err, nerr,
                                "pooled repeat loss exceeds double range");
         if (c == nc)
            w->baseline[o] = (double)sum;
         else {
            double loss = (double)sum;
            w->candidate_loss[(size_t)c * q + o] = loss;
            if (c == 0 || loss < w->loss[o] ||
                (loss == w->loss[o] &&
                 THD_picard_strength(p->full, p->fractions[c]) <
                     THD_picard_strength(p->full, p->fractions[best]))) {
               best = c;
               w->loss[o] = loss;
            }
         }
      }
      w->candidate[o] = best;
      w->fractions[o] = p->fractions[best];
   }
   if (!THD_picard_fit_batch(w->full, y, q, w->fractions, err, nerr))
      return 0;
   w->result.q = q;
   w->result.fit = THD_picard_get_result(w->full);
   w->result.candidate = w->candidate;
   w->result.loss = w->loss;
   w->result.reference_loss = w->baseline;
   w->result.candidate_loss = w->candidate_loss;
   w->result.fold_loss = w->fold_loss;
   w->valid = 1;
   return 1;
}

const THD_picard_repeat_result *
THD_picard_repeat_get_result(const THD_picard_repeat_workspace *w)
{
   return w && w->valid ? &w->result : NULL;
}

/* ========================================================================
   3. TIME-SERIES PREDICTION AND GLOBAL SELECTION
   ======================================================================== */

/* Numeric prediction workflow. Crossval's existing header brings AFNI types
   transitively, but this module neither opens datasets nor uses sample tables.
 */

typedef struct {
   THD_picard_plan *train;
   THD_linproj *test_projection;
   double *x_test, *z_test, *baseline_vector;
} PredictFold;
struct THD_picard_predict_plan {
   const THD_picard_plan *full;
   THD_foldplan *folds;
   THD_picard_predict_info info;
   double *fractions;
   PredictFold *fold;
};
struct THD_picard_predict_workspace {
   const THD_picard_predict_plan *plan;
   THD_picard_workspace *full, **train;
   THD_linproj_workspace **test_projection;
   int maxq, q, valid;
   const double *source;
   double *y_train, *y_test, *fractions, *means;
   double *fold_loss, *fold_baseline, *candidate_loss, *selection;
   double *loss, *baseline, *r2;
   int *candidate, *status;
   THD_picard_predict_result result;
};
static int predict_error(char *err, int nerr, const char *message)
{
   if (err && nerr > 0)
      snprintf(err, (size_t)nerr, "%s", message);
   return 0;
}
static double *predict_array(size_t n, size_t d)
{
   if (!d)
      d = 1;
   return n > SIZE_MAX / sizeof(double) / d ? NULL
                                            : calloc(n * d, sizeof(double));
}
void THD_picard_predict_plan_free(THD_picard_predict_plan *p)
{
   int f;
   if (!p)
      return;
   if (p->fold)
      for (f = 0; f < p->info.nfold; ++f) {
         THD_picard_plan_free(p->fold[f].train);
         THD_linproj_free(p->fold[f].test_projection);
         free(p->fold[f].x_test);
         free(p->fold[f].z_test);
         free(p->fold[f].baseline_vector);
      }
   THD_foldplan_free(p->folds);
   free(p->fold);
   free(p->fractions);
   free(p);
}

/* A training minimum-norm nuisance solution predicts invariantly on test rows
   only if those rows are in its estimable row space. This admits redundant
   nuisance columns while rejecting an unseen held-out-run intercept. */
/* A raw held-out prediction needs nuisance effects that training could
   estimate. Check that each held-out nuisance row lies in the training row
   space. Otherwise different equally good training coefficients would imply
   different test predictions, so reject instead of choosing an arbitrary one.
 */
static int nuisance_estimable(const PredictFold *fold, int nt, char *err,
                              int nerr)
{
   const THD_picard_plan *tr = fold->train;
   THD_numeric_svd *svd;
   double *residual;
   int i, j, r, ok = 1;
   if (!tr->info.k)
      return 1;
   svd = THD_numeric_svd_prepare(tr->weighted_z, tr->info.n, tr->info.k,
                                 tr->info.rank_tol, err, nerr);
   if (!svd)
      return 0;
   residual = predict_array(tr->info.k, 1);
   if (!residual)
      ok = predict_error(err, nerr, "cannot check nuisance estimability");
   for (i = 0; i < nt && ok; ++i) {
      double norm = 0, remaining = 0;
      for (j = 0; j < tr->info.k; ++j) {
         residual[j] = fold->z_test[(size_t)i * tr->info.k + j];
         norm = hypot(norm, residual[j]);
      }
      for (r = 0; r < svd->r; ++r)
         if (svd->s[r] > svd->cutoff) {
            double dot = 0;
            for (j = 0; j < tr->info.k; ++j)
               dot += fold->z_test[(size_t)i * tr->info.k + j] *
                      svd->v[(size_t)j * svd->r + r];
            for (j = 0; j < tr->info.k; ++j)
               residual[j] -= dot * svd->v[(size_t)j * svd->r + r];
         }
      for (j = 0; j < tr->info.k; ++j)
         remaining = hypot(remaining, residual[j]);
      if (remaining > fmax(tr->info.rank_tol, 64 * DBL_EPSILON) * norm)
         ok = predict_error(
             err, nerr,
             "test nuisance is not estimable from training; "
             "use shared nuisance or explicitly request project scoring");
   }
   free(residual);
   THD_numeric_svd_free(svd);
   return ok;
}

/* Preserve the original column roles and scales. Only the row set changes;
   no training-only standardization or implicit intercept is introduced. */
static int prepare_prediction_fold(THD_picard_predict_plan *p, int f, char *err,
                                   int nerr)
{
   const THD_picard_plan *full = p->full;
   PredictFold *fold = p->fold + f;
   int n = p->folds->ntrain[f], nt = p->folds->ntest[f], d = full->info.ncol;
   int *roles = calloc((size_t)d, sizeof(int)),
       *rows = calloc((size_t)n, sizeof(int));
   int *runs = calloc((size_t)n, sizeof(int)), i, j, ok = 1;
   double *matrix = predict_array(n, d);
   THD_linproj_workspace *project = NULL;
   fold->x_test = predict_array(nt, full->info.p);
   fold->z_test = predict_array(nt, full->info.k);
   fold->baseline_vector = predict_array(nt, 1);
   if (!roles || !rows || !runs || !matrix || !fold->x_test || !fold->z_test ||
       !fold->baseline_vector)
      ok = predict_error(err, nerr, "cannot allocate prediction fold");
   if (ok) {
      for (j = 0; j < full->info.p; ++j)
         roles[full->pen[j]] = 1;
      for (i = 0; i < n; ++i) {
         int source = p->folds->train[f][i];
         rows[i] = full->rows[source];
         runs[i] = full->runs[source];
         for (j = 0; j < full->info.p; ++j)
            matrix[(size_t)i * d + full->pen[j]] =
                full->x[(size_t)source * full->info.p + j];
         for (j = 0; j < full->info.k; ++j)
            matrix[(size_t)i * d + full->nuis[j]] =
                full->z[(size_t)source * full->info.k + j];
      }
      for (i = 0; i < nt; ++i) {
         int source = p->folds->test[f][i];
         for (j = 0; j < full->info.p; ++j)
            fold->x_test[(size_t)i * full->info.p + j] =
                full->x[(size_t)source * full->info.p + j];
         for (j = 0; j < full->info.k; ++j)
            fold->z_test[(size_t)i * full->info.k + j] =
                full->z[(size_t)source * full->info.k + j];
         fold->baseline_vector[i] = 1;
      }
      /* Unique-trial columns are not shared predictors and cannot be tuned
         this way. Repeat-condition tuning has a different contract for them. */
      for (j = 0; j < full->info.p && ok; ++j) {
         double norm = 0;
         for (i = 0; i < n; ++i)
            norm = hypot(norm, matrix[(size_t)i * d + full->pen[j]]);
         if (norm == 0)
            ok = predict_error(err, nerr,
                               "each penalized predictor needs "
                               "training support in every fold; trial columns "
                               "need repeat tuning");
      }
      if (ok)
         fold->train = THD_picard_plan_like(full, matrix, n, d, roles, rows,
                                            NULL, runs, full->info.nfull,
                                            full->info.rank_tol, err, nerr);
      ok = ok && fold->train != NULL;
   }
   free(matrix);
   free(roles);
   free(rows);
   free(runs);
   if (!ok)
      return 0;
   if (p->info.nuisance_mode == THD_PICARD_PREDICT_RAW)
      return nuisance_estimable(fold, nt, err, nerr);
   fold->test_projection = THD_linproj_prepare(fold->z_test, nt, full->info.k,
                                               full->info.rank_tol, err, nerr);
   if (!fold->test_projection)
      return 0;
   if (THD_linproj_rank(fold->test_projection) >= nt)
      return predict_error(
          err, nerr,
          "held-out nuisance span removes every test degree of freedom");
   project = THD_linproj_workspace_alloc(fold->test_projection, full->info.p,
                                         err, nerr);
   ok = project &&
        THD_linproj_project(project, fold->x_test, full->info.p, fold->x_test,
                            err, nerr) &&
        THD_linproj_project(project, fold->baseline_vector, 1,
                            fold->baseline_vector, err, nerr);
   THD_linproj_workspace_free(project);
   return ok;
}
THD_picard_predict_plan *
THD_picard_predict_plan_create(const THD_picard_plan *full,
                               const THD_foldplan *folds,
                               const double *fractions, int count, int mode,
                               int one_se, char *err, int nerr)
{
   THD_picard_predict_plan *p;
   int i;
   if (!full || !folds || !fractions || count < 2 ||
       (mode != THD_PICARD_PREDICT_RAW && mode != THD_PICARD_PREDICT_PROJECT) ||
       (one_se != 0 && one_se != 1) || folds->nsample != full->info.n ||
       folds->nfold < 2 || folds->kind == THD_CV_FIXED ||
       !THD_foldplan_validate(folds, err, nerr)) {
      predict_error(
          err, nerr,
          "prediction requires a valid full-coverage CV plan and grid");
      return NULL;
   }
   if (!THD_picard_grid_valid(full, fractions, count, err, nerr))
      return NULL;
   p = calloc(1, sizeof(*p));
   if (!p) {
      predict_error(err, nerr, "cannot allocate prediction plan");
      return NULL;
   }
   p->full = full;
   p->info.nfold = folds->nfold;
   p->info.ncandidate = count;
   p->info.nuisance_mode = mode;
   p->info.one_se = one_se;
   p->fractions = predict_array(count, 1);
   p->fold = calloc((size_t)folds->nfold, sizeof(*p->fold));
   p->folds = THD_foldplan_copy(folds, err, nerr);
   if (!p->fractions || !p->fold || !p->folds) {
      THD_picard_predict_plan_free(p);
      predict_error(err, nerr, "cannot allocate prediction plans");
      return NULL;
   }
   memcpy(p->fractions, fractions, (size_t)count * sizeof(double));
   p->info.fractions = p->fractions;
   for (i = 0; i < folds->nfold; ++i)
      if (!prepare_prediction_fold(p, i, err, nerr)) {
         THD_picard_predict_plan_free(p);
         return NULL;
      }
   return p;
}
const THD_picard_predict_info *
THD_picard_predict_plan_info(const THD_picard_predict_plan *p)
{
   return p ? &p->info : NULL;
}
const THD_foldplan *THD_picard_predict_folds(const THD_picard_predict_plan *p)
{
   return p ? p->folds : NULL;
}

/* This selector also accepts voxel-normalized global fold means. Weighting is
   always by test rows; it never gives a small fold the same mass as a large
   one. */
int THD_picard_predict_select(const THD_picard_predict_plan *p,
                              const double *losses, int *winner,
                              double *best_se, char *err, int nerr)
{
   int c, f, best = 0, nc;
   long double best_mean = 0, sumw2 = 0, se;
   if (!p || !losses || !winner || !best_se)
      return predict_error(err, nerr, "invalid prediction selection inputs");
   nc = p->info.ncandidate;
   for (f = 0; f < p->info.nfold; ++f) {
      long double weight = (long double)p->folds->ntest[f] / p->full->info.n;
      sumw2 += weight * weight;
   }
   for (c = 0; c < nc; ++c) {
      long double mean = 0;
      for (f = 0; f < p->info.nfold; ++f) {
         double value = losses[(size_t)f * nc + c];
         if (!isfinite(value) || value < 0)
            return predict_error(err, nerr, "invalid fold loss");
         mean += (long double)p->folds->ntest[f] / p->full->info.n * value;
      }
      if (!c || mean < best_mean ||
          (mean == best_mean &&
           THD_picard_strength(p->full, p->fractions[c]) <
               THD_picard_strength(p->full, p->fractions[best]))) {
         best = c;
         best_mean = mean;
      }
   }
   se = 0;
   for (f = 0; f < p->info.nfold; ++f) {
      long double delta = losses[(size_t)f * nc + best] - best_mean;
      se += (long double)p->folds->ntest[f] / p->full->info.n * delta * delta;
   }
   se = sqrtl(se * sumw2 / (1 - sumw2));
   if (!isfinite(se) || se > DBL_MAX)
      return predict_error(err, nerr, "selection SE exceeds double range");
   if (p->info.one_se)
      for (c = 0; c < nc; ++c) {
         long double mean = 0;
         for (f = 0; f < p->info.nfold; ++f)
            mean += (long double)p->folds->ntest[f] / p->full->info.n *
                    losses[(size_t)f * nc + c];
         if (mean <= best_mean + se &&
             THD_picard_strength(p->full, p->fractions[c]) >
                 THD_picard_strength(p->full, p->fractions[best]))
            best = c;
      }
   *winner = best;
   *best_se = (double)se;
   return 1;
}

void THD_picard_predict_workspace_free(THD_picard_predict_workspace *w)
{
   int f;
   if (!w)
      return;
   for (f = 0; f < w->plan->info.nfold; ++f) {
      if (w->train)
         THD_picard_workspace_free(w->train[f]);
      if (w->test_projection)
         THD_linproj_workspace_free(w->test_projection[f]);
   }
   THD_picard_workspace_free(w->full);
   free(w->train);
   free(w->test_projection);
   free(w->y_train);
   free(w->y_test);
   free(w->fractions);
   free(w->means);
   free(w->fold_loss);
   free(w->fold_baseline);
   free(w->candidate_loss);
   free(w->selection);
   free(w->loss);
   free(w->baseline);
   free(w->r2);
   free(w->candidate);
   free(w->status);
   free(w);
}
THD_picard_predict_workspace *
THD_picard_predict_workspace_alloc(const THD_picard_predict_plan *p, int q,
                                   char *err, int nerr)
{
   THD_picard_predict_workspace *w;
   int f, ok;
   if (!p || q < 1) {
      predict_error(err, nerr, "invalid prediction worker size");
      return NULL;
   }
   w = calloc(1, sizeof(*w));
   if (!w) {
      predict_error(err, nerr, "cannot allocate prediction worker");
      return NULL;
   }
   w->plan = p;
   w->maxq = q;
   w->full = THD_picard_workspace_alloc(p->full, q, err, nerr);
   w->train = calloc((size_t)p->info.nfold, sizeof(*w->train));
   w->test_projection =
       calloc((size_t)p->info.nfold, sizeof(*w->test_projection));
   w->y_train = predict_array(p->full->info.n, q);
   w->y_test = predict_array(p->full->info.n, q);
   w->fractions = predict_array(q, 1);
   w->means = predict_array(q, 1);
   w->fold_loss = predict_array((size_t)p->info.nfold * p->info.ncandidate, q);
   w->fold_baseline = predict_array(p->info.nfold, q);
   w->candidate_loss = predict_array(p->info.ncandidate, q);
   w->selection = predict_array(p->info.nfold, p->info.ncandidate);
   w->loss = predict_array(q, 1);
   w->baseline = predict_array(q, 1);
   w->r2 = predict_array(q, 1);
   w->candidate = calloc((size_t)q, sizeof(int));
   w->status = calloc((size_t)q, sizeof(int));
   ok = w->full && w->train && w->test_projection && w->y_train && w->y_test &&
        w->fractions && w->means && w->fold_loss && w->fold_baseline &&
        w->candidate_loss && w->selection && w->loss && w->baseline && w->r2 &&
        w->candidate && w->status;
   for (f = 0; f < p->info.nfold && ok; ++f) {
      w->train[f] = THD_picard_workspace_alloc(p->fold[f].train, q, err, nerr);
      ok = w->train[f] != NULL;
      if (ok && p->fold[f].test_projection) {
         w->test_projection[f] = THD_linproj_workspace_alloc(
             p->fold[f].test_projection, q, err, nerr);
         ok = w->test_projection[f] != NULL;
      }
   }
   if (!ok) {
      THD_picard_predict_workspace_free(w);
      predict_error(err, nerr, "cannot allocate prediction scratch");
      return NULL;
   }
   return w;
}

/* Runtime preparation gathers response rows only. SVDs were prepared once in
   plan creation, before seeing any response data, and shared by all workers. */
static int prediction_prepare(void *context, int f, const int *train, int n,
                              const int *test, int nt, char *err, int nerr)
{
   THD_picard_predict_workspace *w = context;
   int i, o, q = w->q;
   for (o = 0; o < q; ++o) {
      long double mean = 0;
      for (i = 0; i < n; ++i) {
         double value = w->source[(size_t)train[i] * q + o];
         w->y_train[(size_t)i * q + o] = value;
         mean += (long double)value / n;
      }
      w->means[o] = (double)mean;
   }
   for (i = 0; i < nt; ++i)
      for (o = 0; o < q; ++o)
         w->y_test[(size_t)i * q + o] = w->source[(size_t)test[i] * q + o];
   /* Conditional scoring removes the supplied test nuisance span from Y.
      Raw scoring leaves Y alone and uses nuisance coefficients learned on
      training rows. This branch is an explicit scoring policy, not a fit
      of new penalized coefficients to held-out responses. */
   if (w->test_projection[f] &&
       !THD_linproj_project(w->test_projection[f], w->y_test, q, w->y_test, err,
                            nerr))
      return 0;
   for (o = 0; o < q; ++o) {
      long double sum = 0;
      for (i = 0; i < nt; ++i) {
         long double delta =
             w->y_test[(size_t)i * q + o] -
             (long double)w->plan->fold[f].baseline_vector[i] * w->means[o];
         sum += delta * delta / nt;
      }
      if (!isfinite(sum) || sum > DBL_MAX)
         return predict_error(err, nerr,
                              "prediction baseline exceeds double range");
      w->fold_baseline[(size_t)f * q + o] = (double)sum;
   }
   return 1;
}
/* Fit the parameter grid using training Y, then predict held-out rows.
   Each fold starts cold; later elastic-net grid points may reuse the previous
   coefficients in that same fold. Held-out Y enters only the error score. */
static int prediction_path(void *context, int f, char *err, int nerr)
{
   THD_picard_predict_workspace *w = context;
   const THD_picard_predict_plan *p = w->plan;
   const PredictFold *fold = p->fold + f;
   int c, i, j, o, q = w->q, nt = p->folds->ntest[f];
   for (c = 0; c < p->info.ncandidate; ++c) {
      const THD_picard_result *r;
      for (o = 0; o < q; ++o)
         w->fractions[o] = p->fractions[c];
      if (!THD_picard_fit_step(w->train[f], w->y_train, q, w->fractions, c > 0,
                               err, nerr))
         return 0;
      r = THD_picard_get_result(w->train[f]);
      for (o = 0; o < q; ++o) {
         long double sum = 0;
         for (i = 0; i < nt; ++i) {
            long double prediction = 0, delta;
            for (j = 0; j < r->p; ++j)
               prediction += (long double)fold->x_test[(size_t)i * r->p + j] *
                             r->beta[(size_t)j * q + o];
            if (p->info.nuisance_mode == THD_PICARD_PREDICT_RAW)
               for (j = 0; j < r->k; ++j)
                  prediction +=
                      (long double)fold->z_test[(size_t)i * r->k + j] *
                      r->gamma[(size_t)j * q + o];
            delta = w->y_test[(size_t)i * q + o] - prediction;
            sum += delta * delta / nt;
         }
         if (!isfinite(sum) || sum > DBL_MAX)
            return predict_error(err, nerr,
                                 "prediction loss exceeds double range");
         w->fold_loss[((size_t)f * p->info.ncandidate + c) * q + o] =
             (double)sum;
      }
   }
   return 1;
}
/* Pool held-out errors for one response at a time. Weight each fold by its
   number of test rows so each held-out sample contributes equally. Preserve
   fold-level losses as well: the one-SE rule needs their variability. */
static int prediction_finish(void *context, char *err, int nerr)
{
   THD_picard_predict_workspace *w = context;
   const THD_picard_predict_plan *p = w->plan;
   int o, f, c, q = w->q, nc = p->info.ncandidate;
   for (o = 0; o < q; ++o) {
      long double baseline = 0;
      double se;
      for (c = 0; c < nc; ++c) {
         long double mean = 0;
         for (f = 0; f < p->info.nfold; ++f) {
            double value = w->fold_loss[((size_t)f * nc + c) * q + o];
            w->selection[(size_t)f * nc + c] = value;
            mean += (long double)p->folds->ntest[f] / p->full->info.n * value;
         }
         w->candidate_loss[(size_t)c * q + o] = (double)mean;
      }
      for (f = 0; f < p->info.nfold; ++f)
         baseline += (long double)p->folds->ntest[f] / p->full->info.n *
                     w->fold_baseline[(size_t)f * q + o];
      w->baseline[o] = (double)baseline;
      if (!THD_picard_predict_select(p, w->selection, w->candidate + o, &se,
                                     err, nerr))
         return 0;
   }
   return 1;
}
/* Run prepare/path/finish callbacks over the explicit folds, then refit the
   selected parameter on all retained rows. score_only skips that last fit
   for a global-selection prepass; forced>=0 applies the previously chosen
   global index while still computing this batch's prediction diagnostics. */
int THD_picard_predict_fit(THD_picard_predict_workspace *w, const double *y,
                           int q, int score_only, int forced, char *err,
                           int nerr)
{
   const THD_cv_path_ops ops = {prediction_prepare, prediction_path,
                                prediction_finish};
   int i, o;
   if (w) {
      w->valid = 0;
      w->result.fit = NULL;
   }
   if (!w || !y || q < 1 || q > w->maxq || forced < -1 ||
       forced >= w->plan->info.ncandidate)
      return predict_error(err, nerr,
                           "invalid prediction batch/forced candidate");
   for (i = 0; i < w->plan->full->info.n; ++i)
      for (o = 0; o < q; ++o)
         if (!isfinite(y[(size_t)i * q + o]))
            return predict_error(err, nerr, "nonfinite prediction response");
   w->source = y;
   w->q = q;
   if (!THD_cv_path_search(w->plan->folds, &ops, w, err, nerr))
      return 0;
   for (o = 0; o < q; ++o) {
      if (forced >= 0)
         w->candidate[o] = forced;
      w->fractions[o] = w->plan->fractions[w->candidate[o]];
      w->loss[o] = w->candidate_loss[(size_t)w->candidate[o] * q + o];
      w->status[o] = w->baseline[o] == 0;
      w->r2[o] = w->status[o] ? NAN : 1 - w->loss[o] / w->baseline[o];
      if (!w->status[o] && !isfinite(w->r2[o]))
         return predict_error(err, nerr, "prediction R2 exceeds double range");
   }
   /* The exported coefficients are from a new full-data fit, not an average
      of fold coefficients. CV scores still describe held-out fold fits. */
   if (!score_only) {
      if (!THD_picard_fit_batch(w->full, y, q, w->fractions, err, nerr))
         return 0;
      w->result.fit = THD_picard_get_result(w->full);
   }
   w->result.q = q;
   w->result.candidate = w->candidate;
   w->result.status = w->status;
   w->result.loss = w->loss;
   w->result.baseline = w->baseline;
   w->result.r2 = w->r2;
   w->result.candidate_loss = w->candidate_loss;
   w->result.fold_loss = w->fold_loss;
   w->result.fold_baseline = w->fold_baseline;
   w->valid = 1;
   return 1;
}
const THD_picard_predict_result *
THD_picard_predict_get_result(const THD_picard_predict_workspace *w)
{
   return w && w->valid ? &w->result : NULL;
}

struct THD_picard_predict_global {
   const THD_picard_predict_plan *plan;
   long double *sum;
   double *mean;
   int nvalid, finished, failed;
   THD_picard_predict_global_result result;
};
void THD_picard_predict_global_free(THD_picard_predict_global *g)
{
   if (!g)
      return;
   free(g->sum);
   free(g->mean);
   free(g);
}
THD_picard_predict_global *
THD_picard_predict_global_alloc(const THD_picard_predict_plan *p, char *err,
                                int nerr)
{
   THD_picard_predict_global *g;
   size_t cells;
   if (!p || (size_t)p->info.nfold >
                 SIZE_MAX / sizeof(long double) / p->info.ncandidate) {
      predict_error(err, nerr, "invalid global prediction dimensions");
      return NULL;
   }
   cells = (size_t)p->info.nfold * p->info.ncandidate;
   g = calloc(1, sizeof(*g));
   if (!g) {
      predict_error(err, nerr, "cannot allocate global prediction");
      return NULL;
   }
   g->plan = p;
   g->sum = calloc(cells, sizeof(long double));
   g->mean = calloc(cells, sizeof(double));
   if (!g->sum || !g->mean) {
      THD_picard_predict_global_free(g);
      predict_error(err, nerr, "cannot allocate global score buffers");
      return NULL;
   }
   return g;
}

/* Retain no voxel-sized score table. The explicit response loop fixes the
   summation order even when callers change the batch capacity. */
int THD_picard_predict_global_add(THD_picard_predict_global *g,
                                  const THD_picard_predict_workspace *w,
                                  char *err, int nerr)
{
   int o;
   size_t cell, cells;
   if (!g || !w || !w->valid || g->plan != w->plan || g->finished ||
       g->failed) {
      if (g)
         g->failed = 1;
      return predict_error(err, nerr, "invalid global prediction batch/state");
   }
   /* Divide by each voxel's baseline error before pooling. Otherwise a
      voxel with large signal units could dominate global selection solely
      because its squared errors are larger. A zero baseline is undefined
      for this normalization and contributes no vote. */
   cells = (size_t)g->plan->info.nfold * g->plan->info.ncandidate;
   for (o = 0; o < w->q; ++o)
      if (w->baseline[o] > 0) {
         if (g->nvalid == INT_MAX) {
            g->failed = 1;
            return predict_error(err, nerr, "too many global responses");
         }
         ++g->nvalid;
         for (cell = 0; cell < cells; ++cell)
            g->sum[cell] +=
                (long double)w->fold_loss[cell * w->q + o] / w->baseline[o];
      }
   return 1;
}
const THD_picard_predict_global_result *
THD_picard_predict_global_finish(THD_picard_predict_global *g, char *err,
                                 int nerr)
{
   size_t cell, cells;
   int c;
   if (!g || g->finished || g->failed) {
      predict_error(err, nerr, "global prediction already finished or invalid");
      return NULL;
   }
   g->finished = 1;
   cells = (size_t)g->plan->info.nfold * g->plan->info.ncandidate;
   if (g->nvalid) {
      for (cell = 0; cell < cells; ++cell) {
         long double value = g->sum[cell] / g->nvalid;
         if (!isfinite(value) || value > DBL_MAX) {
            predict_error(err, nerr,
                          "normalized global loss exceeds double range");
            return NULL;
         }
         g->mean[cell] = (double)value;
      }
      if (!THD_picard_predict_select(g->plan, g->mean, &g->result.candidate,
                                     &g->result.best_se, err, nerr))
         return NULL;
   } else {
      for (c = 1; c < g->plan->info.ncandidate; ++c)
         if (THD_picard_strength(g->plan->full, g->plan->fractions[c]) <
             THD_picard_strength(g->plan->full,
                                 g->plan->fractions[g->result.candidate]))
            g->result.candidate = c;
   }
   g->result.nvalid = g->nvalid;
   g->result.fold_loss = g->mean;
   return &g->result;
}

/* ========================================================================
   4. CANDIDATE DESIGNS AND LOCKED HRF CHOICES
   ======================================================================== */

struct THD_picard_candidates {
   int count, *ids;
   const THD_picard_plan **design;
   THD_picard_candidates_config config;
   THD_picard_repeat_plan **repeat;
   THD_picard_predict_plan **predict;
};
struct THD_picard_candidates_workspace {
   const THD_picard_candidates *plan;
   int maxq, valid;
   THD_picard_workspace **fixed;
   THD_picard_repeat_workspace **repeat;
   THD_picard_predict_workspace **predict;
   double *parameter, *selected_y, *beta, *gamma, *fitted, *residual;
   int *selected_response;
   THD_fraction_diagnostic *fraction_diagnostic;
   THD_elastic_diagnostic *elastic_diagnostic;
   THD_banded_diagnostic *banded_diagnostic;
   int *support_size, *support_rank;
   int *index, *id, *rank, *parameter_index, *status;
   double *selection_loss, *design_loss, *loss, *baseline, *r2, *parameter_loss;
   THD_picard_result fit;
   THD_picard_candidates_result result;
};

static int candidate_error(char *err, int nerr, const char *message)
{
   if (err && nerr > 0)
      snprintf(err, (size_t)nerr, "%s", message);
   return 0;
}
/* Allocate at least one element for legitimate zero-column/grid buffers. */
static double *candidate_array(size_t rows, size_t cols)
{
   if (!rows)
      rows = 1;
   if (!cols)
      cols = 1;
   return rows > SIZE_MAX / sizeof(double) / cols
              ? NULL
              : calloc(rows * cols, sizeof(double));
}

/* Compare everything that gives numeric columns and observations their role.
   Semantic string/event identities belong to the file adapter, not a solver. */
/* Compare scientific identity requirements before allocating a collection.
   Task values may differ between candidates, but response sampling, nuisance
   values, coefficient roles and estimator settings must have the same meaning
   so the winning coefficients can share one output layout. */
static int candidate_compatible(const THD_picard_plan *a,
                                const THD_picard_plan *b)
{
   const THD_picard_design_info *x = &a->info, *y = &b->info;
   return x->n == y->n && x->ncol == y->ncol && x->p == y->p && x->k == y->k &&
          x->nfull == y->nfull && x->rank_tol == y->rank_tol &&
          x->elastic == y->elastic && x->banded == y->banded &&
          x->ar1 == y->ar1 && x->support_refit == y->support_refit &&
          x->covariance_enabled == y->covariance_enabled &&
          x->covariance.model == y->covariance.model &&
          x->covariance.reset_gaps == y->covariance.reset_gaps &&
          x->covariance.a == y->covariance.a &&
          x->covariance.b == y->covariance.b &&
          (!x->banded ||
           (x->ngroup == y->ngroup && x->ncandidate == y->ncandidate &&
            !memcmp(a->groups, b->groups, sizeof(int) * (size_t)x->p) &&
            !memcmp(a->penalties, b->penalties,
                    sizeof(double) * (size_t)x->ngroup * x->ncandidate))) &&
          (!x->elastic ||
           (a->config.rho == b->config.rho &&
            a->config.tolerance == b->config.tolerance &&
            a->config.max_iterations == b->config.max_iterations &&
            a->config.standardize == b->config.standardize)) &&
          !memcmp(a->pen, b->pen, (size_t)x->p * sizeof(int)) &&
          !memcmp(a->nuis, b->nuis, (size_t)x->k * sizeof(int)) &&
          !memcmp(a->rows, b->rows, (size_t)x->n * sizeof(int)) &&
          !memcmp(a->runs, b->runs, (size_t)x->n * sizeof(int)) &&
          !memcmp(a->z, b->z, (size_t)x->n * x->k * sizeof(double));
}

/* Fixed designs are borrowed; only the collection's own tuning plans die. */
void THD_picard_candidates_free(THD_picard_candidates *p)
{
   int c;
   if (!p)
      return;
   for (c = 0; c < p->count; c++) {
      if (p->repeat)
         THD_picard_repeat_plan_free(p->repeat[c]);
      if (p->predict)
         THD_picard_predict_plan_free(p->predict[c]);
   }
   free(p->repeat);
   free(p->predict);
   free(p->design);
   free(p->ids);
   free(p);
}

THD_picard_candidates *THD_picard_candidates_create(
    const THD_picard_plan *const *designs, const int *ids, int count,
    const THD_picard_candidates_config *cfg, char *err, int nerr)
{
   THD_picard_candidates *p;
   int c, j, ok = 1;
   if (!designs || !ids || count < 1 || !cfg || !designs[0] ||
       (cfg->policy != THD_PICARD_SELECT_HRF &&
        cfg->policy != THD_PICARD_SELECT_PREDICT) ||
       cfg->nparameter < 0 ||
       (cfg->nparameter && (!cfg->parameters || cfg->nparameter < 2)) ||
       (cfg->policy == THD_PICARD_SELECT_HRF && cfg->folds) ||
       (cfg->policy == THD_PICARD_SELECT_HRF && cfg->nparameter &&
        !cfg->trials) ||
       (cfg->policy == THD_PICARD_SELECT_PREDICT &&
        (!cfg->folds || !cfg->nparameter || cfg->trials)) ||
       (!cfg->nparameter &&
        (!isfinite(cfg->fixed_parameter) ||
         (designs[0]->banded ? 0
          : designs[0]->elastic
              ? cfg->fixed_parameter <= 0
              : cfg->fixed_parameter < 0 || cfg->fixed_parameter > 1)))) {
      candidate_error(err, nerr, "invalid candidate collection configuration");
      return NULL;
   }
   if (!cfg->nparameter &&
       !THD_picard_grid_valid(designs[0], &cfg->fixed_parameter, 1, err, nerr))
      return NULL;
   for (c = 0; c < count; c++) {
      if (!designs[c] || ids[c] < 0 || ids[c] > 16777215 ||
          !candidate_compatible(designs[0], designs[c])) {
         candidate_error(
             err, nerr,
             "candidate numeric identities/nuisance/settings differ");
         return NULL;
      }
      if (cfg->policy == THD_PICARD_SELECT_HRF &&
          designs[c]->info.xrank != designs[0]->info.xrank) {
         candidate_error(err, nerr,
                         "HRF candidates must have equal projected rank");
         return NULL;
      }
      for (j = 0; j < c; j++)
         if (ids[c] == ids[j]) {
            candidate_error(err, nerr, "duplicate stable candidate ID");
            return NULL;
         }
   }
   p = calloc(1, sizeof(*p));
   if (!p) {
      candidate_error(err, nerr, "cannot allocate candidate collection");
      return NULL;
   }
   p->count = count;
   p->config = *cfg;
   /* Pointer fields in config are construction inputs only, not retained. */
   p->config.parameters = NULL;
   p->config.trials = NULL;
   p->config.folds = NULL;
   p->ids = calloc((size_t)count, sizeof(int));
   p->design = calloc((size_t)count, sizeof(*p->design));
   p->repeat = calloc((size_t)count, sizeof(*p->repeat));
   p->predict = calloc((size_t)count, sizeof(*p->predict));
   ok = p->ids && p->design && p->repeat && p->predict;
   if (!ok)
      candidate_error(err, nerr, "cannot allocate candidate plan arrays");
   for (c = 0; c < count && ok; c++) {
      p->ids[c] = ids[c];
      p->design[c] = designs[c];
      if (cfg->policy == THD_PICARD_SELECT_PREDICT) {
         p->predict[c] = THD_picard_predict_plan_create(
             designs[c], cfg->folds, cfg->parameters, cfg->nparameter,
             cfg->nuisance_mode, cfg->one_se, err, nerr);
         ok = p->predict[c] != NULL;
      } else if (cfg->nparameter) {
         p->repeat[c] = THD_picard_repeat_plan_create(
             designs[c], cfg->trials, cfg->parameters, cfg->nparameter, err,
             nerr);
         ok = p->repeat[c] != NULL;
      }
   }
   if (!ok) {
      THD_picard_candidates_free(p);
      return NULL;
   }
   return p;
}

void THD_picard_candidates_workspace_free(THD_picard_candidates_workspace *w)
{
   int c;
   if (!w)
      return;
   for (c = 0; c < w->plan->count; c++) {
      if (w->fixed)
         THD_picard_workspace_free(w->fixed[c]);
      if (w->repeat)
         THD_picard_repeat_workspace_free(w->repeat[c]);
      if (w->predict)
         THD_picard_predict_workspace_free(w->predict[c]);
   }
   free(w->fixed);
   free(w->repeat);
   free(w->predict);
   free(w->selected_y);
   free(w->selected_response);
   free(w->parameter);
   free(w->beta);
   free(w->gamma);
   free(w->fitted);
   free(w->residual);
   free(w->fraction_diagnostic);
   free(w->elastic_diagnostic);
   free(w->banded_diagnostic);
   free(w->support_size);
   free(w->support_rank);
   free(w->index);
   free(w->id);
   free(w->rank);
   free(w->parameter_index);
   free(w->status);
   free(w->selection_loss);
   free(w->design_loss);
   free(w->loss);
   free(w->baseline);
   free(w->r2);
   free(w->parameter_loss);
   free(w);
}

/* Plans share immutable factors. Each candidate retains batch-bounded solver
   scratch; export buffers retain one winning fit per response plus scalar
   design/path score grids for the requested batch. */
THD_picard_candidates_workspace *
THD_picard_candidates_workspace_alloc(const THD_picard_candidates *p, int q,
                                      char *err, int nerr)
{
   THD_picard_candidates_workspace *w;
   const THD_picard_design_info *d;
   int c, ok;
   if (!p || q < 1) {
      candidate_error(err, nerr, "invalid candidate worker capacity");
      return NULL;
   }
   w = calloc(1, sizeof(*w));
   if (!w) {
      candidate_error(err, nerr, "cannot allocate candidate worker");
      return NULL;
   }
   w->plan = p;
   w->maxq = q;
   d = &p->design[0]->info;
   w->fixed = calloc((size_t)p->count, sizeof(*w->fixed));
   w->repeat = calloc((size_t)p->count, sizeof(*w->repeat));
   w->predict = calloc((size_t)p->count, sizeof(*w->predict));
   w->selected_y = candidate_array(d->n, q);
   w->selected_response = calloc((size_t)q, sizeof(int));
   w->parameter = candidate_array(q, 1);
   w->beta = candidate_array(d->p, q);
   w->gamma = candidate_array(d->k, q);
   w->fitted = candidate_array(d->n, q);
   w->residual = candidate_array(d->n, q);
   w->fraction_diagnostic = calloc((size_t)q, sizeof(*w->fraction_diagnostic));
   w->elastic_diagnostic = calloc((size_t)q, sizeof(*w->elastic_diagnostic));
   w->banded_diagnostic = calloc((size_t)q, sizeof(*w->banded_diagnostic));
   w->support_size = calloc((size_t)q, sizeof(int));
   w->support_rank = calloc((size_t)q, sizeof(int));
   w->index = calloc((size_t)q, sizeof(int));
   w->id = calloc((size_t)q, sizeof(int));
   w->rank = calloc((size_t)q, sizeof(int));
   w->parameter_index = calloc((size_t)q, sizeof(int));
   w->status = calloc((size_t)q, sizeof(int));
   w->selection_loss = candidate_array(q, 1);
   w->design_loss = candidate_array(p->count, q);
   w->loss = candidate_array(q, 1);
   w->baseline = candidate_array(q, 1);
   w->r2 = candidate_array(q, 1);
   w->parameter_loss = candidate_array(p->config.nparameter, q);
   ok = w->banded_diagnostic && w->support_size && w->support_rank &&
        w->selected_y && w->selected_response && w->fixed && w->repeat &&
        w->predict && w->parameter && w->beta && w->gamma && w->fitted &&
        w->residual && w->fraction_diagnostic && w->elastic_diagnostic &&
        w->index && w->id && w->rank && w->parameter_index && w->status &&
        w->selection_loss && w->design_loss && w->loss && w->baseline &&
        w->r2 && w->parameter_loss;
   if (!ok)
      candidate_error(err, nerr, "cannot allocate candidate result buffers");
   for (c = 0; c < p->count && ok; c++) {
      if (p->predict[c]) {
         w->predict[c] =
             THD_picard_predict_workspace_alloc(p->predict[c], q, err, nerr);
         ok = w->predict[c] != NULL;
      } else {
         w->fixed[c] = THD_picard_workspace_alloc(p->design[c], q, err, nerr);
         ok = w->fixed[c] != NULL;
         if (ok && p->repeat[c]) {
            w->repeat[c] =
                THD_picard_repeat_workspace_alloc(p->repeat[c], q, err, nerr);
            ok = w->repeat[c] != NULL;
         }
      }
   }
   if (!ok) {
      THD_picard_candidates_workspace_free(w);
      return NULL;
   }
   return w;
}

/* Copy one response before reusing the candidate's borrowed solver result.
   source indexes the candidate's possibly compacted batch; o indexes the
   original batch. Their widths (r->q and q) can differ, so each side needs
   its own row stride. This is essential when different voxels choose HRFs. */
static void candidate_retain_fit(THD_picard_candidates_workspace *w,
                                 const THD_picard_result *r, int source, int o,
                                 int q)
{
   int j;
   for (j = 0; j < r->p; j++)
      w->beta[(size_t)j * q + o] = r->beta[(size_t)j * r->q + source];
   for (j = 0; j < r->k; j++)
      w->gamma[(size_t)j * q + o] = r->gamma[(size_t)j * r->q + source];
   for (j = 0; j < r->n; j++) {
      w->fitted[(size_t)j * q + o] = r->fitted[(size_t)j * r->q + source];
      w->residual[(size_t)j * q + o] = r->residual[(size_t)j * r->q + source];
   }
   if (r->diagnostic)
      w->fraction_diagnostic[o] = r->diagnostic[source];
   if (r->elastic_diagnostic)
      w->elastic_diagnostic[o] = r->elastic_diagnostic[source];
   if (r->banded_diagnostic)
      w->banded_diagnostic[o] = r->banded_diagnostic[source];
   if (r->support_size) {
      w->support_size[o] = r->support_size[source];
      w->support_rank[o] = r->support_rank[source];
   }
}

/* Preserve the winning design's fit AND its penalty-selection diagnostics.
   Copying the coefficient maps alone would leave scores from another design. */
static void candidate_retain_prediction(THD_picard_candidates_workspace *w,
                                        const THD_picard_predict_result *r,
                                        int o, int q)
{
   int g;
   candidate_retain_fit(w, r->fit, o, o, q);
   w->parameter_index[o] = r->candidate[o];
   w->status[o] = r->status[o];
   w->loss[o] = r->loss[o];
   w->baseline[o] = r->baseline[o];
   w->r2[o] = r->r2[o];
   for (g = 0; g < w->plan->config.nparameter; g++)
      w->parameter_loss[(size_t)g * q + o] =
          r->candidate_loss[(size_t)g * q + o];
}

/* Select all HRFs before any penalized refit: regularization must not alter
   the HRF policy. Prediction instead scores each design using fold-local fits.
 */
static int candidate_fit(THD_picard_candidates_workspace *w, const double *y,
                         int q, const int *forced, const double *scores,
                         char *err, int nerr)
{
   const THD_picard_candidates *p;
   const THD_picard_design_info *d;
   int c, o, i, g;
   if (w)
      w->valid = 0;
   if (!w || !y || q < 1 || q > w->maxq)
      return candidate_error(err, nerr, "invalid candidate response batch");
   p = w->plan;
   d = &p->design[0]->info;
   if (forced) {
      if (p->config.policy != THD_PICARD_SELECT_HRF || !scores)
         return candidate_error(err, nerr,
                                "locked designs require HRF policy and scores");
      for (o = 0; o < q; o++) {
         c = forced[o];
         if (c < 0 || c >= p->count || !isfinite(scores[o]) || scores[o] < 0)
            return candidate_error(err, nerr, "invalid locked HRF selection");
         w->index[o] = c;
         w->id[o] = p->ids[c];
         w->rank[o] = p->design[c]->info.xrank;
         w->selection_loss[o] = scores[o];
      }
   }
   /* Pass 1: compare all designs. HRF policy uses full-data OLS error;
      prediction policy uses each design's selected cross-validation error.
      Both are selection scores: these same data helped choose the design,
      so neither is an independent assessment of the final chosen workflow. */
   for (c = 0; !forced && c < p->count; c++) {
      const THD_picard_predict_result *prediction = NULL;
      if (w->predict[c]) {
         if (!THD_picard_predict_fit(w->predict[c], y, q, 0, -1, err, nerr))
            return 0;
         prediction = THD_picard_predict_get_result(w->predict[c]);
      } else {
         if (!THD_picard_fit_ols(w->fixed[c], y, q, err, nerr))
            return 0;
         if (!THD_picard_selection_mse(w->fixed[c], w->design_loss+(size_t)c*q,
                                        err, nerr))
            return 0;
      }
      for (o = 0; o < q; o++) {
         long double mse = 0;
         if (prediction)
            mse = prediction->loss[o];
         else
            mse = w->design_loss[(size_t)c*q+o];
         if (!isfinite(mse) || mse > DBL_MAX)
            return candidate_error(err, nerr,
                                   "candidate selection loss exceeds range");
         w->design_loss[(size_t)c * q + o] = (double)mse;
         /* Initialize from the first candidate; exact ties use stable IDs,
            not manifest order. Keep only the best fit per response, while
            retaining scalar losses for optional per-design diagnostics. */
         if (c == 0 || (double)mse < w->selection_loss[o] ||
             ((double)mse == w->selection_loss[o] && p->ids[c] < w->id[o])) {
            w->index[o] = c;
            w->id[o] = p->ids[c];
            w->rank[o] = p->design[c]->info.xrank;
            w->selection_loss[o] = (double)mse;
            if (prediction)
               candidate_retain_prediction(w, prediction, o, q);
         }
      }
   }
   /* Pass 2 (HRF only): the design choices are now fixed. Group the winning
      responses by design before applying the requested penalty/tuning. The
      prediction policy already refitted each design during pass 1. */
   if (p->config.policy == THD_PICARD_SELECT_HRF)
      for (c = 0; c < p->count; c++) {
         const THD_picard_result *fit;
         const THD_picard_repeat_result *repeat = NULL;
         int selected = 0, source;
         for (o = 0; o < q; o++)
            if (w->index[o] == c)
               w->selected_response[selected++] = o;
         if (!selected)
            continue;
         /* Refit only the responses that actually selected this design.
            This also prevents irrelevant sparse fits from failing convergence.
          */
         for (i = 0; i < d->n; i++)
            for (source = 0; source < selected; source++)
               w->selected_y[(size_t)i * selected + source] =
                   y[(size_t)i * q + w->selected_response[source]];
         if (w->repeat[c]) {
            if (!THD_picard_repeat_fit(w->repeat[c], w->selected_y, selected,
                                       err, nerr))
               return 0;
            repeat = THD_picard_repeat_get_result(w->repeat[c]);
            fit = repeat->fit;
         } else {
            for (source = 0; source < selected; source++)
               w->parameter[source] = p->config.fixed_parameter;
            if (!THD_picard_fit_batch(w->fixed[c], w->selected_y, selected,
                                      w->parameter, err, nerr))
               return 0;
            fit = THD_picard_get_result(w->fixed[c]);
         }
         /* Scatter the compact fit back into original batch-column order.
            selected_response is the inverse map built for this design. */
         for (source = 0; source < selected; source++) {
            o = w->selected_response[source];
            candidate_retain_fit(w, fit, source, o, q);
            if (repeat) {
               w->parameter_index[o] = repeat->candidate[source];
               w->loss[o] = repeat->loss[source];
               w->baseline[o] = repeat->reference_loss[source];
               for (g = 0; g < p->config.nparameter; g++)
                  w->parameter_loss[(size_t)g * q + o] =
                      repeat->candidate_loss[(size_t)g * selected + source];
            }
         }
      }
   /* Publish one composite fit: each response may come from a different
      design, but all candidates share coefficient identities and dimensions.
      NULL fields below mark diagnostics that this policy did not compute. */
   w->fit = (THD_picard_result){
       d->n,
       d->p,
       d->k,
       q,
       w->beta,
       w->gamma,
       w->fitted,
       w->residual,
       d->elastic || d->banded ? NULL : w->fraction_diagnostic,
       d->elastic ? w->elastic_diagnostic : NULL,
       d->banded ? w->banded_diagnostic : NULL,
       d->support_refit ? w->support_size : NULL,
       d->support_refit ? w->support_rank : NULL};
   w->result = (THD_picard_candidates_result){
       q,
       &w->fit,
       w->index,
       w->id,
       w->rank,
       w->selection_loss,
       forced ? NULL : w->design_loss,
       p->config.nparameter ? w->parameter_index : NULL,
       p->config.policy == THD_PICARD_SELECT_PREDICT ? w->status : NULL,
       p->config.nparameter ? w->loss : NULL,
       p->config.nparameter ? w->baseline : NULL,
       p->config.policy == THD_PICARD_SELECT_PREDICT ? w->r2 : NULL,
       p->config.nparameter ? w->parameter_loss : NULL};
   w->valid = 1;
   return 1;
}

/* A staged workflow locks HRFs before learning noise PCs. Refitting must
   not silently reselect HRFs after adding those nuisance regressors. */
int THD_picard_candidates_fit_selected(THD_picard_candidates_workspace *w,
                                       const double *y, int q,
                                       const int *indices, const double *scores,
                                       char *err, int nerr)
{
   if (!indices) {
      if (w)
         w->valid = 0;
      return candidate_error(err, nerr, "missing locked HRFs");
   }
   return candidate_fit(w, y, q, indices, scores, err, nerr);
}
int THD_picard_candidates_fit(THD_picard_candidates_workspace *w,
                              const double *y, int q, char *err, int nerr)
{
   return candidate_fit(w, y, q, NULL, NULL, err, nerr);
}

const THD_picard_candidates_result *
THD_picard_candidates_get_result(const THD_picard_candidates_workspace *w)
{
   return w && w->valid ? &w->result : NULL;
}

/* ========================================================================
   5. NOISE PCS AND SINGLE-TRIAL PC SELECTION
   ======================================================================== */

/* Single-trial orchestration primitives. Dataset I/O and voxel-mask policy
   stay in 3dPicard; this file owns only numeric preprocessing and scoring.
   Temporal covariance, rather than spatial loadings, bounds PCA memory by
   run length squared. Use the existing SVD, not a second eigensolver. */

/* first/n index one contiguous run in retained-row space. gram accumulates
   outer products; batch is reused for gathered and residualized responses.
   Each run owns its own projector and scratch, so no cross-run data enter PCA.
 */
typedef struct {
   int first, n, usable, rank;
   double *gram, *batch;
   THD_linproj *projection;
   THD_linproj_workspace *work;
} NoiseRun;
struct THD_picard_noise {
   const THD_picard_plan *base;
   int nrun, maxpc, batch, valid, failed;
   NoiseRun *run;
   double *pcs;
   int *ranks;
};
static int noise_error(char *err, int nerr, const char *message)
{
   if (err && nerr > 0)
      snprintf(err, (size_t)nerr, "%s", message);
   return 0;
}
static double *noise_array(size_t n, size_t p)
{
   if (!p)
      p = 1;
   return n > SIZE_MAX / sizeof(double) / p ? NULL
                                            : calloc(n * p, sizeof(double));
}
void THD_picard_noise_free(THD_picard_noise *s)
{
   int r;
   if (!s)
      return;
   if (s->run)
      for (r = 0; r < s->nrun; r++) {
         free(s->run[r].gram);
         free(s->run[r].batch);
         THD_linproj_workspace_free(s->run[r].work);
         THD_linproj_free(s->run[r].projection);
      }
   free(s->run);
   free(s->pcs);
   free(s->ranks);
   free(s);
}
/* Construct only the run layout and nuisance projectors, before seeing data.
   Include a constant even when the supplied nuisance span lacks one: PCA
   here describes fluctuations, not each pool voxel's baseline intensity.
   Redundant constants are safe because the projector is rank revealing. */
THD_picard_noise *THD_picard_noise_create(const THD_picard_plan *p, int maxpc,
                                          int batch, char *err, int nerr)
{
   THD_picard_noise *s;
   int i, j, r = 0, ok = 1;
   if (!p || maxpc < 0 || batch < 1) {
      noise_error(err, nerr, "invalid noise PCA configuration");
      return NULL;
   }
   s = calloc(1, sizeof(*s));
   if (!s) {
      noise_error(err, nerr, "cannot allocate noise PCA");
      return NULL;
   }
   s->base = p;
   s->maxpc = maxpc;
   s->batch = batch;
   for (i = 0; i < p->info.n; i++)
      s->nrun += !i || p->runs[i] != p->runs[i - 1];
   s->run = calloc((size_t)s->nrun, sizeof(*s->run));
   s->ranks = calloc((size_t)s->nrun, sizeof(int));
   s->pcs = noise_array(p->info.n, maxpc);
   if (!s->run || !s->ranks || !s->pcs)
      ok = 0;
   for (i = 0; i < p->info.n && ok; i++) {
      if (i && p->runs[i] != p->runs[i - 1])
         r++;
      if (!s->run[r].n)
         s->run[r].first = i;
      s->run[r].n++;
   }
   for (r = 0; r < s->nrun && ok; r++) {
      NoiseRun *run = s->run + r;
      int k = p->info.k + 1;
      double *z = noise_array(run->n, k);
      run->gram = noise_array(run->n, run->n);
      run->batch = noise_array(run->n, batch);
      if (!z || !run->gram || !run->batch)
         ok = 0;
      if (ok) {
         for (i = 0; i < run->n; i++) {
            z[(size_t)i * k] = 1;
            for (j = 0; j < p->info.k; j++)
               z[(size_t)i * k + j + 1] =
                   p->z[(size_t)(run->first + i) * p->info.k + j];
         }
         run->projection =
             THD_linproj_prepare(z, run->n, k, p->info.rank_tol, err, nerr);
         if (run->projection)
            run->work =
                THD_linproj_workspace_alloc(run->projection, batch, err, nerr);
         ok = run->work != NULL &&
              run->n - THD_linproj_rank(run->projection) >= maxpc;
      }
      free(z);
   }
   if (!ok) {
      THD_picard_noise_free(s);
      noise_error(
          err, nerr,
          "noise PCA allocation or residual degrees of freedom insufficient");
      return NULL;
   }
   return s;
}
/* Normalize each pool voxel separately within each run. Zero residuals cast
   no vote. Reconstructing covariance in voxel order makes results independent
   of the adapter's batch capacity. Failed input permanently invalidates state.
 */
int THD_picard_noise_add(THD_picard_noise *s, const double *y, int q, char *err,
                         int nerr)
{
   int r, i, j, o;
   if (!s || !y || q < 1 || q > s->batch || s->valid || s->failed)
      return noise_error(err, nerr, "invalid noise PCA accumulation state");
   s->failed = 1;
   for (r = 0; r < s->nrun; r++) {
      NoiseRun *run = s->run + r;
      for (i = 0; i < run->n; i++)
         for (o = 0; o < q; o++)
            run->batch[(size_t)i * q + o] = y[(size_t)(run->first + i) * q + o];
      if (!THD_linproj_project(run->work, run->batch, q, run->batch, err, nerr))
         return 0;
      for (o = 0; o < q; o++) {
         double norm = 0;
         for (i = 0; i < run->n; i++)
            norm = hypot(norm, run->batch[(size_t)i * q + o]);
         if (!isfinite(norm))
            return noise_error(err, nerr, "noise norm exceeds range");
         if (!norm)
            continue;
         if (run->usable == INT_MAX)
            return noise_error(err, nerr, "too many noise voxels");
         run->usable++;
         for (i = 0; i < run->n; i++)
            for (j = 0; j <= i; j++)
               run->gram[(size_t)i * run->n + j] +=
                   (run->batch[(size_t)i * q + o] / norm) *
                   (run->batch[(size_t)j * q + o] / norm);
      }
   }
   s->failed = 0;
   return 1;
}
/* Convert accumulated covariance into ordered temporal components. Gram
   eigenvalues are squared singular values of normalized pool responses, so
   square the design tolerance. A roundoff floor prevents tiny Gram artifacts
   from masquerading as additional usable PCs. SVD order is not guaranteed by
   the underlying routine; sort an index array instead of moving its columns. */
int THD_picard_noise_finish(THD_picard_noise *s, char *err, int nerr)
{
   int r, i, j, c;
   if (!s || s->valid || s->failed)
      return noise_error(err, nerr, "invalid PCA finish state");
   s->failed = 1;
   for (r = 0; r < s->nrun; r++) {
      NoiseRun *run = s->run + r;
      THD_numeric_svd *svd;
      int *order;
      double tol = fmax(s->base->info.rank_tol * s->base->info.rank_tol,
                        64 * DBL_EPSILON * run->n);
      if (!run->usable && s->maxpc)
         return noise_error(err, nerr, "noise pool has zero run rank");
      for (i = 0; i < run->n; i++)
         for (j = 0; j < i; j++)
            run->gram[(size_t)j * run->n + i] =
                run->gram[(size_t)i * run->n + j];
      svd = THD_numeric_svd_prepare(run->gram, run->n, run->n, tol, err, nerr);
      if (!svd)
         return 0;
      order = calloc((size_t)run->n, sizeof(int));
      if (!order) {
         THD_numeric_svd_free(svd);
         return noise_error(err, nerr, "cannot sort PCs");
      }
      for (i = 0; i < run->n; i++)
         order[i] = i;
      for (i = 0; i < run->n; i++)
         for (j = i + 1; j < run->n; j++)
            if (svd->s[order[j]] > svd->s[order[i]]) {
               int swap = order[i];
               order[i] = order[j];
               order[j] = swap;
            }
      run->rank = s->ranks[r] = svd->rank;
      if (run->rank < s->maxpc) {
         free(order);
         THD_numeric_svd_free(svd);
         return noise_error(err, nerr,
                            "requested PCs exceed noise rank in a run");
      }
      /* An SVD vector's sign is arbitrary. A deterministic pivot makes saved
         PC provenance and nuisance coefficients comparable across runs of the
         program; the fitted subspace itself does not depend on this sign. */
      for (c = 0; c < s->maxpc; c++) {
         int index = order[c], pivot = 0;
         double mean = 0, sd = 0, scale;
         for (i = 0; i < run->n; i++) {
            double value = svd->u[(size_t)i * run->n + index];
            mean += value / run->n;
            if (fabs(value) > fabs(svd->u[(size_t)pivot * run->n + index]))
               pivot = i;
         }
         for (i = 0; i < run->n; i++)
            sd = hypot(sd, svd->u[(size_t)i * run->n + index] - mean);
         scale = copysign(sqrt(run->n - 1) / sd,
                          svd->u[(size_t)pivot * run->n + index]);
         if (!isfinite(scale)) {
            free(order);
            THD_numeric_svd_free(svd);
            return noise_error(err, nerr, "undefined noise PC scale");
         }
         for (i = 0; i < run->n; i++)
            s->pcs[(size_t)(run->first + i) * s->maxpc + c] =
                svd->u[(size_t)i * run->n + index] * scale;
      }
      free(order);
      THD_numeric_svd_free(svd);
   }
   s->failed = 0;
   s->valid = 1;
   return 1;
}
const double *THD_picard_noise_pcs(const THD_picard_noise *s)
{
   return s && s->valid ? s->pcs : NULL;
}
int THD_picard_noise_runs(const THD_picard_noise *s)
{
   return s ? s->nrun : 0;
}
const int *THD_picard_noise_ranks(const THD_picard_noise *s)
{
   return s && s->valid ? s->ranks : NULL;
}

/* Rebuild through the ordinary plan constructor, rather than altering an
   existing factorization. Original X and Z are scattered back to their source
   columns. Appended columns are zero outside their own run and have role 0,
   so penalties still address exactly the original task coefficients. */
THD_picard_plan *THD_picard_noise_plan(const THD_picard_noise *s,
                                       const THD_picard_plan *p, int count,
                                       int ols, char *err, int nerr)
{
   THD_picard_plan *result = NULL;
   double *matrix;
   int *roles, i, j, r, d;
   THD_picard_config cfg;
   if (!s || !s->valid || !p || count < 0 || count > s->maxpc ||
       p->info.n != s->base->info.n ||
       memcmp(p->rows, s->base->rows, sizeof(int) * (size_t)p->info.n) ||
       memcmp(p->runs, s->base->runs, sizeof(int) * (size_t)p->info.n) ||
       count > (INT_MAX - p->info.ncol) / s->nrun) {
      noise_error(err, nerr, "invalid noise-augmented design");
      return NULL;
   }
   d = p->info.ncol + s->nrun * count;
   matrix = noise_array(p->info.n, d);
   roles = calloc((size_t)d, sizeof(int));
   if (!matrix || !roles) {
      free(matrix);
      free(roles);
      noise_error(err, nerr, "cannot allocate PC design");
      return NULL;
   }
   for (j = 0; j < p->info.p; j++)
      roles[p->pen[j]] = 1;
   for (i = 0, r = 0; i < p->info.n; i++) {
      if (i && p->runs[i] != p->runs[i - 1])
         r++;
      for (j = 0; j < p->info.p; j++)
         matrix[(size_t)i * d + p->pen[j]] = p->x[(size_t)i * p->info.p + j];
      for (j = 0; j < p->info.k; j++)
         matrix[(size_t)i * d + p->nuis[j]] = p->z[(size_t)i * p->info.k + j];
      for (j = 0; j < count; j++)
         matrix[(size_t)i * d + p->info.ncol + r * count + j] =
             s->pcs[(size_t)i * s->maxpc + j];
   }
   cfg = (THD_picard_config){!ols && p->elastic ? &p->config : NULL,
                             !ols && p->banded ? &p->band_config : NULL,
                             p->info.ar1, !ols && p->info.support_refit,
                             p->info.covariance_enabled ?
                                 &p->info.covariance : NULL};
   result = THD_picard_plan_create_config(matrix, p->info.n, d, roles, p->rows,
                                          p->runs, p->info.nfull,
                                          p->info.rank_tol, &cfg, err, nerr);
   free(matrix);
   free(roles);
   return result;
}
/* The run-local contract makes full-data unpenalized fits separable by run.
   Training betas therefore do not require held-out responses. Noise/HRF
   choices are nevertheless learned upstream, so the score is conditional. */
int THD_picard_single_trials(const THD_picard_plan *p,
                             const THD_picard_trial *trials, char *err,
                             int nerr)
{
   int i, j, k, run;
   if (!p || !trials)
      return noise_error(err, nerr, "missing single-trial identities");
   for (j = 0; j < p->info.p; j++) {
      if (trials[j].column != p->pen[j] || trials[j].trial < 0 ||
          trials[j].condition < 0 || trials[j].run < 0)
         return noise_error(err, nerr, "invalid single-trial identity");
      for (k = 0; k < j; k++)
         if (trials[j].trial == trials[k].trial)
            return noise_error(err, nerr, "duplicate trial identity");
      run = -1;
      for (i = 0; i < p->info.n; i++)
         if (p->x[(size_t)i * p->info.p + j] != 0) {
            if (p->runs[i] != trials[j].run)
               return noise_error(err, nerr, "trial not run-local");
            run = p->runs[i];
         }
      if (run < 0)
         return noise_error(err, nerr, "trial has no retained support");
   }
   for (j = 0; j < p->info.k; j++) {
      run = -1;
      for (i = 0; i < p->info.n; i++)
         if (p->z[(size_t)i * p->info.k + j] != 0) {
            if (run >= 0 && run != p->runs[i])
               return noise_error(err, nerr,
                                  "single-trial nuisance must be run-local");
            run = p->runs[i];
         }
   }
   return 1;
}
/* Compare every ordered repeated-condition pair. Both the reference mean
   and SD are fixed across candidate counts; subtracting that same mean from
   both betas cancels in their difference. Long-double accumulation reduces
   cancellation/overflow while the public result remains ordinary double. */
int THD_picard_single_loss(const THD_picard_trial *t, int p,
                           const double *reference, const double *beta, int q,
                           double *loss, char *err, int nerr)
{
   int i, j, o;
   long double pairs = 0;
   if (!t || !reference || !beta || !loss || p < 2 || q < 1)
      return noise_error(err, nerr, "invalid PC scoring inputs");
   for (i = 0; i < p; i++)
      for (j = 0; j < p; j++)
         pairs += t[i].run != t[j].run && t[i].condition == t[j].condition;
   if (!pairs)
      return noise_error(err, nerr,
                         "PC selection requires cross-run condition repeats");
   for (o = 0; o < q; o++) {
      long double mean = 0, variance = 0, sum = 0;
      for (i = 0; i < p; i++) {
         if (!isfinite(reference[(size_t)i * q + o]) ||
             !isfinite(beta[(size_t)i * q + o]))
            return noise_error(err, nerr, "nonfinite PC scoring coefficients");
         mean += (long double)reference[(size_t)i * q + o] / p;
      }
      for (i = 0; i < p; i++) {
         long double delta = reference[(size_t)i * q + o] - mean;
         variance += delta * delta / (p - 1);
         for (j = 0; j < p; j++)
            if (t[i].run != t[j].run && t[i].condition == t[j].condition) {
               delta = (long double)beta[(size_t)i * q + o] -
                       reference[(size_t)j * q + o];
               sum += delta * delta;
            }
      }
      loss[o] = variance > 0 ? (double)(sum / pairs / variance) : -1;
      if (!isfinite(loss[o]))
         return noise_error(err, nerr, "PC selection loss exceeds range");
   }
   return 1;
}
/* Increasing counts encode the simplicity preference. Stop at the first
   count sufficiently close to the best improvement, not simply its argmin.
   In particular, a flat or worsening curve chooses the zero-PC baseline. */
int THD_picard_single_select(const double *loss, int count, double stop)
{
   int c;
   double best = 0;
   if (!loss || count < 1 || !isfinite(stop) || stop < 1)
      return -1;
   for (c = 0; c < count; c++) {
      if (!isfinite(loss[c]) || loss[c] < 0)
         return -1;
      best = fmax(best, loss[0] - loss[c]);
   }
   for (c = 0; c < count; c++)
      if ((loss[0] - loss[c]) * stop >= best)
         return c;
   return 0;
}
