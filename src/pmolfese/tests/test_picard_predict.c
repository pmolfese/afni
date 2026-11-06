#include "thd_crossval.h"
#include "thd_picard.h"
#include <math.h>

#define CHECK(c, msg)                                                          \
   do {                                                                        \
      if (!(c)) {                                                              \
         fprintf(stderr, "line %d: %s\n", __LINE__, msg);                      \
         return 1;                                                             \
      }                                                                        \
   } while (0)

/* Count generic traversal events and deliberately fail a path callback. */
typedef struct {
   int prepared, paths, finished, reject;
} Visits;
static int visit_prepare(void *ctx, int f, const int *train, int n,
                         const int *test, int nt, char *err, int ne)
{
   Visits *v = ctx;
   (void)f;
   (void)train;
   (void)test;
   (void)err;
   (void)ne;
   if (!n || !nt)
      return 0;
   ++v->prepared;
   return 1;
}
static int visit_path(void *ctx, int f, char *err, int ne)
{
   Visits *v = ctx;
   (void)f;
   (void)err;
   (void)ne;
   ++v->paths;
   return !v->reject;
}
static int visit_finish(void *ctx, char *err, int ne)
{
   Visits *v = ctx;
   (void)err;
   (void)ne;
   ++v->finished;
   return 1;
}
static int numeric_folds(void)
{
   int rows[] = {1, 4, 5, 8, 11, 13, 17, 19}, runs[] = {0, 0, 0, 0, 1, 1, 1, 1};
   int starts[] = {0, 10}, groups[] = {1, 1, 0, 0, 1, 1, 0, 0};
   char err[256];
   THD_foldplan *p = THD_foldplan_time_series(8, rows, runs, 2, starts, NULL, 4,
                                              1, err, sizeof(err));
   THD_foldplan *copy;
   Visits v = {0};
   const THD_cv_path_ops ops = {visit_prepare, visit_path, visit_finish};
   CHECK(p, err);
   CHECK(p->nfold == 6 && p->test[0][0] == 0,
         "blocks use original run-relative times");
   CHECK(p->train[0][0] == 1,
         "gap is not measured in compressed retained rows");
   copy = THD_foldplan_copy(p, err, sizeof(err));
   CHECK(copy, err);
   p->train[0][0] = 7;
   CHECK(copy->train[0][0] == 1, "copy owns excluded/training indices");
   THD_foldplan_free(p);
   CHECK(THD_cv_path_search(copy, &ops, &v, err, sizeof(err)), err);
   CHECK(v.prepared == 6 && v.paths == 6 && v.finished == 1,
         "one path callback per fold");
   v = (Visits){0, 0, 0, 1};
   CHECK(!THD_cv_path_search(copy, &ops, &v, err, sizeof(err)) &&
            v.prepared == 1 && v.paths == 1 && v.finished == 0,
         "path abort propagation");
   THD_foldplan_free(copy);
   CHECK(!THD_foldplan_time_series(8, rows, runs, 2, starts, groups, 0, 100,
                                   err, sizeof(err)),
         "reject gaps that empty training sets");
   return 0;
}
static int known_prediction(void)
{
   double x[] = {1, 1, 1, 1}, y[] = {2, 0, 2, 0, 4, 0, 4, 0},
          grid[] = {1, .5, .8, 0};
   int roles[] = {1}, rows[] = {0, 1, 2, 3}, runs[] = {0, 0, 1, 1},
       starts[] = {0, 2};
   char err[256];
   THD_foldplan *folds = THD_foldplan_time_series(4, rows, runs, 2, starts,
                                                  NULL, 0, 0, err, sizeof(err));
   THD_picard_plan *full = THD_picard_plan_create(x, 4, 1, roles, rows, runs, 4,
                                                  1e-10, err, sizeof(err));
   THD_picard_predict_plan *p;
   THD_picard_predict_workspace *w;
   const THD_picard_predict_result *r;
   CHECK(folds && full, err);
   p = THD_picard_predict_plan_create(full, folds, grid, 4, 0, 0, err,
                                      sizeof(err));
   CHECK(p, err);
   THD_foldplan_free(folds); /* Prediction plan copied its indices. */
   w = THD_picard_predict_workspace_alloc(p, 2, err, sizeof(err));
   CHECK(w, err);
   CHECK(THD_picard_predict_fit(w, y, 2, 0, -1, err, sizeof(err)), err);
   r = THD_picard_predict_get_result(w);
   CHECK(r->candidate[0] == 2 && fabs(r->loss[0] - 3.6) < 1e-8,
         "known MSE winner");
   CHECK(fabs(r->r2[0] - .1) < 1e-8 && r->baseline[0] == 4,
         "training mean baseline");
   CHECK(fabs(r->fit->beta[0] - 2.4) < 1e-8,
         "full-design refit, not a fold fit");
   CHECK(r->status[1] == 1 && isnan(r->r2[1]) && r->candidate[1] == 0,
         "zero-baseline policy");
   {
      THD_picard_predict_global *g =
         THD_picard_predict_global_alloc(p, err, sizeof(err));
      const THD_picard_predict_global_result *global;
      CHECK(g, err);
      CHECK(THD_picard_predict_global_add(g, w, err, sizeof(err)), err);
      global = THD_picard_predict_global_finish(g, err, sizeof(err));
      CHECK(global && global->candidate == 2 && global->nvalid == 1,
            "global selection excludes zero-baseline responses");
      CHECK(fabs(global->fold_loss[2] - .36) < 1e-8,
            "global baseline normalization");
      CHECK(!THD_picard_predict_global_add(g, w, err, sizeof(err)),
            "finished global accumulators reject more batches");
      THD_picard_predict_global_free(g);
      g = THD_picard_predict_global_alloc(p, err, sizeof(err));
      CHECK(g, err);
      CHECK(!THD_picard_predict_global_add(g, NULL, err, sizeof(err)) &&
               !THD_picard_predict_global_finish(g, err, sizeof(err)),
            "failed global batch cannot expose a partial selection");
      THD_picard_predict_global_free(g);
   }
   y[0] = 200;
   y[2] = 200;
   CHECK(THD_picard_predict_fit(w, y, 2, 1, -1, err, sizeof(err)), err);
   r = THD_picard_predict_get_result(w);
   CHECK(!r->fit && fabs(r->fold_loss[2] - 198 * 198) < 1e-6,
         "held-out perturbation does not alter training fractional fit");
   y[0] = 2;
   y[2] = 2;
   CHECK(THD_picard_predict_fit(w, y, 2, 0, 1, err, sizeof(err)), err);
   CHECK(fabs(THD_picard_predict_get_result(w)->fit->beta[0] - 1.5) < 1e-8,
         "forced global candidate refit");
   y[0] = NAN;
   CHECK(!THD_picard_predict_fit(w, y, 2, 0, -1, err, sizeof(err)) &&
            !THD_picard_predict_get_result(w),
         "failed prediction invalidates result");
   THD_picard_predict_workspace_free(w);
   THD_picard_predict_plan_free(p);
   folds = THD_foldplan_time_series(4, rows, runs, 2, starts, NULL, 0, 0, err,
                                    sizeof(err));
   p = THD_picard_predict_plan_create(full, folds, grid, 4, 0, 1, err,
                                      sizeof(err));
   CHECK(p, err);
   {
      double loss[] = {0, 5, 4, 8, 8, 5, 4, 8}, se;
      int winner;
      CHECK(THD_picard_predict_select(p, loss, &winner, &se, err, sizeof(err)),
            err);
      CHECK(winner == 3 && fabs(se - 4) < 1e-12,
            "one-SE chooses strongest eligible fraction");
   }
   THD_picard_predict_plan_free(p);
   THD_foldplan_free(folds);
   THD_picard_plan_free(full);
   return 0;
}
static int conditional_prediction(void)
{
   double x[] = {1, 0, -1, 1, 0, 1, 0, 1, -1, 0, 1, 1};
   double y[] = {18, 22, 36, 44}, grid[] = {1, .5, .8, 0};
   int roles[] = {0, 0, 1}, rows[] = {0, 1, 2, 3}, runs[] = {0, 0, 1, 1},
       starts[] = {0, 2};
   char err[256];
   THD_foldplan *folds = THD_foldplan_time_series(4, rows, runs, 2, starts,
                                                  NULL, 0, 0, err, sizeof(err));
   THD_picard_plan *full = THD_picard_plan_create(x, 4, 3, roles, rows, runs, 4,
                                                  1e-10, err, sizeof(err));
   THD_picard_predict_plan *p;
   THD_picard_predict_workspace *w;
   const THD_picard_predict_result *r;
   CHECK(folds && full, err);
   CHECK(!THD_picard_predict_plan_create(full, folds, grid, 4, 0, 0, err,
                                         sizeof(err)) &&
            strstr(err, "not estimable"),
         "raw prediction rejects unseen run intercept");
   p = THD_picard_predict_plan_create(full, folds, grid, 4, 1, 0, err,
                                      sizeof(err));
   CHECK(p, err);
   w = THD_picard_predict_workspace_alloc(p, 1, err, sizeof(err));
   CHECK(w, err);
   CHECK(THD_picard_predict_fit(w, y, 1, 0, -1, err, sizeof(err)), err);
   r = THD_picard_predict_get_result(w);
   CHECK(r->candidate[0] == 2 && fabs(r->loss[0] - 3.6) < 1e-8 &&
            fabs(r->baseline[0] - 10) < 1e-10 && fabs(r->r2[0] - .64) < 1e-8,
         "conditional score projects test signal and baseline consistently");
   THD_picard_predict_workspace_free(w);
   THD_picard_predict_plan_free(p);
   THD_foldplan_free(folds);
   THD_picard_plan_free(full);
   return 0;
}
int main(void)
{
   int failed = numeric_folds() + known_prediction() + conditional_prediction();
   if (!failed)
      puts("PASS numeric folds, traversal, prediction isolation, selection and "
           "nuisance scoring");
   return failed;
}
