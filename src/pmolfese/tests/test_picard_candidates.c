#include "thd_picard.h"
#include "thd_crossval.h"
#include <math.h>
#include <stdio.h>

#define CHECK(c, m)                                                            \
   do {                                                                        \
      if (!(c)) {                                                              \
         fprintf(stderr, "line %d: %s\n", __LINE__, m);                        \
         return 1;                                                             \
      }                                                                        \
   } while (0)

/* Both designs are equally sized but explain different response columns.
   Known slopes verify that selecting and compacting responses preserves
   coefficient/response identities, including ties and a short last batch. */
int main(void)
{
   double first[16], second[16], y[24], grid[] = {1, .5};
   int roles[] = {0, 1}, rows[8], runs[8], starts[] = {0, 4}, ids[] = {7, 2};
   char err[512];
   THD_picard_plan *a, *b;
   const THD_picard_plan *designs[2];
   THD_picard_candidates_config cfg = {0};
   THD_foldplan *folds;
   int i, policy;
   for (i = 0; i < 8; i++) {
      first[2 * i] = second[2 * i] = 1;
      first[2 * i + 1] = i % 4 - 1;
      second[2 * i + 1] = i % 2 ? -1 : 1;
      y[3 * i] = 10 + 2 * first[2 * i + 1];
      y[3 * i + 1] = 5 + 3 * second[2 * i + 1];
      y[3 * i + 2] = 0;
      rows[i] = i;
      runs[i] = i / 4;
   }
   a = THD_picard_plan_create(first, 8, 2, roles, rows, runs, 8, 1e-10, err,
                              512);
   b = THD_picard_plan_create(second, 8, 2, roles, rows, runs, 8, 1e-10, err,
                              512);
   CHECK(a && b, err);
   designs[0] = a;
   designs[1] = b;
   folds =
       THD_foldplan_time_series(8, rows, runs, 2, starts, NULL, 0, 0, err, 512);
   CHECK(folds, err);
   for (policy = 0; policy < 2; policy++) {
      THD_picard_candidates *p;
      THD_picard_candidates_workspace *w;
      const THD_picard_candidates_result *r;
      cfg.policy = policy;
      cfg.fixed_parameter = .5;
      cfg.nparameter = policy ? 2 : 0;
      cfg.parameters = policy ? grid : NULL;
      cfg.folds = policy ? folds : NULL;
      p = THD_picard_candidates_create(designs, ids, 2, &cfg, err, 512);
      CHECK(p, err);
      w = THD_picard_candidates_workspace_alloc(p, 3, err, 512);
      CHECK(w, err);
      CHECK(!THD_picard_candidates_get_result(w), "unfitted result invalid");
      CHECK(THD_picard_candidates_fit(w, y, 3, err, 512), err);
      r = THD_picard_candidates_get_result(w);
      CHECK(r->design_id[0] == 7 && r->design_id[1] == 2 &&
                r->design_id[2] == 2,
            "response-wise winner and stable-ID tie");
      CHECK(fabs(r->fit->beta[0] - (policy ? 2 : 1)) < 1e-8 &&
                fabs(r->fit->beta[1] - (policy ? 3 : 1.5)) < 1e-8,
            "winner coefficients in correct response slots");
      CHECK(r->fit->beta[2] == 0, "zero response coefficient");
      CHECK(r->selection_loss[0] < 1e-20 && r->selection_loss[1] < 1e-20,
            "known perfect OLS/prediction scores");
      if (!policy) {
         int forced[] = {1, 0, 0};
         double scores[] = {4, 5, 6};
         CHECK(THD_picard_candidates_fit_selected(w, y, 3, forced, scores, err,
                                                  512),
               err);
         r = THD_picard_candidates_get_result(w);
         CHECK(r->design_id[0] == 2 && r->design_id[1] == 7 && !r->design_loss,
               "locked HRFs bypass reselection and do not expose stale grid");
         CHECK(r->selection_loss[1] == 5, "locked pre-PC score preserved");
         forced[0] = 2;
         CHECK(!THD_picard_candidates_fit_selected(w, y, 3, forced, scores, err,
                                                   512) &&
                   !THD_picard_candidates_get_result(w),
               "bad lock invalidates result");
      }
      y[0] = NAN;
      CHECK(!THD_picard_candidates_fit(w, y, 3, err, 512) &&
                !THD_picard_candidates_get_result(w),
            "failure invalidates selection");
      y[0] = 8;
      CHECK(THD_picard_candidates_fit(w, y, 3, err, 512), err);
      {
         double short_y[8];
         for (i = 0; i < 8; i++)
            short_y[i] = y[3 * i + 1];
         CHECK(THD_picard_candidates_fit(w, short_y, 1, err, 512), err);
         CHECK(THD_picard_candidates_get_result(w)->design_id[0] == 2,
               "changed active batch stride");
      }
      THD_picard_candidates_workspace_free(w);
      THD_picard_candidates_free(p);
   }
   ids[1] = ids[0];
   CHECK(!THD_picard_candidates_create(designs, ids, 2, &cfg, err, 512),
         "duplicate ID rejected");
   ids[1] = 16777216;
   CHECK(!THD_picard_candidates_create(designs, ids, 2, &cfg, err, 512),
         "inexact float ID rejected");
   ids[1] = 2;
   {
      THD_picard_plan *different;
      second[0] = 2; /* Change the nuisance intercept, not the stimulus. */
      different = THD_picard_plan_create(second, 8, 2, roles, rows, runs, 8,
                                         1e-10, err, 512);
      CHECK(different, err);
      designs[1] = different;
      CHECK(!THD_picard_candidates_create(designs, ids, 2, &cfg, err, 512),
            "different nuisance values rejected");
      THD_picard_plan_free(different);
   }
   THD_foldplan_free(folds);
   THD_picard_plan_free(a);
   THD_picard_plan_free(b);
   puts("PASS numeric candidate selection, response identity, ties and failure "
        "state");
   return 0;
}
