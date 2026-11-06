/* Numerical contracts independent of AFNI I/O: analytic rank-one PCs,
   batching, undefined scores and poisoned/partial workspace cleanup. */
#include "thd_picard.h"
#include <math.h>
#include <stdio.h>
#define CHECK(c)                                                               \
   do {                                                                        \
      if (!(c)) {                                                              \
         fprintf(stderr, "line %d: %s (%s)\n", __LINE__, #c, err);             \
         return 1;                                                             \
      }                                                                        \
   } while (0)
int main(void)
{
   char err[256] = "";
   double x[32] = {0}, y[16], pc[] = {-3, -1, 1, 3};
   int rows[8], runs[8], roles[] = {1, 1, 0, 0}, i, j;
   THD_picard_trial trials[] = {{0, 0, 0, 0}, {1, 1, 0, 1}};
   THD_picard_plan *plan, *aug;
   THD_picard_noise *noise, *bad;
   const double *values;
   for (i = 0; i < 8; i++) {
      rows[i] = i;
      runs[i] = i / 4;
      x[4 * i + i / 4] = i % 2 ? 1 : -1;
      x[4 * i + 2 + i / 4] = 1;
      y[2 * i] = 5 + pc[i % 4];
      y[2 * i + 1] = 9 + 2 * pc[i % 4];
   }
   plan =
       THD_picard_plan_create(x, 8, 4, roles, rows, runs, 8, 1e-10, err, 256);
   CHECK(plan);
   CHECK(THD_picard_single_trials(plan, trials, err, 256));
   noise = THD_picard_noise_create(plan, 1, 2, err, 256);
   CHECK(noise);
   CHECK(!THD_picard_noise_pcs(noise));
   CHECK(THD_picard_noise_add(noise, y, 2, err, 256));
   CHECK(THD_picard_noise_finish(noise, err, 256));
   values = THD_picard_noise_pcs(noise);
   CHECK(values);
   for (i = 0; i < 8; i++)
      CHECK(fabs(fabs(values[i]) - fabs(pc[i % 4]) / sqrt(20.0 / 3)) < 1e-10);
   CHECK(THD_picard_noise_runs(noise) == 2);
   CHECK(THD_picard_noise_ranks(noise)[0] == 1);
   aug = THD_picard_noise_plan(noise, plan, 1, 1, err, 256);
   CHECK(aug);
   CHECK(THD_picard_plan_info(aug)->k == 4);
   CHECK(THD_picard_single_trials(aug, trials, err, 256));
   THD_picard_plan_free(aug);
   bad = THD_picard_noise_create(plan, 2, 2, err, 256);
   CHECK(bad);
   CHECK(THD_picard_noise_add(bad, y, 2, err, 256));
   CHECK(!THD_picard_noise_finish(bad, err, 256));
   CHECK(!THD_picard_noise_pcs(bad));
   THD_picard_noise_free(bad);
   bad = THD_picard_noise_create(plan, 1, 1, err, 256);
   CHECK(bad);
   for (j = 0; j < 2; j++) {
      double batch[8];
      for (i = 0; i < 8; i++)
         batch[i] = y[2 * i + j];
      CHECK(THD_picard_noise_add(bad, batch, 1, err, 256));
   }
   CHECK(THD_picard_noise_finish(bad, err, 256));
   for (i = 0; i < 8; i++)
      CHECK(fabs(values[i] - THD_picard_noise_pcs(bad)[i]) < 1e-12);
   THD_picard_noise_free(bad);
   bad = THD_picard_noise_create(plan, 1, 2, err, 256);
   CHECK(bad);
   y[1] = NAN;
   CHECK(!THD_picard_noise_add(bad, y, 2, err, 256));
   CHECK(!THD_picard_noise_finish(bad, err, 256));
   THD_picard_noise_free(bad);
   {
      double ref[] = {1, 0, 3, 0}, beta[] = {2, 0, 2, 0}, loss[2];
      double trend[] = {10, 6, 5.9, 7}, flat[] = {3, 3, 4};
      CHECK(THD_picard_single_loss(trials, 2, ref, beta, 2, loss, err, 256));
      CHECK(fabs(loss[0] - .5) < 1e-12 && loss[1] == -1);
      CHECK(THD_picard_single_select(trend, 4, 1.05) == 1);
      CHECK(THD_picard_single_select(trend, 4, 1) == 2);
      CHECK(THD_picard_single_select(flat, 3, 1.05) == 0);
      trials[1].condition = 1;
      CHECK(!THD_picard_single_loss(trials, 2, ref, beta, 2, loss, err, 256));
   }
   THD_picard_noise_free(noise);
   THD_picard_plan_free(plan);
   puts("PASS noise PCA, rank, batching, PC loss and selection contracts");
   return 0;
}
