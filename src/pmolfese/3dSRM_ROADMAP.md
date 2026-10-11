# 3dSRM roadmap

Where 3dSRM (Shared Response Model) came from, what it does now, how it was built and checked, and what comes next.
It merges the old `3dSHM_ideas.md` (the motivation) and `3dSRM_next_steps.md` (the handoff notes). Last updated 2026-10-10.

If you are a new instance picking this up: read all of it, then the memory files in
`/Users/molfesepj/.claude/projects/-Users-molfesepj-Documents-Programming-afni/memory/`. The repo is
`/Users/molfesepj/Documents/Programming/afni`, branch `SRMthings`. The user is Peter Molfese (AFNI developer).

---

## Part 1. Background: Shared Response Models

### The problem they solve

Suppose 20 people watch the same movie in the scanner. They all have a "face-responsive" region, but when you line the
brains up anatomically (MNI or Talairach), voxel #1234 in person A isn't doing the same thing as voxel #1234 in person B.
Anatomical alignment matches folds and landmarks, but **functional topography** (where fine-grained information lives)
varies a lot between people. Averaging across people therefore blurs the fine-grained patterns, and a decoder trained on
one person doesn't transfer to the next.

### The core idea

SRM assumes that everyone's brain is running **the same underlying response to the stimulus**, but each person's voxels
"read it out" differently:

```
X_i  ≈  W_i · S  +  noise
```

- `X_i` is subject *i*'s data (voxels × time points).
- `S` is the **shared response** (k × time points). It's a small number of latent time series (k is often 10–100) that
  every subject shares.
- `W_i` is a **subject-specific mapping** (voxels × k) with orthonormal columns. It rotates and selects from that
  person's voxel space to reproduce `S`.

The model learns `S` and every `W_i` at the same time, using the fact that the stimulus was time-locked across
subjects. Think of it as finding a common language, with `W_i` as each person's translator. It descends from Haxby's
**hyperalignment** (2011). The Chen et al. 2015 version adds a probabilistic formulation and dimensionality reduction
(k ≪ voxels). BrainIAK's `SRM`, `DetSRM` and `RSRM` classes implement it.

### Why it's useful in fMRI

- **Cross-subject decoding.** You can train a classifier in the shared space and test it on a new person. This usually
  beats anatomical alignment by a wide margin.
- **Better inter-subject correlation (ISC).** After alignment, ISC and inter-subject pattern correlation rise because
  you've removed the topographic mismatch.
- **Denoising.** Only the part of the signal that is shared across people survives the k-dimensional bottleneck, so
  idiosyncratic noise and artifacts get dropped.
- **Fine-grained comparison across people.** You can compare representations without needing the same voxels.
- **Built-in alignment.** A new subject only needs to watch some of the shared stimulus to estimate their own `W_i`. You
  can then map them into the shared space.

### Caveats

- **It needs a common temporal structure.** The movie, story or task must be time-locked across subjects. It doesn't
  make sense for plain resting-state data, because there's no shared time series to align on. Connectivity-based
  hyperalignment is the resting-state variant. (3dSRM's item axis relaxes "time-locked" to "the same items, in any
  order", for betas.)
- **You must cross-validate properly.** Fit `W_i` on some time segments and evaluate on held-out segments. Otherwise the
  model can fit noise and give you circular, inflated results.
- **k is a free parameter.** Too small throws away signal and too large lets noise back in.
- **It's linear and orthogonal.** It handles rotations of the representation, not arbitrary nonlinear differences.

### Could it help with EEG, or EEG-fMRI?

SRM only needs several "views" of the same time series, with different feature spaces per view. That doesn't have to
mean different voxels in different people, so conceptually it extends. The following is an assessment, not settled
practice.

**EEG alone, across subjects.** This is the most natural fit. EEG has the same problem: electrode patterns differ across
people because of head shape, cap placement and cortical folding. You'd use a stimulus-locked paradigm (movie, speech,
music). SRM would map each person's ~64 channels into a shared k-dimensional space, where k must be below the channel
count, so maybe 5–20. The established EEG analogue is **Correlated Component Analysis** (Parra, Dmochowski), which is
close in spirit.

**Separate EEG and fMRI sessions with the same stimulus.** This is plausible, with two hurdles:
1. The timescales differ. EEG is sampled at hundreds of Hz and fMRI at a TR of 1–2 s. You'd work with EEG band-power
   envelopes at TR resolution, and the HRF has to be handled somewhere (convolve the EEG features, or allow a lag).
2. The two modalities measure different physiology, so only the slow, stimulus-driven portion of the EEG power will
   share structure with the BOLD signal. Treat the EEG group and the fMRI group as different "subjects" feeding one
   shared response. Different participants in each modality is fine, since SRM doesn't need the same people.

**Simultaneous EEG-fMRI.** Here you have two views of the same brain at the same time, and you could treat EEG and fMRI
as views sharing a latent response. The neurovascular lag and the heavy artifacts (gradient, ballistocardiogram) make
this the hardest case. Established alternatives such as EEG-informed fMRI regressors, joint ICA, multiset CCA and PLS
are probably safer baselines to compare against. SRM could still help in a supporting role, for example denoising the
in-scanner EEG using the cross-subject shared component.

**Overall assessment.** Within a modality, SRM is well established for fMRI. For EEG it's reasonable but less mature.
For cross-modal EEG↔fMRI it's an open research idea, and a result there would need careful validation against existing
fusion methods. 3dSRM already takes `.1D` inputs (EEG band power at TR resolution, parcel means) alongside datasets in
one table, and its probabilistic model gives each subject (or modality) its own noise variance.

---

## Part 2. Working rules (from the user; do not break them)

- **Never run `git push`, `git add`, or any staging command.** Ask before committing. Hand over git commands for the
  user to run.
- **Cap threads at about 4 and `nice -n 15` heavy jobs.** Set `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`,
  `MKL_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS` and `NUMEXPR_NUM_THREADS` to 4. Run long jobs in the background.
- **AFNI log style.** The program name appears once, on the banner line (like 3dttest++). No per-line `3dSRM:` prefix.
  Use `INFO_message` / `ININFO_message` / `WARNING_message` / `ERROR_exit`. Progress bars: one per stage, labelled,
  with elapsed time and an estimate of time left. `check_3dSRM.py` asserts the log style.
- Use the scratchpad directory for temporary files, not `/tmp`.
- Write code in the style of the surrounding code: dense AFNI C, same comment density, same naming.
- The help text is the specification. Change it with the code, and keep lines within 80 columns.
- Check every feature against an independent computation (numpy, or the `srm_cli` driver), not against 3dSRM's own
  output. Report what passed, what is still refused, and the git commands.

---

## Part 3. Where things are

| File | What it is |
|---|---|
| `src/thd_srm.h`, `src/thd_srm.c` | The library: fit (prob/det; shared or run-coupled W; voxel or Gram kernels; missing data), project, update, canonicalize, metrics (isc, segment matching with masks, item identification, r2, references), model file I/O. The header documents the API. |
| `src/3dSRM.c` | The program: help, options, data stage (time and item axes), censoring, Xmat handling, CV, scoring (`-score_by`), outputs, `-model_in`, `-update_S`, `-apply_W`. |
| `src/thd_datatable.{c,h}` | The `-dataTable` reader, shared with 3dRSA, 3dInSync and 1dTrdm (compiled into each). |
| `src/thd_xmat.{c,h}` | Design matrices (X.xmat.1D): reading, every attribute as key/value, runs, censoring, per-TR conditions, items, and the table helpers (`THD_datatable_expand_xmat_items`, `THD_datatable_add_xmat_columns`). Unit test `thd_xmat:unit`. |
| `src/pmolfese/tests/test_srm.c` | Library unit tests (ctest `thd_srm:unit`). |
| `src/pmolfese/tests/srm_cli.c` | Library driver with named options (`--obs`, `--kernel`, `--init` ...): the independent reference in the end-to-end checks. |
| `tests/scripts/check_3dSRM.py` | End-to-end checks of 3dSRM against numpy and the driver (about 400 checks). Needs numpy and nibabel; uses the built 3dDeconvolve for the Xmat cases. |
| `tests/CMakeLists.txt` | ctests `3dSRM:check-auto`, `3dSRM:check-gram` (all checks, with the fit forced onto Gram matrices) and `thd_srm:unit`, and the `srm_cli` target. |
| `src/CMakeLists_binaries.txt`, `src/Makefile.INCLUDE` | Build wiring for 3dSRM. |
| `~/Documents/Programming/afni_srm_compare/compare_SRM_BrainIak.py` | The BrainIAK comparison and real-data studies, kept outside the repo. Set `AFNI_REPO` to point it at a checkout. |

## Part 4. Build and test

```bash
cd /Users/molfesepj/Documents/Programming/afni
nice -n 15 ninja -C build -j4 3dSRM srm_cli thd_srm_unit
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 VECLIB_MAXIMUM_THREADS=4 NUMEXPR_NUM_THREADS=4
nice -n 15 ctest --test-dir build -L 3dSRM                 # everything, about 40 s
# one area at a time (the build venv has numpy and nibabel):
nice -n 15 build/venv/bin/python3 tests/scripts/check_3dSRM.py --srm-cli build/targets_built/srm_cli \
     --workdir <scratchpad>/check --only-items      # or --only-stage3, --only-concat, --only-missing, --only-xmat
```

A change to `thd_datatable.c` also touches 3dRSA, 3dInSync and 1dTrdm: run
`ctest --test-dir build -R "RSA|InSync|Trdm|patterns|mapinfer|insync"` too.

---

## Part 5. What 3dSRM does now (all checked)

**Inputs.** A data table (`Subj`, `InputFile`, optional `Group`, `Run`, `Mask`, `Censor`, `Xmat`, covariates).
Datasets and one-row-per-time-point `.1D` files can be mixed. Per-subject masks and grids.

**Two axes.**
- *Time* (`-axis time`, default): time-locked series, runs concatenated. Runs come from a `Run` column, or from inside
  one input per subject (afni_proc.py's errts) via `-concat '1D: 0 300 600'` (as in 3dDeconvolve) or an `Xmat` column
  (its RunStart).
- *Item* (`-axis item`): one row per Subj × Item, each a single volume (`[n]` or `[label]` selector, each file read once
  per subject) or a one-row `.1D`. Items matched across subjects by key; repeats keyed `Item:1, Item:2` (run order, then
  table order) or by a `Repeat` column. Subjects' runs may differ (randomized designs). An all-zero item is an error.
  Time-only options are refused with what to do instead.
- *Items from afni_proc.py output*: one row per subject (`Subj`, `Xmat`, `InputFile` = stats bucket) expands to one item
  per stimulus column of the Xmat (`house`; `face#0, face#1, ...` for `-stim_times_IM` or TENT), reading the
  `LABEL#n_Coef` volumes, with `Condition` and (when each regressor lies in one run) `Run` columns. `-xmat_stims`
  picks conditions; `-show_table` prints the expanded table.

**Tables.** `#` starts a comment only at a word's start (`face#0`, `stats+tlrc[face#0_Coef]` work unquoted); a whole
value may be quoted (`"face#0"`); a quote inside a word is text (`x.1D'`). `-show_table` quotes what needs it.

**Normalization.** `-zscore run` (time default), `all` (item default: betas of one GLM share a baseline), `none`.
Features constant in a run (or overall, for `all`) are dropped.

**Censoring.** Global `-censor` (drops for everyone). Per subject: `Censor` column, Xmat `GoodList`, `-zcensor`
(all-zero volumes). `-cenmode interp` (default: linear within the run, as 3dTproject NTRP; keeps every subject's
projected series complete for 3dInSync/3dRSA), `drop` (for everyone; warned about), or `missing` (that subject only).

**Fit.** Probabilistic (EM) or deterministic SRM, `-nstart`, `-run_coupling` (a W per run tied to an anchor).
Iterates on per-subject Gram matrices of Z = [X | W_start] when that is cheaper (9× on 22 subjects of story data, same
fit to 2e-8; `AFNI_SRM_KERNEL=VOXEL|GRAM` forces one). Missing data (`SRM_params.obs`) in the shared-W fit: DET averages
the observed subjects at each t; PROB uses their precision only. Results do not depend on the thread count.

**Validation.** `-cv run` (leave a run, a block, or a `-fold_column` level out; each subject's own folds when items were
randomized, via missing data), `-cv loso` (leave a subject out), `-cv none`. Metrics: isc, segment matching (`match`) or
item identification (`ident`), r2, each with its chance level; `-save_r2` maps. `-score_by` scores each level of a
factor within each subject: table columns, `-time_labels`, or `Condition` from the Xmat (each TR's strongest stimulus
regressor at `-xmat_label_frac` of its peak). `-nfeatures` lists give a comparison table. `-cv_fit within` (a model
per group) or `reference:GROUP` (others scored as new subjects in one group's model), and `-balance B` (fits on
subsets the size of the smallest group, scores averaged). Reliability: `-reliability split` (stability of S across halves of each
group) and `within` (each subject's own split-half reliability, from repeats; `-repeated_stimulus` averages repeated
runs). `-nperm` tests group and condition differences in the scores by permutation (`PREFIX.tests.txt`).

**Models and new data.** `PREFIX.srm` (plain text; S, noise model, provenance, item keys). `-model_in` maps new
subjects into a saved model (with `-cv run`); `-update_S` (EXPERIMENTAL) lets S absorb them; `-apply_W` carries a
subject's W to other data. All three work on either axis, and `-apply_W` also goes from a time-axis model to item
betas (the movie-to-task design).

**Outputs.** `S.1D` (item keys in its header on the item axis), `proj.SUBJ` (labelled with item keys), `fit.txt`
(features, rho2, censored count), `objective.1D`, `rho2.1D`, `scores.txt` and `scores_folds.txt` (Factor/Level rows),
`coupling.txt`, `nfeat.txt`, and optional W, recon and r2 maps.

### Decisions made with the user (keep them)

- Item axis: subjects may have different run assignments; repeats are their own columns; `-zscore all` is the item
  default; `#` and quotes in tables as above; Xmat item keys keep AFNI's names.
- Censoring: per subject, `interp` the default (complete projected series); `missing` available; `drop` warned.
- `-cv run` on randomized items holds out each subject's own runs (missing data); `-cv loso` likewise.
- `-score_by` matching is among the level's own samples (windows wholly inside it, or its items), with its chance.
- `-update_S` is labelled EXPERIMENTAL.

### Pitfalls already paid for

- AFNI writes ny=nz=1 datasets as 1D text unless `AFNI_WRITE_1D_AS_PREFIX=YES` (set before writing k×1×1 projections).
- `mri_read_1D` is transposed: nx = rows (time), ny = columns (features).
- `THD_open_dataset` needs the full name for AFNI-format datasets: use `srm_dset_find`.
- Library reductions are ordered on purpose so results are bit-identical across threads. Keep it that way.
- Printed tables carry 7 significant digits: tolerances of about 1e-6 for fractions.
- Gauge: S is defined up to a rotation. Compare projections up to Procrustes; canonicalize before comparing S. Held-out
  scores are rotation-invariant; held-out projections are aligned to the full fit's gauge for output.
- `#` in a table is a comment only at a word's start; a trailing apostrophe is AFNI's transpose (`x.1D'`), not a quote.
- Repeats and runs are ordered by first appearance in the table; changing which run label appears first renumbers
  repeats. A `Repeat` column avoids that.
- Per-run offsets in betas need `-zscore run`; with `all` they swamp item identification (a check caught this).
- The TR in a NIfTI header matters to 3dDeconvolve: a wrong TR silently empties late regressors.
- numpy's `slogdet` on macOS Accelerate raises spurious floating-point warnings: take log-determinants from a Cholesky
  factor in references.
- Converged EM wobbles at about 1e-8 relative with float32 W (voxel kernel); the Gram kernel keeps W in doubles.
- BSD `sed -i` needs an argument on macOS. Use python for multi-line edits.

---

## Part 6. Roadmap

### Done 2026-10-10: a shared Xmat module

`src/thd_xmat.{c,h}` (unit test `thd_xmat:unit`, `src/pmolfese/tests/test_xmat.c`); 3dSRM uses it for runs and
censoring from one input, items from the stats bucket, per-TR conditions, and the `Xmat.*` covariate columns (on
whenever there is an Xmat column). Like `thd_datatable.c`, it is compiled into each program that uses it (not libmri);
the table-level helpers live in `thd_xmat.c` so programs without design matrices need not link it. The plan as built:

1. **`src/thd_xmat.{c,h}`**: read an X.xmat.1D once into a struct (matrix, column labels, TR, RunStart, NRowFull,
   GoodList, stimuli with label, columns, basis option and formula, GLTs, command line), with **every header attribute
   available as a key/value pair** (`THD_xmat_attr(xm,"CommandLine")`). Derived helpers: the run of a time point, the
   censored mask on the full timeline, per-TR condition labels (strongest regressor at a fraction of its peak), and the
   item list (key, condition, `_Coef` label, run when the regressor lies in one run). Unit tests against a real
   3dDeconvolve Xmat.
2. **`thd_datatable` aware of it** (opt-in, so existing programs see nothing new):
   `THD_datatable_expand_xmat_items()` (one row per subject → one per stimulus column, as 3dSRM does now) and
   `THD_datatable_add_xmat_columns()` (per-subject facts as ordinary columns: `Xmat.TR`, `Xmat.nrun`, `Xmat.ntime`,
   `Xmat.ncensored`, `Xmat.frac_censored`, `Xmat.nstim`; numeric, so they travel with scores as covariates).
3. **3dSRM switches to the shared code** with no change in behavior (the Xmat checks against real 3dDeconvolve output
   guard that).

Still open: further facts as columns (events per condition, censored fraction per run).

### Then: ask the user before doing any of these

**Other programs (only when asked):**
- *3dInSync*: per-subject censoring and run boundaries from one errts (`GoodList`, `RunStart`) instead of `-zcensor`'s
  union; run-aware temporal nulls (no circular shift across runs); condition ISC from one continuous file (per-TR
  condition labels make the Subj × Condition cells). Caveat to document: an errts from a model with task regressors has
  the task signal removed, so task ISC wants a nuisance-only model (the task Xmat can still give the labels).
- *3dRSA*: conditions chosen and named from `StimLabels` (no hand-written selectors, order guaranteed); per-run condition
  estimates for crossnobis from one afni_proc GLM with trial-wise betas (average each condition's trials within each
  run); exemplar (trial) RSA with the same-order caveat; possibly noise normalization from errts residuals (check how
  3dRSA estimates noise first).
- Whether 3dRSA and 3dInSync should move to the one-time program-name logging style.

**Back to 3dSRM: still refused ("not available in this build yet"), with the help already written as the spec:**
1. (Done 2026-10-10: `-cv_fit within|reference:GROUP` and `-balance B`; the written model stays the pooled fit,
   held-out projections are rotated into its gauge, `reference` needs `-cv loso`, draws in `PREFIX.balance.txt`,
   a subject in no subset gets NA with a warning.)
2. (Done 2026-10-10: `-reliability split` (`PREFIX.reliability.txt`, `_splits.txt`; halves warm-started from the full
   fit, mean cos^2 of the principal angles), `-reliability within` (odd vs even repeats or runs, Spearman-Brown, in
   scores.txt with rel_ratio = isc/rel_pattern, per level too), `-repeated_stimulus` (runs averaged for the fit; not yet
   with -score_by, -cenmode missing, -run_coupling), `-nperm` (group / level / interaction tests on the subject scores
   with thd_permute: Welch or one-sample t, exact when the relabelings are few, family-wise p over a test's rows).)
3. (Done 2026-10-10: `-cv_adapt FRAC` with `-run_coupling` under `-cv run` and `-cv loso`. Each held-out run's W
   is estimated from its first FRAC (`SRM_fit_view_coupled`, tied to the anchor with the fit's lambda) and the rest
   is scored. The main columns hold the run W's scores, and `isc_anchor`, `match_anchor` and `r2_anchor` hold the
   anchor's on the same samples. Under `-cv run`, S of the first samples is the mean of every subject's anchor
   projection, and every subject gets a run W. Under `-cv loso`, S is the group fit's. `PREFIX.proj` stays the
   anchor's. With `-cv_fit`/`-balance` both columns are averaged over the fits, and under `-cv_fit reference` the
   other groups adapt against the reference S, as in `-cv loso`. Checked against numpy on the driver's coupled
   fold fits for `-cv run` (and for each group under `-cv_fit within`). For `-cv loso` (its warm-started refit is not reproducible outside), only the limits are checked:
   lambda 1e9 makes the run W the anchor.)
4. (Done 2026-10-10: `-model_in` and `-apply_W` on the item axis. With `-model_in`, new subjects must have exactly
   the model's keys, in any table order; the columns are put in the model's order, and a mismatch names the first
   missing and extra key. `-cv run` uses item folds and may use each subject's own randomized runs, since S is
   fixed. `-update_S` writes an item model. With `-apply_W`, each subject is read from its own rows, so subjects may
   have different items (at least 3 each). A model z-scored within runs, applied to items without a Run column,
   z-scores over all of the subject's items and says so.)
5. (Done 2026-10-10: `-concat`/Xmat with `-apply_W` (per subject, so subjects' Xmats may differ). `-cenmode missing`
   with `-run_coupling`: the coupled fit now takes masks, with the same E-step and noise updates as the shared fit,
   and lambda 0 needs k observed points per run and view. With `-update_S`: the model carries `nsubj_t` and, for
   PROB, `rho_inv_t` (per-time sums, saved as data blocks), so an update is right whether the old fit, the new
   subjects, or both have missing points. Before this, a model fit with missing data and updated without them used
   the wrong denominators; now it cannot. `SRM_model_add_view_obs` places a subject from observed columns only. With
   `-apply_W`: z-scored over observed samples, 0 at missing ones. `-cv_adapt` with missing data uses observed samples
   only.)
6. (Done 2026-10-10: `-repeated_stimulus` with `-cenmode missing` (each time point averaged over the runs observed
   there; 'within' over time points both halves observed) and with `-score_by` (labels read before averaging, must
   agree across repeats; `-time_labels` may give one label per stimulus time point).)

Nothing is refused as "not available in this build yet" any more. The combinations that stay errors are errors by
design, each saying why: for example, `-run_coupling` on the item axis, `-repeated_stimulus` with `-run_coupling`,
`-cv loso` with `-model_in`, and `-reliability` with `-model_in`.

**Memory (found on the Finn paranoia data, 2026-10-10):** `-cv run` used to allocate each fold's training and test
copies anew, at a different size per fold. macOS does not give such freed blocks back promptly, so the peak climbed
by most of a copy of the data over the folds (48.3 GiB measured against a 35.8 GiB estimate, 22 subjects × 161k voxels
× 1310 TRs). Each subject's pair of fold buffers is now allocated once, at the largest fold's size, and reused
(38.0 GiB, scores bit-identical). The preflight now counts what is held: the data, plus the larger of the largest
fold's copy ((T − shortest run) + longest run, under `-cv run`; one subject's data for `-cv loso` and `-model_in`)
and one input on its full grid while it is read. Estimate 38.4 GiB, measured 38.0.

**`-compress` (2026-10-11, opt-in while verified):** each subject is read, censored and z-scored one at a time,
then replaced by R (rank × T, at most T × T) with R'R = X'X, a pivoted Cholesky of X'X formed in double
(`SRM_compress`). Everything the fit and CV do with X goes through X'X, and W stays in the span of X, so a fit of R
with `SRM_params.ndim` = voxels (the noise variance is per voxel) is the fit of X with W = U W_R. It starts where
the default would, `SRM_compress_start` = U'W0 carried in as `SRM_params.w0`, so S, the scores, projections, rho2
and the model agree with the default to about 1e-7, in every configuration tested (unit tests; check_3dSRM
`case_compress`; the paranoia data). Voxel maps (`-save_W`, `-save_recon`) read each subject again: W = X M,
M = `SRM_compress_map(R, W_R)` (T × k), with `SRM_compress_expand` threaded. To make fits of subsets agree too,
every fit now gives each subject the random-start stream of its whole-table index (`SRM_params.stream`), in both
modes. A subject therefore starts the same way in every fit; this changed the default's `-cv_fit`/`-balance`
starts. The remaining gaps are closed:
- `-save_r2`: each scored fold notes its W and means in R, its reference, and its held-out columns. At the re-read,
  W = X M is scored on the voxels with `SRM_reconstruction_r2`, averaged over folds and then over fits, as in the
  default.
- Item-axis `-reliability within`: made from each subject's voxels before they are compressed. The `-score_by`
  labels come from the table first.
- `-cenmode drop -zcensor`: a first read of every subject finds the all-zero volumes.

All agree with the default to about 1e-7. Not with `-model_in` or `-apply_W`, which hold only the new subjects. On the full Shen mask (22 × 161k voxels × 1310 TRs, `-cv run`, k = 50): 38.0 GiB / 6.6 min →
3.7 GiB / 3.9 min.
`-float16` was considered and dropped: the fit and CV hand the library every subject's data at once as floats, so
16-bit storage would need half-precision kernels throughout; with `-compress` it would save nothing.

**Performance, if wanted:** reuse each subject's Gram matrix across CV folds and `-nstart` starts (also removes the CV
copy of the data); tuned plain-C kernels for the passes that still read the voxels (forming G, final W, projections, r2);
form G in parallel within a subject; count the full-run read (all grid voxels) in the memory preflight, or read volume
by volume.

**Real-data validation (2026-10-10, Finn et al. paranoia story, 22 subjects, 3 runs of 526/452/332 TRs, fmriprep
MNI 2 mm, motion censoring from FD > 0.5 mm: 0-108 TRs per subject):** separate runs and the concatenated files with
`-concat` give bit-identical S and scores (full Shen mask, 161k voxels). On a 27-parcel temporal ROI (18.7k voxels),
every combination ran: censoring interp/missing/drop, `-run_coupling` + missing + `-cv_adapt`, `-cv_fit` within and
reference with `-balance`, `-score_by Run`, `-nperm`, `-reliability split`, a k scan, `-model_in`/`-update_S` with
missing data, and `-apply_W` with an Xmat and with missing data. Checks against numpy on the most-censored subject
(fit projections; `-apply_W -cenmode missing`) agree to 3e-8, and `-apply_W` on the concatenated file with an Xmat
reproduces the fit's projections exactly. Two memory problems found and fixed there (see Memory above, and: the
library no longer copies a view with missing points when its missing values are already 0 and it is not centered,
7.1 → 4.8 GiB on the ROI). Not exercised on real data: `-axis item` and `-repeated_stimulus` (no betas or repeats).

**Other open items:** the PR description for the branch (the user keeps `*_PR.md` drafts in the repo root).

### History

- Sep–Oct 2026: library and program (fit, CV, model file, `-model_in`, `-update_S`, `-apply_W`, run coupling).
- 2026-10-08: Gram-matrix fit (9× on real data); ctests; `-update_S` labelled experimental; BrainIAK comparison moved out.
- 2026-10-09: item axis (data stage, CV, `ident`); `-zscore all`; repeats; `-concat` and Xmat runs; per-subject censoring
  (`interp`, `drop`).
- 2026-10-10: missing data (library and program; `-cv run` on randomized items); items from Xmat + stats bucket;
  `-score_by` with Condition from the Xmat; quotes and word-start comments in data tables; this roadmap (merging the
  motivation notes and the handoff); the shared `thd_xmat` module; `-cv_fit`, `-balance`, `-reliability`,
  `-repeated_stimulus`, `-nperm`; item-axis `-model_in`/`-update_S`/`-apply_W`; `-cv_adapt`; the last refused
  combinations (missing data with coupling, update and apply_W; Xmat/-concat with apply_W; repeated stimulus with
  missing data and -score_by; -cv_adapt with -cv_fit/-balance).
