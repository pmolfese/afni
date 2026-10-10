#include "mrilib.h"
#include "thd_datatable.h"
#include "thd_mapinfer.h"
#include "thd_srm.h"

#ifdef USE_OMP
#include <omp.h>
#endif

#define PROGRAM_NAME "3dSRM"


static char *srm_options[] = {
   "-prefix", "-nfeatures", "-algorithm", "-niter", "-tol", "-seed", "-nstart",
   "-mask", "-run_column", "-group_column", "-zscore",
   "-censor", "-zcensor", "-zcensor_warn",
   "-axis", "-item_column", "-repeated_stimulus", "-time_labels",
   "-model_in", "-update_S", "-run_coupling", "-cv_adapt", "-apply_W", "-score_by", "-cv_fit", "-balance", "-reliability",
   "-nsplit", "-nperm", "-cv", "-cv_blocks", "-cv_metrics", "-match_window", "-loso_iter",
   "-save_W", "-save_recon", "-save_r2", "-proj_format", "-no_projected",
   "-memory_limit", "-memory_override", "-quiet", "-progress",
   "-dataTable", "-dataTableFile", "-show_table", "-help", "-h", NULL
} ;

static void usage_3dSRM(void)
{
   printf(
"\n"
"Usage: 3dSRM [options]                                         ~1~\n"
"\n"
"Fit a Shared Response Model (SRM; Chen et al., NIPS 2015) to time-locked\n"
"data from several subjects, and measure how well it generalizes.\n"
"\n"
"Each subject i contributes a matrix X_i (features x time). SRM learns one\n"
"shared response S (k x time) and, for each subject, a mapping W_i\n"
"(features x k, orthonormal columns) so that\n"
"\n"
"      X_i  ~=  W_i S  +  noise\n"
"\n"
"The subjects need not have the same number or arrangement of features.\n"
"Subject i's voxels (or EEG channels) are read out into the same k-dimensional\n"
"space, so no voxel-to-voxel correspondence across subjects is required.\n"
"Only the columns must correspond. The columns can be time points of a\n"
"stimulus every subject experienced at the same time (a movie, a story), or\n"
"items (trial or condition betas) matched across subjects by a label, which\n"
"is the right choice when trial order was randomized. See -axis. SRM is\n"
"meaningless for resting-state data.\n"
"\n"
"Typical uses:\n"
"  * Subject-wise k x time series in a common space (the 'projected' output)\n"
"    for 3dRSA, 3dInSync, decoding, or a group t-test of a held-out score.\n"
"  * Denoising: only activity shared across subjects survives the bottleneck.\n"
"  * Placing new subjects into an existing space, now or months later, from a\n"
"    saved model and without the original subjects' data (see -model_in).\n"
"\n"
"Required input:                                             ~1~\n"
"\n"
"  -dataTable @TABLE\n"
"  -dataTableFile TABLE\n"
"  -dataTable Subj Group InputFile s01 control s01_movie+tlrc ...\n"
"\n"
"      TABLE needs Subj and InputFile columns, plus optional columns:\n"
"        Run    : with one row per Subj x Run, each subject's runs are\n"
"                 concatenated in time in the order the runs first appear.\n"
"                 Every subject must have the same set of runs, each of the\n"
"                 same length. See -run_column.\n"
"        Item   : item label, with -axis item. See -item_column.\n"
"        Xmat   : the subject's X.xmat.1D (from afni_proc.py or 3dDeconvolve),\n"
"                 when one InputFile holds all of the subject's runs, such as\n"
"                 afni_proc.py's errts. Its RunStart splits the input into\n"
"                 runs, and its NRowFull must equal the input's length. See\n"
"                 -concat.\n"
"        Mask   : a mask dataset for that subject, on that subject's own\n"
"                 grid. Subjects may then have different voxel counts and\n"
"                 different grids. A subject's Mask must be the same on all\n"
"                 their rows.\n"
"        Group  : subject group, used to label and compare scores. See\n"
"                 -group_column.\n"
"        others : any other column can be named in -score_by. Numeric\n"
"                 columns are carried to the output tables, so covariates\n"
"                 such as motion travel with the scores.\n"
"\n"
"      InputFile is an AFNI dataset, or a .1D file. A .1D file is read as\n"
"      time x features (one row per time point, one column per channel or\n"
"      region), the usual AFNI orientation. Datasets and .1D files can be\n"
"      mixed in one table, which is how non-voxel data (EEG band power at TR\n"
"      resolution, parcel means) enter the same model. Every input must have\n"
"      the same number of time points.\n"
"\n"
"  -nfeatures K[,K2,...]\n"
"                       Dimension k of the shared space. Required, except with\n"
"                       -model_in, which takes it from the model. K must not\n"
"                       exceed the smallest feature count or the number of\n"
"                       time points. Several values require -cv and produce a\n"
"                       comparison table (PREFIX.nfeat.txt) and no model\n"
"                       outputs; rerun with the chosen K for those.\n"
"\n"
"Model options:                                              ~1~\n"
"\n"
"  -prefix PREFIX       Output prefix [SRM].\n"
"\n"
"  -algorithm ALG       prob : probabilistic SRM, fitted by EM [default].\n"
"                              Estimates a noise variance for every subject,\n"
"                              so noisier subjects (or noisier modalities)\n"
"                              count for less in S.\n"
"                       det  : deterministic SRM, alternating Procrustes\n"
"                              steps. Subjects are weighted by their feature\n"
"                              count and variance only.\n"
"\n"
"  -niter N             Maximum iterations [500].\n"
"  -tol T               Stop when the relative change in the objective falls\n"
"                       below T [1e-7].\n"
"\n"
"                       SRM converges slowly on real data, which is why the\n"
"                       defaults are generous. For 22 subjects hearing a story\n"
"                       (1000 voxels each, 15 components), the shared response\n"
"                       after 10 iterations, BrainIAK's default, differed from\n"
"                       the converged one by 38%%, after 20 by 30%%, after 100 by\n"
"                       15%%, and after 200 by 4%%. Full convergence (1e-13) took\n"
"                       about 550 iterations. -tol 1e-7 stops near 170 and\n"
"                       -tol 1e-8 near 250. A fit stopped earlier is not\n"
"                       wrong, only unconverged, and it depends on the seed.\n"
"                       Whether a converged or an earlier fit predicts held-out\n"
"                       data better is an empirical question: -cv answers it.\n"
"                       For a quick look, -niter 20 is enough.\n"
"\n"
"  -seed S              Seed for the random starting W_i [1234567]. Until a fit\n"
"                       has converged, its subspaces depend on the seed.\n"
"                       Compare fits with rotation-invariant quantities, not\n"
"                       component by component.\n"
"  -nstart N            Fit from N different random starts and keep the one\n"
"                       with the best objective [1]. Starts run one after\n"
"                       another, so the cost is N times as much. The\n"
"                       objective of every start is in PREFIX.fit.txt: if they\n"
"                       agree to the tolerance you ran at, they found the same\n"
"                       optimum. If they differ by more, run longer first.\n"
"                       Only then is a difference evidence of a local optimum.\n"
"                       In tests, both models reached the same optimum from\n"
"                       every start on real data. One hard synthetic problem\n"
"                       (noise differing 6-fold between subjects) stalled the\n"
"                       deterministic model at a slightly worse stationary\n"
"                       point from one start, so N of 3 to 5 is cheap\n"
"                       insurance when the data are that uneven.\n"
"\n"
"  -run_coupling L      shared : one W per subject for all of their runs [default]\n"
"                       L >= 0 : a separate W for each subject and run, each\n"
"                                tied to a per-subject anchor W. The tie is a\n"
"                                prior worth L time points of data: L=0 fits\n"
"                                every run on its own, L equal to a run's\n"
"                                length weighs the anchor as much as that\n"
"                                run's data, and L much larger than the runs\n"
"                                is nearly 'shared'. Try values from about\n"
"                                0.2 to 2 times a typical run length. Requires\n"
"                                runs, and -axis time. The fit starts from the shared-W fit\n"
"                                (so it takes about twice as long), and a very\n"
"                                large L stays close to it: use 'shared'\n"
"                                instead.\n"
"\n"
"                                Runs not in the fit (held-out runs, new data)\n"
"                                use the anchor. See -cv_adapt. In-sample, a\n"
"                                smaller L always fits better, because a W per\n"
"                                run has more freedom; whether it generalizes\n"
"                                is what -cv shows. With L=0 the runs can fall\n"
"                                into different orientations of S.\n"
"\n"
"Data options:                                               ~1~\n"
"\n"
"  -mask MASK           A mask for every subject lacking a Mask column. All\n"
"                       such subjects then must share MASK's grid. Datasets\n"
"                       need a mask from one source or the other; .1D inputs\n"
"                       never do. Voxels with zero variance in any run (over\n"
"                       all samples, with -zscore all) are dropped for that\n"
"                       subject and counted in the log.\n"
"\n"
"  -run_column NAME     Run column in the table [Run]. If absent, each subject\n"
"                       is a single run, unless -concat or an Xmat column\n"
"                       says where runs start inside each input.\n"
"  -concat SPEC         Run starts within each subject's one input, as in\n"
"                       3dDeconvolve and 3dTproject: the time-point index at\n"
"                       which each run begins, starting with 0, as a 1D file or\n"
"                       '1D: 0 300 600'. The same for every subject. For one\n"
"                       input per subject that holds several runs (afni_proc.py's\n"
"                       errts.SUBJ+tlrc, say); without it such an input is one\n"
"                       long run, so z-scoring and -cv would ignore the\n"
"                       boundaries. An Xmat column does the same per subject,\n"
"                       from RunStart (every subject must still have the same\n"
"                       runs). Not with a Run column.\n"
"  -group_column NAME   Group column in the table [Group]. Only used to label\n"
"                       the output tables.\n"
"\n"
"  -zscore MODE         run  : z-score every feature separately within each run\n"
"                              before fitting [default with -axis time]\n"
"                       all  : z-score every feature over all of the subject's\n"
"                              samples together [default with -axis item:\n"
"                              betas from one GLM share a baseline]\n"
"                       none : use the data as given (it should already be\n"
"                              centered and scaled)\n"
"\n"
"  -censor FILE         One common 0/1 vector with one value per input time\n"
"                       point. Retained time points are identical for all\n"
"                       subjects. Applied before z-scoring. -axis time only.\n"
"  -zcensor             Omit a time point when any subject's dataset volume is\n"
"                       all zero in that subject's mask. Common to everyone.\n"
"                       -axis time only (an all-zero item is always an error).\n"
"  -zcensor_warn PCT    Warn for each subject with more than PCT percent\n"
"                       all-zero volumes [20].\n"
"\n"
"  3dSRM does not regress nuisance signals, filter, or change temporal\n"
"  alignment. Inputs should already be preprocessed and aligned to the\n"
"  stimulus. For EEG-derived inputs, the HRF and any lag must be handled\n"
"  before the data reach this program.\n"
"\n"
"Time-locked or item-indexed:                                ~1~\n"
"\n"
"  -axis time           Columns are the time points of a stimulus presented\n"
"                       identically to everyone [default]. With one row per\n"
"                       Subj x Run, runs are concatenated.\n"
"  -axis item           Columns are items. The table has one row per Subj x\n"
"                       Item, whose InputFile is a single volume (or a\n"
"                       sub-brick selector such as betas+tlrc'[7]' or\n"
"                       betas+tlrc'[face#0_Coef]'), or a .1D file with one\n"
"                       row. Every subject must have every item exactly once,\n"
"                       in any row order; items are matched by label, not\n"
"                       position. Missing, repeated, or extra items are errors\n"
"                       that name the Subj and Item. A dataset that holds many\n"
"                       items is read once per subject, however many rows\n"
"                       select from it. An item that is zero everywhere in the\n"
"                       mask (usually a failed or missing estimate) is an error.\n"
"\n"
"  -item_column NAME    Item-label column [Item]. Labels may not contain ':'.\n"
"  -repeat_column NAME  Repeat-label column [Repeat], used if present: it names\n"
"                       an item's repeats (1, 2, ... or early, late).\n"
"\n"
"  Normalization. Betas from one GLM over all runs (the usual afni_proc.py\n"
"  analysis, including trial-wise betas from -stim_times_IM) share one\n"
"  baseline, so the default, -zscore all, z-scores each feature over all of\n"
"  the subject's items together. With betas from a separate GLM per run, use\n"
"  -zscore run, which works over the items of each of the subject's own runs\n"
"  (and needs a Run column). Each run then needs at least 3 items, and\n"
"  z-scores over fewer than about 10 are noisy (a warning says so).\n"
"\n"
"  Runs. A Run column orders repeats and can define the -cv folds. Subjects\n"
"  need not share runs: when trial order was randomized, an item may be in run\n"
"  2 for one subject and run 3 for another (the folds then come from\n"
"  -fold_column).\n"
"\n"
"  Repeats. An item measured more than once (separate betas: trial-wise, or\n"
"  one per run) gives a column per repeat, keyed Item:1, Item:2, ... in run\n"
"  order (runs ordered as they first appear in the table, as with -axis time)\n"
"  and then table order, so two repeats in the same run are numbered by their\n"
"  rows. A Repeat column names them instead (Item:R; see -repeat_column),\n"
"  which does not depend on row order. The model is free to use the repeats,\n"
"  and they are what within-subject reliability is computed from (see\n"
"  -reliability within). Every subject must have the same repeats of each\n"
"  item.\n"
"\n"
"  What changes with -axis item: runs are not concatenated in time, there is\n"
"  no temporal adjacency, so 'match' becomes item identification (the held-out\n"
"  item must be the best correlate of its own counterpart among all items),\n"
"  and the folds of -cv are runs (when every subject has the same items in each\n"
"  run), or the levels of -fold_column. Options about time are errors with\n"
"  -axis item, each saying what to do instead: -censor, -zcensor, -cv_blocks,\n"
"  -match_window, -time_labels, -repeated_stimulus, and -run_coupling. The\n"
"  outputs list the item keys: in the header of PREFIX.S.1D (one row per item)\n"
"  and as the volume labels of PREFIX.proj.SUBJ.\n"
"\n"
"  -fold_column NAME    The table column whose levels are the folds of -cv run\n"
"                       [Run]. Naming a category column gives\n"
"                       leave-one-category-out: the model must then generalize\n"
"                       to a kind of item it never saw, a stricter test.\n"
"\n"
"  -time_labels NAME FILE\n"
"                       With -axis time, one numeric label per retained time\n"
"                       point (0 = unlabelled), to score by stimulus content\n"
"                       such as speech/no speech or scene type. NAME is the\n"
"                       factor name for -score_by. May be given repeatedly.\n"
"\n"
"  Which one? If every subject saw the same stimuli in the same order, time\n"
"  mode uses the full time series. If trial order differed, only item mode is\n"
"  valid. Time mode cannot align subjects whose events are not simultaneous.\n"
"\n""Mapping new subjects into a saved model:                   ~1~\n"
"\n"
"  Every fit writes PREFIX.srm. It holds the shared response S, the noise\n"
"  variances, and exactly how the data were put on the model's timeline (the\n"
"  censored time points, the run lengths, and the z-scoring). It does not hold\n"
"  any subject's W_i or data, so it is small and can be passed between sites.\n"
"  It is plain text (about 1.2 MB for 30 components and 1310 time points),\n"
"  written to a temporary name and renamed, so an interrupted run never leaves\n"
"  half a model. The numbers are stored with enough digits to be read back\n"
"  exactly, and a model that is read and then updated gives exactly the result\n"
"  that updating the original would have. The format is described in thd_srm.h;\n"
"  keys and sections this version does not know are skipped, and a file from a\n"
"  newer version is refused with a message, not misread. A damaged or\n"
"  truncated file is reported as such.\n"
"\n"
"  -model_in MODEL.srm  Do not fit S. Instead estimate a W for every subject in\n"
"                       the table, against the S in MODEL.srm. Each subject\n"
"                       needs the same stimulus and the same number of time\n"
"                       points and runs the model was built from; this is\n"
"                       checked first, and the data are censored and z-scored\n"
"                       as the model's were. Their feature counts, masks, and\n"
"                       grids can be anything with at least k features.\n"
"\n"
"                       The model defines -nfeatures, -algorithm, -seed,\n"
"                       -run_coupling, -censor, -zcensor, and -zscore; giving\n"
"                       any of them with -model_in is an error, and so is\n"
"                       -niter or -tol unless -update_S is given (they then\n"
"                       control that update). -cv loso is also an error,\n"
"                       since there is no group fit to leave anyone out of.\n"
"                       -cv none or -cv run are allowed, and mean something\n"
"                       different here: see below.\n"
"\n"
"                       A subject who is already in the model is a warning\n"
"                       (their scores are optimistic: S contains them), and an\n"
"                       error with -update_S. A model made with -run_coupling\n"
"                       places new subjects with one W for all their runs,\n"
"                       unless -update_S, which gives them a W per run.\n"
"\n"
"  A new subject's W is fitted to their own data, so its fit to that data is\n"
"  flattering by construction. Use -cv run with -model_in: W is then estimated\n"
"  from all of the subject's runs but one (or blocks but one), and the subject\n"
"  is scored on the part left out. The scores are then directly\n"
"  comparable to those from -cv loso on the original subjects.\n"
"\n"
"  -update_S            ** EXPERIMENTAL **\n"
"                       Let S absorb the new subjects, and write a new\n"
"                       PREFIX.srm that includes them. Without this option the\n"
"                       saved model is untouched and S does not move.\n"
"                       It is an approximation to a refit (see below), and its\n"
"                       behavior and the model file's sums may change.\n"
"\n"
"                       The model keeps running sums of its subjects' projected\n"
"                       data (weighted by their noise variances), so S can be\n"
"                       updated from the new subjects alone, with no access to\n"
"                       the original data. This is the point of the model file:\n"
"                       a site can add its subjects to a shared model without\n"
"                       receiving or sending any subject's data.\n"
"\n"
"                       What it is, and is not. The old subjects' W stay as\n"
"                       they were fitted; only S and the new subjects' W move.\n"
"                       That is an approximation to refitting everyone. The\n"
"                       subjects added together in one run are treated jointly,\n"
"                       so their order within the table cannot matter; subjects\n"
"                       added in different runs can matter, slightly (S differed\n"
"                       by about 3%% between arrival orders in tests).\n"
"\n"
"                       What it costs, measured on 22 subjects hearing a story\n"
"                       (30 components, fit on two runs, scored on a third by\n"
"                       the correlation of each subject's projected data with the\n"
"                       shared response from everyone else). Adding 1, 3, and 6\n"
"                       subjects to 21, 19, and 16 gave a drift of 0.003, 0.015,\n"
"                       and 0.038. Against refitting everyone from scratch, the\n"
"                       EXISTING subjects' scores were lower by 1.6%%, 1.2%%, and\n"
"                       5.0%% (their W had not been refitted); the added subjects'\n"
"                       scores were no lower. Updating S brought it 16-25%% closer\n"
"                       to the refit S than leaving it frozen, but did not raise\n"
"                       held-out scores above those of simply mapping the new\n"
"                       subjects into the frozen S. So the update keeps the\n"
"                       model's S current and the model growing; it is not a\n"
"                       substitute for a refit. The loss was 1-2%% up to adding\n"
"                       about a sixth of the original subjects and 5%% at a third,\n"
"                       so refit from the data when drift reaches a few percent.\n"
"                       These are one data set and one split: differences of 0.001\n"
"                       score of 0.13 are within noise.\n"
"\n"
"                       The gauge of S is never changed by an update, so the\n"
"                       components and the projected series of earlier\n"
"                       analyses stay valid. Scores for the new subjects are\n"
"                       always computed BEFORE the update, from the model they\n"
"                       were not part of (-cv with -model_in).\n"
"\n"
"                       The model file contains sums over subjects of\n"
"                       projected data. With only a few subjects those sums are\n"
"                       not anonymous, so treat a model built from fewer than\n"
"                       about five as if it contained their data.\n"
"\n"
"                       An updated model is written to PREFIX.srm; it will not\n"
"                       overwrite the -model_in file unless it is the same name\n"
"                       and -overwrite is given.\n"
"\n"
"  Without -update_S, adding subjects never changes old results.\n"
"\n"
"Applying a subject's W to other data from the same subject: ~1~\n"
"\n"
"  The cleanest way to analyze a task is to align on data that are independent\n"
"  of what you test: fit on a movie (or on some runs), then carry each subject\n"
"  over to the task data with that subject's own W. Nothing about the task\n"
"  labels then touches the alignment.\n"
"\n"
"  -apply_W PREFIX      Use the W datasets from an earlier -save_W run (the\n"
"                       PREFIX.W.SUBJ datasets and PREFIX.srm). For each table\n"
"                       subject, find their W, take the new InputFile data on\n"
"                       W's own grid and voxels (no mask needed or accepted),\n"
"                       z-score as the model did, and write the data in the\n"
"                       shared space (PREFIX.proj.SUBJ). The new data may have\n"
"                       any length, order, and -axis; nothing is fitted. A\n"
"                       subject in the table but not in the model is an error.\n"
"                       Subjects may differ in length from one another. A\n"
"                       feature with no variance in a run of the new data is\n"
"                       zeroed. With -zscore none in the earlier run, each\n"
"                       feature is centered on the new data, since the means\n"
"                       of the old data are not saved. -mask, -cv and -save_W\n"
"                       are errors here.\n"
"\n"
"                       The projected series of different subjects can then go\n"
"                       to 3dRSA, 3dInSync, or a classifier. Since W is\n"
"                       fitted on other data, an effect found there cannot come\n"
"                       from the alignment.\n"
"\n"
"                       Requires that the earlier run used -save_W.\n"
"\n""Validation:                                                 ~1~\n"
"\n"
"  An SRM fits every subject's data with nfeatures x voxels free parameters.\n"
"  Anything measured on the data it was fit to is optimistic. These options\n"
"  evaluate on data the fit never saw.\n"
"\n"
"  -cv MODE             none : fit once to everything. Reported scores are\n"
"                              in-sample and are labelled as such [default].\n"
"                       run  : leave one run out. For each run, fit on the\n"
"                              other runs, then score every subject on the\n"
"                              held-out run (or -cv_blocks block). With\n"
"                              -model_in, S is fixed, and only each subject's\n"
"                              own W is refitted without the held-out run.\n"
"                       loso : leave one subject out. For each subject, fit S\n"
"                              on all OTHER subjects and every time point.\n"
"                              Then estimate that subject's W from their own\n"
"                              training runs and score them on a held-out run.\n"
"                              The held-out data touch neither S nor any W, so\n"
"                              the score is a clean estimate of how well a new\n"
"                              person maps into an existing space. Folds\n"
"                              within the subject are leave-one-run-out, or\n"
"                              -cv_blocks contiguous blocks if there is one run.\n"
"                              Costs one group fit per subject.\n"
"\n"
"  -cv_blocks B         Contiguous time blocks used as folds in place of runs\n"
"                       when the data have a single run [2]. Time points near a\n"
"                       block edge are autocorrelated with the other side, so\n"
"                       scores are a little optimistic; separate runs are\n"
"                       better.\n"
"  -loso_iter M         Iterations for each leave-one-subject-out group refit.\n"
"                       The refit starts from the full fit minus that subject,\n"
"                       so a few iterations are enough [5].\n"
"\n"
"  -cv_metrics LIST     Comma list from: isc, match, r2 [all three].\n"
"                       isc   : correlation of the subject's projected held-out\n"
"                               data with the shared response estimated from the\n"
"                               other subjects, per component, and their mean\n"
"                               in Fisher-z space (returned as r).\n"
"                       match : time-segment matching accuracy (Chen et al.),\n"
"                               as in BrainIAK's SRM example. Each window of\n"
"                               -match_window TRs of the subject's projected\n"
"                               data is matched to the best-correlated window\n"
"                               of the others' shared response, among its own\n"
"                               counterpart and every window that does not\n"
"                               overlap it. With -axis item it becomes item\n"
"                               identification: each held-out item must be the\n"
"                               best correlate of its own counterpart among\n"
"                               all items. Chance is reported alongside.\n"
"                       r2    : fraction of the subject's held-out variance\n"
"                               reconstructed from W_i and that shared\n"
"                               response. On voxel data this is usually tiny\n"
"                               or slightly negative, because variance that is\n"
"                               specific to a subject dominates (about -0.001\n"
"                               for 5000 voxels of story data, where segment\n"
"                               matching reached 12%% against 0.3%% chance).\n"
"                               Use it to compare methods or settings on the\n"
"                               same data, not as a measure of fit; isc and\n"
"                               match are the informative scores.\n"
"  -match_window N      Segment length in TRs for 'match' [9].\n"
"\n"
"  -cv_adapt FRAC       With -run_coupling L and -cv: the held-out run has no W\n"
"                       of its own, so it uses the anchor. FRAC in (0,1) instead\n"
"                       estimates that run's W from its first FRAC of samples,\n"
"                       tied to the anchor with the same L, and scores only the\n"
"                       rest. The anchor-only score on those same samples is\n"
"                       reported beside it (r2_anchor, r2_adapt, and so on), so\n"
"                       the comparison shows whether a run-specific W helps.\n"
"                       The samples used to estimate W are never scored [0].\n"
"\n"
"Reliability, and comparing groups or conditions:            ~1~\n"
"\n"
"  A low score for a group or a condition has three possible causes: its\n"
"  signal is unreliable, its signal is reliable but not shared with the others,\n"
"  or the model was built mostly from someone else. These options separate\n"
"  them.\n"
"\n"
"  -score_by NAMES      Comma list of table columns (or -time_labels names)\n"
"                       whose levels define subsets of the columns. Every\n"
"                       -cv metric is computed separately for each level, for\n"
"                       every subject, so you get faces vs houses within each\n"
"                       subject and then within each group. The model is\n"
"                       still fitted to all samples.\n"
"\n"
"  -cv_fit MODE         pooled    : one model from everyone [default]\n"
"                       within    : one model per Group, fit only to that\n"
"                                   group. Asks: does this group have a shared\n"
"                                   response at all?\n"
"                       reference:LABEL\n"
"                                 : fit on group LABEL only, then place the\n"
"                                   other groups' subjects in that model (with their\n"
"                                   W estimated from their own training folds).\n"
"                                   Asks: do they use the reference group's\n"
"                                   space? Reference subjects are scored\n"
"                                   leave-one-subject-out, the others as new\n"
"                                   subjects, so the scores are comparable.\n"
"  -balance B           Groups of unequal size bias scores: a model built from\n"
"                       40 people maps its members better than one built from\n"
"                       12. With -balance, every group fit is repeated on B\n"
"                       random subsets the size of the smallest group, and the\n"
"                       scores averaged [off].\n"
"\n"
"  -reliability LIST    Comma list from: within, split.\n"
"                       within : a subject's own reliability, computed from\n"
"                                repeats of the same stimulus (the same Item in\n"
"                                different runs; or with -axis time and\n"
"                                -repeated_stimulus, the runs). It is the\n"
"                                correlation of the two halves' patterns\n"
"                                (r_pattern across voxels, r_series across\n"
"                                time or items), with the Spearman-Brown\n"
"                                correction, in the subject's own voxel space\n"
"                                with no model involved. It bounds what the\n"
"                                shared-space score can reach.\n"
"                       split  : group-level. Randomly split each group's\n"
"                                subjects in two, fit S to each half, and\n"
"                                compare them with a rotation-invariant\n"
"                                subspace similarity. A group with no stable\n"
"                                shared response gives halves that disagree.\n"
"  -nsplit N            Number of random splits for 'split' [100]. Each is two\n"
"                       group fits, warm-started from the full fit.\n"
"  -repeated_stimulus   With -axis time: every run shows the same stimulus. The\n"
"                       runs are averaged for the fit, and their agreement\n"
"                       gives 'within'.\n"
"\n"
"  How to read them together (all subjects scored on held-out data):\n"
"    low within, low shared  : the group or condition is simply noisy.\n"
"    high within, low shared : each person is reliable but they do not\n"
"                              share a response, or the model does not\n"
"                              speak their language. Compare -cv_fit within\n"
"                              against reference.\n"
"    high within, high shared: SRM works for them.\n"
"  'Within' needs real repeats. Without them it is not estimated, and a low\n"
"  score cannot be attributed.\n"
"\n"
"  -nperm P             Test group and condition differences in the scores by\n"
"                       permutation, written to PREFIX.tests.txt [0 = none]:\n"
"                         two groups : group labels permuted over subjects\n"
"                         two levels of a -score_by factor : level labels\n"
"                                      swapped within each subject\n"
"                         both : the group x level interaction, by permuting\n"
"                                the within-subject differences between groups\n"
"                       P must be at least 20. -seed applies. Subjects are the\n"
"                       exchangeable units; the folds of one subject are never\n"
"                       split. For anything more elaborate (covariates, three\n"
"                       groups) take PREFIX.scores.txt to R.\n"
"\n""Outputs:                                                    ~1~\n"
"\n"
"  PREFIX.srm           The saved model (see above). With -model_in no new model\n"
"                       is written.\n"
"  PREFIX.S.1D          Shared response, time x k (one column per component).\n"
"                       Components are ordered by variance and sign-fixed, so\n"
"                       repeated fits are comparable up to noise.\n"
"  PREFIX.rho2.1D       Noise variance per subject (-algorithm prob).\n"
"  PREFIX.objective.1D  Objective per iteration. For prob it is the marginal\n"
"                       log-likelihood and rises; for det it is the summed\n"
"                       squared error and falls.\n"
"  PREFIX.fit.txt       One row per subject: features used, features dropped,\n"
"                       rho2, and the log of any censoring. With -update_S it\n"
"                       also gives the model's drift and number of updates.\n"
"  PREFIX.coupling.txt  With -run_coupling L: per subject and run, how far that\n"
"                       run's W is from the anchor (0 = identical, 2 =\n"
"                       orthogonal, as |W_ir - W_i|^2 / nfeatures) and the\n"
"                       gain in reconstruction R^2 over the anchor in the fit\n"
"                       data. Large deviations show where the subject's\n"
"                       mapping really moved.\n"
"\n"
"  Unless -no_projected:\n"
"  PREFIX.proj.SUBJ     Each subject's data in the shared space, W_i' X_i.\n"
"                       With -cv run or loso these are the held-out\n"
"                       projections concatenated across folds; with -cv none\n"
"                       they are in-sample.\n"
"  -proj_format FMT     dset : a dataset with k x 1 x 1 voxels and one volume\n"
"                              per time point, directly usable as input by\n"
"                              3dInSync (per-component ISC) and 3dRSA [default]\n"
"                       1D   : time x k text file\n"
"\n"
"  With -cv run or loso:\n"
"  PREFIX.scores.txt    One row per Subj x level of each -score_by factor (and\n"
"                       'all'), averaged over folds. Columns: Subj, Group,\n"
"                       Factor, Level, nsamp, isc, match (or ident), r2, and\n"
"                       with -reliability within: rel_pattern, rel_series, and\n"
"                       the ratio of the shared score to rel_pattern. Numeric\n"
"                       table columns are carried along. One score per subject\n"
"                       per condition is the right input for a group\n"
"                       comparison.\n"
"  PREFIX.scores_folds.txt\n"
"                       The same, one row per Subj x fold.\n"
"  PREFIX.reliability.txt\n"
"                       With -reliability split: per Group, the subspace\n"
"                       similarity of the two halves over the splits (mean,\n"
"                       and 2.5 and 97.5 percentiles).\n"
"  PREFIX.tests.txt     With -nperm: the tests above.\n"
"\n"
"  Optional, since they are large:\n"
"  -save_W              PREFIX.W.SUBJ: k volumes per subject, on that subject's\n"
"                       grid. These are the spatial maps of the components.\n"
"                       A .1D subject gets one row per column of the input;\n"
"                       the rows of dropped (constant) features are zero.\n"
"  -save_recon          PREFIX.recon.SUBJ: W_i S, the denoised reconstruction of\n"
"                       the subject's time series.\n"
"  -save_r2             PREFIX.r2.SUBJ: voxelwise reconstruction R^2 for each\n"
"                       subject (held-out under -cv). Maps on a common grid\n"
"                       can go to 3dttest++.\n"
"\n"
"Resources:                                                  ~1~\n"
"\n"
"  -memory_limit G      Refuse an estimated peak above G GiB. Without this\n"
"                       option the limit is 80%% of detected physical RAM.\n"
"  -memory_override     Continue after the memory estimate exceeds its limit.\n"
"  -progress M          auto, bar, line, or off [auto].\n"
"  -quiet               Suppress progress and informational diagnostics.\n"
"  -show_table          Print the parsed table and exit.\n"
"\n"
"  Memory is roughly the sum over subjects of voxels x time x 4 bytes for the\n"
"  data, plus voxels x k x 4 bytes per subject for W. Most of the time goes to\n"
"  matrix products of that size, run subject-parallel under OpenMP.\n"
"\n"
"  When every subject has at least T+k voxels (T time points), the fit forms\n"
"  each subject's (T+k) x (T+k) Gram matrix once and iterates on that, so an\n"
"  iteration costs about (T+k) x T x k per subject instead of voxels x T x k,\n"
"  and W is formed from the data at the end. This adds (T+k)^2 x 8 bytes per\n"
"  subject. The fit is the same either way, up to rounding. Setting\n"
"  AFNI_SRM_KERNEL to VOXEL or GRAM forces one path (for testing).\n"
"\n"
"Examples:                                                   ~1~\n"
"\n"
"  # Fit, with held-out run scores\n"
"  3dSRM -prefix movie_srm -nfeatures 50 -cv run \134\n"
"        -mask brain+tlrc -dataTableFile subjects.txt\n"
"\n"
"  # Subjects in native space, each with their own mask\n"
"  3dSRM -prefix movie_srm -nfeatures 30 -cv loso -dataTableFile native.txt\n"
"\n"
"  # Months later, place new subjects in that model and score them\n"
"  3dSRM -prefix new_subj -model_in movie_srm.srm -cv run \134\n"
"        -mask brain+tlrc -dataTableFile new_subjects.txt\n"
"\n"
"  # Item-indexed betas, scored by condition, controls as the reference\n"
"  3dSRM -prefix faces -axis item -nfeatures 20 -cv run \134\n"
"        -score_by Category -cv_fit reference:control -balance 50 \134\n"
"        -reliability within,split -nperm 5000 -dataTableFile betas.txt\n"
"\n"
"  # Align on the movie, then carry each subject over to the task\n"
"  3dSRM -prefix movie -nfeatures 50 -save_W -mask brain+tlrc \134\n"
"        -dataTableFile movie.txt\n"
"  3dSRM -prefix task -axis item -apply_W movie -dataTableFile task.txt\n"
"\n"
"  # Add this site's subjects to the shared model; the data stay local\n"
"  3dSRM -prefix movie_v2 -model_in movie_srm.srm -update_S -cv run \134\n"
"        -mask brain+tlrc -dataTableFile site2.txt\n"
"\n"
"  # Allow each run its own W, tied loosely to the subject's anchor\n"
"  3dSRM -prefix movie_rc -nfeatures 50 -run_coupling 200 -cv run \134\n"
"        -cv_adapt 0.25 -mask brain+tlrc -dataTableFile subjects.txt\n"
"\n"
"  # Choose k\n"
"  3dSRM -prefix scan -nfeatures 10,20,30,50,80 -cv run \134\n"
"        -mask brain+tlrc -dataTableFile subjects.txt\n"
"\n"
"  # subjects.txt                       # native.txt\n"
"  Subj  Group    Run  InputFile        Subj  Run  Mask             InputFile\n"
"  s01   control  r1   s01_r1+tlrc      s01   r1   s01_mask+orig    s01_r1+orig\n"
"  s01   control  r2   s01_r2+tlrc      s01   r2   s01_mask+orig    s01_r2+orig\n"
"  s02   control  r1   s02_r1+tlrc      s02   r1   s02_mask+orig    s02_r1+orig\n"
"  s02   control  r2   s02_r2+tlrc      s02   r2   s02_mask+orig    s02_r2+orig\n"
"\n"
"Reading the results:                                        ~1~\n"
"\n"
"  * The components have no fixed meaning. Any rotation of S, applied with the\n"
"    opposite rotation to every W_i, fits equally well, and other seeds find\n"
"    others. Do not test component 3's map across subjects, or compare\n"
"    component numbers between runs of the program.\n"
"  * Quantities that do not depend on that rotation are safe: the reconstruction\n"
"    R^2, the mean projected-space correlation, and segment-matching accuracy.\n"
"    Use their held-out versions for inference.\n"
"  * Under -cv loso every subject's score comes from a model that never saw\n"
"    them, so the scores are fair to compare between subjects and groups.\n"
"    Under -cv run the subjects all contributed to the model that scores them\n"
"    on other runs, which suits comparing methods, k, or ROIs, but is less\n"
"    clean for between-subject inference.\n"
"  * To compare groups, use -cv_fit and -balance. A pooled model is built more\n"
"    from the larger or less noisy group and maps that group better, with\n"
"    nothing to do with the groups differing in function.\n"
"  * The model fits samples in proportion to their variance, so a condition\n"
"    with a large response can dominate and a weak one be fitted poorly even\n"
"    in everyone. Score by condition, and check that a weak condition is weak\n"
"    in 'within' too before reading it as a failure of SRM.\n"
"  * If you will RSA the aligned task data, align on other data (-apply_W),\n"
"    or score with held-out folds. SRM maximizes cross-subject agreement of\n"
"    its own columns, which inflates their shared structure.\n"
"  * Larger k always fits the training data better. Choose k by held-out\n"
"    score, not by training fit.\n"
"\n"
"Scope boundary:                                             ~1~\n"
"\n"
"  Version 1 fits and validates SRM on time-locked or item-indexed columns,\n"
"  maps new subjects into a saved model and optionally updates the model's S\n"
"  from them, carries a fitted W to other data, allows a separate run-wise W\n"
"  tied to a subject anchor, and scores reliability by subject, group, and\n"
"  condition. Not yet provided: a prior on S that depends on condition,\n"
"  exact (data-based) refitting from a model alone, searchlight or\n"
"  atlas-parcel looping, robust SRM, and time-lagged or convolved fusion\n"
"  of EEG with fMRI. The fitting library\n"
"  (thd_srm) is written so that these can be added without changing the\n"
"  options above.\n"
"\n") ;
   PRINT_AFNI_OMP_USAGE(PROGRAM_NAME,NULL) ;
   PRINT_COMPILE_DATE ;
}


/*============================================================================
  Options
==============================================================================*/

#define GIB 1073741824.0

#define PROG_AUTO 0
#define PROG_BAR  1
#define PROG_LINE 2
#define PROG_OFF  3

typedef struct {
   char   *prefix ;
   int     nk , *klist ;             /* -nfeatures */
   SRM_params par ;                  /* algorithm, niter, tol, seed, coupling, ... */
   int     nstart ;
   char   *maskname , *run_col , *group_col , *censor_name ;
   int     zscore ;                  /* SRM_ZSCORE_* */
   int     zcensor ; float zcensor_warn ; int zcensor_warn_given ;
   int     axis_item ; char *item_col ; int item_col_given ; char *repeat_col ;
   int     repeated_stimulus ;
   char   *cvmode ;                  /* none | run | loso */
   int     cv_blocks , match_window , loso_iter ; int cv_blocks_given , match_window_given ;
   char   *cv_metrics , *fold_col , *score_by , *cv_fit ;
   int     balance , nsplit , nperm ;
   char   *reliability ;
   double  cv_adapt ; int cv_adapt_given ;
   char   *model_in ; int update_S ;
   char   *apply_W ;
   int     save_W , save_recon , save_r2 , no_projected ;
   int     proj_dset ;               /* 1: dataset, 0: 1D */
   double  memory_limit_gib ; int memory_limit_given , memory_override ;
   int     quiet , progress ;
   int     show_table ;
   int     ntl ; char **tl_name , **tl_file ;
   int     free_runs ;               /* run lengths may differ across subjects (applying a saved W) */
   char   *concat ;                  /* -concat: run starts within each subject's one input */
   int     given_algo , given_niter , given_tol , given_seed , given_zscore , given_coupling , given_mask ;
} srm_opts ;

static void srm_opts_default( srm_opts *o )
{
   memset(o,0,sizeof(*o)) ;
   SRM_params_default(&o->par) ;
   { char *ek=my_getenv("AFNI_SRM_KERNEL") ;      /* see Resources in the help */
     if( ek!=NULL && strcasecmp(ek,"VOXEL")==0 ) o->par.kernel=SRM_KERNEL_VOXEL ;
     else if( ek!=NULL && strcasecmp(ek,"GRAM")==0 ) o->par.kernel=SRM_KERNEL_GRAM ; }
   o->prefix=strdup("SRM") ;
   o->nstart=1 ; o->zscore=SRM_ZSCORE_RUN ; o->zcensor_warn=20.0f ;
   o->item_col=strdup("Item") ; o->repeat_col=strdup("Repeat") ; o->run_col=strdup("Run") ; o->group_col=strdup("Group") ;
   o->cvmode=strdup("none") ; o->cv_blocks=2 ; o->match_window=9 ; o->loso_iter=5 ;
   o->proj_dset=1 ; o->progress=PROG_AUTO ;
}

/*! Comma list of positive integers. */
static int * srm_parse_intlist( const char *s , int *n )
{
   char *copy=strdup(s) , *tok ; int *v=NULL , nn=0 ;
   for( tok=strtok(copy,",") ; tok!=NULL ; tok=strtok(NULL,",") ){
     char *end ; long x=strtol(tok,&end,10) ;
     if( *tok=='\0' || *end!='\0' || x<1 || x>100000 ){ free(copy) ; free(v) ; return NULL ; }
     v=(int *)realloc(v,sizeof(int)*(nn+1)) ; v[nn++]=(int)x ;
   }
   free(copy) ; *n=nn ;
   return v ;
}

/*============================================================================
  Progress: a bar like 3dRSA's and 3dInSync's, labelled by what is going on
==============================================================================*/

typedef struct {
   int enabled , mode , is_tty , total , done , last_bucket ;
   double start , last ;
   char label[24] , unit[24] ;
} srm_prog ;

static void srm_prog_begin( srm_prog *p , int mode , int quiet , const char *label ,
                            const char *unit , int total )
{
   memset(p,0,sizeof(*p)) ;
   p->enabled=!quiet && mode!=PROG_OFF && total>0 ;
   p->is_tty=isatty(fileno(stderr)) ; p->total=total ;
   p->mode=(mode==PROG_AUTO) ? (p->is_tty ? PROG_BAR : PROG_LINE) : mode ;
   p->start=p->last=0.001*(double)NI_clock_time() ;
   p->last_bucket=-1 ;
   snprintf(p->label,sizeof(p->label),"%s",label) ; snprintf(p->unit,sizeof(p->unit),"%s",unit) ;
   if( p->enabled && p->mode==PROG_BAR ){
     fprintf(stderr,"++ %-8s [--------------------] 0/%d %s",p->label,total,p->unit) ;
     fflush(stderr) ;
   }
}

/*! Report done (absolute, 0..total).  Bars repaint at most once a second;
    line mode reports every tenth. */
static void srm_prog_set( srm_prog *p , int done , const char *extra )
{
   double now ; int bucket ;
   if( !p->enabled ) return ;
   if( done>p->total ) done=p->total ;
   p->done=done ;
   now=0.001*(double)NI_clock_time() ; bucket=(10*done)/p->total ;
   if( done==p->total || (p->mode==PROG_BAR && now-p->last>=1.0) ||
       (p->mode==PROG_LINE && bucket>p->last_bucket) ){
     int ii , fill=(20*done)/p->total ; char tm[48] ; double el=now-p->start ;
     if( done<p->total && 10*done>=p->total && done>0 )
       snprintf(tm,sizeof(tm)," [%.0fs, about %.0fs to go]",el,el*(p->total-done)/done) ;
     else snprintf(tm,sizeof(tm)," [%.0fs]",el) ;
     if( p->mode==PROG_BAR ){
       fprintf(stderr,"\r++ %-8s [",p->label) ;
       for( ii=0 ; ii<20 ; ii++ ) fputc(ii<fill ? '#' : '-',stderr) ;
       fprintf(stderr,"] %d/%d %s%s%s%s",done,p->total,p->unit,extra?", ":"",extra?extra:"",tm) ;
       if( p->is_tty ) fputs("\033[K",stderr) ;
       if( done==p->total ) fputc('\n',stderr) ;
       fflush(stderr) ;
     } else {
       INFO_message("%s: %d/%d %s%s%s%s",p->label,done,p->total,p->unit,extra?", ":"",extra?extra:"",tm) ;
     }
     p->last=now ; p->last_bucket=bucket ;
   }
}

/*! End a bar that stopped short of its total (a fit that converged early). */
static void srm_prog_end( srm_prog *p )
{
   if( !p->enabled || p->done>=p->total ) return ;
   if( p->mode==PROG_BAR ){
     if( p->is_tty ) fputs("\033[K",stderr) ;
     fputc('\n',stderr) ; fflush(stderr) ;
   }
}

/* progress for the fit: iterations of the current start */
typedef struct { srm_prog *bar ; int istart , nstart ; double first_obj ; const char *tag ; } fit_prog ;

static void srm_fit_progress( int it , int niter , double obj , void *vctx )
{
   fit_prog *c=(fit_prog *)vctx ; char extra[96] ;
   if( it==0 ){ c->istart++ ; c->first_obj=obj ; }
   if( c->nstart>1 ) snprintf(extra,sizeof(extra),"%s%sstart %d/%d, objective %.8g",c->tag?c->tag:"",c->tag?", ":"",c->istart,c->nstart,obj) ;
   else              snprintf(extra,sizeof(extra),"%s%sobjective %.8g",c->tag?c->tag:"",c->tag?", ":"",obj) ;
   srm_prog_set(c->bar,it,extra) ;
}

/*============================================================================
  Subjects and their data
==============================================================================*/

typedef struct {
   char   *label ;
   int     group ;
   int     nrun ;
   int    *row ;                     /* table row of each run, in run order */
   char   *mask_name ;               /* this subject's Mask, or NULL */
   int     is1D ;
   THD_3dim_dataset *geom ;          /* first dataset: geometry for outputs */
   int     nvox_mask , nvox , ndrop ;
   int    *vox ;                     /* grid index of each feature (datasets) */
   float  *x ;                       /* nvox x T, retained timeline, normalized */
   float  *xt ;                      /* spare */
   /* -axis item: column j is the model's item j */
   int    *irow , *irun ;            /* [nitem] its table row, and the run it is in (0 with no Run column) */
   int    *ifile , *ibrick ;         /* [nitem] which of file[] holds it, and the sub-brick there */
   char  **file ; int nfile ;        /* the distinct datasets (selector removed), each read once */
   /* one input holding all runs (-concat, or an Xmat column) */
   char   *xmat ;                    /* this subject's X.xmat.1D, or NULL */
   int    *good , ngood ;            /* its GoodList: the uncensored time points (NULL: all) */
} srm_subj ;

static int srm_is1D( const char *name )
{
   const char *br=strchr(name,'[') ; size_t n=br ? (size_t)(br-name) : strlen(name) ;
   if( strncmp(name,"1D:",3)==0 ) return 1 ;
   return n>3 && strncasecmp(name+n-3,".1D",3)==0 ;
}

/*! File-name-safe form of a label. */
static void srm_safe( char *dst , size_t nd , const char *src )
{
   size_t i ;
   for( i=0 ; i+1<nd && src[i] ; i++ ){
     char c=src[i] ;
     dst[i]=(isalnum((unsigned char)c) || c=='_' || c=='-' || c=='.') ? c : '_' ;
   }
   dst[i]='\0' ;
}

/*! Does a dataset with this prefix already exist, under any of AFNI's names? */
/*! Does a dataset with this prefix exist?  If so and found!=NULL, its openable name goes in found. */
static int srm_dset_find( const char *prefix , char *found , size_t nfound )
{
   static const char *suf[]={"+orig.HEAD","+tlrc.HEAD","+acpc.HEAD",".nii",".nii.gz",NULL} ;
   char fn[THD_MAX_NAME] ; int ii ;
   for( ii=0 ; suf[ii]!=NULL ; ii++ ){
     snprintf(fn,sizeof(fn),"%s%s",prefix,suf[ii]) ;
     if( THD_is_file(fn) ){ if( found!=NULL ) snprintf(found,nfound,"%s",fn) ; return 1 ; }
   }
   if( THD_is_file((char *)prefix) ){ if( found!=NULL ) snprintf(found,nfound,"%s",prefix) ; return 1 ; }
   return 0 ;
}

static int srm_dset_exists( const char *prefix ){ return srm_dset_find(prefix,NULL,0) ; }

typedef struct { char *name ; byte *mask ; int nvox , nmask ; THD_3dim_dataset *ds ; } mask_entry ;

static mask_entry * srm_mask_get( mask_entry **cache , int *ncache , const char *name )
{
   int ii ; mask_entry *e ;
   for( ii=0 ; ii<*ncache ; ii++ ) if( strcmp((*cache)[ii].name,name)==0 ) return (*cache)+ii ;
   *cache=(mask_entry *)realloc(*cache,sizeof(mask_entry)*(*ncache+1)) ;
   e=(*cache)+*ncache ; (*ncache)++ ;
   e->name=strdup(name) ;
   e->ds=THD_open_dataset((char *)name) ; CHECK_OPEN_ERROR(e->ds,name) ;
   DSET_load(e->ds) ; CHECK_LOAD_ERROR(e->ds) ;
   e->nvox=DSET_NVOX(e->ds) ;
   e->mask=THD_makemask(e->ds,0,0.5f,0.0f) ;
   if( e->mask==NULL ) ERROR_exit("cannot make a mask from '%s'",name) ;
   e->nmask=THD_countmask(e->nvox,e->mask) ;
   if( e->nmask<1 ) ERROR_exit("mask '%s' has no nonzero voxels",name) ;
   return e ;
}

/*! Read dataset ds (masked voxels vox[0..nvox-1]) into x, a nvox x ldx array,
    columns col0 ... col0+nt-1.  Also flags, for -zcensor, volumes that are all
    zero in the mask.  Returns the number of non-finite values met (and the
    first one's location) so the caller can name the culprit. */
static int srm_read_dset( THD_3dim_dataset *ds , const int *vox , int nvox , float *x ,
                          int ldx , int col0 , byte *zero_flag , int *bad_t )
{
   int nt=DSET_NVALS(ds) , tt , nbad=0 ;
   DSET_load(ds) ; CHECK_LOAD_ERROR(ds) ;
   *bad_t=-1 ;
   for( tt=0 ; tt<nt ; tt++ ){
     MRI_IMAGE *fim=THD_extract_float_brick(tt,ds) ;
     float *far ; int vv , nz=0 ;
     if( fim==NULL ) ERROR_exit("cannot read volume %d of '%s'",tt,DSET_BRIKNAME(ds)) ;
     far=MRI_FLOAT_PTR(fim) ;
     for( vv=0 ; vv<nvox ; vv++ ){
       float f=far[vox[vv]] ;
       if( !isfinite(f) ){ nbad++ ; if( *bad_t<0 ) *bad_t=tt ; f=0.0f ; }
       if( f!=0.0f ) nz=1 ;
       x[(size_t)vv*ldx+col0+tt]=f ;
     }
     if( zero_flag!=NULL ) zero_flag[col0+tt]=(byte)!nz ;
     mri_free(fim) ;
   }
   DSET_unload(ds) ;
   return nbad ;
}

/*============================================================================
  Output helpers
==============================================================================*/

static void srm_write_matrix_1D( const char *fname , const float *a , int nrow , int ncol ,
                                 const char *comment )
{
   FILE *fp=fopen(fname,"w") ; int i , j ;
   if( fp==NULL ) ERROR_exit("cannot write '%s'",fname) ;
   if( comment!=NULL ) fprintf(fp,"# %s\n",comment) ;
   for( i=0 ; i<nrow ; i++ ){
     for( j=0 ; j<ncol ; j++ ) fprintf(fp,"%s%.7g",j?" ":"",a[(size_t)i*ncol+j]) ;
     fputc('\n',fp) ;
   }
   fclose(fp) ;
}

/*! Transposed write: a is ncol x nrow row-major (components by time); the file
    has one row per time point. */
static void srm_write_transposed_1D( const char *fname , const float *a , int ncomp , int ntime ,
                                     const char *comment )
{
   FILE *fp=fopen(fname,"w") ; int i , j ;
   if( fp==NULL ) ERROR_exit("cannot write '%s'",fname) ;
   if( comment!=NULL ) fprintf(fp,"# %s\n",comment) ;
   for( j=0 ; j<ntime ; j++ ){
     for( i=0 ; i<ncomp ; i++ ) fprintf(fp,"%s%.7g",i?" ":"",a[(size_t)i*ntime+j]) ;
     fputc('\n',fp) ;
   }
   fclose(fp) ;
}

/*! Dataset of nbrick volumes on ds's grid; vol(b) returns the nvox-sized
    buffer for brick b (the dataset takes ownership). */
static void srm_write_grid_dset( THD_3dim_dataset *geom , const char *prefix , int nbrick ,
                                 float **vols , const char *label_stem , int argc , char **argv )
{
   THD_3dim_dataset *out=EDIT_empty_copy(geom) ; int b ; char lab[64] ;
   EDIT_dset_items(out,ADN_prefix,prefix,ADN_nvals,nbrick,ADN_ntt,0,ADN_brick_fac,NULL,
                   ADN_type,HEAD_FUNC_TYPE,ADN_func_type,FUNC_BUCK_TYPE,ADN_datum_all,MRI_float,ADN_none) ;
   tross_Copy_History(geom,out) ; tross_Make_History("3dSRM",argc,argv,out) ;
   for( b=0 ; b<nbrick ; b++ ){
     EDIT_substitute_brick(out,b,MRI_float,vols[b]) ;
     snprintf(lab,sizeof(lab),"%s%d",label_stem,b+1) ;
     EDIT_BRICK_LABEL(out,b,lab) ; EDIT_BRICK_TO_NOSTAT(out,b) ;
   }
   DSET_write(out) ; WROTE_DSET(out) ; DSET_delete(out) ;
}

/*============================================================================
  Cross-validation
==============================================================================*/

#define M_ISC    0
#define M_MATCH  1
#define M_CHANCE 2
#define M_R2     3
#define M_NSAMP  4
#define NMET     5

typedef struct { int isc , match , r2 , ident ; } srm_want ;   /* ident: 'match' is item identification */

typedef struct {
   srm_opts *o ; srm_subj *sb ; int nsub , nrun , T , K ;
   int *seglen ; char **runname ; char **glabel ; int ngroup ;
   SRM_params par ;                       /* algorithm, K, center, coupling ... */
   int quiet , progress ;
} srm_ctx ;

typedef struct {
   int nfold , by_blocks , items ;
   int **col , *len ;                     /* held-out columns of each fold, ascending */
   char **name ;
   const char *unit ;                     /* "time points" or "items" */
} srm_folds ;

typedef struct {
   int nsub , nfold ;
   double *score ;                        /* [nsub][nfold][NMET] */
   float **proj ;                         /* [nsub] K x T held-out projections, or NULL */
   float **vr2 ;                          /* [nsub] mean voxelwise held-out R^2, or NULL */
} srm_cvres ;

static srm_want srm_parse_metrics( const char *s )
{
   srm_want w ; char *copy , *tok ;
   memset(&w,0,sizeof(w)) ;
   if( s==NULL ){ w.isc=w.match=w.r2=1 ; return w ; }
   copy=strdup(s) ;
   for( tok=strtok(copy,",") ; tok!=NULL ; tok=strtok(NULL,",") ){
     if( strcasecmp(tok,"isc")==0 ) w.isc=1 ;
     else if( strcasecmp(tok,"match")==0 || strcasecmp(tok,"ident")==0 ) w.match=1 ;
     else if( strcasecmp(tok,"r2")==0 ) w.r2=1 ;
     else ERROR_exit("-cv_metrics: '%s' is not isc, match, or r2",tok) ;
   }
   free(copy) ;
   if( !w.isc && !w.match && !w.r2 ) ERROR_exit("-cv_metrics names no metric") ;
   return w ;
}

/*! Folds of the time axis: the runs, or with one run, contiguous blocks of time. */
static void srm_make_folds( srm_folds *f , int nrun , const int *seglen , char **runname , int T , int blocks )
{
   int ii , jj , acc=0 ;
   memset(f,0,sizeof(*f)) ;
   f->unit="time points" ;
   if( nrun>=2 ){
     f->nfold=nrun ; f->by_blocks=0 ;
     f->len=(int *)malloc(sizeof(int)*(size_t)nrun) ; f->name=(char **)malloc(sizeof(char *)*(size_t)nrun) ;
     for( ii=0 ; ii<nrun ; ii++ ){ f->len[ii]=seglen[ii] ; f->name[ii]=strdup(runname[ii]) ; }
   } else {
     int base=T/blocks , rem=T%blocks ; char nm[24] ;
     f->nfold=blocks ; f->by_blocks=1 ;
     f->len=(int *)malloc(sizeof(int)*(size_t)blocks) ; f->name=(char **)malloc(sizeof(char *)*(size_t)blocks) ;
     for( ii=0 ; ii<blocks ; ii++ ){
       f->len[ii]=base+(ii<rem) ;
       snprintf(nm,sizeof(nm),"block%d",ii+1) ; f->name[ii]=strdup(nm) ;
     }
   }
   f->col=(int **)malloc(sizeof(int *)*(size_t)f->nfold) ;
   for( ii=0 ; ii<f->nfold ; ii++ ){
     f->col[ii]=(int *)malloc(sizeof(int)*(size_t)f->len[ii]) ;
     for( jj=0 ; jj<f->len[ii] ; jj++ ) f->col[ii][jj]=acc+jj ;
     acc+=f->len[ii] ;
   }
}

static void srm_free_folds( srm_folds *f )
{
   int ii ;
   for( ii=0 ; ii<f->nfold ; ii++ ){ free(f->name[ii]) ; free(f->col[ii]) ; }
   free(f->name) ; free(f->col) ; free(f->len) ;
}

/*! Columns col[0..len-1] of x (V x T), in that order, as a new V x len array. */
static float * srm_cols_in( const float *x , int V , int T , const int *col , int len )
{
   float *y=(float *)malloc(sizeof(float)*(size_t)V*(size_t)len) ; int f , t ;
   if( y==NULL ) ERROR_exit("out of memory building a fold (%.2f GiB)",4.0*V*(double)len/GIB) ;
   for( f=0 ; f<V ; f++ ){
     const float *src=x+(size_t)f*T ; float *dst=y+(size_t)f*len ;
     for( t=0 ; t<len ; t++ ) dst[t]=src[col[t]] ;
   }
   return y ;
}

/*! All columns of x except col[0..len-1], in their order, as a new V x (T-len) array. */
static float * srm_cols_out( const float *x , int V , int T , const int *col , int len )
{
   int Tr=T-len , f , t , q ; byte *out=(byte *)calloc((size_t)T,1) ;
   float *y=(float *)malloc(sizeof(float)*(size_t)V*(size_t)Tr) ;
   if( y==NULL || out==NULL ) ERROR_exit("out of memory building a fold (%.2f GiB)",4.0*V*(double)Tr/GIB) ;
   for( t=0 ; t<len ; t++ ) out[col[t]]=1 ;
   for( f=0 ; f<V ; f++ ){
     const float *src=x+(size_t)f*T ; float *dst=y+(size_t)f*Tr ;
     for( t=0,q=0 ; t<T ; t++ ) if( !out[t] ) dst[q++]=src[t] ;
   }
   free(out) ;
   return y ;
}

/*! Same for a K x T matrix. */
static float * srm_mat_cols_out( const float *s , int K , int T , const int *col , int len )
{ return srm_cols_out(s,K,T,col,len) ; }

static void srm_cvres_alloc( srm_cvres *r , int nsub , int nfold , int want_proj , int want_vr2 , const srm_ctx *cx )
{
   int ii ;
   memset(r,0,sizeof(*r)) ;
   r->nsub=nsub ; r->nfold=nfold ;
   r->score=(double *)calloc((size_t)nsub*nfold*NMET,sizeof(double)) ;
   if( want_proj ){
     r->proj=(float **)calloc((size_t)nsub,sizeof(float *)) ;
     for( ii=0 ; ii<nsub ; ii++ ) r->proj[ii]=(float *)calloc((size_t)cx->K*cx->T,sizeof(float)) ;
   }
   if( want_vr2 ){
     r->vr2=(float **)calloc((size_t)nsub,sizeof(float *)) ;
     for( ii=0 ; ii<nsub ; ii++ ) r->vr2[ii]=(float *)calloc((size_t)cx->sb[ii].nvox,sizeof(float)) ;
   }
}

static void srm_cvres_free( srm_cvres *r )
{
   int ii ;
   if( r->proj!=NULL ){ for( ii=0 ; ii<r->nsub ; ii++ ) free(r->proj[ii]) ; free(r->proj) ; }
   if( r->vr2!=NULL ){ for( ii=0 ; ii<r->nsub ; ii++ ) free(r->vr2[ii]) ; free(r->vr2) ; }
   free(r->score) ;
}

/*! Score one subject's held-out columns.  P is the subject's projected test
    data (K x len); ref_u the reference shared response for correlation and
    matching; ref_r the one used to reconstruct (the posterior mean for a PROB model).
    model/iview locate the subject's W for the reconstruction. */
static int srm_score_one( const SRM_model *model , int iview , const float *xte , const float *P ,
                          const float *ref_u , const float *ref_r , int K , int len , int win ,
                          srm_want w , double *out , float *vox_r2 )
{
   int rc=SRM_OK ; double a , b , r2 ;
   out[M_NSAMP]=len ; out[M_ISC]=out[M_MATCH]=out[M_CHANCE]=out[M_R2]=NAN ;
   if( w.isc ){
     rc=SRM_component_corr(P,ref_u,K,len,NULL,&a) ;
     if( rc!=SRM_OK ) return rc ;
     out[M_ISC]=a ;
   }
   if( w.match ){
     rc = w.ident ? SRM_item_identify(P,ref_u,K,len,&a,&b) : SRM_segment_match(P,ref_u,K,len,win,&a,&b) ;
     if( rc!=SRM_OK ) return rc ;
     out[M_MATCH]=a ; out[M_CHANCE]=b ;
   }
   if( w.r2 ){
     rc=SRM_reconstruction_r2(model,iview,xte,ref_r,len,&r2,vox_r2) ;
     if( rc!=SRM_OK ) return rc ;
     out[M_R2]=r2 ;
   }
   return SRM_OK ;
}

static void srm_check_fold_lengths( const srm_folds *f , srm_want w , int win )
{
   int ii , minlen=1<<30 ;
   for( ii=0 ; ii<f->nfold ; ii++ ) if( f->len[ii]<minlen ) minlen=f->len[ii] ;
   if( minlen<3 ) ERROR_exit("the smallest fold (%s) has %d %s; need at least 3",f->name[0],minlen,f->unit) ;
   if( f->items ) return ;                       /* item identification has no window */
   if( w.match && minlen<2*win+1 )
     ERROR_exit("time-segment matching with a %d-TR window needs folds of at least %d time points; the smallest has %d "
                "(use a smaller -match_window, or drop 'match' from -cv_metrics)",win,2*win+1,minlen) ;
}

/*! -cv run: leave one run (or block) out, refitting the group each time. */
static void srm_cv_by_runs( const srm_ctx *cx , const SRM_model *mfull , const srm_folds *f , srm_want w ,
                            srm_cvres *res )
{
   int nsub=cx->nsub , K=cx->K , T=cx->T , ff , ii , jj ;
   const srm_opts *o=cx->o ;
   float **S_tr=(float **)calloc((size_t)f->nfold,sizeof(float *)) ;      /* each fold's training S, for the gauge */
   float ***Pf=(float ***)calloc((size_t)f->nfold,sizeof(float **)) ;     /* [fold][subject] projections */
   int *nvox=(int *)malloc(sizeof(int)*(size_t)nsub) ;

   for( ii=0 ; ii<nsub ; ii++ ) nvox[ii]=cx->sb[ii].nvox ;
   for( ff=0 ; ff<f->nfold ; ff++ ){
     int len=f->len[ff] , *col=f->col[ff] , Ttr=T-len ;
     float **xtr=(float **)malloc(sizeof(float *)*(size_t)nsub) , **xte=(float **)malloc(sizeof(float *)*(size_t)nsub) ;
     SRM_params par=cx->par ; SRM_model *mf=NULL ; SRM_fitinfo fi ; srm_prog bar ; fit_prog fp ;
     int segtr[256] , rc , *all=(int *)malloc(sizeof(int)*(size_t)nsub) ;
     float **P=(float **)malloc(sizeof(float *)*(size_t)nsub) ;
     float *ref_u=(float *)malloc(sizeof(float)*(size_t)K*len) , *ref_r=(float *)malloc(sizeof(float)*(size_t)K*len) ;
     char tag[96] ;

     for( ii=0 ; ii<nsub ; ii++ ){
       xtr[ii]=srm_cols_out(cx->sb[ii].x,nvox[ii],T,col,len) ;
       xte[ii]=srm_cols_in(cx->sb[ii].x,nvox[ii],T,col,len) ;
       all[ii]=ii ;
     }
     if( par.coupling!=SRM_COUPLING_SHARED ){        /* the runs that remain, with their lengths */
       int q=0 ;
       if( f->by_blocks ) ERROR_exit("-run_coupling needs runs; with one run -cv uses time blocks") ;
       for( jj=0 ; jj<cx->nrun ; jj++ ) if( jj!=ff ) segtr[q++]=cx->seglen[jj] ;
       if( q<2 ) ERROR_exit("-run_coupling with -cv run needs at least 3 runs (2 remain when one is held out)") ;
       par.nseg=q ; par.seglen=segtr ;
     }
     snprintf(tag,sizeof(tag),"%s %d/%d",f->items?"fold":f->by_blocks?"block":"run",ff+1,f->nfold) ;
     srm_prog_begin(&bar,cx->progress,cx->quiet,"fold","iterations",par.niter) ;
     fp.bar=&bar ; fp.istart=0 ; fp.nstart=o->nstart ; fp.first_obj=0.0 ; fp.tag=tag ;
     par.progress=srm_fit_progress ; par.progress_ctx=&fp ;
     rc=SRM_fit_multistart(nsub,Ttr,nvox,xtr,&par,o->nstart,&mf,&fi,NULL,NULL) ;
     srm_prog_end(&bar) ;
     if( rc!=SRM_OK ) ERROR_exit("the fit for %s failed: %s",tag,SRM_strerror(rc)) ;
     if( !cx->quiet ) ININFO_message("%s (%s): %s after %d iterations",tag,f->name[ff],fi.converged?"converged":"iteration limit",fi.niter_done) ;
     rc=SRM_canonicalize(mf) ;
     if( rc!=SRM_OK ) ERROR_exit("could not canonicalize the fit for %s: %s",tag,SRM_strerror(rc)) ;

     for( ii=0 ; ii<nsub ; ii++ ){                    /* each subject's held-out data in the shared space, once */
       P[ii]=(float *)malloc(sizeof(float)*(size_t)K*len) ;
       rc=SRM_project(mf,ii,xte[ii],len,P[ii]) ;
       if( rc!=SRM_OK ) ERROR_exit("projection failed for Subj %s: %s",cx->sb[ii].label,SRM_strerror(rc)) ;
     }
     for( ii=0 ; ii<nsub ; ii++ ){                    /* the reference for subject ii: everyone else */
       int *oth=(int *)malloc(sizeof(int)*(size_t)nsub) , no=0 ; float **po=(float **)malloc(sizeof(float *)*(size_t)nsub) ;
       double *sc=res->score+((size_t)ii*f->nfold+ff)*NMET ; float *vr=(res->vr2!=NULL) ? (float *)malloc(sizeof(float)*(size_t)nvox[ii]) : NULL ;
       for( jj=0 ; jj<nsub ; jj++ ) if( jj!=ii ){ oth[no]=jj ; po[no]=P[jj] ; no++ ; }
       rc=SRM_shared_from_projections(mf,no,oth,po,len,0,ref_u) ;
       if( rc==SRM_OK ) rc=SRM_shared_from_projections(mf,no,oth,po,len,1,ref_r) ;
       if( rc==SRM_OK ) rc=srm_score_one(mf,ii,xte[ii],P[ii],ref_u,ref_r,K,len,o->match_window,w,sc,vr) ;
       if( rc!=SRM_OK ) ERROR_exit("scoring Subj %s on %s failed: %s",cx->sb[ii].label,tag,SRM_strerror(rc)) ;
       if( vr!=NULL ){ int q ; for( q=0 ; q<nvox[ii] ; q++ ) res->vr2[ii][q]+=vr[q]/(float)f->nfold ; free(vr) ; }
       free(oth) ; free(po) ;
     }
     if( res->proj!=NULL ){
       Pf[ff]=P ; S_tr[ff]=(float *)malloc(sizeof(float)*(size_t)K*Ttr) ; memcpy(S_tr[ff],mf->s,sizeof(float)*(size_t)K*Ttr) ;
     } else { for( ii=0 ; ii<nsub ; ii++ ) free(P[ii]) ; free(P) ; }
     for( ii=0 ; ii<nsub ; ii++ ){ free(xtr[ii]) ; free(xte[ii]) ; }
     free(xtr) ; free(xte) ; free(all) ; free(ref_u) ; free(ref_r) ;
     SRM_free_fitinfo(&fi) ; SRM_free_model(mf) ;
   }

   /* held-out projections in one gauge: each fold's S aligned to the full fit's, same rotation for all subjects.
      This changes labels of components, never the scores (they were computed before this and are rotation-invariant). */
   if( res->proj!=NULL ){
     for( ff=0 ; ff<f->nfold ; ff++ ){
       int len=f->len[ff] , *col=f->col[ff] , Ttr=T-len ; float R[1]={0} ;
       float *Sfull_tr=srm_mat_cols_out(mfull->s,K,T,col,len) ; float *Rm=(float *)malloc(sizeof(float)*(size_t)K*K) ;
       int rc=SRM_procrustes_rotation(S_tr[ff],Sfull_tr,K,Ttr,Rm) ;
       (void)R ;
       if( rc!=SRM_OK ) ERROR_exit("could not align the gauge of %s: %s",f->name[ff],SRM_strerror(rc)) ;
       for( ii=0 ; ii<nsub ; ii++ ){
         int c , d , t ;
         for( c=0 ; c<K ; c++ ) for( t=0 ; t<len ; t++ ){
           double v=0.0 ; for( d=0 ; d<K ; d++ ) v+=(double)Rm[c*K+d]*Pf[ff][ii][d*len+t] ;
           res->proj[ii][(size_t)c*T+col[t]]=(float)v ;
         }
         free(Pf[ff][ii]) ;
       }
       free(Pf[ff]) ; free(S_tr[ff]) ; free(Rm) ; free(Sfull_tr) ;
     }
   }
   free(S_tr) ; free(Pf) ; free(nvox) ;
}

/*! -cv loso: fit the group to everyone else, then estimate this subject's W from their own
    training folds and score them on the held-out fold.  Neither S nor any other W ever saw the
    subject's test data. */
static void srm_cv_loso( const srm_ctx *cx , const SRM_model *mfull , const srm_folds *f , srm_want w ,
                         srm_cvres *res )
{
   int nsub=cx->nsub , K=cx->K , T=cx->T , ii , jj , ff ; srm_prog bar ;
   const srm_opts *o=cx->o ;
   srm_prog_begin(&bar,cx->progress,cx->quiet,"loso","subjects",nsub) ;
   for( ii=0 ; ii<nsub ; ii++ ){
     SRM_params par=cx->par ; SRM_model *mdrop , *ml=NULL ; SRM_fitinfo fi ;
     int *nv=(int *)malloc(sizeof(int)*(size_t)(nsub-1)) ; float **xl=(float **)malloc(sizeof(float *)*(size_t)(nsub-1)) ;
     int q=0 , rc ; float *Ralign=(float *)malloc(sizeof(float)*(size_t)K*K) ;
     for( jj=0 ; jj<nsub ; jj++ ) if( jj!=ii ){ nv[q]=cx->sb[jj].nvox ; xl[q]=cx->sb[jj].x ; q++ ; }
     mdrop=SRM_model_copy(mfull,ii) ;
     if( mdrop==NULL ) ERROR_exit("out of memory copying the model for leave-one-out") ;
     par.init=SRM_INIT_MODEL ; par.init_model=mdrop ; par.niter=o->loso_iter ; par.progress=NULL ; par.progress_ctx=NULL ;
     if( par.coupling!=SRM_COUPLING_SHARED ){ par.nseg=cx->nrun ; par.seglen=cx->seglen ; }
     rc=SRM_fit(nsub-1,T,nv,xl,&par,&ml,&fi) ;
     if( rc!=SRM_OK ) ERROR_exit("the leave-one-out fit without Subj %s failed: %s",cx->sb[ii].label,SRM_strerror(rc)) ;
     SRM_free_fitinfo(&fi) ; SRM_free_model(mdrop) ;
     rc=SRM_procrustes_rotation(ml->s,mfull->s,K,T,Ralign) ;
     if( rc!=SRM_OK ) ERROR_exit("could not align the gauge for Subj %s: %s",cx->sb[ii].label,SRM_strerror(rc)) ;

     for( ff=0 ; ff<f->nfold ; ff++ ){                 /* folds within this subject */
       int len=f->len[ff] , *col=f->col[ff] , Ttr=T-len , V=cx->sb[ii].nvox ;
       float *xtr=srm_cols_out(cx->sb[ii].x,V,T,col,len) , *xte=srm_cols_in(cx->sb[ii].x,V,T,col,len) ;
       float *str=srm_mat_cols_out(ml->s,K,T,col,len) , *ste=srm_cols_in(ml->s,K,T,col,len) ;
       float *Wf=(float *)malloc(sizeof(float)*(size_t)V*K) , *P=(float *)malloc(sizeof(float)*(size_t)K*len) ;
       float *mu=cx->par.center ? (float *)malloc(sizeof(float)*(size_t)V) : NULL ; SRM_model *mv ;
       double *sc=res->score+((size_t)ii*f->nfold+ff)*NMET ; float *vr=(res->vr2!=NULL) ? (float *)malloc(sizeof(float)*(size_t)V) : NULL ;
       rc=SRM_fit_view(str,K,Ttr,xtr,V,cx->par.center,Wf,NULL,mu) ;     /* W from the training columns only */
       if( rc!=SRM_OK ) ERROR_exit("Subj %s, %s: could not estimate W: %s",cx->sb[ii].label,f->name[ff],SRM_strerror(rc)) ;
       mv=SRM_model_from_w(V,K,Wf,mu) ;
       rc=SRM_project(mv,0,xte,len,P) ;
       if( rc==SRM_OK ) rc=srm_score_one(mv,0,xte,P,ste,ste,K,len,o->match_window,w,sc,vr) ;
       if( rc!=SRM_OK ) ERROR_exit("Subj %s, %s: scoring failed: %s",cx->sb[ii].label,f->name[ff],SRM_strerror(rc)) ;
       if( vr!=NULL ){ int z ; for( z=0 ; z<V ; z++ ) res->vr2[ii][z]+=vr[z]/(float)f->nfold ; free(vr) ; }
       if( res->proj!=NULL ){
         int c , d , t ;
         for( c=0 ; c<K ; c++ ) for( t=0 ; t<len ; t++ ){
           double v=0.0 ; for( d=0 ; d<K ; d++ ) v+=(double)Ralign[c*K+d]*P[d*len+t] ;
           res->proj[ii][(size_t)c*T+col[t]]=(float)v ;
         }
       }
       SRM_free_model(mv) ; free(xtr) ; free(xte) ; free(str) ; free(ste) ; free(Wf) ; free(P) ; free(mu) ;
     }
     SRM_free_model(ml) ; free(nv) ; free(xl) ; free(Ralign) ;
     { char extra[128] ; snprintf(extra,sizeof(extra),"left out Subj %s, %d folds scored",cx->sb[ii].label,f->nfold) ;
       srm_prog_set(&bar,ii+1,extra) ; }
   }
}

/*! -cv with -model_in: S is fixed.  Each subject's W is estimated from the training columns only
    and the subject is scored on the held-out columns against the model's own S there.  The model
    never saw the subject, and the held-out columns touch neither W nor S. */
static void srm_cv_frozen( const srm_ctx *cx , const SRM_model *m , const srm_folds *f , srm_want w ,
                           srm_cvres *res )
{
   int nsub=cx->nsub , K=cx->K , T=cx->T , ii , ff ; srm_prog bar ;
   const srm_opts *o=cx->o ;
   srm_prog_begin(&bar,cx->progress,cx->quiet,"folds","subjects",nsub) ;
   for( ii=0 ; ii<nsub ; ii++ ){
     for( ff=0 ; ff<f->nfold ; ff++ ){
       int len=f->len[ff] , *col=f->col[ff] , Ttr=T-len , V=cx->sb[ii].nvox , rc ;
       float *xtr=srm_cols_out(cx->sb[ii].x,V,T,col,len) , *xte=srm_cols_in(cx->sb[ii].x,V,T,col,len) ;
       float *str=srm_mat_cols_out(m->s,K,T,col,len) , *ste=srm_cols_in(m->s,K,T,col,len) ;
       float *Wf=(float *)malloc(sizeof(float)*(size_t)V*K) , *P=(float *)malloc(sizeof(float)*(size_t)K*len) ;
       float *mu=cx->par.center ? (float *)malloc(sizeof(float)*(size_t)V) : NULL ; SRM_model *mv ;
       double *sc=res->score+((size_t)ii*f->nfold+ff)*NMET ; float *vr=(res->vr2!=NULL) ? (float *)malloc(sizeof(float)*(size_t)V) : NULL ;
       rc=SRM_fit_view(str,K,Ttr,xtr,V,cx->par.center,Wf,NULL,mu) ;
       if( rc!=SRM_OK ) ERROR_exit("Subj %s, %s: could not estimate W: %s",cx->sb[ii].label,f->name[ff],SRM_strerror(rc)) ;
       mv=SRM_model_from_w(V,K,Wf,mu) ;
       if( mv==NULL ) ERROR_exit("out of memory scoring Subj %s",cx->sb[ii].label) ;
       rc=SRM_project(mv,0,xte,len,P) ;
       if( rc==SRM_OK ) rc=srm_score_one(mv,0,xte,P,ste,ste,K,len,o->match_window,w,sc,vr) ;
       if( rc!=SRM_OK ) ERROR_exit("Subj %s, %s: scoring failed: %s",cx->sb[ii].label,f->name[ff],SRM_strerror(rc)) ;
       if( vr!=NULL ){ int z ; for( z=0 ; z<V ; z++ ) res->vr2[ii][z]+=vr[z]/(float)f->nfold ; free(vr) ; }
       if( res->proj!=NULL ){
         int c , t ;
         for( c=0 ; c<K ; c++ ) for( t=0 ; t<len ; t++ ) res->proj[ii][(size_t)c*T+col[t]]=P[(size_t)c*len+t] ;
       }
       SRM_free_model(mv) ; free(xtr) ; free(xte) ; free(str) ; free(ste) ; free(Wf) ; free(P) ; free(mu) ;
     }
     srm_prog_set(&bar,ii+1,NULL) ;
   }
}

/*! Mean over folds of one metric for one subject (NaN if not computed). */
static double srm_subj_metric( const srm_cvres *r , int ii , int met )
{
   int ff , n=0 ; double s=0.0 ;
   for( ff=0 ; ff<r->nfold ; ff++ ){
     double v=r->score[((size_t)ii*r->nfold+ff)*NMET+met] ;
     if( isfinite(v) ){ s+=v ; n++ ; }
   }
   return n>0 ? s/n : NAN ;
}

static void srm_write_scores( const srm_ctx *cx , const srm_folds *f , const srm_cvres *r , srm_want w ,
                              const char *mode , THD_datatable *tab , int *carry , int ncarry )
{
   char fn[THD_MAX_NAME] ; FILE *fp ; int ii , ff , cc ; const char *mname = w.ident ? "ident" : "match" ;
   snprintf(fn,sizeof(fn),"%s.scores_folds.txt",cx->o->prefix) ;
   fp=fopen(fn,"w") ; if( fp==NULL ) ERROR_exit("cannot write '%s'",fn) ;
   fprintf(fp,"# held-out scores, one row per subject and fold (%s; %d folds)\n",mode,f->nfold) ;
   fprintf(fp,"Subj\tGroup\tFold\tnsamp\tisc\t%s\tchance\tr2\n",mname) ;
   for( ii=0 ; ii<cx->nsub ; ii++ ) for( ff=0 ; ff<f->nfold ; ff++ ){
     const double *sc=r->score+((size_t)ii*f->nfold+ff)*NMET ;
     fprintf(fp,"%s\t%s\t%s\t%d",cx->sb[ii].label,cx->glabel[cx->sb[ii].group],f->name[ff],(int)sc[M_NSAMP]) ;
     for( cc=0 ; cc<4 ; cc++ ){ int mi[4]={M_ISC,M_MATCH,M_CHANCE,M_R2} ;
       if( isfinite(sc[mi[cc]]) ) fprintf(fp,"\t%.7g",sc[mi[cc]]) ; else fprintf(fp,"\tNA") ; }
     fputc('\n',fp) ;
   }
   fclose(fp) ;
   snprintf(fn,sizeof(fn),"%s.scores.txt",cx->o->prefix) ;
   fp=fopen(fn,"w") ; if( fp==NULL ) ERROR_exit("cannot write '%s'",fn) ;
   fprintf(fp,"# held-out scores per subject, the mean over folds (%s; %d folds). isc is the mean over components of the "
              "correlation with the others' shared response, in Fisher-z space; %s (chance in the next column); r2 is "
              "reconstruction R^2\n",mode,f->nfold,
              w.ident ? "ident is item identification accuracy: each held-out item must correlate best with its own "
                        "counterpart among the fold's items" : "match is time-segment matching accuracy") ;
   fprintf(fp,"Subj\tGroup\tFactor\tLevel\tnsamp\tisc\t%s\tchance\tr2",mname) ;
   for( cc=0 ; cc<ncarry ; cc++ ) fprintf(fp,"\t%s",tab->cname[carry[cc]]) ;
   fputc('\n',fp) ;
   for( ii=0 ; ii<cx->nsub ; ii++ ){
     int ns=0 ;
     for( ff=0 ; ff<f->nfold ; ff++ ) ns+=(int)r->score[((size_t)ii*f->nfold+ff)*NMET+M_NSAMP] ;
     fprintf(fp,"%s\t%s\tall\tall\t%d",cx->sb[ii].label,cx->glabel[cx->sb[ii].group],ns) ;
     { int mi[4]={M_ISC,M_MATCH,M_CHANCE,M_R2} ;
       for( cc=0 ; cc<4 ; cc++ ){ double v=srm_subj_metric(r,ii,mi[cc]) ; if( isfinite(v) ) fprintf(fp,"\t%.7g",v) ; else fprintf(fp,"\tNA") ; } }
     for( cc=0 ; cc<ncarry ; cc++ ) fprintf(fp,"\t%s",DT_CELL(tab,cx->sb[ii].row[0],carry[cc])) ;
     fputc('\n',fp) ;
   }
   fclose(fp) ;
}

/*! Group means (and standard errors across subjects) for the log, and for -nfeatures. */
static void srm_cv_summary( const srm_ctx *cx , const srm_cvres *r , srm_want w , const char *mode , int nfold )
{
   int g , ii , met ; const char *nm[4]={"isc",w.ident?"ident":"match","chance","r2"} ; int mi[4]={M_ISC,M_MATCH,M_CHANCE,M_R2} ;
   int sel[4]={w.isc,w.match,w.match,w.r2} ;
   if( cx->quiet ) return ;
   INFO_message("held-out scores (%s, %d folds), mean over subjects with its standard error:",mode,nfold) ;
   for( g=-1 ; g<cx->ngroup ; g++ ){
     char line[512] ; int pos=0 , n=0 ;
     if( g>=0 && cx->ngroup<2 ) break ;
     for( ii=0 ; ii<cx->nsub ; ii++ ) n+=(g<0 || cx->sb[ii].group==g) ;
     pos+=snprintf(line,sizeof(line),"%s (%d subjects):",g<0?"all":cx->glabel[g],n) ;
     for( met=0 ; met<4 ; met++ ){
       double s=0.0 , s2=0.0 ; int cnt=0 ;
       if( !sel[met] ) continue ;
       for( ii=0 ; ii<cx->nsub ; ii++ ) if( g<0 || cx->sb[ii].group==g ){
         double v=srm_subj_metric(r,ii,mi[met]) ; if( isfinite(v) ){ s+=v ; s2+=v*v ; cnt++ ; } }
       if( cnt>0 ){
         double mean=s/cnt , sd=(cnt>1) ? sqrt((s2-cnt*mean*mean)/(cnt-1)) : 0.0 ;
         if( met==2 ) pos+=snprintf(line+pos,sizeof(line)-pos,"  %s %.4g",nm[met],mean) ;
         else pos+=snprintf(line+pos,sizeof(line)-pos,"  %s %.4g (%.2g)",nm[met],mean,cnt>1?sd/sqrt((double)cnt):0.0) ;
       }
     }
     ININFO_message("%s",line) ;
   }
}

/*! A feature to drop: variance zero in some run, for zscore none. */
static int srm_flag_constant( float *x , int nvox , int ntime , int nseg , const int *seglen ,
                              unsigned char *bad )
{
   int f , nb=0 ;
   for( f=0 ; f<nvox ; f++ ){
     const float *row=x+(size_t)f*ntime ; int rr , c0=0 , isb=0 ;
     for( rr=0 ; rr<nseg && !isb ; c0+=seglen[rr],rr++ ){
       double m=0.0 , v=0.0 ; int t , L=seglen[rr] ;
       for( t=0 ; t<L ; t++ ) m+=row[c0+t] ;
       m/=L ;
       for( t=0 ; t<L ; t++ ){ double d=row[c0+t]-m ; v+=d*d ; }
       if( !(v>0.0) ) isb=1 ;
     }
     bad[f]=(unsigned char)isb ; nb+=isb ;
   }
   return nb ;
}

/*! Keep only the rows with bad[f]==0, packed to the front. */
static int srm_compact_rows( float *x , int nvox , int ntime , const unsigned char *bad , int *vox )
{
   int f , q=0 ;
   for( f=0 ; f<nvox ; f++ ) if( !bad[f] ){
     if( q!=f ) memmove(x+(size_t)q*ntime,x+(size_t)f*ntime,sizeof(float)*(size_t)ntime) ;
     if( vox!=NULL ) vox[q]=vox[f] ;
     q++ ;
   }
   return q ;
}

/*============================================================================
  The data stage, as functions: table -> subjects -> reading -> timeline -> normalized.
  Fitting, applying a saved model, and applying a saved W all go through it.
==============================================================================*/

typedef struct {
   THD_datatable *tab ; THD_datatable_index *dx ;
   int nsub , nrun , have_run , ngroup ;
   char **glabel , **runname ;
   srm_subj *sb ;
   int *runlen , Tin , *c0in ;
   int icol_mask , icol_group , icol_run ;
   mask_entry *mcache ; int nmcache ; mask_entry *gmask ;
   byte *keep ;                        /* [Tin] retained time points */
   byte *zero_any ; int *nzero ;       /* -zcensor bookkeeping */
   int T ; int *seglen , *time_index ; /* retained timeline */
   int axis_item , nitem ;             /* -axis item: T == nitem columns, one per key */
   int concat , *cstart ;              /* runs within one input: their starts (from -concat or Xmat) */
   int icol_xmat ;
   char **item_key ;                   /* [nitem] Item, or Item:n for the n-th repeat */
} srm_data ;

/*! Subjects, their runs, groups and masks, from the table. */
static void srm_data_table( const srm_opts *o , THD_datatable *tab , srm_data *d , int min_subjects )
{
   char *tab_cols[2] ; int ss , rr , jj ;
   d->tab=tab ;
   if( tab->icol_input<0 ) ERROR_exit("the data table needs an InputFile column") ;
   d->icol_run=THD_datatable_column(tab,o->run_col) ;
   d->have_run=(d->icol_run>=0) ;
   if( d->have_run && tab->icol_subj<0 ) ERROR_exit("the table needs a Subj column when it has a '%s' column",o->run_col) ;
   if( d->have_run ){
     tab_cols[0]=tab->cname[tab->icol_subj] ; tab_cols[1]=o->run_col ;
     d->dx=THD_datatable_index_columns(tab,2,tab_cols,NULL,NULL) ;
     d->nsub=d->dx->nlevel[0] ; d->nrun=d->dx->nlevel[1] ;
   } else {
     if( tab->icol_subj>=0 ){
       tab_cols[0]=tab->cname[tab->icol_subj] ;
       d->dx=THD_datatable_index_columns(tab,1,tab_cols,NULL,NULL) ;   /* duplicate subjects are fatal */
       d->nsub=d->dx->nlevel[0] ;
     } else d->nsub=tab->nrow ;
     d->nrun=1 ;
   }
   if( d->nsub<min_subjects ){
     if( min_subjects==2 ) ERROR_exit("SRM needs at least two subjects; the table has %d",d->nsub) ;
     ERROR_exit("the table has %d subject%s; need at least %d",d->nsub,d->nsub==1?"":"s",min_subjects) ;
   }

   d->sb=(srm_subj *)calloc((size_t)d->nsub,sizeof(srm_subj)) ;
   d->icol_xmat=THD_datatable_column(tab,"Xmat") ;
   if( (d->icol_xmat>=0 || o->concat!=NULL) && d->have_run )
     ERROR_exit("%s splits each subject's one input into runs; the table also has a '%s' column. Use one or the other",
                d->icol_xmat>=0 ? "The Xmat column" : "-concat",o->run_col) ;
   if( d->icol_xmat>=0 && o->concat!=NULL ) ERROR_exit("give the runs once: -concat, or the Xmat column, not both") ;
   d->icol_mask=THD_datatable_column(tab,"Mask") ;
   d->icol_group=THD_datatable_column(tab,o->group_col) ;
   for( ss=0 ; ss<d->nsub ; ss++ ){
     srm_subj *s=d->sb+ss ;
     s->nrun=d->nrun ; s->row=(int *)calloc((size_t)d->nrun,sizeof(int)) ;
     s->label=strdup(d->dx ? d->dx->level[0][ss] : tab->subj[ss]) ;
     for( rr=0 ; rr<d->nrun ; rr++ )
       s->row[rr]=d->have_run ? d->dx->row_of[(size_t)ss*d->nrun+rr] : (d->dx ? d->dx->row_of[ss] : ss) ;
     if( d->icol_mask>=0 ){                             /* one Mask per subject, on all their rows */
       const char *mn=DT_CELL(tab,s->row[0],d->icol_mask) ;
       for( rr=1 ; rr<d->nrun ; rr++ )
         if( strcmp(DT_CELL(tab,s->row[rr],d->icol_mask),mn)!=0 )
           ERROR_exit("Subj %s has different Mask entries on different rows",s->label) ;
       if( mn[0]!='\0' && strcmp(mn,"-")!=0 && strcmp(mn,"NA")!=0 ) s->mask_name=strdup(mn) ;
     }
     if( d->icol_xmat>=0 ){
       const char *xn=DT_CELL(tab,s->row[0],d->icol_xmat) ;
       if( xn[0]!='\0' && strcmp(xn,"-")!=0 && strcmp(xn,"NA")!=0 ) s->xmat=strdup(xn) ;
     }
     if( d->icol_group>=0 ){
       const char *gn=DT_CELL(tab,s->row[0],d->icol_group) ;
       for( rr=1 ; rr<d->nrun ; rr++ )
         if( strcmp(DT_CELL(tab,s->row[rr],d->icol_group),gn)!=0 )
           ERROR_exit("Subj %s has different %s entries on different rows",s->label,o->group_col) ;
       for( jj=0 ; jj<d->ngroup ; jj++ ) if( strcmp(d->glabel[jj],gn)==0 ) break ;
       if( jj==d->ngroup ){ d->glabel=(char **)realloc(d->glabel,sizeof(char *)*(size_t)(d->ngroup+1)) ; d->glabel[d->ngroup++]=strdup(gn) ; }
       s->group=jj ;
     }
   }
   if( d->ngroup==0 ){ d->glabel=(char **)malloc(sizeof(char *)) ; d->glabel[0]=strdup("All") ; d->ngroup=1 ; }
}

/*============================================================================
  Several runs in one input (afni_proc.py's errts): their starts come from -concat, as in
  3dDeconvolve and 3dTproject, or from the subject's X.xmat.1D (RunStart, NRowFull, GoodList).
==============================================================================*/

static void srm_int_list( char *dst , size_t nd , const int *a , int n ) ;

/*! -concat: run starts, as a 1D list ('1D: 0 300 600' or a file); *n receives the count. */
static int * srm_read_concat( const char *spec , int *n )
{
   MRI_IMAGE *cim=mri_read_1D((char *)spec) ; float *car ; int *st , ii ;
   if( cim==NULL ) ERROR_exit("cannot read -concat '%s'",spec) ;
   *n=cim->nvox ; car=MRI_FLOAT_PTR(cim) ;
   st=(int *)malloc(sizeof(int)*(size_t)(*n)) ;
   for( ii=0 ; ii<*n ; ii++ ){
     if( !isfinite(car[ii]) || car[ii]!=rintf(car[ii]) || car[ii]<0.0f )
       ERROR_exit("-concat entry %d is %g; run starts are time-point indices (0, then each run's first time point)",ii+1,car[ii]) ;
     st[ii]=(int)rintf(car[ii]) ;
   }
   mri_free(cim) ;
   if( st[0]!=0 ) ERROR_exit("-concat must start at 0 (it lists where each run starts), not %d",st[0]) ;
   for( ii=1 ; ii<*n ; ii++ ) if( st[ii]<=st[ii-1] )
     ERROR_exit("-concat entries must increase: entry %d (%d) is not after entry %d (%d)",ii+1,st[ii],ii,st[ii-1]) ;
   return st ;
}

/*! A subject's X.xmat.1D: run starts (RunStart; one run if absent), the full length (NRowFull),
    and the uncensored time points (GoodList).  *nfull, *nrun, *start are set; s->good/ngood too. */
static int * srm_read_xmat( srm_subj *s , int *nrun , int *nfull )
{
   NI_element *nel=(NI_element *)NI_read_element_fromfile(s->xmat) ; char *a ; int *st , ii ;
   if( nel==NULL || nel->type!=NI_ELEMENT_TYPE )
     ERROR_exit("Subj %s: cannot read '%s' as a design matrix (an X.xmat.1D from 3dDeconvolve or afni_proc.py)",s->label,s->xmat) ;
   a=NI_get_attribute_nocase(nel,"NRowFull") ;
   if( a==NULL ) ERROR_exit("Subj %s: '%s' has no NRowFull attribute; is it an X.xmat.1D?",s->label,s->xmat) ;
   *nfull=(int)strtod(a,NULL) ;
   a=NI_get_attribute_nocase(nel,"RunStart") ;
   if( a!=NULL ){
     NI_int_array *ia=NI_decode_int_list(a,";,") ;
     if( ia==NULL || ia->num<1 ) ERROR_exit("Subj %s: RunStart in '%s' is badly formatted",s->label,s->xmat) ;
     *nrun=ia->num ; st=(int *)malloc(sizeof(int)*(size_t)ia->num) ; memcpy(st,ia->ar,sizeof(int)*(size_t)ia->num) ;
     NI_delete_int_array(ia) ;
   } else { *nrun=1 ; st=(int *)calloc(1,sizeof(int)) ; }
   for( ii=0 ; ii<*nrun ; ii++ ) if( st[ii]<0 || st[ii]>=*nfull || (ii>0 && st[ii]<=st[ii-1]) || (ii==0 && st[ii]!=0) )
     ERROR_exit("Subj %s: RunStart in '%s' has bad values",s->label,s->xmat) ;
   a=NI_get_attribute_nocase(nel,"GoodList") ;
   if( a!=NULL ){
     NI_int_array *ia=NI_decode_int_list(a,";,") ;
     if( ia==NULL ) ERROR_exit("Subj %s: GoodList in '%s' is badly formatted",s->label,s->xmat) ;
     for( ii=0 ; ii<ia->num ; ii++ ) if( ia->ar[ii]<0 || ia->ar[ii]>=*nfull || (ii>0 && ia->ar[ii]<=ia->ar[ii-1]) )
       ERROR_exit("Subj %s: GoodList in '%s' has bad values",s->label,s->xmat) ;
     if( ia->num<*nfull ){
       s->ngood=ia->num ; s->good=(int *)malloc(sizeof(int)*(size_t)ia->num) ; memcpy(s->good,ia->ar,sizeof(int)*(size_t)ia->num) ;
     }
     NI_delete_int_array(ia) ;
   }
   NI_free_element(nel) ;
   return st ;
}

/*! After srm_data_shape has read each subject's one input: the runs inside it, the same for everyone. */
static void srm_data_concat_runs( const srm_opts *o , srm_data *d )
{
   int ss , rr , nr=0 , nt=d->runlen[0] , *st=NULL , ncens=0 ;
   if( o->concat!=NULL ) st=srm_read_concat(o->concat,&nr) ;
   for( ss=0 ; ss<d->nsub ; ss++ ){
     srm_subj *s=d->sb+ss ; int n2 , *st2 , nfull ;
     if( s->xmat==NULL ) continue ;
     st2=srm_read_xmat(s,&n2,&nfull) ;
     if( nfull!=nt )
       ERROR_exit("Subj %s: '%s' describes %d time points (NRowFull), but the input has %d",s->label,s->xmat,nfull,nt) ;
     if( st==NULL ){ st=st2 ; nr=n2 ; continue ; }
     for( rr=0 ; rr<nr && n2==nr ; rr++ ) if( st2[rr]!=st[rr] ) break ;
     if( n2!=nr || rr<nr ){
       char a1[256] , a2[256] ;
       srm_int_list(a1,sizeof(a1),st,nr) ; srm_int_list(a2,sizeof(a2),st2,n2) ;
       ERROR_exit("Subj %s: the runs in '%s' start at %s, not at %s as for %s; on the time axis every subject needs the "
                  "same runs",s->label,s->xmat,a2,a1,o->concat!=NULL?"-concat":"the first subject with an Xmat") ;
     }
     free(st2) ;
   }
   if( st==NULL ) return ;
   if( st[nr-1]>=nt ) ERROR_exit("the last run starts at time point %d, but the inputs have %d",st[nr-1],nt) ;
   d->concat=1 ; d->cstart=st ; d->nrun=nr ;
   free(d->runlen) ; d->runlen=(int *)malloc(sizeof(int)*(size_t)nr) ;
   for( rr=0 ; rr<nr ; rr++ ) d->runlen[rr]=((rr<nr-1) ? st[rr+1] : nt)-st[rr] ;
   if( !o->quiet ){
     char ls[256] ; int pos=0 ;
     for( rr=0 ; rr<nr && pos<220 ; rr++ ) pos+=snprintf(ls+pos,sizeof(ls)-pos,"%s%d",rr?",":"",d->runlen[rr]) ;
     INFO_message("each input holds %d runs (%s, from %s): lengths %s",nr,o->concat!=NULL?"-concat":"the Xmat column",
                  o->concat!=NULL?o->concat:"RunStart",ls) ;
   }
   for( ss=0 ; ss<d->nsub ; ss++ ) if( d->sb[ss].good!=NULL ) ncens++ ;
   if( ncens>0 )
     WARNING_message("%d subject%s' Xmat list%s censored time points (GoodList); censoring from the Xmat is not applied in this "
                     "build: those time points are used as data. Give -censor, or -zcensor if they are zero in the input",
                     ncens,ncens==1?"":"s",ncens==1?"s":"") ;
}

/*! Open every input (header only for datasets): lengths, grids, masks.  With need_mask == 0
    (applying a saved W) no mask is read; the features come from W. */
static void srm_data_shape( const srm_opts *o , srm_data *d , int need_mask , int Kmin )
{
   int ss , rr , ii , jj ; THD_datatable *tab=d->tab ;
   if( need_mask && o->maskname!=NULL ) d->gmask=srm_mask_get(&d->mcache,&d->nmcache,o->maskname) ;
   d->runlen=(int *)calloc((size_t)d->nrun,sizeof(int)) ;
   for( ss=0 ; ss<d->nsub ; ss++ ){
     srm_subj *s=d->sb+ss ;
     for( rr=0 ; rr<d->nrun ; rr++ ){
       const char *fn=DT_CELL(tab,s->row[rr],tab->icol_input) ; int nt , is1 ;
       is1=srm_is1D(fn) ;
       if( rr==0 ) s->is1D=is1 ;
       else if( is1!=s->is1D ) ERROR_exit("Subj %s mixes .1D files and datasets across runs",s->label) ;
       if( is1 ){
         MRI_IMAGE *im=mri_read_1D((char *)fn) ;
         if( im==NULL ) ERROR_exit("cannot read 1D file '%s' (Subj %s)",fn,s->label) ;
         nt=im->nx ;
         if( rr==0 ){
           s->nvox_mask=im->ny ;
           if( need_mask ){ s->vox=(int *)malloc(sizeof(int)*(size_t)im->ny) ; for( jj=0 ; jj<im->ny ; jj++ ) s->vox[jj]=jj ; }
         }
         else if( im->ny!=s->nvox_mask ) ERROR_exit("Subj %s: '%s' has %d columns; the first run had %d",s->label,fn,im->ny,s->nvox_mask) ;
         mri_free(im) ;
       } else {
         THD_3dim_dataset *ds=THD_open_dataset((char *)fn) ;
         CHECK_OPEN_ERROR(ds,fn) ;
         nt=DSET_NVALS(ds) ;
         if( rr==0 ){
           s->geom=ds ;
           if( need_mask ){
             mask_entry *me ;
             if( s->mask_name!=NULL ) me=srm_mask_get(&d->mcache,&d->nmcache,s->mask_name) ;
             else if( d->gmask!=NULL ) me=d->gmask ;
             else ERROR_exit("Subj %s: datasets need a mask: give -mask or a Mask column in the table",s->label) ;
             if( !EQUIV_GRIDS(ds,me->ds) )
               ERROR_exit("Subj %s: '%s' is not on the grid of mask '%s'",s->label,fn,me->name) ;
             s->nvox_mask=me->nmask ;
             s->vox=(int *)malloc(sizeof(int)*(size_t)me->nmask) ;
             for( jj=0,ii=0 ; jj<me->nvox ; jj++ ) if( me->mask[jj] ) s->vox[ii++]=jj ;
           }
         } else {
           if( !EQUIV_GRIDS(ds,s->geom) ) ERROR_exit("Subj %s: run %s is off the grid of run 1",s->label,d->have_run?d->dx->level[1][rr]:"?") ;
           DSET_delete(ds) ;
         }
       }
       if( ss==0 ) d->runlen[rr]=nt ;
       else if( !o->free_runs && nt!=d->runlen[rr] )
         ERROR_exit("Subj %s, run %s has %d time points; the first subject's has %d",s->label,d->have_run?d->dx->level[1][rr]:"1",nt,d->runlen[rr]) ;
     }
     if( need_mask && s->nvox_mask<Kmin ) ERROR_exit("Subj %s has %d features; -nfeatures is %d",s->label,s->nvox_mask,Kmin) ;
   }
   if( !d->have_run && (o->concat!=NULL || d->icol_xmat>=0) ){
     if( o->free_runs ) ERROR_exit("an Xmat column with -apply_W is not available in this build yet") ;
     srm_data_concat_runs(o,d) ;
   }
   d->Tin=0 ;
   for( rr=0 ; rr<d->nrun ; rr++ ){ d->Tin+=d->runlen[rr] ; if( d->runlen[rr]<2 ) ERROR_exit("run %d has %d time points; need at least 2",rr+1,d->runlen[rr]) ; }
   d->c0in=(int *)calloc((size_t)d->nrun,sizeof(int)) ;
   for( rr=1 ; rr<d->nrun ; rr++ ) d->c0in[rr]=d->c0in[rr-1]+d->runlen[rr-1] ;
   d->runname=(char **)malloc(sizeof(char *)*(size_t)d->nrun) ;
   for( rr=0 ; rr<d->nrun ; rr++ ){ char nm[24] ; snprintf(nm,sizeof(nm),"%d",rr+1) ; d->runname[rr]=strdup(d->have_run ? d->dx->level[1][rr] : nm) ; }
   if( !o->quiet ){
     INFO_message("%d subjects, %d group%s, %d run%s, %d time points in all",d->nsub,d->ngroup,d->ngroup>1?"s":"",d->nrun,d->nrun>1?"s":"",d->Tin) ;
     if( d->nrun==1 && !d->have_run )
       ININFO_message("one run: if each input holds several runs (an afni_proc.py errts, say), give -concat or an Xmat column, "
                      "so that z-scoring and -cv respect the run boundaries") ;
     if( d->ngroup>1 ){ for( jj=0 ; jj<d->ngroup ; jj++ ){ int cnt=0 ; for( ss=0 ; ss<d->nsub ; ss++ ) cnt+=(d->sb[ss].group==jj) ;
         ININFO_message("group %s: %d subjects",d->glabel[jj],cnt) ; } }
   }
}

/*! The retained time points: from -censor, or from a saved model's mask (given_keep). */
static void srm_data_censor( const srm_opts *o , srm_data *d , const unsigned char *given_keep )
{
   int ii ;
   d->keep=(byte *)malloc((size_t)d->Tin) ; memset(d->keep,1,(size_t)d->Tin) ;
   if( given_keep!=NULL ) memcpy(d->keep,given_keep,(size_t)d->Tin) ;
   else if( o->censor_name!=NULL ){
     MRI_IMAGE *cim=mri_read_1D(o->censor_name) , *fim ; float *car ;
     if( cim==NULL ) ERROR_exit("cannot read -censor file '%s'",o->censor_name) ;
     fim=mri_to_float(cim) ; if( fim!=cim ) mri_free(cim) ;
     if( fim->nvox!=d->Tin ) ERROR_exit("-censor has %d values; the data have %d time points",fim->nvox,d->Tin) ;
     car=MRI_FLOAT_PTR(fim) ;
     for( ii=0 ; ii<d->Tin ; ii++ ){
       if( !isfinite(car[ii]) || (car[ii]!=0.0f && car[ii]!=1.0f) )
         ERROR_exit("-censor value %d is %.7g; expected 0 or 1",ii+1,car[ii]) ;
       d->keep[ii]=(byte)(car[ii]!=0.0f) ;
     }
     mri_free(fim) ;
   }
   if( o->zcensor && given_keep==NULL ){
     d->zero_any=(byte *)calloc((size_t)d->Tin,1) ; d->nzero=(int *)calloc((size_t)d->nsub,sizeof(int)) ;
   }
}

/*! Read every subject's masked time series, all runs, onto the input timeline. */
static void srm_data_read( const srm_opts *o , srm_data *d )
{
   srm_prog bar ; int ss , rr , ii ; THD_datatable *tab=d->tab ; int Tin=d->Tin ;
   srm_prog_begin(&bar,o->progress,o->quiet,"reading","subjects",d->nsub) ;
   for( ss=0 ; ss<d->nsub ; ss++ ){
     srm_subj *s=d->sb+ss ;
     byte *zf = d->zero_any ? (byte *)calloc((size_t)Tin,1) : NULL ;
     s->x=(float *)malloc(sizeof(float)*(size_t)s->nvox_mask*(size_t)Tin) ;
     if( s->x==NULL ) ERROR_exit("cannot allocate %.2f GiB for Subj %s",4.0*s->nvox_mask*(double)Tin/GIB,s->label) ;
     for( rr=0 ; rr<(d->concat ? 1 : d->nrun) ; rr++ ){    /* -concat/Xmat: one input holds every run */
       const char *fn=DT_CELL(tab,s->row[rr],tab->icol_input) ;
       if( s->is1D ){
         MRI_IMAGE *im=mri_read_1D((char *)fn) , *fm ; float *ar ; int nt , f , t ;
         fm=mri_to_float(im) ; ar=MRI_FLOAT_PTR(fm) ; nt=fm->nx ;
         for( f=0 ; f<s->nvox_mask ; f++ ) for( t=0 ; t<nt ; t++ ){
           float v=ar[(size_t)f*nt+t] ;
           if( !isfinite(v) ) ERROR_exit("Subj %s: '%s' has a non-finite value at row %d, column %d",s->label,fn,t+1,f+1) ;
           s->x[(size_t)f*Tin+d->c0in[rr]+t]=v ;
         }
         if( zf!=NULL ) for( t=0 ; t<nt ; t++ ){ int nz=0 ; for( f=0 ; f<s->nvox_mask && !nz ; f++ ) nz=(ar[(size_t)f*nt+t]!=0.0f) ; zf[d->c0in[rr]+t]=(byte)!nz ; }
         if( fm!=im ) mri_free(fm) ;
         mri_free(im) ;
       } else {
         THD_3dim_dataset *ds=(rr==0) ? s->geom : THD_open_dataset((char *)fn) ; int badt ;
         if( ds==NULL ) ERROR_exit("cannot open '%s'",fn) ;
         if( srm_read_dset(ds,s->vox,s->nvox_mask,s->x,Tin,d->c0in[rr],zf,&badt)>0 )
           ERROR_exit("Subj %s: '%s' has non-finite values in the mask (first at volume %d)",s->label,fn,badt+1) ;
         if( rr>0 ) DSET_delete(ds) ;
       }
     }
     if( zf!=NULL ){
       for( ii=0 ; ii<Tin ; ii++ ) if( zf[ii] ){ d->nzero[ss]++ ; d->zero_any[ii]=1 ; }
       free(zf) ;
     }
     srm_prog_set(&bar,ss+1,NULL) ;
   }
}

/*! Apply -zcensor and the censor mask: the retained timeline, its run lengths, packed data. */
static void srm_data_timeline( const srm_opts *o , srm_data *d , int Kmin )
{
   int ii , rr , ss ; int Tin=d->Tin ;
   if( o->zcensor && d->zero_any!=NULL ){
     int nz=0 ;
     for( ii=0 ; ii<Tin ; ii++ ) if( d->zero_any[ii] ){ nz++ ; d->keep[ii]=0 ; }
     if( !o->quiet ) INFO_message("-zcensor found %d all-zero time points; they are dropped for everyone",nz) ;
     for( ss=0 ; ss<d->nsub ; ss++ ){
       float pct=100.0f*d->nzero[ss]/(float)Tin ;
       if( pct>o->zcensor_warn )
         WARNING_message("Subj %s has %d of %d time points (%.1f%%) that are all zero in the mask; threshold is %.1f%%",
                         d->sb[ss].label,d->nzero[ss],Tin,pct,o->zcensor_warn) ;
     }
   }
   d->T=0 ;
   d->seglen=(int *)calloc((size_t)d->nrun,sizeof(int)) ; d->time_index=(int *)malloc(sizeof(int)*(size_t)Tin) ;
   for( rr=0 ; rr<d->nrun ; rr++ ) for( ii=d->c0in[rr] ; ii<d->c0in[rr]+d->runlen[rr] ; ii++ ) if( d->keep[ii] ){ d->time_index[d->T++]=ii ; d->seglen[rr]++ ; }
   for( rr=0 ; rr<d->nrun ; rr++ ) if( d->seglen[rr]<3 )
     ERROR_exit("run %d keeps only %d of %d time points after censoring; need at least 3",rr+1,d->seglen[rr],d->runlen[rr]) ;
   if( d->T<Kmin ) ERROR_exit("%d time points are left after censoring; -nfeatures is %d",d->T,Kmin) ;
   if( d->T<Tin && !o->quiet ) INFO_message("%d of %d time points retained after censoring",d->T,Tin) ;
   if( d->T<Tin ){                                    /* pack the columns we keep */
     for( ss=0 ; ss<d->nsub ; ss++ ){
       srm_subj *s=d->sb+ss ;
       float *packed=(float *)malloc(sizeof(float)*(size_t)s->nvox_mask*(size_t)d->T) ; int f , t ;
       if( packed==NULL ) ERROR_exit("cannot allocate the censored data for Subj %s",s->label) ;
       for( f=0 ; f<s->nvox_mask ; f++ ){
         const float *src=s->x+(size_t)f*Tin ; float *dst=packed+(size_t)f*d->T ;
         for( t=0 ; t<d->T ; t++ ) dst[t]=src[d->time_index[t]] ;
       }
       free(s->x) ; s->x=packed ;
     }
   }
}

/*============================================================================
  -axis item: one table row per Subj x Item (x Run), each a single volume or a one-row .1D.
  Columns are keyed by Item, or Item:n for the n-th occurrence (in run order) of an item that a
  subject has more than once; every subject must have the same keys.  A subject's runs need not
  match anyone else's (randomized designs): they only group that subject's items for z-scoring
  and number the repeats.
==============================================================================*/

typedef struct { int row , run , pos ; const char *item ; } srm_irow ;

static int srm_irow_by_run( const void *a , const void *b )      /* run order, then table order */
{
   const srm_irow *x=(const srm_irow *)a , *y=(const srm_irow *)b ;
   if( x->run!=y->run ) return (x->run<y->run) ? -1 : 1 ;
   return (x->row<y->row) ? -1 : (x->row>y->row) ;
}

static int srm_irow_by_item( const void *a , const void *b )     /* item, then position in run order */
{
   const srm_irow *x=(const srm_irow *)a , *y=(const srm_irow *)b ; int c=strcmp(x->item,y->item) ;
   if( c!=0 ) return c ;
   return (x->pos<y->pos) ? -1 : (x->pos>y->pos) ;
}

typedef struct { const char *key ; int idx ; } srm_ikey ;

static int srm_ikey_cmp( const void *a , const void *b )
{
   return strcmp(((const srm_ikey *)a)->key,((const srm_ikey *)b)->key) ;
}

/*! Subjects, runs, groups, masks, and the item keys of every subject, matched to the first
    subject's order. */
static void srm_data_items( const srm_opts *o , THD_datatable *tab , srm_data *d , int min_subjects )
{
   int icol_item , icol_rep , r , ss , jj , nsub=0 , *subj_of , *run_of , nrep=0 , minrun=-1 ; char **slab=NULL ;
   srm_ikey *canon=NULL ;

   d->tab=tab ; d->axis_item=1 ;
   if( tab->icol_input<0 ) ERROR_exit("the data table needs an InputFile column") ;
   if( tab->icol_subj<0 ) ERROR_exit("-axis item needs a Subj column") ;
   icol_item=THD_datatable_column(tab,o->item_col) ;
   if( icol_item<0 ) ERROR_exit("-axis item needs an item column: the table has no '%s' column (see -item_column)",o->item_col) ;
   d->icol_run=THD_datatable_column(tab,o->run_col) ; d->have_run=(d->icol_run>=0) ;
   d->icol_mask=THD_datatable_column(tab,"Mask") ; d->icol_group=THD_datatable_column(tab,o->group_col) ;
   icol_rep=THD_datatable_column(tab,o->repeat_col) ;
   if( o->zscore==SRM_ZSCORE_RUN && !d->have_run )
     ERROR_exit("-zscore run works over the items of each run: the table has no '%s' column (the default on the item axis is "
                "-zscore all)",o->run_col) ;

   /* subjects and runs, in order of first appearance */
   subj_of=(int *)malloc(sizeof(int)*(size_t)tab->nrow) ; run_of=(int *)calloc((size_t)tab->nrow,sizeof(int)) ;
   for( r=0 ; r<tab->nrow ; r++ ){
     const char *sl=DT_CELL(tab,r,tab->icol_subj) , *il=DT_CELL(tab,r,icol_item) , *cc ;
     if( il[0]=='\0' ) ERROR_exit("row %d of the table has an empty %s",r+1,o->item_col) ;
     for( cc=il ; *cc!='\0' ; cc++ )
       if( isspace((unsigned char)*cc) || *cc==':' ) ERROR_exit("%s '%s' (row %d): item labels may not contain ':' or spaces",o->item_col,il,r+1) ;
     if( icol_rep>=0 ) for( cc=DT_CELL(tab,r,icol_rep) ; *cc!='\0' ; cc++ )
       if( isspace((unsigned char)*cc) || *cc==':' ) ERROR_exit("%s '%s' (row %d): repeat labels may not contain ':' or spaces",o->repeat_col,DT_CELL(tab,r,icol_rep),r+1) ;
     for( ss=0 ; ss<nsub ; ss++ ) if( strcmp(slab[ss],sl)==0 ) break ;
     if( ss==nsub ){ slab=(char **)realloc(slab,sizeof(char *)*(size_t)(nsub+1)) ; slab[nsub++]=strdup(sl) ; }
     subj_of[r]=ss ;
     if( d->have_run ){
       const char *rl=DT_CELL(tab,r,d->icol_run) ;
       for( jj=0 ; jj<d->nrun ; jj++ ) if( strcmp(d->runname[jj],rl)==0 ) break ;
       if( jj==d->nrun ){ d->runname=(char **)realloc(d->runname,sizeof(char *)*(size_t)(d->nrun+1)) ; d->runname[d->nrun++]=strdup(rl) ; }
       run_of[r]=jj ;
     }
   }
   if( !d->have_run ){ d->nrun=1 ; d->runname=(char **)malloc(sizeof(char *)) ; d->runname[0]=strdup("1") ; }
   d->nsub=nsub ;
   if( nsub<min_subjects ){
     if( min_subjects==2 ) ERROR_exit("SRM needs at least two subjects; the table has %d",nsub) ;
     ERROR_exit("the table has %d subject%s; need at least %d",nsub,nsub==1?"":"s",min_subjects) ;
   }

   d->sb=(srm_subj *)calloc((size_t)nsub,sizeof(srm_subj)) ;
   for( ss=0 ; ss<nsub ; ss++ ){
     srm_subj *s=d->sb+ss ; srm_irow *ir ; int n=0 , k , k0 , *seen=NULL , *perrun ;
     char **key ;
     s->label=slab[ss] ;
     for( r=0 ; r<tab->nrow ; r++ ) n+=(subj_of[r]==ss) ;
     ir=(srm_irow *)malloc(sizeof(srm_irow)*(size_t)n) ; key=(char **)calloc((size_t)n,sizeof(char *)) ;
     for( r=0,k=0 ; r<tab->nrow ; r++ ) if( subj_of[r]==ss ){ ir[k].row=r ; ir[k].run=run_of[r] ; ir[k].item=DT_CELL(tab,r,icol_item) ; k++ ; }
     qsort(ir,(size_t)n,sizeof(srm_irow),srm_irow_by_run) ;
     for( k=0 ; k<n ; k++ ) ir[k].pos=k ;

     /* repeats: an item's occurrences are Item:1, Item:2, ... in run order, then table order (or
        Item:R, R from the repeat column) */
     qsort(ir,(size_t)n,sizeof(srm_irow),srm_irow_by_item) ;
     for( k0=0 ; k0<n ; ){
       int k1=k0+1 , q , q2 ;
       while( k1<n && strcmp(ir[k1].item,ir[k0].item)==0 ) k1++ ;
       if( k1-k0>1 && icol_rep>=0 )
         for( q=k0 ; q<k1 ; q++ ) for( q2=q+1 ; q2<k1 ; q2++ )
           if( strcmp(DT_CELL(tab,ir[q].row,icol_rep),DT_CELL(tab,ir[q2].row,icol_rep))==0 )
             ERROR_exit("Subj %s has %s %s with %s %s twice",s->label,o->item_col,ir[q].item,o->repeat_col,DT_CELL(tab,ir[q].row,icol_rep)) ;
       for( q=k0 ; q<k1 ; q++ ){
         char kb[THD_MAX_NAME] ;
         if( k1-k0>1 && icol_rep>=0 ) snprintf(kb,sizeof(kb),"%s:%s",ir[q].item,DT_CELL(tab,ir[q].row,icol_rep)) ;
         else if( k1-k0>1 )           snprintf(kb,sizeof(kb),"%s:%d",ir[q].item,q-k0+1) ;
         else                         snprintf(kb,sizeof(kb),"%s",ir[q].item) ;
         key[ir[q].pos]=strdup(kb) ;
       }
       k0=k1 ;
     }
     qsort(ir,(size_t)n,sizeof(srm_irow),srm_irow_by_run) ;        /* back to run order: ir[k] has key[k] */

     if( ss==0 ){                                                    /* the first subject sets the columns */
       d->nitem=n ; d->item_key=key ; key=NULL ;
       canon=(srm_ikey *)malloc(sizeof(srm_ikey)*(size_t)n) ;
       for( k=0 ; k<n ; k++ ){ canon[k].key=d->item_key[k] ; canon[k].idx=k ; if( strchr(d->item_key[k],':')!=NULL ) nrep++ ; }
       qsort(canon,(size_t)n,sizeof(srm_ikey),srm_ikey_cmp) ;
     }
     s->irow=(int *)malloc(sizeof(int)*(size_t)d->nitem) ; s->irun=(int *)malloc(sizeof(int)*(size_t)d->nitem) ;
     seen=(int *)calloc((size_t)d->nitem,sizeof(int)) ;
     for( k=0 ; k<n ; k++ ){
       srm_ikey want , *hit ; int jc ;
       want.key=(ss==0) ? d->item_key[k] : key[k] ;
       hit=(srm_ikey *)bsearch(&want,canon,(size_t)d->nitem,sizeof(srm_ikey),srm_ikey_cmp) ;
       if( hit==NULL ){
         int mine=0 , theirs=0 , q ; size_t L=strlen(ir[k].item) ;
         for( q=0 ; q<n ; q++ ) mine+=(strcmp(ir[q].item,ir[k].item)==0) ;
         for( q=0 ; q<d->nitem ; q++ )
           theirs+=(strncmp(d->item_key[q],ir[k].item,L)==0 && (d->item_key[q][L]=='\0' || d->item_key[q][L]==':')) ;
         if( theirs==0 ) ERROR_exit("Subj %s has %s %s, which Subj %s does not have",s->label,o->item_col,ir[k].item,d->sb[0].label) ;
         ERROR_exit("Subj %s has %s %s %d time%s, Subj %s has it %d time%s: a repeated item must recur the same "
                    "number of times for every subject",s->label,o->item_col,ir[k].item,mine,mine==1?"":"s",
                    d->sb[0].label,theirs,theirs==1?"":"s") ;
       }
       jc=hit->idx ; seen[jc]=1 ; s->irow[jc]=ir[k].row ; s->irun[jc]=ir[k].run ;
     }
     for( k=0 ; k<d->nitem ; k++ ) if( !seen[k] )
       ERROR_exit("Subj %s is missing %s %s, which Subj %s has%s",s->label,o->item_col,d->item_key[k],d->sb[0].label,
                  strchr(d->item_key[k],':')!=NULL ? " (an item numbered :n recurs; every subject needs the same number of repeats)" : "") ;
     if( key!=NULL ){ for( k=0 ; k<n ; k++ ) free(key[k]) ; free(key) ; }
     free(seen) ; free(ir) ;

     /* the subject's own runs: z-scoring works over the items of each one */
     perrun=(int *)calloc((size_t)d->nrun,sizeof(int)) ;
     for( k=0 ; k<d->nitem ; k++ ) perrun[s->irun[k]]++ ;
     for( k=0 ; k<d->nrun ; k++ ) if( perrun[k]>0 ){
       if( o->zscore==SRM_ZSCORE_RUN && perrun[k]<3 )
         ERROR_exit("Subj %s has %d item%s in run %s: z-scoring within a run needs at least 3 (or use -zscore none)",
                    s->label,perrun[k],perrun[k]==1?"":"s",d->runname[k]) ;
       if( minrun<0 || perrun[k]<minrun ) minrun=perrun[k] ;
     }
     free(perrun) ;

     /* one Mask and one Group per subject */
     s->nrun=1 ; s->row=(int *)malloc(sizeof(int)) ; s->row[0]=s->irow[0] ;
     if( d->icol_mask>=0 ){
       const char *mn=DT_CELL(tab,s->irow[0],d->icol_mask) ;
       for( k=1 ; k<d->nitem ; k++ ) if( strcmp(DT_CELL(tab,s->irow[k],d->icol_mask),mn)!=0 )
         ERROR_exit("Subj %s has different Mask entries on different rows",s->label) ;
       if( mn[0]!='\0' && strcmp(mn,"-")!=0 && strcmp(mn,"NA")!=0 ) s->mask_name=strdup(mn) ;
     }
     if( d->icol_group>=0 ){
       const char *gn=DT_CELL(tab,s->irow[0],d->icol_group) ;
       for( k=1 ; k<d->nitem ; k++ ) if( strcmp(DT_CELL(tab,s->irow[k],d->icol_group),gn)!=0 )
         ERROR_exit("Subj %s has different %s entries on different rows",s->label,o->group_col) ;
       for( jj=0 ; jj<d->ngroup ; jj++ ) if( strcmp(d->glabel[jj],gn)==0 ) break ;
       if( jj==d->ngroup ){ d->glabel=(char **)realloc(d->glabel,sizeof(char *)*(size_t)(d->ngroup+1)) ; d->glabel[d->ngroup++]=strdup(gn) ; }
       s->group=jj ;
     }
   }
   if( d->ngroup==0 ){ d->glabel=(char **)malloc(sizeof(char *)) ; d->glabel[0]=strdup("All") ; d->ngroup=1 ; }
   if( !o->quiet ){
     INFO_message("%d subjects, %d group%s, %d item columns each%s",nsub,d->ngroup,d->ngroup>1?"s":"",d->nitem,
                  nrep>0 ? " (counting each repeat of an item as its own column)" : "") ;
     if( nrep>0 ) ININFO_message("%d columns are repeats, keyed Item:%s",nrep,
                                 icol_rep>=0 ? "R by the repeat column" : "n in run order, then table order") ;
     if( d->have_run ) ININFO_message("%d run%s; %s",d->nrun,d->nrun>1?"s":"",
                                      o->zscore==SRM_ZSCORE_RUN ? "each subject's items are z-scored within their own runs"
                                                                : "they order repeats (and can be the -cv folds)") ;
     ININFO_message("-zscore %s",o->zscore==SRM_ZSCORE_ALL ? "all: each feature over all of the subject's items"
                                 : o->zscore==SRM_ZSCORE_RUN ? "run: each feature over the items of each of the subject's runs"
                                 : "none: the betas as given (centered by the fit)") ;
     if( d->ngroup>1 ){ for( jj=0 ; jj<d->ngroup ; jj++ ){ int cnt=0 ; for( ss=0 ; ss<nsub ; ss++ ) cnt+=(d->sb[ss].group==jj) ;
         ININFO_message("group %s: %d subjects",d->glabel[jj],cnt) ; } }
   }
   if( o->zscore==SRM_ZSCORE_RUN && minrun>=0 && minrun<10 && !o->quiet )
     WARNING_message("a run has only %d items: z-scores over so few items are noisy; consider fewer, larger runs, or -zscore none",minrun) ;
   free(subj_of) ; free(run_of) ; free(slab) ; free(canon) ;
}

/*! Open every item's input once per file (header only): grids, masks, and the sub-brick each
    item selects.  A selector like betas+tlrc'[face#0_Coef]' must pick exactly one volume. */
static void srm_data_shape_items( const srm_opts *o , srm_data *d , int need_mask , int Kmin )
{
   int ss , j , f , ii , n=d->nitem ; THD_datatable *tab=d->tab ;
   if( need_mask && o->maskname!=NULL ) d->gmask=srm_mask_get(&d->mcache,&d->nmcache,o->maskname) ;
   for( ss=0 ; ss<d->nsub ; ss++ ){
     srm_subj *s=d->sb+ss ;
     s->ifile=(int *)malloc(sizeof(int)*(size_t)n) ; s->ibrick=(int *)calloc((size_t)n,sizeof(int)) ;
     for( j=0 ; j<n ; j++ ){
       const char *fn=DT_CELL(tab,s->irow[j],tab->icol_input) ; int is1=srm_is1D(fn) ;
       if( j==0 ) s->is1D=is1 ;
       else if( is1!=s->is1D ) ERROR_exit("Subj %s mixes .1D files and datasets",s->label) ;
       if( is1 ){
         MRI_IMAGE *im=mri_read_1D((char *)fn) ;
         if( im==NULL ) ERROR_exit("cannot read 1D file '%s' (Subj %s, item %s)",fn,s->label,d->item_key[j]) ;
         if( im->nx!=1 )
           ERROR_exit("Subj %s, item %s: '%s' has %d rows; with -axis item a .1D input is one row (the item's features)",
                      s->label,d->item_key[j],fn,im->nx) ;
         if( j==0 ){
           s->nvox_mask=im->ny ;
           if( need_mask ){ s->vox=(int *)malloc(sizeof(int)*(size_t)im->ny) ; for( ii=0 ; ii<im->ny ; ii++ ) s->vox[ii]=ii ; }
         } else if( im->ny!=s->nvox_mask )
           ERROR_exit("Subj %s, item %s: '%s' has %d columns; the first item had %d",s->label,d->item_key[j],fn,im->ny,s->nvox_mask) ;
         mri_free(im) ; s->ifile[j]=-1 ;
       } else {
         char base[THD_MAX_NAME] ; const char *sel=strchr(fn,'[') ; size_t len=sel!=NULL ? (size_t)(sel-fn) : strlen(fn) ;
         if( strchr(fn,'<')!=NULL || strchr(fn,'{')!=NULL )
           ERROR_exit("Subj %s, item %s: '%s': select the item's volume with [..]; <..> and {..} selectors are not used with -axis item",
                      s->label,d->item_key[j],fn) ;
         if( len>=sizeof(base) ) ERROR_exit("input name too long: '%s'",fn) ;
         memcpy(base,fn,len) ; base[len]='\0' ;
         for( f=0 ; f<s->nfile ; f++ ) if( strcmp(s->file[f],base)==0 ) break ;
         if( f==s->nfile ){ s->file=(char **)realloc(s->file,sizeof(char *)*(size_t)(s->nfile+1)) ; s->file[s->nfile++]=strdup(base) ; }
         s->ifile[j]=f ;
       }
     }
     for( f=0 ; f<s->nfile ; f++ ){
       THD_3dim_dataset *ds=THD_open_dataset(s->file[f]) ;
       CHECK_OPEN_ERROR(ds,s->file[f]) ;
       if( f==0 ){
         s->geom=ds ;
         if( need_mask ){
           mask_entry *me ;
           if( s->mask_name!=NULL ) me=srm_mask_get(&d->mcache,&d->nmcache,s->mask_name) ;
           else if( d->gmask!=NULL ) me=d->gmask ;
           else ERROR_exit("Subj %s: datasets need a mask: give -mask or a Mask column in the table",s->label) ;
           if( !EQUIV_GRIDS(ds,me->ds) ) ERROR_exit("Subj %s: '%s' is not on the grid of mask '%s'",s->label,s->file[f],me->name) ;
           s->nvox_mask=me->nmask ;
           s->vox=(int *)malloc(sizeof(int)*(size_t)me->nmask) ;
           for( j=0,ii=0 ; j<me->nvox ; j++ ) if( me->mask[j] ) s->vox[ii++]=j ;
         }
       } else if( !EQUIV_GRIDS(ds,s->geom) )
         ERROR_exit("Subj %s: '%s' is not on the grid of '%s'",s->label,s->file[f],s->file[0]) ;
       for( j=0 ; j<n ; j++ ) if( s->ifile[j]==f ){
         const char *fn=DT_CELL(tab,s->irow[j],tab->icol_input) , *sel=strchr(fn,'[') ;
         if( sel!=NULL ){
           int *iv=MCW_get_thd_intlist(ds,(char *)sel) ;
           if( iv==NULL || iv[0]!=1 )
             ERROR_exit("Subj %s, item %s: '%s' selects %d volumes; each item needs exactly one",
                        s->label,d->item_key[j],fn,iv==NULL ? 0 : iv[0]) ;
           s->ibrick[j]=iv[1] ; free(iv) ;
         } else if( DSET_NVALS(ds)!=1 )
           ERROR_exit("Subj %s, item %s: '%s' has %d volumes; select the item's one, as in %s'[3]'",
                      s->label,d->item_key[j],fn,DSET_NVALS(ds),fn) ;
       }
       if( f>0 ) DSET_delete(ds) ;
     }
     if( need_mask && s->nvox_mask<Kmin ) ERROR_exit("Subj %s has %d features; -nfeatures is %d",s->label,s->nvox_mask,Kmin) ;
   }
   d->Tin=n ;
}

/*! Read every item, one load per file.  A non-finite value, or an item that is zero everywhere
    in the mask (a failed or missing estimate, usually), is an error naming the Subj and item. */
static void srm_data_read_items( const srm_opts *o , srm_data *d )
{
   srm_prog bar ; int ss , j , f , v , n=d->nitem , done=0 ; THD_datatable *tab=d->tab ; char extra[200] ;
   srm_prog_begin(&bar,o->progress,o->quiet,"reading","items",d->nsub*n) ;
   for( ss=0 ; ss<d->nsub ; ss++ ){
     srm_subj *s=d->sb+ss ;
     s->x=(float *)malloc(sizeof(float)*(size_t)s->nvox_mask*(size_t)n) ;
     if( s->x==NULL ) ERROR_exit("cannot allocate %.2f GiB for Subj %s",4.0*s->nvox_mask*(double)n/GIB,s->label) ;
     if( s->is1D ){
       for( j=0 ; j<n ; j++ ){
         const char *fn=DT_CELL(tab,s->irow[j],tab->icol_input) ;
         MRI_IMAGE *im=mri_read_1D((char *)fn) , *fm ; float *ar ; int nz=0 ;
         if( im==NULL ) ERROR_exit("cannot read 1D file '%s' (Subj %s, item %s)",fn,s->label,d->item_key[j]) ;
         fm=mri_to_float(im) ; ar=MRI_FLOAT_PTR(fm) ;
         for( v=0 ; v<s->nvox_mask ; v++ ){
           float a=ar[v] ;
           if( !isfinite(a) ) ERROR_exit("Subj %s, item %s: '%s' has a non-finite value in column %d",s->label,d->item_key[j],fn,v+1) ;
           if( a!=0.0f ) nz=1 ;
           s->x[(size_t)v*n+j]=a ;
         }
         if( fm!=im ) mri_free(fm) ;
         mri_free(im) ;
         if( !nz ) ERROR_exit("Subj %s, item %s: '%s' is all zero (a failed or missing estimate?)",s->label,d->item_key[j],fn) ;
         snprintf(extra,sizeof(extra),"Subj %s (%d/%d), %s",s->label,ss+1,d->nsub,THD_trailname((char *)fn,0)) ;
         srm_prog_set(&bar,++done,extra) ;
       }
     } else for( f=0 ; f<s->nfile ; f++ ){
       THD_3dim_dataset *ds=THD_open_dataset(s->file[f]) ;
       CHECK_OPEN_ERROR(ds,s->file[f]) ;
       DSET_load(ds) ; CHECK_LOAD_ERROR(ds) ;
       for( j=0 ; j<n ; j++ ) if( s->ifile[j]==f ){
         MRI_IMAGE *fim=THD_extract_float_brick(s->ibrick[j],ds) ; float *far ; int nz=0 ;
         if( fim==NULL ) ERROR_exit("cannot read volume %d of '%s' (Subj %s, item %s)",s->ibrick[j],s->file[f],s->label,d->item_key[j]) ;
         far=MRI_FLOAT_PTR(fim) ;
         for( v=0 ; v<s->nvox_mask ; v++ ){
           float a=far[s->vox[v]] ;
           if( !isfinite(a) ) ERROR_exit("Subj %s, item %s: volume %d of '%s' has non-finite values in the mask",s->label,d->item_key[j],s->ibrick[j],s->file[f]) ;
           if( a!=0.0f ) nz=1 ;
           s->x[(size_t)v*n+j]=a ;
         }
         mri_free(fim) ;
         if( !nz ) ERROR_exit("Subj %s, item %s: volume %d of '%s' is all zero in the mask (a failed or missing estimate?)",
                              s->label,d->item_key[j],s->ibrick[j],s->file[f]) ;
         snprintf(extra,sizeof(extra),"Subj %s (%d/%d), %s",s->label,ss+1,d->nsub,THD_trailname(s->file[f],0)) ;
         srm_prog_set(&bar,++done,extra) ;
       }
       DSET_delete(ds) ;
     }
   }
}

/*! The item "timeline": every column kept, in key order. */
static void srm_data_items_timeline( srm_data *d )
{
   int j ;
   d->T=d->nitem ;
   d->keep=(byte *)malloc((size_t)d->T) ; memset(d->keep,1,(size_t)d->T) ;
   d->time_index=(int *)malloc(sizeof(int)*(size_t)d->T) ;
   for( j=0 ; j<d->T ; j++ ) d->time_index[j]=j ;
}

/*! -axis item: z-score every feature over all of the subject's items (-zscore all), or over the
    items of each of the subject's own runs (-zscore run); with -zscore none, only check for features
    that never vary.  A run's items need not be adjacent columns, so they are gathered by run, done
    as segments, and put back. */
static int srm_item_normalize( srm_subj *s , int nitem , int nrun , int zscore , unsigned char *bad )
{
   int *perm=(int *)malloc(sizeof(int)*(size_t)nitem) , *seg=(int *)calloc((size_t)nrun+1,sizeof(int)) ;
   int nseg=0 , rr , j , q=0 , f , nb ; float *y ;
   if( zscore==SRM_ZSCORE_RUN ){
     for( rr=0 ; rr<nrun ; rr++ ){
       int c=0 ;
       for( j=0 ; j<nitem ; j++ ) if( s->irun[j]==rr ){ perm[q++]=j ; c++ ; }
       if( c>0 ) seg[nseg++]=c ;
     }
   } else {
     for( j=0 ; j<nitem ; j++ ) perm[j]=j ;
     seg[0]=nitem ; nseg=1 ;
   }
   y=(float *)malloc(sizeof(float)*(size_t)s->nvox_mask*(size_t)nitem) ;
   if( y==NULL ) ERROR_exit("out of memory z-scoring Subj %s",s->label) ;
   for( f=0 ; f<s->nvox_mask ; f++ ){
     const float *xr=s->x+(size_t)f*nitem ; float *yr=y+(size_t)f*nitem ;
     for( j=0 ; j<nitem ; j++ ) yr[j]=xr[perm[j]] ;
   }
   if( zscore!=SRM_ZSCORE_NONE ){
     nb=SRM_zscore_segments(y,s->nvox_mask,nitem,nseg,seg,bad) ;
     if( nb<0 ) ERROR_exit("Subj %s: %s",s->label,SRM_strerror(nb)) ;
   } else nb=srm_flag_constant(y,s->nvox_mask,nitem,nseg,seg,bad) ;
   for( f=0 ; f<s->nvox_mask ; f++ ){
     float *xr=s->x+(size_t)f*nitem ; const float *yr=y+(size_t)f*nitem ;
     for( j=0 ; j<nitem ; j++ ) xr[perm[j]]=yr[j] ;
   }
   free(y) ; free(perm) ; free(seg) ;
   return nb ;
}

/*! Folds of the item axis: the levels of -fold_column (the Run column by default), in order of first
    appearance.  A fold is a set of item columns, the same for every subject, so the column must give
    an item the same value in every subject; runs do that only when the run assignment was not
    randomized. */
static void srm_make_item_folds( const srm_opts *o , const srm_data *d , srm_folds *f )
{
   THD_datatable *tab=d->tab ; const char *fc=(o->fold_col!=NULL) ? o->fold_col : o->run_col ;
   int icol=THD_datatable_column(tab,(char *)fc) , j , ss , q , n=d->nitem , *lev ;
   memset(f,0,sizeof(*f)) ; f->items=1 ; f->unit="items" ;
   if( icol<0 ){
     if( o->fold_col!=NULL ) ERROR_exit("-fold_column: the table has no '%s' column",fc) ;
     ERROR_exit("-cv with -axis item needs folds: a '%s' column (when every subject has the same items in each run), "
                "or -fold_column",o->run_col) ;
   }
   lev=(int *)malloc(sizeof(int)*(size_t)n) ;
   for( j=0 ; j<n ; j++ ){
     const char *v0=DT_CELL(tab,d->sb[0].irow[j],icol) ;
     for( ss=1 ; ss<d->nsub ; ss++ ){
       const char *v=DT_CELL(tab,d->sb[ss].irow[j],icol) ;
       if( strcmp(v,v0)==0 ) continue ;
       if( o->fold_col==NULL )
         ERROR_exit("-cv: item %s is in run %s for Subj %s but in run %s for Subj %s, so the runs cannot be folds (a fold "
                    "must hold out the same items for everyone). Give -fold_column a column whose value is the same for "
                    "an item in every subject: a fold number, or a category for leave-one-category-out",
                    d->item_key[j],v0,d->sb[0].label,v,d->sb[ss].label) ;
       ERROR_exit("-fold_column %s: item %s is '%s' for Subj %s but '%s' for Subj %s; a fold must hold out the same "
                  "items for everyone",fc,d->item_key[j],v0,d->sb[0].label,v,d->sb[ss].label) ;
     }
     for( q=0 ; q<f->nfold ; q++ ) if( strcmp(f->name[q],v0)==0 ) break ;
     if( q==f->nfold ){ f->name=(char **)realloc(f->name,sizeof(char *)*(size_t)(f->nfold+1)) ; f->name[f->nfold++]=strdup(v0) ; }
     lev[j]=q ;
   }
   if( f->nfold<2 ) ERROR_exit("%s '%s' has one level (%s); -cv needs at least 2 folds",
                               o->fold_col!=NULL?"-fold_column":"the run column",fc,f->name[0]) ;
   f->len=(int *)calloc((size_t)f->nfold,sizeof(int)) ; f->col=(int **)malloc(sizeof(int *)*(size_t)f->nfold) ;
   for( j=0 ; j<n ; j++ ) f->len[lev[j]]++ ;
   for( q=0 ; q<f->nfold ; q++ ){
     int c=0 ;
     f->col[q]=(int *)malloc(sizeof(int)*(size_t)f->len[q]) ;
     for( j=0 ; j<n ; j++ ) if( lev[j]==q ) f->col[q][c++]=j ;
   }
   free(lev) ;
}

/*! z-score within runs and deal with features that never vary: drop them (fitting), or zero them
    and keep their place (applying a saved W, whose rows are fixed). */
static void srm_data_normalize( const srm_opts *o , srm_data *d , int zscore , int Kmin , int keep_rows )
{
   int ss , total_drop=0 ; srm_prog bar ; char extra[160] ;
   srm_prog_begin(&bar,o->progress,o->quiet,zscore!=SRM_ZSCORE_NONE ? "z-score" : "checking","subjects",d->nsub) ;
   for( ss=0 ; ss<d->nsub ; ss++ ){
     srm_subj *s=d->sb+ss ;
     unsigned char *bad=(unsigned char *)calloc((size_t)s->nvox_mask,1) ; int nb ;
     if( d->axis_item ) nb=srm_item_normalize(s,d->T,d->nrun,zscore,bad) ;
     else if( zscore!=SRM_ZSCORE_NONE ){
       int one=d->T ;
       nb = (zscore==SRM_ZSCORE_ALL) ? SRM_zscore_segments(s->x,s->nvox_mask,d->T,1,&one,bad)
                                     : SRM_zscore_segments(s->x,s->nvox_mask,d->T,d->nrun,d->seglen,bad) ;
       if( nb<0 ) ERROR_exit("Subj %s: %s",s->label,SRM_strerror(nb)) ;
     } else nb=srm_flag_constant(s->x,s->nvox_mask,d->T,d->nrun,d->seglen,bad) ;
     if( keep_rows ){
       int f ;
       for( f=0 ; f<s->nvox_mask ; f++ ) if( bad[f] ) memset(s->x+(size_t)f*d->T,0,sizeof(float)*(size_t)d->T) ;
       s->nvox=s->nvox_mask ; s->ndrop=nb ; total_drop+=nb ;
     } else {
       s->nvox=srm_compact_rows(s->x,s->nvox_mask,d->T,bad,s->vox) ;
       s->ndrop=s->nvox_mask-s->nvox ; total_drop+=s->ndrop ;
       if( s->nvox<Kmin ) ERROR_exit("Subj %s has %d features with any variance; -nfeatures is %d",s->label,s->nvox,Kmin) ;
       if( s->nvox<s->nvox_mask ){
         float *shrunk=(float *)realloc(s->x,sizeof(float)*(size_t)s->nvox*(size_t)d->T) ;
         if( shrunk!=NULL ) s->x=shrunk ;
       }
     }
     free(bad) ;
     snprintf(extra,sizeof(extra),"Subj %s, %d feature%s without variance so far",s->label,total_drop,total_drop==1?"":"s") ;
     srm_prog_set(&bar,ss+1,extra) ;
   }
   if( !o->quiet && total_drop>0 ){
     if( keep_rows ) INFO_message("%d features have no variance in some run of the new data; they are zeroed",total_drop) ;
     else INFO_message("%d feature%s with no variance in some run %s dropped (see %s.fit.txt)",total_drop,
                       total_drop==1?"":"s",total_drop==1?"was":"were",o->prefix) ;
   }
}

/*============================================================================
  Pieces shared by fitting, -model_in, and -apply_W
==============================================================================*/

/*! Rows of a kept-feature matrix (nkept x ncol) put back on the input's own columns,
    zeros for the features that were dropped. */
static float * srm_expand_rows( const float *a , int nkept , int ncol , const int *vox , int nfull )
{
   float *y=(float *)calloc((size_t)nfull*(size_t)ncol,sizeof(float)) ; int f ;
   if( y==NULL ) ERROR_exit("out of memory writing a matrix (%.2f GiB)",4.0*nfull*(double)ncol/GIB) ;
   for( f=0 ; f<nkept ; f++ ) memcpy(y+(size_t)vox[f]*ncol,a+(size_t)f*ncol,sizeof(float)*(size_t)ncol) ;
   return y ;
}

/*! A 1D comment: head, then the nkey item keys (the rows, in order) on a '# items:' line.  Free it. */
static char * srm_keys_comment( const char *head , char * const *key , int nkey )
{
   size_t len=strlen(head)+16 ; int j ; char *c ;
   for( j=0 ; j<nkey ; j++ ) len+=strlen(key[j])+1 ;
   c=(char *)malloc(len) ;
   strcpy(c,head) ; strcat(c,"\n# items:") ;
   for( j=0 ; j<nkey ; j++ ){ strcat(c," ") ; strcat(c,key[j]) ; }
   return c ;
}

/*! One subject's series in the shared space, K x T: a k x 1 x 1 dataset, or a text file. */
static void srm_write_proj( const srm_opts *o , const srm_subj *s , const float *P , int K , int T ,
                            const char *comment , char * const *collab , int argc , char **argv )
{
   char sn[128] , fn[THD_MAX_NAME] ; srm_safe(sn,sizeof(sn),s->label) ;
   if( o->proj_dset ){
     THD_3dim_dataset *out=EDIT_empty_copy(NULL) ; THD_ivec3 ixyz ; THD_fvec3 dxyz , oxyz ; float TR=1.0f ; int t ;
     AFNI_setenv("AFNI_WRITE_1D_AS_PREFIX=YES") ;   /* a k x 1 x 1 dataset is otherwise written as a 1D text file */
     if( s->geom!=NULL && DSET_TR(s->geom)>0.0f ) TR=DSET_TR(s->geom) ;
     LOAD_IVEC3(ixyz,K,1,1) ; LOAD_FVEC3(dxyz,1.0,1.0,1.0) ; LOAD_FVEC3(oxyz,0.0,0.0,0.0) ;
     snprintf(fn,sizeof(fn),"%s.proj.%s",o->prefix,sn) ;
     if( collab!=NULL )                             /* items: a bucket, each volume labelled with its key */
       EDIT_dset_items(out,ADN_nxyz,ixyz,ADN_xyzdel,dxyz,ADN_xyzorg,oxyz,ADN_prefix,fn,ADN_nvals,T,ADN_ntt,0,
                       ADN_datum_all,MRI_float,ADN_none) ;
     else
       EDIT_dset_items(out,ADN_nxyz,ixyz,ADN_xyzdel,dxyz,ADN_xyzorg,oxyz,ADN_prefix,fn,ADN_nvals,T,ADN_ntt,T,
                       ADN_ttdel,TR,ADN_ttorg,0.0,ADN_tunits,UNITS_SEC_TYPE,ADN_datum_all,MRI_float,ADN_none) ;
     tross_Make_History("3dSRM",argc,argv,out) ;
     for( t=0 ; t<T ; t++ ){
       float *col=(float *)malloc(sizeof(float)*(size_t)K) ; int c ;
       for( c=0 ; c<K ; c++ ) col[c]=P[(size_t)c*T+t] ;
       EDIT_substitute_brick(out,t,MRI_float,col) ;
       if( collab!=NULL ) EDIT_BRICK_LABEL(out,t,collab[t]) ;
     }
     DSET_write(out) ; DSET_delete(out) ;
   } else {
     snprintf(fn,sizeof(fn),"%s.proj.%s.1D",o->prefix,sn) ;
     if( collab!=NULL ){
       char *ck=srm_keys_comment(comment,collab,T) ;
       srm_write_transposed_1D(fn,P,K,T,ck) ; free(ck) ;
     } else srm_write_transposed_1D(fn,P,K,T,comment) ;
   }
}

/*! The outputs that exist for every subject: the data in the shared space (held out, when res
    has them), W, voxelwise held-out R^2, and the reconstruction.  Subject ss is view voff+ss of m. */
static void srm_write_subject_outputs( const srm_opts *o , const srm_data *d , const SRM_model *m , int voff ,
                                       const srm_cvres *res , int K , int argc , char **argv )
{
   srm_prog pbar ; int nwork=(!o->no_projected)+(o->save_W)+(o->save_recon)+(o->save_r2) , ss , rc , T=d->T ;
   srm_prog_begin(&pbar,o->progress,o->quiet,"writing","subjects",nwork>0 ? d->nsub : 0) ;
   for( ss=0 ; ss<d->nsub ; ss++ ){
     const srm_subj *s=d->sb+ss ; char sn[128] , fn[THD_MAX_NAME] ; srm_safe(sn,sizeof(sn),s->label) ;
     if( !o->no_projected ){
       float *P=(float *)malloc(sizeof(float)*(size_t)K*(size_t)T) ; int held=(res!=NULL && res->proj!=NULL) ;
       if( held ) memcpy(P,res->proj[ss],sizeof(float)*(size_t)K*(size_t)T) ;
       else {
         rc=SRM_project(m,voff+ss,s->x,T,P) ;
         if( rc!=SRM_OK ) ERROR_exit("Subj %s: projection failed: %s",s->label,SRM_strerror(rc)) ;
       }
       if( d->axis_item )
         srm_write_proj(o,s,P,K,T,held ? "subject's HELD-OUT items in the shared space (each from a model that never saw it), aligned to the full fit: one row per item, one column per component"
                                       : "subject's items in the shared space: one row per item, one column per component",d->item_key,argc,argv) ;
       else
         srm_write_proj(o,s,P,K,T,held ? "subject's HELD-OUT data in the shared space (each stretch from a model that never saw it), aligned to the full fit: one row per retained time point, one column per component"
                                       : "subject's data in the shared space: one row per retained time point, one column per component",NULL,argc,argv) ;
       free(P) ;
     }
     if( o->save_W ){
       if( s->is1D ){
         float *Wf=srm_expand_rows(m->w[voff+ss],s->nvox,K,s->vox,s->nvox_mask) ;
         snprintf(fn,sizeof(fn),"%s.W.%s.1D",o->prefix,sn) ;
         srm_write_matrix_1D(fn,Wf,s->nvox_mask,K,"mapping W: one row per feature (column of the input; zeros for dropped features), one column per component") ;
         free(Wf) ;
       } else {
         int nvg=DSET_NVOX(s->geom) , c , f ; float **vols=(float **)malloc(sizeof(float *)*(size_t)K) ;
         for( c=0 ; c<K ; c++ ){
           vols[c]=(float *)calloc((size_t)nvg,sizeof(float)) ;
           for( f=0 ; f<s->nvox ; f++ ) vols[c][s->vox[f]]=m->w[voff+ss][(size_t)f*K+c] ;
         }
         snprintf(fn,sizeof(fn),"%s.W.%s",o->prefix,sn) ;
         srm_write_grid_dset(s->geom,fn,K,vols,"comp",argc,argv) ;
         free(vols) ;
       }
     }
     if( o->save_r2 && res!=NULL && res->vr2!=NULL ){
       if( s->is1D ){
         float *rf=srm_expand_rows(res->vr2[ss],s->nvox,1,s->vox,s->nvox_mask) ;
         snprintf(fn,sizeof(fn),"%s.r2.%s.1D",o->prefix,sn) ;
         srm_write_matrix_1D(fn,rf,s->nvox_mask,1,"held-out R^2 of each feature, the mean over folds (0 for dropped features)") ;
         free(rf) ;
       } else {
         int nvg=DSET_NVOX(s->geom) , f ; float **vols=(float **)malloc(sizeof(float *)) ;
         vols[0]=(float *)calloc((size_t)nvg,sizeof(float)) ;
         for( f=0 ; f<s->nvox ; f++ ) vols[0][s->vox[f]]=res->vr2[ss][f] ;
         snprintf(fn,sizeof(fn),"%s.r2.%s",o->prefix,sn) ;
         srm_write_grid_dset(s->geom,fn,1,vols,"r2_",argc,argv) ;
         free(vols) ;
       }
     }
     if( o->save_recon ){
       float *xh=(float *)malloc(sizeof(float)*(size_t)s->nvox*(size_t)T) ;
       rc=SRM_reconstruct(m,voff+ss,m->s,T,xh) ;
       if( rc!=SRM_OK ) ERROR_exit("Subj %s: reconstruction failed: %s",s->label,SRM_strerror(rc)) ;
       if( s->is1D ){
         float *xf=srm_expand_rows(xh,s->nvox,T,s->vox,s->nvox_mask) ;
         snprintf(fn,sizeof(fn),"%s.recon.%s.1D",o->prefix,sn) ;
         srm_write_transposed_1D(fn,xf,s->nvox_mask,T,"W S: one row per retained time point, one column per feature (zeros for dropped features)") ;
         free(xf) ;
       } else {
         int nvg=DSET_NVOX(s->geom) , t , f ; float **vols=(float **)malloc(sizeof(float *)*(size_t)T) ;
         THD_3dim_dataset *out ;
         for( t=0 ; t<T ; t++ ){
           vols[t]=(float *)calloc((size_t)nvg,sizeof(float)) ;
           for( f=0 ; f<s->nvox ; f++ ) vols[t][s->vox[f]]=xh[(size_t)f*T+t] ;
         }
         snprintf(fn,sizeof(fn),"%s.recon.%s",o->prefix,sn) ;
         out=EDIT_empty_copy(s->geom) ;
         EDIT_dset_items(out,ADN_prefix,fn,ADN_nvals,T,ADN_ntt,T,ADN_datum_all,MRI_float,ADN_brick_fac,NULL,ADN_none) ;
         tross_Copy_History(s->geom,out) ; tross_Make_History("3dSRM",argc,argv,out) ;
         for( t=0 ; t<T ; t++ ) EDIT_substitute_brick(out,t,MRI_float,vols[t]) ;
         DSET_write(out) ; DSET_delete(out) ; free(vols) ;
       }
       free(xh) ;
     }
     srm_prog_set(&pbar,ss+1,NULL) ;
   }
}

/*! Refuse to overwrite any per-subject output, before the work is done. */
static void srm_check_subject_outputs( const srm_opts *o , const srm_data *d )
{
   char fn[THD_MAX_NAME] ; int ss ;
   if( THD_ok_overwrite() ) return ;
   for( ss=0 ; ss<d->nsub ; ss++ ){
     const srm_subj *s=d->sb+ss ; char sn[128] ; srm_safe(sn,sizeof(sn),s->label) ;
     if( !o->no_projected ){
       if( o->proj_dset ){ snprintf(fn,sizeof(fn),"%s.proj.%s",o->prefix,sn) ; if( srm_dset_exists(fn) ) ERROR_exit("a dataset '%s' already exists (use -overwrite)",fn) ; }
       else { snprintf(fn,sizeof(fn),"%s.proj.%s.1D",o->prefix,sn) ; if( THD_is_file(fn) ) ERROR_exit("'%s' already exists (use -overwrite)",fn) ; }
     }
     if( o->save_W ){ snprintf(fn,sizeof(fn),s->is1D?"%s.W.%s.1D":"%s.W.%s",o->prefix,sn) ;
       if( s->is1D ? THD_is_file(fn) : srm_dset_exists(fn) ) ERROR_exit("'%s' already exists (use -overwrite)",fn) ; }
     if( o->save_r2 ){ snprintf(fn,sizeof(fn),s->is1D?"%s.r2.%s.1D":"%s.r2.%s",o->prefix,sn) ;
       if( s->is1D ? THD_is_file(fn) : srm_dset_exists(fn) ) ERROR_exit("'%s' already exists (use -overwrite)",fn) ; }
     if( o->save_recon ){ snprintf(fn,sizeof(fn),s->is1D?"%s.recon.%s.1D":"%s.recon.%s",o->prefix,sn) ;
       if( s->is1D ? THD_is_file(fn) : srm_dset_exists(fn) ) ERROR_exit("'%s' already exists (use -overwrite)",fn) ; }
   }
}

static void srm_check_new_file( const srm_opts *o , const char *suffix )
{
   char fn[THD_MAX_NAME] ;
   if( THD_ok_overwrite() ) return ;
   snprintf(fn,sizeof(fn),"%s%s",o->prefix,suffix) ;
   if( THD_is_file(fn) ) ERROR_exit("'%s' already exists (use -overwrite)",fn) ;
}

/*! The memory contract, before the big allocations. */
static void srm_memory_preflight( const srm_opts *o , const srm_data *d , int K , int do_cv , int nsegfit )
{
   THD_memory_plan me ; double inb=0.0 , Wb=0.0 ; int *nv0=(int *)malloc(sizeof(int)*(size_t)d->nsub) , ss ;
   for( ss=0 ; ss<d->nsub ; ss++ ){ nv0[ss]=d->sb[ss].nvox_mask ; inb+=4.0*d->sb[ss].nvox_mask*(double)d->Tin ; Wb+=4.0*d->sb[ss].nvox_mask*K ; }
   memset(&me,0,sizeof(me)) ;
   me.input=inb*(do_cv ? 2.0 : 1.0) ;        /* a cross-validation fold holds a second copy of the data */
   me.shared=SRM_estimate_bytes(d->nsub,nv0,d->Tin,K,o->par.algo,0)
            + Wb*(nsegfit>0 ? (double)nsegfit : 0.0) + (o->nstart>1 ? Wb : 0.0) ;
   me.system=(double)AFNI_get_memsize() ;
   me.limit=o->memory_limit_given ? o->memory_limit_gib*GIB : ((me.system>0.0) ? 0.80*me.system : 0.0) ;
   me.nthread=1 ; THD_memory_plan_finish(&me) ;
   if( !o->quiet && (me.total>0.50*me.system || o->memory_limit_given) ){
     if( me.limit>0.0 )
       INFO_message("memory preflight: about %.2f GiB peak (%.2f GiB of data + %.2f GiB for the fit); limit %.2f GiB",
                    me.total/GIB,me.input/GIB,me.shared/GIB,me.limit/GIB) ;
     else INFO_message("memory preflight: about %.2f GiB peak; physical memory is unavailable",me.total/GIB) ;
   }
   if( me.limit>0.0 && me.total>me.limit ){
     if( !o->memory_override )
       ERROR_exit("estimated peak %.2f GiB exceeds the %.2f GiB limit: use fewer voxels or subjects, a smaller -nfeatures, "
                  "set -memory_limit, or acknowledge with -memory_override",me.total/GIB,me.limit/GIB) ;
     WARNING_message("estimated peak exceeds the memory limit; continuing because -memory_override was given") ;
   }
   free(nv0) ;
}

/*! The numeric table columns that travel with the scores (age, motion, ...). */
static int * srm_carry_columns( const THD_datatable *tab , const srm_data *d , int *ncarry )
{
   int *carry=NULL , jj ;
   *ncarry=0 ;
   for( jj=0 ; jj<tab->ncol ; jj++ )
     if( tab->isnum[jj] && jj!=tab->icol_subj && jj!=tab->icol_input && jj!=d->icol_run && jj!=d->icol_mask && jj!=d->icol_group ){
       if( d->axis_item ){                         /* a column that varies over a subject's items is about items */
         int ss , k , same=1 ;
         for( ss=0 ; ss<d->nsub && same ; ss++ ) for( k=1 ; k<d->nitem && same ; k++ )
           same=(strcmp(DT_CELL(tab,d->sb[ss].irow[k],jj),DT_CELL(tab,d->sb[ss].irow[0],jj))==0) ;
         if( !same ) continue ;
       }
       carry=(int *)realloc(carry,sizeof(int)*(size_t)(*ncarry+1)) ; carry[(*ncarry)++]=jj ;
     }
   return carry ;
}

static void srm_int_list( char *dst , size_t nd , const int *a , int n )
{
   size_t pos=0 ; int ii ;
   dst[0]='\0' ;
   for( ii=0 ; ii<n && pos+12<nd ; ii++ ) pos+=snprintf(dst+pos,nd-pos,"%s%d",ii?",":"",a[ii]) ;
   if( ii<n && nd>4 ) snprintf(dst+nd-4,4,"...") ;
}

/*============================================================================
  -model_in: place subjects in a saved model; optionally let S absorb them
==============================================================================*/

static void srm_run_model_in( srm_opts *o , THD_datatable *tab , int argc , char **argv )
{
   srm_data d ; SRM_model *m=NULL ; SRM_provenance *pv=NULL ; SRM_fitinfo fi ; SRM_params par ;
   srm_ctx cx ; srm_want want ; srm_folds folds ; srm_cvres res ;
   int rc , ss , rr , q , K , T , nold , nsub , do_cv=(strcasecmp(o->cvmode,"none")!=0) , have_res=0 ;
   int *carry=NULL , ncarry=0 , nupd0 ; double drift0 ;
   memset(&d,0,sizeof(d)) ; memset(&fi,0,sizeof(fi)) ; memset(&want,0,sizeof(want)) ; memset(&folds,0,sizeof(folds)) ; memset(&res,0,sizeof(res)) ;

   rc=SRM_load_model(o->model_in,&m,&pv) ;
   if( rc!=SRM_OK ) ERROR_exit("cannot use the model '%s': %s",o->model_in,SRM_strerror(rc)) ;
   if( pv->axis!=SRM_AXIS_TIME ) ERROR_exit("'%s' is a model on the item axis; this build can only place subjects in a time-axis model",o->model_in) ;
   K=m->nfeat ; nold=m->nview ;
   o->par.algo=m->algo ; o->zscore=pv->zscore ;
   if( !o->quiet ){
     INFO_message("model %s: %s, %d components, %d subjects, %d time points in %d run%s%s",o->model_in,
                  m->algo==SRM_ALGO_PROB?"probabilistic":"deterministic",K,nold,m->ntime,pv->nseg,pv->nseg>1?"s":"",
                  m->coupled?", a W per run":"") ;
     if( m->nupdate>0 ) ININFO_message("it has absorbed new subjects %d time%s; drift %.4g",m->nupdate,m->nupdate>1?"s":"",m->drift) ;
   }

   srm_data_table(o,tab,&d,1) ;
   srm_data_shape(o,&d,1,K) ;
   nsub=d.nsub ;

   /*-- the new data must be the model's own stimulus, in the model's own runs --*/
   {
     int *sl=(int *)calloc((size_t)d.nrun,sizeof(int)) , t ;
     rc=SRM_OK ;
     if( d.Tin!=pv->ntime_in || d.nrun!=pv->nseg ) rc=SRM_ERR_MISMATCH ;
     else {
       for( rr=0 ; rr<d.nrun ; rr++ ) for( t=0 ; t<d.runlen[rr] ; t++ ) sl[rr]+=(pv->keep[d.c0in[rr]+t]!=0) ;
       rc=SRM_provenance_check(pv,d.Tin,d.nrun,sl) ;
     }
     if( rc!=SRM_OK ){
       char a[256] , b[256] ;
       srm_int_list(a,sizeof(a),d.runlen,d.nrun) ;
       srm_int_list(b,sizeof(b),pv->seglen,pv->nseg) ;
       ERROR_exit("the new data do not fit the model: they have %d time points in %d run%s (lengths %s); the model was built from "
                  "%d time points in %d run%s (%d kept, per run %s). The stimulus, runs, and time points must be the same",
                  d.Tin,d.nrun,d.nrun>1?"s":"",a,pv->ntime_in,pv->nseg,pv->nseg>1?"s":"",m->ntime,b) ;
     }
     free(sl) ;
   }
   srm_data_censor(o,&d,pv->keep) ;

   /*-- subjects already in the model --*/
   for( ss=0 ; ss<nsub ; ss++ ){
     char sn[128] ; srm_safe(sn,sizeof(sn),d.sb[ss].label) ;
     for( q=0 ; q<pv->nsubj ; q++ ) if( strcmp(sn,pv->label[q])==0 ){
       if( o->update_S ) ERROR_exit("Subj %s is already one of the %d subjects in the model; -update_S would count them twice",d.sb[ss].label,pv->nsubj) ;
       WARNING_message("Subj %s is also a subject the model was fit to: their S already contains this person, so their scores are optimistic",d.sb[ss].label) ;
       break ;
     }
   }

   /*-- do not clobber anything --*/
   if( o->update_S ){ srm_check_new_file(o,".srm") ; srm_check_new_file(o,".S.1D") ; srm_check_new_file(o,".objective.1D") ;
                      if( m->coupled ) srm_check_new_file(o,".coupling.txt") ; }
   srm_check_new_file(o,".fit.txt") ;
   if( m->algo==SRM_ALGO_PROB ) srm_check_new_file(o,".rho2.1D") ;
   if( do_cv ){ srm_check_new_file(o,".scores.txt") ; srm_check_new_file(o,".scores_folds.txt") ; }
   srm_check_subject_outputs(o,&d) ;
   srm_memory_preflight(o,&d,K,do_cv,0) ;

   srm_data_read(o,&d) ;
   srm_data_timeline(o,&d,K) ;
   srm_data_normalize(o,&d,pv->zscore,K,0) ;
   T=d.T ;
   if( T!=m->ntime ) ERROR_exit("internal error: %d time points retained; the model has %d",T,m->ntime) ;

   /*-- scores first, from the model as it stands: the new subjects have not touched it --*/
   memset(&cx,0,sizeof(cx)) ;
   cx.o=o ; cx.sb=d.sb ; cx.nsub=nsub ; cx.nrun=d.nrun ; cx.T=T ; cx.K=K ; cx.seglen=d.seglen ; cx.runname=d.runname ;
   cx.glabel=d.glabel ; cx.ngroup=d.ngroup ; cx.quiet=o->quiet ; cx.progress=o->progress ;
   par=o->par ; par.nfeat=K ; par.center=m->center ; cx.par=par ;
   if( do_cv ){
     want=srm_parse_metrics(o->cv_metrics) ;
     srm_make_folds(&folds,d.nrun,d.seglen,d.runname,T,o->cv_blocks) ;
     srm_check_fold_lengths(&folds,want,o->match_window) ;
     { int ff , big=0 ;
       for( ff=1 ; ff<folds.nfold ; ff++ ) if( folds.len[ff]>folds.len[big] ) big=ff ;
       if( T-folds.len[big]<K ) ERROR_exit("holding out fold %s leaves %d time points to fit W; the model has %d components",
                                           folds.name[big],T-folds.len[big],K) ; }
     if( !o->quiet ){
       INFO_message("-cv %s with -model_in: S stays fixed; each subject's W is estimated without each of %d %s and scored on it",
                    o->cvmode,folds.nfold,folds.by_blocks?"time blocks":"runs") ;
       if( folds.by_blocks ) WARNING_message("one run only: the folds are contiguous time blocks. Time points next to a block's edge "
                                             "resemble the other side, so scores are somewhat optimistic; separate runs are better") ;
     }
     carry=srm_carry_columns(tab,&d,&ncarry) ;
     srm_cvres_alloc(&res,nsub,folds.nfold,!o->no_projected,o->save_r2,&cx) ; have_res=1 ;
     srm_cv_frozen(&cx,m,&folds,want,&res) ;
     srm_cv_summary(&cx,&res,want,o->cvmode,folds.nfold) ;
   }

   /*-- then the new subjects join: a W each, and (with -update_S) S moves --*/
   drift0=m->drift ; nupd0=m->nupdate ;
   if( o->update_S ){
     int *nv=(int *)malloc(sizeof(int)*(size_t)nsub) ; float **xv=(float **)malloc(sizeof(float *)*(size_t)nsub) ;
     srm_prog bar ; fit_prog fp ;
     for( ss=0 ; ss<nsub ; ss++ ){ nv[ss]=d.sb[ss].nvox ; xv[ss]=d.sb[ss].x ; }
     if( m->coupled ){ par.coupling=m->coupling ; par.nseg=pv->nseg ; par.seglen=pv->seglen ; }
     if( !o->quiet ) INFO_message("updating S with %d new subject%s (up to %d iterations)",nsub,nsub>1?"s":"",par.niter) ;
     srm_prog_begin(&bar,o->progress,o->quiet,"update","iterations",par.niter) ;
     fp.bar=&bar ; fp.istart=0 ; fp.nstart=1 ; fp.first_obj=0.0 ; fp.tag=NULL ;
     par.progress=srm_fit_progress ; par.progress_ctx=&fp ;
     rc=SRM_update(m,nsub,nv,xv,&par,&fi) ;
     srm_prog_end(&bar) ;
     if( rc!=SRM_OK ) ERROR_exit("the update failed: %s",SRM_strerror(rc)) ;
     if( !o->quiet ){
       INFO_message("%s after %d iterations; this update moved S by %.4g (the model's drift is now %.4g after %d update%s)",
                    fi.converged?"converged":"stopped at the iteration limit",fi.niter_done,m->drift-drift0,m->drift,m->nupdate,m->nupdate>1?"s":"") ;
       if( m->drift>0.03 ) WARNING_message("the model's drift is %.3g: the oldest subjects' W were fitted to a noticeably different S. "
                                           "Refit from the data when you can",m->drift) ;
     }
     free(nv) ; free(xv) ;
   } else {
     srm_prog bar ;
     srm_prog_begin(&bar,o->progress,o->quiet,"mapping","subjects",nsub) ;
     for( ss=0 ; ss<nsub ; ss++ ){
       rc=SRM_model_add_view(m,d.sb[ss].x,d.sb[ss].nvox) ;
       if( rc!=SRM_OK ) ERROR_exit("Subj %s: could not estimate W against the model: %s",d.sb[ss].label,SRM_strerror(rc)) ;
       srm_prog_set(&bar,ss+1,NULL) ;
     }
     srm_prog_end(&bar) ;
   }

   /*-- outputs --*/
   {
     char fn[THD_MAX_NAME] ; FILE *fq ; int ii ;
     if( o->update_S ){
       SRM_provenance *pn=pv ;
       snprintf(fn,sizeof(fn),"%s.S.1D",o->prefix) ;
       { char c2[256] ; snprintf(c2,sizeof(c2),"shared response after the update: %d time points (retained) x %d components",T,K) ;
         srm_write_transposed_1D(fn,m->s,K,T,c2) ; }
       if( !o->quiet ) INFO_message("wrote %s",fn) ;
       snprintf(fn,sizeof(fn),"%s.objective.1D",o->prefix) ;
       fq=fopen(fn,"w") ; if( fq==NULL ) ERROR_exit("cannot write '%s'",fn) ;
       fprintf(fq,"# objective of the restricted problem at the start of each iteration, and at the end (%s); comparable within this file only\n",
               m->algo==SRM_ALGO_PROB?"log-likelihood, rises":"summed squared error, falls") ;
       for( ii=0 ; ii<=fi.niter_done ; ii++ ) fprintf(fq,"%.12g\n",fi.objective[ii]) ;
       fclose(fq) ;
       pn->label=(char **)realloc(pn->label,sizeof(char *)*(size_t)(pn->nsubj+nsub)) ;
       for( ss=0 ; ss<nsub ; ss++ ){ char sn[128] ; srm_safe(sn,sizeof(sn),d.sb[ss].label) ; pn->label[pn->nsubj+ss]=strdup(sn) ; }
       pn->nsubj+=nsub ;
       snprintf(fn,sizeof(fn),"%s.srm",o->prefix) ;
       rc=SRM_save_model(fn,m,pn) ;
       if( rc!=SRM_OK ) ERROR_exit("cannot write the model '%s': %s",fn,SRM_strerror(rc)) ;
       if( !o->quiet ) INFO_message("wrote %s (%d subjects)",fn,m->nview) ;
     }
     if( m->algo==SRM_ALGO_PROB ){
       snprintf(fn,sizeof(fn),"%s.rho2.1D",o->prefix) ;
       fq=fopen(fn,"w") ; if( fq==NULL ) ERROR_exit("cannot write '%s'",fn) ;
       fprintf(fq,"# noise variance per new subject, in the order of the table:") ;
       for( ss=0 ; ss<nsub ; ss++ ) fprintf(fq," %s",d.sb[ss].label) ;
       fputc('\n',fq) ;
       for( ss=0 ; ss<nsub ; ss++ ) fprintf(fq,"%.9g\n",m->rho2[nold+ss]) ;
       fclose(fq) ;
     }
     snprintf(fn,sizeof(fn),"%s.fit.txt",o->prefix) ;
     fq=fopen(fn,"w") ; if( fq==NULL ) ERROR_exit("cannot write '%s'",fn) ;
     fprintf(fq,"# 3dSRM -model_in %s: %s, %d components, %d new subjects, %d time points (%d runs)\n",
             o->model_in,m->algo==SRM_ALGO_PROB?"probabilistic":"deterministic",K,nsub,T,d.nrun) ;
     if( o->update_S )
       fprintf(fq,"# -update_S: iterations %d, converged %d, this update moved S by %.6g; drift %.6g after %d update%s; model now has %d subjects\n",
               fi.niter_done,fi.converged,m->drift-drift0,m->drift,m->nupdate,m->nupdate==1?"":"s",m->nview) ;
     else
       fprintf(fq,"# S was not changed (no -update_S); the model has %d subjects and %d update%s so far\n",nold,nupd0,nupd0==1?"":"s") ;
     fprintf(fq,"# zscore %s, the model's censoring applied: %d of %d input time points dropped\n",
             pv->zscore==SRM_ZSCORE_RUN?"run":pv->zscore==SRM_ZSCORE_ALL?"all":"none",d.Tin-T,d.Tin) ;
     fprintf(fq,"Subj\tGroup\tFeatures_in_mask\tFeatures_used\tFeatures_dropped\trho2\n") ;
     for( ss=0 ; ss<nsub ; ss++ ){
       char r2s[40] ;
       if( m->algo==SRM_ALGO_PROB ) snprintf(r2s,sizeof(r2s),"%.9g",m->rho2[nold+ss]) ; else snprintf(r2s,sizeof(r2s),"NA") ;
       fprintf(fq,"%s\t%s\t%d\t%d\t%d\t%s\n",d.sb[ss].label,d.glabel[d.sb[ss].group],d.sb[ss].nvox_mask,d.sb[ss].nvox,d.sb[ss].ndrop,r2s) ;
     }
     fclose(fq) ;
     if( o->update_S && m->coupled && m->wrun!=NULL ){
       snprintf(fn,sizeof(fn),"%s.coupling.txt",o->prefix) ;
       fq=fopen(fn,"w") ; if( fq==NULL ) ERROR_exit("cannot write '%s'",fn) ;
       fprintf(fq,"# deviation of each run's W from the subject's anchor: |W_ir - W_i|^2 / %d (0 identical, 2 unrelated); lambda %.6g\n",K,m->coupling) ;
       fprintf(fq,"Subj\tGroup\tRun\tDeviation\n") ;
       for( ss=0 ; ss<nsub ; ss++ ) for( rr=0 ; rr<d.nrun ; rr++ ){
         double dv ; SRM_run_deviation(m,nold+ss,rr,&dv) ;
         fprintf(fq,"%s\t%s\t%s\t%.6g\n",d.sb[ss].label,d.glabel[d.sb[ss].group],d.runname[rr],dv) ;
       }
       fclose(fq) ;
     }
     if( have_res ){
       srm_write_scores(&cx,&folds,&res,want,o->cvmode,tab,carry,ncarry) ;
       if( !o->quiet ) INFO_message("wrote %s.scores.txt and %s.scores_folds.txt",o->prefix,o->prefix) ;
     }
   }
   srm_write_subject_outputs(o,&d,m,nold,have_res?&res:NULL,K,argc,argv) ;
   if( have_res ){ srm_cvres_free(&res) ; srm_free_folds(&folds) ; }
   free(carry) ;
   SRM_free_fitinfo(&fi) ; SRM_free_model(m) ; SRM_provenance_free(pv) ;
}

/*============================================================================
  -apply_W: carry each subject's W from an earlier run over to other data
==============================================================================*/

static void srm_run_apply_W( srm_opts *o , THD_datatable *tab , int argc , char **argv )
{
   srm_data d ; SRM_model *m0=NULL ; SRM_provenance *pv=NULL ; int rc , ss , q , K ;
   char fn[THD_MAX_NAME] ; srm_prog bar ;
   memset(&d,0,sizeof(d)) ;
   snprintf(fn,sizeof(fn),"%s.srm",o->apply_W) ;
   rc=SRM_load_model(fn,&m0,&pv) ;
   if( rc!=SRM_OK ) ERROR_exit("-apply_W needs '%s' from the earlier run: %s",fn,SRM_strerror(rc)) ;
   K=m0->nfeat ; o->zscore=pv->zscore ;
   srm_data_table(o,tab,&d,1) ;

   /*-- every subject must be in the model and have a W, before any work --*/
   for( ss=0 ; ss<d.nsub ; ss++ ){
     const srm_subj *s=d.sb+ss ; char sn[128] ; const char *first=DT_CELL(tab,s->row[0],tab->icol_input) ; int is1=srm_is1D(first) , ok ;
     srm_safe(sn,sizeof(sn),s->label) ;
     for( q=0 ; q<pv->nsubj ; q++ ) if( strcmp(sn,pv->label[q])==0 ) break ;
     if( q==pv->nsubj ) ERROR_exit("Subj %s is not in the model '%s' (it has %d subjects)",s->label,fn,pv->nsubj) ;
     snprintf(fn,sizeof(fn),is1?"%s.W.%s.1D":"%s.W.%s",o->apply_W,sn) ;
     ok = is1 ? THD_is_file(fn) : srm_dset_exists(fn) ;
     if( !ok ) ERROR_exit("no W for Subj %s: '%s' was not found (did the earlier run use -save_W?)",s->label,fn) ;
     if( !THD_ok_overwrite() ){
       if( !o->no_projected ){
         if( o->proj_dset ){ snprintf(fn,sizeof(fn),"%s.proj.%s",o->prefix,sn) ; if( srm_dset_exists(fn) ) ERROR_exit("a dataset '%s' already exists (use -overwrite)",fn) ; }
         else { snprintf(fn,sizeof(fn),"%s.proj.%s.1D",o->prefix,sn) ; if( THD_is_file(fn) ) ERROR_exit("'%s' already exists (use -overwrite)",fn) ; }
       }
     }
   }
   if( !o->quiet )
     INFO_message("%s: %d components, %s; applying each of %d subject%s' own W to the new data",o->apply_W,K,
                  pv->zscore==SRM_ZSCORE_RUN?"z-scored within runs":pv->zscore==SRM_ZSCORE_ALL?"z-scored over all samples":
                  "not z-scored (each feature is centered on the new data)",d.nsub,"s") ;

   srm_prog_begin(&bar,o->progress,o->quiet,"applying","subjects",d.nsub) ;
   for( ss=0 ; ss<d.nsub ; ss++ ){
     srm_data d1 ; srm_opts oq=*o ; srm_subj *s=d.sb+ss ; char sn[128] ;
     float *W=NULL , *mu=NULL , *P ; int V , c , f , t ; SRM_model *mv ;
     srm_safe(sn,sizeof(sn),s->label) ;
     d1=d ; d1.nsub=1 ; d1.sb=s ; d1.mcache=NULL ; d1.nmcache=0 ; d1.gmask=NULL ;
     d1.runlen=NULL ; d1.c0in=NULL ; d1.runname=NULL ; d1.keep=NULL ; d1.zero_any=NULL ; d1.nzero=NULL ;
     d1.seglen=NULL ; d1.time_index=NULL ; d1.T=0 ; d1.Tin=0 ;
     oq.quiet=1 ; oq.maskname=NULL ;
     srm_data_shape(&oq,&d1,0,0) ;                      /* this subject's lengths and grid; no mask */
     if( s->is1D ){
       MRI_IMAGE *im , *fm ; float *ar ; int nx ;
       snprintf(fn,sizeof(fn),"%s.W.%s.1D",o->apply_W,sn) ;
       im=mri_read_1D(fn) ; if( im==NULL ) ERROR_exit("cannot read '%s'",fn) ;
       fm=mri_to_float(im) ; ar=MRI_FLOAT_PTR(fm) ; nx=fm->nx ;
       if( fm->ny!=K ) ERROR_exit("'%s' has %d columns; the model has %d components",fn,fm->ny,K) ;
       if( nx!=s->nvox_mask ) ERROR_exit("Subj %s: W has %d rows but the new data have %d columns",s->label,nx,s->nvox_mask) ;
       V=nx ; W=(float *)malloc(sizeof(float)*(size_t)V*K) ;
       for( f=0 ; f<V ; f++ ) for( c=0 ; c<K ; c++ ) W[(size_t)f*K+c]=ar[(size_t)c*nx+f] ;
       if( fm!=im ) mri_free(fm) ;
       mri_free(im) ;
     } else {
       char wfound[THD_MAX_NAME] ; THD_3dim_dataset *wds ; MRI_IMAGE **bim=(MRI_IMAGE **)malloc(sizeof(MRI_IMAGE *)*(size_t)K) ; float **fb=(float **)malloc(sizeof(float *)*(size_t)K) ; int nvg ; byte *in ;
       snprintf(fn,sizeof(fn),"%s.W.%s",o->apply_W,sn) ;
       if( !srm_dset_find(fn,wfound,sizeof(wfound)) ) ERROR_exit("cannot find the W dataset '%s'",fn) ;
       wds=THD_open_dataset(wfound) ; CHECK_OPEN_ERROR(wds,wfound) ;
       if( DSET_NVALS(wds)!=K ) ERROR_exit("'%s' has %d volumes; the model has %d components",fn,DSET_NVALS(wds),K) ;
       if( !EQUIV_GRIDS(wds,s->geom) )
         ERROR_exit("Subj %s: the new data are not on the grid of '%s' (W is applied on its own grid and voxels)",s->label,fn) ;
       DSET_load(wds) ; nvg=DSET_NVOX(wds) ; in=(byte *)calloc((size_t)nvg,1) ;
       for( c=0 ; c<K ; c++ ){
         bim[c]=THD_extract_float_brick(c,wds) ; fb[c]=MRI_FLOAT_PTR(bim[c]) ;
         for( f=0 ; f<nvg ; f++ ) if( fb[c][f]!=0.0f ) in[f]=1 ;
       }
       for( V=0,f=0 ; f<nvg ; f++ ) V+=in[f] ;
       if( V<K ) ERROR_exit("'%s' has only %d voxels with a nonzero W",fn,V) ;
       s->vox=(int *)malloc(sizeof(int)*(size_t)V) ; s->nvox_mask=V ; W=(float *)malloc(sizeof(float)*(size_t)V*K) ;
       for( q=0,f=0 ; f<nvg ; f++ ) if( in[f] ){
         s->vox[q]=f ; for( c=0 ; c<K ; c++ ) W[(size_t)q*K+c]=fb[c][f] ;
         q++ ;
       }
       for( c=0 ; c<K ; c++ ) mri_free(bim[c]) ;
       free(bim) ; free(fb) ; free(in) ; DSET_delete(wds) ;
     }
     srm_data_censor(&oq,&d1,NULL) ;
     srm_data_read(&oq,&d1) ;
     srm_data_timeline(&oq,&d1,0) ;
     srm_data_normalize(&oq,&d1,pv->zscore,0,1) ;
     t=d1.T ;
     if( pv->zscore==SRM_ZSCORE_NONE ){                 /* no mu in the files: center on the new data */
       mu=(float *)malloc(sizeof(float)*(size_t)V) ;
       for( f=0 ; f<V ; f++ ){ double a=0.0 ; int z ; for( z=0 ; z<t ; z++ ) a+=s->x[(size_t)f*t+z] ; mu[f]=(float)(a/t) ; }
     }
     mv=SRM_model_from_w(V,K,W,mu) ;
     if( mv==NULL ) ERROR_exit("out of memory applying W for Subj %s",s->label) ;
     P=(float *)malloc(sizeof(float)*(size_t)K*(size_t)t) ;
     rc=SRM_project(mv,0,s->x,t,P) ;
     if( rc!=SRM_OK ) ERROR_exit("Subj %s: projection failed: %s",s->label,SRM_strerror(rc)) ;
     s->nvox=V ;
     if( !o->no_projected )
       srm_write_proj(o,s,P,K,t,"subject's NEW data in the shared space, through this subject's own W from the earlier run: one row per retained time point, one column per component",NULL,argc,argv) ;
     SRM_free_model(mv) ; free(W) ; free(mu) ; free(P) ; free(s->x) ; s->x=NULL ;
     srm_prog_set(&bar,ss+1,NULL) ;
   }
   srm_prog_end(&bar) ;
   SRM_free_model(m0) ; SRM_provenance_free(pv) ;
}

int main( int argc , char **argv )
{
   srm_opts o ;
   int nopt=1 , ii , jj , rr , ss ;
   THD_datatable *tab=NULL ; THD_datatable_index *dx=NULL ;
   srm_data d ;
   int nsub , nrun , have_run , ngroup ; char **glabel=NULL ;
   srm_subj *sb ;
   int *runlen=NULL , Tin=0 , *c0in=NULL , T=0 , *seglen=NULL ;
   int icol_mask , icol_group ;
   byte *censor_keep=NULL ;
   int *time_index=NULL ; char **runname=NULL ;
   int K=0 , do_cv=0 ;

   /*-- startup, the AFNI way: the program is named once, here --*/

   if( argc<2 ){ usage_3dSRM() ; return 0 ; }
   mainENTRY("3dSRM main") ; machdep() ; AFNI_SETUP_OMP(0) ;
   AFNI_logger("3dSRM",argc,argv) ; PRINT_VERSION("3dSRM") ; AUTHOR("P Molfese") ;
   THD_check_AFNI_version("3dSRM") ;

   srm_opts_default(&o) ; memset(&d,0,sizeof(d)) ;

   /*-- options --*/

   while( nopt<argc && argv[nopt][0]=='-' ){
#define NEEDARG(what) do{ if( ++nopt>=argc ) ERROR_exit("need %s after '%s'",(what),argv[nopt-1]) ; }while(0)
#define IS(s) (strcasecmp(argv[nopt],(s))==0)
     if( IS("-help") || IS("-h") ){ usage_3dSRM() ; return 0 ; }
     else if( IS("-prefix") ){ NEEDARG("a prefix") ; free(o.prefix) ; o.prefix=strdup(argv[nopt++]) ; }
     else if( IS("-nfeatures") ){
       NEEDARG("an integer or comma list") ;
       o.klist=srm_parse_intlist(argv[nopt],&o.nk) ;
       if( o.klist==NULL ) ERROR_exit("-nfeatures needs positive integers separated by commas, not '%s'",argv[nopt]) ;
       nopt++ ;
     }
     else if( IS("-algorithm") ){
       NEEDARG("prob or det") ;
       if( strcasecmp(argv[nopt],"prob")==0 ) o.par.algo=SRM_ALGO_PROB ;
       else if( strcasecmp(argv[nopt],"det")==0 ) o.par.algo=SRM_ALGO_DET ;
       else ERROR_exit("-algorithm must be prob or det, not '%s'",argv[nopt]) ;
       o.given_algo=1 ; nopt++ ;
     }
     else if( IS("-niter") ){ NEEDARG("an integer") ; o.par.niter=atoi(argv[nopt++]) ; o.given_niter=1 ;
       if( o.par.niter<1 ) ERROR_exit("-niter must be at least 1") ; }
     else if( IS("-tol") ){ NEEDARG("a number") ; o.par.tol=atof(argv[nopt++]) ; o.given_tol=1 ;
       if( !(o.par.tol>=0.0) ) ERROR_exit("-tol must be a nonnegative number") ; }
     else if( IS("-seed") ){ NEEDARG("an integer") ; o.par.seed=(unsigned int)strtoul(argv[nopt++],NULL,10) ; o.given_seed=1 ; }
     else if( IS("-nstart") ){ NEEDARG("an integer") ; o.nstart=atoi(argv[nopt++]) ;
       if( o.nstart<1 ) ERROR_exit("-nstart must be at least 1") ; }
     else if( IS("-mask") ){ NEEDARG("a dataset") ; o.maskname=argv[nopt++] ; o.given_mask=1 ; }
     else if( IS("-concat") ){ NEEDARG("a 1D list of run starts") ; o.concat=argv[nopt++] ; }
     else if( IS("-run_column") ){ NEEDARG("a name") ; free(o.run_col) ; o.run_col=strdup(argv[nopt++]) ; }
     else if( IS("-group_column") ){ NEEDARG("a name") ; free(o.group_col) ; o.group_col=strdup(argv[nopt++]) ; }
     else if( IS("-zscore") ){
       NEEDARG("run, all, or none") ;
       if( strcasecmp(argv[nopt],"run")==0 ) o.zscore=SRM_ZSCORE_RUN ;
       else if( strcasecmp(argv[nopt],"all")==0 ) o.zscore=SRM_ZSCORE_ALL ;
       else if( strcasecmp(argv[nopt],"none")==0 ) o.zscore=SRM_ZSCORE_NONE ;
       else ERROR_exit("-zscore must be run, all, or none, not '%s'",argv[nopt]) ;
       o.given_zscore=1 ; nopt++ ;
     }
     else if( IS("-censor") ){ NEEDARG("a 1D file") ; o.censor_name=argv[nopt++] ; }
     else if( IS("-zcensor") ){ o.zcensor=1 ; nopt++ ; }
     else if( IS("-zcensor_warn") ){ NEEDARG("a percentage") ; o.zcensor_warn=(float)atof(argv[nopt++]) ; o.zcensor_warn_given=1 ; }
     else if( IS("-run_coupling") ){
       NEEDARG("shared or a nonnegative number") ;
       if( strcasecmp(argv[nopt],"shared")==0 ) o.par.coupling=SRM_COUPLING_SHARED ;
       else { char *end ; double v=strtod(argv[nopt],&end) ;
         if( *end!='\0' || !isfinite(v) || v<0.0 ) ERROR_exit("-run_coupling must be shared or a number >= 0, not '%s'",argv[nopt]) ;
         o.par.coupling=v ; }
       o.given_coupling=1 ; nopt++ ;
     }
     else if( IS("-axis") ){
       NEEDARG("time or item") ;
       if( strcasecmp(argv[nopt],"item")==0 ) o.axis_item=1 ;
       else if( strcasecmp(argv[nopt],"time")==0 ) o.axis_item=0 ;
       else ERROR_exit("-axis must be time or item, not '%s'",argv[nopt]) ;
       nopt++ ;
     }
     else if( IS("-repeat_column") ){ NEEDARG("a name") ; free(o.repeat_col) ; o.repeat_col=strdup(argv[nopt++]) ; o.item_col_given=1 ; }
     else if( IS("-item_column") ){ NEEDARG("a name") ; free(o.item_col) ; o.item_col=strdup(argv[nopt++]) ; o.item_col_given=1 ; }
     else if( IS("-repeated_stimulus") ){ o.repeated_stimulus=1 ; nopt++ ; }
     else if( IS("-time_labels") ){
       if( nopt+2>=argc ) ERROR_exit("need a name and a file after -time_labels") ;
       o.tl_name=(char **)realloc(o.tl_name,sizeof(char *)*(o.ntl+1)) ; o.tl_file=(char **)realloc(o.tl_file,sizeof(char *)*(o.ntl+1)) ;
       o.tl_name[o.ntl]=argv[nopt+1] ; o.tl_file[o.ntl]=argv[nopt+2] ; o.ntl++ ; nopt+=3 ;
     }
     else if( IS("-model_in") ){ NEEDARG("a model file") ; o.model_in=argv[nopt++] ; }
     else if( IS("-update_S") ){ o.update_S=1 ; nopt++ ; }
     else if( IS("-apply_W") ){ NEEDARG("a prefix") ; o.apply_W=argv[nopt++] ; }
     else if( IS("-cv") ){
       NEEDARG("none, run, or loso") ;
       if( strcasecmp(argv[nopt],"none")!=0 && strcasecmp(argv[nopt],"run")!=0 && strcasecmp(argv[nopt],"loso")!=0 )
         ERROR_exit("-cv must be none, run, or loso, not '%s'",argv[nopt]) ;
       free(o.cvmode) ; o.cvmode=strdup(argv[nopt++]) ;
     }
     else if( IS("-cv_blocks") ){ NEEDARG("an integer") ; o.cv_blocks=atoi(argv[nopt++]) ; o.cv_blocks_given=1 ;
       if( o.cv_blocks<2 ) ERROR_exit("-cv_blocks must be at least 2") ; }
     else if( IS("-cv_metrics") ){ NEEDARG("a comma list") ; o.cv_metrics=argv[nopt++] ; }
     else if( IS("-match_window") ){ NEEDARG("an integer") ; o.match_window=atoi(argv[nopt++]) ; o.match_window_given=1 ;
       if( o.match_window<1 ) ERROR_exit("-match_window must be at least 1") ; }
     else if( IS("-loso_iter") ){ NEEDARG("an integer") ; o.loso_iter=atoi(argv[nopt++]) ;
       if( o.loso_iter<1 ) ERROR_exit("-loso_iter must be at least 1") ; }
     else if( IS("-fold_column") ){ NEEDARG("a name") ; o.fold_col=argv[nopt++] ; }
     else if( IS("-cv_adapt") ){ NEEDARG("a fraction") ; o.cv_adapt=atof(argv[nopt++]) ; o.cv_adapt_given=1 ;
       if( !(o.cv_adapt>0.0 && o.cv_adapt<1.0) ) ERROR_exit("-cv_adapt must be between 0 and 1") ; }
     else if( IS("-score_by") ){ NEEDARG("a comma list") ; o.score_by=argv[nopt++] ; }
     else if( IS("-cv_fit") ){ NEEDARG("pooled, within, or reference:LABEL") ; o.cv_fit=argv[nopt++] ; }
     else if( IS("-balance") ){ NEEDARG("an integer") ; o.balance=atoi(argv[nopt++]) ; }
     else if( IS("-reliability") ){ NEEDARG("a comma list") ; o.reliability=argv[nopt++] ; }
     else if( IS("-nsplit") ){ NEEDARG("an integer") ; o.nsplit=atoi(argv[nopt++]) ; }
     else if( IS("-nperm") ){ NEEDARG("an integer") ; o.nperm=atoi(argv[nopt++]) ; }
     else if( IS("-save_W") ){ o.save_W=1 ; nopt++ ; }
     else if( IS("-save_recon") ){ o.save_recon=1 ; nopt++ ; }
     else if( IS("-save_r2") ){ o.save_r2=1 ; nopt++ ; }
     else if( IS("-no_projected") ){ o.no_projected=1 ; nopt++ ; }
     else if( IS("-proj_format") ){
       NEEDARG("dset or 1D") ;
       if( strcasecmp(argv[nopt],"dset")==0 ) o.proj_dset=1 ;
       else if( strcasecmp(argv[nopt],"1D")==0 ) o.proj_dset=0 ;
       else ERROR_exit("-proj_format must be dset or 1D, not '%s'",argv[nopt]) ;
       nopt++ ;
     }
     else if( IS("-memory_limit") ){
       NEEDARG("GiB") ; o.memory_limit_gib=atof(argv[nopt++]) ; o.memory_limit_given=1 ;
       if( !(o.memory_limit_gib>0.0) || !isfinite(o.memory_limit_gib) ) ERROR_exit("-memory_limit must be a positive number of GiB") ;
     }
     else if( IS("-memory_override") ){ o.memory_override=1 ; nopt++ ; }
     else if( IS("-quiet") ){ o.quiet=1 ; nopt++ ; }
     else if( IS("-progress") ){
       NEEDARG("auto, bar, line, or off") ;
       if( strcasecmp(argv[nopt],"auto")==0 ) o.progress=PROG_AUTO ;
       else if( strcasecmp(argv[nopt],"bar")==0 ) o.progress=PROG_BAR ;
       else if( strcasecmp(argv[nopt],"line")==0 ) o.progress=PROG_LINE ;
       else if( strcasecmp(argv[nopt],"off")==0 ) o.progress=PROG_OFF ;
       else ERROR_exit("-progress must be auto, bar, line, or off, not '%s'",argv[nopt]) ;
       nopt++ ;
     }
     else if( IS("-show_table") ){ o.show_table=1 ; nopt++ ; }
     else if( IS("-dataTableFile") ){
       NEEDARG("a file") ;
       if( tab!=NULL ) ERROR_exit("the data table was given twice") ;
       tab=THD_read_datatable_file(argv[nopt++]) ;
     }
     else if( IS("-dataTable") ){
       int nused=0 ;
       NEEDARG("a table") ;
       if( tab!=NULL ) ERROR_exit("the data table was given twice") ;
       tab=THD_read_datatable_args(argc,argv,nopt,srm_options,&nused) ;
       nopt+=nused ;
     }
     else ERROR_exit("unknown option '%s'",argv[nopt]) ;
   }
   if( nopt<argc ) ERROR_exit("unexpected argument '%s'",argv[nopt]) ;
   if( tab==NULL ) ERROR_exit("a data table is required: -dataTable or -dataTableFile") ;
   if( o.show_table ){ THD_datatable_print(tab,stdout) ; return 0 ; }

   /*-- options that belong to one axis --*/
   if( o.axis_item && !o.given_zscore ) o.zscore=SRM_ZSCORE_ALL ;   /* betas of one GLM share a baseline */
   if( o.axis_item ){
     if( o.censor_name!=NULL )
       ERROR_exit("-censor drops time points; with -axis item, leave an item out by removing its rows (for every subject)") ;
     if( o.zcensor || o.zcensor_warn_given )
       ERROR_exit("-zcensor is for time points; with -axis item an all-zero item is always an error, named by Subj and item") ;
     if( o.cv_blocks_given )
       ERROR_exit("-cv_blocks cuts time into contiguous blocks; items have no order: the folds of -axis item are runs, or the levels of -fold_column") ;
     if( o.match_window_given )
       ERROR_exit("-match_window is the window of segment matching in time; with -axis item each item is identified on its own") ;
     if( o.ntl>0 )
       ERROR_exit("-time_labels label time points; with -axis item put the labels in a table column (for -score_by)") ;
     if( o.repeated_stimulus )
       ERROR_exit("-repeated_stimulus is for -axis time; with -axis item the repeats are the items that recur across runs") ;
     if( o.par.coupling!=SRM_COUPLING_SHARED )
       ERROR_exit("-run_coupling needs the same runs in the same columns for every subject, so it is for -axis time only") ;
   } else if( o.item_col_given ) ERROR_exit("-item_column and -repeat_column need -axis item") ;

   /*-- what this build does not do yet: say so, plainly --*/
   if( o.repeated_stimulus || o.ntl>0 )
     ERROR_exit("-repeated_stimulus and -time_labels are not available in this build yet") ;
   if( o.axis_item && (o.model_in!=NULL || o.apply_W!=NULL) )
     ERROR_exit("with -axis item, -model_in and -apply_W are not available in this build yet") ;
   if( o.fold_col!=NULL && !o.axis_item )
     ERROR_exit("-fold_column is for -axis item; with -axis time the folds are the runs (or -cv_blocks)") ;
   if( o.fold_col!=NULL && strcasecmp(o.cvmode,"none")==0 ) ERROR_exit("-fold_column needs -cv run or -cv loso") ;
   if( o.concat!=NULL && o.axis_item ) ERROR_exit("-concat splits a time series into runs; it is for -axis time") ;
   if( o.concat!=NULL && o.apply_W!=NULL ) ERROR_exit("-concat with -apply_W is not available in this build yet") ;
   if( o.cv_adapt_given ) ERROR_exit("-cv_adapt is not available in this build yet") ;
   if( o.score_by!=NULL || o.cv_fit!=NULL || o.balance>0 || o.reliability!=NULL || o.nsplit>0 || o.nperm>0 )
     ERROR_exit("-score_by, -cv_fit, -balance, -reliability, -nsplit and -nperm are not available in this build yet") ;
   if( o.update_S && o.model_in==NULL ) ERROR_exit("-update_S needs -model_in") ;
   if( o.model_in!=NULL && o.apply_W!=NULL ) ERROR_exit("-model_in and -apply_W cannot be combined") ;
   if( o.model_in!=NULL || o.apply_W!=NULL ){
     const char *who = o.model_in!=NULL ? "-model_in" : "-apply_W" ;
     const char *from = o.model_in!=NULL ? "the model" : "the earlier run's model" ;
     if( o.nk>0 )             ERROR_exit("-nfeatures is defined by %s; do not give it with %s",from,who) ;
     if( o.given_algo )       ERROR_exit("-algorithm is defined by %s; do not give it with %s",from,who) ;
     if( o.given_zscore )     ERROR_exit("-zscore is defined by %s; do not give it with %s",from,who) ;
     if( o.given_coupling )   ERROR_exit("-run_coupling is defined by %s; do not give it with %s",from,who) ;
     if( o.given_seed )       ERROR_exit("-seed is defined by %s; do not give it with %s",from,who) ;
     if( o.nstart!=1 )        ERROR_exit("-nstart is not used with %s: nothing is fitted from scratch",who) ;
     if( (o.given_niter || o.given_tol) && !o.update_S )
                              ERROR_exit("-niter and -tol are defined by %s; with %s they apply only to -update_S",from,who) ;
     if( o.model_in!=NULL ){
       if( o.censor_name!=NULL ) ERROR_exit("-censor is defined by the model; do not give it with -model_in") ;
       if( o.zcensor )           ERROR_exit("-zcensor is defined by the model; do not give it with -model_in") ;
       if( strcasecmp(o.cvmode,"loso")==0 )
         ERROR_exit("-cv loso needs a group fit to leave subjects out of; with -model_in use -cv run, which scores each subject on runs the model and their W never saw") ;
     } else {
       if( o.maskname!=NULL )    ERROR_exit("-apply_W takes the voxels from W: do not give -mask") ;
       if( strcasecmp(o.cvmode,"none")!=0 || o.cv_metrics!=NULL || o.save_r2 )
                                 ERROR_exit("-apply_W only projects: it takes no -cv, -cv_metrics or -save_r2") ;
       if( o.save_W || o.save_recon ) ERROR_exit("-apply_W only projects: it takes no -save_W or -save_recon") ;
     }
     if( o.model_in!=NULL ) srm_run_model_in(&o,tab,argc,argv) ; else srm_run_apply_W(&o,tab,argc,argv) ;
     if( !o.quiet ) INFO_message("done") ;
     return 0 ;
   }
   if( o.nk<1 ) ERROR_exit("-nfeatures is required") ;
   do_cv=(strcasecmp(o.cvmode,"none")!=0) ;
   if( o.nk>1 && !do_cv ) ERROR_exit("several -nfeatures values need -cv: they are compared on held-out data") ;
   if( !do_cv && (o.cv_metrics!=NULL || o.save_r2) ) ERROR_exit("-cv_metrics and -save_r2 need -cv run or -cv loso") ;
   for( ii=0 ; ii<o.nk ; ii++ ) if( o.klist[ii]>K ) K=o.klist[ii] ;       /* K is the largest until the loop below */
   if( !isfinite(o.zcensor_warn) || o.zcensor_warn<0.0f || o.zcensor_warn>100.0f ) ERROR_exit("-zcensor_warn must be between 0 and 100") ;
   if( o.zcensor_warn_given && !o.zcensor ) ERROR_exit("-zcensor_warn requires -zcensor") ;
   if( o.par.coupling!=SRM_COUPLING_SHARED && THD_datatable_column(tab,o.run_col)<0 )
     ERROR_exit("-run_coupling needs runs: the table has no '%s' column",o.run_col) ;

   /*-- subjects, runs, masks, and the censor mask --*/

   if( o.axis_item ){
     srm_data_items(&o,tab,&d,2) ;
     srm_data_shape_items(&o,&d,1,K) ;
   } else {
     srm_data_table(&o,tab,&d,2) ;
     srm_data_shape(&o,&d,1,K) ;
     srm_data_censor(&o,&d,NULL) ;
   }
   nsub=d.nsub ; nrun=d.nrun ; sb=d.sb ; Tin=d.Tin ; glabel=d.glabel ; ngroup=d.ngroup ; runname=d.runname ;
   have_run=d.have_run ; dx=d.dx ; icol_mask=d.icol_mask ; icol_group=d.icol_group ; runlen=d.runlen ; c0in=d.c0in ;

   /*-- do not clobber anything: check every output now, not after the work --*/

   if( !THD_ok_overwrite() ){
     char fn[THD_MAX_NAME] ;
     snprintf(fn,sizeof(fn),"%s.S.1D",o.prefix) ;   if( THD_is_file(fn) ) ERROR_exit("'%s' already exists (use -overwrite)",fn) ;
     snprintf(fn,sizeof(fn),"%s.srm",o.prefix) ;    if( THD_is_file(fn) ) ERROR_exit("'%s' already exists (use -overwrite)",fn) ;
     snprintf(fn,sizeof(fn),"%s.fit.txt",o.prefix) ; if( THD_is_file(fn) ) ERROR_exit("'%s' already exists (use -overwrite)",fn) ;
     snprintf(fn,sizeof(fn),"%s.objective.1D",o.prefix) ; if( THD_is_file(fn) ) ERROR_exit("'%s' already exists (use -overwrite)",fn) ;
     if( o.par.algo==SRM_ALGO_PROB ){ snprintf(fn,sizeof(fn),"%s.rho2.1D",o.prefix) ; if( THD_is_file(fn) ) ERROR_exit("'%s' already exists (use -overwrite)",fn) ; }
     if( do_cv ){
       snprintf(fn,sizeof(fn),"%s.scores.txt",o.prefix) ; if( THD_is_file(fn) ) ERROR_exit("'%s' already exists (use -overwrite)",fn) ;
       snprintf(fn,sizeof(fn),"%s.scores_folds.txt",o.prefix) ; if( THD_is_file(fn) ) ERROR_exit("'%s' already exists (use -overwrite)",fn) ;
       if( o.nk>1 ){ snprintf(fn,sizeof(fn),"%s.nfeat.txt",o.prefix) ; if( THD_is_file(fn) ) ERROR_exit("'%s' already exists (use -overwrite)",fn) ; }
     }
     if( o.par.coupling!=SRM_COUPLING_SHARED ){ snprintf(fn,sizeof(fn),"%s.coupling.txt",o.prefix) ; if( THD_is_file(fn) ) ERROR_exit("'%s' already exists (use -overwrite)",fn) ; }
   }
   srm_check_subject_outputs(&o,&d) ;

   /*-- the memory contract, before the big allocations --*/

   srm_memory_preflight(&o,&d,K,do_cv,(o.par.coupling!=SRM_COUPLING_SHARED) ? nrun : 0) ;

   /*-- read everything, then the retained timeline, then normalize --*/

   if( o.axis_item ){
     srm_data_read_items(&o,&d) ;
     srm_data_items_timeline(&d) ;
     if( d.T<K ) ERROR_exit("%d items; -nfeatures is %d",d.T,K) ;
   } else {
     srm_data_read(&o,&d) ;
     srm_data_timeline(&o,&d,K) ;
   }
   srm_data_normalize(&o,&d,o.zscore,K,0) ;
   T=d.T ; seglen=d.seglen ; time_index=d.time_index ; censor_keep=d.keep ;

   /*-- the fits, and cross-validation --*/
   {
     srm_ctx cx ; srm_want want ; srm_folds folds ; int ik , rc ;
     int *carry=NULL , ncarry=0 , icol_run=THD_datatable_column(tab,o.run_col) ;
     int *nvox_v=(int *)malloc(sizeof(int)*(size_t)nsub) ; float **xv=(float **)malloc(sizeof(float *)*(size_t)nsub) ;
     double *kbl=NULL ; FILE *kfp=NULL ; char **gl=glabel ;
     memset(&cx,0,sizeof(cx)) ;
     cx.o=&o ; cx.sb=sb ; cx.nsub=nsub ; cx.nrun=nrun ; cx.T=T ; cx.seglen=seglen ; cx.runname=runname ;
     cx.glabel=gl ; cx.ngroup=ngroup ; cx.quiet=o.quiet ; cx.progress=o.progress ;
     for( ss=0 ; ss<nsub ; ss++ ){ nvox_v[ss]=sb[ss].nvox ; xv[ss]=sb[ss].x ; }
     memset(&want,0,sizeof(want)) ; memset(&folds,0,sizeof(folds)) ;

     if( do_cv ){
       want=srm_parse_metrics(o.cv_metrics) ; want.ident=o.axis_item ;
       if( o.axis_item ) srm_make_item_folds(&o,&d,&folds) ;
       else              srm_make_folds(&folds,nrun,seglen,runname,T,o.cv_blocks) ;
       srm_check_fold_lengths(&folds,want,o.match_window) ;
       { int ff , big=0 ;
         for( ff=1 ; ff<folds.nfold ; ff++ ) if( folds.len[ff]>folds.len[big] ) big=ff ;
         if( T-folds.len[big]<K ) ERROR_exit("holding out fold %s leaves %d %s to fit; -nfeatures is %d",
                                             folds.name[big],T-folds.len[big],folds.unit,K) ; }
       if( !o.quiet ){
         if( o.axis_item ){
           char lv[256] ; int ff , pos=0 ;
           for( ff=0 ; ff<folds.nfold && pos<200 ; ff++ ) pos+=snprintf(lv+pos,sizeof(lv)-pos,"%s%s (%d)",ff?", ":"",folds.name[ff],folds.len[ff]) ;
           if( ff<folds.nfold ) snprintf(lv+pos,sizeof(lv)-pos,", ...") ;
           INFO_message("-cv %s: %d folds of the items, from column %s: %s",o.cvmode,folds.nfold,
                        o.fold_col!=NULL?o.fold_col:o.run_col,lv) ;
           ININFO_message("%s; scored by item identification (chance is 1 over the items in the fold)",
                          strcasecmp(o.cvmode,"loso")==0 ? "each subject scored by a model that never saw them" : "the group refitted without each fold") ;
         } else
         INFO_message("-cv %s: %d folds of the %s, %s",o.cvmode,folds.nfold,folds.by_blocks?"time points in blocks":"runs",
                      strcasecmp(o.cvmode,"loso")==0 ? "each subject scored by a model that never saw them" : "the group refitted without each fold") ;
         if( folds.by_blocks ) WARNING_message("one run only: the folds are contiguous time blocks. Time points next to a block's edge "
                                               "resemble the other side, so scores are somewhat optimistic; separate runs are better") ;
       }
       if( o.nk>1 ){
         kbl=(double *)calloc((size_t)o.nk*8,sizeof(double)) ;
       }
       carry=srm_carry_columns(tab,&d,&ncarry) ;
     }

     for( ik=0 ; ik<o.nk ; ik++ ){
       SRM_model *m=NULL ; SRM_fitinfo fi ; srm_prog bar ; fit_prog fp ; SRM_params par=o.par ;
       double *start_obj=(double *)calloc((size_t)o.nstart,sizeof(double)) ; int best=0 ;
       srm_cvres res ; int have_res=0 ;
       K=o.klist[ik] ; cx.K=K ;
       par.nfeat=K ; par.center=(o.zscore==SRM_ZSCORE_NONE) ;
       if( par.coupling!=SRM_COUPLING_SHARED ){ par.nseg=nrun ; par.seglen=seglen ; }
       cx.par=par ;
       if( !o.quiet ){
         INFO_message("%sfitting a %s SRM: %d subjects, %d components, %d %s, up to %d iterations%s",
                      o.nk>1?"[-nfeatures ] ":"",par.algo==SRM_ALGO_PROB?"probabilistic":"deterministic",nsub,K,T,
                      d.axis_item?"items":"time points",par.niter,
                      par.coupling!=SRM_COUPLING_SHARED?", a W per run tied to the subject's anchor":"") ;
         if( par.coupling!=SRM_COUPLING_SHARED ) ININFO_message("it starts from the shared-W fit, then adds the per-run mappings") ;
       }
       srm_prog_begin(&bar,o.progress,o.quiet,"fit","iterations",par.niter) ;
       fp.bar=&bar ; fp.istart=0 ; fp.nstart=o.nstart ; fp.first_obj=0.0 ; fp.tag=NULL ;
       par.progress=srm_fit_progress ; par.progress_ctx=&fp ;
       rc=SRM_fit_multistart(nsub,T,nvox_v,xv,&par,o.nstart,&m,&fi,start_obj,&best) ;
       srm_prog_end(&bar) ;
       if( rc!=SRM_OK ) ERROR_exit("the fit failed: %s",SRM_strerror(rc)) ;
       if( !o.quiet ){
         INFO_message("%s after %d iterations (objective %.10g)",fi.converged?"converged":"stopped at the iteration limit",
                      fi.niter_done,fi.objective[fi.niter_done]) ;
         if( !fi.converged ) ININFO_message("not converged to -tol %.3g: the fit still depends on the seed; raise -niter or relax -tol",par.tol) ;
         if( o.nstart>1 ){
           double lo=start_obj[0] , hi=start_obj[0] ;
           for( ii=1 ; ii<o.nstart ; ii++ ){ if( start_obj[ii]<lo ) lo=start_obj[ii] ; if( start_obj[ii]>hi ) hi=start_obj[ii] ; }
           ININFO_message("best of %d starts: start %d; objectives span %.3g relative",o.nstart,best+1,
                          fabs(hi-lo)/(fabs(start_obj[best])>0.0?fabs(start_obj[best]):1.0)) ;
         }
       }
       rc=SRM_canonicalize(m) ;
       if( rc!=SRM_OK ) ERROR_exit("could not put the model in its canonical gauge: %s",SRM_strerror(rc)) ;

       /*-- held-out scoring --*/
       if( do_cv ){
         have_res=1 ;
         srm_cvres_alloc(&res,nsub,folds.nfold,!o.no_projected && o.nk==1,o.save_r2 && o.nk==1,&cx) ;
         if( strcasecmp(o.cvmode,"loso")==0 ) srm_cv_loso(&cx,m,&folds,want,&res) ;
         else                                srm_cv_by_runs(&cx,m,&folds,want,&res) ;
         srm_cv_summary(&cx,&res,want,o.cvmode,folds.nfold) ;
         if( o.nk>1 ){                                    /* one row per K for the comparison table */
           int met , mi[4]={M_ISC,M_MATCH,M_CHANCE,M_R2} ;
           for( met=0 ; met<4 ; met++ ){
             double s=0.0 , s2=0.0 ; int cnt=0 ;
             for( ss=0 ; ss<nsub ; ss++ ){ double v=srm_subj_metric(&res,ss,mi[met]) ; if( isfinite(v) ){ s+=v ; s2+=v*v ; cnt++ ; } }
             kbl[(size_t)ik*8+2*met]=cnt>0 ? s/cnt : NAN ;
             kbl[(size_t)ik*8+2*met+1]=(cnt>1) ? sqrt(fmax(0.0,(s2-s*s/cnt)/(cnt-1))/cnt) : NAN ;
           }
         }
       }
       if( o.nk>1 ){                                     /* comparing K only: no model outputs */
         if( have_res ) srm_cvres_free(&res) ;
         SRM_free_fitinfo(&fi) ; SRM_free_model(m) ; free(start_obj) ;
         continue ;
       }

       /*-- outputs --*/
       {
         char fn[THD_MAX_NAME] , cm[256] ; FILE *fq ;
         snprintf(fn,sizeof(fn),"%s.S.1D",o.prefix) ;
         if( d.axis_item ){
           char *ck ;
           snprintf(cm,sizeof(cm),"shared response: %d items x %d components, ordered by variance",T,K) ;
           ck=srm_keys_comment(cm,d.item_key,T) ;
           srm_write_transposed_1D(fn,m->s,K,T,ck) ; free(ck) ;
         } else {
           snprintf(cm,sizeof(cm),"shared response: %d time points (retained) x %d components, ordered by variance",T,K) ;
           srm_write_transposed_1D(fn,m->s,K,T,cm) ;
         }
         if( !o.quiet ) INFO_message("wrote %s",fn) ;

         snprintf(fn,sizeof(fn),"%s.objective.1D",o.prefix) ;
         fq=fopen(fn,"w") ; if( fq==NULL ) ERROR_exit("cannot write '%s'",fn) ;
         fprintf(fq,"# objective at the start of each iteration, and at the end (%s)\n",
                 par.algo==SRM_ALGO_PROB?"log-likelihood, rises":"summed squared error, falls") ;
         for( ii=0 ; ii<=fi.niter_done ; ii++ ) fprintf(fq,"%.12g\n",fi.objective[ii]) ;
         fclose(fq) ;

         if( par.algo==SRM_ALGO_PROB ){
           snprintf(fn,sizeof(fn),"%s.rho2.1D",o.prefix) ;
           fq=fopen(fn,"w") ; if( fq==NULL ) ERROR_exit("cannot write '%s'",fn) ;
           fprintf(fq,"# noise variance per subject, in the order of the table:") ;
           for( ss=0 ; ss<nsub ; ss++ ) fprintf(fq," %s",sb[ss].label) ;
           fputc('\n',fq) ;
           for( ss=0 ; ss<nsub ; ss++ ) fprintf(fq,"%.9g\n",m->rho2[ss]) ;
           fclose(fq) ;
         }

         snprintf(fn,sizeof(fn),"%s.fit.txt",o.prefix) ;
         fq=fopen(fn,"w") ; if( fq==NULL ) ERROR_exit("cannot write '%s'",fn) ;
         fprintf(fq,"# 3dSRM fit: %s, %d components, %d subjects, %d %s (%d runs), seed %u\n",
                 par.algo==SRM_ALGO_PROB?"probabilistic":"deterministic",K,nsub,T,d.axis_item?"items":"time points",nrun,par.seed) ;
         fprintf(fq,"# iterations %d, converged %d, tol %.3g, final objective %.12g\n",fi.niter_done,fi.converged,par.tol,fi.objective[fi.niter_done]) ;
         fprintf(fq,"# starts %d, kept start %d, objectives",o.nstart,best+1) ;
         for( ii=0 ; ii<o.nstart ; ii++ ) fprintf(fq," %.12g",start_obj[ii]) ;
         fputc('\n',fq) ;
         if( d.axis_item )
           fprintf(fq,"# axis item: %d item columns; zscore %s\n",T,o.zscore==SRM_ZSCORE_RUN?"run (over the items of each subject's own runs)":
                                                                        o.zscore==SRM_ZSCORE_ALL?"all (over all of each subject's items)":"none") ;
         else
           fprintf(fq,"# zscore %s, censored %d of %d input time points\n",o.zscore==SRM_ZSCORE_RUN?"run":o.zscore==SRM_ZSCORE_ALL?"all":"none",Tin-T,Tin) ;
         fprintf(fq,"Subj\tGroup\tFeatures_in_mask\tFeatures_used\tFeatures_dropped\trho2\n") ;
         for( ss=0 ; ss<nsub ; ss++ )
           fprintf(fq,"%s\t%s\t%d\t%d\t%d\t%s\n",sb[ss].label,glabel[sb[ss].group],sb[ss].nvox_mask,sb[ss].nvox,sb[ss].ndrop,
                   par.algo==SRM_ALGO_PROB ? (snprintf(cm,sizeof(cm),"%.9g",m->rho2[ss]),cm) : "NA") ;
         fclose(fq) ;

         if( par.coupling!=SRM_COUPLING_SHARED ){
           snprintf(fn,sizeof(fn),"%s.coupling.txt",o.prefix) ;
           fq=fopen(fn,"w") ; if( fq==NULL ) ERROR_exit("cannot write '%s'",fn) ;
           fprintf(fq,"# deviation of each run's W from the subject's anchor: |W_ir - W_i|^2 / %d (0 identical, 2 unrelated); lambda %.6g\n",K,par.coupling) ;
           fprintf(fq,"Subj\tGroup\tRun\tDeviation\n") ;
           for( ss=0 ; ss<nsub ; ss++ ) for( rr=0 ; rr<nrun ; rr++ ){
             double dv ; SRM_run_deviation(m,ss,rr,&dv) ;
             fprintf(fq,"%s\t%s\t%s\t%.6g\n",sb[ss].label,glabel[sb[ss].group],runname[rr],dv) ;
           }
           fclose(fq) ;
         }
         if( have_res ){
           srm_write_scores(&cx,&folds,&res,want,o.cvmode,tab,carry,ncarry) ;
           if( !o.quiet ) INFO_message("wrote %s.scores.txt and %s.scores_folds.txt",o.prefix,o.prefix) ;
         }
       }

       /* the model: what a later run needs to take in new subjects */
       {
         char fn[THD_MAX_NAME] ; SRM_provenance *pv=SRM_provenance_new() ;
         snprintf(fn,sizeof(fn),"%s.srm",o.prefix) ;
         pv->zscore=o.zscore ;
         if( d.axis_item ){
           pv->axis=SRM_AXIS_ITEM ; pv->nitem=T ; pv->item_key=(char **)malloc(sizeof(char *)*(size_t)T) ;
           for( ii=0 ; ii<T ; ii++ ) pv->item_key[ii]=strdup(d.item_key[ii]) ;
         } else {
           pv->axis=SRM_AXIS_TIME ; pv->ntime_in=Tin ; pv->nseg=nrun ;
           pv->seglen=(int *)malloc(sizeof(int)*(size_t)nrun) ; memcpy(pv->seglen,seglen,sizeof(int)*(size_t)nrun) ;
           pv->keep=(unsigned char *)malloc((size_t)Tin) ; memcpy(pv->keep,censor_keep,(size_t)Tin) ;
         }
         pv->nsubj=nsub ; pv->label=(char **)malloc(sizeof(char *)*(size_t)nsub) ;
         for( ss=0 ; ss<nsub ; ss++ ){
           char sn[128] ; srm_safe(sn,sizeof(sn),sb[ss].label) ; pv->label[ss]=strdup(sn) ;
         }
         rc=SRM_save_model(fn,m,pv) ;
         if( rc!=SRM_OK ) ERROR_exit("cannot write the model '%s': %s",fn,SRM_strerror(rc)) ;
         SRM_provenance_free(pv) ;
         if( !o.quiet ) INFO_message("wrote %s",fn) ;
       }

       /* per subject: the data in the shared space (held out, under -cv), and the maps */
       srm_write_subject_outputs(&o,&d,m,0,have_res?&res:NULL,K,argc,argv) ;
       if( have_res ) srm_cvres_free(&res) ;
       SRM_free_fitinfo(&fi) ; SRM_free_model(m) ; free(start_obj) ;
     }

     /*-- comparing numbers of components --*/
     if( o.nk>1 ){
       char fn[THD_MAX_NAME] ; int best_isc=-1 , best_match=-1 , best_r2=-1 ;
       snprintf(fn,sizeof(fn),"%s.nfeat.txt",o.prefix) ;
       kfp=fopen(fn,"w") ; if( kfp==NULL ) ERROR_exit("cannot write '%s'",fn) ;
       fprintf(kfp,"# held-out scores by number of components (%s; mean over %d subjects, with the standard error across subjects)\n",o.cvmode,nsub) ;
       fprintf(kfp,"nfeatures\tisc\tisc_se\t%s\t%s_se\tchance\tchance_se\tr2\tr2_se\n",want.ident?"ident":"match",want.ident?"ident":"match") ;
       for( ik=0 ; ik<o.nk ; ik++ ){
         int q ;
         fprintf(kfp,"%d",o.klist[ik]) ;
         for( q=0 ; q<8 ; q++ ){ double v=kbl[(size_t)ik*8+q] ; if( isfinite(v) ) fprintf(kfp,"\t%.7g",v) ; else fprintf(kfp,"\tNA") ; }
         fputc('\n',kfp) ;
         if( want.isc   && (best_isc<0   || kbl[(size_t)ik*8+0]>kbl[(size_t)best_isc*8+0]) )   best_isc=ik ;
         if( want.match && (best_match<0 || kbl[(size_t)ik*8+2]>kbl[(size_t)best_match*8+2]) ) best_match=ik ;
         if( want.r2    && (best_r2<0    || kbl[(size_t)ik*8+6]>kbl[(size_t)best_r2*8+6]) )    best_r2=ik ;
       }
       fclose(kfp) ;
       if( !o.quiet ){
         INFO_message("wrote %s",fn) ;
         if( best_isc>=0 )   INFO_message("best -nfeatures by held-out isc: %d",o.klist[best_isc]) ;
         if( best_match>=0 ) INFO_message("best -nfeatures by %s: %d",want.ident?"item identification":"segment matching",o.klist[best_match]) ;
         if( best_r2>=0 )    INFO_message("best -nfeatures by held-out r2: %d",o.klist[best_r2]) ;
         INFO_message("no model was written; rerun with the -nfeatures you choose for the model and its maps") ;
       }
     }
     free(nvox_v) ; free(xv) ;
     if( do_cv ) srm_free_folds(&folds) ;
   }
   if( !o.quiet ) INFO_message("done") ;
   return 0 ;
}
