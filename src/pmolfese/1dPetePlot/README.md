# 1dPetePlot — native and headless RSA plotting prototype

A Rust prototype for AFNI-style plotting with native desktop interaction and
CPU-rendered PNG/SVG export. The built-in gallery contains **synthetic data**.

## Try it

From this directory, with Rust/Cargo 1.95 or newer:

```sh
cargo run --locked --release -- -demo
```

Drag inside a panel to pan. Hold **Shift** and drag a rectangle to zoom; release
to apply it. Double-click a panel or press **R** to reset. Scrolling does not zoom.
Hover for sample or matrix-cell values. A scrollable left sidebar shows thumbnails
and titles; select one graph or **All graphs**. Drag thumbnails onto another card
to move a graph to that position (the others shift, preserving their relative order).
Right-click a graph or thumbnail to **Copy this graph**, **Copy all graphs**, or
save either scope. Clipboard images include axes and labels, at export resolution.
**Titles & arrangement** edits the figure title, graph titles, and column count;
it also provides Earlier/Later buttons for reordering without dragging. The export controls save the current view by default; select
**Full plot** to restore the original data ranges for export. Export dimensions
are independent of the window and display scaling. **Save Image…** opens the OS
file chooser. Select PNG or SVG in the inspector before opening it. The dialog
remembers the last destination during this session. Cancelling leaves files alone.
A replacement confirmation may follow the platform's own overwrite prompt.

In full screen, an **Exit full screen** button remains visible in the inspector.
**Escape** leaves full screen without closing the app. The View menu also toggles
full screen: **Control–Command–F** on macOS, **F11** on Linux. Standard native
window decorations remain enabled outside full screen.

An output filename saves and exits. Add `-show` to also open the native window:

```sh
cargo run --locked --release -- -demo -prefix gallery.png
cargo run --locked --release -- -demo -prefix gallery.svg -overwrite
cargo run --locked --release -- -demo heat -show
```

The demo scenes have fixed data/drawing settings, but accept titles and panel layout options. Use the example `.1D` files to
exercise the CLI drawing options. `-help` lists the implemented options.

A local macOS app can be built with `./scripts/package_macos.sh`; its output
is `target/1dPetePlot.app`. It is a local development bundle, not a signed or
notarized distribution. Opening it starts an empty import workspace; `-demo` still opens the synthetic gallery. The ordinary
`target/release/1dPetePlot` binary accepts data files and works from the shell.

## Open / Drop → Preview → Add

Launch with no arguments for an empty workspace, or import into an existing figure.
Choose **Open Data…** (Command–I / Ctrl–I) or drop local `.1D` / headered `.tsv`
files onto the window. Import review fills the main window, temporarily replacing
the figure workspace. Try `examples/import-demo.tsv` (explicitly synthetic).

1. Inspect column names, numeric/text classification, missing counts, and the first
   eight source rows. Text columns remain visible; they are not coerced into numbers.
2. Choose **Lines**, **Distribution**, **Scatter**, or **Heatmap**. The explanation
   states what each view means. Select numeric X/Y columns; Lines can use row index.
   Heatmap requires a square numeric selection and does not imply an RDM.
3. Choose the missing-value policy: reject (default), allow gaps/omit missing from
   counts, or explicitly impute column means. X coordinates must remain finite.
4. Click **Add graph to figure** in the fixed footer. Previewing does not alter
   the figure. **Cancel import** or Escape discards the pending previews;
   **Reload file / Retry** rereads a changed source, and **Skip file** advances a queue.
   The preview and column controls scroll while the action buttons stay visible.
   Cancel or Escape returns to the figure; the review page resizes with the window.

Up to 16 files can be queued, each with its own preview and Add decision. Import
another view of the same file to compare different plot types. Graphs from different
files can coexist; overlay/pairing within this slice uses columns of one table,
without guessing cross-file joins. Headers become series labels. A single worker
per window parses and renders previews in the background; obsolete results are
discarded. Fingerprints are checked before/after reading and again before Add.

Add Graph is undoable, including returning to an empty figure. Sources and original
column indices have stable identities through reordering and Undo/Redo. Save Project
writes version 4 recipes with independent source recipes, explicit column roles,
style overrides, and duplicate identities; versions 1–3 still open.
Saving after Undo excludes undone imports. PNG/SVG replay of a mixed figure works in a headless build. Adding graphs grows the export
height when needed; that export preference remains separate from presentation Undo.
An empty workspace can be saved/reopened but cannot produce an image until it has a graph.

Interactive limits: regular files up to 16 MiB, 100,000 rows, 1,000,000 cells, and
256 columns; 16 selected series per line/scatter/histogram, 50,000 selected points
per line/scatter preview, and square heatmaps up to 256 × 256. The session retains
at most 64 graphs including Undo history. No preview downsampling is performed;
histogram counts use all selected observations. Exceeding limits gives an explicit
message; CLI use remains available for larger tables within its existing limits.

Import does not infer units, ROI identities, condition meanings, run/TR timing,
or cross-file alignment. Scientific roles and attachments are explicit (below).
CSV, automatic attachment suggestions, cross-file joins, and a simultaneous
multi-thumbnail candidate gallery remain later work. Wide and long RSA tables
are supported through the explicit scientific families below.

The included `examples/imported-figure.peteplot.json` demonstrates three plot types
from a named TSV, with editable labels/captions and reproducible export:

```sh
./target/release/1dPetePlot -recipe examples/imported-figure.peteplot.json
./target/release/1dPetePlot -recipe examples/imported-figure.peteplot.json -prefix imported.png
```

## Scientific roles and attachments

In Open/Drop preview, use **Scientific roles & attachments** after choosing X/Y.
For new imports, the target is **this preview, added as a new graph**. For an
existing imported line/scatter/distribution graph, select it and choose
**Edit → Edit Scientific Attachments** (also in its right-click menu). This opens
the same full-window role/attachment review with **Apply changes to graph** in
the fixed footer. Nothing changes until Add or Apply.

- **Group column:** any original text or numeric category column, distinct from
  X/Y. Groups appear in first-seen order with consistent category colors across Y
  series. Lines break between groups; scatter pairs retain their row positions;
  grouped histograms use common bin edges. No aggregation or interpolation occurs.
  Missing group labels and groups without finite observations are rejected. The
  limit is 16 total series (Y columns × groups), with 50,000 stored line/scatter
  points in interactive previews. No mean imputation across groups.
- **Supplied uncertainty (lines):** choose absolute ± width columns, a width file,
  or lower/upper endpoint columns for each Y. A shared width column may be assigned
  deliberately to multiple series. External width files need one numeric column
  per selected Y, in the same order. Endpoints must bracket the supplied estimate;
  asymmetry is retained. No SE/CI estimation or confidence level is inferred.
  Include bands in automatic Y limits is on by default in import. Gap mode breaks
  bands where bounds are incomplete, preserving observed estimates; mean imputation
  is unavailable with these roles.
- **Censor column or file (lines):** exactly 0 = censored and 1 = retained, with
  one value per original sample. Censoring adds shading and hover metadata; it
  does not delete observations. Choose a column or an external single-column mask.
- **AFNI X matrix (lines):** choose an unselected X.mat.1D and optionally set a
  threshold percentage and condition colors (space-separated names or hex colors).
  Positive condition groups supply labels/colors; baseline and motion are excluded.
  NRowFull and GoodList must match the full plotted timeline. Censored matrix rows
  stay unshaded. This depicts modeled response above a fraction of each condition's
  peak magnitude, **not stimulus onset timing**. Censor and condition shading can
  coexist in multiple colors.

Sidecars must be regular files up to 16 MiB. The preview checks counts and metadata;
you must also confirm that sidecars use this graph's full row order/timeline.
Matching lengths alone cannot establish alignment. X must be finite and increasing
for timeline shading. There is no resampling or inferred cross-file join.
Changing a sidecar path clears that confirmation. Clear scientific roles resets
attachments; changing plot families preserves choices and reports incompatibilities.
Reload/Retry preserves mappings when the original file path and column names/order
are unchanged; a changed schema starts a new mapping.

Source paths/selectors, original column indices, grouping and attachment options
are retained in the version 4 recipe. All auxiliary files are fingerprinted,
checked before/after preview and again on Add/save/reopen, and rebased with the
project. Add/Undo, duplication, styling, and headless replay retain the attachments.
Editing attachments preserves the graph identity, titles, captions, axes, and
compatible styling. Use Reset styling if changed series roles make old styles
incompatible; review retained legend names when the scientific meaning changes.
The edited graph fits its new data bounds. Sidecar timeline confirmation is
required again. Apply is one Undo action, including the data and source recipe;
existing duplicates remain independent. Removing a sidecar removes its replay
dependency unless another graph still uses it. Save/reopen and headless replay
retain the edited recipe. No dependencies were added.

This editor currently requires the stored column roles of an ordinary graph added
through Open Data/drop. CLI-created base graphs and RSA/RDM graphs do not yet have
this attachment editor; import an ordinary table to assign editable roles. Separate
external lower/upper files are not editable here; use same-table endpoint columns
or a symmetric width file. Changing source column names/order requires reimport.

Open `examples/attachments/signals.tsv` for a small synthetic example: X = time,
Y = effect, lower/upper = lower/upper, Group = condition, and Censor = keep.
Attach `examples/attachments/X.mat.1D` for condition shading. A complete replayable
example is included:

```sh
./target/release/1dPetePlot -recipe examples/attachments/attached-figure.peteplot.json
./target/release/1dPetePlot -recipe examples/attachments/attached-figure.peteplot.json -prefix attached.png

# The same scientific mechanics are available directly from the CLI:
./target/release/1dPetePlot -infiles 'examples/attachments/signals.tsv[effect]' \
  -xfile 'examples/attachments/signals.tsv[time]' -one -group-by condition \
  -yfiles_lo 'examples/attachments/signals.tsv[lower]' \
  -yfiles_hi 'examples/attachments/signals.tsv[upper]' -ylim_use_pm \
  -censor examples/attachments/censor.1D \
  -shade-xmat examples/attachments/X.mat.1D -prefix attached.svg
```

`-group-by COLUMN` accepts an original name or zero-based index and requires one
regular data file and `-one`; row selectors and transpose are rejected for grouped
inputs. `-yfiles_lo` / `-yfiles_hi` accept matching file/selector lists, one endpoint
column per Y; they cannot be combined with `-yfiles_pm`. CLI grouping also rejects
mean imputation. The existing ungrouped `-yfiles_pm` behavior is unchanged.

## RSA presets, confidence intervals, and semantic import (M7)

`-rsa` reads an ordinary ROI-wide `3dRSA` results table directly: the commented
`#ROI label nvox ...` header, ROI IDs/labels, and named model effects. It does not
need a `grep`/`awk` extraction step. With no output, the summary opens in the
native window; an output path renders without a display/GPU.

```sh
./target/release/1dPetePlot -rsa results.rsa.1D
./target/release/1dPetePlot -rsa results.rsa.1D -prefix rsa-summary.png
./target/release/1dPetePlot -rsa results.rsa.1D \
  -rsa-models semantic,visual -rsa-rois 3,12,27 \
  -rsa-plots scatter,heatmap -rsa-pair semantic,visual -prefix comparison.svg
./target/release/1dPetePlot -rsa regression.rsa.1D -rsa-effect pr -prefix partial-r.png
./target/release/1dPetePlot -rsa results.rsa.1D \
  -rsa-plots intervals -rsa-ci boot -prefix confidence.svg
```

Default composition is one distribution per selected model, one scatter of the
first two models, and a **ROI × model effect heatmap** for up to 40 selected ROIs.
All models are selected by default, up to eight; larger model sets require an
explicit `-rsa-models` selection. Every selected ROI contributes to distributions
and paired scatter when finite. Above 40 ROIs, default composition omits the
heatmap and says so in the caption; it never samples distributions or silently
truncates rows. Explicit heatmap requests above that limit fail with instructions
to select ROIs. Export height grows to fit ROI labels when size is automatic.
Explicit sizes that cannot fit the labels fail before writing a file.

| Option | Meaning |
|---|---|
| `-rsa-models A,B` | Model names, in requested order; up to 8 |
| `-rsa-effect auto\|r\|b\|pr` | Default r for ordinary comparison; b for regression; pr selects partial r explicitly |
| `-rsa-rois 3,12,27` | Unique ROI IDs in requested order; default all, sorted by ID |
| `-rsa-plots histogram,scatter,heatmap,intervals` | Plot families, in requested order; intervals are opt-in |
| `-rsa-ci auto\|boot\|cboot\|dual` | Interval method; auto requires exactly one method common to the selected models |
| `-rsa-pair A,B` | Scatter X/Y models, both in the selected model set |
| `-missing gap` | Omit nonfinite histogram values and incomplete scatter pairs; gray heatmap cells; incomplete interval triplets retain a labeled missing row |

Missing selected effects fail by default. RSA presets reject imputation. Scatter
pairs use unique ROI IDs; repeated or absent requested IDs fail explicitly. Hover
shows the ROI ID, full label, voxel count, and effect values. Labels containing
spaces and model names containing underscores are supported. r, beta, and partial r
are distinct roles; a plot cannot silently compare them as equivalent effects.
`r` is the source's RSA comparison coefficient, not a promise of Pearson correlation;
the supplied metric description is retained in the default caption.

The reader recognizes p/q/FWE, bootstrap/condition-bootstrap/dual-bootstrap bounds,
fitted metrics/weights, contrasts, leave-one-out metrics, reliability, and noise-ceiling
columns as auxiliary data. They are **not** treated as effects or automatically
plotted. Probability ranges and supplied interval ordering are validated. There is
no significance filtering, inference recomputation, RDM reconstruction, or assumption
that an effect heatmap is an RDM. Fitted/commonality quantities require their
explicit families below. Unknown schema columns and malformed inputs fail with
an explanation. Long tables have the separate
identity-aware workflow below.
Input limits are 16 MiB, 10,000 ROIs, 512 columns, and 1,000,000 cells.

**Confidence-interval plots:** `-rsa-plots intervals` makes one horizontal ROI
effect plot per selected model, with points, asymmetric whiskers, a zero reference,
and a common x range covering all selected effects and endpoints. Combine families
with a comma-separated list, e.g. `-rsa-plots intervals,heatmap`. Up to 40 selected
ROIs fit each interval panel; automatic height grows to keep labels readable.
ROI order follows `-rsa-rois`, or sorted ROI IDs by default. Long labels are
shortened on the axis; hover retains the full label, ROI ID, nvox, estimate,
endpoints, method and effect scale.

`boot` reads `MODEL_bootLo/bootHi` (subject bootstrap), `cboot` reads
`MODEL_cbootLo/cbootHi` (condition bootstrap), and `dual` reads
`MODEL_dualLo/dualHi` (dual bootstrap). If more than one method is available,
choose one explicitly. Every selected model must have the requested pair of
columns. Confidence level and percentile/corrected-variance t construction are
read from matching source comments when recorded. Missing level metadata is
labeled **CI level not recorded**; 95% is never assumed.

Endpoints are already on the plotted effect scale: no symmetrization,
back-transformation, clamping, or inference recomputation. Regression bounds
belong to **beta**, so `-rsa-effect pr` with intervals fails rather than borrowing
beta bounds. Finite lower bounds must not exceed upper bounds; an estimate need
not lie inside its interval. Missing effects or endpoints fail by default.
With `-missing gap`, any incomplete triplet keeps its ROI row marked missing,
with neither a point nor a whisker; an entirely incomplete panel fails.

**Open Data / drag-drop also recognizes RSA tables.** The full-window review offers
one model distribution, a same-kind two-model scatter, an effect heatmap, or one
model's **Confidence intervals** per Add. Interval review includes the CI-method
selector and uses the same validation as the CLI.
Choose the scientific family, its compatible quantity columns, and optionally
enter comma-separated ROI IDs. Primary effects use r/b/pr; probability columns
remain metadata and cannot become plotted effects accidentally.
The ordinary `.1D` line/matrix workflow remains available for non-RSA inputs.
Interactive RSA tables currently allow at most 256 columns. All previews run through
the existing background worker, confirmation, source-change checks, and Undo.

Saving a project freezes **RSA preset version 1**, model/effect/ROI/plot choices,
the comparison pair and resolved CI method alongside source fingerprints and presentation edits.
Reopening or headless `-recipe` replay does not reselect future defaults. RSA graphs
can coexist with ordinary imported graphs, reorder, undo, and save like other panels.
No new runtime dependencies or external services were added.

Try the clearly synthetic fixtures `examples/rsa/ordinary.rsa.1D` and
`examples/rsa/regression.rsa.1D`. Saved matrices now have the `-rdm` workflow below;
long time/run schemas and fitted/commonality quantities are supported below.

## Long RSA tables (M7)

`-rsa` also recognizes the two long schemas emitted by `3dRSA.c`:

- **Model series:** `ROI label nvox time_index time_label effect stat ...`.
  Default: one curve per ROI, using numeric `time_index` order. Tick labels include
  the original `time_label`; these labels are never interpreted as seconds or
  sorted alphabetically. The logical model name is `series`. `stat` is auxiliary.
- **Run results:** `ROI label nvox summary model effect ...`, or joint-regression
  `beta partial_r ...`. Default: one categorical point plot per ROI, with a color
  for each selected model. `MEAN`, run labels, and named contrasts remain distinct
  summaries; no lines connect them and no averaging is performed. Source contrast
  descriptions are retained in the default caption. Beta and partial r stay separate.

```sh
# Curves with the source's supplied subject-bootstrap band
./target/release/1dPetePlot -rsa examples/rsa/series.rsa.1D \
  -rsa-ci boot -prefix series.png
# Select by time_index, retaining source index order; horizontal CI alternative
./target/release/1dPetePlot -rsa examples/rsa/series.rsa.1D \
  -rsa-rois 3 -rsa-times 0,2,5,10 -rsa-plots intervals -rsa-ci boot -prefix time-ci.svg
# Categorical run estimates; explicit summary and model order
./target/release/1dPetePlot -rsa examples/rsa/runs-joint.rsa.1D \
  -rsa-models visual,semantic -rsa-summaries task,rest,task_minus_rest \
  -rsa-effect pr -prefix runs.png
```

Long time tables support `-rsa-plots curves,heatmap,intervals`; run tables support
`points,heatmap`. Each requested mode makes one panel per selected ROI. Heatmaps
have models as rows and time indices/summaries as columns; they are **not RDMs**.
Axes and heatmap ranges are independently fitted per ROI. Curves omit CI bands
unless a supplied method is explicitly requested with `-rsa-ci boot|cboot|dual`.
Interval plots resolve `auto` only when exactly one supplied method exists. CI
level/construction comes from source metadata, as with ROI-wide intervals.
`-rsa-effect auto` resolves to `effect`, or `b` for joint-run tables. `effect` keeps
the source scale without promising a particular correlation metric.

`-rsa-times` selects existing nonnegative indices in index order;
`-rsa-summaries` selects literal names in requested order. Models default to their
first appearance in the file, summaries to first appearance, and ROIs to ID order.
Rows are keyed by ROI/model/time or ROI/model/summary. Duplicate keys, conflicting
ROI labels/voxel counts, conflicting time labels, invalid probabilities, reversed
CI endpoints, and unsupported schemas fail explicitly. A missing selected
combination or value fails by default. With `-missing gap`, absent combinations
remain gaps/gray cells; incomplete CI rows retain a missing marker. For curves,
an available estimate remains visible when only its band is missing. No
imputation, pooling across runs, significance filtering, or inference is added.

Open or drop the same `.1D` files to preview these plots. Select **one ROI per
Add**, an effect kind, model names (run tables), time indices or summary names,
and optionally a supplied CI band. The table inspector retains source row order;
full identities and auxiliary values are available on curve/point hover. Add,
Undo/Redo, styling, duplication, Save Project, and headless project replay use the
same data model. Recipes freeze resolved selections and fingerprint the source;
changed data is rejected rather than silently substituted.

Limits: 8 selected models, 40 ROIs, 2,000 time indices/summaries, 50,000 selected
observations, and 64 panels. Heatmaps allow 40 time/summary columns; interval panels
allow 40 rows. These are explicit selection limits, not automatic sampling.
The ordinary input bounds still apply. Synthetic fixtures are
`examples/rsa/{series,runs,runs-joint}.rsa.1D`. Cross-file ROI joins and additional
fixed-model contrast/leave-one-out views remain future work.

## Fitted models and commonality (M7)

Select a scientific family explicitly for wide `3dRSA` results:

| `-rsa-family` | Plotted quantities | Meaning |
|---|---|---|
| `effects` (default) | `MODEL_r`, `MODEL_b`, `MODEL_pr` | Existing primary-effect workflow |
| `fit` | `FIT_cvR` | Supplied cross-validated accuracy, on its source scale |
| `weights` | `FIT_w_COMPONENT` | Descriptive mean normalized fold weights |
| `fit-contrast` | `CONTRAST_cvDiff` | Supplied paired held-fold Fisher-z accuracy difference |
| `commonality` | `uniq_*`, `common_*` | Raw R² decomposition components, including negative shared components |
| `partial-r2` | `partialR2_*` | Conditional partial R², kept separate from raw decomposition |

Use `-rsa-measures COLUMN,COLUMN` to select 1–8 literal column names; omitted
selection uses all available quantities in that family within this limit. A
ROI-by-quantity heatmap is the default. `-rsa-plots histogram`, `scatter` (two
quantities), or `intervals` selects the existing geometry. Select ROIs with
`-rsa-rois`. Do not combine these families with primary-effect selectors
`-rsa-models`, `-rsa-effect`, or `-rsa-pair`.

No refitting, conversion, normalization, significance filtering, aggregation, or
clipping is performed. Negative shared components remain visible. Conditional
partial R² values are not additive raw decomposition regions; mean fold weights
are descriptive and are not inferential evidence. Source method metadata remains
in default captions. Probability columns are validated auxiliary data.

Commonality and partial-R² interval plots require supplied matching endpoint
columns and use the existing explicit `-rsa-ci boot|cboot|dual` selection and
recorded confidence-level rules. Fit accuracy, weights, and fit contrasts have no
CI endpoints in the writer and reject interval requests rather than estimating
uncertainty. Pairwise and triple-decomposition column identities remain intact.

```sh
./target/release/1dPetePlot -rsa examples/rsa/fitted-commonality.rsa.1D \
    -rsa-family fit -prefix accuracy.png
./target/release/1dPetePlot -rsa examples/rsa/fitted-commonality.rsa.1D \
    -rsa-family commonality -rsa-plots intervals -rsa-ci boot -prefix components.svg
./target/release/1dPetePlot -rsa examples/rsa/commonality-triple.rsa.1D \
    -rsa-family partial-r2 -prefix conditional.png
```

Open/drop these synthetic fixtures to select the same families and compatible
quantities in preview. Add, Undo, styling, duplication and project replay retain
explicit frozen selections. Ordinary fixed-model contrasts (`_diff`, `_rDiff`,
`_zDiff`), leave-one-out views, and cross-file ROI joins remain later work.

## Saved RDMs and condition labels

Use `-rdm` for actual saved square matrices, including numeric `.1D` files from
`3dRSA -save_rdm`. Provide a separate label file with **one unique condition name
per line, in the matrix's original row/column order**. Blank lines and `#` comments
are ignored; names may contain spaces. Names are case sensitive and leading/trailing
whitespace is trimmed. Numeric IDs can be names, but are never interpreted as indices.

```sh
# One matrix with labels; opens in the native window.
1dPetePlot -rdm model.1D -rdm-labels conditions.txt

# Same source condition order in both matrices: one shared label file.
1dPetePlot -rdm model.1D neural.1D -rdm-labels conditions.txt \
  -rdm-kind dissimilarity -rdm-scale shared -prefix comparison.svg

# Different source orders: repeat labels once per matrix. Align by name,
# then display both axes in the order given in display-order.txt.
1dPetePlot -rdm model.1D neural.1D \
  -rdm-labels model-conditions.txt -rdm-labels neural-conditions.txt \
  -rdm-order display-order.txt -rdm-scale shared -prefix aligned.png
```

The first label file determines the displayed order unless `-rdm-order FILE` is
given. Every matrix must have exactly the same unique condition names, even with
independent color scales. A display-order file is a full permutation of those
names, not a subset. Both axes and values move together; hover reports the full
row/column names and value. A single shared label file asserts that all input
matrices already have the same source order—unlabeled numbers cannot verify that
assertion. Label counts, duplicates, missing names, and mismatched sets fail before
export. No row-position join or RDM reconstruction from RSA effect tables occurs.

| Option | Meaning |
|---|---|
| `-rdm FILE ...` | 1–8 numeric square symmetric matrices, 2–256 conditions each |
| `-rdm-labels FILE` | Required once for shared source order, or repeat once per input matrix |
| `-rdm-order FILE` | Optional condition names in desired display order, one per line |
| `-rdm-kind unspecified\|similarity\|dissimilarity` | Declared matrix sense; default unspecified; changes labeling only |
| `-rdm-scale independent\|shared` | Independent per-matrix color limits by default; shared derives one range from all inputs |
| `-clim MIN:MAX` | Explicit fixed color limits applied to all input matrices |
| `-zerocenter` | Symmetric color limits around zero; cannot combine with `-clim` |
| `-cbar NAME` | Existing viridis/diverging palettes; default viridis |

AFNI's saved matrices can contain **similarities or dissimilarities**; neither
shape nor filename establishes their sense or units. Choose the sense explicitly
when known. Shared scales are opt-in: use them only for values you intend to
compare in the same units. They do not normalize, change signs, or convert `r`
to `1-r`. Automatic limits include every finite cell, including the diagonal.
Explicit limits can clip displayed colors while hover retains the original value.

Negative distances (including crossnobis estimates) and the original diagonal
are preserved. The reader checks symmetry with tolerance
`1e-6 × max(1, |a|, |b|)` without averaging values. Missing values fail by default;
`-missing gap` permits paired missing cells, displayed gray. One-sided missing
cells, infinities, entirely missing matrices, and imputation are rejected.
Numeric matrices may contain whitespace and comments; triangular vectors,
headered/text-bearing TSV matrices, AFNI selectors, and stdin are outside this
slice. Each matrix/label/order file is limited to 16 MiB; labels to 256 bytes.

**Open Data / drop:** open a numeric square `.1D` file, select **RDM / named
matrix**, then choose its condition-label file with the native file chooser.
Optionally attach a display-order file, select matrix sense, and set fixed color
limits. Preview validates the attachments before Add. Reuse identical fixed
limits across comparable imports; automatic shared ranges across multiple inputs
are currently a CLI operation. Imported graphs retain independent recipes.
Changing label/order files after preview prevents Add until you refresh the preview.

Save Project fingerprints the matrices **and their label/order files**, freezes
RDM preset version 1, scale/sense/palette choices, and retains edits and views.
Relative sources can move with the project; changed attachments are detected on
replay. Named axes get more room; long names are abbreviated and dense labels are
thinned to avoid overlap. Hover retains every full identity. PNG/SVG and native
views use the same scene, with no new dependencies.

Synthetic examples live in `examples/rdm/`: `model.1D` and `neural.1D` have
different source orders, described by `conditions.txt` and `neural-conditions.txt`.
`display-order.txt` provides a common grouped order. These are demonstration
values, not experimental findings.

## Editable figures and projects

Click a figure title, graph title, axis label, or existing caption to edit it.
**Enter** commits a single-line field; **Done** commits a multiline caption;
**Escape** cancels. **Titles & arrangement** also provides keyboard-accessible
editing buttons, including adding an empty figure/graph caption. Text uses the
same renderer in the window, PNG, and SVG. Titles/captions wrap using the bundled
font's metrics; if the text leaves insufficient plotting space, the edit stays
uncommitted with an explanation. Increase export height or shorten the text.
Long axis labels scale to fit. Margin-free QC strips omit all editable text.

**Undo/Redo** covers committed text/style edits, panel order, column layout, Add Graph,
and Duplicate Graph. It
preserves each panel's identity and current view when reordering. Pan/zoom and
export settings are separate from presentation undo. History is session-local
(up to 100 presentation edits). An unsaved indicator appears in the window title
and inspector; opening another project or closing the window offers Save/Discard/
Cancel when necessary.

**Style graph…** opens a full-window workspace with a live rendered preview and
fixed **Apply / Cancel** footer. Choose **This graph** for series colors, legend
names, line widths, marker radii, series-legend visibility, or heatmap palettes.
Expand each series to edit it independently. Unchecked settings retain the original
plot defaults; **Restore original graph styling** resets the draft. Line markers
represent isolated observations; changing their radius does not add markers along
continuous curves. Histograms expose series colors/names; confidence intervals
expose color, whisker width, and estimate radius. Sizes are logical drawing units
that scale with export resolution. Heatmap palette changes retain color limits.

**All graphs · shared appearance** has separate settings, initially unchecked.
Checked widths/radii replace corresponding per-series overrides in compatible
graphs; legend visibility and heatmap palettes also apply only where supported.
Series colors and names stay individual. Apply commits the whole change as one
Undo action. Cancel/Escape leaves the figure unchanged. Data, censor/shading keys,
reference labels, CI endpoints, and condition order are preserved.

**Duplicate graph** inserts an independent copy immediately after the selected
graph, keeping its data, annotations, styling, and current view. The copy is
selected for editing. Use the inspector, right-click a graph, or choose the Edit
menu (Command–D / Ctrl–D). Duplicates use stable IDs; copying a copy, reordering,
Undo/Redo, and headless replay retain the intended source and presentation. Image
height grows when needed to fit the additional graph; export size is outside
Undo history, as with Add Graph. Margin-free QC strips cannot be duplicated into
a multi-graph figure. The session limit remains 64 retained graphs including Undo.

**Save Project…** writes a versioned `.peteplot.json` recipe containing the original
plot options, edits, panel identities/order, all panels' view bounds, export size,
style overrides, copy identities, and Current view/Full plot preference. Version 4
is written; versions 1–3 remain readable with their original styling.
**Open Project…** reconstructs the figure.
Selecting a single graph in the sidebar does not remove the other graphs from
the saved project. **Save Image…** remains a separate operation.

| Action | macOS | Linux |
|---|---|---|
| Save Image | Command–S | Ctrl–S |
| Save Project | Command–Shift–S | Ctrl–Shift–S |
| Open Project | Command–O | Ctrl–O |
| Open Data | Command–I | Ctrl–I |
| Duplicate Graph | Command–D | Ctrl–D |
| Undo / Redo | Command–Z / Command–Shift–Z | Ctrl–Z / Ctrl–Shift–Z |

Project files reference inputs with paths relative to the project where possible;
keep their relative arrangement when moving the project/data folder. Byte
fingerprints detect changed inputs, including masks, bands, x coordinates, order
files, and X matrices. Normal replay and saving reject missing/changed inputs.
They never silently substitute observations. The fingerprint is a change detector,
not a security signature. Portable bundles containing data remain future work. Inline `1D:` inputs and
demos are self-contained; stdin must be saved to a file and reopened before it
can be saved as a project. Recipes use the current executable's plotting algorithms.

**File → Reload / Relink Sources** shows source status and replacement paths.
Keep a path to reload a rerun analysis, or use the native file chooser to relink a
moved file. **Preview changes** rebuilds in the background and reports point/cell
counts, ranges and changed series labels. Imported column names/order must match;
changed graph/series roles require reimport. Existing selectors and frozen RSA
selections are retained. Titles, captions, styles and panel identity/order survive
compatible reloads. New bounds are fitted by default; Keep current views is explicit.
**Apply changes** rechecks all files and commits one Undo action. A source changed
after preview requires a new preview. Cancel leaves the figure untouched.

If Open Project cannot reconstruct a source, it offers the same recovery workspace;
the currently open document remains intact until the replacement passes validation
and is applied. Save Project after recovery records the new fingerprints. GUI
recovery does not weaken strict headless `-recipe` replay.

Data Undo/Redo retains up to eight prior data versions in memory alongside the
ordinary 100-entry presentation history; older history requiring evicted versions
is pruned. Undo restores the old observed values and source recipe, not disk files.
You may still export that cached figure, but saving its recipe requires those old
source bytes to be available and to match. Undo history is not serialized. Reload
accepts regular sources up to 64 MiB each; narrower table/sidecar limits still apply.

Replay a saved project without a display, or open the included editable example:

```sh
./target/release/1dPetePlot -recipe examples/edited-figure.peteplot.json
./target/release/1dPetePlot -recipe examples/edited-figure.peteplot.json -prefix figure.png
./target/release/1dPetePlot -recipe examples/edited-figure.peteplot.json -prefix figure.svg
./target/release/1dPetePlot -recipe examples/styled-copies.peteplot.json
```

Replay exports **all panels**, with saved views (or full bounds if Full plot was
saved). Output/size options and `-show` are accepted after `-recipe FILE`; drawing
options cannot override saved edits. CLI replacements still require `-overwrite`.
**Copy replay command** becomes available after saving the current project.
Captions are also available without a project: `-caption TEXT` and repeatable
`-panel-caption N TEXT`, with original 1-based panel numbering before reordering.

## Real file inputs

```sh
# Overlay two columns with explicit legends.
./target/release/1dPetePlot -one -infiles examples/signals.1D \
  -legend_labels 'Model A' 'Model B' -xlabel 'Time (s)' \
  -title 'Model evidence' -width 1600 -height 1000 -prefix lines.png

# Actual binning of observations, with zero and highlighted reference lines.
./target/release/1dPetePlot -histogram -histbin 16 -density \
  -vline 0 -vline_hl 0.4 'Reference 0.4' \
  -infiles 'examples/effects.1D[0]' -prefix histogram.png

# Matched model effects, equal physical axis scaling, identity reference.
./target/release/1dPetePlot -scatter -diagonal \
  -xfile 'examples/effects.1D[0]' -infiles 'examples/effects.1D[1]' \
  -xlabel 'Model A (r)' -ylabels 'Model B (r)' -prefix scatter.svg

# Apply the same permutation to both matrix axes.
./target/release/1dPetePlot -heat -infiles examples/rdm.1D \
  -reorder examples/order.1D -zerocenter -cbar Reds_and_Blues_Inv \
  -figsize 8 8 -dpi 150 -prefix rdm.png
```

### QC, uncertainty, and missing observations

These commands use synthetic fixtures and follow the Python-style commands
emitted by `afni_proc.py` QC. Use PNG instead of its current JPEG filenames:

```sh
./target/release/1dPetePlot \
  -infiles 'examples/qc.tsv[enorm,outlier]' \
  -yfiles_pm 'examples/qc.tsv[se_enorm,se_outlier]' -ylim_use_pm \
  -censor_files examples/qc_censor.1D -censor_trs '0..2' \
  -censor_hline 0.2 0.05 -patches 60 60 -boxplot_on \
  -yaxis 0:0.6 0:0.1 -ylabels 'enorm (~mm)' 'outlier fraction' \
  -title 'Motion and outlier QC · synthetic data' -prefix qc.png

# Native interactive view of the same six motion parameters used by QC.
./target/release/1dPetePlot -infiles examples/qc_motion.1D \
  -sepscl -ylabels VOLREG -boxplot_on -patches 60 60 \
  -censor_files examples/qc_censor.1D

# Missing observations keep their time positions. Curves and bands break at gaps.
./target/release/1dPetePlot -infiles '1D: 0 1 NaN 3 2' \
  -missing gap -yfiles_pm '1D: 5@0.2' -ylim_use_pm -prefix gaps.svg
```

Run `./scripts/render_qc_gallery.sh` to generate reproducible QC PNG/SVG examples,
including a six-panel motion plot and the 1800 × 75 margin-free strip used above
APQC grayplots. See [M3_COMPATIBILITY.md](M3_COMPATIBILITY.md) for coverage and
remaining differences.

- `-censor_files` combines masks with `-censor_trs` by union. Masks contain
  **0 = censored, 1 = retained**, with exactly the selected data's row count.
  `-censor FILE` is a single-file alias. Curves retain censored samples; shading
  identifies them, and hover reports **CENSORED**. Censoring does not interpolate
  or delete observations. Index lists refer to selected rows, starting at zero.
- `-patches` gives run lengths (sum must match selected rows). Backgrounds
  alternate white/gray; an initial zero-length run starts with gray. Censor and
  run spans use sample midpoints, including explicit increasing x coordinates.
  Non-increasing x values with these annotations are rejected.
- `-censor_RGB` accepts a color with optional opacity, or three/four RGB(A)
  values in 0..1. Quote values together or separate them. Default pink opacity
  is 0.45. Censor and run shading are behind the data.
- `-boxplot_on` adds summaries: BC before censoring and AC after censoring.
  `-bplot_view BC_ONLY|AC_ONLY` chooses one; without censoring, one All box is
  shown. Quartiles use linear percentiles; whiskers end at the outermost observed
  values within 1.5 IQR. Outliers are circles; labels include finite sample counts.
  An empty AC sample is explicitly labeled empty. Boxplots summarize the entire
  selected series, even after viewport pan/zoom. They require separate line panels.
- `-yfiles_pm` supplies absolute symmetric half-widths, one column per selected
  y series with matching rows. Automatic limits include the bands only with
  `-ylim_use_pm`; explicit bounds take precedence. Negative half-widths use their
  absolute values, following Python behavior.
- `-scale` accepts one factor or one per series, `SCALE_TO_HLINE`, or
  `SCALE_TO_MAX` (maximum, **not** maximum absolute value). Curves, band widths,
  and horizontal thresholds scale together. Explicit axis limits do not scale.
  Zero/missing normalization denominators fail. Thresholds may be negative or
  zero; `NONE` explicitly omits a series' threshold.
- `-missing error` is the default. `-missing gap` allows NaN, NA, N/A, and empty
  TSV fields. Lines/bands break, scatter/histograms omit missing observations,
  boxplots exclude them, and heatmaps show gray cells. All-missing plotted
  series/matrices fail. Explicit x values and censor/reorder files remain strict.
  `-missing mean` imputes each selected output column from its finite values,
  after selectors and transpose and before scaling; it reports replacements to
  stderr. Entirely missing columns cannot be imputed. No implicit zero filling.
- Inline `1D:` accepts comma/space-separated values, positive `N@value` repeats,
  and `|` to start another series. A trailing apostrophe transposes after row/column
  selection. Quote the whole specification, e.g. `"matrix.1D[0,1]'"`.
- TSV input uses a required, unique header row, numeric indices or exact header
  names in selectors, and skips text columns when no selector is given. Selected
  text columns fail explicitly. Numeric indices retain original table positions.
  TSV headers do not automatically become legends or heatmap condition labels yet.
- Scientific labels accept Unicode Greek and brace-delimited `r^{2}` / `β_{model}`.
  Ordinary underscores remain literal. This is basic formatting, not a TeX engine.

### Multiple shading colors and AFNI design matrices

Repeat `-shade START:END` for independent x-coordinate intervals. Color, label,
and alpha options modify the **most recent** range; a new range gets a fresh
palette color and default opacity 0.18. End must exceed start. Ranges apply to
all panels, may overlap, and are clipped to the current view without changing
axis limits. They never censor samples or alter BC/AC boxplot statistics.

```sh
./target/release/1dPetePlot -infiles examples/qc_motion.1D \
  -shade 10:20 -shade-color blue -shade-label 'Task A' \
  -shade 35:45 -shade-color orange -shade-label 'Task B' -shade-alpha 0.25

./target/release/1dPetePlot -infiles 'examples/qc_motion.1D[0]' \
  -shade-xmat examples/X.demo.mat.1D -shade-xmat-colors blue orange \
  -censor_files examples/qc_censor.1D -boxplot_on -prefix conditions.png
```

Manual shading works with line, scatter, and histogram data plots. Labeled ranges
get color keys even without `-legend_on`; repeated identical label/color pairs
share one key. Keys wrap rather than silently dropping conditions. Labels are
hidden in `-margin_off` strips. Annotation layers blend in command order, followed
by design conditions in group order, then run/censor backgrounds; data remain
on top. The existing `-censor_RGB` still controls the combined censor mask.

`-shade-xmat` requires an unselected AFNI matrix file and line plots. It reads
`ColumnGroups` and `ColumnLabels`; positive groups become conditions, while
baseline/motion/nuisance groups (<=0) are skipped. `StimLabels`/`StimBots`/`StimTops`
provide condition names when present, otherwise column names supply them. Basis
columns in a group are combined by logical OR. Standard X matrices **do not store
colors**; colors come from the default palette or `-shade-xmat-colors` (cycled).

This is **modeled-response shading, not exact stimulus timing**: a row is active
when any basis column's magnitude exceeds `-shade-xmat-threshold F` times the
condition's maximum absolute matrix value. Default F is 0.05; allowed range is
0 <= F < 1. Negative responses count, zero-only conditions produce no spans,
and overlapping conditions remain overlapping. Legends say `(model)` and a CLI
diagnostic states the threshold. Use explicit ranges for exact event intervals.

The plot must have the full `NRowFull` sample count. Censored matrices use
`GoodList` to place rows without compressing time; missing rows remain unshaded
because their design values are unknown. Prefer `X.nocensor.xmat.1D` to show model
activity through censored volumes. Uncensored matrices containing full rows are
used directly even if they retain a shorter GoodList. Condition boundaries use
midpoints of the plotted x coordinates. `RowTR` does not silently change units:
use `-dx TR` for seconds, or an explicit increasing `-xfile`. User-selected input
rows must still represent that full timeline.

`examples/X.demo.mat.1D` was generated from synthetic timings with:

```sh
3dDeconvolve -nodata 120 1 -polort 1 -num_stimts 2 \
  -stim_times 1 '1D: 12 54' 'BLOCK(8,1)' -stim_label 1 Task_A \
  -stim_times 2 '1D: 30 78' 'BLOCK(10,1)' -stim_label 2 Task_B \
  -x1D X.demo.mat.1D -x1D_stop
```

### Titles and panel arrangement

```sh
./target/release/1dPetePlot -demo \
  -title 'RSA model comparison' \
  -panel-title 1 'Evidence over time' \
  -panel-title 4 'Condition similarity' \
  -panel-order 4,1,2,3 -columns 2
```

`-title` names the whole figure. Repeat `-panel-title N TEXT` for individual
panels. Panel numbers start at **1**, refer to input-generated order, and are
resolved **before** reordering. `-panel-order` must contain each panel once;
`-columns` sets 1–4 columns (capped at the number of panels). Display fills rows
left to right. Do not combine `-panel-order` with `-reverse_order`. Smaller
column counts may need more image height; use `-height` if the layout is too short.
AFNI data selectors and heatmap permutation files remain **zero-based**.

GUI changes affect subsequent saved/copied images. They do not rewrite input
files or persist a session. Mixed plot types from arbitrary files and a reusable
figure recipe are still future work; the demo already illustrates mixed panels.

Run `./scripts/render_gallery.sh` after a release build to create PNG and SVG
examples under `gallery/`. Generated files and compiled binaries are ignored by
Git. The script intentionally overwrites its own named gallery outputs.

### Current behavior

- Input is numeric, whitespace-delimited `.1D`; `#` comments and scientific
  notation are accepted. Use positional files, `-infiles`, or `-yfiles`.
- Selectors: `[columns]` and `{rows}` support comma lists, inclusive ranges,
  reverse ranges, `$`, `$-N`, and positive steps such as `0..$(2)`.
  Selection preserves order and repeated indices. Quote selectors in the shell.
- Use `-stdin` or `-` to read standard input. It may be consumed only once.
- Multiple columns become separate panels; `-one` / `-one_graph` overlays them.
  First input is at the top. `-reverse_order` reverses panel order.
- Line/scatter panels share automatic ranges by default; `-sepscl` permits
  independent y scales. `-xaxis MIN:MAX` and `-yaxis MIN:MAX` set explicit limits
  for all panels. `-yaxis` also accepts one range per separate panel, `MIN:`,
  `:MAX`, and `:` for automatic endpoints. Options bind before panel reordering.
- Histograms use `max(8, floor(sqrt(N)))` bins by default. Internal edges belong
  to the bin on their right; the final edge is included in the final bin.
  Explicit widths produce uniform-width bins covering the data, so the final
  edge may exceed the maximum. Density integrates to one. Constant observations
  get a nonzero binning range. Overlay histograms currently bin each series
  independently; shared-bin overlays remain future work.
- Heatmaps require a square numeric matrix. First row is at the top. A reorder
  file must be a complete zero-based permutation. Axis labels track original
  indices after reordering. `-clim MIN:MAX` fixes the scale; `-zerocenter` makes it
  symmetric. Constant matrices receive a nondegenerate scale.
- The prototype provides sampled viridis and red/blue palettes, including the
  PR's aliases. They are not exact reproductions of every Matplotlib lookup table.
- Ragged rows and infinities fail. Missing observations require `-missing gap`
  or `-missing mean`; the default remains strict.
- PNG dimensions are exact for the whole figure. SVG has matching dimensions and
  scalable vector content. Default size starts at 1500 × 1200 and grows in height
  for additional panels; explicit dimensions remain exact. `-figsize` uses inches and
  defaults to 150 DPI. Pixel sizing cannot be combined with inches/DPI.
- Long headings/axis labels are abbreviated to fit this prototype. Supported
  image bounds are 400–8192 pixels wide, 300–8192 high, and at most 32 million
  pixels. `-margin_off` allows height down to 32 pixels for a single plot-only
  panel. Up to 64 panels are accepted within image/layout limits; many panels
  require multiple columns and explicit sizing.
- CLI outputs require `-overwrite` for replacement. GUI Save always confirms
  replacement. Rendering finishes into a temporary file in the destination
  directory before publishing, with an atomic no-clobber path when replacement
  is not permitted. Failed rendering or cancellation preserves an existing file.

## Renderer and native windows

The prototype chooses a **shared SVG scene** instead of separate Plotters and
GUI plot implementations. A small layout/drawing module emits vector geometry;
`usvg` resolves the bundled Ubuntu font into paths; `resvg`/`tiny-skia` rasterize
that same scene for both PNG and the native preview. Pan/zoom update shared data
bounds and re-render the scene. Hover uses the same coordinate transform.

This establishes one path for fonts, clipping, colors, labels, and annotations.
SVG remains vector, including text outlines, with no embedded plot screenshot.
Outlining avoids system-font dependencies but labels are no longer editable as
SVG text. See `assets/README.md` for font licensing and coverage.

The native shell is egui/eframe with wgpu (Metal on macOS, Vulkan on Linux).
There is no X11 window backend or XQuartz requirement. Linux native viewing
requires a working Wayland session and Vulkan graphics stack. The Linux
clipboard dependency still contains pure-Rust X11 protocol code transitively;
that is distinct from enabling an X11 window backend. Neither is present in the
headless build. Linux runtime portability has not yet been validated.

On macOS, `muda` supplies the native AppKit menu bar (File/Edit/View/Window) and
`rfd` supplies native Save dialogs. Linux uses an egui menu bar **inside the
window** and `rfd`'s XDG desktop portal backend, which selects the desktop's file
chooser. A supported GNOME/GTK/KDE portal backend must be available; rfd may use
Zenity as a fallback. No GTK widget toolkit is linked into our plotting shell.
The image clipboard uses `arboard`, with Wayland data-control support enabled.
Some compositors do not expose that protocol; arboard can fall back to X11
clipboard transport, and failures appear in the status line with a Save fallback.
Keep the app open until pasting on Linux, where clipboard ownership may be transient.

These desktop dependencies are optional and excluded with `--no-default-features`.
Rust itself provides no GUI widgets: egui draws the interior controls, while the
OS provides the outer window and the explicitly integrated menus/dialogs.
See [rfd's platform notes](https://docs.rs/rfd/0.17.2/rfd/) and
[muda's platform notes](https://docs.rs/muda/0.21.2/muda/).

CPU export does not initialize a window, graphics device, or event loop. For
compute nodes or machines without desktop libraries, build a headless binary:

```sh
cargo build --locked --release --no-default-features
./target/release/1dPetePlot -heat -infiles examples/rdm.1D -prefix rdm.svg
```

A headless build reports an actionable error if asked to open a window. Cargo
writes GUI and headless variants to the same target location by default; use
`--target-dir target/headless` to retain both.

## Build wrappers

These are standalone wrappers for the prototype. AFNI's top-level Make/CMake,
program lists, installation components, and `build_afni.py` are unchanged.

```sh
make                          # GUI release build
make GUI=0                    # headless release build
make install PREFIX=/tmp/peteplot-install

cmake -S . -B target/cmake -DPETEPLOT_GUI=OFF \
  -DCMAKE_INSTALL_PREFIX=/tmp/peteplot-install
cmake --build target/cmake
cmake --install target/cmake
```

CMake isolates Cargo outputs in its build tree and installs the executable into
`bin` by default (`PETEPLOT_INSTALL_BINDIR` is configurable). Cargo.lock is
included with the source. The actual development toolchain is Rust 1.99.1 on
macOS ARM64; the dependency baseline is Rust 1.95, not yet separately tested.
`--locked` fixes resolution, not network access; use `cargo fetch --locked` before
an offline build, then add `--offline`. A vendored release workflow is future work.

## Validation

```sh
cargo test --locked
cargo test --locked --no-default-features
cargo clippy --locked --all-targets -- -D warnings
cargo clippy --locked --no-default-features --all-targets -- -D warnings
cargo fmt --check
```

Tests cover selectors, malformed input, histogram counts/density, coordinate
transforms, real-file exports for all four plot types, dimensions, existing-file
protection, panel title/order validation, and egui hover/Shift-box-zoom/pan/reset,
right-click targeting, current-view saving, full-screen escape, and view identity
when rearranging panels. Document tests cover save/open with relocated relative
sources, changed/missing source detection, schema validation, NaN gap preservation,
identical rendered pixels after replay, atomic replacement, history branching,
and headless recipe export. GUI tests cover direct text hit regions at 1×/2×
display scaling, editing/reordering/undo/reopen, and rejecting captions that do not
fit. RSA tests cover ordinary/regression schemas, optional statistics, ROI order/
missing pairs, rectangular heatmap coordinates, asymmetric interval geometry,
CI method/level metadata, beta versus partial-r guards, missing triplets, 40-ROI
layout, semantic import, Add/Undo/reopen, headless exports, and preset freezing
with identical current-view replay pixels. Import tests cover queued files,
retry/cancel, stale results and changed sources,
missing-data choices, full-data counts, mixed-source replay, and clicking Add without
scrolling on an 850 × 500 viewport. A measured 25,001-point inspection + 900 × 600
preview took about 0.36 seconds in the macOS ARM64 debug test; this is a representative
measurement, not a universal performance guarantee. The egui interaction
harness runs in process without a desktop; it is not a substitute for a visual
native-window test.

Current platform evidence: macOS ARM64 builds and headless exports; generated
figures visually inspected. The revised macOS window, native menu bar, and native
Save dialog were inspected through desktop automation. Full-screen entry and
Escape restoration were verified; copied plot images opened successfully in
Preview. The thumbnail selection/reordering harness checks that other panels
retain their relative order. Linux Wayland and other
architectures remain unverified. For the new document/editing slice, desktop
automation timed out while attaching to the launched macOS app; native project
dialogs and text editing still need a hands-on check. In-process GUI tests and
visually inspected captioned exports pass. M6 builds/launches and mixed-source
exports are verified; desktop automation again timed out attaching to the app.
The user confirmed native drag/drop reaches the preview. That feedback exposed
scroll-hidden actions; the footer is now fixed and covered by a click test.
Import review now fills the main window, with tested resizing and Escape return.
Native verification of the revised footer, Open Data chooser, and Wayland drop remains
pending. The M7 confidence-interval and RDM slices pass 85 GUI-build tests and 60
headless-build tests, Clippy in both builds, and formatting checks. Correlation
and regression interval exports were visually inspected; the release app was
rebuilt. RDM tests also cover name-based axis alignment, matrix/label rejection, shared/fixed
color scales, attachment fingerprints, Add/Undo, and relocated pixel replay.
Hands-on native interval/RDM controls and Linux runtime checks remain pending.
The M5 styling/duplication slice passes 92 GUI-build tests and 63 headless-build
tests, Clippy in both builds, and formatting. Coverage includes scoped styling,
one-action Undo, a short-window Apply click, Escape precedence, copies of copies,
imported CI/RDM duplication, history pruning, unchanged legacy recipe pixels,
invalid-style/copy rejection, and exact replay pixels. Native styling controls
and Linux runtime still require hands-on verification. No dependencies were added.
The M6 scientific-attachment slice passes 101 GUI-build tests and 66 headless-build
tests, including combined bounds/group/censor/design previews, stale-sidecar Add
rejection, asymmetric and negatively scaled bounds, shared histogram bins,
Add/Undo/duplicate/replay, relocated source fingerprints, and headless PNG/SVG.
Native attachment choosers and Linux runtime still need hands-on verification.
The 2026-10-30 source/relink, existing-attachment and fitted/commonality slices pass
120 GUI-build tests (102 unit + 18 CLI) and 77 headless-build tests (58 unit + 19
CLI), Clippy with warnings denied in both configurations, and formatting. Coverage
includes stale preview rejection, schema-change rejection, source recovery, data
Undo/Redo with text edits, independent copies, removed sidecar dependencies,
version 4 replay, family separation and supplied 90% commonality bounds. Release
commonality heatmap/interval PNGs were visually inspected, and the macOS ARM64
binary/app were rebuilt offline. Native source/recovery and attachment-edit controls
still need hands-on verification; Linux remains unverified. No dependencies were
added. See the adjacent `../1dPetePlot_PLAN.md` for the milestone status.

## Deliberately outside this prototype

Full AFNI option parity, supplied-height `-hist`, arbitrary expression transforms,
normalization/detrending, multiple x columns, full metadata handling,
PDF/JPEG/PNM/PostScript, configurable fonts, and large-data rendering optimizations.
C's run-relative `-CENSORTR`/`-concat` and multiple censor-color layers remain
future work. Unsupported flags fail explicitly. AFNI build/distribution and QC
HTML integration remain separate from the implemented plotting interfaces.
