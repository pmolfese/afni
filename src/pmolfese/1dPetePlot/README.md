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
writes version 2 recipes with independent source recipes and explicit column roles;
version 1 projects still open. Saving after Undo excludes undone imports. PNG/SVG
replay of a mixed figure works in a headless build. Adding graphs grows the export
height when needed; that export preference remains separate from presentation Undo.
An empty workspace can be saved/reopened but cannot produce an image until it has a graph.

Interactive limits: regular files up to 16 MiB, 100,000 rows, 1,000,000 cells, and
256 columns; 16 selected series per line/scatter/histogram, 50,000 selected points
per line/scatter preview, and square heatmaps up to 256 × 256. The session retains
at most 64 graphs including Undo history. No preview downsampling is performed;
histogram counts use all selected observations. Exceeding limits gives an explicit
message; CLI use remains available for larger tables within its existing limits.

This slice does not infer units, ROI identities, condition meanings, run/TR timing,
uncertainty/group roles, or cross-file alignment. X-matrix attachment suggestions,
CSV, long/advanced RSA schemas, and a simultaneous multi-thumbnail candidate gallery
remain later work. Ordinary ROI-wide RSA import is supported below. Existing CLI X-matrix shading remains available.

The included `examples/imported-figure.peteplot.json` demonstrates three plot types
from a named TSV, with editable labels/captions and reproducible export:

```sh
./target/release/1dPetePlot -recipe examples/imported-figure.peteplot.json
./target/release/1dPetePlot -recipe examples/imported-figure.peteplot.json -prefix imported.png
```

## RSA presets and semantic import (M7 first slice)

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
| `-rsa-plots histogram,scatter,heatmap` | Plot families, in requested order |
| `-rsa-pair A,B` | Scatter X/Y models, both in the selected model set |
| `-missing gap` | Omit nonfinite histogram values and incomplete scatter pairs; missing heatmap cells are gray |

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
that an effect heatmap is an RDM. Unknown schema columns, commonality outputs,
long time/series/run tables, and malformed inputs fail with an explanation.
Input limits are 16 MiB, 10,000 ROIs, 512 columns, and 1,000,000 cells.

**Open Data / drag-drop also recognizes RSA tables.** The full-window review offers
one model distribution, a same-kind two-model scatter, or an effect heatmap per Add.
Choose r/b/pr effect columns and optionally enter comma-separated ROI IDs. Other
columns remain visible as metadata and cannot become plotted effects accidentally.
The ordinary `.1D` line/matrix workflow remains available for non-RSA inputs.
Interactive RSA tables currently allow at most 256 columns. All previews run through
the existing background worker, confirmation, source-change checks, and Undo.

Saving a project freezes **RSA preset version 1**, model/effect/ROI/plot choices,
and the comparison pair alongside source fingerprints and presentation edits.
Reopening or headless `-recipe` replay does not reselect future defaults. RSA graphs
can coexist with ordinary imported graphs, reorder, undo, and save like other panels.
No new runtime dependencies or external services were added.

Try the clearly synthetic fixtures `examples/rsa/ordinary.rsa.1D` and
`examples/rsa/regression.rsa.1D`. Actual saved RDMs still use `-heat`; a dedicated
RDM/condition-label workflow, asymmetric confidence-interval plots, long schemas,
and fitted/contrast effect visualization remain subsequent M7 slices.

## Editable figures and projects

Click a figure title, graph title, axis label, or existing caption to edit it.
**Enter** commits a single-line field; **Done** commits a multiline caption;
**Escape** cancels. **Titles & arrangement** also provides keyboard-accessible
editing buttons, including adding an empty figure/graph caption. Text uses the
same renderer in the window, PNG, and SVG. Titles/captions wrap using the bundled
font's metrics; if the text leaves insufficient plotting space, the edit stays
uncommitted with an explanation. Increase export height or shorten the text.
Long axis labels scale to fit. Margin-free QC strips omit all editable text.

**Undo/Redo** covers committed text edits, panel order, column layout, and Add Graph. It
preserves each panel's identity and current view when reordering. Pan/zoom and
export settings are separate from presentation undo. History is session-local
(up to 100 presentation edits). An unsaved indicator appears in the window title
and inspector; opening another project or closing the window offers Save/Discard/
Cancel when necessary.

**Save Project…** writes a versioned `.peteplot.json` recipe containing the original
plot options, edits, panel identities/order, all panels' view bounds, export size,
and Current view/Full plot preference. **Open Project…** reconstructs the figure.
Selecting a single graph in the sidebar does not remove the other graphs from
the saved project. **Save Image…** remains a separate operation.

| Action | macOS | Linux |
|---|---|---|
| Save Image | Command–S | Ctrl–S |
| Save Project | Command–Shift–S | Ctrl–Shift–S |
| Open Project | Command–O | Ctrl–O |
| Open Data | Command–I | Ctrl–I |
| Undo / Redo | Command–Z / Command–Shift–Z | Ctrl–Z / Ctrl–Shift–Z |

Project files reference inputs with paths relative to the project where possible;
keep their relative arrangement when moving the project/data folder. Byte
fingerprints detect changed inputs, including masks, bands, x coordinates, order
files, and X matrices. Missing or changed inputs stop opening/saving the project
with a source-specific error. They never silently reload different observations.
The fingerprint is a change detector, not a security signature. A GUI relink flow
and portable bundles containing data are future work. Inline `1D:` inputs and
demos are self-contained; stdin must be saved to a file and reopened before it
can be saved as a project. Recipes use the current executable's plotting algorithms.

Replay a saved project without a display, or open the included editable example:

```sh
./target/release/1dPetePlot -recipe examples/edited-figure.peteplot.json
./target/release/1dPetePlot -recipe examples/edited-figure.peteplot.json -prefix figure.png
./target/release/1dPetePlot -recipe examples/edited-figure.peteplot.json -prefix figure.svg
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
missing pairs, rectangular heatmap coordinates, semantic import, Add/Undo/reopen,
headless exports, and preset freezing with identical replay pixels. Import tests cover queued files, retry/cancel, stale results and changed sources,
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
pending. See the adjacent
`../1dPetePlot_PLAN.md` for the milestone status.

## Deliberately outside this prototype

Full AFNI option parity, supplied-height `-hist`, arbitrary expression transforms,
normalization/detrending, multiple x columns, full metadata handling,
PDF/JPEG/PNM/PostScript, configurable fonts, and large-data rendering optimizations.
C's run-relative `-CENSORTR`/`-concat` and multiple censor-color layers remain
future work. Unsupported flags fail explicitly. AFNI build/distribution and QC
HTML integration remain separate from the implemented plotting interfaces.
