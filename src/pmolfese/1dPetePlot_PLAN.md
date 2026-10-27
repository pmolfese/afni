# 1dPetePlot implementation plan

Date: 2026-10-26

Status as of 2026-10-27: **M1 Complete (basic plotting prototype scope)**;
**M0, M2, M3, M4, M5, M6, and M7 Partially complete**. **M5 first slice Complete**
(direct editing, Undo/Redo, projects/replay). **M6 first slice Complete
(implementation and automated checks)**: Open/Drop → Preview → Add, empty launch,
column mapping, background previews, and mixed-source projects/replay. **M7 first
slice Complete (implementation and automated checks)**: ROI-wide RSA reader,
`-rsa` presets, semantic Open/Drop review, and reproducible projects/exports.
Native drag/drop reaches preview per user feedback; the revised fixed action footer
and native Open Data chooser still need hands-on verification.
Detailed completion markers appear
within each milestone below. Complete means the stated work is implemented and
verified to the extent explicitly recorded; it does not imply untested platform
support or full legacy compatibility.

Native GUI and headless builds, PNG/SVG exports, and interaction tests work on
macOS ARM64. The native macOS window, menus, Save dialog, full-screen escape, and
clipboard images in Preview were inspected. Linux runtime checks remain pending.

## Milestone handoff preference — Accepted 2026-10-27

After each implementation milestone or usable slice, update the completion markers
and always prompt the user with concrete options for the next implementation
milestone in this PLAN/ROADMAP. Recommend one option, explain its scope and any
prerequisites briefly, and let the user choose before starting the next milestone.
Do not end a milestone with only a completion report or silently start another.

The user reports that the M5 first slice looks to be working. Record this as user
feedback; detailed native-dialog and cross-platform verification remain separate.

### Next implementation choices — Awaiting user selection after M7 first slice

The user selected M7 on 2026-10-27; the first RSA slice is implemented below.

1. **M7 next slice: ROI effects with supplied confidence intervals (recommended).**
   Add ROI comparison/forest plots using the actual asymmetric bootstrap, condition-
   bootstrap, or dual-bootstrap bounds. Preserve interval method and effect scale;
   do not recompute inference or substitute symmetric bands. This builds on the
   semantic reader and makes the supplied uncertainty useful in figures.
2. **M5 next slice: Figure styling and duplication.** Editable series colors,
   line widths, marker sizes, and legend labels, plus Duplicate Graph. Keep these
   changes undoable and persistent, with deliberate figure/selected-graph scope.
3. **M6 next slice: Scientific attachments.** Add uncertainty/group/label roles
   for generic tables and compatible censor/X-matrix attachment previews. Validate
   timeline alignment and make the target graph explicit before attachment.

These choices do not supersede remaining M2 distribution or M3 parity work;
those remain available priorities if the user prefers them.

## Objective and agreed requirements

Build a Rust plotting program, provisionally named `1dPetePlot`, that covers
the static plotting mechanics of `1dplot` and `1dplot.py`, prioritizing the RSA
plots in [AFNI PR #919](https://github.com/afni/afni/pull/919).

- Preserve familiar AFNI command options where their semantics can be supported.
- Support direct command-line image generation with explicit dimensions, labels,
  axes, colors, reference lines, and other drawing options.
- Support native desktop windows, preferably through egui/eframe or a comparable
  Rust framework. X11 windows and XQuartz are outside the intended design.
- Export must work without a desktop session, display server, GPU, or GUI event
  loop. Saving a plot must not require opening a window or taking a screenshot.
- Keep scientific content consistent between window display and file export.
- Integrate compilation and installation into both AFNI Make and CMake builds.
- Install alongside the existing plotters initially, with a distinct name.
- First-release windows provide pan/zoom, reset view, hover values, and Save.
- First-release file formats are PNG and SVG; PDF follows in a later milestone.
- An output filename selects save-and-exit; no output selects a native window;
  `-show` can request a window as well as file export.
- Follow Python-style defaults with C option aliases and the RSA options.
- Support both pixel dimensions and inches/DPI, rejecting conflicting settings.
- GUI Save defaults to the current view after pan/zoom, with a full-plot option;
  render at the requested export resolution rather than capturing the window.
- CLI replacement of an existing output file requires `-overwrite`; GUI Save
  asks before replacing an existing file.

Visual identity need not duplicate the older plotters. Full compatibility with
every historical option is a staged goal, not a claim for the initial release.

## Current implementation reference points

- `src/1dplot.c`: C options, input processing, transforms, and plotting behavior.
- `src/python_scripts/scripts/1dplot.py`: Python entry point.
- `src/python_scripts/afnipy/lib_apqc_io.py`: Python options and plot models.
- `src/python_scripts/afnipy/lib_plot_1D.py`: Python data preparation and drawing.
- `src/python_scripts/afnipy/lib_apqc_tcsh.py`: `afni_proc.py` QC plotting commands.
- `src/3dRSA.c`: generated example commands near the plotting suggestions block.
- PR #919: histogram, heatmap, and scatter additions, including highlighted
  vertical lines. Head inspected during planning: `aed65af011903d75b9087a83f27403369217d145`.
- `src/Makefile.INCLUDE`: program lists, developer targets, and installation.
- `src/CMakeLists.txt`, `cmake/afni_cmake_build_options.cmake`, and
  `cmake/get_build_macros_and_functions.cmake`: CMake integration and parity checks.
- `src/python_scripts/scripts/build_afni.py`: build orchestration and installation.

Use these as behavioral references. Do not import unrelated development branches
to obtain plotting functionality.

## Prototype results (2026-10-26) — Complete for the macOS prototype

- Implemented native egui/eframe shell with panel focus, pan/zoom, reset, hover,
  current-view/full-plot export, dimension controls, and replacement confirmation.
- Implemented CPU PNG and outlined-vector SVG export, `.1D` inputs/selectors,
  line plots, scatter/identity plots, histograms, and reordered heatmaps.
- Added a four-panel synthetic RSA demonstration, input fixtures, a reproducible
  gallery script, and a local macOS `.app` packaging script.
- Chose a shared SVG scene with usvg/resvg/tiny-skia instead of separate Plotters
  and egui_plot renderers. This gives preview/PNG/SVG a common layout and geometry;
  automated raster comparison verifies outlined SVG against the preview scene.
- Standalone Make and CMake wrappers build/install successfully into temporary
  prefixes. AFNI's top-level build and `build_afni.py` remain unchanged; M2 still
  includes their integration and packaging/program-list work.
- Rust 1.99.1 / macOS ARM64 tested; manifest baseline 1.95 comes from dependencies
  and has not been tested separately. Cargo.lock is included. No global toolchain
  configuration was changed.
- GUI and headless test suites cover numerical behavior, real-file CLI exports,
  dimensions, replacement protection, and in-process egui interaction. Clippy and
  formatting checks are included in the documented verification commands.
- Generated PNG figures and the revised macOS window were visually inspected.
  Native menu and Save dialog checks ran through desktop automation; full-screen
  entry/Escape restoration and clipboard images in Preview were also verified.
  Linux Wayland and other architectures remain unverified.
- UI refinement: no brand banner; Preview-style scrollable thumbnail sidebar,
  graph selection and drag reordering; Shift-drag rectangular zoom, ordinary
  drag pan, double-click reset; right-click single/all graph image copy/save.
- Native Save dialogs (rfd), macOS AppKit menus (muda), in-window Linux menus,
  optional image clipboard support, and a vector P/plot icon with macOS packaging.
- Explicit full-screen escape via a visible button, Escape, and platform menu
  shortcut. Native window decorations remain enabled.
- `-panel-title N TEXT`, `-panel-order 3,1,2`, and `-columns N` supplement the
  existing overall `-title`; one-based panel identifiers resolve before ordering.
  GUI title edits and panel order travel into saved/copied images.
- Fully AppKit-styled interior controls are not required; retain egui and avoid
  adding another full GUI toolkit. New desktop dependencies remain optional.
- X11 window features are disabled. Linux's clipboard dependency still includes
  pure-Rust X11 protocol code; the headless build excludes GUI dependencies. The
  tested headless macOS binary links only libSystem and libiconv.

See `1dPetePlot/README.md` for commands, scope limits, and validation details.
Remaining M0 work is Linux Wayland/headless validation and portability refinement;
the full compatibility and AFNI integration milestones remain staged below.

## Architecture selected for the prototype — Complete (design selection)

Use one data/plot model with separate export and native-window adapters:

```text
AFNI-style CLI -> validated options -> input reader -> numerical preparation
                                                        |
                                                     PlotSpec
                                                        |
                                  shared axes, colors, labels, layout policy
                                         /                         \
                              CPU file exporter              native window
                              PNG / SVG / ...                egui / eframe
```

1. **Input and computation:** parse `.1D` tables and selectors; validate data;
   calculate transforms, histogram bins, boxplot statistics, and matrix ordering.
   Keep these operations independent of graphical libraries.
2. **Plot specification:** represent panels, series, matrix cells, ranges,
   annotations, legends, font choices, dimensions, and color normalization.
   Both frontends consume the same computed values.
3. **Export:** one SVG scene is resolved through usvg, with bundled font outlines,
   and rasterized through resvg/tiny-skia for PNG. SVG export retains the vector
   scene. PDF follows through a separately validated backend/conversion library.
   No runtime shelling out to Python or a screenshot utility.
4. **Window:** egui/eframe with wgpu presents the shared raster scene. A custom
   interaction adapter updates data bounds for pan/zoom, reset, and hover; Save
   re-renders the scene at export resolution. All statistics and scale decisions
   remain in the shared model. Profile large data before deciding whether a more
   specialized interactive drawing path is necessary.
5. **Interaction state:** keep zoom/pan bounds distinct from the original plot
   specification. Save defaults to the current view and offers the full plot as
   an alternative. Capture current axis bounds and re-render at export resolution;
   do not use the displayed image's pixel dimensions as the export size.

The prototype establishes the rendering adapter. Continue assessing text layout,
large-data performance, heatmap labels/colorbars, and native export fidelity
before expanding full compatibility. Avoid independent scientific pipelines.

### Native platform policy

Use macOS native windows and Linux Wayland for the initial desktop targets.
An X11-only Linux session can use file export, but cannot open the proposed native
window without a Wayland environment. Do not silently fall back to X11.

egui is the UI layer; the window backend and Cargo features determine whether
X11 is used. Explicitly select dependencies/features and audit the resolved tree:
eframe defaults currently include both Wayland and X11. Platform selection alone
does not guarantee the whole dependency tree is free of X11 requirements.

Provide a build with GUI dependencies disabled for compute nodes and minimal
installations. Also verify that export in the desktop build never initializes the
GUI/GPU. Whether that desktop binary can load on a machine lacking GUI libraries
must be tested; a separate headless build avoids that dependency assumption.

Windows support can remain a future target unless requested. Native GUI renderer
choice (for example wgpu versus glow) follows the portability prototype, not a
promise that every installed graphics driver will work.

## Command and output contract

Invocation, first-release formats, compatibility defaults, size units, GUI Save
view, and overwrite policy are accepted (D1-D3, D5, D6a-D6b, D6d). Final sizing
flags and default dimensions are implemented for the prototype (D6c):

- An output option such as `-prefix figure.png` saves and exits by default.
- With no output option, open a native window in a GUI-enabled build.
- `-show` explicitly requests a window, including after exporting a file.
- Without an output option or an available native display, report how to request
  an output file. Do not invent an output filename or silently launch another tool.
- **Complete:** `-width W -height H` specifies pixels (default 1500 × 1200).
  `-figsize W H` specifies inches with `-dpi D` (default 150); conflicting size
  specifications are rejected. Rounding and supported sizes are documented.
- Export dimensions describe the entire figure, including labels and colorbars;
  automatic layout must not silently change the requested pixel dimensions.
- Window dimensions and desktop scaling are independent of export resolution.
- GUI Save uses the current view by default, with a full-plot option. Preserve
  axis bounds, labels, annotations, and color normalization when exporting the
  current view. A new export size changes layout, not the scientific data range.
- If a CLI output already exists, fail without changing it unless `-overwrite`
  is present. GUI Save asks before replacement, including when CLI `-overwrite`
  was used for the initial export. Cancellation leaves the existing file intact.
- Plot choices, ranges, titles, labels, colors, annotations, and sizes are fully
  controllable from the CLI. GUI defaults must not be required to complete a plot.
- Support `-help`, useful diagnostics, and nonzero exit status on failure.

Examples of the implemented prototype interface (supply the named input files):

```sh
1dPetePlot -histogram -histbin 30 -vline 0 \
  -vline_hl 0.42 'mean 0.42' -infiles effects.1D \
  -width 1600 -height 1000 -prefix effects.png

1dPetePlot -heat -zerocenter -cbar Reds_and_Blues_Inv \
  -reorder order.1D -infiles rdm.1D -prefix rdm.svg

1dPetePlot -scatter -diagonal -xfile model_a.1D \
  -infiles model_b.1D -show
```

Build a compatibility table before implementing the full parser. Record accepted
aliases, defaults, conflicts, and unsupported features. In particular:

- `-one` / `-one_graph`, `-x` / `-xfile`, positional input / `-infiles` / `-yfiles`.
- Preserve the distinction between `-hist` (draw supplied heights) and
  `-histogram` (bin observations).
- Use Python-style `-yaxis` syntax as the baseline; handle any supported C syntax
  explicitly rather than guessing when arguments are ambiguous.
- Use Python panel order (first input at the top); support `-reverse_order`.
- Parse negative numbers and variadic arguments without confusing them for flags.
- Reject unsupported options rather than accepting them without an effect.
- Scope legacy `-D` environment settings, expression syntax, and label escapes
  individually; familiar spellings do not imply full legacy emulation.

## Feature milestones

**Status key:** **Complete** identifies finished work within the stated scope;
**Remaining** identifies unfinished work; **Deferred** identifies later work that
is intentionally outside the current prototype. A partially complete milestone
lists its completed parts separately rather than marking the whole milestone done.

### M0: decisions and portability prototype — Partially complete

- **Complete:** Initial user decisions D1-D3, D5, D6a-D6b, and D6d are settled;
  renderer and prototype sizing choices D4 and D6c are selected.
- **Complete:** Cargo project under `src/pmolfese/1dPetePlot/`, package
  `afni-peteplot`, and executable `1dPetePlot`.
- **Complete:** Native macOS line plots and labeled heatmaps, plus CPU export.
  Native window, menu, Save dialog, full-screen escape, and clipboard images in
  Preview inspected on macOS ARM64.
- **Complete:** Exact-size PNG, vector SVG, bundled fonts, clipping, and labels;
  exported SVG/CPU-preview comparison and generated-figure inspection.
- **Complete:** GUI-independent export and a GUI-disabled build tested on macOS;
  CLI export tests run without DISPLAY or WAYLAND_DISPLAY. X11 window features
  are disabled and the resolved window-backend feature graph was inspected.
- **Remaining:** Exercise native Wayland and headless execution on Linux;
  macOS results do not establish Linux runtime support.
- **Complete:** Dependency versions resolved in Cargo.lock; core rendering/UI
  dependencies pinned; Rust minimum version declared as 1.95.
- **Remaining:** Validate the declared Rust 1.95 minimum separately from the
  tested Rust 1.99.1 toolchain; resolve remaining deployment/compatibility details
  recorded under D7-D8.
- **Complete:** Shared SVG/usvg/resvg renderer for PNG/SVG and egui/eframe viewing.
  PDF remains a later milestone.

### M1: RSA and basic plotting — Complete (prototype scope)

- **Complete:** Numeric whitespace-delimited `.1D` files, comments, multiple
  columns/files, and stdin. Row/column selectors support ranges, end-of-axis
  notation, ordering, and duplicates, with invalid-input tests.
- **Complete:** Current supported input syntax and deferred conveniences are
  documented in the README. Transpose, inline data, and full metadata handling
  remain explicitly deferred rather than being implied by `.1D` support.
- **Complete:** Line/scatter plots, overlays and separate panels, shared ranges
  or independent y ranges, explicit/generated x values, legends, labels, colors,
  and dimensions. Remaining legacy options belong to M3.
- **Complete:** True histograms with count/density, bin count/width, reference
  lines, and highlighted labeled lines. Supplied-height `-hist` is explicitly
  distinguished and deferred; histogram overlay shared-bin behavior is also
  documented as future work.
- **Complete:** Scatter identity lines with equal x/y limits and physical scaling.
- **Complete:** Square heatmaps with validated row/column permutation, colorbar,
  symmetric zero-centered limits, and documented matrix orientation.
- **Complete:** Selected viridis and red/blue palettes with the PR's
  Reds_and_Blues aliases. Sampled palettes are documented; full Matplotlib
  colormap parity is not claimed.
- **Complete:** CLI PNG/SVG output and native viewing with pan/zoom, reset, hover,
  and current-view/full-plot Save. Invocation, sizing, and overwrite contracts
  are implemented. Linux runtime validation remains tracked in M0/M2.
- **Complete:** RSA-style histogram, scatter, and reordered-heatmap commands
  exercised against fixtures through CLI tests and `scripts/render_gallery.sh`.
  Exhaustive behavioral comparison with the legacy plotters remains later work.

### M2: build and distribution integration — Partially complete

Begin build plumbing during M0/M1; finish release integration here.

- **Complete:** Standalone Make wrapper runs Cargo, supports GUI/headless builds,
  installs through PREFIX/DESTDIR, and provides test, gallery, app, and clean
  targets. Build/install into a temporary prefix was verified.
- **Remaining — AFNI Make integration:** Connect existing developer-directory
  hooks; install into `$(INSTALLDIR)`; register program discovery lists and source
  distribution; validate clean/dependency rules during parallel build/install.
- **Complete:** Standalone CMake wrapper locates Cargo, defines a custom target
  and binary byproduct, puts Cargo output under the build tree, delegates
  incremental rebuilding to Cargo, and installs into a temporary prefix.
- **Remaining — AFNI CMake integration:** Install through
  `${AFNI_INSTALL_RUNTIME_DIR}` and register AFNI components and build-parity
  checks; wire the prototype into AFNI's top-level configuration.
- **Complete:** Optional Cargo GUI feature, Make `GUI=0`, CMake `PETEPLOT_GUI=OFF`,
  and a headless executable without GUI dependencies.
- **Remaining:** Separate top-level AFNI inclusion/GUI switches, initial opt-in
  policy, and prerequisite diagnostics for explicitly requested Rust tools.
- **Remaining:** `build_afni.py` prerequisite checks and option propagation through
  Make/CMake targets. If its CMake installation path is in scope, handle that
  explicitly; it currently builds without running installation.
- **Complete:** Cargo.lock generated and included in the working source; wrappers
  use `cargo build --release --locked`. Cached offline builds and the tested
  toolchain are documented.
- **Remaining:** Commit the source/lockfile when requested, select and pin the
  supported release toolchain, and define a vendor/cache distribution workflow.
  `--locked` by itself does not prohibit network access.
- **Complete:** macOS ARM64 GUI/headless builds and local `.app` packaging,
  including the custom icon.
- **Remaining:** Select deployment baselines with AFNI distribution requirements;
  build/test macOS x86_64 and Linux x86_64/ARM64. Signed/notarized macOS distribution
  and official release packaging remain unfinished.
- **Complete:** Bundled Ubuntu font license included; macOS runtime libraries and
  Cargo window-backend features inspected. The headless macOS executable links
  only libSystem and libiconv; no Python, Matplotlib, or XQuartz runtime is needed.
- **Remaining:** Complete dependency/license and executable-loading checks for
  each distribution target, including machines without desktop libraries.

### M3: broader existing plotting coverage — Partially complete

#### M3a: uncertainty and Python-style QC — Complete

- **Complete:** User-supplied `-yfiles_pm` bands and `-ylim_use_pm`, including
  selectors, matching dimensions, absolute half-widths, and scaling with curves.
- **Complete:** Union of `-censor_files` / `-censor FILE` and `-censor_trs` index
  lists; translucent highlighting retains all samples; hover identifies censoring.
- **Complete:** Per-series `-censor_hline` (including zero/negative/NONE), configurable
  censor RGB(A), and `-patches` run backgrounds. Increasing explicit x coordinates
  use midpoint sample boundaries. Mask/run lengths must match selected rows.
- **Complete:** Before/after-censoring boxplots, BC_ONLY/AC_ONLY selection, linear
  percentile quartiles, 1.5-IQR observed whiskers/outliers, counts and empty AC states.
- **Complete:** Python `-scale` numbers, SCALE_TO_HLINE, and SCALE_TO_MAX; one factor
  or one per series. Curves/bands/thresholds scale together; axis limits do not.
- **Complete:** Per-panel/partial y ranges; VOLREG and ALLINPAR label presets,
  label wrapping; `-margin_off` including the APQC 12 × 0.5 inch strip at 150 DPI.
  Automatic default height accommodates motion panels; explicit sizes stay exact.
- **Complete:** Synthetic acceptance commands drawn from `lib_apqc_tcsh.py`, plus
  existing RSA histogram/scatter/heatmap fixtures. PNG/SVG gallery inspected and
  numerical, CLI, and GUI-state regression tests added. No new Cargo dependencies.

- **Complete (2026-10-27):** Independent repeatable `-shade START:END` annotations
  with per-range color, label, and opacity; labels wrap in the key, ranges clip to
  the viewport, and censor masks/statistics remain unchanged. Line/scatter/histogram
  support; annotations survive native selection/reordering and current-view export.
- **Complete (2026-10-27):** `-shade-xmat` reads AFNI condition groups/labels and
  combines basis columns into modeled-response activity intervals, with a palette
  override and adjustable relative-magnitude threshold (default 5%). Standard
  matrices contain no color metadata. GoodList maps censored rows onto NRowFull;
  absent rows remain unshaded. Full-row matrices are preferred. Plot x coordinates
  determine units; RowTR does not silently change them. This is model activity,
  not reconstructed exact stimulus timing. A real 3dDeconvolve-generated synthetic
  matrix fixture, alignment tests, and PNG/SVG acceptance examples verify behavior.

#### M3b: extended inputs and everyday parity — Partially complete

- **Complete:** Row/column selectors, explicit/generated x, per-series colors,
  vertical reference lines, and highlighted labeled lines from the prototype.
- **Complete:** Inline `1D:` repeats/multiple series and trailing-apostrophe transpose.
  Named-column TSV input moved forward from M4 at the user's request.
- **Complete:** Explicit missing-data policies: strict default, gaps, or selected
  column-mean imputation with diagnostics. No time-index compression; finite-only
  boxplot/histogram statistics. Missing heatmap cells use gray.
- **Complete:** Unicode Greek and explicit brace-delimited sub/superscripts;
  literal bare underscores. No full TeX engine required.
- **Complete:** Reviewed installed C/Python help and command generators; detailed
  coverage/semantic differences in `1dPetePlot/M3_COMPATIBILITY.md`.
- **Remaining:** Normalization, demeaning/detrending, common log transforms,
  `-xvals`, multiple x columns, broader line/marker styling, font controls,
  legend locations, exact tick/range controls, and further compatibility aliases.

#### M3c: extended C compatibility — Remaining / deferred

- **Remaining:** Supplied-height `-hist`, drawing rectangles/line segments,
  run-relative `-CENSORTR`/`-concat`, multiple censor-color layers, and remaining
  legacy input/label options selected by actual use.
- **Deferred:** General expression language and full TeX-like formatting. User
  accepted simpler Python-oriented transforms and basic scientific labels first.
- **Deferred:** JPEG/PNM/PostScript; unsupported output extensions fail explicitly.
  APQC currently writes JPEG and HTML references; changing generators/consumers
  is a separate integration step, not implied by PNG-equivalent acceptance tests.

### M4: RSA improvements and desktop refinement — Partially complete

- **Complete:** Heatmap index labels track reordered conditions; explicit color
  limits and zero-centered limits are available.
- **Remaining:** User-supplied condition labels, automatic shared color limits
  across matrices, optional matrix masks, and improved publication layouts.
- **Complete:** Named-column TSV selection (delivered with M3).
- **Complete in M6 interactive imports:** Table headers become series labels.
  **Complete in M7:** Semantic ordinary ROI-wide RSA tables.
  **Remaining:** General CLI header defaults and advanced scientific table formats.
- **Complete:** Preview-style scrollable thumbnail sidebar, individual/all-graph
  views, drag reordering, and Earlier/Later controls. Graphs retain their view
  bounds when rearranged.
- **Complete:** Figure/per-panel title editing and column layout in the GUI;
  `-title`, repeatable `-panel-title N TEXT`, `-panel-order`, and `-columns` in the
  CLI. Saved/copied images use the edited titles and arrangement.
- **Complete:** Shift-drag box zoom, ordinary drag pan, double-click reset, and
  full-screen escape by visible button, Escape, and platform menu shortcut.
- **Complete:** Right-click copy/save for one graph or all graphs at export
  resolution; native Save dialogs; native macOS menus; custom P/plot icon; banner
  removed. macOS clipboard images verified in Preview.
- **Complete (implementation):** In-window Linux menus and desktop-portal Save
  dialogs, with Wayland clipboard support enabled. **Remaining:** Linux runtime
  checks, including compositor clipboard availability and dialog backends.
- **Complete in the M5 first slice:** Reusable figure recipes and direct text editing.
- **Complete in M6 first slice:** Mixed plot types from files and interactive import.
- **Remaining — expanded in M5–M7:** Richer document/styles, scientific import
  metadata/attachments, and task presets.
- **Deferred:** Direct PDF export through a validated backend/conversion library.
- **Deferred:** Update 3dRSA's suggested executable only after distribution support
  exists; installing a new tool does not update the commands it prints.

## Product direction: an editable scientific figure workspace — Proposed 2026-10-27

The user requested organization and planning for direct title/caption editing,
drag-and-drop data with suitable plot previews, and task-family entry points such
as `-rsa`. M5's first document/editing slice is now implemented; import and
presets remain roadmap additions. Preserve the
existing native-window, no-X11, shared-renderer, headless export, and minimal
additional dependency requirements. Keep M2 distribution and remaining M3 parity
work visible; the new desktop roadmap does not supersede them.

The organizing principle is **one editable figure document shared by the GUI,
CLI, and exporter**. A useful workflow is:

1. Open/drop data, or launch a task preset.
2. Inspect the recognized columns, units, missingness, and data roles.
3. Choose a suggested plot preview or assign x/y/group roles explicitly.
4. Add the plot to the figure; edit labels, captions, colors, and arrangement.
5. Save the editable project, export an image, or copy a reproducible command.

Keep the current left plot-thumbnail navigator. A source-data drawer or separate
sidebar tab can hold imported datasets without confusing datasets with panels.
Preview cards are temporary candidates until the user chooses Add to figure;
opening another dataset should not replace an existing figure silently.

### Feasibility and boundaries

| Idea | Assessment | Main work / boundary |
|---|---|---|
| Click a title, axis label, or caption and edit it | High feasibility; bounded work | Shared renderer reports semantic text hit regions; native overlay text field edits model strings. Captions need new fields/layout; existing inspector title editing is already implemented. |
| Drag local files into the window | High feasibility; moderate workflow work | Existing egui version exposes drop events; add import, cancellation, previews, errors, and Open File fallback. Verify macOS and supported Wayland desktops. |
| Suggest plots appropriate to the data | High feasibility for known schemas; moderate-to-large scope | Deterministic shape/metadata rules and explicit column-role selection. Numeric shape alone cannot identify scientific meaning. |
| `-rsa` / later `-qc` presets | High feasibility once source schemas are represented | Reuse existing plot types; add readers, grouping, joins, mixed-panel construction, and recorded preset defaults. |
| Save a figure and reopen/edit it | High feasibility; foundational moderate work | Versioned document, stable IDs, source references, data mappings, undo history boundaries, and reproducible export. |
| Linked ROI/condition selection across plots | Feasible, later | Shared row/condition identities and joins; equal row counts do not prove correspondence. |
| Instant previews for arbitrarily large inputs | Cannot promise universal instant response | Background work, progress/cancel, cached previews, bounded memory, and declared downsampling. |
| Exact event timing from a convolved X matrix | Not generally recoverable from that matrix alone | Show modeled activity as implemented; exact event shading needs original timings or explicit ranges. |
| Recover a fully editable project from PNG or outlined SVG | Original semantics cannot generally be recovered | Save the project/recipe alongside the image; pixels/outlines do not retain data-role mappings. |
| Direct BRIK/NIfTI/GIFTI/CIFTI plotting | Possible but a separate substantial scope | Format readers, spatial/temporal axes, masks/ROI extraction, and user-selected reductions. Initially import tables and AFNI command outputs. |
| Full vector illustration editor / arbitrary layout canvas | Possible, expensive to maintain | Prioritize scientific panels, captions, guides, and templates before freeform vector authoring. |

Technical evidence: egui 0.36.2 exposes `hovered_files` and `dropped_files` through
[RawInput](https://docs.rs/egui/0.36.2/egui/struct.RawInput.html). This confirms the
input hook, not cross-desktop runtime behavior. Figure/Plot remains the shared
computed model; M5 now stores source arguments and editable presentation separately,
with fingerprints and headless recipe replay. M6 adds independent source layers,
original column indices, and explicit numeric roles; richer scientific identity
and semantic roles for ordinary ROI-wide RSA are now implemented in M7;
long schemas and cross-file scientific joins remain future work. Current
`3dRSA.c` writes mixed ROI/label/nvox/effect/statistic columns and separate long
formats for series/run-resolved output; its `.rsa.1D` requires a semantic reader,
not just the existing numeric `.1D` parser or TSV reader.

### M5: editable documents and direct figure editing — Partially complete

**First usable slice — Complete (implementation and automated verification), 2026-10-27.**
Native project-dialog and direct-editing inspection remains pending: desktop
automation timed out while attaching to the launched macOS app. Earlier native
window/menu/image-dialog checks do not count as verification of these new controls.

**M5a — document foundation: Partially complete**

- **Complete:** Version 1 `.peteplot.json` recipe stores source arguments (including
  selectors, transformations, plot options, and annotations), stable panel IDs,
  editable titles/axis labels/captions, display order, column count, per-panel views,
  export size, and Current view/Full plot preference. Computed data and raster
  caches are separate; reopening reconstructs data and reapplies presentation.
- **Complete:** Undo/Redo for committed text edits, panel reordering, and column
  layout; 100-entry session history, modified indicator, and Save/Discard/Cancel
  on open/close. Selection/view bounds stay with panel identity through reordering
  and undo. Pan/zoom and export preferences are not presentation-history actions.
- **Complete:** Save Project, Save Project As, and Open Project with native file
  dialogs; atomic writes and replacement confirmation. Save Image remains separate.
  Projects always retain all panels, regardless of sidebar selection.
- **Complete:** Relative source paths and deterministic byte fingerprints for data
  and auxiliary files. Missing/changed files fail explicitly; source changes are
  checked again before saving. Moving the project and sources together while
  preserving their relative layout works. Symlink basenames retain legend labels.
- **Complete:** Headless `-recipe FILE` replay and Copy replay command for saved,
  unmodified projects. Export/size overrides are explicit; drawing overrides are
  rejected. Replay exports all panels using saved view/full-range preferences.
- **Complete in M6:** Version 2 mixed-source recipes, stable import-layer IDs and
  original column indices, explicit numeric role mappings, and version 1 read
  compatibility. Undo/Redo includes Add Graph; saved projects exclude undone imports.
- **Remaining:** GUI relink/reload for existing panels, richer scientific roles,
  style actions, and subsequent schema migrations. Import Retry rereads a pending
  source; it does not replace an already-added graph.
- **Deferred:** Portable bundles/data snapshots and persistent undo history. Current
  recipes reference external inputs; demos and inline `1D:` data are self-contained.
  Stdin needs to be saved to a file and reopened before Save Project is available.

**M5b — direct editing: Partially complete**

- **Complete:** Click rendered figure/panel titles, axis labels, and existing
  captions using shared semantic text regions. An egui editor opens by the clicked
  region; Enter commits simple fields, Done commits captions, Escape cancels.
  Inspector buttons support keyboard access and adding initially empty captions.
  Active editing takes priority over navigation/full-screen Escape.
- **Complete:** Separate figure/panel caption fields; generated notes and scientific
  provenance remain distinct. `-caption` and repeatable `-panel-caption N TEXT`
  also work from the CLI, before panel reordering.
- **Complete:** Titles/captions wrap using the bundled font metrics; long axis labels
  scale to fit. PNG/SVG and window share the layout. Oversized edits remain
  uncommitted with an actionable message, rather than silently losing text.
  Margin-free QC strips intentionally omit text and reject caption options.
- **Remaining:** Figure-wide styles, selected-panel style overrides, consistent
  condition color controls, duplicate panel, and broader accessibility/typography
  refinement. Compact navigator previews omit captions; full views/exports retain them.

**Evidence:** Automated edit/type/commit/cancel/undo/redo/reorder/save/reopen tests;
1×/2× semantic click targets and focused-panel targeting; source relocation and
changed/missing-input detection; retained NaN gaps; identical CPU-rendered pixels
before and after project reconstruction; CLI PNG/SVG replay, sizing and overwrite
checks; caption wrapping/axis separation. Captioned output visually inspected.
Example project: `1dPetePlot/examples/edited-figure.peteplot.json`.

### M6: interactive data import and plot suggestions — Partially complete

**First usable slice — Complete (implementation and automated verification), 2026-10-27.**

- **Complete:** No-argument native/Finder launch opens an empty workspace. Open
  Data (Cmd/Ctrl-I) and file-drop both enter an explicit preview workflow; the
  synthetic gallery remains an intentional choice or `-demo`. Headless launch
  still requires an input/project and output destination.
- **Complete:** Numeric `.1D` and headered TSV inspection, column names, numeric/text
  classification, missing counts, and the first eight source rows. Text remains
  inspectable and is not coerced. X/Y and matrix column roles are explicit;
  headers label series, and missing policy is reject/gap/mean, with finite X.
- **Complete:** Deterministic line, distribution, scatter, and square-heatmap
  choices with explanations and a full-data preview of the selected choice.
  A square table is labeled a heatmap, not assumed to be an RDM. No AI/network
  dependency and no preview downsampling.
- **Complete:** Add explicitly commits one graph. Multiple dropped/opened files
  queue separate previews; Cancel/Escape, Skip, and Reload/Retry leave existing
  graphs alone. Import review fills the main window and resizes with it, replacing
  the floating sub-window following user feedback. Cancel/Escape returns to the
  figure workspace. The action footer stays visible while preview/roles scroll,
  following user feedback that Add was too easy to miss below the fold.
- **Complete:** One background worker per window with bounded jobs, cancellation,
  stale-result suppression, status/error display, and fingerprint checks before/
  after reading and before Add. Source changes require explicit retry.
- **Complete:** Stable source-layer/column identities, mixed-source/plot-type
  figures, Undo/Redo Add, version 2 Save/Open Project, relative source paths,
  and headless PNG/SVG replay. Version 1 projects remain readable. Undo back to
  empty is supported; saving excludes undone imports. Export height grows as
  needed on Add; export preferences remain separate from presentation Undo.
- **Complete:** Bounded interactive inputs: 16 MiB, 100,000 rows / 1,000,000 cells,
  256 columns; 16 series and 50,000 selected line/scatter points; 256-square heatmaps;
  16 queued files and 64 retained graphs including Undo. Limits report actionable
  errors and do not silently sample or truncate scientific observations.
- **Remaining:** Units, row/ROI identities, run/TR metadata, uncertainty/label/group
  roles, cross-source joins/alignment, and recognized X-matrix design views or
  censor/condition attachment suggestions. Existing CLI X-matrix shading is unchanged.
  Require matching timelines/metadata and explicit attachment targets.
- **Remaining:** Simultaneous candidate thumbnail gallery (the first slice shows
  one selected preview), CSV, advanced RSA schemas, richer metadata, and
  large-data optimization beyond the explicit interactive limits.
- **Remaining:** Native revised-footer/Open Data dialog and Wayland drop checks.
  The user confirmed native drag/drop reaches the preview. Desktop automation
  could launch the rebuilt macOS app but timed out attaching; in-process tests
  do not substitute for those platform checks.

**Evidence:** Queue/skip/retry/cancel/invalid-file/changed-source tests; an actual
pointer click on Add without scrolling at 850 × 500; full-window review bounds at
850 × 500 and 1440 × 900, and Escape return; preview-to-added-panel pixel
comparison; mixed-source Undo/reorder/save/reopen pixel equality; headless version
2 PNG/SVG replay and malformed identity/role rejection; finite-X/missing-policy
validation. A 25,001-point inspection and 900 × 600 PNG preview took approximately
0.36 seconds in the macOS ARM64 debug test, with all points retained. Oversized
selection rejects explicitly; two histograms count all 25,001 observations each.
This is measured bounded work, not a universal UI responsiveness guarantee.
The three-panel mixed-source example export was visually inspected:
`1dPetePlot/examples/imported-figure.peteplot.json`; try importing its synthetic
`import-demo.tsv` source through the window.

### M7: scientific presets, beginning with RSA — Partially complete

**First usable slice — Complete (implementation and automated verification), 2026-10-27.**

- **Complete:** Semantic reader for ordinary ROI-wide `3dRSA` tables, using the
  current local `3dRSA.c` header/output code as the reference. Supports commented
  headers, labels with spaces, named models, unique ROI identities, and voxel counts.
  Numeric RDM input remains the ordinary `-heat` path.
- **Complete:** Distinct r, beta b, and partial-r roles; default r for ordinary
  comparison and b for regression. Probability/FWE, interval, fit/weight, contrast,
  LOO, reliability and noise-ceiling columns are retained as auxiliary values,
  not plotted as primary effects. Probabilities/interval ordering are validated.
  Unknown/commonality and long series/run schemas fail explicitly.
- **Complete:** `-rsa FILE` generates per-model distributions, one ROI-matched
  model comparison with ROI/label/nvox hover, and a rectangular ROI-by-model effect
  heatmap. Separate row/model labels and symmetric color normalization; a coefficient
  heatmap is explicitly distinguished from an RDM. Selected missing values fail
  by default; `-missing gap` omits incomplete pairs together and marks missing cells.
  Imputation and automatic significance filtering are absent.
- **Complete:** `-rsa-models`, `-rsa-effect`, `-rsa-rois`, `-rsa-plots`, and
  `-rsa-pair` provide deterministic selection. Default all models up to 8; larger
  sets require selection. Default all ROI IDs in ID order, with requested order
  preserved. Heatmap limited to 40 selected ROIs: omitted with a caption explanation
  in larger default summaries, rejected if explicitly requested without subsetting.
  All selected finite observations still contribute to distributions/scatter.
- **Complete:** Open Data/drop recognizes the semantic RSA header and offers one
  distribution, comparison, or heatmap per Add in the full-window review. Model
  effect columns and ROI IDs are selectable; statistics remain non-plottable metadata.
  Existing background work, cancellation, source verification, Undo, and mixed-source
  project handling apply. Multiple distributions can be added separately.
- **Complete:** Native-window or headless PNG/SVG output through the shared model.
  Save Project freezes resolved models/effect/ROIs/modes/pair and preset version 1
  in the source recipe. Fingerprinted relative paths, edits, identities and exact
  replay pixels are retained. No new dependencies.
- **Complete:** Explicit source/plot limits and automatic heatmap label sizing;
  oversized explicit dimensions/layouts fail before export. Input limits are 16 MiB,
  10,000 ROIs / 1,000,000 cells / 512 columns (256 columns for interactive import).
- **Remaining:** Supplied asymmetric confidence-interval ROI plots; actual RDM
  inputs with condition-label matching and explicit shared scales; long series/run,
  fits/contrasts/commonality visualization; multi-file ROI joins. A results table
  cannot reconstruct an RDM. Recognizing auxiliary columns does not complete their
  plotting workflows or scientific method labeling.
- **Remaining:** Native hands-on inspection of the RSA import controls and hover;
  Linux/other architecture checks. In-process GUI tests and CPU export inspection
  are evidence for implementation, not those platform checks.

**Implemented invocation examples:**

```sh
1dPetePlot -rsa results.rsa.1D
1dPetePlot -rsa results.rsa.1D -prefix rsa-summary.png
1dPetePlot -rsa results.rsa.1D -rsa-models semantic,visual \
  -rsa-rois 3,12,27 -rsa-plots scatter,heatmap -prefix selected.svg
1dPetePlot -rsa regression.rsa.1D -rsa-effect pr -prefix partial-r.png
1dPetePlot -recipe figure.peteplot.json -prefix publication.svg
```

**Evidence:** Synthetic ordinary/regression fixtures with optional fitted metrics,
probabilities, bounds and reliability; unordered ROI rows; missing effects/complete
pair omission; unknown/long/duplicate schemas; missing requested ROI/model checks;
rectangular heatmap hover/orientation; semantic Open/Drop review and Add/Undo/reopen;
headless PNG/SVG exports and output preservation on validation failure; frozen
preset/relative-source relocation with identical replay pixels and changed-source
refusal. The ordinary four-panel summary PNG was visually inspected.
Fixtures: `1dPetePlot/examples/rsa/{ordinary,regression}.rsa.1D`.

**Later QC preset:** load a compatible motion/outlier/censor/design set, suggest
existing QC plots, and validate row/run/TR alignment. Recognizing a folder can
suggest related files; show which files will be used rather than guessing silently.

**Remaining validation for later slices:** RDM condition-label mismatch fixtures,
asymmetric interval rendering with explicit method/scale labels, long-schema
recognition, and cross-file ROI joins. Continue preserving GUI/headless equality
and the no-inference/no-undocumented-filtering contract.

### Recommended order and useful later additions

1. **Complete (first slice):** Minimal document/action foundation, direct title/
   caption editing, Undo/Redo, Save/Open Project, and headless replay. Finish the
   new native-control hands-on check and continue the remaining M5 items.
2. **Complete (first slice):** Open/Drop, numeric column preview, deterministic plot
   choices, Add/Undo, and mixed-source replay. Continue scientific attachments later.
3. **Complete (first slice):** Ordinary-table `-rsa` and semantic import. Continue
   supplied intervals, actual RDM inputs, and richer result schemas next.
4. Add publication templates, explicit linked axes/colorbars, linked ROI/condition
   selection, and replace-data-with-same-schema while preserving figure styling.
5. Consider editable shading handles with explicit coordinate readouts and undo,
   explicit Reload Data, batch application of a saved recipe, and an export gallery.
   Live file watching should be opt-in and clearly mark changed source data.

Useful product improvements: show why each plot is suggested; let users inspect
which source column produced a curve; save condition colors across the figure;
provide accessible palettes and readable multi-panel presets; make changing one
figure's layout a reversible action. Keep presentation editing separate from any
new analysis/inference features, which need their own statistical specifications.

Keep estimates qualitative until prototypes and representative large files exist.
The direct-editing, initial import and ordinary RSA slices are delivered; scientific
attachments and broader presets still need their own semantic contracts. No GUI
framework replacement was needed. Serialization uses serde/serde_json; font metrics
reuse ttf-parser, already present in the renderer's dependency tree.

## Validation and completion criteria — Partially complete

**Complete for the current implementation:** 71 GUI-build tests and 48 headless
build tests pass; Clippy and formatting checks pass. Coverage includes selectors,
numeric validation, histogram behavior, heatmap permutations, coordinate mapping,
PNG/SVG exports, output sizing/replacement protection, title/order validation,
hover/pan/box zoom, thumbnail selection/reordering, and full-screen Escape.
M3 adds mask union/sample boundaries, box statistics, scaled bands, per-panel
limits, inline/TSV/transpose input, missing-data policies, APQC-shaped CLI
exports, gap rendering, and QC annotations through GUI selection/reordering.
macOS UI inspection and standalone Make/CMake temporary-prefix installs are also
complete for earlier features. M5 adds document replay/source-fingerprint,
caption, direct-editing, and undo tests; native inspection of the new M5 controls
remains pending after desktop automation timeouts. See the prototype results and
README for the exact platform evidence. M6 adds background-worker cancellation,
queue/retry, bounded large-data previews, visible/clickable footer, Add/Undo, mixed
project pixel equality and headless replay tests. Native drop reaches preview per
user feedback; revised-footer/chooser/Wayland checks remain pending.

**Remaining:** Linux and other architecture/runtime checks, release integration,
full legacy behavioral comparisons, expanded layout/font testing, and tests for
features not yet implemented. The following criteria remain the checklist for
those future additions; completion above does not mark every case below done.

Test scientific computations separately from rendering:

- Histogram edge inclusion, constant inputs, negative values, partial final bins,
  and density integral; document the chosen binning semantics.
- Square-matrix validation, orientation, reordering, labels, zero-centered colors,
  and all-zero/constant matrices.
- Shared scatter limits and equal aspect, multiple series, and empty/singleton data.
- Selector semantics, comments, malformed rows, mismatched lengths, non-finite
  values, censor masks, transformation order, and boxplot statistics.
- Define per-mode missing-value behavior; never silently turn NaN into zero.

Use representative C/Python commands as behavioral references, with deliberate
differences recorded. Compare numerical values and plot specifications rather
than requiring pixel identity with Matplotlib or the old C renderer.

Render and visually inspect a compact gallery: line panels, overlays, RSA
histogram, labeled scatter, heatmap, long labels, error bands, and censor/boxplots
as implemented. Check pixel dimensions, unclipped labels, legend placement, and
consistent colors. Test Unicode with the bundled font's documented coverage.

Verify CLI export on a machine/session without DISPLAY or WAYLAND_DISPLAY and
without GPU access; verify a headless build lacks GUI dependencies. Test native
window opening, resizing, closing, and saving on the supported desktops. Inspect
the resolved Cargo feature graph and runtime dependencies for unintended X11 use.

Verify refusal to overwrite from the CLI without `-overwrite`, explicit CLI
replacement, GUI replacement confirmation/cancellation, and preservation of an
existing output if rendering fails. Ensure the no-overwrite policy also handles
a destination created while rendering, rather than relying only on an early
existence check.

Verify Make and CMake build/install paths independently, incremental rebuilds,
installation into a temporary prefix, and program-list/component registration.
Only mark a platform verified after an actual build/runtime check.

## Decision register

Accepted entries record the user's choices. Remaining recommendations are
proposals; update them as decisions are made before dependent implementation work.

| ID | Decision | Recommendation | Alternative / consequence | Status |
|---|---|---|---|---|
| D1 | Initial native-window interaction | Pan/zoom, reset view, hover values, and Save | Preview-only scope not selected | Accepted 2026-10-26 |
| D2 | Invocation behavior | Output filename means save-and-exit; no output means window; `-show` requests both | No-output behavior in a headless build reports how to export | Accepted 2026-10-26 |
| D3 | First-release export formats | PNG and SVG; add PDF next | Direct PDF deferred beyond the first release | Accepted 2026-10-26 |
| D4 | Renderer implementation | Shared SVG scene, usvg/resvg CPU rendering, egui/eframe + wgpu with custom interaction adapter | Profile large-data interaction before extending this into full feature coverage | Selected for prototype 2026-10-26 |
| D5 | Default compatibility dialect | Python-style defaults and RSA flags, with explicit C aliases | Separate compatibility modes only if real conflicts warrant them | Accepted 2026-10-26 |
| D6a | Output size units | Support both pixels and inches/DPI; reject conflicting specifications | Keep the existing Python sizing convention available | Accepted 2026-10-26 |
| D6b | GUI Save after pan/zoom | Default to current view, with a full-plot choice; re-render at export resolution | Full plot remains explicitly selectable | Accepted 2026-10-26 |
| D6c | Remaining size details | Prototype defaults: 1500 × 1200; `-width`/`-height`, or `-figsize` with DPI defaulting to 150 | Exact explicit dimensions; incompatible sizing options rejected; default height now grows for additional panels | Selected for prototype 2026-10-26 |
| D6d | Existing output files | CLI requires `-overwrite` to replace an existing file; GUI asks before replacing | Cancellation/refusal preserves the existing file | Accepted 2026-10-26 |
| D7 | Deployment baselines and GUI packaging | Native macOS + Wayland desktop build, plus headless build; initially opt in through AFNI | Select exact OS/compiler baselines and official package composition with AFNI maintainers | Open |
| D8 | Scientific defaults and extended input | Python-oriented semantics; inline/transpose, named TSV, explicit gap/mean-imputation modes | Strict default; no silent imputation. Remaining bin/colormap differences documented | Accepted 2026-10-26; extended input implemented |
| D9 | M3 priority and acceptance examples | Uncertainty/QC first; use 3dRSA and afni_proc.py QC command patterns | PNG/SVG equivalents now; command-generator/HTML integration later | Accepted 2026-10-26 |
| D10 | Transform scope | Python scaling plus common normalization/detrending/log operations before a general expression language | General AFNI expression evaluator deferred | Accepted 2026-10-26; scaling implemented |
| D11 | Scientific text | Unicode Greek and basic explicit superscripts/subscripts, literal ordinary underscores | Full TeX-like engine deferred | Accepted 2026-10-26; basic rendering implemented |
| D12 | Editable workspace direction | Direct editing, drop/preview/add workflow, and RSA family entry point | Preserve shared CLI/export model and existing scientific contracts | User-requested direction 2026-10-27; M5/M6/M7 first slices implemented |
| D13 | Document and reproducibility design | Versioned recipe, stable IDs, undoable edits, explicit source links, GUI/headless replay | Version 1 JSON recipe selected; rich source/series identities, relinking, and portable bundles remain open | First document/editing/replay slice implemented 2026-10-27 |
| D14 | Plot suggestions and presets | Deterministic schema-aware candidates; role selection when ambiguous; ordinary RSA table first | Universal scientific intent cannot be inferred from numeric shape | First ordinary-table RSA slice implemented in M7; broader schemas remain |

The no-X11 requirement is accepted, not a pending choice. Linux systems without
Wayland retain command-line export. Any future X11 backend would require a change
to the agreed scope.

## External technical references

Checked during planning; recheck exact APIs when selecting dependency versions.

- [eframe documentation and feature flags](https://docs.rs/eframe/latest/eframe/)
- [winit native platform backends](https://docs.rs/winit/latest/winit/platform/index.html)
- [egui_plot plotting primitives](https://docs.rs/egui_plot/latest/egui_plot/)
- [Plotters backends and features](https://docs.rs/plotters/latest/plotters/)
- [Cargo build options](https://doc.rust-lang.org/cargo/commands/cargo-build.html)
- [Cargo vendoring](https://doc.rust-lang.org/cargo/commands/cargo-vendor.html)
- [Rust platform support](https://doc.rust-lang.org/rustc/platform-support.html)
- [CMake custom targets](https://cmake.org/cmake/help/latest/command/add_custom_target.html)
- [CMake installation](https://cmake.org/cmake/help/latest/command/install.html)
