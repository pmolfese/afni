# M3 compatibility inventory — updated 2026-10-27

References: installed `1dplot -help` (AFNI_26.3.01) and `1dplot.py -help`, local
`src/1dplot.c`, `lib_apqc_io.py`, `lib_plot_1D.py`, and `lib_apqc_tcsh.py`.
The [master 3dRSA plotting suggestions](https://github.com/afni/afni/blob/master/src/3dRSA.c#L9557)
were also read through GitHub (file blob `1537c279181123c320cdb42e3bc1067456d229f5`).
Tests compare scientific behavior and useful syntax, not pixel identity.

| Feature / option family | Status |
|---|---|
| Lines, separate/overlay, shared/independent y ranges, order | Implemented |
| RSA histogram, scatter/identity, heatmap/reorder, highlighted vertical lines | Implemented; sampled colormaps and independent overlay histogram bins remain documented differences |
| `-yfiles_pm`, `-ylim_use_pm` | Implemented; match flattened selected columns/rows |
| `-censor_files`, `-censor_trs`, `-censor` | Implemented; union mask, exact selected length, preserved samples |
| Repeatable `-shade`, color/label/alpha | Implemented independently of censoring; line/scatter/histogram; x-coordinate ranges |
| `-shade-xmat`, condition palettes and threshold | Implemented model-activity shading; positive groups, GoodList alignment, full timeline required |
| `-censor_hline`, `-censor_RGB`, `-patches` | Implemented for line plots |
| `-boxplot_on`, `-bplot_view` | Implemented; separate line panels; BC/AC, finite sample counts |
| `-scale` numbers / SCALE_TO_HLINE / SCALE_TO_MAX | Implemented, including band/threshold scaling |
| Python per-panel/partial `-yaxis` | Implemented; one range broadcasts, or one per separate panel |
| VOLREG, ALLINPAR6/9/12, `-volreg`, `-ylabels_maxlen` | Implemented |
| `-margin_off` and 12 × 0.5 inch QC strip | Implemented; no labels, title, legend, or boxplot; grid remains |
| Inline `1D:`, repeats, multiple series, transpose | Implemented; selectors then transpose |
| TSV header-name selection; gap/mean policies | Implemented; exact header names; no implicit imputation |
| Basic scientific label text | Unicode + `_{...}` / `^{...}`; ordinary underscores literal |
| Numeric grayscale, hex colors, common names | Implemented; not every Matplotlib named color |
| `-xvals`, `-xmulti`, `-xl10`, `-norm*`, repeated `-demean` | Remaining |
| `-box`, `-noline`, `-dashed`, `-thick`, `-THICK` | Remaining; `-scatter` already available |
| `-ynames`, legend locations, font options, background color | Remaining |
| `-ignore`, `-use`, C tick syntax `b:t:n:m`, reversed axes, `-nopush`, `-aspect`, `-naked` | Remaining |
| C `-hist`, `-rbox`, `-Rbox`, `-line` | Remaining; supplied-height bars differ from histogram binning |
| Run-relative `-CENSORTR`, `-concat`, layered censor colors | Remaining |
| `-ytran`, `-xtran`, full TeX-like syntax | Deferred general expression/formatting engines |
| JPEG/TIFF/PNM/PostScript, AFNI environment-variable controls | Unsupported; explicit errors for unsupported flags/formats |
| X11-specific options / rendering | Outside scope; native macOS / Wayland remain the window targets |

## Scientific behavior and deliberate differences

- Follow Python's top-to-bottom ordering and scale meanings. Do not copy help
  typos (Python help reverses row/column selector notation); `[columns]` and
  `{rows}` remain AFNI conventions.
- A zero mask value means censored. Combining masks/index lists marks the union.
  The plotted curve includes censored samples; only AC statistics remove them.
  Masks, run lengths, and index lists refer to data **after selection**. Apply
  matching selectors to auxiliary files yourself; no implicit remapping occurs.
- Intervals use midpoint boundaries and extrapolated half-spacing at endpoints.
  With one sample the span is x ± 0.5. Censoring and patches require increasing x;
  unlike C's `-x` limitation, explicit increasing coordinates are supported.
- Band widths use absolute values as Python does. Scaling applies to both curves
  and bands, correcting Python's current separate-band scaling behavior. Bounds
  include bands only when requested, and explicit ranges take priority.
- Zero/negative horizontal thresholds are valid. Zero scaling denominators are
  rejected. `SCALE_TO_MAX` uses the signed maximum; it is not C's `-normx`.
- Missing tokens are NaN/NA/N/A and empty TSV cells. Gaps do not renumber samples.
  Mean imputation uses each selected/transposed output column's finite values,
  occurs before scaling, and reports replacements. Infinity is unsupported.
  Entirely missing plotted series/matrices or imputation columns fail explicitly.
- Boxplots use linear percentiles, observed whiskers within 1.5 IQR, and outliers.
  They summarize full selected data even when the plot viewport changes. An
  entirely censored AC selection displays empty/n=0, not a fabricated statistic.
- All plot-associated options bind in input order before reverse/custom ordering.
  Native selection, clipboard, and save use the same plot model and renderer.
- Default height grows for larger panel counts. Explicit pixel/inch dimensions
  remain exact. The baseline default is 1500 × 1200, rather than Python's height
  formula. Native Save defaults to current view, as previously agreed.

## Acceptance workflows

`tests/cli.rs` and `scripts/render_qc_gallery.sh` exercise:

1. Enorm/outlier panels: censor mask + index union, thresholds, partial/per-panel
   limits, run patches, uncertainty, BC/AC boxplots, PNG and SVG output.
2. Six motion columns: VOLREG, reverse order, independent scales, censoring,
   boxplots, automatic height.
3. Grayplot top strip: margin_off, overlay, SCALE_TO_HLINE, fixed 0:3 range,
   blue/green curves, run patches, mask, 12 × 0.5 inches at 150 DPI.
4. Inline missing data: broken curve/band, explicit mean imputation, and basic
   scientific labels. Input/statistics tests check actual numbers and semantics.
5. Existing RSA histogram, scatter, and reordered heatmap fixture commands.

The executable name changes to `1dPetePlot`; QC JPEG filenames change to PNG/SVG.
The actual `3dRSA` and APQC command generators, HTML references, AFNI install/build
integration, and official distribution are unchanged. These are plotting
acceptance examples, not a claim that end-to-end AFNI QC integration is complete.
