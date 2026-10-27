#!/bin/sh
# Synthetic acceptance examples patterned after lib_apqc_tcsh.py.
set -eu
PETEPLOT_ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
PETEPLOT_BIN=${1:-"$PETEPLOT_ROOT/target/release/1dPetePlot"}
PETEPLOT_OUT="$PETEPLOT_ROOT/gallery"
mkdir -p "$PETEPLOT_OUT"
for PETEPLOT_EXT in png svg; do
  "$PETEPLOT_BIN" -infiles "$PETEPLOT_ROOT/examples/qc.tsv[enorm,outlier]" \
    -yfiles_pm "$PETEPLOT_ROOT/examples/qc.tsv[se_enorm,se_outlier]" -ylim_use_pm \
    -censor_files "$PETEPLOT_ROOT/examples/qc_censor.1D" -censor_trs '0..2' \
    -censor_hline 0.2 0.05 -patches 60 60 -boxplot_on -sepscl \
    -yaxis 0:0.6 0:0.1 -ylabels 'enorm (~mm)' 'outlier fraction' -xlabel 'Volume index' \
    -title 'Motion and outlier QC · synthetic data' \
    -panel-title 1 'Motion Euclidean norm' -panel-title 2 'Outlier fraction' \
    -prefix "$PETEPLOT_OUT/qc-censor.$PETEPLOT_EXT" -overwrite
  "$PETEPLOT_BIN" -infiles "$PETEPLOT_ROOT/examples/qc.tsv[enorm,outlier]" \
    -margin_off -one_graph -figsize 12 0.5 -dpi 150 -patches 60 60 \
    -scale SCALE_TO_HLINE -censor_hline 0.2 0.05 -yaxis 0:3 -colors blue green \
    -censor_files "$PETEPLOT_ROOT/examples/qc_censor.1D" \
    -prefix "$PETEPLOT_OUT/qc-strip.$PETEPLOT_EXT" -overwrite
  "$PETEPLOT_BIN" -infiles "$PETEPLOT_ROOT/examples/qc_motion.1D" \
    -sepscl -boxplot_on -patches 60 60 -reverse_order -ylabels VOLREG \
    -censor_files "$PETEPLOT_ROOT/examples/qc_censor.1D" -title 'Motion parameters · synthetic data' \
    -prefix "$PETEPLOT_OUT/qc-motion.$PETEPLOT_EXT" -overwrite
  "$PETEPLOT_BIN" -infiles '1D: 0 1 2 NaN NaN 1 0 1 2' -missing gap \
    -yfiles_pm '1D: 9@0.25' -ylim_use_pm -title 'Missing observations remain gaps · synthetic data' \
    -ylabels 'β_{model}' -xlabel 'Sample' -prefix "$PETEPLOT_OUT/qc-gaps.$PETEPLOT_EXT" -overwrite
  "$PETEPLOT_BIN" -infiles "$PETEPLOT_ROOT/examples/qc_motion.1D[0]" \
    -shade 3:11 -shade-color purple -shade-label Practice \
    -shade-xmat "$PETEPLOT_ROOT/examples/X.demo.mat.1D" -shade-xmat-colors blue orange \
    -censor_trs '23..25,82..84' -boxplot_on -title 'Condition shading · synthetic data' \
    -panel-title 1 'Independent annotations plus AFNI model activity' -xlabel 'Volume index' \
    -prefix "$PETEPLOT_OUT/shading.$PETEPLOT_EXT" -overwrite

done
