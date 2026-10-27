#!/bin/sh
# Reproducible synthetic examples. Intentionally replaces only these gallery files.
set -eu
PETEPLOT_ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
PETEPLOT_BIN=${1:-"$PETEPLOT_ROOT/target/release/1dPetePlot"}
PETEPLOT_OUT="$PETEPLOT_ROOT/gallery"
mkdir -p "$PETEPLOT_OUT"
for PETEPLOT_EXT in png svg; do
  "$PETEPLOT_BIN" -demo -prefix "$PETEPLOT_OUT/rsa-gallery.$PETEPLOT_EXT" -overwrite
  "$PETEPLOT_BIN" -one -infiles "$PETEPLOT_ROOT/examples/signals.1D" \
    -title 'Model evidence through time' -xlabel 'Time (s)' \
    -legend_labels 'Model A' 'Model B' -prefix "$PETEPLOT_OUT/lines.$PETEPLOT_EXT" -overwrite
  "$PETEPLOT_BIN" -histogram -histbin 16 -vline 0 -vline_hl 0.40 'Reference 0.40' \
    -infiles "$PETEPLOT_ROOT/examples/effects.1D[0]" -xlabel 'RSA association (r)' \
    -title 'Synthetic effects across regions' -prefix "$PETEPLOT_OUT/histogram.$PETEPLOT_EXT" -overwrite
  "$PETEPLOT_BIN" -scatter -diagonal -xfile "$PETEPLOT_ROOT/examples/effects.1D[0]" \
    -infiles "$PETEPLOT_ROOT/examples/effects.1D[1]" -xlabel 'Model A (r)' -ylabels 'Model B (r)' \
    -title 'Compare model associations' -prefix "$PETEPLOT_OUT/scatter.$PETEPLOT_EXT" -overwrite
  "$PETEPLOT_BIN" -heat -infiles "$PETEPLOT_ROOT/examples/rdm.1D" \
    -reorder "$PETEPLOT_ROOT/examples/order.1D" -zerocenter -cbar Reds_and_Blues_Inv \
    -title 'Reordered representational dissimilarity' -prefix "$PETEPLOT_OUT/heatmap.$PETEPLOT_EXT" -overwrite
done
