#!/usr/bin/env bash
set -euo pipefail

mvpa=$1
createico=$2
convertdset=$3
dinfo=$4
work=$(mktemp -d "/tmp/thd_mvpa_surface.XXXXXX")
trap 'rm -rf "$work"' EXIT

"$createico" -ld 1 -rad 10 -prefix "$work/mesh.gii" >/dev/null 2>&1

cat > "$work/nodes.1D" <<'EOF'
0
2
5
9
EOF
cat > "$work/values.1D" <<'EOF'
0 10 1 9
0 10 1 9
5 5 5 5
5 5 5 5
EOF
cat > "$work/mask.1D" <<'EOF'
1
1
1
1
EOF

"$convertdset" -o_niml -input "$work/values.1D" \
  -node_index_1D "$work/nodes.1D" -dset_labels 'a1 b1 a2 b2' \
  -prefix "$work/samples.niml.dset"
"$convertdset" -o_niml -input "$work/mask.1D" \
  -node_index_1D "$work/nodes.1D" -dset_labels mask \
  -prefix "$work/mask.niml.dset"

cat > "$work/samples.tsv" <<EOF
Subj Run Sample Condition InputFile
s01 r1 a1 A $work/samples.niml.dset[0]
s01 r1 b1 B $work/samples.niml.dset[1]
s01 r2 a2 A $work/samples.niml.dset[2]
s01 r2 b2 B $work/samples.niml.dset[3]
EOF

OMP_NUM_THREADS=2 "$mvpa" -mode classify -dataTableFile "$work/samples.tsv" \
  -target_column Condition -cv_column Run -mask "$work/mask.niml.dset" \
  -surf "$work/mesh.gii" -searchlight 100 \
  -nperm 99 -permute_within Run -seed 42 \
  -prefix "$work/masked.niml.dset" \
  -save_predictions "$work/masked.predictions.tsv" \
  -save_fold_scores "$work/masked.folds.tsv" \
  > "$work/masked.audit.txt" 2>&1

grep -q 'geodesic surface searchlight, radius 100 mm' "$work/masked.audit.txt"
grep -q '4 centers, at most 4 features, 2 workers' "$work/masked.audit.txt"
grep -q '4 synchronized label relabelings (exact; 2 blocks from Run)' \
  "$work/masked.audit.txt"
test "$("$dinfo" -nv "$work/masked.niml.dset")" = 9

"$convertdset" -o_1Dp_stdout -input "$work/masked.niml.dset" \
  -prepend_node_index_1D > "$work/masked.1D"
awk '
  NF {
    n++
    if ($1 != id[n] || $2 != 1 || $3 != .5 || $4 != .5 ||
        $6 != .5 || $7 != .5 || $9 != 4 || $10 != 2) exit 1
  }
  BEGIN { id[1]=0; id[2]=2; id[3]=5; id[4]=9 }
  END { if (n != 4) exit 1 }
' "$work/masked.1D"
awk -F '\t' '
  /^#/ { next }
  $1 == "UnitType" { next }
  {
    n++; seen[$2]=1
    if ($1 != "Node" || $4 != "NA" || $5 != "NA" || $6 != "NA" ||
        $7 != 4 || $16 != $17 || $18 != 1) exit 1
  }
  END { if (n != 16 || !seen[0] || !seen[2] || !seen[5] || !seen[9]) exit 1 }
' "$work/masked.predictions.tsv"
awk -F '\t' '
  /^#/ { next }
  $1 == "UnitType" { next }
  { n++; if ($1 != "Node" || $7 != 4 || $10 != 2 || $11 != 2 || $12 != 1 || $13 != 1) exit 1 }
  END { if (n != 8) exit 1 }
' "$work/masked.folds.tsv"

OMP_NUM_THREADS=1 "$mvpa" -mode classify -dataTableFile "$work/samples.tsv" \
  -target_column Condition -cv_column Run -surf "$work/mesh.gii" \
  -searchlight 100 -prefix "$work/nomask.niml.dset" \
  -nperm 99 -permute_within Run -seed 42 \
  -save_predictions "$work/nomask.predictions.tsv" \
  -save_fold_scores "$work/nomask.folds.tsv" -quiet
"$convertdset" -o_1Dp_stdout -input "$work/nomask.niml.dset" \
  -prepend_node_index_1D > "$work/nomask.1D"
cmp "$work/masked.1D" "$work/nomask.1D"
cmp "$work/masked.predictions.tsv" "$work/nomask.predictions.tsv"
cmp "$work/masked.folds.tsv" "$work/nomask.folds.tsv"

cat > "$work/ridge_class.tsv" <<EOF
Subj Run Sample Condition InputFile
s01 r1 a1 A $work/samples.niml.dset[0]
s01 r1 b1 B $work/samples.niml.dset[1]
s01 r2 a2 A $work/samples.niml.dset[2]
s01 r2 b2 B $work/samples.niml.dset[3]
s01 r3 a3 A $work/samples.niml.dset[0]
s01 r3 b3 B $work/samples.niml.dset[1]
EOF
for workers in 1 2; do
  OMP_NUM_THREADS=$workers "$mvpa" -mode classify -method ridge \
    -dataTableFile "$work/ridge_class.tsv" -target_column Condition -cv_column Run \
    -mask "$work/mask.niml.dset" -surf "$work/mesh.gii" -searchlight 100 \
    -lambda_grid 10,1 -nperm 99 -permute_within Run -seed 42 \
    -prefix "$work/ridge_class_$workers.niml.dset" \
    -save_predictions "$work/ridge_class_$workers.pred.tsv" \
    -save_fold_scores "$work/ridge_class_$workers.fold.tsv"
  "$convertdset" -o_1Dp_stdout -input "$work/ridge_class_$workers.niml.dset" \
    -prepend_node_index_1D > "$work/ridge_class_$workers.1D"
done
cmp "$work/ridge_class_1.1D" "$work/ridge_class_2.1D"
cmp "$work/ridge_class_1.pred.tsv" "$work/ridge_class_2.pred.tsv"
cmp "$work/ridge_class_1.fold.tsv" "$work/ridge_class_2.fold.tsv"
cmp "$work/ridge_class_1.niml.dset.tuning.tsv" "$work/ridge_class_2.niml.dset.tuning.tsv"
awk 'BEGIN{id[1]=0;id[2]=2;id[3]=5;id[4]=9}
  NF{n++;if($1!=id[n] || $2!=1 || $4!=.25 || $7!=.25 || $10!=3)exit 1}
  END{if(n!=4)exit 1}' "$work/ridge_class_1.1D"
"$mvpa" -mode classify -method ridge -dataTableFile "$work/ridge_class.tsv" \
  -target_column Condition -cv_column Run -mask "$work/mask.niml.dset" \
  -surf "$work/mesh.gii" -searchlight 100 -lambda_grid 10,1 \
  -prefix "$work/ridge_observed.niml.dset" -save_predictions "$work/ridge_observed.pred.tsv"
cmp "$work/ridge_observed.pred.tsv" "$work/ridge_class_1.pred.tsv"
cmp "$work/ridge_observed.niml.dset.tuning.tsv" "$work/ridge_class_1.niml.dset.tuning.tsv"

for method in lda_diag lda_shrinkage; do
  for workers in 1 2; do
    OMP_NUM_THREADS=$workers "$mvpa" -mode classify -method "$method" \
      -dataTableFile "$work/ridge_class.tsv" -target_column Condition -cv_column Run \
      -mask "$work/mask.niml.dset" -surf "$work/mesh.gii" -searchlight 100 \
      -nperm 99 -permute_within Run -seed 42 -prefix "$work/${method}_$workers.niml.dset" \
      -save_predictions "$work/${method}_$workers.pred.tsv" -save_fold_scores "$work/${method}_$workers.fold.tsv"
    "$convertdset" -o_1Dp_stdout -input "$work/${method}_$workers.niml.dset" \
      -prepend_node_index_1D > "$work/${method}_$workers.1D"
  done
  cmp "$work/${method}_1.1D" "$work/${method}_2.1D"
  cmp "$work/${method}_1.pred.tsv" "$work/${method}_2.pred.tsv"
  cmp "$work/${method}_1.fold.tsv" "$work/${method}_2.fold.tsv"
  awk 'BEGIN{id[1]=0;id[2]=2;id[3]=5;id[4]=9}
    NF{n++;if($1!=id[n] || $2!=1 || $4!=.25 || $7!=.25 || $10!=3)exit 1}
    END{if(n!=4)exit 1}' "$work/${method}_1.1D"
  "$mvpa" -mode classify -method "$method" -dataTableFile "$work/ridge_class.tsv" \
    -target_column Condition -cv_column Run -mask "$work/mask.niml.dset" \
    -surf "$work/mesh.gii" -searchlight 100 -prefix "$work/${method}_observed.niml.dset" \
    -save_predictions "$work/${method}_observed.pred.tsv" -save_fold_scores "$work/${method}_observed.fold.tsv"
  cmp "$work/${method}_observed.pred.tsv" "$work/${method}_1.pred.tsv"
  cmp "$work/${method}_observed.fold.tsv" "$work/${method}_1.fold.tsv"
done
if "$mvpa" -mode classify -method lda_diag -dataTableFile "$work/samples.tsv" \
  -target_column Condition -cv_column Run -mask "$work/mask.niml.dset" \
  -prefix "$work/lda_insufficient.niml.dset" > "$work/lda_insufficient.txt" 2>&1; then
  echo "FAIL LDA accepted zero residual covariance degrees of freedom" >&2
  exit 1
fi
grep -q 'more training rows than classes' "$work/lda_insufficient.txt"

cat > "$work/continuous.tsv" <<EOF
Subj Run Sample Rating InputFile
s01 r1 a1 0 $work/samples.niml.dset[0]
s01 r1 b1 10 $work/samples.niml.dset[1]
s01 r2 a2 1 $work/samples.niml.dset[2]
s01 r2 b2 9 $work/samples.niml.dset[3]
EOF
OMP_NUM_THREADS=2 "$mvpa" -mode decode -dataTableFile "$work/continuous.tsv" \
  -target_column Rating -cv_column Run -scale none -mask "$work/mask.niml.dset" \
  -surf "$work/mesh.gii" -searchlight 100 -prefix "$work/decode.niml.dset" \
  -save_predictions "$work/decode.pred.tsv" -save_fold_scores "$work/decode.fold.tsv"
OMP_NUM_THREADS=1 "$mvpa" -mode decode -dataTableFile "$work/continuous.tsv" \
  -target_column Rating -cv_column Run -scale none -mask "$work/mask.niml.dset" \
  -surf "$work/mesh.gii" -searchlight 100 -prefix "$work/decode_one.niml.dset" \
  -save_predictions "$work/decode_one.pred.tsv" -save_fold_scores "$work/decode_one.fold.tsv"
cmp "$work/decode.pred.tsv" "$work/decode_one.pred.tsv"
cmp "$work/decode.fold.tsv" "$work/decode_one.fold.tsv"
cmp "$work/decode.niml.dset.continuous.tsv" "$work/decode_one.niml.dset.continuous.tsv"
"$convertdset" -o_1Dp_stdout -input "$work/decode.niml.dset" \
  -prepend_node_index_1D > "$work/decode.1D"
awk 'BEGIN{id[1]=0;id[2]=2;id[3]=5;id[4]=9} NF{n++;if($1!=id[n] || $2<.9 || $8!=1)exit 1} END{if(n!=4)exit 1}' "$work/decode.1D"

"$mvpa" -mode encode -dataTableFile "$work/continuous.tsv" \
  -design_columns Rating -cv_column Run -scale none -mask "$work/mask.niml.dset" \
  -score voxel -prefix "$work/encode.niml.dset" -save_predictions "$work/encode.pred.tsv"
"$convertdset" -o_1Dp_stdout -input "$work/encode.niml.dset" \
  -prepend_node_index_1D > "$work/encode.1D"
awk 'BEGIN{id[1]=0;id[2]=2;id[3]=5;id[4]=9} NF{n++;if($1!=id[n])exit 1; if(n<=2){if($2<.9 || $8!=1)exit 1}else if($2!=0 || $8!=0)exit 1} END{if(n!=4)exit 1}' "$work/encode.1D"
awk -F '\t' 'NR>1 {seen[$2]=1} END{if(!seen["node_0"] || !seen["node_2"] || !seen["node_5"] || !seen["node_9"])exit 1}' "$work/encode.pred.tsv"
OMP_NUM_THREADS=2 "$mvpa" -mode encode -dataTableFile "$work/continuous.tsv" \
  -design_columns Rating -cv_column Run -scale none -mask "$work/mask.niml.dset" \
  -surf "$work/mesh.gii" -searchlight 100 -score pattern -prefix "$work/encode_pattern.niml.dset"
awk -F '\t' '/^#/ || $1=="Unit"{next} {n++;if($3!=4 || $4!=2 || $5<.9)exit 1} END{if(n!=4)exit 1}' "$work/encode_pattern.niml.dset.continuous.tsv"

for workers in 1 2; do
  OMP_NUM_THREADS=$workers "$mvpa" -mode decode -dataTableFile "$work/continuous.tsv" \
    -target_column Rating -cv_column Run -scale none -mask "$work/mask.niml.dset" \
    -surf "$work/mesh.gii" -searchlight 100 -lambda_grid 10,1 -inner_cv_column Rating \
    -prefix "$work/tuned_$workers.niml.dset" -save_fold_scores "$work/tuned_$workers.folds.tsv"
done
cmp "$work/tuned_1.niml.dset.tuning.tsv" "$work/tuned_2.niml.dset.tuning.tsv"
cmp "$work/tuned_1.folds.tsv" "$work/tuned_2.folds.tsv"
awk -F '\t' '/^#/ || $1=="Unit"{next} {seen[$1]=1;n++} END{if(n!=16 || !seen[0] || !seen[2] || !seen[5] || !seen[9])exit 1}' "$work/tuned_1.niml.dset.tuning.tsv"
OMP_NUM_THREADS=2 "$mvpa" -mode encode -dataTableFile "$work/continuous.tsv" \
  -design_columns Rating -cv_column Run -scale none -mask "$work/mask.niml.dset" \
  -surf "$work/mesh.gii" -searchlight 100 -score pattern -lambda auto -inner_cv_column Rating \
  -prefix "$work/tuned_encode.niml.dset"
test -s "$work/tuned_encode.niml.dset.tuning.tsv"

for workers in 1 2; do
  OMP_NUM_THREADS=$workers "$mvpa" -mode decode -dataTableFile "$work/continuous.tsv" \
    -target_column Rating -cv_column Run -scale none -mask "$work/mask.niml.dset" \
    -surf "$work/mesh.gii" -searchlight 100 -lambda_grid 10,1 -inner_cv_column Rating \
    -nperm 99 -permute_within Run -seed 42 -prefix "$work/perm_$workers.niml.dset" \
    -save_fold_scores "$work/perm_$workers.folds.tsv"
done
cmp "$work/perm_1.niml.dset.inference.tsv" "$work/perm_2.niml.dset.inference.tsv"
cmp "$work/perm_1.niml.dset.tuning.tsv" "$work/tuned_1.niml.dset.tuning.tsv"
cmp "$work/perm_1.folds.tsv" "$work/tuned_1.folds.tsv"
test "$("$dinfo" -nv "$work/perm_1.niml.dset")" = 17
"$convertdset" -o_1Dp_stdout -input "$work/perm_1.niml.dset" \
  -prepend_node_index_1D > "$work/perm.1D"
awk 'BEGIN{id[1]=0;id[2]=2;id[3]=5;id[4]=9} NF{n++;if($1!=id[n])exit 1} END{if(n!=4)exit 1}' "$work/perm.1D"
awk -F '\t' '/^#/ || $1=="Unit"{next} {n++;if($6<.25 || $6>1 || $9<$6 || $12!=4)exit 1} END{if(n!=4)exit 1}' \
  "$work/perm_1.niml.dset.inference.tsv"
OMP_NUM_THREADS=2 "$mvpa" -mode encode -dataTableFile "$work/continuous.tsv" \
  -design_columns Rating -cv_column Run -scale none -mask "$work/mask.niml.dset" \
  -surf "$work/mesh.gii" -searchlight 100 -score pattern -metric rmse \
  -nperm 99 -permute_within Run -prefix "$work/perm_encode.niml.dset"
awk -F '\t' '/^#/ || $1=="Unit"{next} {n++;if($5<0 || $6<.25 || $9<$6 || $12!=4)exit 1} END{if(n!=4)exit 1}' \
  "$work/perm_encode.niml.dset.inference.tsv"

cat > "$work/wrong_nodes.1D" <<'EOF'
0
2
5
8
EOF
"$convertdset" -o_niml -input "$work/mask.1D" \
  -node_index_1D "$work/wrong_nodes.1D" -dset_labels mask \
  -prefix "$work/wrong_mask.niml.dset"
if "$mvpa" -mode classify -dataTableFile "$work/samples.tsv" \
  -target_column Condition -cv_column Run -mask "$work/wrong_mask.niml.dset" \
  -surf "$work/mesh.gii" -searchlight 100 -prefix "$work/wrong_domain.niml.dset" \
  > "$work/wrong_domain.txt" 2>&1; then
  echo "FAIL 3dMVPA accepted mismatched surface node-row domains" >&2
  exit 1
fi
grep -q 'surface node-row domains differ at row 4' "$work/wrong_domain.txt"

cat > "$work/outside_nodes.1D" <<'EOF'
0
2
5
99
EOF
"$convertdset" -o_niml -input "$work/values.1D" \
  -node_index_1D "$work/outside_nodes.1D" -dset_labels 'a1 b1 a2 b2' \
  -prefix "$work/outside_samples.niml.dset"
"$convertdset" -o_niml -input "$work/mask.1D" \
  -node_index_1D "$work/outside_nodes.1D" -dset_labels mask \
  -prefix "$work/outside_mask.niml.dset"
sed "s#samples.niml.dset#outside_samples.niml.dset#" \
  "$work/samples.tsv" > "$work/outside.tsv"
if "$mvpa" -mode classify -dataTableFile "$work/outside.tsv" \
  -target_column Condition -cv_column Run -mask "$work/outside_mask.niml.dset" \
  -surf "$work/mesh.gii" -searchlight 100 -prefix "$work/outside.niml.dset" \
  > "$work/outside.txt" 2>&1; then
  echo "FAIL 3dMVPA accepted a node outside the mesh domain" >&2
  exit 1
fi
grep -q 'node 99 outside mesh range 0..11' "$work/outside.txt"

if "$mvpa" -mode classify -dataTableFile "$work/samples.tsv" \
  -target_column Condition -cv_column Run -surf "$work/mesh.gii" \
  -searchlight 'SPHERE(5)' -prefix "$work/bad_radius.niml.dset" \
  > "$work/bad_radius.txt" 2>&1; then
  echo "FAIL 3dMVPA accepted a shaped surface neighborhood" >&2
  exit 1
fi
grep -q 'surface -searchlight radius needs one finite number > 0' "$work/bad_radius.txt"

echo "PASS 3dMVPA sparse geodesic surface searchlights and domain validation"
