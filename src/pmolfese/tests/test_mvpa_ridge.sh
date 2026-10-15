#!/usr/bin/env bash
set -euo pipefail
mvpa=$1
undump=$2
dinfo=$3
dump=$4
work=$(mktemp -d /tmp/thd_mvpa_ridge.XXXXXX)
trap 'rm -rf "$work"' EXIT

make_data() {
  local name=$1
  shift
  local i=0
  for value in "$@"; do
    printf '%d 0 0 %s\n' "$i" "$value"
    i=$((i+1))
  done > "$work/$name.txt"
  "$undump" -dimen 4 2 1 -ijk -datum float -prefix "$work/$name" \
    "$work/$name.txt" >/dev/null 2>&1
}
make_data mask 1 1 2 2
make_data low -1 -2 5 5
make_data mid 0 0 5 5
make_data high 1 2 5 5
cat > "$work/samples.tsv" <<EOF
Subj Run Sample X Y Z InputFile
s1 r1 a -1 8 23 $work/low+orig
s1 r1 b 0 10 20 $work/mid+orig
s1 r1 c 1 12 17 $work/high+orig
s1 r2 d -1 8 23 $work/low+orig
s1 r2 e 0 10 20 $work/mid+orig
s1 r2 f 1 12 17 $work/high+orig
EOF
common=(-dataTableFile "$work/samples.tsv" -cv_column Run -mask "$work/mask+orig" -scale none -lambda 1)
"$mvpa" -mode decode "${common[@]}" -target_columns Y,Z \
  -prefix "$work/decode" -save_predictions "$work/pred.tsv" \
  -save_fold_scores "$work/fold.tsv"
awk -F '\t' '
  function abs(x) {return x<0?-x:x}
  /^#/ || $1=="Unit" {next}
  $1==1 {if ($4!=1 || abs($5-120/121)>1e-6 || $6!=1 || $7!=1) exit 1; signal++}
  $1==2 {if ($5!=0 || $4!=1) exit 1; flat++}
  END {if(signal!=2 || flat!=2) exit 1}
' "$work/decode.continuous.tsv"
awk -F '\t' '
  function abs(x) {return x<0?-x:x}
  NR==1 {next}
  $1==1 && $2=="Y" && $3==1 {
    if($4!="a" || $7!="r1" || abs($9-(10-20/11))>1e-6 || $10!=10) exit 1; found=1
  }
  END {if(NR!=25 || !found) exit 1}
' "$work/pred.tsv"
test "$("$dinfo" -nv "$work/decode+orig")" = 16
test "$("$dinfo" -label2index CVR2_Y "$work/decode+orig")" = 3
test "$("$dinfo" -label2index SCORE_Z "$work/decode+orig")" = 7
awk -F '\t' 'NR>1 {n++;if($5!=3 || $6!=3)exit 1} END {if(n!=8)exit 1}' "$work/fold.tsv"

"$mvpa" -mode encode "${common[@]}" -design_columns X -score pattern \
  -prefix "$work/pattern" -save_predictions "$work/encoding.tsv"
awk -F '\t' '
  function abs(x) {return x<0?-x:x}
  /^#/ || $1=="Unit" {next}
  $1==1 {if($3!=2 || $4!=2 || abs($5-8/9)>1e-6)exit 1;signal=1}
  $1==2 {if($3!=2 || $4!=0 || tolower($5)!="nan")exit 1;flat=1}
  END {if(!signal || !flat)exit 1}
' "$work/pattern.continuous.tsv"
"$mvpa" -mode encode "${common[@]}" -design_columns X -score voxel -prefix "$work/voxel"
"$dump" -noijk -mask "$work/mask+orig" "$work/voxel+orig[0,6]" > "$work/voxel.txt"
awk 'function abs(x){return x<0?-x:x} NR<=2 {if(abs($1-8/9)>1e-6 || $2!=1)exit 1} NR>2 {if($1!=0 || $2!=0)exit 1} END {if(NR!=4)exit 1}' "$work/voxel.txt"
"$mvpa" -mode decode "${common[@]}" -target_column Y -searchlight 'SPHERE(0.1)' -prefix "$work/search"
test "$("$dinfo" -nv "$work/search+orig")" = 9
"$mvpa" -mode encode "${common[@]}" -design_columns X -score pattern \
  -searchlight 'SPHERE(0.1)' -prefix "$work/encode_search"

expect_failure() {
  if "$mvpa" "$@" > "$work/error.txt" 2>&1; then
    echo 'Expected invalid continuous options to fail' >&2; exit 1
  fi
}
expect_failure -mode decode "${common[@]}" -target_column Run -prefix "$work/bad"
grep -q 'finite numeric' "$work/error.txt"
expect_failure -mode decode -dataTableFile "$work/samples.tsv" -cv_column Run \
  -mask "$work/mask+orig" -target_column Y -lambda 0 -prefix "$work/bad"
grep -q 'finite and positive' "$work/error.txt"
expect_failure -mode decode "${common[@]}" -target_column Y -nperm 10 -prefix "$work/bad"
grep -q 'requires -permute_within' "$work/error.txt"
expect_failure -mode encode "${common[@]}" -design_columns X,X -prefix "$work/bad"
grep -q 'duplicate numeric column' "$work/error.txt"

cat "$work/samples.tsv" > "$work/nested.tsv"
cat >> "$work/nested.tsv" <<EOF
s1 r3 g -1 8 23 $work/low+orig
s1 r3 h 0 10 20 $work/mid+orig
s1 r3 i 1 12 17 $work/high+orig
EOF
nested_common=(-dataTableFile "$work/nested.tsv" -cv_column Run -mask "$work/mask+orig" -scale none -lambda_grid 10,1)
"$mvpa" -mode decode "${nested_common[@]}" -target_columns Y,Z \
  -prefix "$work/nested_decode" -save_fold_scores "$work/nested.folds.tsv" -show_folds > "$work/nested.audit.txt"
grep -q $'^OuterFold\tInnerFold\tStratum\tSet\tRow\tSample$' "$work/nested.audit.txt"
awk -F '\t' '
  function abs(x){return x<0?-x:x}
  /^#/ || $1=="Unit"{next}
  {n++;if($4!=2 || $5!=2)exit 1}
  $1==1 && $6==1 {if(abs($7-1/121)>1e-6 || $8!=1)exit 1; signal++}
  $1==1 && $6==10 {if(abs($7-.25)>1e-6 || $8!=0)exit 1}
  $1==2 && $6==10 {if($7!=1 || $8!=1)exit 1; tie++}
  END{if(n!=12 || signal!=3 || tie!=3)exit 1}
' "$work/nested_decode.tuning.tsv"
awk -F '\t' '
  NR==1{if($12!="Lambda" || $13!="InnerNFold" || $14!="InnerNValid" || $15!="InnerLoss")exit 1;next}
  {n++;if($13!=2 || $14!=2)exit 1}
  $1==1{if($12!=1)exit 1}
  $1==2{if($12!=10)exit 1}
  END{if(n!=12)exit 1}
' "$work/nested.folds.tsv"
"$mvpa" -mode encode "${nested_common[@]}" -design_columns X -score pattern \
  -prefix "$work/nested_encode" -inner_cv_k 2
awk -F '\t' '
  /^#/ || $1=="Unit"{next}
  $1==1 && $6==1 {if($5!=2 || $8!=1)exit 1; signal++}
  $1==2 && $6==10 {if($5!=0 || tolower($7)!="nan" || $8!=1)exit 1; constant++}
  END{if(signal!=3 || constant!=3)exit 1}
' "$work/nested_encode.tuning.tsv"
"$mvpa" -mode encode "${nested_common[@]}" -design_columns X -score voxel -prefix "$work/nested_voxel"
OMP_NUM_THREADS=1 "$mvpa" -mode decode "${nested_common[@]}" -target_column Y \
  -searchlight 'SPHERE(0.1)' -prefix "$work/tuned_one" -save_predictions "$work/tuned_one.pred.tsv"
OMP_NUM_THREADS=3 "$mvpa" -mode decode "${nested_common[@]}" -target_column Y \
  -searchlight 'SPHERE(0.1)' -prefix "$work/tuned_three" -save_predictions "$work/tuned_three.pred.tsv"
cmp "$work/tuned_one.tuning.tsv" "$work/tuned_three.tuning.tsv"
cmp "$work/tuned_one.pred.tsv" "$work/tuned_three.pred.tsv"
cmp "$work/tuned_one.continuous.tsv" "$work/tuned_three.continuous.tsv"
# A numeric lambda and a grid cannot be combined, even with custom grouping.
expect_failure -mode decode "${common[@]}" -target_column Y -lambda_grid 1,10 \
  -inner_cv_column X -prefix "$work/custom_inner" -inner_cv_k 2
grep -q 'cannot accompany a numeric' "$work/error.txt"
"$mvpa" -mode decode -dataTableFile "$work/samples.tsv" -cv_column Run \
  -mask "$work/mask+orig" -scale none -target_column Y -lambda auto \
  -inner_cv_column X -prefix "$work/custom_inner" -inner_cv_k 2
test -s "$work/custom_inner.tuning.tsv"
expect_failure -mode decode -dataTableFile "$work/samples.tsv" -cv_column Run \
  -mask "$work/mask+orig" -target_column Y -lambda auto -prefix "$work/bad"
grep -q 'inner groups' "$work/error.txt"
expect_failure -mode decode "${nested_common[@]}" -target_column Y -inner_cv_k 3 -prefix "$work/bad"
grep -q 'need at least 3' "$work/error.txt"
expect_failure -mode decode -dataTableFile "$work/nested.tsv" -cv_column Run \
  -mask "$work/mask+orig" -target_column Y -lambda_grid 1,1 -prefix "$work/bad"
grep -q 'duplicate penalties' "$work/error.txt"

"$mvpa" -mode decode "${common[@]}" -target_columns Y,Z -nperm 99 -permute_within Run \
  -prefix "$work/continuous_perm" -save_predictions "$work/perm.pred.tsv" \
  -save_fold_scores "$work/perm.fold.tsv" > "$work/perm.log" 2>&1
grep -q '36 synchronized row relabelings (exact; 2 blocks from Run)' "$work/perm.log"
cmp "$work/pred.tsv" "$work/perm.pred.tsv"
cmp "$work/fold.tsv" "$work/perm.fold.tsv"
awk -F '\t' '
  function abs(x){return x<0?-x:x}
  /^#/ || $1=="Unit"{next}
  $1==1 {if(abs($6-2/36)>1e-7 || abs($8-4/36)>1e-7 || abs($9-$6)>1e-7 || $5<=0 || $11!=36 || $12!=36)exit 1; signal++}
  $1==2 {if($6!=1 || $8!=1 || $9!=1)exit 1; flat++}
  END{if(signal!=2 || flat!=2)exit 1}
' "$work/continuous_perm.inference.tsv"
test "$("$dinfo" -nv "$work/continuous_perm+orig")" = 32
test "$("$dinfo" -label2index PVAL_Y "$work/continuous_perm+orig")" = 16
test "$("$dinfo" -label2index ZSCR_Y "$work/continuous_perm+orig")" = 17
test "$("$dinfo" -label2index ZFWE_Z "$work/continuous_perm+orig")" = 28
"$dinfo" -verb "$work/continuous_perm+orig" > "$work/perm.info.txt" 2>&1
test "$(grep -c 'statcode = fizt' "$work/perm.info.txt")" = 4
"$mvpa" -mode decode "${common[@]}" -target_column Y -nperm 99 -permute_within Run \
  -metric MAE -prefix "$work/mae_perm"
awk -F '\t' '
  function abs(x){return x<0?-x:x}
  /^#/ || $1=="Unit"{next}
  $1==1{if(abs($6-2/36)>1e-7 || $4<=$3 || $5<=0 || $9<$6)exit 1; found=1}
  END{if(!found)exit 1}
' "$work/mae_perm.inference.tsv"
awk 'NR>1 && $5==10 {$5=8} {print}' "$work/samples.tsv" > "$work/tied_targets.tsv"
"$mvpa" -mode decode -dataTableFile "$work/tied_targets.tsv" -cv_column Run \
  -mask "$work/mask+orig" -scale none -lambda 1 -target_column Y -metric pearson \
  -nperm 99 -permute_within Run -prefix "$work/undefined_null"
awk -F '\t' '
  /^#/ || $1=="Unit"{next}
  $1==1{if(tolower($3)=="nan" || $11<=0 || $11>=36 || $6<(36-$11)/36 || $9<$6 || $12!=36)exit 1;found=1}
  END{if(!found)exit 1}
' "$work/undefined_null.inference.tsv"
"$mvpa" -mode encode "${common[@]}" -design_columns X -score pattern \
  -nperm 99 -permute_within Run -prefix "$work/encoding_perm"
awk -F '\t' '
  function abs(x){return x<0?-x:x}
  /^#/ || $1=="Unit"{next}
  $1==1{if(abs($6-2/36)>1e-7 || abs($9-$6)>1e-7)exit 1;found=1}
  $1==2{if($6!=1 || $7!=0 || $8!=1 || $9!=1 || $10!=0 || $11!=0)exit 1;undefined=1}
  END{if(!found || !undefined)exit 1}
' "$work/encoding_perm.inference.tsv"
for workers in 1 3; do
  OMP_NUM_THREADS=$workers "$mvpa" -mode decode "${nested_common[@]}" -target_column Y \
    -searchlight 'SPHERE(0.1)' -nperm 9 -permute_within Run -seed 42 \
    -prefix "$work/nested_perm_$workers" -save_predictions "$work/nested_perm_$workers.pred.tsv" \
    -save_fold_scores "$work/nested_perm_$workers.fold.tsv"
done
cmp "$work/nested_perm_1.inference.tsv" "$work/nested_perm_3.inference.tsv"
cmp "$work/nested_perm_1.tuning.tsv" "$work/nested_perm_3.tuning.tsv"
cmp "$work/nested_perm_1.pred.tsv" "$work/nested_perm_3.pred.tsv"
cmp "$work/nested_perm_1.fold.tsv" "$work/nested_perm_3.fold.tsv"
cmp "$work/tuned_one.pred.tsv" "$work/nested_perm_1.pred.tsv"
cmp "$work/tuned_one.tuning.tsv" "$work/nested_perm_1.tuning.tsv"
awk -F '\t' '/^#/ || $1=="Unit"{next}{if($6<1/9 || $6>1 || $9<$6 || $9>1 || $12!=9)exit 1;n++}END{if(n!=4)exit 1}' "$work/nested_perm_1.inference.tsv"
"$mvpa" -mode encode "${nested_common[@]}" -design_columns X -score voxel -metric RMSE \
  -nperm 9 -permute_within Run -prefix "$work/nested_voxel_perm"
test -s "$work/nested_voxel_perm.inference.tsv"

awk 'NR==1{print $0,"Class";next}{print $0,($4<0?"A":"B")}' \
  "$work/nested.tsv" > "$work/classify.tsv"
class_common=(-mode classify -method ridge -dataTableFile "$work/classify.tsv" \
  -cv_column Run -mask "$work/mask+orig" -target_column Class -scale none)
"$mvpa" "${class_common[@]}" -lambda 1 -prefix "$work/class_fixed" \
  -save_predictions "$work/class_fixed.pred.tsv" -save_fold_scores "$work/class_fixed.fold.tsv"
awk -F '\t' '/^#/ || $1=="ROI"{next}
  $1==1{if($3!=1 || $4!=1 || $5!=1)exit 1;signal=1}
  $1==2{if($3!=.5 || $5!=.5)exit 1;flat=1}
  END{if(!signal || !flat)exit 1}' "$work/class_fixed.roi.tsv"
test "$("$dinfo" -nv "$work/class_fixed+orig")" = 3
awk -F '\t' '/^#/{next} $1=="UnitType"{if($14!="Lambda")exit 1;next}
  {n++;if($14!=1 || $15!=0)exit 1} END{if(n!=6)exit 1}' "$work/class_fixed.fold.tsv"
"$mvpa" "${class_common[@]}" -lambda_grid 10,1 -prefix "$work/class_tuned" \
  -save_predictions "$work/class_tuned.pred.tsv" -save_fold_scores "$work/class_tuned.fold.tsv"
"$mvpa" "${class_common[@]}" -lambda_grid 10,1 -nperm 9 -permute_within Run -seed 42 \
  -prefix "$work/class_perm" -save_predictions "$work/class_perm.pred.tsv" \
  -save_fold_scores "$work/class_perm.fold.tsv"
cmp "$work/class_tuned.pred.tsv" "$work/class_perm.pred.tsv"
cmp "$work/class_tuned.fold.tsv" "$work/class_perm.fold.tsv"
cmp "$work/class_tuned.tuning.tsv" "$work/class_perm.tuning.tsv"
awk -F '\t' '/^#/ || $1=="Unit"{next}
  $1==1 && $8==1{if($6!=1)exit 1;signal++}
  $1==2 && $8==1{if($6!=10)exit 1;flat++}
  END{if(signal!=3 || flat!=3)exit 1}' "$work/class_tuned.tuning.tsv"
for workers in 1 3; do
  OMP_NUM_THREADS=$workers "$mvpa" "${class_common[@]}" -lambda auto \
    -searchlight 'SPHERE(0.1)' -nperm 9 -permute_within Run -seed 42 \
    -prefix "$work/class_volume_$workers" -save_predictions "$work/class_volume_$workers.pred.tsv" \
    -save_fold_scores "$work/class_volume_$workers.fold.tsv"
  "$dump" -noijk -mask "$work/mask+orig" "$work/class_volume_$workers+orig" \
    > "$work/class_volume_$workers.txt"
done
cmp "$work/class_volume_1.txt" "$work/class_volume_3.txt"
cmp "$work/class_volume_1.pred.tsv" "$work/class_volume_3.pred.tsv"
cmp "$work/class_volume_1.fold.tsv" "$work/class_volume_3.fold.tsv"
cmp "$work/class_volume_1.tuning.tsv" "$work/class_volume_3.tuning.tsv"
awk '{n++;if($3<1/9 || $3>1 || $6<$3 || $6>1)exit 1} END{if(n!=4)exit 1}' "$work/class_volume_1.txt"
expect_failure "${class_common[@]}" -distance euclidean -prefix "$work/bad_class_distance"
grep -q 'only used by centroid' "$work/error.txt"
expect_failure "${class_common[@]}" -method centroid -lambda 1 -prefix "$work/bad_centroid_lambda"
grep -q 'require -method ridge' "$work/error.txt"
expect_failure "${class_common[@]}" -metric pearson -prefix "$work/bad_class_metric"
grep -q 'accuracy or balanced_accuracy' "$work/error.txt"
expect_failure "${class_common[@]}" -lambda_grid 1,10 -inner_cv_column Class -prefix "$work/bad_inner_classes"
grep -q 'inner fold.*has no class' "$work/error.txt"
lda_common=(-mode classify -dataTableFile "$work/classify.tsv" -cv_column Run \
  -mask "$work/mask+orig" -target_column Class -scale none)
for method in lda_diag lda_shrinkage; do
  "$mvpa" "${lda_common[@]}" -method "$method" -prefix "$work/$method" \
    -save_predictions "$work/$method.pred.tsv" -save_fold_scores "$work/$method.fold.tsv"
  "$mvpa" "${lda_common[@]}" -method "$method" -nperm 99 -permute_within Run \
    -prefix "$work/${method}_perm" -save_predictions "$work/${method}_perm.pred.tsv" \
    -save_fold_scores "$work/${method}_perm.fold.tsv"
  cmp "$work/$method.pred.tsv" "$work/${method}_perm.pred.tsv"
  cmp "$work/$method.fold.tsv" "$work/${method}_perm.fold.tsv"
  awk -F '\t' '/^#/ || $1=="ROI"{next}
    $1==1{if($3!=1)exit 1;signal=1} $1==2{if($3!=.5)exit 1;flat=1}
    END{if(!signal || !flat)exit 1}' "$work/$method.roi.tsv"
  awk -F '\t' '/^#/{next} $1=="UnitType"{if($14!="Shrinkage" || $15!="CovarianceScale" || $17!="CovarianceForm")exit 1;next}
    {n++;if($15<=0 || $16<=0 || ($17!="diagonal" && $17!="primal"))exit 1}
    END{if(n!=6)exit 1}' "$work/$method.fold.tsv"
  for workers in 1 3; do
    OMP_NUM_THREADS=$workers "$mvpa" "${lda_common[@]}" -method "$method" \
      -searchlight 'SPHERE(0.1)' -nperm 9 -permute_within Run -seed 42 \
      -prefix "$work/${method}_volume_$workers" \
      -save_predictions "$work/${method}_volume_$workers.pred.tsv" \
      -save_fold_scores "$work/${method}_volume_$workers.fold.tsv"
    "$dump" -noijk -mask "$work/mask+orig" "$work/${method}_volume_$workers+orig" \
      > "$work/${method}_volume_$workers.txt"
  done
  cmp "$work/${method}_volume_1.txt" "$work/${method}_volume_3.txt"
  cmp "$work/${method}_volume_1.pred.tsv" "$work/${method}_volume_3.pred.tsv"
  cmp "$work/${method}_volume_1.fold.tsv" "$work/${method}_volume_3.fold.tsv"
  "$mvpa" "${lda_common[@]}" -method "$method" -class_prior empirical \
    -prefix "$work/${method}_empirical" -save_predictions "$work/${method}_empirical.pred.tsv"
  awk -F '\t' '/^#/ || $1=="UnitType"{next} $2==2{n++;if($17!="B")exit 1}
    END{if(n!=9)exit 1}' "$work/${method}_empirical.pred.tsv"
  expect_failure "${lda_common[@]}" -method "$method" -lambda 1 -prefix "$work/bad_lda_lambda"
  grep -q 'require -method ridge' "$work/error.txt"
  expect_failure "${lda_common[@]}" -method "$method" -distance cosine -prefix "$work/bad_lda_distance"
  grep -q 'only used by centroid' "$work/error.txt"
done
"$mvpa" "${lda_common[@]}" -method lda_shrinkage -shrinkage 1 -variance_floor .01 \
  -prefix "$work/lda_spherical" -save_fold_scores "$work/lda_spherical.fold.tsv"
awk -F '\t' '/^#/ || $1=="UnitType"{next}{n++;if($14!=1 || $17!="diagonal")exit 1} END{if(n!=6)exit 1}' \
  "$work/lda_spherical.fold.tsv"
expect_failure "${lda_common[@]}" -method lda_diag -shrinkage .5 -prefix "$work/bad_lda_shrink"
grep -q 'requires -method lda_shrinkage' "$work/error.txt"
expect_failure "${lda_common[@]}" -method lda_shrinkage -shrinkage 0 -prefix "$work/bad_lda_zero"
grep -q 'must be in' "$work/error.txt"
expect_failure "${lda_common[@]}" -method lda_diag -variance_floor 0 -prefix "$work/bad_lda_floor"
grep -q 'must be in' "$work/error.txt"
expect_failure "${class_common[@]}" -class_prior empirical -prefix "$work/bad_ridge_prior"
grep -q 'require LDA classification' "$work/error.txt"
echo 'PASS ridge decoding/encoding/classification, diagonal/shrinkage LDA, audits, maps, searchlights, and invalid options'
