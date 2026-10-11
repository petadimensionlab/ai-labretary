#!/bin/zsh
# EXP-20261007-001 の事前登録を sha256 で固定する。
# **PI の判断（plan.md §11.1）を反映し、試し実行と評価コードの通しが済んでから走らせる。**
# plan.md / config.yaml を触ったら、変更ログ・チェックリスト・この固定の 3 点をまとめてやり直す。
set -eu
cd "${0:A:h}/../../.."
E=experiments/EXP-20261007-001
OUT=$E/work/SHA256SUMS-preregistered.txt
files=(
  $E/plan.md $E/config.yaml
  $E/analysis/precheck_011.py $E/analysis/make_config_011.py $E/analysis/evaluate_011.py
  experiments/EXP-20261004-001/analysis/evaluate_009.py
  experiments/EXP-20261001-001/analysis/evaluate_005.py
  scripts/analysis/movement_correction.py
  scripts/geneformer/run_isp.py scripts/geneformer/prepare_h5ad.py
  $E/work/prepare_011.yaml $E/work/genes_011.json $E/work/donor_groups_011.json $E/work/precheck_011.json
  $E/work/state_rows-seed81.json $E/work/state_rows-seed82.json $E/work/state_rows-seed83.json
  outputs/design/2026-10-07-A2-feasibility/axis_feasibility_A2.json
  $E/analysis/dtype_resolution_011.py $E/work/dtype_resolution_011.json $E/work/dtype_genes_011.json
  $E/analysis/machine_check_011.py $E/work/machine_check_011.json $E/analysis/merge_results_011.py
  $E/work/config-cf5b.yaml $E/work/config-c083.yaml scripts/remote/watch_run.py
)
for f in $files; do [[ -s $f ]] || { echo "無い・空: $f" >&2; exit 2; }; done
shasum -a 256 $files > $OUT
shasum -a 256 -c $OUT
echo "固定した: $OUT（${#files} 件）"
