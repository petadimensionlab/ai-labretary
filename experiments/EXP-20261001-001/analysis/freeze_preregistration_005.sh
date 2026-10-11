#!/bin/bash
# EXP-20261001-001 の事前登録を sha256 で固定する。
# plan.md §11 の 3 点セット（変更ログ・チェックリスト・ハッシュの固定し直し）のうち (3) を行う。
# **plan.md や config を触ったら必ずこれを走らせ直し、その場で shasum -c が全件一致することを確かめる**
# （guardrails 2026-09-27, 2026-09-28）。
set -euo pipefail
cd "$(dirname "$0")/.."

FILES=(
  plan.md
  config.yaml
  config-seed81.yaml
  config-seed82.yaml
  config-seed83.yaml
  analysis/gene_sets_005.json
  analysis/prepare_005.py
  analysis/state_rows_005.py
  analysis/select_genes_005.py
  analysis/precheck_eligibility_005.py
  analysis/power_005.py
  analysis/make_configs_005.py
  analysis/evaluate_005.py
  analysis/merge_results_005.py
  analysis/freeze_preregistration_005.sh
  work/source_files_005.json
  work/input_states_005.json
  work/symbol_to_ensembl_005.json
  work/depth_match_summary_005.json
  work/state_rows-seed81.json
  work/state_rows-seed82.json
  work/state_rows-seed83.json
  work/axisfeas_centroids-seed81.npz
  work/axisfeas_centroids-seed82.npz
  work/axisfeas_centroids-seed83.npz
  work/gene_stats_005.csv
  work/genes_005.json
  work/null_pool_005.json
  work/subsets_005.json
  work/neutral_rule_005.json
  work/overlap_check_005.json
  work/precheck_coverage_005.json
  work/precheck_eligibility_005-A_300.json
  work/precheck_eligibility_005-B_500.json
  work/power_005.json
  work/axis_feasibility.json
  work/tokenized.dataset/data-00000-of-00001.arrow
  work/tokenized.dataset/dataset_info.json
  work/tokenized.dataset/state.json
  ../../scripts/geneformer/run_isp.py
  ../../scripts/analysis/movement_correction.py
  ../../scripts/analysis/axis_feasibility.py
  ../../outputs/design/2026-10-01-intervention-gene-set/procedure.md
  ../../outputs/design/2026-10-01-intervention-gene-set/counts.json
  ../../outputs/design/2026-10-01-intervention-gene-set/count_sets.py
  ../../outputs/design/2026-10-01-intervention-gene-set/axis_feasibility.json
  ../../outputs/design/2026-10-01-intervention-gene-set/axis_feasibility_config.yaml
  ../../experiments/EXP-20260928-002/work/pathway_set_002.json
  ../../outputs/literature/2026-09-30-independent-goal-state/markers/adaptive_stromal_state_aFIB_plus_MYOF_minimal.csv
  ../../outputs/literature/2026-09-30-independent-goal-state/markers/adaptive_stromal_state_aFIB_plus_MYOF_full_conserved.csv
  ../../outputs/literature/2026-09-30-independent-goal-state/markers/lake2023_fibroblast_state_markers.json
  ../../hypotheses/HYP-20261001-01.md
)
# **含めないもの**（理由つき）:
#   work/config-smoke.yaml     試し実行専用（plan.md §11）
#   results-smoke/ / results*/ 実行の出力
#   work/input.h5ad            この実験では作らない（トークン化済み入力を複製して使う。
#                              複製元の sha256 は work/source_files_005.json に残してある）
#
# 判定に使う値の出所は plan.md §12 でたどってある。**たどった先のコードはすべて上の一覧に入れる**
# （guardrails 2026-09-28）。実験ディレクトリの外のファイルも相対パスで入れてある。
#   - counts.json / pathway_set_002.json は**介入セットそのもの**の出どころなので必須
#   - movement_correction.py は主指標の交絡除去（use_dnorm=False で呼ぶ）
#   - axis_feasibility.py は着手前チェック

shasum -a 256 "${FILES[@]}" > work/SHA256SUMS-preregistered.txt
echo "固定した ${#FILES[@]} ファイル → work/SHA256SUMS-preregistered.txt"
shasum -a 256 -c work/SHA256SUMS-preregistered.txt
