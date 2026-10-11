#!/bin/bash
# EXP-20260930-001 の事前登録を sha256 で固定する。
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
  analysis/gene_sets_004.json
  analysis/prepare_004.py
  analysis/state_rows_004.py
  analysis/centroids_from_embs_004.py
  analysis/select_genes_004.py
  analysis/precheck_eligibility_004.py
  analysis/power_004.py
  analysis/make_configs_004.py
  analysis/evaluate_004.py
  analysis/merge_results_004.py
  analysis/freeze_preregistration_004.sh
  work/source_files_004.json
  work/input_states_004.json
  work/symbol_to_ensembl_004.json
  work/depth_match_summary_004.json
  work/state_rows-seed81.json
  work/state_rows-seed82.json
  work/state_rows-seed83.json
  work/axisfeas_centroids-seed81.npz
  work/axisfeas_centroids-seed82.npz
  work/axisfeas_centroids-seed83.npz
  work/gene_stats_004.csv
  work/genes_004.json
  work/null_pool_004.json
  work/subsets_004.json
  work/overlap_check_004.json
  work/precheck_coverage_004.json
  work/precheck_eligibility_004-A_300.json
  work/precheck_eligibility_004-B_500.json
  work/power_004.json
  work/axis_feasibility.json
  work/tokenized.dataset/data-00000-of-00001.arrow
  work/tokenized.dataset/dataset_info.json
  work/tokenized.dataset/state.json
  ../../scripts/geneformer/run_isp.py
  ../../scripts/analysis/movement_correction.py
  ../../scripts/analysis/axis_feasibility.py
  ../../outputs/design/2026-09-30-annotation-goal-state/procedure-v2.md
  ../../outputs/design/2026-09-30-annotation-goal-state/axis_feasibility_v3.json
  ../../outputs/design/2026-09-30-annotation-goal-state/state_rows_v2-seed81.json
  ../../outputs/design/2026-09-30-annotation-goal-state/state_rows_v2-seed82.json
  ../../outputs/design/2026-09-30-annotation-goal-state/state_rows_v2-seed83.json
  ../../outputs/literature/2026-09-30-independent-goal-state/markers/adaptive_stromal_state_aFIB_plus_MYOF_minimal.csv
  ../../outputs/literature/2026-09-30-independent-goal-state/markers/adaptive_stromal_state_aFIB_plus_MYOF_full_conserved.csv
  ../../outputs/literature/2026-09-30-independent-goal-state/markers/lake2023_fibroblast_state_markers.json
  ../../hypotheses/HYP-20260930-01.md
)
# **含めないもの**（理由つき）:
#   work/config-smoke.yaml     試し実行専用（plan.md §11）
#   results-smoke/ / results*/ 実行の出力
#   work/input.h5ad            1.6 GB。EXP-20260929-001 の build_input_003.py が data/raw から
#                              決定的に作り直せる。sha256 は work/source_files_004.json に残してある
#
# 判定に使う値の出所は plan.md §12 でたどってある。**たどった先のコードはすべて上の一覧に入れる**
# （guardrails 2026-09-28）。実験ディレクトリの外のファイルも相対パスで入れてある。

shasum -a 256 "${FILES[@]}" > work/SHA256SUMS-preregistered.txt
echo "固定した ${#FILES[@]} ファイル → work/SHA256SUMS-preregistered.txt"
shasum -a 256 -c work/SHA256SUMS-preregistered.txt
