#!/bin/bash
# EXP-20260929-001 の事前登録を sha256 で固定する。
# plan.md §11 の 3 点セット（変更ログ・チェックリスト・ハッシュの固定し直し）のうち (3) を行う。
# **plan.md や config を触ったら必ずこれを走らせ直し、その場で shasum -c が全件一致することを確かめる**
# （guardrails 2026-09-27, 2026-09-28）。
set -euo pipefail
cd "$(dirname "$0")/.."

FILES=(
  plan.md
  config.yaml
  config-seed61.yaml
  config-seed62.yaml
  config-seed63.yaml
  analysis/excluded_symbols_003.json
  analysis/build_input_003.py
  analysis/symbol_map_003.py
  analysis/tokenize_003.py
  analysis/precheck_counts_003.py
  analysis/state_rows_003.py
  analysis/precheck_auc_extra_003.py
  analysis/precheck_eligibility_003.py
  analysis/select_genes_003.py
  analysis/make_configs_003.py
  analysis/extra_centroids_003.py
  analysis/power_003.py
  analysis/evaluate_003.py
  analysis/machine_check_003.py
  analysis/merge_results_003.py
  analysis/freeze_preregistration_003.sh
  work/input_states.json
  work/symbol_to_ensembl_003.json
  work/symbol_map_003_meta.json
  work/precheck_counts_003.json
  work/precheck_auc_extra_003.json
  work/precheck_eligibility_003-A_300.json
  work/precheck_eligibility_003-B_500.json
  work/depth_match_summary_003.json
  work/state_rows-seed61.json
  work/state_rows-seed62.json
  work/state_rows-seed63.json
  work/gene_stats_003.csv
  work/genes_003.json
  work/null_pool_003.json
  work/subsets_003.json
  work/precheck_coverage_003.json
  work/power_003.json
  work/config-precheck.yaml
  work/config-machinecheck.yaml
  work/tokenized.dataset/data-00000-of-00001.arrow
  work/tokenized.dataset/dataset_info.json
  work/tokenized.dataset/state.json
)
# work/config-smoke.yaml は試し実行専用なので**含めない**（plan.md §11）。
# results-smoke/ / results-precheck/ / results*/ / work/input.h5ad（1.6 GB）も含めない。
# work/input.h5ad は analysis/build_input_003.py が data/raw から決定的に作り直せる。

shasum -a 256 "${FILES[@]}" > work/SHA256SUMS-preregistered.txt
echo "固定した ${#FILES[@]} ファイル → work/SHA256SUMS-preregistered.txt"
shasum -a 256 -c work/SHA256SUMS-preregistered.txt
