#!/bin/bash
# EXP-20260928-002 の事前登録を sha256 で固定する。
# plan.md §10 の 3 点セット（変更ログ・チェックリスト・ハッシュの固定し直し）のうち (3) を行う。
# **plan.md や config を触ったら必ずこれを走らせ直し、その場で shasum -c が全件一致することを確かめる**
# （guardrails 2026-09-27, 2026-09-28）。
set -euo pipefail
cd "$(dirname "$0")/.."

FILES=(
  plan.md
  config.yaml
  config-seed51.yaml
  config-seed52.yaml
  config-seed53.yaml
  config-seed51-stage1.yaml
  config-seed51-stage2.yaml
  config-seed52-stage1.yaml
  config-seed52-stage2.yaml
  config-seed53-stage1.yaml
  config-seed53-stage2.yaml
  analysis/pathway_candidates_002.tsv
  analysis/select_genes_002.py
  analysis/select_genes_fib.py
  analysis/make_configs_002.py
  analysis/extra_centroids.py
  analysis/power_002.py
  analysis/evaluate_002.py
  analysis/machine_check_002.py
  analysis/merge_results_002.py
  analysis/freeze_preregistration_002.sh
  work/gene_stats_002.csv
  work/pathway_set_002.json
  work/null_pool_002.json
  work/subsets_002.json
  work/genes_002.json
  work/stage_genes_002.json
  work/state_rows-seed51.json
  work/state_rows-seed52.json
  work/state_rows-seed53.json
  work/depth_match_summary_002.json
  work/precheck_coverage_002.json
  work/config_gene_lists_002.yaml
  work/power_002.json
  work/config-machinecheck.yaml
  work/tokenized.dataset/data-00000-of-00001.arrow
  work/tokenized.dataset/dataset_info.json
  work/tokenized.dataset/state.json
)
# work/config-smoke.yaml は試し実行専用なので**含めない**（plan.md §11）。
# results-smoke/ と results*/ も含めない。

shasum -a 256 "${FILES[@]}" > work/SHA256SUMS-preregistered.txt
echo "固定した ${#FILES[@]} ファイル → work/SHA256SUMS-preregistered.txt"
shasum -a 256 -c work/SHA256SUMS-preregistered.txt
