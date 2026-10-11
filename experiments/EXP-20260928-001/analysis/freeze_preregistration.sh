#!/bin/bash
# EXP-20260928-001 の事前登録を sha256 で固定する。
# plan.md §12 の 3 点セット（変更ログ・チェックリスト・ハッシュの固定し直し）のうち (3) を行う。
# **plan.md を触ったら必ずこれを走らせ直し、その場で shasum -c が全件一致することを確かめる**
# （guardrails 2026-09-27, 2026-09-28）。
set -euo pipefail
cd "$(dirname "$0")/.."

FILES=(
  plan.md
  config.yaml
  analysis/select_genes_001.py
  analysis/select_genes_fib.py
  analysis/extra_centroids.py
  analysis/evaluate_001.py
  analysis/power_001.py
  analysis/confound_probe_001.py
  work/genes_001.json
  work/marker_sets_001.csv
  work/null_pool_001.json
  work/gene_stats_001.csv
  work/precheck_coverage_001.json
  work/depth_match_summary_001.json
  work/state_rows-seed41.json
  work/state_rows-seed42.json
  work/state_rows-seed43.json
  work/power_001.json
  work/confound_probe_001.json
  work/gene_direction_001.csv
  work/config_gene_lists_001.yaml
  work/tokenized.dataset/data-00000-of-00001.arrow
  work/tokenized.dataset/dataset_info.json
  work/tokenized.dataset/state.json
)
# work/config-smoke.yaml は試し実行専用なので**含めない**（plan.md §11.1）。

shasum -a 256 "${FILES[@]}" > work/SHA256SUMS-preregistered.txt
echo "固定した ${#FILES[@]} ファイル → work/SHA256SUMS-preregistered.txt"
shasum -a 256 -c work/SHA256SUMS-preregistered.txt
