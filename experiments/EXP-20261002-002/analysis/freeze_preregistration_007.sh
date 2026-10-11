#!/bin/bash
# EXP-20261002-002 の事前登録を sha256 で固定する。
# **plan.md や config.yaml を触ったら必ずこれを走らせ直し、その場で shasum -c が
# 全件一致することを確かめる**（guardrails 2026-09-27, 2026-09-28）。
#
# **EXP-20261001-001 と EXP-20261002-001 の事前登録済みファイルは書き換えない。**
# ここでは「判定に使う値の出所」として**読むだけ**で一覧に入れ、ハッシュを控える。
set -euo pipefail
cd "$(dirname "$0")/.."

FILES=(
  # ---- この実験そのもの
  plan.md
  config.yaml
  config-seed81.yaml
  config-seed82.yaml
  config-seed83.yaml
  analysis/select_extra_controls_007.py
  analysis/make_configs_007.py
  analysis/evaluate_007.py
  analysis/machine_check_007.py
  analysis/stage_remote_007.sh
  analysis/freeze_preregistration_007.sh
  # ---- 追加する対照の決定と、その確認（**結果を見る前に固定する**）
  work/extra_null_007.json
  work/overlap_check_007.json
  work/precheck_coverage_007.json
  # ---- 状態の重心の行（EXP-20261001-001 からの複製。sha256 が同一であること）
  work/state_rows-seed81.json
  work/state_rows-seed82.json
  work/state_rows-seed83.json
  # ---- 健全性の検査の実装（**EXP-20261002-001 で固定したものをそのまま使う**）
  ../../scripts/analysis/bias_budget.py
  ../../scripts/analysis/movement_correction.py
  ../EXP-20261002-001/plan.md
  ../EXP-20261002-001/config.yaml
  ../EXP-20261002-001/work/selftest_bias_budget_006.json
  ../EXP-20261002-001/work/bias_budget_validation_006.json
  ../EXP-20261002-001/results/criteria_006.json
  # ---- 既存の結果とセットの定義（**読むだけ。書き換えない**）
  ../EXP-20261001-001/plan.md
  ../EXP-20261001-001/config.yaml
  ../EXP-20261001-001/analysis/evaluate_005.py
  ../EXP-20261001-001/analysis/select_genes_005.py
  ../EXP-20261001-001/analysis/gene_sets_005.json
  ../EXP-20261001-001/work/genes_005.json
  ../EXP-20261001-001/work/subsets_005.json
  ../EXP-20261001-001/work/null_pool_005.json
  ../EXP-20261001-001/work/power_005.json
  ../EXP-20261001-001/work/gene_stats_005.csv
  ../EXP-20261001-001/work/symbol_to_ensembl_005.json
  ../EXP-20261001-001/work/state_rows-seed81.json
  ../EXP-20261001-001/work/state_rows-seed82.json
  ../EXP-20261001-001/work/state_rows-seed83.json
  ../EXP-20261001-001/results/eval/gene_seed_005.csv
  ../EXP-20261001-001/results/eval/criteria_005.json
  ../EXP-20261001-001/work/SHA256SUMS-preregistered.txt
  # ---- 仮説
  ../../hypotheses/HYP-20261001-01.md
)
# **含めないもの**（理由つき）:
#   work/config-smoke.yaml / work/config-machinecheck.yaml  試し実行・機体確認専用（事前登録の対象外）
#   results-smoke/*                                          試し実行の出力（判定ではない）
#   results/criteria_007*.json                               これから作る出力（固定したあとに回す）
#   results/env*.json                                        実行時に作る

shasum -a 256 "${FILES[@]}" > work/SHA256SUMS-preregistered.txt
echo "固定した ${#FILES[@]} ファイル → work/SHA256SUMS-preregistered.txt"
shasum -a 256 -c work/SHA256SUMS-preregistered.txt
