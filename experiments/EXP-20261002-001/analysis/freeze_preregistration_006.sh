#!/bin/bash
# EXP-20261002-001 の事前登録を sha256 で固定する。
# **plan.md や config.yaml を触ったら必ずこれを走らせ直し、その場で shasum -c が
# 全件一致することを確かめる**（guardrails 2026-09-27, 2026-09-28）。
#
# **EXP-20261001-001 の事前登録済みファイルは書き換えない。** ここでは
# 「判定に使う値の出所」として**読むだけ**で一覧に入れ、ハッシュを控える。
set -euo pipefail
cd "$(dirname "$0")/.."

FILES=(
  # ---- この再評価そのもの
  plan.md
  config.yaml
  analysis/reevaluate_006.py
  analysis/validate_bias_budget_006.py
  analysis/freeze_preregistration_006.sh
  # ---- 新しい検査の実装と、その点検の結果（**結果を見る前に固定する**）
  ../../scripts/analysis/bias_budget.py
  work/selftest_bias_budget_006.json
  work/bias_budget_validation_006.json
  # ---- 入力（EXP-20261001-001 の確定済みの出力。**読むだけ**）
  ../EXP-20261001-001/results/eval/gene_seed_005.csv
  ../EXP-20261001-001/results/eval/criteria_005.json
  # ---- セットと対照の定義（**読むだけ。書き換えない**）
  ../EXP-20261001-001/plan.md
  ../EXP-20261001-001/config.yaml
  ../EXP-20261001-001/analysis/evaluate_005.py
  ../EXP-20261001-001/work/genes_005.json
  ../EXP-20261001-001/work/subsets_005.json
  ../EXP-20261001-001/work/null_pool_005.json
  ../EXP-20261001-001/work/power_005.json
  ../EXP-20261001-001/work/SHA256SUMS-preregistered.txt
  # ---- 検査を作り直す根拠になった点検（2026-10-02）
  ../EXP-20261001-001/posthoc-m5/audit_part1.json
  ../EXP-20261001-001/posthoc-m5/audit_part3.json
  ../EXP-20261001-001/posthoc-m5/audit_part4.json
  # ---- 仮説
  ../../hypotheses/HYP-20261001-01.md
)
# **含めないもの**（理由つき）:
#   results/criteria_006.json   これから作る出力（固定したあとに回す）
#   results/env.json            実行時に作る
#   posthoc-m5/audit_part2.json 掃引の生データ（part1/3/4 に要約が入っている。大きいだけ）

shasum -a 256 "${FILES[@]}" > work/SHA256SUMS-preregistered.txt
echo "固定した ${#FILES[@]} ファイル → work/SHA256SUMS-preregistered.txt"
shasum -a 256 -c work/SHA256SUMS-preregistered.txt
