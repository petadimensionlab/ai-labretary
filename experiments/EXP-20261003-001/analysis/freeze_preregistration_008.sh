#!/bin/bash
# EXP-20261003-001 の事前登録を sha256 で固定する。
# **plan.md や config.yaml を触ったら必ずこれを走らせ直し、その場で shasum -c が
# 全件一致することを確かめる**（guardrails 2026-09-27, 2026-09-28）。
#
# **EXP-20261001-001 / EXP-20261002-001 / EXP-20261002-002 の事前登録済みファイルは書き換えない。**
# ここでは「判定に使う値の出所」として**読むだけ**で一覧に入れ、ハッシュを控える。
set -euo pipefail
cd "$(dirname "$0")/.."

FILES=(
  # ---- この実験そのもの
  plan.md
  config.yaml
  analysis/precheck_008.py
  analysis/make_configs_008.py
  analysis/evaluate_008.py
  analysis/machine_check_008.py
  analysis/stage_remote_008.sh
  analysis/freeze_preregistration_008.sh
  # ---- 着手前チェックの出力（**結果を見る前に固定する**）
  work/precheck_008.json
  # ---- 状態の重心の行（EXP-20261001-001 からの複製。sha256 が同一であること）
  work/state_rows-seed81.json
  work/state_rows-seed82.json
  work/state_rows-seed83.json
  # ---- 指標と補正の実装（**既存で固定したものをそのまま使う**）
  ../../scripts/analysis/movement_correction.py
  ../../scripts/geneformer/run_isp.py
  # ---- 既存の結果とセットの定義（**読むだけ。書き換えない**）
  ../EXP-20261001-001/plan.md
  ../EXP-20261001-001/config.yaml
  ../EXP-20261001-001/analysis/evaluate_005.py
  ../EXP-20261001-001/work/genes_005.json
  ../EXP-20261001-001/work/subsets_005.json
  ../EXP-20261001-001/work/null_pool_005.json
  ../EXP-20261001-001/work/power_005.json
  ../EXP-20261001-001/work/symbol_to_ensembl_005.json
  ../EXP-20261001-001/work/state_rows-seed81.json
  ../EXP-20261001-001/work/state_rows-seed82.json
  ../EXP-20261001-001/work/state_rows-seed83.json
  ../EXP-20261001-001/work/SHA256SUMS-preregistered.txt
  ../EXP-20261002-002/plan.md
  ../EXP-20261002-002/config.yaml
  ../EXP-20261002-002/analysis/evaluate_007.py
  ../EXP-20261002-002/work/extra_null_007.json
  ../EXP-20261002-002/work/SHA256SUMS-preregistered.txt
  ../EXP-20261002-002/results/criteria_007.json
  ../EXP-20261002-002/results/gene_seed_007-81.csv
  ../EXP-20261002-002/results/gene_seed_007-82.csv
  ../EXP-20261002-002/results/gene_seed_007-83.csv
  # ---- 検出下限 2.73 の出どころ
  ../EXP-20260928-001/work/power_001.json
  # ---- 事実関係の出どころ（目標の状態での発現）
  ../EXP-20261001-001/posthoc-eng/q1_expression.csv
  ../EXP-20261001-001/posthoc-eng/validation.md
  # ---- 仮説
  ../../hypotheses/HYP-20261001-01.md
)
# **含めないもの**（理由つき）:
#   work/config-smoke.yaml / work/config-machinecheck.yaml  試し実行・機体確認専用（事前登録の対象外）
#   work/machine_check_008.json                             機体確認の出力（事前登録の対象外。本番前に通す門）
#   results-smoke/* / results-machinecheck/*                 試し実行・機体確認の出力（判定ではない）
#   results/criteria_008*.json                               これから作る出力（固定したあとに回す）
#   results/env*.json / results/device.json                  実行時に作る

shasum -a 256 "${FILES[@]}" > work/SHA256SUMS-preregistered.txt
echo "固定した ${#FILES[@]} ファイル → work/SHA256SUMS-preregistered.txt"
shasum -a 256 -c work/SHA256SUMS-preregistered.txt
