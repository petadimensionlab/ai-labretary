#!/bin/bash
# EXP-20261004-001 の事前登録を sha256 で固定する。
# **plan.md や config.yaml を触ったら必ずこれを走らせ直し、その場で shasum -c が
# 全件一致することを確かめる**（guardrails 2026-09-27, 2026-09-28）。
#
# **C-9（試し実行）・C-10（評価コードを試し実行の出力で最後まで通す）・C-11（機体一致確認）が
# 全部合格してから走らせること。** 固定したあとに評価コードを直すことになると、
# 事前登録の意味がなくなる（plan.md §12「固定してから摂動を回す」）。
#
# **EXP-20260928-001 / EXP-20261001-001 / EXP-20261002-002 の事前登録済みファイルは書き換えない。**
# ここでは「判定に使う値の出所」として**読むだけ**で一覧に入れ、ハッシュを控える。
set -euo pipefail
cd "$(dirname "$0")/.."

FILES=(
  # ---- この実験そのもの
  plan.md
  config.yaml
  analysis/precheck_009.py
  analysis/make_config_009.py
  analysis/dtype_resolution_009.py
  analysis/check_common_cells_009.py
  analysis/evaluate_009.py
  analysis/machine_check_009.py
  analysis/stage_remote_009.sh
  analysis/freeze_preregistration_009.sh
  # ---- 着手前チェック・着手前計算の出力（**結果を見る前に固定する**）
  work/precheck_009.json
  work/axis_feasibility_009.json
  work/dtype_resolution_009.json
  work/genes_009.json
  work/config-precheck.yaml
  # ---- 機体一致確認（C-11）の出力。**事前登録の対象ではないが、合格したことの記録として控える**
  work/machine_check_009.json
  # ---- 状態の重心の行（EXP-20261001-001 からの複製。sha256 が同一であること）
  work/state_rows-seed81.json
  work/state_rows-seed82.json
  work/state_rows-seed83.json
  # ---- 指標と補正の実装（**既存で固定したものをそのまま使う**）
  ../../scripts/analysis/movement_correction.py
  ../../scripts/analysis/compare_conditions.py
  ../../scripts/geneformer/run_isp.py
  # ---- 既存の結果とセットの定義（**読むだけ。書き換えない**）
  ../EXP-20261001-001/plan.md
  ../EXP-20261001-001/config.yaml
  ../EXP-20261001-001/analysis/evaluate_005.py
  ../EXP-20261001-001/work/genes_005.json
  ../EXP-20261001-001/work/subsets_005.json
  ../EXP-20261001-001/work/null_pool_005.json
  ../EXP-20261001-001/work/power_005.json
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
)
# **含めないもの**（理由つき）:
#   work/config-smoke.yaml / work/config-machinecheck.yaml   試し実行・機体確認専用（事前登録の対象外）
#   work/common_cells_check_009*.json                        A-2 チェックの出力（都度の実行で作る。事前登録の対象外）
#   results-smoke/* / results-machinecheck-{cf5b,c083}/*      試し実行・機体確認の出力（判定ではない）
#   results/criteria_009*.json / results/gene_*_009*.csv      これから作る出力（固定したあとに回す）
#   results/env*.json / results/device.json                   実行時に作る
#   hypothesis_id                                             この実験は null（G1 を立てない。plan.md §0.3-1）

shasum -a 256 "${FILES[@]}" > work/SHA256SUMS-preregistered.txt
echo "固定した ${#FILES[@]} ファイル → work/SHA256SUMS-preregistered.txt"
shasum -a 256 -c work/SHA256SUMS-preregistered.txt
