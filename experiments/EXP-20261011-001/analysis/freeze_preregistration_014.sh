#!/usr/bin/env bash
# EXP-20261011-001 の事前登録を sha256 で固定する（plan.md §13）。
# **本番の計算の前に走らせる。** plan.md / config.yaml / 評価コード / 着手前チェックなどの事前登録物を固定する。
# 固定の前に plan.md / config.yaml を触ったら、変更ログ・チェックリスト・この固定の 3 点をまとめてやり直す
# （guardrail 2026-09-27）。
# 実行: bash experiments/EXP-20261011-001/analysis/freeze_preregistration_014.sh
set -eu
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
cd "$SCRIPT_DIR/../../.."
E=experiments/EXP-20261011-001
OUT=$E/analysis/SHA256SUMS-preregistered.txt

# 固定するファイル（plan.md・config.yaml・評価コード・事前確認のコードと出力）
FILES=(
  "$E/plan.md"
  "$E/config.yaml"
  "$E/analysis/evaluate_014.py"
  "$E/analysis/precheck_014.py"
  "$E/analysis/geometry_014.py"
  "$E/analysis/power_014.py"
  "$E/work/axis-A2-ada.yaml"
  "$E/work/axis-A2-deg.yaml"
  "$E/work/axis-A1-ada.yaml"
  "$E/work/axis-A1-deg.yaml"
  "$E/work/axis_feasibility_A2_ada.json"
  "$E/work/axis_feasibility_A2_deg.json"
  "$E/work/axis_feasibility_A1_ada.json"
  "$E/work/axis_feasibility_A1_deg.json"
  "$E/work/geometry_014.json"
  "$E/work/precheck_014.json"
  "$E/work/power_014.json"
)

missing=0
for f in "${FILES[@]}"; do
  [ -s "$f" ] || { echo "無い・空: $f" >&2; missing=1; }
done
if [ "$missing" -ne 0 ]; then
  exit 2
fi

shasum -a 256 "${FILES[@]}" > "$OUT"
shasum -a 256 -c "$OUT"
echo "固定した: ${OUT} (${#FILES[@]} 件)"