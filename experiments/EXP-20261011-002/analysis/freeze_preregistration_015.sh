#!/usr/bin/env bash
# EXP-20261011-002 の事前登録を sha256 で固定する（plan.md §13）。
# **本番の評価（evaluate_015.py）の前に走らせる。** plan.md / config.yaml / 評価コード / 事前確認の
# コードと出力を固定する。
# 固定の前に plan.md / config.yaml を触ったら、変更ログ・チェックリスト・この固定の 3 点をまとめてやり直す
# （guardrail 2026-09-27）。
# 実行: bash experiments/EXP-20261011-002/analysis/freeze_preregistration_015.sh
set -eu
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
cd "$SCRIPT_DIR/../../.."
E=experiments/EXP-20261011-002
OUT=$E/analysis/SHA256SUMS-preregistered.txt

# 固定するファイル（plan.md・config.yaml・評価コード・事前確認のコードと出力）
FILES=(
  "$E/plan.md"
  "$E/config.yaml"
  "$E/analysis/evaluate_015.py"
  "$E/analysis/precheck_015.py"
  "$E/analysis/power_015.py"
  "$E/work/axis-A1-deg23-015.yaml"
  "$E/work/axis-A2-deg23-015.yaml"
  "$E/work/axis_feasibility_A1_deg23.json"
  "$E/work/axis_feasibility_A2_deg23.json"
  "$E/work/precheck_015.json"
  "$E/work/power_015.json"
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