#!/usr/bin/env bash
# EXP-20261011-003 の事前登録を sha256 で固定する（plan.md §13）。
# **本番（埋め込み embed_a6.py）の前に走らせる。** plan.md / config.yaml / 事前確認のコードと出力を固定する。
# 固定の前に plan.md / config.yaml を触ったら、変更ログ・チェックリスト・この固定の 3 点をまとめてやり直す
# （guardrail 2026-09-27）。
# 実行: bash experiments/EXP-20261011-003/analysis/freeze_preregistration_a6.sh
set -eu
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
cd "$SCRIPT_DIR/../../.."
E=experiments/EXP-20261011-003
OUT=$E/analysis/SHA256SUMS-preregistered.txt

# 固定するファイル（plan.md・config.yaml・事前確認のコードと出力・着手前チェックの設定・陽性対照の材料・
# 埋め込みと評価のコード。**結果を見る前（2026-10-11）に smoke を通した後で追加した。実行指示どおり**）
FILES=(
  "$E/plan.md"
  "$E/config.yaml"
  "$E/analysis/precheck_a6.py"
  "$E/analysis/embed_a6.py"
  "$E/analysis/evaluate_a6.py"
  "$E/work/precheck_a6.json"
  "$E/work/axis-A6-interstitial.yaml"
  "$E/work/axis-A6-tubular.yaml"
  "outputs/literature/2026-10-11-tubular-degenerative-positive-control/candidate_markers.csv"
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