#!/usr/bin/env bash
# EXP-20261007-003 段階 2 の事前登録を sha256 で固定する（plan.md §13）。
# **本番の計算の前に走らせる。** plan.md / config.yaml / analysis/evaluate_013.py を固定する。
# 固定の前に plan.md / config.yaml を触ったら、変更ログ・チェックリスト・この固定の 3 点をまとめてやり直す
# （guardrail 2026-09-27）。
# 実行: bash experiments/EXP-20261007-003/analysis/freeze_preregistration_013.sh
set -eu
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
cd "$SCRIPT_DIR/../../.."
E=experiments/EXP-20261007-003
OUT=$E/analysis/SHA256SUMS-preregistered.txt

# 固定するファイル（plan.md・config.yaml・実行スクリプト。plan.md §7・§13）
F1=$E/plan.md
F2=$E/config.yaml
F3=$E/analysis/evaluate_013.py

for f in "$F1" "$F2" "$F3"; do
  [ -s "$f" ] || { echo "無い・空: $f" >&2; exit 2; }
done

shasum -a 256 "$F1" "$F2" "$F3" > "$OUT"
shasum -a 256 -c "$OUT"
echo "固定した: ${OUT} (3 件)"
