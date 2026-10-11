#!/bin/bash
# EXP-20261003-001 — 1 台の機体に**必要なものを全部まとめて先に送る**。
#
# ## なぜこのスクリプトがあるか（2026-10-02 の教訓）
#
# `isp_remote.sh push` が送るのは run_isp.py・capture_env.sh・**その 1 本の config**・
# トークン化済み入力・モデルだけで、**`state_embs.rows_file` は送らない**。
# この実験は rows ファイルを使う（状態の重心を EXP-20261001-001 と 1 ビットも違わないものに
# するため）ので、送り忘れると起動できない。
#
# また、**実験の途中で push を走らせるとモデルへのリンクを壊す**（2026-10-02 に実際に壊した）。
# → **本番の config・機体確認の config・rows ファイル 3 本を、開始前に 1 回で送る。**
#   送った後は実験が終わるまで push を走らせない。
#
# 使い方:
#   bash experiments/EXP-20261003-001/analysis/stage_remote_008.sh HOST
set -euo pipefail
cd "$(dirname "$0")/../../.."

HOST="${1:?usage: stage_remote_008.sh HOST}"
RROOT="labretary-remote"
ED="experiments/EXP-20261003-001"

echo "=== $HOST に EXP-20261003-001 の一式を送る ==="

# 1) 本番の config で push（run_isp.py・capture_env.sh・トークン化済み入力・モデル）
#    **トークン化済み入力は EXP-20261001-001 のもの**（作り直すと行の並びが変わる）
bash scripts/remote/isp_remote.sh push "$HOST" "$ED/config.yaml"

# 2) **機体確認の config と rows ファイル 3 本を、同じ呼び出しのうちに全部送る**
FILES=("$ED/work/config-machinecheck.yaml")
for s in 81 82 83; do FILES+=("$ED/work/state_rows-seed${s}.json"); done
for f in "${FILES[@]}"; do
  ssh -o BatchMode=yes "$HOST" "mkdir -p $RROOT/$(dirname "$f") && cat > $RROOT/$f.tmp && mv $RROOT/$f.tmp $RROOT/$f" < "$f"
  want="$(shasum -a 256 "$f" | awk '{print $1}')"
  got="$(ssh -o BatchMode=yes "$HOST" "sha256sum $RROOT/$f" | awk '{print $1}')"
  if [[ "$want" == "$got" ]]; then echo "OK: $f"; else echo "NG: $f の sha256 が一致しません"; exit 1; fi
done

# 3) 送れたことを機体の側で確かめる（ls で実在を見る。guardrails）
echo "--- $HOST で実在を確認 ---"
ssh -o BatchMode=yes "$HOST" "cd $RROOT && ls -d scripts/geneformer/run_isp.py $ED/config.yaml ${FILES[*]} \
  experiments/EXP-20261001-001/work/tokenized.dataset models/Geneformer/Geneformer-V2-316M/config.json"

echo
echo "=== 完了。**この後、実験が終わるまで isp_remote.sh push を走らせないこと** ==="
echo "次: bash scripts/remote/isp_remote.sh run $HOST $ED/work/config-machinecheck.yaml   # C-4"
echo "    （合格したら） bash scripts/remote/isp_remote.sh run $HOST $ED/config.yaml       # 本番"
