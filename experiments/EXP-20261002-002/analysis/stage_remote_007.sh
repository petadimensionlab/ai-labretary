#!/bin/bash
# EXP-20261002-002 — 1 台の機体に**必要なものを全部まとめて先に送る**。
#
# ## なぜこのスクリプトがあるか（2026-10-02 の教訓）
#
# EXP-20261001-001 の seed83 は起動に 2 回失敗した。
#   1 回目: その seed の **config を送っていなかった**ので FileNotFoundError
#   2 回目: 実験の途中で `isp_remote.sh push` を走らせ、**モデルへのリンクを壊した**ので OSError
#
# `isp_remote.sh push` が送るのは run_isp.py・capture_env.sh・**その 1 本の config**・
# トークン化済み入力・モデルだけである。**`state_embs.rows_file` は送らない。**
# この実験は rows_file を使う（状態の重心を EXP-20261001-001 と 1 ビットも違わないものにするため）
# ので、送り忘れると起動できない。
#
# → **1 台が担当する seed の config を全部と、rows ファイル 3 本を、開始前に 1 回で送る。**
#   送った後は実験が終わるまで push を走らせない。
#
# 使い方:
#   bash experiments/EXP-20261002-002/analysis/stage_remote_007.sh HOST SEED [SEED ...]
# 例:
#   bash .../stage_remote_007.sh petadimensionlab@100.104.171.109 81
#   bash .../stage_remote_007.sh coi-next1@100.94.57.17 82 83     # ← 2 seed 分を 1 回で
set -euo pipefail
cd "$(dirname "$0")/../../.."

HOST="${1:?usage: stage_remote_007.sh HOST SEED [SEED ...]}"; shift
SEEDS=("$@")
[[ ${#SEEDS[@]} -ge 1 ]] || { echo "seed を 1 つ以上指定すること"; exit 1; }
RROOT="labretary-remote"
ED="experiments/EXP-20261002-002"

echo "=== $HOST に seed ${SEEDS[*]} の一式を送る ==="

# 1) 1 本目の config で push（run_isp.py・capture_env.sh・トークン化済み入力・モデル）
#    **トークン化済み入力は EXP-20261001-001 のもの**（作り直すと行の並びが変わる危険がある）
bash scripts/remote/isp_remote.sh push "$HOST" "$ED/config-seed${SEEDS[0]}.yaml"

# 2) **残りの seed の config と rows ファイルを、同じ呼び出しのうちに全部送る**
FILES=()
for s in "${SEEDS[@]}"; do
  FILES+=("$ED/config-seed${s}.yaml" "$ED/work/state_rows-seed${s}.json")
done
for f in "${FILES[@]}"; do
  ssh -o BatchMode=yes "$HOST" "mkdir -p $RROOT/$(dirname "$f") && cat > $RROOT/$f.tmp && mv $RROOT/$f.tmp $RROOT/$f" < "$f"
  want="$(shasum -a 256 "$f" | awk '{print $1}')"
  got="$(ssh -o BatchMode=yes "$HOST" "sha256sum $RROOT/$f" | awk '{print $1}')"
  if [[ "$want" == "$got" ]]; then echo "OK: $f"; else echo "NG: $f の sha256 が一致しません"; exit 1; fi
done

# 3) 送れたことを機体の側で確かめる（ls で実在を見る。guardrails）
echo "--- $HOST で実在を確認 ---"
ssh -o BatchMode=yes "$HOST" "cd $RROOT && ls -d scripts/geneformer/run_isp.py ${FILES[*]} \
  experiments/EXP-20261001-001/work/tokenized.dataset models/Geneformer/Geneformer-V2-316M/config.json"

echo
echo "=== 完了。**この後、実験が終わるまで isp_remote.sh push を走らせないこと** ==="
echo "次: bash scripts/remote/isp_remote.sh run $HOST $ED/config-seed${SEEDS[0]}.yaml"
