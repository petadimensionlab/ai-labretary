#!/bin/bash
# EXP-20261003-001 — seed 83 を回し直すため、**spark-c083 に必要なものを全部まとめて先に送る**。
#
# ## spark-c083 に固有の注意（2026-09-26 / 2026-10-02 の記録）
#
#   REMOTE_GF_HOME       = $HOME/workspace/research/Geneformer
#   REMOTE_PRETRAINED_DIR= $HOME/workspace/research/Geneformer/geneformer_hf/Geneformer-V2-316M
#
# **この 2 つを指定せずに push すると、モデルへのリンクが宙に浮いて壊れる**
# （spark-c083 では重みが geneformer_hf/ 配下にあるため。2026-10-02 に実際に壊した）。
# `isp_remote.sh push` を素で叩かず、**必ずこのスクリプト経由で送る。**
# push の中で `model.safetensors` の sha256 を Mac と照合するので、
# **リンクが正しく張れていなければここで止まる。**
#
# また `isp_remote.sh push` は `state_embs.rows_file` を送らない。送り忘れると起動できない。
# → **設定と rows ファイルを 1 回でまとめて送り、sha256 で照合する。**
#   送った後は、実験が終わるまで push を走らせない。
#
# 使い方:
#   bash experiments/EXP-20261003-001/analysis/stage_remote_008_c083.sh coi-next1@100.94.57.17
set -euo pipefail
cd "$(dirname "$0")/../../.."

HOST="${1:?usage: stage_remote_008_c083.sh HOST}"
RROOT="labretary-remote"
ED="experiments/EXP-20261003-001"
CFG="$ED/work/config-seed83-c083.yaml"

export REMOTE_GF_HOME='$HOME/workspace/research/Geneformer'
export REMOTE_PRETRAINED_DIR='$HOME/workspace/research/Geneformer/geneformer_hf/Geneformer-V2-316M'

echo "=== $HOST に seed 83 の一式を送る ==="
echo "REMOTE_GF_HOME=$REMOTE_GF_HOME"
echo "REMOTE_PRETRAINED_DIR=$REMOTE_PRETRAINED_DIR"

# 1) push（run_isp.py・capture_env.sh・この設定・トークン化済み入力・モデルのリンク）
#    **トークン化済み入力は EXP-20261001-001 のもの**（作り直すと行の並びが変わる）
bash scripts/remote/isp_remote.sh push "$HOST" "$CFG"

# 2) **rows ファイルを同じ呼び出しのうちに送る**（push は送らない）
f="$ED/work/state_rows-seed83.json"
ssh -o BatchMode=yes "$HOST" "mkdir -p $RROOT/$(dirname "$f") && cat > $RROOT/$f.tmp && mv $RROOT/$f.tmp $RROOT/$f" < "$f"
want="$(shasum -a 256 "$f" | awk '{print $1}')"
got="$(ssh -o BatchMode=yes "$HOST" "sha256sum $RROOT/$f" | awk '{print $1}')"
if [[ "$want" == "$got" ]]; then echo "OK: $f"; else echo "NG: $f の sha256 が一致しません"; exit 1; fi

# 3) **モデルのリンクが生きていることを、重みの実体まで見て確かめる**（2026-10-02 の教訓）
echo "--- モデルのリンクの確認 ---"
ssh -o BatchMode=yes "$HOST" "cd $RROOT && ls -la models/Geneformer/ && ls -l models/Geneformer/Geneformer-V2-316M/model.safetensors"

# 4) 送れたことを機体の側で確かめる（ls で実在を見る。guardrails）
echo "--- $HOST で実在を確認 ---"
ssh -o BatchMode=yes "$HOST" "cd $RROOT && ls -d scripts/geneformer/run_isp.py $CFG $f \
  experiments/EXP-20261001-001/work/tokenized.dataset models/Geneformer/Geneformer-V2-316M/config.json"

echo
echo "=== 完了。**この後、実験が終わるまで isp_remote.sh push を走らせないこと** ==="
echo "次: bash scripts/remote/isp_remote.sh run $HOST $CFG"
