#!/bin/bash
# EXP-20261011-001 — 1 台の機体に、その config が参照する外部ファイルを**数え上げて 1 回でまとめて送る**。
# EXP-20261004-001/analysis/stage_remote_009.sh の踏襲（plan.md §7・長時間ジョブ表の本番の行）。
#
# この実験の config が参照する外部ファイル（数え上げ）:
#   1. config 本体（isp_remote.sh push が送る）
#   2. トークン化済み入力 data.tokenized_dataset（isp_remote.sh push が送る）
#   3. モデル model.directory（isp_remote.sh push がリンク・sha256 照合する）
#   4. state_embs.rows_file（{seed} 展開後。この config は使うので SEED 指定が必須）
#   5. run_isp.py / capture_env.sh（isp_remote.sh push が送る）
#
# 使い方:
#   bash experiments/EXP-20261011-001/analysis/stage_remote_014.sh HOST CONFIG SEED ...
# 例:
#   bash .../stage_remote_014.sh petadimensionlab@100.104.171.109 experiments/EXP-20261011-001/config-seed81.yaml 81
#   bash .../stage_remote_014.sh coi-next1@100.94.57.17 experiments/EXP-20261011-001/config-seed82.yaml 82
#   bash .../stage_remote_014.sh coi-next1@100.94.57.17 experiments/EXP-20261011-001/config-seed83.yaml 83
set -euo pipefail
cd "$(dirname "$0")/../../.."

HOST="${1:?usage: stage_remote_014.sh HOST CONFIG SEED ...}"; shift
CONFIG="${1:?CONFIG}"; shift
SEEDS=("$@")
RROOT="labretary-remote"
PY=".venv/bin/python"

cfgval() { "$PY" -c "import yaml,sys;c=yaml.safe_load(open(sys.argv[1]));v=c
for k in sys.argv[2].split('.'): v=(v or {}).get(k)
print(v if v is not None else '')" "$CONFIG" "$1"; }

echo "=== $HOST に $CONFIG の一式を送る ==="
echo "--- 1) isp_remote.sh push（run_isp.py・capture_env.sh・config・トークン化済み入力・モデル） ---"
bash scripts/remote/isp_remote.sh push "$HOST" "$CONFIG"

echo "--- 2) state_embs.rows_file（${#SEEDS[@]} 本）を送って sha256 で照合 ---"
for s in "${SEEDS[@]}"; do
  f="$(cfgval state_embs.rows_file)"; f="${f//\{seed\}/$s}"
  [[ -f "$f" ]] || { echo "rows_file がありません: $f"; exit 1; }
  ssh -o BatchMode=yes "$HOST" "mkdir -p $RROOT/$(dirname "$f") && cat > $RROOT/$f.tmp && mv $RROOT/$f.tmp $RROOT/$f" < "$f"
  want="$(shasum -a 256 "$f" | awk '{print $1}')"
  got="$(ssh -o BatchMode=yes "$HOST" "sha256sum $RROOT/$f" | awk '{print $1}')"
  if [[ "$want" == "$got" ]]; then echo "OK: $f"; else echo "NG: $f の sha256 が一致しません"; exit 1; fi
done

echo "--- 3) $HOST で実在を確認（guardrails: rsync が無転送で成功することがある） ---"
TOK="$(cfgval data.tokenized_dataset)"; MODEL="$(cfgval model.directory)"
for s in "${SEEDS[@]}"; do
  ROWS="${ROWS:-} $(cfgval state_embs.rows_file | sed "s/{seed}/$s/g")"
done
ssh -o BatchMode=yes "$HOST" "cd $RROOT && ls -d scripts/geneformer/run_isp.py $CONFIG $TOK $MODEL/config.json ${ROWS:-}"

echo
echo "=== 完了。**この後、このジョブが終わるまで isp_remote.sh push を走らせないこと** ==="
echo "次: bash scripts/remote/isp_remote.sh run $HOST $CONFIG"