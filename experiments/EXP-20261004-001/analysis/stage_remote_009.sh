#!/bin/bash
# EXP-20261004-001 — 1 台の機体に、ある config が参照する外部ファイルを**数え上げて 1 回でまとめて送る**。
#
# ## なぜこのスクリプトがあるか（guardrails 2026-09-28 / 2026-10-02）
#
# `isp_remote.sh push` が送るのは run_isp.py・capture_env.sh・**その 1 本の config**・
# トークン化済み入力・モデルだけである。**`state_embs.rows_file` は送らない。**
# 送り忘れると FileNotFoundError で起動に失敗する（EXP-20261001-001 の seed83 で実際に起きた）。
#
# さらに **実験の途中で `isp_remote.sh push` を走らせるとモデルへのリンクを壊す**（guardrails 2026-10-02）。
# → **開始前に、その config が参照する外部ファイルを全部洗い出し、1 回でまとめて送る。**
#   送った後は、その機体でのジョブが終わるまで push を走らせない。
#
# この実験で config が参照しうる外部ファイル（**数え上げ**）:
#   1. config 本体（isp_remote.sh push が送る）
#   2. トークン化済み入力 data.tokenized_dataset（isp_remote.sh push が送る）
#   3. モデル model.directory（isp_remote.sh push がリンク/sha256照合する）
#   4. state_embs.rows_file（{seed} 展開後。**config-machinecheck.yaml だけが使う。
#      config-smoke.yaml・config.yaml（本番）は rows_file を持たないので対象外**）
#   5. run_isp.py / capture_env.sh（isp_remote.sh push が送る。ここでは重複して送らない）
#
# 使い方:
#   bash experiments/EXP-20261004-001/analysis/stage_remote_009.sh HOST CONFIG [SEED ...]
# 例（機体一致確認。rows_file を使うので SEED が要る）:
#   bash .../stage_remote_009.sh petadimensionlab@100.104.171.109 experiments/EXP-20261004-001/work/config-machinecheck.yaml 81
#   bash .../stage_remote_009.sh coi-next1@100.94.57.17 experiments/EXP-20261004-001/work/config-machinecheck.yaml 81
# 例（試し実行。rows_file を使わないので SEED 無指定で可）:
#   bash .../stage_remote_009.sh petadimensionlab@100.104.171.109 experiments/EXP-20261004-001/work/config-smoke.yaml
set -euo pipefail
cd "$(dirname "$0")/../../.."

HOST="${1:?usage: stage_remote_009.sh HOST CONFIG [SEED ...]}"; shift
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

HAS_ROWS="$(cfgval state_embs.rows_file)"
FILES=()
if [[ -n "$HAS_ROWS" ]]; then
  [[ ${#SEEDS[@]} -ge 1 ]] || { echo "この config は state_embs.rows_file を使う。SEED を 1 つ以上指定すること"; exit 1; }
  for s in "${SEEDS[@]}"; do
    f="${HAS_ROWS//\{seed\}/$s}"
    [[ -f "$f" ]] || { echo "rows_file がありません: $f"; exit 1; }
    FILES+=("$f")
  done
  echo "--- 2) state_embs.rows_file（${#FILES[@]} 本）を送って sha256 で照合 ---"
  for f in "${FILES[@]}"; do
    ssh -o BatchMode=yes "$HOST" "mkdir -p $RROOT/$(dirname "$f") && cat > $RROOT/$f.tmp && mv $RROOT/$f.tmp $RROOT/$f" < "$f"
    want="$(shasum -a 256 "$f" | awk '{print $1}')"
    got="$(ssh -o BatchMode=yes "$HOST" "sha256sum $RROOT/$f" | awk '{print $1}')"
    if [[ "$want" == "$got" ]]; then echo "OK: $f"; else echo "NG: $f の sha256 が一致しません"; exit 1; fi
  done
else
  echo "--- 2) この config は state_embs.rows_file を使わない（送るものなし） ---"
fi

echo "--- 3) $HOST で実在を確認（guardrails: rsync が無転送で成功することがある） ---"
TOK="$(cfgval data.tokenized_dataset)"; MODEL="$(cfgval model.directory)"
ssh -o BatchMode=yes "$HOST" "cd $RROOT && ls -d scripts/geneformer/run_isp.py $CONFIG $TOK $MODEL/config.json ${FILES[*]:-}"

echo
echo "=== 完了。**この後、このジョブが終わるまで isp_remote.sh push を走らせないこと** ==="
echo "次: bash scripts/remote/isp_remote.sh run $HOST $CONFIG"
