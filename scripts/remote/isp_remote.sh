#!/usr/bin/env bash
# Labretary: run_isp.py をリモート GPU（Tailscale 経由の DGX Spark 系など）で実行する。
# config は Mac と同じものを使い、実行先だけ LABRETARY_DEVICE=cuda で切り替える（device.json に記録）。
# リモートには config が参照するファイルだけを、プロジェクトと同じ相対パスで置く（~/labretary-remote/）。
#
#   bash scripts/remote/isp_remote.sh push   HOST CONFIG   # 必要なファイルを送る（tokenized dataset, モデル, run_isp.py, geneformer パッケージ）
#   bash scripts/remote/isp_remote.sh run    HOST CONFIG   # tmux で実行（切断しても続く。同じコマンドで遺伝子単位の再開）
#   bash scripts/remote/isp_remote.sh status HOST CONFIG   # ログの末尾と tmux の有無
#   bash scripts/remote/isp_remote.sh pull   HOST CONFIG   # 結果ディレクトリとログを Mac に戻す
#   bash scripts/remote/isp_remote.sh clean  HOST          # リモートの ~/labretary-remote を削除（ヒト由来データの複製を残さない）
#
# HOST 例: petadimensionlab@100.104.171.109。資格情報はここに書かない（SSH 鍵を使う）。
# 環境変数: REMOTE_GF_HOME（既定 ~/workspace/Geneformer。venv と Geneformer-V2-316M がある場所）
#           REMOTE_PRETRAINED_DIR（事前学習モデルの場所が $REMOTE_GF_HOME/Geneformer-V2-316M でないとき）
# ホスト別の設定（2026-09-26）:
#   spark-cf5b: HOST=petadimensionlab@100.104.171.109（既定のままで可）
#   spark-c083: HOST=coi-next1@100.94.57.17  REMOTE_GF_HOME='$HOME/workspace/research/Geneformer'
#               REMOTE_PRETRAINED_DIR='$HOME/workspace/research/Geneformer/geneformer_hf/Geneformer-V2-316M'
# tmux のセッション名は config 名から作るので、同じホストで別 config を同時に回すと GPU を取り合う（1 ホスト 1 ジョブ推奨）
#           GF_PKG_LOCAL（送る geneformer パッケージ。既定 ~/workspace/research/Geneformer/geneformer_hf/geneformer）
set -euo pipefail
MODE="${1:?usage: isp_remote.sh push|run|status|pull|clean HOST [CONFIG]}"
HOST="${2:?HOST}"
CONFIG="${3:-}"
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
RROOT="labretary-remote"                       # リモートの作業ルート（$HOME 相対）
GF_HOME="${REMOTE_GF_HOME:-\$HOME/workspace/Geneformer}"
GF_PKG_LOCAL="${GF_PKG_LOCAL:-$HOME/workspace/research/Geneformer/geneformer_hf/geneformer}"
PY="$ROOT/.venv/bin/python"
cd "$ROOT"

cfgval() { "$PY" -c "import yaml,sys;c=yaml.safe_load(open(sys.argv[1]));v=c
for k in sys.argv[2].split('.'): v=(v or {}).get(k)
print(v if v is not None else '')" "$CONFIG" "$1"; }

need_config() { [[ -n "$CONFIG" && -f "$CONFIG" ]] || { echo "CONFIG が見つかりません: $CONFIG"; exit 1; }; }

exp_dir() { dirname "$CONFIG"; }
results_dir() { local r; r="$(cfgval output.results_dir)"; echo "$(exp_dir)/${r:-results}"; }
tag() { basename "$CONFIG" .yaml | tr -c 'A-Za-z0-9_-' '_'; }

case "$MODE" in
  push)
    need_config
    TOK="$(cfgval data.tokenized_dataset)"; MODEL="$(cfgval model.directory)"
    [[ -d "$TOK" ]] || { echo "tokenized dataset がありません: $TOK"; exit 1; }
    ssh -o BatchMode=yes "$HOST" "mkdir -p $RROOT/geneformer_pkg $RROOT/models/Geneformer"
    # openrsync は -a だけだと無転送で成功することがある → -v で送り、最後にリモートで確認する（guardrails）
    # 小さいファイル（コード・config）は ssh 経由で書き、sha256 を照合する（既存ファイルの更新が rsync/scp で反映されないことがあった。guardrails）
    for f in scripts/geneformer/run_isp.py scripts/capture_env.sh "$CONFIG"; do
      ssh -o BatchMode=yes "$HOST" "mkdir -p $RROOT/$(dirname "$f") && cat > $RROOT/$f.tmp && mv $RROOT/$f.tmp $RROOT/$f" < "$f"
      want="$(shasum -a 256 "$f" | awk '{print $1}')"; got="$(ssh -o BatchMode=yes "$HOST" "sha256sum $RROOT/$f" | awk '{print $1}')"
      [[ "$want" == "$got" ]] && echo "OK: $f" || { echo "NG: $f の sha256 が一致しません"; exit 1; }
    done
    rsync -av --relative "$TOK" "$HOST:$RROOT/" | tail -2
    if [[ "$MODEL" == models/Geneformer/* ]]; then
      # 事前学習モデルはリモートの既存のものを使う（sha256 を照合）
      # REMOTE_PRETRAINED_DIR: 事前学習モデルの場所がホストで違う場合（coi-next1 は $GF_HOME/geneformer_hf/Geneformer-V2-316M）
      ssh -o BatchMode=yes "$HOST" "ln -sfn ${REMOTE_PRETRAINED_DIR:-$GF_HOME/$(basename "$MODEL")} $RROOT/$MODEL"
      want="$(shasum -a 256 "$MODEL/model.safetensors" | awk '{print $1}')"
      got="$(ssh -o BatchMode=yes "$HOST" "sha256sum $RROOT/$MODEL/model.safetensors" | awk '{print $1}')"
      [[ "$want" == "$got" ]] && echo "OK: $MODEL の重みは Mac と同一" || { echo "NG: $MODEL の重みが一致しません"; exit 1; }
    else
      rsync -av --relative "$MODEL" "$HOST:$RROOT/" | tail -2
    fi
    rsync -av --exclude '__pycache__' "$GF_PKG_LOCAL" "$HOST:$RROOT/geneformer_pkg/" | tail -2
    ssh -o BatchMode=yes "$HOST" "cd $RROOT && ls -d scripts/geneformer/run_isp.py $CONFIG $TOK $MODEL/config.json geneformer_pkg/geneformer/device.py"
    ;;
  run)
    need_config
    LOG="$(exp_dir)/logs/remote-$(tag).log"
    ssh -o BatchMode=yes "$HOST" "cd $RROOT && mkdir -p $(exp_dir)/logs && \
      tmux new -d -s isp-$(tag) 'export PYTHONPATH=\$HOME/$RROOT/geneformer_pkg LABRETARY_DEVICE=cuda TOKENIZERS_PARALLELISM=false; \
        PY=$GF_HOME/.venv/bin/python; \
        echo \"=== start \$(date +%F\ %T) \$(hostname)\" >> $LOG; \
        PY=\$PY bash scripts/capture_env.sh $(exp_dir) >> $LOG 2>&1; \
        \$PY scripts/geneformer/run_isp.py --config $CONFIG --approved >> $LOG 2>&1; \
        echo \"=== exit=\$? \$(date +%F\ %T)\" >> $LOG' && sleep 5 && tmux ls"
    ;;
  status)
    need_config
    ssh -o BatchMode=yes "$HOST" "tmux has-session -t isp-$(tag) 2>/dev/null && echo RUNNING || echo NOT_RUNNING; \
      grep -E '^\[|^===|Error|Traceback' $RROOT/$(exp_dir)/logs/remote-$(tag).log | tail -6"
    ;;
  pull)
    need_config
    R="$(results_dir)"
    mkdir -p "$R" "$(exp_dir)/logs"
    rsync -av "$HOST:$RROOT/$R/" "$R/" | tail -2
    rsync -av "$HOST:$RROOT/$(exp_dir)/logs/remote-$(tag).log" "$(exp_dir)/logs/" | tail -1
    rsync -av "$HOST:$RROOT/$(exp_dir)/results/env.json" "$R/env-remote.json" 2>/dev/null | tail -1 || true
    ls "$R" | head
    ;;
  clean)
    ssh -o BatchMode=yes "$HOST" "du -sh $RROOT 2>/dev/null; rm -rf -- \$HOME/$RROOT && echo removed"
    ;;
  *) echo "MODE は push | run | status | pull | clean"; exit 1 ;;
esac
