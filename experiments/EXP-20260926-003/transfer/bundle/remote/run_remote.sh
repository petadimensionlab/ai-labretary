#!/usr/bin/env bash
# EXP-20260926-003: リモート NVIDIA GPU での fine-tuning（bundle/remote/ にコピーされて使われる）
# 前提: リモートに Geneformer がセットアップ済み（既定 $HOME/workspace/Geneformer。venv と Geneformer-V2-316M を含む）
#
#   bash run_remote.sh check            # 環境の確認だけ（torch / CUDA / bf16 / geneformer / モデル / bundle の sha256）
#   bash run_remote.sh probe  [PROFILE] # 20 optimizer step だけ回して s/step とピークメモリを測る
#   bash run_remote.sh train  [PROFILE] # 本番（中断しても同じコマンドで checkpoint から再開）
#
# PROFILE: vram16 | vram24 | vram48 | vram80（既定 vram80 = DGX Spark / GB10 128 GB 統合メモリ向け。実効バッチ 16 は共通）
# 環境変数（必要なときだけ上書き。ホスト名や資格情報はここに書かない）:
#   GF_HOME   リモートの Geneformer（既定: $HOME/workspace/Geneformer。spark-cf5b の配置）
#   GF_VENV   venv（既定: $GF_HOME/.venv。なければ PATH の python3）
#   GF_BASE_MODEL_DIR  Geneformer-V2-316M（既定: $GF_HOME/Geneformer-V2-316M、なければ $GF_HOME/geneformer_hf/Geneformer-V2-316M）
#   GF_PKG    geneformer パッケージの親ディレクトリ（既定: bundle の隣の geneformer_pkg。あれば PYTHONPATH に入れる。venv は変更しない）
#   GF_OUT    出力（既定: bundle の隣の out/ft-main）
set -euo pipefail
MODE="${1:?usage: run_remote.sh check|probe|train [PROFILE]}"
PROFILE="${2:-vram80}"
HERE="$(cd "$(dirname "$0")" && pwd)"
BUNDLE="$(cd "$HERE/.." && pwd)"
GF_HOME="${GF_HOME:-$HOME/workspace/Geneformer}"
if [[ -n "${GF_VENV:-}" ]]; then PY="$GF_VENV/bin/python"
elif [[ -x "$GF_HOME/.venv/bin/python" ]]; then PY="$GF_HOME/.venv/bin/python"
else PY="$(command -v python3)"; fi
if [[ -z "${GF_BASE_MODEL_DIR:-}" ]]; then
  if [[ -f "$GF_HOME/Geneformer-V2-316M/config.json" ]]; then GF_BASE_MODEL_DIR="$GF_HOME/Geneformer-V2-316M"
  else GF_BASE_MODEL_DIR="$GF_HOME/geneformer_hf/Geneformer-V2-316M"; fi
fi
export GF_BASE_MODEL_DIR
GF_PKG="${GF_PKG:-$(dirname "$BUNDLE")/geneformer_pkg}"
[[ -d "$GF_PKG/geneformer" ]] && export PYTHONPATH="$GF_PKG${PYTHONPATH:+:$PYTHONPATH}"
OUT="${GF_OUT:-$(dirname "$BUNDLE")/out/ft-main}"
export CUDA_VISIBLE_DEVICES="${GPU_ID:-0}"
export GENEFORMER_DEVICE=cuda
export TOKENIZERS_PARALLELISM=false

[[ -x "$PY" ]] || { echo "NG: python が見つかりません（GF_VENV を指定）"; exit 1; }
[[ -f "$GF_BASE_MODEL_DIR/config.json" ]] || { echo "NG: モデルが見つかりません: $GF_BASE_MODEL_DIR（GF_BASE_MODEL_DIR を指定）"; exit 1; }

echo "== 環境 =="
echo "python: $PY"
"$PY" - <<'EOF'
import torch, transformers, datasets, geneformer, os
print("torch", torch.__version__, "cuda", torch.version.cuda, "available", torch.cuda.is_available())
if torch.cuda.is_available():
    print("gpu", torch.cuda.get_device_name(0), "bf16", torch.cuda.is_bf16_supported())
print("transformers", transformers.__version__, "datasets", datasets.__version__)
d = os.path.dirname(geneformer.__file__)
print("geneformer", d, "| device パッチ:", os.path.exists(os.path.join(d, "device.py")))
EOF
nvidia-smi --query-gpu=index,name,driver_version,memory.total --format=csv,noheader || true

echo "== bundle の sha256 照合 =="
"$PY" "$HERE/verify_sha256.py" check "$BUNDLE"
echo "== モデル重みの sha256 照合 =="
want="$(awk '{print $1}' "$HERE/MODEL_SHA256")"
got="$("$PY" -c "import hashlib,sys;h=hashlib.sha256();f=open(sys.argv[1],'rb');[h.update(b) for b in iter(lambda:f.read(1<<20),b'')];print(h.hexdigest())" "$GF_BASE_MODEL_DIR/model.safetensors")"
if [[ "$want" == "$got" ]]; then echo "OK: model.safetensors（Mac と同一）"
else echo "NG: model.safetensors が Mac の重みと一致しません"; exit 1; fi

ARGS=(--config "$BUNDLE/config-finetune.yaml" --stage train --profile "$PROFILE" --bundle "$BUNDLE" --out "$OUT")
mkdir -p "$(dirname "$OUT")"
case "$MODE" in
  check) "$PY" "$BUNDLE/finetune_classifier.py" "${ARGS[@]}" --dry-run ;;
  probe) "$PY" "$BUNDLE/finetune_classifier.py" "${ARGS[@]}" --probe 20 2>&1 | tee -a "$(dirname "$OUT")/probe.log" ;;
  train)
    "$PY" "$BUNDLE/finetune_classifier.py" "${ARGS[@]}" 2>&1 | tee -a "$(dirname "$OUT")/train.log"
    cp "$(dirname "$OUT")/train.log" "$OUT/train.log"
    "$PY" "$HERE/verify_sha256.py" check "$OUT/final"
    echo "完了: $OUT"
    ;;
  *) echo "MODE は check | probe | train"; exit 1 ;;
esac
