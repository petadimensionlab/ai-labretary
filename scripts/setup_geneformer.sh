#!/usr/bin/env bash
# 既存の Geneformer インストールをこのプロジェクトに接続する（clone / pip install はしない）
# GENEFORMER_HOME: Geneformer ワークスペース（.venv と geneformer_hf/ を含む）。既定は ../Geneformer
set -euo pipefail
cd "$(dirname "$0")/.."
GENEFORMER_HOME="${GENEFORMER_HOME:-$(cd .. && pwd)/Geneformer}"
[[ -d "$GENEFORMER_HOME/geneformer_hf" ]] || { echo "Geneformer が見つかりません: $GENEFORMER_HOME（GENEFORMER_HOME で指定してください）"; exit 1; }
[[ -x "$GENEFORMER_HOME/.venv/bin/python" ]] || { echo "venv が見つかりません: $GENEFORMER_HOME/.venv"; exit 1; }

# 既存の実体（ディレクトリ・ファイル）は上書きしない。シンボリックリンクのみ作成／更新
link() {
  local target="$1" name="$2"
  if [[ -e "$name" && ! -L "$name" ]]; then
    echo "$name は既に存在します（リンクではない）。スキップ"
  else
    ln -sfn "$target" "$name" && echo "$name -> $target"
  fi
}
link "$GENEFORMER_HOME/.venv" .venv
link "$GENEFORMER_HOME/geneformer_hf" models/Geneformer

.venv/bin/python - <<'PY'
import inspect, geneformer
from geneformer import InSilicoPerturber, InSilicoPerturberStats, EmbExtractor, TranscriptomeTokenizer
import yaml, pandas, scipy, matplotlib, scanpy, anndata
print("geneformer OK:", geneformer.__file__)
for c in (TranscriptomeTokenizer, EmbExtractor, InSilicoPerturber, InSilicoPerturberStats):
    print(c.__name__, list(inspect.signature(c.__init__).parameters)[1:])
PY
echo "利用可能なモデル（config.yaml の model.directory に指定）:"
find -L models/Geneformer -maxdepth 2 -name config.json -exec dirname {} \;
