#!/usr/bin/env bash
# fp32 再計算（run 1 の境界遺伝子のみ）を seed 0 → 2 → 1 の順に直列実行（MPS は 1 ジョブずつ）
set -euo pipefail
cd "$(dirname "$0")/../../.."
E=experiments/EXP-20260926-001
for s in 0 2 1; do
  .venv/bin/python scripts/geneformer/run_isp.py --config $E/recheck-fp32/config-seed$s.yaml --run tls_high
done
