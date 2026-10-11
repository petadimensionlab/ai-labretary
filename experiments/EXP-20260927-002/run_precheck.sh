#!/usr/bin/env bash
# EXP-20260927-002 — HYP-20260925-02 の事前チェック（摂動なし）
# 実行するのは run_isp.py --precheck だけ。ISP は G1・G2 の後、別の plan で事前登録してから。
# 使い方: bash experiments/EXP-20260927-002/run_precheck.sh 2>&1 | tee logs/EXP-20260927-002-precheck.log
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"
EXP=experiments/EXP-20260927-002
CFG=$EXP/config.yaml
PY=.venv/bin/python

echo "[1/6] A1 の確認（生カウントか）"
$PY scripts/inspect_h5ad.py data/raw/kidney_tls/A1_KPMP_v1.5_snRNA_Lake2023.h5ad > $EXP/work/inspect_A1.json
$PY -c "import json,sys; r=list(json.load(open('$EXP/work/inspect_A1.json')).values())[0]; \
print(r['shape'], 'has_raw', r['has_raw'], 'counts_integer_like', r['counts_integer_like']); \
sys.exit(0 if r['counts_integer_like'] and r['has_raw'] else 1)"

echo "[2/6] prepare_h5ad（data/raw は読み取りのみ）"
$PY scripts/geneformer/prepare_h5ad.py --config $CFG

echo "[3/6] ラベル付与（tokenize の前に生カウントを再確認）"
$PY $EXP/analysis/label_fib_states.py --config $CFG

echo "[4/6] precheck（seed 0,1,2。MPS / bf16 / nproc=1 / emb_mode=cls）"
for S in 0 1 2; do
  $PY - "$CFG" "$S" <<'PYGEN'
import sys, yaml, pathlib
cfg = yaml.safe_load(open(sys.argv[1])); s = int(sys.argv[2])
cfg["precheck"]["seed"] = s
cfg["output"]["results_dir"] = f'{cfg["output"]["results_dir"]}/seed{s}'
out = pathlib.Path(cfg["output"]["experiment_dir"]) / "work" / f"config-seed{s}.yaml"
out.write_text(yaml.safe_dump(cfg, allow_unicode=True, sort_keys=False))
print(out)
PYGEN
  $PY scripts/geneformer/run_isp.py --config $EXP/work/config-seed$S.yaml --precheck
done

echo "[5/6] 事前登録した解析（深さマッチ・CI・並べ替え・感度分析）"
for S in 0 1 2; do $PY $EXP/analysis/precheck_ci.py --config $CFG --seed $S; done
$PY $EXP/analysis/precheck_ci.py --config $CFG --summary

echo "[6/6] 環境の記録（<EXP>/results/env.json に保存される）"
bash scripts/capture_env.sh $EXP
cp $EXP/results/env.json $EXP/results-precheck/env.json
ls -la $EXP/results-precheck
