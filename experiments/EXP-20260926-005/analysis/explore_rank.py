"""探索的（事前登録外）: overexpress/delete の a2 と、摂動遺伝子の元の順位（input_ids 内の相対位置 pos/length）の関係。
全遺伝子（null 含む）で seed ごとに遺伝子単位 Spearman を取る。出力: <results>/eval/explore_rank.json
実行: .venv/bin/python experiments/EXP-20260926-005/analysis/explore_rank.py --config experiments/EXP-20260926-005/config-overexpress.yaml
"""
import argparse, json, sys
from pathlib import Path
import numpy as np, pandas as pd, yaml
from scipy import stats
sys.path.insert(0, str(Path(__file__).resolve().parent))
import evaluate_positives as ep  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--config", required=True); a = ap.parse_args()
cfg = yaml.safe_load(open(a.config)); res = ep.E / cfg["output"]["results_dir"]
out = {}
for run in sorted(pd.read_csv(res / "cells.csv", usecols=["run"]).run.unique()):
    cells, E0, E1, dP, cen = ep.load_run(res, run, False)
    cells, _ = ep.vp.vector_metrics(cells, E0, E1, cen)
    cells["alt"] = cells[[x for x in cells.columns if x.startswith("shift_alt")][0]]
    cells["pos"] = ep.token_positions(cells)
    cells = ep.rc.corrected(cells)
    cells["relpos"] = cells.pos / cells.length
    g = cells.groupby(["gene", "role"]).agg(m=("m_a2", "mean"), rp=("relpos", "mean")).reset_index()
    nul = g[g.role == "negative"]
    r_all = stats.spearmanr(g.rp, g.m); r_null = stats.spearmanr(nul.rp, nul.m)
    rc_cell = [stats.spearmanr(s.relpos, s.m_a2).statistic for _, s in cells.groupby("gene")]
    out[run] = {"n_genes": len(g), "rho_gene_all": round(float(r_all.statistic), 3), "p_all": float(r_all.pvalue),
                "rho_gene_null20": round(float(r_null.statistic), 3), "p_null20": float(r_null.pvalue),
                "rho_cell_within_gene_median": round(float(np.nanmedian(rc_cell)), 3)}
(res / "eval/explore_rank.json").write_text(json.dumps(out, indent=1) + "\n")
print(json.dumps(out, indent=1))
