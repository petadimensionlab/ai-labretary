#!/usr/bin/env python3
"""
bf16 vs fp32 の matched pair（同じ細胞・同じ遺伝子・同じ重心の細胞で dtype だけ違う）を比較する。
判定基準は ../Geneformer/docs/quantization/EXPERIMENT.md §12:
  Spearman ≥ 0.98、上位 N の一致、符号一致は |shift| の中央値より上（またはノイズ床 3×median|Δ| より上）の遺伝子に限る。
usage: compare_precision.py --ref <fp32 results> --test <bf16 results> [--topn 5] [--out x.json]
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def load(res: Path):
    s = pd.read_csv(res / "summary.csv")
    c = pd.read_csv(res / "cells.csv", dtype={"donor_id": str}) if (res / "cells.csv").exists() else None
    t = [json.loads(x) for x in (res / "timing.jsonl").read_text().splitlines()] if (res / "timing.jsonl").exists() else []
    return s, c, t


def timing(t: list[dict]) -> dict:
    g = [x for x in t if "gene" in x and x.get("n_cells")]
    if not g:
        return {}
    cells = sum(x["n_cells"] for x in g)
    isp = sum(x["isp_seconds"] for x in g)
    load = sum(x["model_load_seconds"] for x in g)
    emb = sum(x["seconds"] for x in t if x.get("step") == "state_embs")
    return {"n_genes": len(g), "cell_gene_pairs": cells, "isp_seconds_total": round(isp, 1),
            "model_load_seconds_in_isp": round(load, 1), "state_embs_seconds": round(emb, 1),
            "sec_per_cell_gene_excl_load": round((isp - load) / cells, 4),
            "sec_per_cell_gene_incl_load": round(isp / cells, 4)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True, help="fp32 の results ディレクトリ")
    ap.add_argument("--test", required=True, help="bf16 の results ディレクトリ")
    ap.add_argument("--topn", type=int, default=5)
    ap.add_argument("--min-spearman", type=float, default=0.98)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    ra, ca, ta = load(Path(a.ref))
    rb, cb, tb = load(Path(a.test))
    col = "shift_goal_mean" if "shift_goal_mean" in ra.columns else "shift"
    m = ra[["gene", col]].merge(rb[["gene", col]], on="gene", suffixes=("_ref", "_test"))
    x, y = m[f"{col}_ref"].to_numpy(), m[f"{col}_test"].to_numpy()
    d = np.abs(x - y)
    med = np.median(np.abs(x))
    above = np.abs(x) > med
    floor = 3 * np.median(d)
    above_floor = np.abs(x) > floor
    topn = min(a.topn, len(m))
    top_ref = set(m.nlargest(topn, f"{col}_ref")["gene"])
    top_test = set(m.nlargest(topn, f"{col}_test")["gene"])
    res = {
        "ref": a.ref, "test": a.test, "metric": col, "n_genes": len(m),
        "spearman": round(float(stats.spearmanr(x, y).statistic), 4),
        "pearson": round(float(stats.pearsonr(x, y).statistic), 4),
        "topn": topn, "topn_overlap": len(top_ref & top_test),
        "sign_agreement_all": f"{int((np.sign(x) == np.sign(y)).sum())}/{len(m)}",
        "sign_agreement_above_median_abs_shift": f"{int((np.sign(x[above]) == np.sign(y[above])).sum())}/{int(above.sum())}",
        "noise_floor_3x_median_abs_diff": float(floor),
        "sign_agreement_above_noise_floor": f"{int((np.sign(x[above_floor]) == np.sign(y[above_floor])).sum())}/{int(above_floor.sum())}",
        "mean_abs_diff": float(d.mean()), "max_abs_diff": float(d.max()), "mean_abs_shift_ref": float(np.abs(x).mean()),
        "mean_abs_diff_over_mean_abs_shift": round(float(d.mean() / np.abs(x).mean()), 4),
        "genes": m.round(8).to_dict("records"),
    }
    if ca is not None and cb is not None:
        k = ["gene", "cell_uid"]
        cc = ca[k + ["shift_goal"]].merge(cb[k + ["shift_goal"]], on=k, suffixes=("_ref", "_test"))
        res["cell_level"] = {"n_pairs": len(cc),
                             "spearman": round(float(stats.spearmanr(cc.shift_goal_ref, cc.shift_goal_test).statistic), 4),
                             "pearson": round(float(stats.pearsonr(cc.shift_goal_ref, cc.shift_goal_test).statistic), 4),
                             "mean_abs_diff": float((cc.shift_goal_ref - cc.shift_goal_test).abs().mean()),
                             "mean_abs_shift_ref": float(cc.shift_goal_ref.abs().mean())}
    res["timing_ref"], res["timing_test"] = timing(ta), timing(tb)
    if res["timing_ref"] and res["timing_test"]:
        res["speedup_excl_load"] = round(res["timing_ref"]["sec_per_cell_gene_excl_load"]
                                         / res["timing_test"]["sec_per_cell_gene_excl_load"], 2)
    res["gate"] = {"spearman_ge_min": res["spearman"] >= a.min_spearman,
                   "topn_overlap_ge_n_minus_1": res["topn_overlap"] >= topn - 1}
    res["gate"]["pass"] = all(res["gate"].values())
    out = Path(a.out) if a.out else Path(a.test) / "precision_vs_ref.json"
    out.write_text(json.dumps(res, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({k: v for k, v in res.items() if k != "genes"}, indent=2, ensure_ascii=False))
    print(f"saved: {out}")


if __name__ == "__main__":
    main()
