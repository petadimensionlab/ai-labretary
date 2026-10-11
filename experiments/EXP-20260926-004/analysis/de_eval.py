#!/usr/bin/env python3
"""EXP-20260926-004 (a): DE マッチ比較遺伝子の評価（plan.md §A.3 / §A.4 に事前登録）。CPU のみ。
1) 統合の妥当性: 再実行した VCAM1・GAPDH（全 seed）と seed 3 の null 1 遺伝子で、細胞・重心・shift・logits を EXP-003 と照合
2) 条件を満たせば EXP-003 results-ft に DE マッチ 10 遺伝子を加えて a2 / raw / dP_goal_a2 を計算し、A-P, A-S1〜S3 を出す
usage: .venv/bin/python experiments/EXP-20260926-004/analysis/de_eval.py
出力: experiments/EXP-20260926-004/de_matched/eval/{merge_check.json, gene_seed.csv, decision.json}
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as cm  # noqa: E402

BASE = cm.ROOT / "experiments/EXP-20260926-003/results-ft"
DERES = cm.ROOT / "experiments/EXP-20260926-004/de_matched/results-de"
OUT = cm.ROOT / "experiments/EXP-20260926-004/de_matched/eval"
RECHECK_NULL = ("run-3-L-1", "ENSG00000112297")
TOL = 1e-6
METRICS = {"m_a2": "a2（主）", "m_raw": "raw shift_goal", "dP_goal_a2": "ΔP(PT_healthy) 長さビン補正"}


def merge_check(base: pd.DataFrame, de: pd.DataFrame) -> dict:
    chk = {"genes": {}, "centroids": {}}
    pairs = [(r, g) for r in sorted(de.run.unique()) for g in ["ENSG00000162692", "ENSG00000111640"]] + [RECHECK_NULL]
    ok = True
    for r, g in pairs:
        a = base[(base.run == r) & (base.gene == g)].set_index("cell_uid")
        b = de[(de.run == r) & (de.gene == g)].set_index("cell_uid")
        same_cells = list(a.index) == list(b.index) or sorted(a.index) == sorted(b.index)
        d = {"n": [len(a), len(b)], "same_cells": bool(same_cells and len(a) == len(b))}
        if d["same_cells"]:
            b = b.loc[a.index]
            for col in ["shift_goal", "shift_start", "dP"]:
                d[f"max_abs_diff_{col}"] = float((a[col] - b[col]).abs().max())
            d["within_tol"] = d["max_abs_diff_shift_goal"] <= TOL
        else:
            d["within_tol"] = False
        ok &= d["within_tol"]
        chk["genes"][f"{cm.SYM.get(g, g)}@{r}"] = d
    for r in sorted(de.run.unique()):
        A = np.load(BASE / r / "state_emb_vectors.npz"); B = np.load(DERES / r / "state_emb_vectors.npz")
        chk["centroids"][r] = {k: float(np.abs(A[k].astype(np.float64) - B[k].astype(np.float64)).max()) for k in A.files}
    chk["merge_ok"] = bool(ok)
    return chk


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    base = cm.load_cells(BASE, logits=True)
    de = cm.load_cells(DERES, logits=True)
    chk = merge_check(base, de)
    (OUT / "merge_check.json").write_text(json.dumps(chk, ensure_ascii=False, indent=1) + "\n")
    print(json.dumps(chk, ensure_ascii=False, indent=1))
    if not chk["merge_ok"]:
        raise SystemExit("統合の条件（plan §A.4）を満たさないため評価しない")
    add = de[de.gene.isin(cm.DE)].copy()
    add["role"] = "de_matched"
    assert add.groupby("run").gene.nunique().min() == len(cm.DE), "DE 遺伝子が欠けている seed がある"
    merged = pd.concat([base, add], ignore_index=True)
    c = cm.add_metrics(merged)
    # 既存遺伝子の a2 は変わらない（補正は null 細胞だけで推定）ことを確認
    gsm = pd.read_csv(BASE / "eval/gene_seed_metrics.csv").set_index(["run", "gene"]).m_a2
    mine = c.groupby(["run", "gene"]).m_a2.mean()
    assert float((mine[gsm.index] - gsm).abs().max()) < 1e-12
    g = c.groupby(["run", "gene", "role"])[list(METRICS)].mean().reset_index()
    g["symbol"] = g.gene.map(cm.SYM)
    g["n_cells"] = c.groupby(["run", "gene", "role"]).size().values
    dec = {"merge_ok": True}
    rows = []
    for m in METRICS:
        per_seed = {}
        for run, s in g.groupby("run"):
            v = s.set_index("gene")[m]
            null = s[s.role == "negative"][m].to_numpy()
            dev = v[list(cm.DE)]
            pv, ph = v["ENSG00000162692"], v["ENSG00000113249"]
            band_v = dev[[k for k, b in cm.DE_BAND.items() if b == "VCAM1"]]
            band_h = dev[[k for k, b in cm.DE_BAND.items() if b == "HAVCR1"]]
            de_null_pct = {cm.SYM[k]: cm.null_pct(x, null) for k, x in dev.items()}
            per_seed[run] = {
                "VCAM1_pct_in_DE": float((dev < pv).mean()), "HAVCR1_pct_in_DE": float((dev < ph).mean()),
                "VCAM1_gt_all_band": bool((band_v < pv).all()), "HAVCR1_gt_all_band": bool((band_h < ph).all()),
                "VCAM1_null_pct": cm.null_pct(pv, null), "HAVCR1_null_pct": cm.null_pct(ph, null),
                "DE_null_pct": de_null_pct, "DE_mean_null_pct": float(np.mean(list(de_null_pct.values()))),
                "VCAM1_value": float(pv), "DE_values": {cm.SYM[k]: float(x) for k, x in dev.items()},
                "null_median": float(np.median(null))}
            for k, x in de_null_pct.items():
                rows.append({"metric": m, "run": run, "symbol": k, "band": cm.DE_BAND[[e for e, s_ in cm.SYM.items() if s_ == k][0]],
                             "null_pct": x})
        runs = sorted(per_seed)
        c1_de = {k: sum(per_seed[r]["DE_null_pct"][k] >= 0.90 for r in runs) for k in cm.DE.values()}
        mean_pct_gene = {k: float(np.mean([per_seed[r]["DE_null_pct"][k] for r in runs])) for k in cm.DE.values()}
        w = stats.wilcoxon(np.array(list(mean_pct_gene.values())) - 0.5, alternative="greater")
        n_vp = sum(per_seed[r]["VCAM1_pct_in_DE"] >= 0.90 for r in runs)
        dec[m] = {
            "VCAM1_pct_in_DE": [per_seed[r]["VCAM1_pct_in_DE"] for r in runs],
            "HAVCR1_pct_in_DE": [per_seed[r]["HAVCR1_pct_in_DE"] for r in runs],
            "VCAM1_n_seeds_pct_ge_0.90": int(n_vp),
            "A_P_verdict": ("VCAM1 は DE マッチより大きい（マーカー除去だけでは説明できない）" if n_vp >= 2
                            else "VCAM1 の C1 合格はマーカー除去と区別できない"),
            "A_S1_VCAM1_gt_all_band_seeds": int(sum(per_seed[r]["VCAM1_gt_all_band"] for r in runs)),
            "A_S1_HAVCR1_gt_all_band_seeds": int(sum(per_seed[r]["HAVCR1_gt_all_band"] for r in runs)),
            "A_S2_DE_pass_C1": {k: v for k, v in c1_de.items()}, "A_S2_n_DE_pass_C1": int(sum(v >= 2 for v in c1_de.values())),
            "A_S3_DE_mean_null_pct_per_seed": [round(per_seed[r]["DE_mean_null_pct"], 3) for r in runs],
            "A_S3_gene_mean_pct": {k: round(v, 3) for k, v in mean_pct_gene.items()},
            "A_S3_wilcoxon_p_greater_0.5": float(w.pvalue),
            "VCAM1_null_pct": [per_seed[r]["VCAM1_null_pct"] for r in runs],
            "HAVCR1_null_pct": [per_seed[r]["HAVCR1_null_pct"] for r in runs],
            "per_seed": per_seed}
    g.to_csv(OUT / "gene_seed.csv", index=False)
    pd.DataFrame(rows).to_csv(OUT / "de_null_pct.csv", index=False)
    (OUT / "decision.json").write_text(json.dumps(dec, ensure_ascii=False, indent=1, default=float) + "\n")
    for m, lab in METRICS.items():
        d = dec[m]
        print(f"== {lab}: VCAM1 pct in DE {d['VCAM1_pct_in_DE']} → {d['A_P_verdict']}")
        print(f"   HAVCR1 pct in DE {d['HAVCR1_pct_in_DE']}; band: VCAM1>{d['A_S1_VCAM1_gt_all_band_seeds']}/3, HAVCR1>{d['A_S1_HAVCR1_gt_all_band_seeds']}/3")
        print(f"   DE pass C1: {d['A_S2_n_DE_pass_C1']}/10 {d['A_S2_DE_pass_C1']}")
        print(f"   DE mean null pct per seed {d['A_S3_DE_mean_null_pct_per_seed']}, gene means {d['A_S3_gene_mean_pct']}, Wilcoxon p={d['A_S3_wilcoxon_p_greater_0.5']:.3g}")
        print(f"   VCAM1 null pct {d['VCAM1_null_pct']}, HAVCR1 {d['HAVCR1_null_pct']}")


if __name__ == "__main__":
    main()
