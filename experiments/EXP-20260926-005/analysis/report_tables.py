"""EXP-20260926-005: 報告用の表（事前登録の判定には使わない。判定は evaluate_positives.py の criteria.json）。

evaluate_positives.py と同じ読み込み・指標計算を再利用し、次を出す:
  - per_seed_pct.csv   : 陽性・HK・パネル・同梱対照の seed ごとの null 分位（a2 / raw / dP_a2 / dP_raw）と値
  - reference_C4_C5.json: 参考基準 C4（donor split-half 20 回・seed 間の固定遺伝子順位）と C5（長さとの相関）
  - boundary.json      : a2 で分位 0.85〜0.95 に入った陽性（plan.md §7 の fp32 確認対象）
盲検: evaluate_positives.load_run が読み込み直後に rc.BLIND を除外（assert）。
実行: .venv/bin/python experiments/EXP-20260926-005/analysis/report_tables.py --config experiments/EXP-20260926-005/config-delete.yaml
"""
from __future__ import annotations

import argparse
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import evaluate_positives as ep  # noqa: E402

rc, vp = ep.rc, ep.vp
METRICS = ["m_a2", "m_raw", "dP_a2", "dP_raw", "dM_raw"]  # dM_raw: 探索的（事前登録外）。logit 差 (l1-l0) の変化


def logit_margin_delta(rd: Path, cells: pd.DataFrame) -> np.ndarray:
    """探索的指標（事前登録外）: 細胞ごとの Δ[(logit_PT_healthy − logit_aPT)]。分類器が飽和（P≈0）して ΔP が 0 に潰れるため併記する。"""
    v = np.full(len(cells), np.nan)
    for gene, s in cells.groupby("gene", sort=False):
        lo = np.load(rd / "genes" / gene / "emb" / "logits_orig.npy").astype(np.float64)
        lp = np.load(rd / "genes" / gene / "emb" / "logits_pert.npy").astype(np.float64)
        per = pd.read_csv(rd / "genes" / gene / "cells.csv", dtype={"donor_id": str})
        m = dict(zip(per.cell_uid, (lp[:, 1] - lp[:, 0]) - (lo[:, 1] - lo[:, 0])))
        v[s.index.to_numpy()] = s.cell_uid.map(m).to_numpy()
    return v


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    a = ap.parse_args()
    cfg = yaml.safe_load(open(a.config))
    res = ep.E / cfg["output"]["results_dir"]
    out = res / "eval"
    out.mkdir(parents=True, exist_ok=True)
    from geneformer import TOKEN_DICTIONARY_FILE
    with open(Path(TOKEN_DICTIONARY_FILE).with_name("gene_name_id_dict_gc104M.pkl"), "rb") as f:
        names = {v: k for k, v in pickle.load(f).items()}

    parts = []
    for run in sorted(pd.read_csv(res / "cells.csv", usecols=["run"]).run.unique()):
        cells, E0, E1, dP, cen = ep.load_run(res, run, True)
        cells, _ = vp.vector_metrics(cells, E0, E1, cen)
        alt = [x for x in cells.columns if x.startswith("shift_alt")][0]
        cells["alt"] = cells[alt]
        cells["pos"] = ep.token_positions(cells)
        cells = rc.corrected(cells)
        cells["dP_raw"] = dP
        ep.bin_correct(cells, "dP_raw", "dP_a2", cfg["analysis"]["length_bins"])
        cells["dM_raw"] = logit_margin_delta(res / run, cells)
        parts.append(cells)
    cells = pd.concat(parts, ignore_index=True)
    runs = sorted(cells.run.unique())

    g = cells.groupby(["run", "gene", "role"]).agg(
        **{m: (m, "mean") for m in METRICS}, L=("length", "mean"), n=("cell_uid", "size")).reset_index()
    g["gene_name"] = g.gene.map(names)
    rows = []
    for run, s in g.groupby("run"):
        null = s[s.role == "negative"]
        for _, r in s[s.role != "negative"].iterrows():
            d = {"seed": run, "gene": r.gene_name, "role": r.role, "n_cells": int(r.n)}
            for m in METRICS:
                d[f"{m}_pct"] = round(float(rc.null_pct(r[m], null[m].to_numpy())), 3)
                d[m] = float(r[m])
            rows.append(d)
        rows.append({"seed": run, "gene": "(null 20 mean)", "role": "negative", "n_cells": int(null.n.sum()),
                     **{m: float(null[m].mean()) for m in METRICS},
                     **{f"{m}_null_sd": float(null[m].std()) for m in METRICS}})
    pt = pd.DataFrame(rows)
    pt.to_csv(out / "per_seed_pct.csv", index=False)

    # 参考 C4 / C5（rc.evaluate と同じ計算。陽性リストだけ config から）
    fixed = [x for x in cfg["controls"]["positive"]] + list(cfg["controls"]["negative_fixed"]) + \
        [x for grp in cfg["controls"].get("extra_fixed", []) or [] for x in grp["genes"]]
    ref = {}
    rng = np.random.default_rng(0)
    for m in ["m_a2", "dP_a2"]:
        rhos = []
        for run, s in cells.groupby("run"):
            donors = np.array(sorted(s.donor_id.unique()))
            for _ in range(cfg["analysis"].get("split_half_n", 20)):
                half = set(rng.permutation(donors)[: len(donors) // 2])
                h = s.donor_id.isin(half)
                x = s[h].groupby("gene")[m].mean(); y = s[~h].groupby("gene")[m].mean()
                idx = x.index.intersection(y.index)
                rhos.append(stats.spearmanr(x[idx], y[idx]).statistic)
        piv = g[g.gene.isin(fixed)].pivot(index="gene", columns="run", values=m)
        rr = [stats.spearmanr(piv[p], piv[q]).statistic for i, p in enumerate(runs) for q in runs[i + 1:]]
        rcell = [stats.spearmanr(s.length, s[m]).statistic for _, s in cells.groupby(["run", "gene"])
                 if s.length.nunique() > 5]
        ref[m] = {"C4i_splithalf_rho_mean": round(float(np.nanmean(rhos)), 3),
                  "C4ii_seed_rho_fixed_genes": round(float(np.mean(rr)), 3), "n_fixed_genes": len(fixed),
                  "C4ii_pairs": [round(float(x), 3) for x in rr],
                  "C5_rho_len_gene": round(float(stats.spearmanr(g.L, g[m]).statistic), 3),
                  "C5_rho_ncells_gene": round(float(stats.spearmanr(g.n, g[m]).statistic), 3),
                  "C5_rho_len_cell_median": round(float(np.nanmedian(rcell)), 3)}
    (out / "reference_C4_C5.json").write_text(json.dumps(ref, indent=1) + "\n")

    pos_names = {names[x] for x in cfg["controls"]["positive"]}
    b = pt[pt.gene.isin(pos_names) & pt.m_a2_pct.between(0.85, 0.95)]
    (out / "boundary.json").write_text(json.dumps(b[["seed", "gene", "m_a2_pct"]].to_dict("records"), indent=1) + "\n")

    pd.set_option("display.width", 250); pd.set_option("display.max_rows", 200)
    show = pt[pt.role != "negative"].pivot_table(index=["role", "gene"], columns="seed",
                                                 values=["m_a2_pct", "dP_a2_pct"]).round(2)
    print(show.to_string())
    print(json.dumps(ref, indent=1))
    print("boundary:", b[["seed", "gene", "m_a2_pct"]].to_dict("records"))


if __name__ == "__main__":
    main()
