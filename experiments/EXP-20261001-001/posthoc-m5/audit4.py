#!/usr/bin/env python3
"""検査そのものの点検のまとめ —
(1) 5 つの統計量の測定精度（対照 500 個を引き直したときの揺れ）
(2) 判断から導いた閾値の案と、その下できれいな指標が通る割合
(3) 性質どうしの重なり（移動量と元の順位は同じものか）
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
import diag  # noqa: E402

HERE = Path(__file__).parent
NULLSET = set(diag.NULLS_B500)


def main():
    cells = pd.read_parquet(diag.CACHE / "cells2.parquet")
    elig = diag.get_elig(cells)
    out = {}

    # ---- 性質どうしの重なり（対照 500 個の遺伝子レベル）
    ov = {}
    for s in diag.SEEDS:
        cs = cells[(cells.seed == s) & cells.gene.isin(NULLSET)]
        g = cs.groupby("gene").agg(ld=("dnorm", "median"), lp=("lp", "median"),
                                   ll=("ll", "median"), rp=("relpos", "median"),
                                   un=("m_uniform", "mean"))
        r = lambda a, b: round(float(stats.spearmanr(g[a], g[b]).statistic), 3)
        ov[s] = {"移動量と元の順位": r("ld", "lp"), "移動量と入力の長さ": r("ld", "ll"),
                 "移動量と相対位置": r("ld", "rp"), "移動量と一様成分": r("ld", "un"),
                 "一様成分と元の順位": r("un", "lp")}
    out["E_性質どうしの重なり（対照 500 個）"] = ov

    # ---- 測定精度（対照 500 個の引き直し）
    tabs = {}
    for s in diag.SEEDS:
        cs = cells[cells.seed == s]
        nl = {h: hd for h, hd in cs[cs.gene.isin(NULLSET)].groupby("gene")}
        f = diag.gene_stats(cs, nl, elig[s], "v_base", ["lp", "ll"])
        fu = diag.gene_stats(cs, nl, elig[s], "m_uniform", ["lp", "ll"])
        N = f[f.gene.isin(NULLSET)].copy()
        N["uni"] = fu.set_index("gene").reindex(N.gene).D_raw.to_numpy()
        tabs[s] = N
    rg = np.random.default_rng(99)
    keys = ("M5", "M4", "M2_lp", "M2_ll", "M2_relpos")
    B = {k: [] for k in keys}
    for _ in range(2000):
        ii = rg.integers(0, 500, 500)
        for k in keys:
            row = []
            for s in diag.SEEDS:
                N = tabs[s]
                y = N.D_raw.to_numpy()[ii]
                x = {"M5": N.dnorm_med, "M4": N.uni, "M2_lp": N.lp_med,
                     "M2_ll": N.ll_med, "M2_relpos": N.relpos_med}[k].to_numpy()[ii]
                row.append(float(stats.spearmanr(x, y).statistic))
            B[k].append(row)
    out["測定精度（対照 500 個を引き直したときの標準誤差）"] = {
        k: [round(float(x), 4) for x in np.array(v).std(axis=0, ddof=1)] for k, v in B.items()}
    out["実測値"] = {k: [round(float(stats.spearmanr(
        {"M5": tabs[s].dnorm_med, "M4": tabs[s].uni, "M2_lp": tabs[s].lp_med,
         "M2_ll": tabs[s].ll_med, "M2_relpos": tabs[s].relpos_med}[k],
        tabs[s].D_raw).statistic), 4) for s in diag.SEEDS] for k in keys}

    # ---- 判断から導いた閾値の案
    sl = json.loads((HERE / "audit_part3.json").read_text())[
        "統計量 1 単位あたりのセット平均 z の傾き"]
    dec = 0.41
    prop = {}
    for stat, cov in (("M5", "ld"), ("M4", "un"), ("M2_relpos", "rp")):
        s75 = sl[stat]["露出 p=0.75（実測 0.75）"] - sl[stat]["露出 p=0.5（実測 0.498）"]
        s60 = sl[stat]["露出 p=0.6（実測 0.598）"] - sl[stat]["露出 p=0.5（実測 0.498）"]
        prop[stat] = {
            "露出 p=0.75 を想定し歪みを判断の線の 25%（0.10 標準偏差）までに抑える閾値":
                round(float(0.25 * dec / abs(s75)), 3),
            "露出 p=0.6 を想定した同じ閾値": round(float(0.25 * dec / abs(s60)), 3),
            "いまの閾値": {"M5": 0.10, "M4": 0.40, "M2_relpos": 0.25}[stat]}
    out["C_判断から導いた閾値の案"] = prop

    # ---- きれいな指標（無作為な向き 60 本）が、案の閾値で通る割合
    ens = json.loads((HERE / "ensemble_random_dirs.json").read_text())
    def vals(key, sub=None):
        return np.array([[e[key][str(s)][sub] if sub else e[key][str(s)]
                          for s in diag.SEEDS] for e in ens])
    rates = {}
    for stat, arr in (("M5", vals("M5")), ("M4", vals("M4")),
                      ("M2_relpos", vals("M2", "relpos")), ("M2_ll", vals("M2", "ll")),
                      ("M2_lp", vals("M2", "lp"))):
        rates[stat] = {"いまの閾値での合格率": None, "統計量の標準偏差":
                       round(float(arr.std(ddof=1)), 3)}
        t = {"M5": 0.10, "M4": 0.40, "M2_relpos": 0.25, "M2_ll": 0.25, "M2_lp": 0.25}[stat]
        rates[stat]["いまの閾値での合格率"] = round(float((np.abs(arr) <= t).all(axis=1).mean()), 3)
        for cand in (0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5):
            rates[stat][f"閾値 {cand} での合格率"] = round(
                float((np.abs(arr) <= cand).all(axis=1).mean()), 3)
    out["B_きれいな指標（無作為な向き 60 本）の合格率"] = rates

    (HERE / "audit_part4.json").write_text(json.dumps(out, ensure_ascii=False, indent=1,
                                                      default=float) + "\n")
    print(json.dumps(out, ensure_ascii=False, indent=1, default=float))


if __name__ == "__main__":
    main()
