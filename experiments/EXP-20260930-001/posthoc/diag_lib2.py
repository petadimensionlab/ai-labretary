#!/usr/bin/env python3
"""事後診断の追加部品（事前登録の外）: 移動量を**遺伝子ごとの代表値**で入れる版。

細胞ごとの移動量 ld を説明変数に入れると、ld の散らばりの 62% が「遺伝子の中」の
散らばりになり、遺伝子の中の傾きが遺伝子間の傾きに化けて効く（生態学的誤謬）。
遺伝子ごとの中央値だけを使えば、補正は遺伝子間の関係だけを使う。
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def _d(lp, ll, ld, form: str) -> np.ndarray:
    one = np.ones_like(lp)
    base = [one, lp, ll, lp ** 2, ll ** 2, lp * ll]
    if form == "none":
        return np.column_stack(base)
    if form == "ld1":
        return np.column_stack(base + [ld])
    if form == "ld2":
        return np.column_stack(base + [ld, ld ** 2])
    if form == "full":
        return np.column_stack(base + [ld, ld ** 2, lp * ld, ll * ld])
    raise ValueError(form)


def gene_stats_genemed(cs: pd.DataFrame, nulls: dict, elig: dict, col: str, min_elig: int,
                       form: str) -> pd.DataFrame:
    """移動量は**その遺伝子の中央値**を細胞すべてに当てる（遺伝子の中の散らばりを使わない）。"""
    med = cs.groupby("gene").dnorm.median()

    def X(df, gene):
        n = len(df)
        return _d(df.lp.to_numpy(float), df.ll.to_numpy(float),
                  np.full(n, np.log2(med[gene])), form)

    hs = list(nulls)
    Xn = {h: X(hd, h) for h, hd in nulls.items()}
    yn = {h: nulls[h][col].to_numpy() for h in hs}
    G = {h: Xn[h].T @ Xn[h] for h in hs}
    B = {h: Xn[h].T @ yn[h] for h in hs}
    Xbar = np.vstack([Xn[h].mean(axis=0) for h in hs])
    ybar = np.array([yn[h].mean() for h in hs])
    hi = {h: i for i, h in enumerate(hs)}
    G_all, B_all = sum(G.values()), sum(B.values())
    rows = []
    for gene, tg in cs.groupby("gene"):
        Gm, Bm = (G_all - G[gene], B_all - B[gene]) if gene in G else (G_all, B_all)
        beta = np.linalg.lstsq(Gm, Bm, rcond=None)[0]
        Xg = X(tg, gene)
        Dg = float(tg[col].to_numpy().mean() - Xg.mean(axis=0) @ beta)
        e = elig.get(gene, [])
        if e:
            ii = np.array([hi[h] for h in e])
            Dh = ybar[ii] - Xbar[ii] @ beta
        else:
            Dh = np.zeros(0)
        q = float((Dh < Dg).mean()) if len(Dh) >= min_elig else np.nan
        rows.append({"gene": gene, "D_raw": Dg, "q": q, "n_elig": len(Dh),
                     "dnorm_med": float(tg.dnorm.median())})
    return pd.DataFrame(rows)
