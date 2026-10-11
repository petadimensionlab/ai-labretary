#!/usr/bin/env python3
"""事後診断の共通部品（事前登録の外）。

`scripts/analysis/movement_correction.py` の `gene_stats` と同じ計算を、
**設計行列だけ差し替えられる**ようにしたもの（元のモジュールは書き換えない）。
leave-self-out の正規方程式で解く手順は元の実装と同一。
"""
from __future__ import annotations

import numpy as np
import pandas as pd


# ---- 設計行列の版 ------------------------------------------------------------
def _d(lp, ll, ld, form: str) -> np.ndarray:
    one = np.ones_like(lp)
    base = [one, lp, ll, lp ** 2, ll ** 2, lp * ll]          # 従来（lp・ll の完全二次）
    if form == "none":                 # 移動量を入れない（EXP-20260927-004 以降の従来版）
        return np.column_stack(base)
    if form == "ld1":                  # ld の 1 次だけ足す
        return np.column_stack(base + [ld])
    if form == "ld2":                  # ld と ld² を足す
        return np.column_stack(base + [ld, ld ** 2])
    if form == "ld_cross":             # ld と交差項（ld² なし）
        return np.column_stack(base + [ld, lp * ld, ll * ld])
    if form == "full":                 # 事前登録の完全二次式（ld・ld²・lp·ld・ll·ld）
        return np.column_stack(base + [ld, ld ** 2, lp * ld, ll * ld])
    if form == "linear3":              # lp・ll・ld すべて 1 次だけ
        return np.column_stack([one, lp, ll, ld])
    if form == "lplll1_ld1":           # lp・ll 1 次 + ld 1 次
        return np.column_stack([one, lp, ll, ld])
    if form == "intercept":            # 補正を一切しない（定数のみ）
        return np.column_stack([one])
    raise ValueError(form)


FORMS = ("intercept", "none", "ld1", "ld2", "ld_cross", "full", "linear3")

FORM_LABEL = {
    "intercept": "補正なし（定数のみ）",
    "none": "従来（元の位置・入力の長さの完全二次、移動量を入れない）",
    "ld1": "従来 ＋ 移動量の対数の 1 次（1 項追加）",
    "ld2": "従来 ＋ 移動量の対数の 1 次と 2 乗（2 項追加）",
    "ld_cross": "従来 ＋ 移動量の対数の 1 次と交差項（2 乗なし、3 項追加）",
    "full": "事前登録の完全二次式（移動量の 1 次・2 乗・交差項 2 本、4 項追加）",
    "linear3": "元の位置・入力の長さ・移動量をすべて 1 次だけ",
}


def gene_stats(cs: pd.DataFrame, nulls: dict, elig: dict, col: str, min_elig: int,
               form: str) -> pd.DataFrame:
    """遺伝子ごとに、対照遺伝子（自分自身は除く）で当てた回帰からの残差平均 D_raw と
    対照の中での順位 q を返す。`movement_correction.gene_stats` と同じ手順。"""
    def X(df):
        return _d(df.lp.to_numpy(float), df.ll.to_numpy(float),
                  np.log2(df.dnorm.to_numpy(float)), form)

    hs = list(nulls)
    Xn = {h: X(hd) for h, hd in nulls.items()}
    yn = {h: nulls[h][col].to_numpy() for h in hs}
    G = {h: Xn[h].T @ Xn[h] for h in hs}
    B = {h: Xn[h].T @ yn[h] for h in hs}
    Xbar = np.vstack([Xn[h].mean(axis=0) for h in hs])
    ybar = np.array([yn[h].mean() for h in hs])
    hi = {h: i for i, h in enumerate(hs)}
    G_all = sum(G.values())
    B_all = sum(B.values())
    rows = []
    for gene, tg in cs.groupby("gene"):
        Gm, Bm = (G_all - G[gene], B_all - B[gene]) if gene in G else (G_all, B_all)
        beta = np.linalg.solve(Gm, Bm)
        Xg = X(tg)
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


def metric_columns(c: pd.DataFrame) -> None:
    """事前登録の `evaluate_004.metric_columns` と同じ主指標・一様成分を付ける。"""
    dn = c["dnorm"]
    c["m_goal"] = (c.s_goal - c.s_start) / dn           # 主指標
    c["m_uniform"] = (c.s_goal + c.s_start) / 2.0       # 一様成分（M4）
    c["m_A11"] = c.s_goal - c.s_start                   # 移動量で割らない版
    c["m_uniform_norm"] = (c.s_goal + c.s_start) / 2.0 / dn   # 一様成分も移動量で割った版（診断）
