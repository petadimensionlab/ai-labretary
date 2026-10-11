#!/usr/bin/env python3
"""M5 の事後診断 — 指標の候補を並べて M1〜M5 と目標の軸の陽性対照 71 個を測る。

**主判定の値は計算しない**（介入セット 119 個とその入れ子の遺伝子は控えに入っていない。
build_cache.py が読み込み直後に捨てている）。

閾値は一切動かさない（M1 0.2 / M2 0.25 / M3 [0.1,0.9] / M4 0.40 / M5 0.10 / 陽性対照 0.53）。
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
import diag  # noqa: E402

HERE = Path(__file__).parent


# --------------------------------------------------------------- 共通の向き
def build_common(cells: pd.DataFrame, disp: np.ndarray) -> dict:
    """対照 500 個から「摂動すると誰もが動く向き」ĉ を seed ごとに作る。
    遺伝子ごとの寄与を控えておき、その遺伝子を除いた ĉ を使えるようにする。"""
    info = {}
    for s in diag.SEEDS:
        m = ((cells.seed == s) & cells.gene.isin(diag.NULLS_B500)).to_numpy()
        idx = np.flatnonzero(m)
        gl = cells.gene.to_numpy()[idx]
        per, tot, n = {}, np.zeros(disp.shape[1]), 0
        for g in np.unique(gl):
            ii = idx[gl == g]
            d = disp[ii].astype(np.float64)
            d = d / np.linalg.norm(d, axis=1, keepdims=True)
            per[g] = (d.sum(axis=0), len(ii))
            tot += per[g][0]
            n += len(ii)
        info[s] = {"per": per, "tot": tot, "n": n}
    return info


def loo_c(info: dict, gene: str) -> np.ndarray:
    tot, n = info["tot"], info["n"]
    if gene in info["per"]:
        sv, k = info["per"][gene]
        tot, n = tot - sv, n - k
    v = tot / n
    return v / np.linalg.norm(v)


def add_columns(cells: pd.DataFrame, disp: np.ndarray, states: dict) -> pd.DataFrame:
    c = cells
    c["m_uniform"] = c[["s_dFIB", "s_dMFIB"]].mean(axis=1)
    c["ld"] = np.log2(c.dnorm)
    info = build_common(c, disp)
    U = {}
    for s in diag.SEEDS:
        st = {k: np.asarray(v) for k, v in states[f"run-{s}-L-1"].items()}
        vg, vs = st["goal"], st["start"]
        u = vg / np.linalg.norm(vg) - vs / np.linalg.norm(vs)
        U[s] = u / np.linalg.norm(u)
    n = len(c)
    dc = np.full(n, np.nan)          # d · ĉ
    dup = np.full(n, np.nan)         # d · û⊥（û から ĉ 成分を抜いた向き）
    dnp = np.full(n, np.nan)         # ‖d − (d·ĉ)ĉ‖
    uc = {}
    for (s, g), part in c.groupby(["seed", "gene"], sort=False):
        ch = loo_c(info[s], g)
        u = U[s]
        uo = u - (u @ ch) * ch
        uo = uo / np.linalg.norm(uo)
        uc[(s, g)] = float(u @ ch)
        ii = part.index.to_numpy()
        d = disp[ii].astype(np.float64)
        a = d @ ch
        dc[ii] = a
        dup[ii] = d @ uo
        dnp[ii] = np.sqrt(np.maximum(np.einsum("ij,ij->i", d, d) - a * a, 0))
    c["dc"] = dc
    c["duperp"] = dup
    c["dnperp"] = dnp
    c["f_common"] = c.dc / c.dnorm          # 変位のうち共通の向きが占める割合
    # ---- 指標の候補（細胞ごと）
    c["v_base"] = c.du / c.dnorm                                  # 事前登録
    c["v_uperp"] = c.duperp / c.dnorm                             # 向きから共通成分を抜く
    c["v_dperp"] = c.duperp / c.dnperp                            # 変位からも共通成分を抜く
    c["v_nodiv"] = c.du                                           # 割らない
    c["v_div_o"] = c.du / c.onorm                                 # 元の埋め込みの長さで割る
    c["v_div_perp"] = c.du / c.dnperp                             # 共通成分を除いた長さで割る
    (HERE / "u_dot_c.json").write_text(json.dumps(
        {str(s): round(float(np.mean([v for (ss, _), v in uc.items() if ss == s])), 5)
         for s in diag.SEEDS}, ensure_ascii=False))
    return c


# --------------------------------------------------------------- 遺伝子単位の後処理
def gene_level_residual(f: pd.DataFrame, nullset: set, mode: str) -> pd.DataFrame:
    """遺伝子単位での追加補正。mode:
       'none'     なにもしない
       'lpll'     遺伝子ごとの中央値 lp・ll の完全二次で、対照だけを使って当てた回帰の残差
       'ld'       さらに log2(移動量の中央値) の 2 次を足す
       'ld_rank'  log2(移動量の中央値) の順位を 10 分位に分け、対照の分位平均を引く（単調でなくてよい）
    """
    f = f.copy()
    if mode == "none":
        return f
    nul = f[f.gene.isin(nullset)]
    lp, ll, ld = f.lp_med.to_numpy(), f.ll_med.to_numpy(), np.log2(f.dnorm_med.to_numpy())
    if mode == "ld_rank":
        q = np.log2(nul.dnorm_med.to_numpy())
        edges = np.quantile(q, np.linspace(0, 1, 11))
        edges[0], edges[-1] = -np.inf, np.inf
        b_all = np.digitize(ld, edges[1:-1])
        b_nul = np.digitize(q, edges[1:-1])
        adj = np.zeros(len(f))
        for b in range(10):
            mu = nul.D_raw.to_numpy()[b_nul == b].mean()
            adj[b_all == b] = mu
        f["D_raw"] = f.D_raw - adj
        return f
    cols = [lp, ll] + ([ld] if mode == "ld" else [])
    X = diag.design(cols, "full_quadratic")
    mN = f.gene.isin(nullset).to_numpy()
    beta, *_ = np.linalg.lstsq(X[mN], f.D_raw.to_numpy()[mN], rcond=None)
    f["D_raw"] = f.D_raw - X @ beta
    return f


# --------------------------------------------------------------- 候補の定義
VARIANTS = [
    # (名前, 細胞ごとの指標の列, 回帰の説明変数, 回帰の形, 集計, 遺伝子単位の後処理)
    ("0 事前登録（基準）",                 "v_base",   ["lp", "ll"], "full_quadratic", "mean", "none"),
    ("1 向きから共通成分を抜く",           "v_uperp",  ["lp", "ll"], "full_quadratic", "mean", "none"),
    ("2 変位からも共通成分を抜く",         "v_dperp",  ["lp", "ll"], "full_quadratic", "mean", "none"),
    ("3 移動量で割らない",                 "v_nodiv",  ["lp", "ll"], "full_quadratic", "mean", "none"),
    ("4 元の埋め込みの長さで割る",         "v_div_o",  ["lp", "ll"], "full_quadratic", "mean", "none"),
    ("5 共通成分を除いた長さで割る",       "v_div_perp", ["lp", "ll"], "full_quadratic", "mean", "none"),
    ("6 細胞ごとの回帰に移動量を足す",     "v_base",   ["lp", "ll", "ld"], "full_quadratic", "mean", "none"),
    ("7 残差を中央値で集計",               "v_base",   ["lp", "ll"], "full_quadratic", "median", "none"),
    ("8 遺伝子単位でやり直す（位置と長さ）", "v_base",  [],           "intercept",      "mean", "lpll"),
    ("9 遺伝子単位で移動量も補正（二次）",  "v_base",   [],           "intercept",      "mean", "ld"),
    ("10 遺伝子単位で移動量も補正（十分位）", "v_base",  [],           "intercept",      "mean", "ld_rank"),
    ("11 基準＋遺伝子単位で移動量（二次）",  "v_base",   ["lp", "ll"], "full_quadratic", "mean", "ld"),
    ("12 基準＋遺伝子単位で移動量（十分位）", "v_base",  ["lp", "ll"], "full_quadratic", "mean", "ld_rank"),
    ("13 向きから共通成分＋遺伝子単位の移動量", "v_uperp", ["lp", "ll"], "full_quadratic", "mean", "ld_rank"),
]


def run(cells: pd.DataFrame, elig: dict, n_perm: int = 200000) -> dict:
    nullset = set(diag.NULLS_B500)
    res = {}
    for name, col, xcols, form, agg, post in VARIANTS:
        t0 = time.time()
        tg, tu = {}, {}
        for s in diag.SEEDS:
            cs = cells[cells.seed == s]
            nl = {h: hd for h, hd in cs[cs.gene.isin(nullset)].groupby("gene")}
            kw = dict(xcols=xcols, form=form, agg=agg)
            if not xcols:
                kw["xcols"] = []
            f = diag.gene_stats(cs, nl, elig[s], col, **kw)
            fu = diag.gene_stats(cs, nl, elig[s], "m_uniform", **kw)
            if post != "none":
                f = gene_level_residual(f, nullset, post)
                fu = gene_level_residual(fu, nullset, post)
                # 順位 q は補正後の残差で引き直す
                for t in (f, fu):
                    nn = t[t.gene.isin(nullset)].set_index("gene").D_raw
                    qs = []
                    for g, d in zip(t.gene, t.D_raw):
                        e = elig[s].get(g, [])
                        qs.append(float((nn.reindex(e).to_numpy() < d).mean())
                                  if len(e) >= 20 else np.nan)
                    t["q"] = qs
            tg[s], tu[s] = f, fu
        res[name] = diag.evaluate(tg, tu, cells, n_perm=n_perm)
        res[name]["秒"] = round(time.time() - t0, 1)
        M = res[name]
        print(f"{name:42s} M5 {[round(M['M5'][s],4) for s in diag.SEEDS]} "
              f"M2max {[round(M['M2_max_abs'][s],3) for s in diag.SEEDS]} "
              f"M4 {[round(M['M4'][s],3) for s in diag.SEEDS]} "
              f"M3 {M['M3_n_in_range']}/4 "
              f"陽性 {M['positive_control71']['z']} "
              f"M1 {'o' if M['M1_pass'] else 'x'} 全 {'合格' if M['all_pass'] else '不合格'}",
              flush=True)
    return res


def main() -> None:
    cells = diag.load_cells()
    disp = diag.load_disp()
    states = json.loads((diag.CACHE / "states.json").read_text())
    elig = diag.get_elig(cells)
    p = diag.CACHE / "cells2.parquet"
    if p.exists():
        cells = pd.read_parquet(p)
    else:
        cells = add_columns(cells, disp, states)
        cells.to_parquet(p)
    res = run(cells, elig)
    (HERE / "variants.json").write_text(json.dumps(res, ensure_ascii=False, indent=1,
                                                   default=float) + "\n")


if __name__ == "__main__":
    main()
