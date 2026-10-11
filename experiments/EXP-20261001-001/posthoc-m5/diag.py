#!/usr/bin/env python3
"""M5 の事後診断 — 指標の候補ごとに M1〜M5 と目標の軸の陽性対照 71 個を計算する。

**主判定（介入セット 119 個）とその入れ子は、控えを作る段階で捨ててある**ので、
ここでは計算できない（build_cache.py を参照）。

事前登録済みのファイルは読むだけで、一切書き換えない。
判定の閾値（M1 0.2 / M2 0.25 / M3 [0.1,0.9] / M4 0.40 / M5 0.10 / 陽性対照 0.53）は動かさない。
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

ROOT = Path("/Users/petadimensionlab/workspace/research/ai-labretary")
E = ROOT / "experiments/EXP-20261001-001"
CACHE = Path("/private/tmp/claude-503/-Users-petadimensionlab-workspace-research-ai-labretary/"
             "3a2e253b-7169-44cc-ba54-2207da1ae321/scratchpad/cache")

CFG = yaml.safe_load(open(E / "config.yaml"))
CR = CFG["analysis"]["criteria"]
PM = CFG["analysis"]["pos_matched_permutation"]
ABS = CFG["analysis"]["set_abs_min"]
GJ = json.loads((E / "work/genes_005.json").read_text())["genes"]
ENS = {g["symbol"]: g["ensembl"] for g in GJ}
SYM = {v: k for k, v in ENS.items()}
SUBS = json.loads((E / "work/subsets_005.json").read_text())["subsets"]
SETS = {k: [ENS[x["symbol"]] for x in v] for k, v in SUBS.items()
        if k in ("positive_control71", "positive_control65", "positive_control39",
                 "neutral_adopted", "housekeeping", "negctrl_panel")}
NULLS_B500 = [ENS[x["symbol"]] for x in
              json.loads((E / "work/null_pool_005.json").read_text())["versions"]["B_500"]]
SEEDS = [81, 82, 83]


# ---------------------------------------------------------------- 設計行列と遺伝子ごとの残差
def design(cols: list[np.ndarray], form: str = "full_quadratic") -> np.ndarray:
    cols = [np.asarray(x, dtype=float) for x in cols]
    if not cols:
        raise ValueError("説明変数が空（定数のみは form='intercept' を使う）")
    out = [np.ones_like(cols[0])] + list(cols)
    if form == "linear":
        return np.column_stack(out)
    for i in range(len(cols)):
        for j in range(i, len(cols)):
            out.append(cols[i] * cols[j])
    return np.column_stack(out)


def gene_stats(cs: pd.DataFrame, nulls: dict, elig: dict, col: str, xcols: list[str],
               min_elig: int = 20, form: str = "full_quadratic",
               agg: str = "mean") -> pd.DataFrame:
    """遺伝子ごとに、対照（自分を除く）で当てた回帰からの残差 D_raw と対照内順位 q。
    evaluate_005 が呼ぶ movement_correction.gene_stats と同じ計算（xcols=['lp','ll'] のとき）。"""
    hs = list(nulls)
    mk = (lambda t: np.ones((len(t), 1))) if not xcols else \
         (lambda t: design([t[x].to_numpy() for x in xcols], form))
    Xn = {h: mk(hd) for h, hd in nulls.items()}
    yn = {h: hd[col].to_numpy() for h, hd in nulls.items()}
    G = {h: Xn[h].T @ Xn[h] for h in hs}
    B = {h: Xn[h].T @ yn[h] for h in hs}
    Xbar = np.vstack([Xn[h].mean(axis=0) for h in hs])
    ybar = np.array([yn[h].mean() for h in hs])
    hi = {h: i for i, h in enumerate(hs)}
    G_all, B_all = sum(G.values()), sum(B.values())
    rows = []
    for gene, tg in cs.groupby("gene"):
        Gm, Bm = (G_all - G[gene], B_all - B[gene]) if gene in G else (G_all, B_all)
        beta = np.linalg.solve(Gm, Bm)
        Xg = mk(tg)
        if agg == "mean":
            Dg = float(tg[col].to_numpy().mean() - Xg.mean(axis=0) @ beta)
            Dh = ybar - Xbar @ beta
        else:   # 残差の中央値（外れ値に強い集計）
            Dg = float(np.median(tg[col].to_numpy() - Xg @ beta))
            Dh = np.array([float(np.median(yn[h] - Xn[h] @ beta)) for h in hs])
        e = elig[gene]
        ii = np.array([hi[h] for h in e], dtype=int) if e else np.zeros(0, int)
        Dhe = Dh[ii]
        q = float((Dhe < Dg).mean()) if len(Dhe) >= min_elig else np.nan
        rows.append({"gene": gene, "D_raw": Dg, "q": q, "n_elig": len(Dhe), "n_cells": len(tg),
                     "lp_med": float(tg.lp.median()), "ll_med": float(tg.ll.median()),
                     "relpos_med": float(tg.relpos.median()),
                     "dnorm_med": float(tg.dnorm.median())})
    return pd.DataFrame(rows)


def eligibility(cs: pd.DataFrame, nulls: dict, tol: float, cover: float) -> dict:
    out = {}
    for gene, tg in cs.groupby("gene"):
        lp, ll = tg.lp.to_numpy(np.float32), tg.ll.to_numpy(np.float32)
        elig = []
        for h, hd in nulls.items():
            if h == gene:
                continue
            D = np.maximum(np.abs(lp[:, None] - hd.lp.to_numpy(np.float32)[None, :]),
                           np.abs(ll[:, None] - hd.ll.to_numpy(np.float32)[None, :]))
            if (D <= tol).any(axis=1).mean() >= cover:
                elig.append(h)
        out[gene] = elig
    return out


# ---------------------------------------------------------------- 判定
def set_null(z: np.ndarray, k: int, n_draw: int, rng) -> np.ndarray:
    return z[rng.integers(0, z.shape[0], size=(n_draw, k))].mean(axis=1)


def evaluate(tabs_goal: dict, tabs_uni: dict, cells: pd.DataFrame,
             n_perm: int = 200000) -> dict:
    """tabs_goal/tabs_uni: seed -> gene_stats の表。M1〜M5 と陽性対照のセット平均 z を返す。"""
    nullset = set(NULLS_B500)
    M = {}
    # ---- M1 / M2 / M5
    M1, M2, M5 = {}, {}, {}
    for s in SEEDS:
        P = tabs_goal[s]
        N = P[P.gene.isin(nullset)]
        q = N.q.dropna()
        M1[s] = {"n": int(len(q)), "hi": float((q >= 0.90).mean()), "lo": float((q < 0.10).mean())}
        M2[s] = {k: float(stats.spearmanr(N[f"{k}_med"], N.D_raw).statistic)
                 for k in ("lp", "ll", "relpos")}
        M5[s] = float(stats.spearmanr(N.dnorm_med, N.D_raw).statistic)
    M["M1"] = M1
    M["M1_pass"] = all(v["n"] >= CR["M1_min_n"] and v["hi"] <= CR["M1_tail_frac_max"]
                       and v["lo"] <= CR["M1_tail_frac_max"] for v in M1.values())
    M["M2"] = M2
    M["M2_max_abs"] = {s: max(abs(x) for x in v.values()) for s, v in M2.items()}
    M["M2_pass"] = sum(all(abs(x) <= CR["M2_abs_rho_max"] for x in v.values())
                       for v in M2.values()) >= CR["M2_min_seeds"]
    M["M5"] = M5
    M["M5_pass"] = all(abs(v) <= CR["M5_abs_rho_movement_max"] for v in M5.values())
    # ---- M3（採用した中立対照 4 個）
    lo, hi = CR["M3_range"]
    M3 = {}
    for g in SETS["neutral_adopted"]:
        qs = []
        for s in SEEDS:
            P = tabs_goal[s]
            r = P[P.gene == g]
            qs.append(float(r.q.iloc[0]) if len(r) else np.nan)
        M3[SYM[g]] = qs
    M["M3"] = M3
    n_in = sum(sum(1 for x in v if not np.isnan(x) and lo <= x <= hi) >= CR["M3_min_seeds"]
               for v in M3.values())
    M["M3_n_in_range"] = int(n_in)
    M["M3_pass"] = bool(n_in >= min(CR["M3_adopted_min_genes"], len(M3)))
    # ---- M4
    M4 = {}
    for s in SEEDS:
        P = tabs_goal[s]
        N = P[P.gene.isin(nullset)].set_index("gene")
        U = tabs_uni[s].set_index("gene").reindex(N.index)
        M4[s] = float(stats.spearmanr(N.D_raw, U.D_raw).statistic)
    M["M4"] = M4
    M["M4_pass"] = all(abs(v) <= CR["M4_abs_rho_uniform_max"] for v in M4.values())
    # ---- 標準化と陽性対照
    Z = {}
    for s in SEEDS:
        P = tabs_goal[s].copy()
        nul = P[P.gene.isin(nullset)]
        P["z"] = (P.D_raw - nul.D_raw.mean()) / nul.D_raw.std(ddof=1)
        Z[s] = P.set_index("gene").z
    rng = np.random.default_rng(CR["null_draw_seed"])
    znull = np.column_stack([Z[s].reindex(sorted(nullset)).to_numpy() for s in SEEDS])
    for name in ("positive_control71", "positive_control65", "positive_control39",
                 "neutral_adopted", "housekeeping", "negctrl_panel"):
        genes = [g for g in SETS[name] if g in Z[SEEDS[0]].index]
        obs = np.array([float(np.mean([Z[s][g] for g in genes])) for s in SEEDS])
        nd = set_null(znull, len(genes), n_perm, rng)
        thr = np.array([np.quantile(nd[:, i], CR["set_percentile"]) for i in range(len(SEEDS))])
        pct = [float((nd[:, i] < obs[i]).mean()) for i in range(len(SEEDS))]
        am = ABS.get(name, -1e9)
        M[name] = {"k": len(genes), "z": [round(float(x), 3) for x in obs],
                   "thr": [round(float(x), 3) for x in thr],
                   "pct": [round(x, 4) for x in pct],
                   "pass": bool(all(o > t for o, t in zip(obs, thr)) and all(o >= am for o in obs))}
    M["P4_pass"] = M["positive_control71"]["pass"]
    M["P5_pass"] = bool(sum(CR["neutral_range"][0] <= x <= CR["neutral_range"][1]
                            for x in M["neutral_adopted"]["pct"]) >= CR["neutral_min_seeds"])
    M["all_pass"] = bool(M["M1_pass"] and M["M2_pass"] and M["M3_pass"] and M["M4_pass"]
                         and M["M5_pass"] and M["P4_pass"] and M["P5_pass"])
    return M


def load_cells() -> pd.DataFrame:
    return pd.read_parquet(CACHE / "cells.parquet")


def load_disp() -> np.ndarray:
    return np.load(CACHE / "disp_kept.npy", mmap_mode="r")


def get_elig(cells: pd.DataFrame) -> dict:
    f = CACHE / "elig.json"
    if f.exists():
        return {int(k): v for k, v in json.loads(f.read_text()).items()}
    out = {}
    for s in SEEDS:
        cs = cells[cells.seed == s]
        nl = {h: hd for h, hd in cs[cs.gene.isin(NULLS_B500)].groupby("gene")}
        out[s] = eligibility(cs, nl, PM["tolerance_log2"], PM["min_cover_frac"])
        print(f"  比べられる対照 seed{s} 完了", flush=True)
    f.write_text(json.dumps(out))
    return out
