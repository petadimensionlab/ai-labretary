#!/usr/bin/env python3
"""Q1 — 旧データ（EXP-20260930-001）と本番（EXP-20261001-001）で、
向きへの射影の M5 がなぜ違うのかを対照遺伝子だけで比べる。

**どちらの実験でも、対照 500 個・目標の軸の陽性対照・中立対照・恒常発現しか触らない。**
旧実験の T1・T2、本番の介入セットは読み込まない。
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

ROOT = Path("/Users/petadimensionlab/workspace/research/ai-labretary")
OLD = ROOT / "experiments/EXP-20260930-001"
NEW = ROOT / "experiments/EXP-20261001-001"


def old_tables():
    npool = json.loads((OLD / "work/null_pool_004.json").read_text())["versions"]["B_500"]
    gj = json.loads((OLD / "work/genes_004.json").read_text())["genes"]
    ens = {g["symbol"]: g["ensembl"] for g in gj}
    nulls = [ens[x["symbol"]] for x in npool]
    subs = json.loads((OLD / "work/subsets_004.json").read_text())["subsets"]
    pc = [ens[x["symbol"]] for x in subs["positive_control71"]]
    c = pd.read_parquet(OLD / "posthoc/cells_cache.parquet")
    keep = set(nulls) | set(pc)
    c = c[c.gene.isin(keep)].reset_index(drop=True)
    el = {}
    for s in diag.SEEDS:
        e = json.loads((OLD / f"posthoc/eligibility-seed{s}.json").read_text())
        el[s] = {g: [h for h in v if h in set(nulls)] for g, v in e.items() if g in keep}
    return c, nulls, pc, el


def m5_for(cells, nulls, el, col, xcols, form="full_quadratic"):
    out = {}
    for s in diag.SEEDS:
        cs = cells[cells.seed == s]
        nl = {h: hd for h, hd in cs[cs.gene.isin(nulls)].groupby("gene")}
        f = diag.gene_stats(cs, nl, el[s], col, xcols, form=form)
        N = f[f.gene.isin(nulls)]
        out[s] = round(float(stats.spearmanr(N.dnorm_med, N.D_raw).statistic), 4)
    return out


def gene_level(cells, nulls):
    """対照 500 個の遺伝子ごと: 移動量の中央値と、各状態への寄り・射影の平均の順位相関。"""
    res = {}
    for s in diag.SEEDS:
        cs = cells[(cells.seed == s) & cells.gene.isin(nulls)]
        g = cs.groupby("gene").agg(dnorm_med=("dnorm", "median"), m=("m_proj", "mean"),
                                   sg=("s_goal", "mean"), ss=("s_start", "mean"),
                                   sf=("s_dFIB", "mean"), sm=("s_dMFIB", "mean"),
                                   lp=("lp", "median"), ll=("ll", "median"),
                                   rp=("relpos", "median"), n=("dnorm", "size"))
        r = lambda a, b: round(float(stats.spearmanr(a, b).statistic), 3)
        res[s] = {"n_genes": len(g),
                  "移動量中央値の中央値": round(float(g.dnorm_med.median()), 4),
                  "移動量中央値の四分位": [round(float(g.dnorm_med.quantile(q)), 4)
                                           for q in (0.25, 0.75)],
                  "log2移動量の遺伝子間SD": round(float(np.log2(g.dnorm_med).std(ddof=1)), 4),
                  "rho(移動量, 射影の平均)": r(g.dnorm_med, g.m),
                  "rho(移動量, 目標への寄り)": r(g.dnorm_med, g.sg),
                  "rho(移動量, 出発点への寄り)": r(g.dnorm_med, g.ss),
                  "rho(移動量, 退化FIBへの寄り)": r(g.dnorm_med, g.sf),
                  "rho(移動量, 退化MFIBへの寄り)": r(g.dnorm_med, g.sm),
                  "rho(移動量, 元の位置)": r(g.dnorm_med, g.lp),
                  "rho(移動量, 入力の長さ)": r(g.dnorm_med, g.ll),
                  "rho(移動量, 相対位置)": r(g.dnorm_med, g.rp),
                  "細胞数の中央値": int(g.n.median())}
    return res


def common_dir(resdir: Path, nulls: list[str], cells: pd.DataFrame, seed: int) -> dict:
    run = f"run-{seed}-L-1"
    sz = np.load(resdir / run / "state_emb_vectors.npz")
    vg = sz["aStr"].astype(np.float64)
    vs = sz["FIB_MFIB"].astype(np.float64)
    u = vg / np.linalg.norm(vg) - vs / np.linalg.norm(vs)
    u = u / np.linalg.norm(u)
    tot = None
    n = 0
    for g in nulls:
        gd = resdir / run / "genes" / g / "emb"
        if not gd.exists():
            continue
        o = np.load(gd / "orig.npy").astype(np.float64)
        p = np.load(gd / "pert.npy").astype(np.float64)
        d = p - o
        dn = np.linalg.norm(d, axis=1)
        d = d[dn > 0] / dn[dn > 0, None]
        tot = d.sum(axis=0) if tot is None else tot + d.sum(axis=0)
        n += len(d)
    cbar = tot / n
    ch = cbar / np.linalg.norm(cbar)
    return {"n_cells": n, "‖共通の向きの平均‖（単位変位の平均の長さ）": round(float(np.linalg.norm(cbar)), 4),
            "û·ĉ": round(float(u @ ch), 4),
            "ĉ·目標の単位ベクトル": round(float(ch @ (vg / np.linalg.norm(vg))), 4),
            "ĉ·出発点の単位ベクトル": round(float(ch @ (vs / np.linalg.norm(vs))), 4),
            "目標と出発点のなす角 cos": round(float((vg @ vs) / np.linalg.norm(vg) / np.linalg.norm(vs)), 4)}


def main():
    out = {}
    # ---- 新（本番）
    cn = diag.load_cells()
    eln = diag.get_elig(cn)
    nulls_new = diag.NULLS_B500
    pc_new = diag.SETS["positive_control71"]
    # ---- 旧
    co, nulls_old, pc_old, elo = old_tables()
    out["対照プールの顔ぶれ"] = {
        "旧 B_500": len(nulls_old), "本番 B_500": len(nulls_new),
        "共通": len(set(nulls_old) & set(nulls_new)),
        "入れ替わった割合": round(1 - len(set(nulls_old) & set(nulls_new)) / 500, 3)}
    out["目標の軸の陽性対照 71 個"] = {"共通": len(set(pc_old) & set(pc_new))}
    out["M5（向きへの射影）"] = {
        "旧・補正なし": m5_for(co, nulls_old, elo, "m_proj", [], form="linear"),
        "旧・元の位置と入力の長さの完全二次": m5_for(co, nulls_old, elo, "m_proj", ["lp", "ll"]),
        "本番・補正なし": m5_for(cn, nulls_new, eln, "m_proj", [], form="linear"),
        "本番・元の位置と入力の長さの完全二次": m5_for(cn, nulls_new, eln, "m_proj", ["lp", "ll"])}
    out["対照 500 個の遺伝子ごとの様子（旧）"] = gene_level(co, nulls_old)
    out["対照 500 個の遺伝子ごとの様子（本番）"] = gene_level(cn, nulls_new)
    # 共通の 500 個だけで比べる
    both = sorted(set(nulls_old) & set(nulls_new))
    out["共通の対照だけで見た遺伝子ごとの様子（旧）"] = gene_level(co, both)
    out["共通の対照だけで見た遺伝子ごとの様子（本番）"] = gene_level(cn, both)
    out["共通の向き ĉ"] = {
        "旧": {s: common_dir(OLD / "results", nulls_old, co, s) for s in diag.SEEDS},
        "本番": {s: common_dir(NEW / "results", nulls_new, cn, s) for s in diag.SEEDS}}
    p = Path(__file__).parent / "q1_compare.json"
    p.write_text(json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
    print(json.dumps(out, ensure_ascii=False, indent=1, default=float))


if __name__ == "__main__":
    main()
