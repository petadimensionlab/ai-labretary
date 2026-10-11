#!/usr/bin/env python3
"""指標の健全性 M1〜M5 という検査そのものの点検（PI 指示 2026-10-02）。

B 感度と特異度 / C 判断との関係 / E 重複 を実測する。

**主判定（介入 119 個）とその部分集合は一切読み込まない**（控えの段階で捨ててある）。
セットの歪みは、対照 500 個から選んだ 119 個のセットと、目標の軸の陽性対照 71 個で測る。
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
NULLSET = set(diag.NULLS_B500)
CR = diag.CR
K = CR["K_set"]          # 119


# --------------------------------------------------------------- 1 回ぶんの評価
def one(cells: pd.DataFrame, elig: dict, col: str, cov_blocks: dict) -> dict:
    """細胞ごとの指標 col について、M1〜M5 と、対照から選んだ 119 個のセットの平均 z を返す。"""
    M1, M2, M4, M5, Zs = {}, {}, {}, {}, {}
    q_of = {}
    for s in diag.SEEDS:
        cs = cells[cells.seed == s]
        nl = {h: hd for h, hd in cs[cs.gene.isin(NULLSET)].groupby("gene")}
        f = diag.gene_stats(cs, nl, elig[s], col, ["lp", "ll"])
        fu = diag.gene_stats(cs, nl, elig[s], "m_uniform", ["lp", "ll"])
        N = f[f.gene.isin(NULLSET)]
        q = N.q.dropna()
        M1[s] = {"n": int(len(q)), "hi": float((q >= 0.90).mean()),
                 "lo": float((q < 0.10).mean())}
        M2[s] = {k: float(stats.spearmanr(N[f"{k}_med"], N.D_raw).statistic)
                 for k in ("lp", "ll", "relpos")}
        M5[s] = float(stats.spearmanr(N.dnorm_med, N.D_raw).statistic)
        U = fu.set_index("gene").reindex(N.gene).D_raw
        M4[s] = float(stats.spearmanr(N.D_raw, U).statistic)
        z = (f.set_index("gene").D_raw - N.D_raw.mean()) / N.D_raw.std(ddof=1)
        Zs[s] = z
        q_of[s] = f.set_index("gene").q
    out = {"M1": M1, "M2": M2, "M4": M4, "M5": M5}
    out["M1_pass"] = all(v["n"] >= CR["M1_min_n"] and v["hi"] <= CR["M1_tail_frac_max"]
                         and v["lo"] <= CR["M1_tail_frac_max"] for v in M1.values())
    out["M2_pass"] = sum(all(abs(x) <= CR["M2_abs_rho_max"] for x in v.values())
                         for v in M2.values()) >= CR["M2_min_seeds"]
    out["M4_pass"] = all(abs(v) <= CR["M4_abs_rho_uniform_max"] for v in M4.values())
    out["M5_pass"] = all(abs(v) <= CR["M5_abs_rho_movement_max"] for v in M5.values())
    lo, hi = CR["M3_range"]
    M3 = {diag.SYM[g]: [float(q_of[s].get(g, np.nan)) for s in diag.SEEDS]
          for g in diag.SETS["neutral_adopted"]}
    out["M3"] = M3
    out["M3_n_in_range"] = int(sum(
        sum(1 for x in v if not np.isnan(x) and lo <= x <= hi) >= CR["M3_min_seeds"]
        for v in M3.values()))
    out["M3_pass"] = bool(out["M3_n_in_range"] >= CR["M3_adopted_min_genes"])
    # ---- セットの平均 z（対照から選んだ 119 個。主判定のセットではない）
    sets = {}
    for name, genes in cov_blocks.items():
        sets[name] = [round(float(np.mean([Zs[s][g] for g in genes])), 4) for s in diag.SEEDS]
    sets["目標の軸の陽性対照71"] = [round(float(np.mean(
        [Zs[s][g] for g in diag.SETS["positive_control71"]])), 4) for s in diag.SEEDS]
    out["セット平均z"] = sets
    return out


def cov_blocks(cells: pd.DataFrame) -> dict:
    """対照 500 個を、ある性質の順位で並べて 119 個の連続した塊を取る（平均分位 p のセット）。
    **主判定のセットではない。対照だけで作った、性質が偏ったセットである。**"""
    s = diag.SEEDS[0]
    cs = cells[(cells.seed == s) & cells.gene.isin(NULLSET)]
    g = cs.groupby("gene").agg(ld=("dnorm", "median"), lp=("lp", "median"),
                               ll=("ll", "median"), un=("m_uniform", "mean"))
    out = {}
    for cov in ("ld", "lp", "ll", "un"):
        order = g[cov].rank(method="first").sort_values().index.to_numpy()
        for p in (0.5, 0.6, 0.75, 0.9):
            start = int(round(p * 500 - K / 2))
            start = max(0, min(500 - K, start))
            out[f"{cov}_p{p}"] = list(order[start:start + K])
    rg = np.random.default_rng(3)
    out["無作為119"] = list(rg.choice(sorted(NULLSET), K, replace=False))
    return out


def main() -> None:
    cells = pd.read_parquet(diag.CACHE / "cells2.parquet")
    disp = diag.load_disp()
    elig = diag.get_elig(cells)
    blocks = cov_blocks(cells)
    out = {}

    # ---- セットの偏りの実例（対照の中での分位）
    s = diag.SEEDS[0]
    cs = cells[cells.seed == s]
    gg = cs.groupby("gene").agg(ld=("dnorm", "median"), lp=("lp", "median"),
                                ll=("ll", "median"), un=("m_uniform", "mean"))
    nul = gg.loc[sorted(NULLSET)]
    def pct(genes, cov):
        v = nul[cov].to_numpy()
        return round(float(np.mean([(v < gg[cov][g]).mean() for g in genes])), 3)
    out["セットの性質の偏り（対照の中での平均分位）"] = {
        "目標の軸の陽性対照71": {c: pct(diag.SETS["positive_control71"], c)
                                 for c in ("ld", "lp", "ll", "un")},
        "恒常発現3": {c: pct(diag.SETS["housekeeping"], c) for c in ("ld", "lp", "ll", "un")},
        "中立対照パネル7": {c: pct(diag.SETS["negctrl_panel"], c) for c in ("ld", "lp", "ll", "un")},
        "無作為な119個の分位の標準偏差": round(float(0.2887 / np.sqrt(K)), 4)}

    # ---- B 特異度: 「きれいな指標」= 無作為な向きへの射影
    t0 = time.time()
    rg = np.random.default_rng(2026)
    ens = []
    nseed_dirs = 60
    for i in range(nseed_dirs):
        r = rg.normal(size=1152)
        r /= np.linalg.norm(r)
        v = np.full(len(cells), np.nan)
        for a in range(0, len(cells), 40000):
            b = min(a + 40000, len(cells))
            v[a:b] = (disp[a:b].astype(np.float64) @ r) / cells.dnorm.to_numpy()[a:b]
        cells["v_rnd"] = v
        ens.append(one(cells, elig, "v_rnd", blocks))
        if (i + 1) % 10 == 0:
            print(f"  無作為な向き {i+1}/{nseed_dirs}  {time.time()-t0:.0f}s", flush=True)
    (HERE / "ensemble_random_dirs.json").write_text(
        json.dumps(ens, ensure_ascii=False, default=float))

    def arr(key, sub=None):
        if sub:
            return np.array([[e[key][s][sub] for s in diag.SEEDS] for e in ens])
        return np.array([[e[key][s] for s in diag.SEEDS] for e in ens])
    spec = {
        "本数": nseed_dirs,
        "M5 の分布": {"平均": round(float(arr("M5").mean()), 4),
                      "標準偏差": round(float(arr("M5").std(ddof=1)), 4),
                      "|値|>0.10 の割合": round(float((np.abs(arr("M5")) > 0.10).mean()), 3)},
        "M4 の分布": {"平均": round(float(arr("M4").mean()), 4),
                      "標準偏差": round(float(arr("M4").std(ddof=1)), 4),
                      "|値|>0.40 の割合": round(float((np.abs(arr("M4")) > 0.40).mean()), 3)},
        "M2(入力の長さ) の分布": {"平均": round(float(arr("M2", "ll").mean()), 4),
                                  "標準偏差": round(float(arr("M2", "ll").std(ddof=1)), 4),
                                  "|値|>0.25 の割合": round(float((np.abs(arr("M2", "ll")) > 0.25).mean()), 3)},
        "M2(元の位置) の分布": {"平均": round(float(arr("M2", "lp").mean()), 4),
                                "標準偏差": round(float(arr("M2", "lp").std(ddof=1)), 4)},
        "M2(相対位置) の分布": {"平均": round(float(arr("M2", "relpos").mean()), 4),
                                "標準偏差": round(float(arr("M2", "relpos").std(ddof=1)), 4)},
        "各検査に合格した本数": {k: int(sum(e[k] for e in ens))
                                 for k in ("M1_pass", "M2_pass", "M3_pass", "M4_pass", "M5_pass")},
        "5 つすべてに合格した本数": int(sum(all(e[k] for k in
                                               ("M1_pass", "M2_pass", "M3_pass", "M4_pass", "M5_pass"))
                                           for e in ens)),
        "M3 の範囲内の個数の分布": {str(k): int(sum(e["M3_n_in_range"] == k for e in ens))
                                    for k in range(5)}}
    out["B_特異度_無作為な向きへの射影"] = spec

    # ---- E 重複: 5 つの統計量どうしの相関（無作為な向きの集団で）
    cols = {"M5": arr("M5").ravel(), "M4": arr("M4").ravel(),
            "M2_lp": arr("M2", "lp").ravel(), "M2_ll": arr("M2", "ll").ravel(),
            "M2_relpos": arr("M2", "relpos").ravel(),
            "M1_hi": arr("M1", "hi").ravel(), "M1_lo": arr("M1", "lo").ravel()}
    names = list(cols)
    C = np.corrcoef(np.vstack([cols[n] for n in names]))
    out["E_検査どうしの相関（無作為な向き 60 本 × 3 seed）"] = {
        names[i]: {names[j]: round(float(C[i, j]), 3) for j in range(len(names))}
        for i in range(len(names))}

    (HERE / "audit_part1.json").write_text(json.dumps(out, ensure_ascii=False, indent=1,
                                                      default=float) + "\n")
    print(json.dumps(out, ensure_ascii=False, indent=1, default=float))


if __name__ == "__main__":
    main()
