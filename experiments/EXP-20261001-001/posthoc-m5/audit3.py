#!/usr/bin/env python3
"""C の続き — 相対位置の不具合を混ぜた掃引（M2 用）と、
各検査が閾値ちょうどのときのセット平均 z の歪みの計算。

元の順位 lp と入力の長さ ll をそのまま混ぜても何も起きない。
交絡を除く回帰が lp・ll の完全二次式なので、**定義により厳密に吸収される**。
M2 が実際に反応しうるのは、その二次式で表せない形（相対位置など）の不具合である。
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
from audit import one  # noqa: E402

HERE = Path(__file__).parent
NULLSET = set(diag.NULLS_B500)
K = diag.CR["K_set"]
GRID = [-0.30, -0.20, -0.10, -0.05, 0.0, 0.05, 0.10, 0.20, 0.30]
# 露出の段階（対照 500 個を、その性質の順位で並べた 119 個の塊）
PS = [0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9]


def blocks_for(cells: pd.DataFrame) -> dict:
    s = diag.SEEDS[0]
    cs = cells[(cells.seed == s) & cells.gene.isin(NULLSET)]
    g = cs.groupby("gene").agg(ld=("dnorm", "median"), un=("m_uniform", "mean"),
                               rp=("relpos", "median"))
    out, expo = {}, {}
    for cov in ("ld", "un", "rp"):
        order = g[cov].rank(method="first").sort_values().index.to_numpy()
        for p in PS:
            st = max(0, min(500 - K, int(round(p * 500 - K / 2))))
            genes = list(order[st:st + K])
            out[f"{cov}_p{p}"] = genes
            expo[f"{cov}_p{p}"] = round(float(np.mean(
                [(g[cov].to_numpy() < g[cov][x]).mean() for x in genes])), 3)
    rg = np.random.default_rng(3)
    out["無作為119"] = list(rg.choice(sorted(NULLSET), K, replace=False))
    return out, expo


def sweep(cells, elig, blocks, X, label):
    base = cells.v_base.to_numpy()
    sd = float(np.std(base))
    X = (X - X.mean()) / X.std()
    rows = []
    for a in GRID:
        cells["v_inj"] = base + a * sd * X
        r = one(cells, elig, "v_inj", blocks)
        rows.append({"a": a, "M5": [r["M5"][s] for s in diag.SEEDS],
                     "M4": [r["M4"][s] for s in diag.SEEDS],
                     "M2_lp": [r["M2"][s]["lp"] for s in diag.SEEDS],
                     "M2_ll": [r["M2"][s]["ll"] for s in diag.SEEDS],
                     "M2_relpos": [r["M2"][s]["relpos"] for s in diag.SEEDS],
                     "M1_hi": [r["M1"][s]["hi"] for s in diag.SEEDS],
                     "M1_lo": [r["M1"][s]["lo"] for s in diag.SEEDS],
                     "M3": r["M3_n_in_range"],
                     "pass": {k: r[k] for k in ("M1_pass", "M2_pass", "M3_pass",
                                                "M4_pass", "M5_pass")},
                     "sets": r["セット平均z"]})
        print(f"  {label} a={a:+.2f}", flush=True)
    return rows


def slope(rows, stat, setname):
    """統計量 1 単位あたり、そのセットの平均 z が何標準偏差動くか（最小二乗の傾き）。"""
    x = np.array([np.mean(r[stat]) for r in rows])
    y = np.array([np.mean(r["sets"][setname]) for r in rows])
    b = np.polyfit(x, y, 1)[0]
    return round(float(b), 3)


def main():
    cells = pd.read_parquet(diag.CACHE / "cells2.parquet")
    elig = diag.get_elig(cells)
    blocks, expo = blocks_for(cells)
    out = {"露出の実測（塊の平均分位）": expo}

    sw = {}
    sw["移動量"] = sweep(cells, elig, blocks, np.log2(cells.dnorm.to_numpy()), "移動量")
    sw["一様成分"] = sweep(cells, elig, blocks, cells.m_uniform.to_numpy(), "一様成分")
    sw["相対位置"] = sweep(cells, elig, blocks, cells.relpos.to_numpy(), "相対位置")
    (HERE / "audit_part3_sweeps.json").write_text(
        json.dumps(sw, ensure_ascii=False, default=float))

    # --- 傾き: 統計量 1 単位あたりのセット平均 z の動き
    tab = {}
    for stat, inj, cov in (("M5", "移動量", "ld"), ("M4", "一様成分", "un"),
                           ("M2_relpos", "相対位置", "rp")):
        rows = sw[inj]
        tab[stat] = {f"露出 p={p}（実測 {expo[f'{cov}_p{p}']}）": slope(rows, stat, f"{cov}_p{p}")
                     for p in PS}
        tab[stat]["無作為な 119 個"] = slope(rows, stat, "無作為119")
        tab[stat]["目標の軸の陽性対照71"] = slope(rows, stat, "目標の軸の陽性対照71")
    out["統計量 1 単位あたりのセット平均 z の傾き"] = tab

    # --- 閾値ちょうどのときの歪み
    thr = {"M5": diag.CR["M5_abs_rho_movement_max"], "M4": diag.CR["M4_abs_rho_uniform_max"],
           "M2_relpos": diag.CR["M2_abs_rho_max"]}
    dec = {"セットの絶対値条件": diag.ABS.get("INT"), "検出下限": None}
    pw = json.loads((diag.E / "work/power_005.json").read_text())["sets"]
    dec["検出下限"] = pw["INT"]["detection_floor_with_gene_sd"]
    out["守りたい判断"] = dec
    bias = {}
    for stat, t in thr.items():
        cov = {"M5": "ld", "M4": "un", "M2_relpos": "rp"}[stat]
        # 露出 p=0.5 の塊を基準にした差（純粋な露出の効果）
        s50 = tab[stat][f"露出 p=0.5（実測 {expo[f'{cov}_p0.5']}）"]
        bias[stat] = {"閾値": t}
        for p in PS:
            sl = tab[stat][f"露出 p={p}（実測 {expo[f'{cov}_p{p}']}）"] - s50
            bias[stat][f"露出 p={p} のセットの歪み"] = round(float(abs(sl) * t), 4)
            bias[stat][f"露出 p={p} の歪み / 判断の線 {dec['セットの絶対値条件']}"] = \
                round(float(abs(sl) * t / dec["セットの絶対値条件"]), 3)
        sl_pc = tab[stat]["目標の軸の陽性対照71"] - s50
        bias[stat]["実在のセット（陽性対照71）の歪み"] = round(float(abs(sl_pc) * t), 4)
    out["閾値ちょうどのときのセット平均 z の歪み"] = bias

    # --- 検査の感度: どの大きさの不具合で鳴るか
    sens = {}
    for stat, inj in (("M5", "移動量"), ("M4", "一様成分"), ("M2_relpos", "相対位置")):
        rows = sw[inj]
        sens[inj] = [{"a": r["a"], stat: round(float(np.mean(r[stat])), 4),
                      "M1_上裾": round(float(np.mean(r["M1_hi"])), 3),
                      "M3_範囲内": r["M3"],
                      "合否": "".join("o" if v else "x" for v in r["pass"].values())}
                     for r in rows]
    out["感度（混ぜた大きさと統計量・合否）"] = sens

    (HERE / "audit_part3.json").write_text(json.dumps(out, ensure_ascii=False, indent=1,
                                                      default=float) + "\n")
    print(json.dumps(out, ensure_ascii=False, indent=1, default=float))


if __name__ == "__main__":
    main()
