#!/usr/bin/env python3
"""EXP-20261002-001 — EXP-20261001-001 のデータの再評価（**摂動は 1 回も回さない**）。

**この再評価で主判定の値は変わらない。** 変えるのは健全性の検査だけである。
主判定（介入セット 119 個）のセット平均 z は 0.208 / 0.098 / 0.142 のままで、
絶対値の条件 0.41 にも検出下限 0.604 にも届かない。**この値はすでに見えている。**
決まるのは「**その値を解釈してよいか**」だけである（plan.md §2）。

入力は EXP-20261001-001 が残した遺伝子ごとの表 `results/eval/gene_seed_005.csv` だけで、
**摂動も、埋め込みの読み直しも行わない**。最初に主判定・陽性対照・中立性を計算し直し、
**事前登録済みの `criteria_005.json` と一致しなければ止まる**（再現の検査）。

走らせ方:
    .venv/bin/python experiments/EXP-20261002-001/analysis/reevaluate_006.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

ROOT = Path(__file__).resolve().parents[3]
E = Path(__file__).resolve().parents[1]
SRC = ROOT / "experiments/EXP-20261001-001"
sys.path.insert(0, str(ROOT / "scripts/analysis"))
sys.path.insert(0, str(SRC / "analysis"))
from bias_budget import bias_budget  # noqa: E402
from evaluate_005 import sham_sets, set_null  # noqa: E402

BLOCKED = {"ENSG00000141510", "ENSG00000118785", "ENSG00000148346"}   # TP53 / SPP1 / LCN2


def main() -> None:
    cfg = yaml.safe_load(open(E / "config.yaml"))
    rc = cfg["reevaluation"]
    bb = cfg["bias_budget"]
    src_cfg = yaml.safe_load(open(SRC / "config.yaml"))
    a = src_cfg["analysis"]
    cr, ABS = a["criteria"], a["set_abs_min"]
    seeds = list(rc["seeds"])

    gj = json.loads((SRC / "work/genes_005.json").read_text())["genes"]
    ens = {g["symbol"]: g["ensembl"] for g in gj}
    sym = {v: k for k, v in ens.items()}
    bias_goal = {ens[g["symbol"]]: g["bias_goal"] for g in gj}
    subs = json.loads((SRC / "work/subsets_005.json").read_text())["subsets"]
    sets_def = {k: [ens[x["symbol"]] for x in v] for k, v in subs.items()}
    nulls = [ens[x["symbol"]] for x in
             json.loads((SRC / "work/null_pool_005.json").read_text())["versions"]["B_500"]]

    g = pd.read_csv(SRC / "results/eval/gene_seed_005.csv")
    assert not (set(g.gene) & BLOCKED), "盲検中の遺伝子が入っている"

    # ---------------------------------------------------------- 遺伝子ごとの表（z と交絡の候補）
    tabs, Z = {}, {}
    for s in seeds:
        x = g[(g.seed == s) & (g.col == "m_goal")].set_index("gene")
        u = g[(g.seed == s) & (g.col == "m_uniform")].set_index("gene").D_raw
        n = x.loc[[k for k in nulls if k in x.index]]
        t = pd.DataFrame({
            "gene": x.index,
            "z": ((x.D_raw - n.D_raw.mean()) / n.D_raw.std(ddof=1)).to_numpy(),
            "q": x.q.to_numpy(), "n_cells": x.n_cells.to_numpy(),
            "dnorm_med": x.dnorm_med.to_numpy(), "lp_med": x.lp_med.to_numpy(),
            "ll_med": x.ll_med.to_numpy(), "relpos_med": x.relpos_med.to_numpy(),
            "uniform": u.reindex(x.index).to_numpy()}).reset_index(drop=True)
        tabs[s] = t
        Z[s] = t.set_index("gene").z
    nulls = [k for k in nulls if k in set(tabs[seeds[0]].gene)]
    out = {"実験": "EXP-20261002-001", "種類": "再評価（摂動は回していない）",
           "入力": "experiments/EXP-20261001-001/results/eval/gene_seed_005.csv",
           "seed": seeds, "対照": len(nulls)}

    # ---------------------------------------------------------- セットの規則（事前登録と同じ）
    rng = np.random.default_rng(cr["null_draw_seed"])
    znull = np.column_stack([Z[s].reindex(sorted(nulls)).to_numpy() for s in seeds])

    def set_test(genes, name):
        gs = [x for x in genes if x in Z[seeds[0]].index]
        obs = np.array([float(np.mean([Z[s][x] for x in gs])) for s in seeds])
        nd = set_null(znull, len(gs), cr["n_permutation"], rng)
        thr = np.array([float(np.quantile(nd[:, i], cr["set_percentile"]))
                        for i in range(len(seeds))])
        am = ABS.get(name, -1e9)
        return {"k": len(gs), "set_mean_z": [round(float(x), 3) for x in obs],
                "threshold": [round(float(x), 3) for x in thr], "abs_min": am,
                "pass": bool(all(o > t for o, t in zip(obs, thr)) and all(o >= am for o in obs))}

    P1 = set_test(sets_def["INT"], "INT")
    P4 = set_test(sets_def["positive_control71"], "positive_control71")
    neu = [x for x in sets_def["neutral_adopted"] if x in Z[seeds[0]].index]
    obs_n = np.array([float(np.mean([Z[s][x] for x in neu])) for s in seeds])
    nd_n = set_null(znull, len(neu), cr["n_permutation"], rng)
    pct_n = [float((nd_n[:, i] < obs_n[i]).mean()) for i in range(len(seeds))]
    P5 = {"k": len(neu), "set_mean_z": [round(float(x), 3) for x in obs_n],
          "null_percentile": [round(x, 4) for x in pct_n],
          "n_seeds_in_range": int(sum(cr["neutral_range"][0] <= x <= cr["neutral_range"][1]
                                      for x in pct_n))}
    P5["pass"] = bool(P5["n_seeds_in_range"] >= cr["neutral_min_seeds"])

    # ---------------------------------------------------------- 再現の検査（**ここで止まりうる**）
    ref = json.loads((SRC / "results/eval/criteria_005.json").read_text())
    rep = {"主判定 P-1": {"今回": P1["set_mean_z"], "前回": ref["P1_primary"]["set_mean_z"]},
           "陽性対照 P-4": {"今回": P4["set_mean_z"],
                            "前回": ref["P4_positive_control"]["set_mean_z"]},
           "中立性 P-5": {"今回": P5["set_mean_z"],
                          "前回": ref["P5_neutral"]["neutral_adopted"]["set_mean_z"]}}
    tol = rc["reproduction_tolerance"]
    bad = [k for k, v in rep.items()
           if max(abs(np.array(v["今回"]) - np.array(v["前回"]))) > tol]
    rep["許容差"] = tol
    rep["一致"] = not bad
    out["再現の検査"] = rep
    if bad:
        (E / "results/criteria_006.json").write_text(
            json.dumps(out, ensure_ascii=False, indent=1, default=float))
        raise SystemExit(f"前回の値を再現できなかった: {bad}。**ここで止める**")

    out["主判定 P-1（**この再評価で変わらない**）"] = P1
    out["検出下限"] = json.loads((SRC / "work/power_005.json").read_text())["sets"]["INT"][
        "detection_floor_with_gene_sd"]

    # ---------------------------------------------------------- S-1 歪みの予算（**新しい本体**）
    S1 = {str(s): bias_budget(tabs[s], sets_def["INT"], nulls,
                              covariates=bb["covariates"], budget=bb["budget"],
                              decision_line=bb["decision_line"]) for s in seeds}
    out["S-1 介入セットの歪みの予算"] = {
        "per_seed": S1,
        "上側": [S1[str(s)]["upper_joint"] for s in seeds],
        "B": [S1[str(s)]["B_joint"] for s in seeds],
        "予算": bb["budget"],
        "pass": bool(all(S1[str(s)]["pass"] for s in seeds))}

    # ---------------------------------------------------------- S-2 偽のセット 200 組
    tgt = tabs[seeds[0]].set_index("gene").loc[
        [x for x in sets_def["INT"] if x in set(tabs[seeds[0]].gene)]].copy()
    tgt["symbol"] = [sym[x] for x in tgt.index]
    tgt["bias_goal"] = [bias_goal[x] for x in tgt.index]
    pool = tabs[seeds[0]].set_index("gene").loc[nulls].copy()
    pool["symbol"] = [sym[x] for x in pool.index]
    pool["bias_goal"] = [bias_goal[x] for x in pool.index]
    rng2 = np.random.default_rng(cr["null_draw_seed"] + 1)
    shams = sham_sets(pool.reset_index(), tgt.reset_index(), len(tgt),
                      cr["n_sham_sets"], rng2, cr["sham_match_bias"])
    k_pass = 0
    sham_z = []
    for ss in shams:
        gs = [ens[x] for x in ss]
        r = set_test(gs, "INT")
        sham_z.append(r["set_mean_z"])
        k_pass += int(r["pass"])
    p_chance = float(ref["chance_pass_rate"]["INT"]["rate"])
    N = len(shams)
    pval = float(stats.binomtest(k_pass, N, max(p_chance, 1e-9), alternative="greater").pvalue)
    out["S-2 偽のセットの二項検定"] = {
        "N": N, "p（単一の偶然合格率。対照で実測）": p_chance,
        "期待合格数": round(N * p_chance, 4), "k（実際の合格数）": k_pass,
        "片側 p 値": round(pval, 4), "偽セットの平均 z": np.round(
            np.array(sham_z).mean(axis=0), 4).tolist(),
        "pass": bool(pval >= bb["sham_pvalue_min"])}

    # ---------------------------------------------------------- S-3 陽性対照（歪みを引いても立つか）
    S3b = {str(s): bias_budget(tabs[s], sets_def["positive_control71"], nulls,
                               covariates=bb["covariates"], budget=bb["budget"],
                               decision_line=bb["decision_line"]) for s in seeds}
    adj = [round(abs(P4["set_mean_z"][i]) - S3b[str(s)]["upper_joint"], 3)
           for i, s in enumerate(seeds)]
    out["S-3 陽性対照 P-4"] = {
        "セット平均 z": P4["set_mean_z"], "規則の合否": P4["pass"],
        "歪みの上側": [S3b[str(s)]["upper_joint"] for s in seeds],
        "歪みを引いた値": adj, "必要な値": ABS["positive_control71"],
        "pass": bool(P4["pass"] and all(x >= ABS["positive_control71"] for x in adj))}
    out["S-4 中立性 P-5"] = P5

    # ---------------------------------------------------------- 記述のみ（判定に使わない）
    desc = {}
    for s in seeds:
        n = tabs[s][tabs[s].gene.isin(set(nulls))]
        desc[str(s)] = {
            "M1 較正（上裾 / 下裾）": [round(float((n.q >= 0.90).mean()), 4),
                                       round(float((n.q < 0.10).mean()), 4)],
            "M2 元の順位": round(float(stats.spearmanr(n.lp_med, n.z).statistic), 4),
            "M2 入力の長さ": round(float(stats.spearmanr(n.ll_med, n.z).statistic), 4),
            "M2 相対位置": round(float(stats.spearmanr(n.relpos_med, n.z).statistic), 4),
            "M4 一様成分": round(float(stats.spearmanr(n.uniform, n.z).statistic), 4),
            "M5 移動量": round(float(stats.spearmanr(n.dnorm_med, n.z).statistic), 4)}
    desc["M3 個々の中立対照の順位"] = {
        sym[x]: [round(float(tabs[s].set_index("gene").q.get(x, np.nan)), 4) for s in seeds]
        for x in sets_def["neutral_adopted"]}
    out["記述のみ（判定に使わない）"] = desc

    # ---------------------------------------------------------- 判定の提案
    health = all(out[k]["pass"] for k in ("S-1 介入セットの歪みの予算",
                                          "S-2 偽のセットの二項検定",
                                          "S-3 陽性対照 P-4", "S-4 中立性 P-5"))
    out["健全性"] = {"S-1": out["S-1 介入セットの歪みの予算"]["pass"],
                     "S-2": out["S-2 偽のセットの二項検定"]["pass"],
                     "S-3": out["S-3 陽性対照 P-4"]["pass"],
                     "S-4": out["S-4 中立性 P-5"]["pass"], "すべて合格": health}
    if not health:
        prop = "判定不能（健全性の検査が不合格）"
    elif P1["pass"]:
        prop = "支持"
    else:
        prop = f"否定（検出下限 {out['検出下限']} 標準偏差以上の効果は認められなかった）"
    out["判定の提案"] = prop
    out["注"] = ("判定は PI（関門 G3）が行う。この道具は提案までにとどまる。"
                 "また、これはコンピュータ上の遺伝子摂動の結果であり、"
                 "実験的に確認された生物学的事実ではない。")

    (E / "results").mkdir(exist_ok=True)
    (E / "results/criteria_006.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
    print(json.dumps({k: out[k] for k in
                      ("再現の検査", "主判定 P-1（**この再評価で変わらない**）", "健全性",
                       "判定の提案")}, ensure_ascii=False, indent=1, default=float))


if __name__ == "__main__":
    main()
