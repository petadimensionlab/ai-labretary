#!/usr/bin/env python3
"""EXP-20261005-001 着手前チェック — 特異性の統計量を、**対照遺伝子だけ**で較正する。

摂動は 1 回も回さない。EXP-20261004-001 が保存した埋め込み（orig.npy / pert.npy）と
状態の重心（state_emb_vectors.npz）を読むだけ。

**候補 8 個と陽性対照 71 個には一切触れない。** 読むのは対照 549 個だけである。
事前登録の閾値を、結果を見る前に実測から決めるための計算である。

出力: experiments/EXP-20261005-001/work/precheck_010.json
"""
from __future__ import annotations

import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "experiments/EXP-20261004-001"
OUT = ROOT / "experiments/EXP-20261005-001/work/precheck_010.json"
sys.path.insert(0, str(ROOT / "scripts/analysis"))
sys.path.insert(0, str(SRC / "analysis"))
from evaluate_009 import design_reduced, resolve_gene_sets  # noqa: E402
from movement_correction import design as design_full       # noqa: E402

SEEDS = [81, 82, 83]
MDE = 2.73  # 単一遺伝子の最小検出可能効果量（EXP-20260928-002 の実測）


def unit(v):
    return v / np.linalg.norm(v)


def directions(npz) -> dict[str, np.ndarray]:
    """重心から向きを作る。recipe は evaluate_009.load_cells と同一
    （重心をそれぞれ単位長にしてから差を取り、もう一度単位長にする）。"""
    h = {k: unit(npz[k].astype(np.float64)) for k in npz.files}
    s = h["FIB_MFIB"]
    d = {"goal": unit(h["aStr"] - s),
         "dFIB": unit(h["dFIB"] - s),
         "dMFIB": unit(h["dMFIB"] - s)}
    C = np.column_stack([d["dFIB"], d["dMFIB"]])
    Q, _ = np.linalg.qr(C)
    r = d["goal"] - Q @ (Q.T @ d["goal"])
    d["spec"] = unit(r)
    d["_resid_norm"] = np.linalg.norm(r)
    return d


def main() -> None:
    cfg = yaml.safe_load((SRC / "config.yaml").read_text())
    sets = resolve_gene_sets(cfg)
    nulls = set(sets["対照549"])
    assert not (nulls & set(sets["候補"])), "対照に候補が混ざっている"
    assert not (nulls & set(sets["陽性対照71"])), "対照に陽性対照が混ざっている"

    res = SRC / "results"
    cells = pd.read_csv(res / "cells.csv", dtype={"donor_id": str})
    cells = cells[cells.gene.isin(nulls)].reset_index(drop=True)
    cells["seed"] = cells.run.astype(str).str.extract(r"run-(\d+)", expand=False).astype(int)

    # lp（元の順位）・ll（入力の長さ）
    from datasets import load_from_disk
    from geneformer import TOKEN_DICTIONARY_FILE
    tokd = pickle.load(open(TOKEN_DICTIONARY_FILE, "rb"))
    ds = load_from_disk(str(ROOT / cfg["data"]["tokenized_dataset"]))
    uids = sorted(cells.cell_uid.unique())
    ids = {u: np.asarray(x) for u, x in zip(uids, ds.select(uids)["input_ids"])}
    pos = np.full(len(cells), np.nan)
    expressed = np.zeros(len(cells), dtype=bool)
    for gene, part in cells.groupby("gene", sort=False):
        t = tokd.get(gene)
        for j, u in zip(part.index.to_numpy(), part.cell_uid):
            hit = np.flatnonzero(ids[u] == t) if t is not None else np.empty(0, dtype=int)
            if len(hit):
                pos[j], expressed[j] = hit[0], True
    cells["expressed"] = expressed
    with np.errstate(divide="ignore", invalid="ignore"):
        cells["lp"] = np.log2(1 + pos)
    cells["ll"] = np.log2(cells.length.astype(float))

    # 向きごとの細胞ごとの射影
    names = ["goal", "dFIB", "dMFIB", "spec"]
    cols = {n: np.full(len(cells), np.nan) for n in names}
    dnorm = np.full(len(cells), np.nan)
    geom = {}
    for run, part_run in cells.groupby("run", sort=False):
        npz = np.load(res / run / "state_emb_vectors.npz")
        d = directions(npz)
        geom[run] = {"resid_norm": float(d.pop("_resid_norm")),
                     "cos(goal,dFIB)": float(d["goal"] @ d["dFIB"]),
                     "cos(goal,dMFIB)": float(d["goal"] @ d["dMFIB"]),
                     "cos(goal,spec)": float(d["goal"] @ d["spec"])}
        for gene, part in part_run.groupby("gene", sort=False):
            gd = res / run / "genes" / gene / "emb"
            o = np.load(gd / "orig.npy").astype(np.float64)
            p = np.load(gd / "pert.npy").astype(np.float64)
            idx = part.index.to_numpy()
            dd = p - o
            dn = np.linalg.norm(dd, axis=1)
            dnorm[idx] = dn
            with np.errstate(divide="ignore", invalid="ignore"):
                for n in names:
                    cols[n][idx] = (dd @ d[n]) / dn
    for n in names:
        cells[f"m_{n}"] = cols[n]
    cells["dnorm"] = dnorm
    zero = cells.dnorm <= 0
    frac_zero = float(zero.mean())
    cells = cells[~zero].reset_index(drop=True)

    # 層別補正 → 遺伝子ごとの D_raw → 対照 549 個で標準化
    from evaluate_009 import gene_stats_stratified
    z = {}
    for seed in SEEDS:
        cs = cells[cells.seed == seed]
        nd = {g: d for g, d in cs.groupby("gene")}
        for n in names:
            gs = gene_stats_stratified(cs, nd, f"m_{n}").set_index("gene")["D_raw"]
            z[(seed, n)] = (gs - gs.mean()) / gs.std(ddof=1)
    zdf = pd.DataFrame({f"{n}_s{s}": z[(s, n)] for s in SEEDS for n in names})

    def rate(col_prefix):
        hit = np.ones(len(zdf), dtype=bool)
        for s in SEEDS:
            hit &= (zdf[f"{col_prefix}_s{s}"] > MDE).to_numpy()
        return int(hit.sum()), float(hit.mean())

    out = {
        "実験": "EXP-20261005-001",
        "注": "摂動は回していない。EXP-20261004-001 の保存済み埋め込みを読んだだけ。"
              "読んだのは対照 549 個のみで、候補 8 個・陽性対照 71 個には触れていない。",
        "対照の個数": int(zdf.shape[0]),
        "向きの幾何（run ごと）": geom,
        "移動量 0 の割合": round(frac_zero, 6),
        "3 seed とも標準得点 > 2.73 の割合（偶然の水準）": {
            n: dict(zip(["k", "rate"], rate(n))) for n in names},
        "遺伝子ごとの標準得点の相関（3 seed 平均、対照 549）": {},
        "統計量の seed 間ばらつき（対照 549、3 seed の標準偏差の中央）": {},
    }
    mean3 = {n: zdf[[f"{n}_s{s}" for s in SEEDS]].mean(axis=1) for n in names}
    sd3 = {n: zdf[[f"{n}_s{s}" for s in SEEDS]].std(axis=1, ddof=1) for n in names}
    for a in names:
        out["統計量の seed 間ばらつき（対照 549、3 seed の標準偏差の中央）"][a] = round(float(sd3[a].median()), 4)
        for b in names:
            if a < b:
                out["遺伝子ごとの標準得点の相関（3 seed 平均、対照 549）"][f"{a} vs {b}"] = {
                    "pearson": round(float(mean3[a].corr(mean3[b])), 4),
                    "spearman": round(float(mean3[a].corr(mean3[b], method="spearman")), 4)}

    # 差の統計量（goal − 比べる相手）を対照で標準化したときの偶然の水準
    for k in ("dFIB", "dMFIB"):
        diff = {s: zdf[f"goal_s{s}"] - zdf[f"{k}_s{s}"] for s in SEEDS}
        zd = {s: (diff[s] - diff[s].mean()) / diff[s].std(ddof=1) for s in SEEDS}
        hit = np.ones(len(zdf), dtype=bool)
        for s in SEEDS:
            hit &= (zd[s] > MDE).to_numpy()
        out[f"差の統計量 goal−{k} を対照で標準化したときの偶然の水準"] = {
            "k": int(hit.sum()), "rate": round(float(hit.mean()), 5),
            "生の差の標準偏差（seed ごと）": [round(float(diff[s].std(ddof=1)), 4) for s in SEEDS]}

    zdf.to_csv(OUT.with_suffix(".csv"))
    OUT.write_text(json.dumps(out, ensure_ascii=False, indent=2))
    print(json.dumps(out, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
