#!/usr/bin/env python3
"""EXP-20261003-001 — 機体の一致確認（**本番の前に通す。不合格なら本番を開始しない**）。

## この実験でとくに重要な理由

新しい 6 個の z は、**EXP-20261001-001 / EXP-20261002-002 で測った対照 1,000 個の分布**で
標準化する。既存の 1,217 遺伝子はすべて **NVIDIA GB10（CUDA / bf16）**で回したもので、
この実験は **Mac（MPS / bf16）**で回す。**機体が違えば、6 個だけが系統的にずれる**。
ずれは「新しい 6 個が動いた／動かない」という結論とまったく同じ形で現れ、見分けられない。

→ **既存の対照 12 個を Mac で回し直し、EXP-20261001-001 の保存済みの出力と突き合わせる。**

| # | 見るもの | 合格条件 |
|---|---|---|
| H-1 | 摂動前後の埋め込みの差（細胞ごと・成分ごと） | 最大の絶対差 <= 1e-6 |
| H-2 | 細胞の顔ぶれ（cell_uid の集合と順序） | 完全一致 |
| H-3 | 状態の重心 4 本 | 最大の絶対差 = 0（rows_file が同じなので同一になるはず） |
| H-4 | 主指標（向きへの射影）の遺伝子ごとの中央値 | 最大の絶対差 <= 1e-6 |

**不合格なら Mac では回さない。** その場合は既存と同じ機体（spark-cf5b / spark-c083）で回す。

使い方:
    .venv/bin/python experiments/EXP-20261003-001/analysis/machine_check_008.py \
        --results experiments/EXP-20261003-001/results-machinecheck
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20261003-001"
P01 = ROOT / "experiments/EXP-20261001-001"

TOL_EMB = 1e-6
TOL_PROJ = 1e-6


def main(results: str, seed: int = 81) -> None:
    new = ROOT / results
    old = P01 / "results"
    run = f"run-{seed}-L-1"
    out: dict = {"機体の出力": str(new), "比べる相手": str(old / run), "seed": seed,
                 "目的": "追加の対照を、既存の対照と同じ尺度で測れることの確認"}

    cn = pd.read_csv(new / "cells.csv", dtype={"donor_id": str})
    co = pd.read_csv(old / "cells.csv", dtype={"donor_id": str})
    co = co[co.run == run]
    genes = sorted(cn.gene.unique())
    out["突き合わせた遺伝子"] = {"n": len(genes), "genes": genes}

    # ---- H-3 状態の重心
    zn, zo = np.load(new / run / "state_emb_vectors.npz"), np.load(old / run / "state_emb_vectors.npz")
    d3 = {k: float(np.max(np.abs(zn[k].astype(float) - zo[k].astype(float)))) for k in zo.files}
    out["H-3 状態の重心"] = {"重心ごとの最大差": d3, "最大差": max(d3.values()),
                             "pass": bool(max(d3.values()) == 0.0)}

    # ---- H-1 / H-2 / H-4
    h1, h2, h4 = [], [], []
    for g in genes:
        pn = cn[cn.gene == g].reset_index(drop=True)
        po = co[co.gene == g].reset_index(drop=True)
        h2.append({"gene": g, "n_new": len(pn), "n_old": len(po),
                   "cell_uid 一致": bool(len(pn) == len(po)
                                          and (pn.cell_uid.to_numpy() == po.cell_uid.to_numpy()).all())})
        on = np.load(new / run / "genes" / g / "emb/orig.npy").astype(np.float64)
        qn = np.load(new / run / "genes" / g / "emb/pert.npy").astype(np.float64)
        oo = np.load(old / run / "genes" / g / "emb/orig.npy").astype(np.float64)
        qo = np.load(old / run / "genes" / g / "emb/pert.npy").astype(np.float64)
        if on.shape != oo.shape:
            h1.append({"gene": g, "形が違う": [list(on.shape), list(oo.shape)]})
            continue
        h1.append({"gene": g,
                   "orig 最大差": float(np.max(np.abs(on - oo))),
                   "pert 最大差": float(np.max(np.abs(qn - qo)))})

        # 主指標（出発点 → 目標の向きへの射影 ÷ 移動量）。**重心は既存のものを使う**
        vg = zo["aStr"].astype(np.float64)
        vs = zo["FIB_MFIB"].astype(np.float64)
        u = vg / np.linalg.norm(vg) - vs / np.linalg.norm(vs)
        u = u / np.linalg.norm(u)
        mn = (qn - on) @ u / np.linalg.norm(qn - on, axis=1)
        mo = (qo - oo) @ u / np.linalg.norm(qo - oo, axis=1)
        h4.append({"gene": g, "主指標の中央値の差": float(abs(np.median(mn) - np.median(mo)))})

    out["H-1 埋め込み"] = {
        "遺伝子ごと": h1,
        "最大差": max([max(x.get("orig 最大差", 0), x.get("pert 最大差", 0)) for x in h1] or [np.inf]),
        "pass": bool(all("形が違う" not in x and max(x["orig 最大差"], x["pert 最大差"]) <= TOL_EMB
                         for x in h1))}
    out["H-2 細胞の顔ぶれ"] = {"遺伝子ごと": h2, "pass": bool(all(x["cell_uid 一致"] for x in h2))}
    out["H-4 主指標"] = {"遺伝子ごと": h4,
                         "最大差": max([x["主指標の中央値の差"] for x in h4] or [np.inf]),
                         "pass": bool(all(x["主指標の中央値の差"] <= TOL_PROJ for x in h4))}
    out["pass"] = bool(all(out[k]["pass"] for k in
                           ("H-1 埋め込み", "H-2 細胞の顔ぶれ", "H-3 状態の重心", "H-4 主指標")))
    out["判定"] = ("合格。追加の対照を既存の対照と同じ尺度で測れる" if out["pass"]
                   else "**不合格。本番を開始しない。** Labretary へ報告して指示を待つ")

    name = f"machine_check_008-{Path(results).name}.json"
    (E / "work" / name).write_text(json.dumps(out, ensure_ascii=False, indent=1, default=float) + "\n")
    print(json.dumps({k: out[k] for k in ("機体の出力", "pass", "判定")}
                     | {k: {"pass": out[k]["pass"], "最大差": out[k].get("最大差")}
                        for k in ("H-1 埋め込み", "H-2 細胞の顔ぶれ", "H-3 状態の重心", "H-4 主指標")},
                     ensure_ascii=False, indent=1, default=float))
    if not out["pass"]:
        raise SystemExit("**不合格。本番を開始しない。**")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--seed", type=int, default=81)
    main(**vars(ap.parse_args()))
