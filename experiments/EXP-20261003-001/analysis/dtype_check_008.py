#!/usr/bin/env python3
"""EXP-20261003-001 — **数値型（bf16 / fp32）の切り分け**。

機体照合が Mac で落ちた原因は MPS ではなく bf16 だった（2026-10-04）。
ここでは同じ 12 遺伝子を 4 つの条件で突き合わせ、**どの条件の違いがどれだけの差を生むか**を
主指標の水準まで降ろして数値で出す。

| 条件 | 機体 | 装置 | 数値型 |
|---|---|---|---|
| `P01`  | NVIDIA GB10 | CUDA | bf16 |（既存 1,217 遺伝子の基準）
| `cf5b` | spark-cf5b  | CUDA | bf16 |
| `mac`  | Mac         | MPS  | bf16 |
| `fp32` | Mac         | MPS  | **fp32** |

**主指標の差は、対照 12 遺伝子の主指標の中央値のばらつき（標準偏差）で割って報告する。**
絶対値ではなく「対照何個ぶんのずれか」が、結論を動かすかどうかを決めるため。

使い方: `.venv/bin/python experiments/EXP-20261003-001/analysis/dtype_check_008.py`
"""
from __future__ import annotations

import itertools, json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20261003-001"
RUN, SEED = "run-81-L-1", 81

CONDS = {
    "P01  (GB10/CUDA/bf16)":  ROOT / "experiments/EXP-20261001-001/results",
    "cf5b (CUDA/bf16)":       E / "results-machinecheck-cf5b",
    "mac  (MPS/bf16)":        E / "results-machinecheck-mac",
    "fp32 (MPS/fp32)":        E / "results-machinecheck-fp32",
}


def load(root: Path):
    cells = pd.read_csv(root / "cells.csv", dtype={"donor_id": str})
    if "run" in cells:
        cells = cells[cells.run == RUN]
    z = np.load(root / RUN / "state_emb_vectors.npz")
    return cells, z


def main() -> None:
    have = {k: v for k, v in CONDS.items() if (v / RUN / "state_emb_vectors.npz").exists()}
    missing = [k for k in CONDS if k not in have]
    if missing:
        print("見つからない条件:", missing)

    data = {k: load(v) for k, v in have.items()}
    ref = ROOT / "experiments/EXP-20261001-001/results"
    zref = np.load(ref / RUN / "state_emb_vectors.npz")
    # 主指標の向きは**どの条件でも同じ重心**（P01）から作る。向きまで動かすと差の出どころが混ざる
    u = zref["aStr"].astype(np.float64) / np.linalg.norm(zref["aStr"].astype(np.float64)) \
        - zref["FIB_MFIB"].astype(np.float64) / np.linalg.norm(zref["FIB_MFIB"].astype(np.float64))
    u /= np.linalg.norm(u)

    genes = sorted(set.intersection(*[set(c.gene.unique()) for c, _ in data.values()]))
    print(f"突き合わせた遺伝子: {len(genes)} 個 / 条件: {len(have)}\n")

    med = {}      # 条件 -> 遺伝子 -> 主指標の中央値
    emb = {}      # 条件 -> 遺伝子 -> (orig, pert)
    for k, root in have.items():
        med[k], emb[k] = {}, {}
        for g in genes:
            o = np.load(root / RUN / "genes" / g / "emb/orig.npy").astype(np.float64)
            p = np.load(root / RUN / "genes" / g / "emb/pert.npy").astype(np.float64)
            emb[k][g] = (o, p)
            med[k][g] = float(np.median((p - o) @ u / np.linalg.norm(p - o, axis=1)))

    base = list(have)[0]
    sd = float(np.std(list(med[base].values()), ddof=1))
    print(f"対照 {len(genes)} 個の主指標の中央値のばらつき（{base} での SD）= {sd:.4f}")
    print("→ 以下の「対照 SD 比」が 1 を超えると、対照 1 個ぶん以上の系統的なずれになる\n")

    print("%-24s %-24s %12s %12s %10s" % ("条件 A", "条件 B", "埋め込み最大差", "主指標最大差", "対照 SD 比"))
    rows = []
    for a, b in itertools.combinations(have, 2):
        de = max(max(float(np.abs(emb[a][g][i] - emb[b][g][i]).max()) for i in (0, 1)) for g in genes)
        dm = max(abs(med[a][g] - med[b][g]) for g in genes)
        print("%-24s %-24s %12.6f %12.6f %10.2f" % (a, b, de, dm, dm / sd))
        rows.append({"A": a, "B": b, "埋め込み最大差": de, "主指標最大差": dm, "対照SD比": dm / sd})

    out = {"目的": "bf16 / fp32 と機体の違いが、主指標でどれだけの差になるかを切り分ける",
           "seed": SEED, "突き合わせた遺伝子": len(genes),
           "主指標の向き": "P01 の重心から作った固定の向き（条件ごとに作り直さない）",
           f"対照 {len(genes)} 個の主指標の中央値の SD": sd,
           "条件どうしの差": rows,
           "見つからない条件": missing}
    (E / "work/dtype_check_008.json").write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n")
    print(f"\n書き出し: {E / 'work/dtype_check_008.json'}")


if __name__ == "__main__":
    main()
