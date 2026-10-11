#!/usr/bin/env python3
"""EXP-20261003-001 — seed 83 を spark-c083 の結果に差し替える（Labretary の決定、案 A）。

## なぜ差し替えるか

評価が plan.md §6.1 の「状態の重心の検査」で止まった（中断条件 A-2）。
**seed 83 の重心 dMFIB だけ**が EXP-20261001-001 と一致しなかった
（最大差 0.0009765625 = bf16 の 1 刻み）。原因は機体の違いで、
**既存の seed 83 は spark-c083、最初に回したのは spark-cf5b** だった。

**事前登録を 1 文字も変えずに通すため、seed 83 だけを spark-c083 で回し直す。**
seed 81・82（spark-cf5b、照合済み）は**そのまま使う。回し直さない。**

## このスクリプトがすること（統合は機械的な作業で、判定基準に関わらない）

 1. 新しい seed 83 の中身を確かめる（6 遺伝子・埋め込みの欠け 0・盲検の混入 0）
 2. **重心が EXP-20261001-001 と最大差 0 になったか**を先に確かめる（ならなければ**止める**）
 3. 最初に回した seed 83（spark-cf5b）を `results-seed83-cf5b/` へ**退避する。消さない**
 4. 新しい seed 83 を `results/run-83-L-1/` に置く
 5. `cells.csv` と `summary.csv` を、seed 83 の行だけ差し替えて作り直す
 6. 差し替えの記録を `results/seed83_swap_report.json` に残す

実行: .venv/bin/python experiments/EXP-20261003-001/analysis/merge_seed83_008.py
"""
from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20261003-001"
P01 = ROOT / "experiments/EXP-20261001-001"
RUN = "run-83-L-1"
BLOCKED = {"ENSG00000141510": "TP53", "ENSG00000118785": "SPP1", "ENSG00000148346": "LCN2"}
TARGETS = {"ENSG00000107796": "ACTA2", "ENSG00000149591": "TAGLN", "ENSG00000130176": "CNN1",
           "ENSG00000133392": "MYH11", "ENSG00000133110": "POSTN", "ENSG00000113721": "PDGFRB"}


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main() -> None:
    new_root = E / "results-seed83-c083"
    res = E / "results"
    new_run = new_root / RUN
    if not new_run.exists():
        raise SystemExit(f"{new_run} がありません。seed 83 の回し直しが回収されていない")

    rep: dict = {"決定": "Labretary 2026-10-03、案 A（事前登録を 1 文字も変えない）",
                 "差し替えるもの": f"results/{RUN}",
                 "新しい出どころ": "spark-c083（既存の EXP-20261001-001 seed83 と同じ機体）",
                 "退避先": "results-seed83-cf5b/（**消さない**）",
                 "回し直していないもの": ["run-81-L-1（spark-cf5b）", "run-82-L-1（spark-cf5b）"]}

    # ---- 1. 中身の確認
    genes = sorted(p.name for p in (new_run / "genes").iterdir() if p.is_dir())
    miss = [g for g in genes if not ((new_run / "genes" / g / "emb/orig.npy").exists()
                                     and (new_run / "genes" / g / "emb/pert.npy").exists())]
    rep["1. 新しい seed 83 の中身"] = {
        "遺伝子数": len(genes), "顔ぶれが対象 6 個と一致": sorted(genes) == sorted(TARGETS),
        "埋め込みの欠け": len(miss), "盲検の混入": sorted(set(genes) & set(BLOCKED)),
        "重心ファイル": (new_run / "state_emb_vectors.npz").exists()}
    if not (len(genes) == 6 and sorted(genes) == sorted(TARGETS) and not miss
            and not (set(genes) & set(BLOCKED)) and (new_run / "state_emb_vectors.npz").exists()):
        raise SystemExit(f"**止める**: 新しい seed 83 の中身が揃っていない: {rep['1. 新しい seed 83 の中身']}")

    # ---- 2. 重心（**ここで止まりうる**）
    a = np.load(P01 / "results" / RUN / "state_emb_vectors.npz")
    b = np.load(new_run / "state_emb_vectors.npz")
    d = {k: float(np.max(np.abs(a[k].astype(float) - b[k].astype(float)))) for k in a.files}
    old_run = res / RUN
    d_old = {}
    if old_run.exists():
        c = np.load(old_run / "state_emb_vectors.npz")
        d_old = {k: float(np.max(np.abs(a[k].astype(float) - c[k].astype(float)))) for k in a.files}
    rep["2. 状態の重心（EXP-20261001-001 との最大の絶対差）"] = {
        "新しい seed 83（spark-c083）": d, "最大": max(d.values()),
        "最初の seed 83（spark-cf5b。参考）": d_old,
        "pass": max(d.values()) == 0.0}
    if max(d.values()) != 0.0:
        (E / "work/merge_seed83_008.json").write_text(
            json.dumps(rep, ensure_ascii=False, indent=1, default=float) + "\n")
        raise SystemExit(f"**止める**: 回し直しても重心が一致しない {d}")

    # ---- 3. 退避（**消さない**）
    keep = E / "results-seed83-cf5b"
    if old_run.exists():
        keep.mkdir(exist_ok=True)
        if (keep / RUN).exists():
            shutil.rmtree(keep / RUN)
        shutil.move(str(old_run), str(keep / RUN))
        rep["3. 退避"] = {"から": f"results/{RUN}", "へ": f"results-seed83-cf5b/{RUN}",
                          "注": "最初に spark-cf5b で回した seed 83。**消していない**"}

    # ---- 4. 置き換え
    shutil.copytree(new_run, old_run, copy_function=shutil.copy2)
    rep["4. 置き換え"] = {"遺伝子数": len(list((old_run / "genes").iterdir()))}

    # ---- 5. cells.csv / summary.csv の作り直し
    cells = pd.read_csv(res / "cells.csv", dtype={"donor_id": str})
    n_before = len(cells)
    cells = cells[cells.run != RUN]
    rows = []
    for g in genes:
        c = pd.read_csv(old_run / "genes" / g / "cells.csv", dtype={"donor_id": str})
        c.insert(0, "gene", g)
        c.insert(0, "run", RUN)
        c["role"] = "target"
        rows.append(c)
    new_cells = pd.concat(rows, ignore_index=True)
    cells = pd.concat([cells, new_cells[cells.columns]], ignore_index=True)
    cells.to_csv(res / "cells.csv", index=False)

    summ = pd.read_csv(res / "summary.csv")
    cols = list(summ.columns)
    summ = summ[summ.run != RUN]
    new_sm = pd.read_csv(new_root / "summary.csv")
    summ = pd.concat([summ, new_sm[cols]], ignore_index=True)
    summ.to_csv(res / "summary.csv", index=False)

    rep["5. 表の作り直し"] = {
        "cells.csv の行数": {"前": int(n_before), "後": int(len(cells))},
        "seed 83 の行数": int((cells.run == RUN).sum()),
        "seed・遺伝子ごとの細胞数":
            {f"{r}/{TARGETS[g]}": int(n) for (r, g), n in cells.groupby(["run", "gene"]).size().items()},
        "summary.csv の行数": int(len(summ)),
        "盲検の混入": sorted(set(cells.gene) & set(BLOCKED)),
        "cells.csv の sha256": sha256(res / "cells.csv"),
        "summary.csv の sha256": sha256(res / "summary.csv")}

    (E / "work/merge_seed83_008.json").write_text(
        json.dumps(rep, ensure_ascii=False, indent=1, default=float) + "\n")
    (res / "seed83_swap_report.json").write_text(
        json.dumps(rep, ensure_ascii=False, indent=1, default=float) + "\n")
    print(json.dumps(rep, ensure_ascii=False, indent=1, default=float))


if __name__ == "__main__":
    main()
