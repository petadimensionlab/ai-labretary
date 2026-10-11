#!/usr/bin/env python3
"""EXP-20261007-001 — 2 台の本番の出力を results/ にまとめる（plan.md §5.4）。**値は変えない。**

  results-cf5b/run-81-L-1          → results/run-81-L-1
  results-c083/run-{82,83}-L-1     → results/run-{82,83}-L-1
  cells.csv を結合（行数が元の和と一致することを確かめる）
  device.json を results/device-{cf5b,c083}.json として残す

止める条件: 群の欠け・余計な群・遺伝子数が 587 でない群・results/ が既にある。

走らせ方:
    .venv/bin/python experiments/EXP-20261007-001/analysis/merge_results_011.py
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20261007-001"
PLAN = {"cf5b": [81], "c083": [82, 83]}
N_GENES = 587


def main() -> None:
    dst = E / "results"
    if dst.exists():
        raise SystemExit(f"**止める**: {dst} が既にある（試し実行の出力が混ざっていないか確かめる）")
    frames, rep = [], {}
    for host, seeds in PLAN.items():
        src = E / f"results-{host}"
        runs = sorted(p.name for p in src.glob("run-*-L-1"))
        want = [f"run-{s}-L-1" for s in seeds]
        if runs != want:
            raise SystemExit(f"**止める**: {src} の群が想定と違う（{runs} / 想定 {want}）")
        for r in runs:
            n = len([p for p in (src / r / "genes").iterdir() if p.is_dir()])
            if n != N_GENES:
                raise SystemExit(f"**止める**: {src / r} の遺伝子数が {n}（想定 {N_GENES}）")
            rep[r] = {"機体": host, "遺伝子数": n}
        c = pd.read_csv(src / "cells.csv", dtype={"donor_id": str})
        frames.append(c)
    dst.mkdir()
    for host, seeds in PLAN.items():
        src = E / f"results-{host}"
        for s in seeds:
            shutil.copytree(src / f"run-{s}-L-1", dst / f"run-{s}-L-1")
        shutil.copy2(src / "device.json", dst / f"device-{host}.json")
        if (src / "timing.jsonl").exists():
            shutil.copy2(src / "timing.jsonl", dst / f"timing-{host}.jsonl")
    cells = pd.concat(frames, ignore_index=True)
    if len(cells) != sum(len(f) for f in frames):
        raise SystemExit("**止める**: cells.csv の結合で行数が合わない")
    cells.to_csv(dst / "cells.csv", index=False)
    rep["cells.csv の行数"] = int(len(cells))
    (E / "work/merge_results_011.json").write_text(json.dumps(rep, ensure_ascii=False, indent=1) + "\n")
    print(json.dumps(rep, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
