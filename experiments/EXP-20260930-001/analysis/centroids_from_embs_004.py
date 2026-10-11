#!/usr/bin/env python3
"""EXP-20260930-001 — 着手前チェック用の状態重心を、**既に計算してある埋め込みから**作り直す。

2026-09-30 の着手前チェック（第 3 版, `axis_feasibility_v3.json`）は、
`outputs/design/2026-09-30-annotation-goal-state/work-v2/axisfeas-embs-seed{81,82,83}/cell_embs.csv`
に保存された埋め込みから重心を作っている。**この埋め込みは摂動を 1 回も使っていない**（前向き計算だけ）。

陽性対照のセットを替えて着手前チェックを回し直すとき、**重心は marker の選び方に一切依存しない**ので、
同じ埋め込みから作り直せば、モデルの前向き計算をやり直さずに済む
（やり直すと Mac で約 25 分かかるうえ、bf16 の丸めで完全一致しない恐れがある）。

本スクリプトは `scripts/analysis/axis_feasibility.py` の `extract_centroids()` と
**同じ計算**（指定した行の埋め込みの平均）を行い、`centroids.mode: precomputed` で読める .npz にする。
出力後、第 3 版の重心の指標（軸の長さの比）が再現することを呼び出し側で確かめる。

出力: work/axisfeas_centroids-seed{81,82,83}.npz
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/Users/petadimensionlab/workspace/research/ai-labretary")
E = ROOT / "experiments/EXP-20260930-001"
EMB = ROOT / "outputs/design/2026-09-30-annotation-goal-state/work-v2"
SEEDS = [81, 82, 83]


def main() -> None:
    for s in SEEDS:
        f = EMB / f"axisfeas-embs-seed{s}/cell_embs.csv"
        if not f.exists():
            raise SystemExit(f"埋め込みがありません: {f}")
        df = pd.read_csv(f, index_col=0)   # 先頭の無名列は行番号なので落とす
        rows = json.loads((E / f"work/state_rows-seed{s}.json").read_text())["rows"]
        X = df.drop(columns=["row_id"]).to_numpy(dtype=np.float64)
        by_row = dict(zip(df["row_id"].astype(int).to_numpy(), range(len(df))))
        out = {}
        for name, keep in rows.items():
            idx = [by_row[int(i)] for i in keep]
            if len(idx) != len(keep):
                raise SystemExit(f"seed {s}: {name} の行が埋め込みに足りない")
            out[name] = X[idx].mean(axis=0)
        np.savez(E / f"work/axisfeas_centroids-seed{s}.npz", **out)
        print(f"seed {s}: " + ", ".join(f"{k} {len(v)}" for k, v in rows.items())
              + f" → work/axisfeas_centroids-seed{s}.npz（{len(out)} 本 × {X.shape[1]} 次元）")


if __name__ == "__main__":
    main()
