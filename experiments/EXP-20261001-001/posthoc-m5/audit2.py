#!/usr/bin/env python3
"""C 判断との関係 — 検査が閾値ちょうどのとき、119 個のセットの平均 z がどれだけ歪むか。

既知の大きさの不具合を細胞ごとの指標に混ぜ、
(1) どの検査がどの大きさで鳴るか（感度）
(2) そのとき 119 個のセットの平均 z がどれだけ動くか（判断との関係）
を実測する。

**主判定のセットは使わない。** 歪みは、対照 500 個から「ある性質の順位で選んだ 119 個」と、
目標の軸の陽性対照 71 個（実在する、生物学で選ばれたセット）で測る。
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
from audit import one, cov_blocks  # noqa: E402

HERE = Path(__file__).parent
NULLSET = set(diag.NULLS_B500)


def main() -> None:
    cells = pd.read_parquet(diag.CACHE / "cells2.parquet")
    elig = diag.get_elig(cells)
    blocks = cov_blocks(cells)
    base = cells.v_base.to_numpy()
    sd = float(np.std(base))
    inj = {
        "移動量（M5 が狙う不具合）": np.log2(cells.dnorm.to_numpy()),
        "一様成分（M4 が狙う不具合）": cells.m_uniform.to_numpy(),
        "元の順位（M2 が狙う不具合）": cells.lp.to_numpy(),
        "入力の長さ（M2 が狙う不具合）": cells.ll.to_numpy(),
    }
    inj = {k: (v - v.mean()) / v.std() for k, v in inj.items()}
    grid = [-0.30, -0.20, -0.10, -0.05, -0.02, 0.0, 0.02, 0.05, 0.10, 0.20, 0.30]
    out = {}
    t0 = time.time()
    for name, X in inj.items():
        rows = []
        for a in grid:
            cells["v_inj"] = base + a * sd * X
            r = one(cells, elig, "v_inj", blocks)
            rows.append({
                "混ぜた大きさ（指標の標準偏差に対する比）": a,
                "M5": [round(r["M5"][s], 4) for s in diag.SEEDS],
                "M4": [round(r["M4"][s], 4) for s in diag.SEEDS],
                "M2_元の位置": [round(r["M2"][s]["lp"], 4) for s in diag.SEEDS],
                "M2_入力の長さ": [round(r["M2"][s]["ll"], 4) for s in diag.SEEDS],
                "M2_相対位置": [round(r["M2"][s]["relpos"], 4) for s in diag.SEEDS],
                "M1_上裾": [round(r["M1"][s]["hi"], 4) for s in diag.SEEDS],
                "M1_下裾": [round(r["M1"][s]["lo"], 4) for s in diag.SEEDS],
                "M3_範囲内": r["M3_n_in_range"],
                "合否": {k: r[k] for k in ("M1_pass", "M2_pass", "M3_pass",
                                           "M4_pass", "M5_pass")},
                "セット平均z": r["セット平均z"]})
            print(f"  {name} a={a:+.2f} M5={rows[-1]['M5']} "
                  f"{time.time()-t0:.0f}s", flush=True)
        out[name] = rows
    (HERE / "audit_part2.json").write_text(json.dumps(out, ensure_ascii=False, indent=1,
                                                      default=float) + "\n")


if __name__ == "__main__":
    main()
