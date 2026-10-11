#!/usr/bin/env python3
"""EXP-20261011-001 — 幾何の事前確認（HYP-20261011-01 §3.2・§9.2）。摂動は 1 回も使わない。

保存済みの状態重心（A2: EXP-20261007-001 / A1: EXP-20261004-001）だけから、
2 本の向き（u_ada・u_deg）のなす角と、差の向き（u_diff）・共通の向き（u_common）の直交性を計算する。

- `u_ada = unit(goal_ada) − unit(start)`、`u_deg = unit(goal_deg) − unit(start)`
- `u_diff = (u_ada − u_deg) / ‖u_ada − u_deg‖`、`u_common = (u_ada + u_deg) / ‖u_ada + u_deg‖`
- **追加の着手条件（仮説 §9.2）**: `sin(なす角) < 0.31`（角 < 18 度）なら着手しない。
  差の向きが残す成分は sin(角) であり、小さすぎると差の向きが雑音に埋もれるため。

**これは結果ではない。** 状態の重心は設定（著者ラベル）の関数であり、摂動の結果を 1 つも使わない。
出力は事前登録のハッシュに含める（仮説 §9・タスク 3）。
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "experiments/EXP-20261011-001/work/geometry_014.json"

# データごとの設定。angle の判定（sin ≥ 0.31）はどちらも満たす（仮説 §3.2 の実測: A2 0.822 / A1 0.980）。
DATASETS = {
    "A2": {
        "seed_files": {s: ROOT / f"experiments/EXP-20261007-001/results/run-{s}-L-1/state_emb_vectors.npz"
                       for s in (81, 82, 83)},
        "start": "FIB", "goal_ada": "aStr", "goal_deg": "dFIB",
    },
    "A1": {
        "seed_files": {s: ROOT / f"experiments/EXP-20261004-001/results/run-{s}-L-1/state_emb_vectors.npz"
                       for s in (81, 82, 83)},
        "start": "FIB_MFIB", "goal_ada": "aStr", "goal_deg": "dFIB",
    },
}


def unit(v: np.ndarray) -> np.ndarray:
    return v / np.linalg.norm(v)


def main() -> None:
    res: dict = {
        "tool": "experiments/EXP-20261011-001/analysis/geometry_014.py",
        "注": "保存済みの状態重心だけから計算。摂動の結果は 1 つも使っていない。事前登録の一部（仮説 §9.2）。",
        "判定の規則": "sin(2 本の向きのなす角) < 0.31（角 < 18 度）なら着手しない。",
        "データごと": {},
    }
    ok_all = True
    for name, d in DATASETS.items():
        per_seed = {}
        for s, f in d["seed_files"].items():
            if not f.exists():
                raise SystemExit(f"重心がありません: {f}")
            z = np.load(str(f))
            ua = unit(z[d["goal_ada"]].astype(np.float64)) - unit(z[d["start"]].astype(np.float64))
            ua = ua / np.linalg.norm(ua)
            ud = unit(z[d["goal_deg"]].astype(np.float64)) - unit(z[d["start"]].astype(np.float64))
            ud = ud / np.linalg.norm(ud)
            udiff = (ua - ud) / np.linalg.norm(ua - ud)
            ucom = (ua + ud) / np.linalg.norm(ua + ud)
            cos_a = float(np.clip(np.dot(ua, ud), -1.0, 1.0))
            ang_deg = float(np.degrees(np.arccos(cos_a)))
            per_seed[str(s)] = {
                "cos(なす角)": round(cos_a, 6),
                "なす角（度）": round(ang_deg, 4),
                "sin(なす角)": round(float(np.sqrt(max(0.0, 1.0 - cos_a ** 2))), 6),
                "u_diff・u_common の内積（直交なら 0）": round(float(np.dot(udiff, ucom)), 12),
            }
        sin_min = min(v["sin(なす角)"] for v in per_seed.values())
        ang_min = min(v["なす角（度）"] for v in per_seed.values())
        ok = bool(sin_min >= 0.31)
        ok_all &= ok
        res["データごと"][name] = {
            "重心の出どころ": {k: str(v.relative_to(ROOT)) for k, v in d["seed_files"].items()},
            "seed ごと": per_seed,
            "sin の最小（3 群）": round(sin_min, 6),
            "なす角の最小（3 群）": round(ang_min, 4),
            "着手条件 sin >= 0.31": ok,
        }
    res["全体の着手条件"] = ok_all
    res["判定"] = "着手してよい（不合格なし）" if ok_all else "着手しない（2 本の向きが近すぎる）"
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, ensure_ascii=False, indent=1) + "\n")
    print(json.dumps(res, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()