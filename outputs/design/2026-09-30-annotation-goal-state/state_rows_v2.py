#!/usr/bin/env python3
"""候補 (I) の 4 本の重心に使う細胞を seed ごとに決める（深さをそろえる）。

手順は experiments/EXP-20260929-001/analysis/state_rows_003.py と同一
（log(トークン長)の 1:1 最近傍マッチ、許容 0.10。guardrails 2026-09-27）。
入力（tokenized.dataset）は読み取りのみ。摂動は行わない。

  aFIB      目標
  FIB_MFIB  出発点（FIB + M-FIB）
  MYOF      参照（軸の長さの比の分母）
  dFIB      第 4 の重心
"""
from __future__ import annotations

import json
import random
from pathlib import Path

import numpy as np
from datasets import load_from_disk

ROOT = Path("/Users/petadimensionlab/workspace/research/ai-labretary")
OUT = ROOT / "outputs/design/2026-09-30-annotation-goal-state"
DS = ROOT / "experiments/EXP-20260929-001/work/tokenized.dataset"

SEEDS = [81, 82, 83]
DEPTH_TOL = 0.10
GOAL_MAX, GOAL_PER_DONOR = 950, 60
START_POOL_MAX = 3000
EXTRA_MAX, EXTRA_PER_DONOR = 900, 60

GOAL_LABELS = ["aFIB", "MYOF"]
START_LABELS = ["FIB", "M-FIB"]
EXTRA = {"dFIB": ["dFIB"], "dMFIB": ["dM-FIB"]}
GOAL, START = "aStr", "FIB_MFIB"


def balanced(idx, donors, n_max, per_donor, rng):
    order = list(idx)
    rng.shuffle(order)
    used, keep = {}, []
    for i in order:
        dn = donors[i]
        if per_donor is not None and used.get(dn, 0) >= per_donor:
            continue
        used[dn] = used.get(dn, 0) + 1
        keep.append(i)
        if n_max is not None and len(keep) >= n_max:
            break
    return sorted(keep)


def depth_match(log_len, pos, neg, tol):
    pi, ni = np.where(pos)[0], np.where(neg)[0]
    ni = ni[np.argsort(log_len[ni], kind="stable")]
    nv = log_len[ni]
    used = np.zeros(len(ni), bool)
    keep = []
    for i in pi:
        v = log_len[i]
        cand = np.where((~used) & (np.abs(nv - v) <= tol))[0]
        if len(cand):
            j = cand[np.argmin(np.abs(nv[cand] - v))]
            used[j] = True
            keep += [int(i), int(ni[j])]
    m = np.zeros(len(log_len), bool)
    m[keep] = True
    return m


def main() -> None:
    ds = load_from_disk(str(DS))
    sub = np.array([str(x) for x in ds["subclass_l2"]])
    donors = [str(x) for x in ds["donor_id"]]
    log_len = np.log(np.array(ds["length"], dtype=float))

    goal_idx = np.flatnonzero(np.isin(sub, GOAL_LABELS)).tolist()
    start_idx = np.flatnonzero(np.isin(sub, START_LABELS)).tolist()
    summary = {"seeds": SEEDS, "depth_tol": DEPTH_TOL,
               "pool": {"goal": len(goal_idx), "start": len(start_idx),
                        **{k: int(np.isin(sub, v).sum()) for k, v in EXTRA.items()}},
               "per_seed": {}}

    for s in SEEDS:
        rng = random.Random(s)
        gcap = balanced(goal_idx, donors, GOAL_MAX, GOAL_PER_DONOR, rng)
        scap = balanced(start_idx, donors, START_POOL_MAX, GOAL_PER_DONOR, rng)
        g0 = np.zeros(len(sub), bool); g0[gcap] = True
        s0 = np.zeros(len(sub), bool); s0[scap] = True
        m = depth_match(log_len, g0, s0, DEPTH_TOL)
        rows = {GOAL: sorted(np.flatnonzero(m & g0).tolist()),
                START: sorted(np.flatnonzero(m & s0).tolist())}

        sref = np.zeros(len(sub), bool); sref[rows[START]] = True
        extra = {}
        for name, labs in EXTRA.items():
            pool = np.flatnonzero(np.isin(sub, labs)).tolist()
            cap = balanced(pool, donors, EXTRA_MAX, EXTRA_PER_DONOR, random.Random(f"{s}-{name}"))
            p0 = np.zeros(len(sub), bool); p0[cap] = True
            mm = depth_match(log_len, p0, sref, DEPTH_TOL)
            extra[name] = sorted(np.flatnonzero(mm & p0).tolist())

        allrows = {**rows, **extra}
        assert len({i for v in allrows.values() for i in v}) == sum(len(v) for v in allrows.values()), \
            "4 本の重心で細胞が重なっている"
        (OUT / f"state_rows_v2-seed{s}.json").write_text(json.dumps(
            {"seed": s, "rows": rows, "extra_rows": extra}, ensure_ascii=False) + "\n")
        summary["per_seed"][str(s)] = {
            k: {"n_cells": len(v), "n_donors": len({donors[i] for i in v}),
                "median_log_len": float(np.median(log_len[v]))} for k, v in allrows.items()}
        print(s, {k: len(v) for k, v in allrows.items()}, flush=True)

    (OUT / "depth_match_summary_v2.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=1) + "\n")
    print("→", OUT / "depth_match_summary_v2.json")


if __name__ == "__main__":
    main()
