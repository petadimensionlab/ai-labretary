#!/usr/bin/env python3
"""EXP-20261001-001 — 4 本の状態重心に使う細胞（行番号）を seed ごとに決める。

手順は EXP-20260930-001 の `analysis/state_rows_004.py` と**完全に同一**
（log(トークン長)の 1:1 最近傍マッチ、許容 0.10。guardrails 2026-09-27「状態の重心は
深さをそろえた後の細胞から作る」）。**摂動は 1 回も行わない。**

**seed は 81 / 82 / 83 で、EXP-20260930-001 と同じものを使う**（plan.md §2 の「seed を
引き継ぐ理由」）。したがって行番号は前実験と 1 行も違わないはずで、違えば止める。
前実験の出力 `experiments/EXP-20260930-001/work/state_rows-seed{S}.json` と照合する。

出力:
  work/state_rows-seed{81,82,83}.json   {"rows": {4 状態: [行番号]}}
  work/depth_match_summary_005.json     状態ごとの細胞数・提供者数・入力の長さの中央値
"""
from __future__ import annotations

import json
import random
import sys
from pathlib import Path

import numpy as np

ROOT = Path("/Users/petadimensionlab/workspace/research/ai-labretary")
E = ROOT / "experiments/EXP-20261001-001"
PREV = ROOT / "experiments/EXP-20260930-001/work"      # **読み取りのみ**

SEEDS = [81, 82, 83]
DEPTH_TOL = 0.10
GOAL_MAX, GOAL_PER_DONOR = 950, 60
START_POOL_MAX = 3000
EXTRA_MAX, EXTRA_PER_DONOR = 900, 60

GOAL, START = "aStr", "FIB_MFIB"
EXTRA = ["dFIB", "dMFIB"]


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
    from datasets import load_from_disk

    ds = load_from_disk(str(E / "work/tokenized.dataset"))
    st = np.array([str(x) for x in ds["fib_anno_state"]])
    donors = [str(x) for x in ds["donor_id"]]
    log_len = np.log(np.array(ds["length"], dtype=float))

    goal_idx = np.flatnonzero(st == GOAL).tolist()
    start_idx = np.flatnonzero(st == START).tolist()
    summary = {"seeds": SEEDS, "depth_tol": DEPTH_TOL,
               "pool": {GOAL: len(goal_idx), START: len(start_idx),
                        **{k: int((st == k).sum()) for k in EXTRA}},
               "per_seed": {}}

    for s in SEEDS:
        rng = random.Random(s)
        gcap = balanced(goal_idx, donors, GOAL_MAX, GOAL_PER_DONOR, rng)
        scap = balanced(start_idx, donors, START_POOL_MAX, GOAL_PER_DONOR, rng)
        g0 = np.zeros(len(st), bool); g0[gcap] = True
        s0 = np.zeros(len(st), bool); s0[scap] = True
        m = depth_match(log_len, g0, s0, DEPTH_TOL)
        rows = {GOAL: sorted(np.flatnonzero(m & g0).tolist()),
                START: sorted(np.flatnonzero(m & s0).tolist())}

        sref = np.zeros(len(st), bool); sref[rows[START]] = True
        for name in EXTRA:
            pool = np.flatnonzero(st == name).tolist()
            cap = balanced(pool, donors, EXTRA_MAX, EXTRA_PER_DONOR, random.Random(f"{s}-{name}"))
            p0 = np.zeros(len(st), bool); p0[cap] = True
            mm = depth_match(log_len, p0, sref, DEPTH_TOL)
            rows[name] = sorted(np.flatnonzero(mm & p0).tolist())

        n_all = sum(len(v) for v in rows.values())
        if len({i for v in rows.values() for i in v}) != n_all:
            sys.exit("4 本の重心で細胞が重なっている")

        # **seed を引き継いでいるので、前実験と 1 行も違わないはず。違えば止める**
        ref = json.loads((PREV / f"state_rows-seed{s}.json").read_text())["rows"]
        if {k: list(v) for k, v in ref.items()} != {k: list(v) for k, v in rows.items()}:
            sys.exit(f"seed {s}: EXP-20260930-001 の state_rows-seed{s}.json と一致しない")

        (E / f"work/state_rows-seed{s}.json").write_text(json.dumps(
            {"seed": s, "rows": rows,
             "note": "EXP-20260930-001/work/state_rows-seed%d.json と同一内容"
                     "（このスクリプトが再現して照合した）" % s},
            ensure_ascii=False) + "\n")
        summary["per_seed"][str(s)] = {
            k: {"n_cells": len(v), "n_donors": len({donors[i] for i in v}),
                "median_log_len": float(np.median(log_len[v]))} for k, v in rows.items()}
        print(s, {k: len(v) for k, v in rows.items()}, "← EXP-20260930-001 と一致", flush=True)

    (E / "work/depth_match_summary_005.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=1) + "\n")
    print("→", E / "work/depth_match_summary_005.json")


if __name__ == "__main__":
    main()
