#!/usr/bin/env python3
"""EXP-20260929-001: 状態の重心に使う細胞（dataset の行番号）を seed ごとに決める。

guardrails 2026-09-27「状態の重心は、深さ（読み取りの深さ）をそろえた後の細胞から作る」。
深さのそろえ方は EXP-20260927-002 / 003 / EXP-20260928-001 / 002 と同じ
（log(トークン長)の 1:1 最近傍マッチ、許容 0.10）。

作る重心は 4 本:
  FIB_orgHigh  目標              … run_isp.py が作る（rows_file の "rows"）
  FIB_orgLow   出発点            … 同上。**目標とマッチした相手の細胞**
  FIB_ecmHigh  対照の状態（第 3 の軸） … extra_centroids_003.py が作る（"extra_rows"）
  MYOF_fib     筋線維芽細胞（第 2 の軸） … 同上

対照の状態 FIB_ecmHigh は**状態ラベルではなく obs の ecm_high 列**で指定する。
こうすると出発点の細胞プール（18,099）が削られず、CXCL13 の余裕（127 細胞）を減らさない。
ただし**出発点の重心を作る細胞からは ecm_high を外す**ので、重心どうしの細胞は重ならない。

出力: work/state_rows-seed{S}.json, work/depth_match_summary_003.json
"""
from __future__ import annotations

import json
import random
import sys
from pathlib import Path

import numpy as np

ROOT = Path("/Users/petadimensionlab/workspace/research/ai-labretary")
E = ROOT / "experiments/EXP-20260929-001"
sys.path.insert(0, str(ROOT / "scripts/geneformer"))

SEEDS = [61, 62, 63]
DEPTH_TOL = 0.10          # |Δlog(トークン長)| ≤ 0.10（EXP-20260927-002 以来 同じ）
GOAL_MAX, GOAL_PER_DONOR = 950, 60
START_POOL_MAX = 3000     # マッチの相手を広く取る（1:1 マッチの質を上げる。採るのはマッチした分だけ）
EXTRA_MAX, EXTRA_PER_DONOR = 900, 60

GOAL, START = "FIB_orgHigh", "FIB_orgLow"


def balanced(idx, donors, n_max, per_donor, rng):
    """run_isp.balanced_subset と同じ規則。"""
    order = list(idx)
    rng.shuffle(order)
    used: dict[str, int] = {}
    keep = []
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
    """log(トークン長)の最近傍 1:1 マッチ（置換なし）。EXP-20260927-002 と同じ実装。"""
    pi, ni = np.where(pos)[0], np.where(neg)[0]
    ni = ni[np.argsort(log_len[ni], kind="stable")]
    nv = log_len[ni]
    used = np.zeros(len(ni), bool)
    keep: list[int] = []
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
    states = np.array([str(x) for x in ds["fib_org_state"]])
    ecm = np.array([str(x) for x in ds["ecm_high"]]) == "yes"
    donors = [str(x) for x in ds["donor_id"]]
    length = np.array(ds["length"], dtype=float)
    log_len = np.log(length)

    goal_idx = np.flatnonzero(states == GOAL).tolist()
    start_idx = np.flatnonzero((states == START) & ~ecm).tolist()   # 重心用（ecm_high を外す）
    ecm_idx = np.flatnonzero(ecm).tolist()
    myof_idx = np.flatnonzero(states == "MYOF_fib").tolist()

    summary = {"seeds": SEEDS, "depth_tol": DEPTH_TOL,
               "pool": {"goal": len(goal_idx), "start_excl_ecm": len(start_idx),
                        "ecm_high": len(ecm_idx), "myof": len(myof_idx)}, "per_seed": {}}

    for s in SEEDS:
        rng = random.Random(s)
        gcap = balanced(goal_idx, donors, GOAL_MAX, GOAL_PER_DONOR, rng)
        scap = balanced(start_idx, donors, START_POOL_MAX, GOAL_PER_DONOR, rng)
        g0 = np.zeros(len(states), bool); g0[gcap] = True
        s0 = np.zeros(len(states), bool); s0[scap] = True
        m = depth_match(log_len, g0, s0, DEPTH_TOL)
        rows = {GOAL: sorted(np.flatnonzero(m & g0).tolist()),
                START: sorted(np.flatnonzero(m & s0).tolist())}

        extra = {}
        sref = np.zeros(len(states), bool); sref[rows[START]] = True
        for name, pool in [("FIB_ecmHigh", ecm_idx), ("MYOF_fib", myof_idx)]:
            cap = balanced(pool, donors, EXTRA_MAX, EXTRA_PER_DONOR, random.Random(f"{s}-{name}"))
            p0 = np.zeros(len(states), bool); p0[cap] = True
            mm = depth_match(log_len, p0, sref, DEPTH_TOL)
            extra[name] = sorted(np.flatnonzero(mm & p0).tolist())

        allrows = {**rows, **extra}
        assert len({i for v in allrows.values() for i in v}) == sum(len(v) for v in allrows.values()), \
            "4 本の重心で細胞が重なっている"
        (E / f"work/state_rows-seed{s}.json").write_text(json.dumps(
            {"seed": s, "rows": rows, "extra_rows": extra,
             "unmatched_goal_rows": sorted(gcap)}, ensure_ascii=False) + "\n")
        summary["per_seed"][str(s)] = {
            n: {"n_cells": len(v), "n_donors": len({donors[i] for i in v}),
                "median_length": float(np.median(length[v])),
                "median_n_genes": float(np.median(np.array(ds["n_genes"], dtype=float)[v]))}
            for n, v in allrows.items()}
        summary["per_seed"][str(s)]["unmatched_goal"] = {"n_cells": len(gcap)}

    (E / "work/depth_match_summary_003.json").write_text(
        json.dumps(summary, indent=1, ensure_ascii=False) + "\n")
    print(json.dumps(summary, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
