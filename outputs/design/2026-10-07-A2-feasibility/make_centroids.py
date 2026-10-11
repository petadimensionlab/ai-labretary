#!/usr/bin/env python3
"""A2 の状態の重心を作り、軸の長さの比を複数の参照状態の取り方で比べる。

- 細胞の抽出: 提供者ごとに上限 20 細胞（seed ごとに乱数の鍵を変える）
- 深さをそろえた版: トークン長の分位ビンごとに、比べる状態の最小個数に合わせて引き直す
- 出力: work/centroids-seed{seed}.npz（主の版）と ratios.json（全版の比較）
"""
import json, sys
from pathlib import Path
import numpy as np
from datasets import load_from_disk

ROOT = Path(__file__).resolve().parents[3]
W = ROOT / "outputs/design/2026-10-07-A2-feasibility/work"

def donor_cap(idx, donors, cap, rng):
    by = {}
    for i in idx:
        by.setdefault(donors[i], []).append(i)
    keep = []
    for d in sorted(by):
        v = sorted(by[d])
        if len(v) > cap:
            v = sorted(rng.choice(v, size=cap, replace=False).tolist())
        keep += v
    return sorted(keep)

def depth_match(groups, lengths, nbins=8, rng=None):
    """groups: {name: [row ids]} をトークン長の分位ビンで個数をそろえる。"""
    allrows = sorted({i for v in groups.values() for i in v})
    edges = np.quantile(lengths[allrows], np.linspace(0, 1, nbins + 1))
    edges[0], edges[-1] = -np.inf, np.inf
    out = {k: [] for k in groups}
    for b in range(nbins):
        pools = {k: [i for i in v if edges[b] <= lengths[i] < edges[b + 1]] for k, v in groups.items()}
        n = min(len(p) for p in pools.values())
        if n == 0:
            continue
        for k, p in pools.items():
            out[k] += sorted(rng.choice(sorted(p), size=n, replace=False).tolist()) if len(p) > n else sorted(p)
    return {k: sorted(v) for k, v in out.items()}

def ratios(cent, start, goal, ref, others=()):
    need = [start, goal, ref] + list(others)
    U = {n: cent[n] / np.linalg.norm(cent[n]) for n in need}
    dg = float(np.linalg.norm(U[start] - U[goal]))
    dr = float(np.linalg.norm(U[start] - U[ref]))
    def direction(n):
        rest = [U[m] for m in need if m != n]
        return U[n] - np.mean(rest, axis=0)
    lg = float(np.linalg.norm(direction(goal))); lr = float(np.linalg.norm(direction(ref)))
    cos_sg = float(np.dot(U[start], U[goal]))
    cos_sr = float(np.dot(U[start], U[ref]))
    return {"dist_start_to_goal": dg, "dist_start_to_reference": dr, "distance_ratio": dg / dr,
            "dirlen_goal": lg, "dirlen_reference": lr, "direction_ratio": lg / lr,
            "cos_start_goal": cos_sg, "angle_start_goal_deg": float(np.degrees(np.arccos(np.clip(cos_sg, -1, 1)))),
            "cos_start_reference": cos_sr,
            "angle_start_reference_deg": float(np.degrees(np.arccos(np.clip(cos_sr, -1, 1))))}

def main():
    z = np.load(W / "cell_embs.npz")
    X, rid = z["X"], z["row_id"].astype(int)
    pos = {int(r): i for i, r in enumerate(rid)}
    ds = load_from_disk(str(W / "tokenized.dataset"))
    st = np.array(ds["isp_state"]); don = np.array(ds["donor_id"]); ln = np.array(ds["length"])
    have = set(pos)
    states = {s: [i for i in np.where(st == s)[0].tolist() if i in have] for s in np.unique(st)}
    CAP = 20
    res = {"n_embedded": len(rid), "dim": int(X.shape[1]), "cap_per_donor": CAP, "per_seed": {}}
    for seed in (81, 82, 83):
        rng = np.random.default_rng(seed)
        capped = {s: donor_cap(v, don, CAP, rng) for s, v in states.items()}
        cent = {s: X[[pos[i] for i in v]].mean(axis=0) for s, v in capped.items()}
        r = {"n_cells": {s: len(v) for s, v in capped.items()},
             "len_median": {s: float(np.median(ln[v])) for s, v in capped.items()},
             "variants": {}}
        r["variants"]["A_dFIB_A1と同じ構成"] = ratios(cent, "FIB", "aStr", "dFIB")
        for ref in ("VSMC", "MC", "VSMCP"):
            r["variants"][f"C_参照={ref}"] = ratios(cent, "FIB", "aStr", ref)
        # 深さをそろえた版（FIB / aStr / dFIB の 3 状態で長さをそろえる）
        g = {k: capped[k] for k in ("FIB", "aStr", "dFIB")}
        dm = depth_match(g, ln, nbins=8, rng=np.random.default_rng(seed + 1000))
        cent_dm = {k: X[[pos[i] for i in v]].mean(axis=0) for k, v in dm.items()}
        r["depth_matched_n"] = {k: len(v) for k, v in dm.items()}
        r["depth_matched_len_median"] = {k: float(np.median(ln[v])) for k, v in dm.items()}
        r["variants"]["B_深さをそろえた_参照=dFIB"] = ratios(cent_dm, "FIB", "aStr", "dFIB")
        # 深さをそろえた版（参照 = VSMC/P）
        g2 = {k: capped[k] for k in ("FIB", "aStr", "VSMCP")}
        dm2 = depth_match(g2, ln, nbins=8, rng=np.random.default_rng(seed + 2000))
        cent_dm2 = {k: X[[pos[i] for i in v]].mean(axis=0) for k, v in dm2.items()}
        r["depth_matched2_n"] = {k: len(v) for k, v in dm2.items()}
        r["variants"]["D_深さをそろえた_参照=VSMCP"] = ratios(cent_dm2, "FIB", "aStr", "VSMCP")
        allc = dict(cent)
        allc.update({f"{k}_dm": v for k, v in cent_dm.items()})
        allc.update({f"{k}_dm2": v for k, v in cent_dm2.items()})
        np.savez(W / f"centroids-seed{seed}.npz", **allc)
        res["per_seed"][str(seed)] = r
    out = ROOT / "outputs/design/2026-10-07-A2-feasibility/centroid_ratios.json"
    out.write_text(json.dumps(res, ensure_ascii=False, indent=1) + "\n")
    names = sorted(res["per_seed"]["81"]["variants"])
    print(f"{'版':<30} {'seed81':>9} {'seed82':>9} {'seed83':>9} {'最小':>9}")
    for n in names:
        v = [res["per_seed"][s]["variants"][n]["distance_ratio"] for s in ("81", "82", "83")]
        print(f"{n:<30} {v[0]:9.4f} {v[1]:9.4f} {v[2]:9.4f} {min(v):9.4f}")
    print()
    for n in names:
        a = res["per_seed"]["81"]["variants"][n]
        print(f"{n}: 出発点→目標 {a['dist_start_to_goal']:.4f} / 出発点→参照 {a['dist_start_to_reference']:.4f} "
              f"/ 出発点と目標の角度 {a['angle_start_goal_deg']:.2f} 度 / 出発点と参照の角度 {a['angle_start_reference_deg']:.2f} 度")
    print(json.dumps({"depth_matched_n": res["per_seed"]["81"]["depth_matched_n"],
                      "depth_matched_len_median": res["per_seed"]["81"]["depth_matched_len_median"],
                      "n_cells": res["per_seed"]["81"]["n_cells"],
                      "len_median": res["per_seed"]["81"]["len_median"]}, ensure_ascii=False, indent=1))

if __name__ == "__main__":
    main()
