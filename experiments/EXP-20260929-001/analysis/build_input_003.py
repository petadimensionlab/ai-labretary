#!/usr/bin/env python3
"""EXP-20260929-001: トークン化の入力 h5ad を作る（data/raw は読み取りのみ）。

HYP-20260929-01 §2.1 の定義をそのまま実装する。
  - 線維芽細胞 = subclass.l2 ∈ {FIB, M-FIB, aFIB, dFIB, dM-FIB}（筋線維芽細胞は含めない）
  - マーカー得点 = オーガナイザー様マーカー 7 個の log1p(CP10K) の和
  - 目標 FIB_orgHigh = 検出遺伝子数の十分位の中で得点 上位 5%
  - 出発点 FIB_orgLow = 残りの線維芽細胞すべて
  - 筋線維芽細胞 MYOF_fib = subclass.l2 ∈ {MYOF, cycMYOF}（第 2 の軸。摂動しない）
  - 対照の状態 ecm_high = 出発点の中で、同じ作り方（十分位内 上位 5%）を
    細胞外基質プログラム 7 個に当てたもの（第 3 の軸。**状態ラベルにはしない**ので
    摂動する細胞プール（出発点 18,099）は削られない）

十分位と上位の取り方は outputs/literature/2026-09-29-goal-redefinition/adopted.py の
strat_goal と完全に同一（目標 950 細胞 / 18,099 細胞を再現することを検算する）。

出力: work/input.h5ad, work/input_states.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path("/Users/petadimensionlab/workspace/research/ai-labretary")
RAW = ROOT / "data/raw/kidney_tls/A1_KPMP_v1.5_snRNA_Lake2023.h5ad"
E = ROOT / "experiments/EXP-20260929-001"

FIB_LABELS = ["FIB", "M-FIB", "aFIB", "dFIB", "dM-FIB"]
MYOF_LABELS = ["MYOF", "cycMYOF"]

# 目標の定義に使うマーカー 7 個（HYP-20260929-01 §2.1。CXCL13/CCL19/CCL21 を含まない）
ORG7 = ["PDPN", "NGFR", "ICAM1", "VCAM1", "TNFSF13B", "MADCAM1", "ENPP2"]

# 対照の状態（細胞外基質プログラム）の候補。**結果を見る前にこの順番で固定した**。
# 先頭から、除外集合（目標マーカー 7 / 摂動 25 / 筋線維芽細胞マーカー 24 / 恒常発現 3 /
# 中立対照パネル 7 / 参照 2 / 盲検 5）に当たらず、A1 に存在するものを 7 個採る。
ECM_CANDIDATES = ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A3", "FBN1", "SPARC",
                  "LUM", "DCN", "FN1", "POSTN", "TNC", "ELN", "MMP2", "TIMP1"]

TOP_FRAC = 0.05
N_BINS = 10


def rd(h, c):
    g = h["obs"][c]
    if isinstance(g, h5py.Group):
        cats = np.array(g["categories"].asstr()[...])
        codes = g["codes"][...]
        return np.where(codes >= 0, cats[np.clip(codes, 0, None)], "nan")
    v = g[...]
    return np.array(g.asstr()[...]) if v.dtype.kind in "OS" else v


def strat_top(n_genes: np.ndarray, score: np.ndarray, top_frac: float, nbin: int) -> np.ndarray:
    """adopted.py の strat_goal と同一。検出遺伝子数の十分位の中で得点 上位 top_frac を True にする。"""
    dq = pd.qcut(pd.Series(n_genes).rank(method="first"), nbin, labels=False).values
    g = np.zeros(len(n_genes), bool)
    for b in range(nbin):
        idx = np.where(dq == b)[0]
        k = max(1, int(round(top_frac * len(idx))))
        g[idx[np.lexsort((idx, -score[idx]))[:k]]] = True
    return g


def main() -> None:
    import anndata as ad

    excluded = set(json.loads((E / "analysis/excluded_symbols_003.json").read_text())["symbols"])

    with h5py.File(RAW, "r") as h:
        ct = rd(h, "subclass.l2")
        don = rd(h, "donor_id")
        dis = rd(h, "disease")
        reg = rd(h, "region")
        n = len(ct)
        fib = np.isin(ct, FIB_LABELS)
        myof = np.isin(ct, MYOF_LABELS)
        keep_mask = fib | myof

        v = h["raw/var"]
        ens = np.array(v["_index"].asstr()[...])
        fn = v["feature_name"]
        sym = (np.array(fn["categories"].asstr()[...])[fn["codes"][...]]
               if isinstance(fn, h5py.Group) else np.array(fn.asstr()[...]))
        s2i: dict[str, int] = {}
        for i, s in enumerate(sym):
            s2i.setdefault(s, i)

        ecm7 = [g for g in ECM_CANDIDATES if g in s2i and g not in excluded][:7]
        if len(ecm7) != 7:
            sys.exit(f"細胞外基質プログラムが 7 個そろわない: {ecm7}")
        miss = [g for g in ORG7 if g not in s2i]
        if miss:
            sys.exit(f"目標マーカーが A1 に無い: {miss}")

        Xg = h["raw/X"]
        indptr = Xg["indptr"][...]
        chk = Xg["data"][indptr[0]:indptr[min(200, n)]][:100000]
        counts_integer_like = bool(np.allclose(chk, np.round(chk)))
        if not counts_integer_like:
            sys.exit("raw/X が整数カウントではありません（guardrails 2026-09-25）")

        blocks, rows = [], []
        step = 20000
        for s in range(0, n, step):
            e = min(n, s + step)
            if not keep_mask[s:e].any():
                continue
            a, b = indptr[s], indptr[e]
            m = sp.csr_matrix((Xg["data"][a:b], Xg["indices"][a:b], indptr[s:e + 1] - a),
                              shape=(e - s, len(ens)))
            loc = np.where(keep_mask[s:e])[0]
            blocks.append(m[loc])
            rows.append(loc + s)
            print(f"  {e}/{n}", file=sys.stderr, end="\r")
        rows = np.concatenate(rows)
        X = sp.vstack(blocks, format="csr").astype(np.float32)

    n_counts = np.asarray(X.sum(1)).ravel()
    n_genes = np.diff(X.indptr)
    is_fib = fib[rows]
    is_myof = myof[rows]
    print(f"\n抽出 {X.shape[0]} 細胞（線維芽細胞 {int(is_fib.sum())} / 筋線維芽細胞 {int(is_myof.sum())}）",
          file=sys.stderr)

    def score(symbols: np.ndarray | list[str], sub: np.ndarray) -> np.ndarray:
        cols = [s2i[g] for g in symbols]
        raw = X[sub][:, cols].toarray().astype(np.float64)
        cp = raw / n_counts[sub][:, None] * 1e4
        return np.log1p(cp).sum(1)

    fi = np.flatnonzero(is_fib)
    goal_local = strat_top(n_genes[fi], score(ORG7, fi), TOP_FRAC, N_BINS)
    state = np.array(["MYOF_fib"] * X.shape[0], dtype=object)
    state[fi] = np.where(goal_local, "FIB_orgHigh", "FIB_orgLow")

    si = fi[~goal_local]                               # 出発点の細胞（A1 の行順のまま）
    ecm_local = strat_top(n_genes[si], score(ecm7, si), TOP_FRAC, N_BINS)
    ecm_high = np.zeros(X.shape[0], bool)
    ecm_high[si[ecm_local]] = True

    obs = pd.DataFrame({
        "donor_id": don[rows].astype(str), "subclass_l2": ct[rows].astype(str),
        "disease": dis[rows].astype(str), "region": reg[rows].astype(str),
        "fib_org_state": state.astype(str), "ecm_high": np.where(ecm_high, "yes", "no"),
        "n_counts": n_counts.astype(np.float32), "n_genes": n_genes.astype(np.int32),
        "a1_row": rows.astype(np.int64)},
        index=[f"c{int(r)}" for r in rows])

    adata = ad.AnnData(X=X, obs=obs)
    adata.var_names = ens.astype(str)
    adata.var["ensembl_id"] = ens.astype(str)
    adata.var["gene_name"] = sym.astype(str)
    out = E / "work/input.h5ad"
    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out)

    def cnt(m):
        return {"n_cells": int(m.sum()), "n_donors": int(pd.Series(don[rows][m]).nunique())}

    info = {
        "source": str(RAW.relative_to(ROOT)), "counts_integer_like": counts_integer_like,
        "shape": list(adata.shape), "org_markers": ORG7,
        "ecm_program": ecm7, "ecm_candidates_order": ECM_CANDIDATES,
        "top_frac": TOP_FRAC, "n_depth_bins": N_BINS,
        "states": {s: cnt(state == s) for s in ["FIB_orgHigh", "FIB_orgLow", "MYOF_fib"]},
        "ecm_high": cnt(ecm_high),
        "ecm_high_overlap_goal": int((ecm_high & (state == "FIB_orgHigh")).sum()),
        "median_n_genes": {s: float(np.median(n_genes[state == s]))
                           for s in ["FIB_orgHigh", "FIB_orgLow", "MYOF_fib"]},
        "median_n_counts": {s: float(np.median(n_counts[state == s]))
                            for s in ["FIB_orgHigh", "FIB_orgLow", "MYOF_fib"]},
        "goal_by_subclass": pd.Series(ct[rows][state == "FIB_orgHigh"]).value_counts().to_dict(),
        "goal_by_disease": pd.Series(dis[rows][state == "FIB_orgHigh"]).value_counts().to_dict(),
        "expected": {"FIB_orgHigh": 950, "FIB_orgLow": 18099, "goal_donors": 70, "start_donors": 79},
    }
    (E / "work/input_states.json").write_text(json.dumps(info, indent=1, ensure_ascii=False) + "\n")
    print(json.dumps(info, ensure_ascii=False, indent=1))

    ok = (info["states"]["FIB_orgHigh"]["n_cells"] == 950
          and info["states"]["FIB_orgLow"]["n_cells"] == 18099
          and info["states"]["FIB_orgHigh"]["n_donors"] == 70
          and info["states"]["FIB_orgLow"]["n_donors"] == 79)
    print("検算（HYP §2.2 の実測と一致）:", "OK" if ok else "不一致", file=sys.stderr)
    if not ok:
        sys.exit("目標・出発点の細胞数が HYP-20260929-01 §2.2 と一致しません")


if __name__ == "__main__":
    main()
