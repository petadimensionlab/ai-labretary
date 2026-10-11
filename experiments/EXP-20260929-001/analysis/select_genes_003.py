#!/usr/bin/env python3
"""EXP-20260929-001: 摂動する遺伝子・対照遺伝子・部分セットを決める（**結果を見る前に固定する**）。

使うのは (i) 仮説・文献で先に決めた遺伝子の一覧と、(ii) 出発点の細胞だけから計算できる
実行可能性の量（トークン化済み入力での検出細胞数・提供者数・入力中の位置・入力長）だけである。
**摂動の結果は一切使わない**（この時点で摂動を 1 回も実行していない）。

出力:
  work/gene_stats_003.csv        全トークンの統計（検出率・偏り・位置・長さ）
  work/genes_003.json            役割ごとの遺伝子（摂動する全遺伝子）
  work/null_pool_003.json        対照遺伝子 500（先頭 300 が主判定。入れ子）
  work/subsets_003.json          部分セット（結果を見る前に固定）
  work/precheck_coverage_003.json 遺伝子ごとの抽出細胞数・提供者数（3 seed）
"""
from __future__ import annotations

import json
import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path("/Users/petadimensionlab/workspace/research/ai-labretary")
E = ROOT / "experiments/EXP-20260929-001"
sys.path.insert(0, str(ROOT / "scripts/geneformer"))
import run_isp as ri  # noqa: E402

SEEDS = [61, 62, 63]
MIN_CELLS, MAX_CELLS, MAX_PER_DONOR = 120, 200, 20
GOAL, START, MYOF = "FIB_orgHigh", "FIB_orgLow", "MYOF_fib"

# 対照遺伝子の条件（EXP-20260928-002 と同じ形。数は出発点が 4 倍になったぶん増える）
NULL_ABS_BIAS_MAX = 0.5          # |bias_goal| < 0.5 かつ |bias_myof| < 0.5 かつ |bias_ecm| < 0.5
NULL_MIN_DETECT, NULL_MIN_DONORS = 250, 12
N_NULL_PRIMARY, N_NULL_TOTAL = 300, 500
N_RESERVE_PER_JUDGED = 2         # 判定遺伝子 1 個につき、近い対照を 2 個ずつ先に確保する
NULL_ORDER_SEED = 20260930
TOL_N = {"lncells": float(np.log(1.5)), "log2pos": 0.25, "log2len": 0.25}

# 共起オッズ比（HYP-20260929-01 §3.1。**目標の定義と摂動の対象の共起。結果を見る前に固定**）
OR_WITH_GOAL = {
    "CXCL13": 0.92, "CCL19": 6.29, "CCL21": 2.96, "LTB": 2.41, "TNFSF14": 5.34, "CXCL16": 1.10,
    "CXCL1": 3.20, "CXCL14": 1.56, "CCL2": 1.92, "CCL11": 1.23, "CCL28": 0.93, "CSF1": 1.40,
    "IL7": 1.51, "IL15": 1.13, "IL33": 2.11, "IL34": 1.52, "IL11": 0.75, "IL18": 1.34,
    "IL32": 1.71, "TNFSF10": 1.70, "TNFSF15": 0.70, "CXCL12": 1.33, "VEGFA": 1.05,
    "PGF": 1.80, "ANGPT1": 0.64}
OR_SPLIT = 2.0                   # HYP §4.2 が事前登録した境目。**変えない**
TLS_CHEMOKINE = ["CXCL13", "CCL19", "CCL21", "LTB", "TNFSF14"]   # HYP §5 P-4


def main() -> None:
    from datasets import load_from_disk

    groups = json.loads((E / "analysis/excluded_symbols_003.json").read_text())["groups"]
    sym2ens = json.loads((E / "work/symbol_to_ensembl_003.json").read_text())
    pc = json.loads((E / "work/precheck_counts_003.json").read_text())
    ecm7 = json.loads((E / "work/input_states.json").read_text())["ecm_program"]
    ecm_cand = json.loads((E / "work/input_states.json").read_text())["ecm_candidates_order"]
    tokd = ri.gene_token_dict()
    ens2sym: dict[str, str] = {}
    for s, e in sym2ens.items():
        ens2sym.setdefault(e, s)

    ds = load_from_disk(str(E / "work/tokenized.dataset"))
    st = np.array([str(x) for x in ds["fib_org_state"]])
    ecm = np.array([str(x) for x in ds["ecm_high"]]) == "yes"
    donors = np.array([str(x) for x in ds["donor_id"]])
    ids = [np.asarray(x) for x in ds["input_ids"]]
    lens = np.array([len(x) for x in ids], float)

    indptr = np.cumsum([0] + [len(x) for x in ids])
    col = np.concatenate(ids)
    ntok = max(int(col.max()) + 1, max(tokd.values()) + 1)
    posdata = np.concatenate([np.arange(1, len(x) + 1) for x in ids]).astype(np.int32)
    Mb = sp.csr_matrix((np.ones(len(col), bool), col, indptr), shape=(len(ids), ntok)).tocsc()
    Mp = sp.csr_matrix((posdata, col, indptr), shape=(len(ids), ntok)).tocsc()

    sets = {"start": st == START, "goal": st == GOAL, "ecm": ecm, "myof": st == MYOF}
    d = pd.DataFrame({"token": np.arange(ntok)})
    for k, v in sets.items():
        n = np.asarray(Mb[v].sum(0)).ravel()
        d[f"n_{k}"] = n
        d[f"det_{k}"] = n / v.sum()
    eps = 1.0 / int(sets["goal"].sum())
    for k in ("goal", "ecm", "myof"):
        d[f"bias_{k}"] = np.log2((d[f"det_{k}"] + eps) / (d.det_start + eps))

    # 出発点での提供者数、上限を当てた後の細胞数（上限後の「数」は seed に依らない）
    sidx = np.flatnonzero(sets["start"])
    Ms = Mb[sidx].tocsc()
    Mps = Mp[sidx].tocsc()
    dcodes = pd.Categorical(donors[sidx]).codes
    ndon = np.zeros(ntok, int)
    ncap = np.zeros(ntok, int)
    lp_med = np.full(ntok, np.nan)
    ll_med = np.full(ntok, np.nan)
    slen = lens[sidx]
    for t in range(ntok):
        r = Ms.indices[Ms.indptr[t]:Ms.indptr[t + 1]]
        if len(r) == 0:
            continue
        dd = dcodes[r]
        cnt = np.bincount(dd)
        ndon[t] = int((cnt > 0).sum())
        ncap[t] = int(min(MAX_CELLS, np.minimum(cnt, MAX_PER_DONOR).sum()))
        p = Mps.data[Mps.indptr[t]:Mps.indptr[t + 1]].astype(float)
        lp_med[t] = float(np.median(np.log2(p)))          # pos は 1 始まり → log2(1+index)
        ll_med[t] = float(np.median(np.log2(slen[r])))
    d["n_donor_start"], d["n_cap20"] = ndon, ncap
    d["lp_med"], d["ll_med"] = lp_med, ll_med
    tok2ens = {v: k for k, v in tokd.items()}
    d["ensembl"] = [tok2ens.get(int(t), "") for t in d.token]
    d["symbol"] = [ens2sym.get(e, "") for e in d.ensembl]
    d = d[(d.ensembl != "")].reset_index(drop=True)
    d.to_csv(E / "work/gene_stats_003.csv", index=False)

    ref = d.set_index("symbol")
    ref = ref[~ref.index.duplicated()]

    # ---------------------------------------------------------------- 役割ごとの遺伝子
    sec = pc["verdict"]["secretory_pass"]               # トークン化後の数え直しに合格した分泌プログラム
    org7 = pc["verdict"]["org_markers7_pass"]
    myof24 = pc["verdict"]["myof24_pass"]
    hk3, panel7, ref2 = groups["housekeeping3"], groups["negctrl_panel7"], groups["reference2"]
    judged = sec + org7 + myof24 + hk3 + panel7 + ref2
    assert len(judged) == len(set(judged)), "役割の間で遺伝子が重複している"

    # ---------------------------------------------------------------- 対照遺伝子
    blocked = set(groups["blinded5"]) | set(ecm_cand) | set(ecm7) | set(judged) \
        | set(groups["secretory25"]) | set(groups["org_markers7"])
    cand = d[(~d.symbol.isin(blocked)) & (d.symbol != "")
             & (d.n_start >= NULL_MIN_DETECT) & (d.n_donor_start >= NULL_MIN_DONORS)
             & (d.n_cap20 >= MIN_CELLS)
             & (d.bias_goal.abs() < NULL_ABS_BIAS_MAX)
             & (d.bias_myof.abs() < NULL_ABS_BIAS_MAX)
             & (d.bias_ecm.abs() < NULL_ABS_BIAS_MAX)].copy()
    cand = cand.set_index("symbol")
    cand = cand[~cand.index.duplicated()]

    order: list[str] = []
    taken: set[str] = set()
    for sym in judged:                                  # 判定遺伝子に近い対照を 2 個ずつ確保する
        r = ref.loc[sym]
        av = cand[~cand.index.isin(taken)]
        dist = pd.DataFrame({
            "a": (np.log(av.n_cap20) - np.log(r.n_cap20)).abs() / TOL_N["lncells"],
            "b": (av.lp_med - r.lp_med).abs() / TOL_N["log2pos"],
            "c": (av.ll_med - r.ll_med).abs() / TOL_N["log2len"]}).max(axis=1)
        for s2 in dist.nsmallest(N_RESERVE_PER_JUDGED).index:
            order.append(str(s2)); taken.add(str(s2))
    rest = [s for s in cand.index if s not in taken]
    random.Random(NULL_ORDER_SEED).shuffle(rest)
    order += rest
    if len(order) < N_NULL_TOTAL:
        sys.exit(f"対照遺伝子が足りない（{len(order)} < {N_NULL_TOTAL}）")
    null500 = order[:N_NULL_TOTAL]
    null300 = null500[:N_NULL_PRIMARY]

    # ---------------------------------------------------------------- 部分セット
    or_high = [g for g in sec if OR_WITH_GOAL[g] >= OR_SPLIT]
    or_low = [g for g in sec if OR_WITH_GOAL[g] < OR_SPLIT]
    tls = [g for g in sec if g in TLS_CHEMOKINE]
    subsets = {
        "secretory_all": sec,                  # 主判定
        "co_or_low": or_low,                   # HYP P-6。**循環でない証拠の要**
        "co_or_high": or_high,                 # 記述値（個数が少ないので判定に使わない場合がある）
        "tls_chemokine": tls,                  # HYP P-4（記述値）
        "non_tls": [g for g in sec if g not in TLS_CHEMOKINE],
        "org_markers7": org7,                  # HYP P-3。**健全性チェック。支持の証拠にしない**
        "myof_markers24": myof24,
        "housekeeping3": hk3,
        "negctrl_panel7": panel7,
    }

    roles = ({g: "target" for g in sec} | {g: "marker_org" for g in org7}
             | {g: "marker_myof" for g in myof24} | {g: "housekeeping" for g in hk3}
             | {g: "negctrl_panel" for g in panel7} | {g: "reference_myof" for g in ref2}
             | {g: "negative" for g in null500})
    allg = sec + org7 + myof24 + hk3 + panel7 + ref2 + null500
    assert len(allg) == len(set(allg)), "遺伝子が重複している"
    assert not (set(allg) & set(groups["blinded5"])), "盲検中の遺伝子が入っている"

    def rec(sym: str) -> dict:
        r = ref.loc[sym]
        return {"symbol": sym, "ensembl": sym2ens[sym], "role": roles[sym],
                "bias_goal": round(float(r.bias_goal), 3), "bias_ecm": round(float(r.bias_ecm), 3),
                "bias_myof": round(float(r.bias_myof), 3), "n_start": int(r.n_start),
                "n_donor_start": int(r.n_donor_start), "n_cap20": int(r.n_cap20),
                "lp_med": round(float(r.lp_med), 4), "ll_med": round(float(r.ll_med), 4)}

    (E / "work/genes_003.json").write_text(json.dumps(
        {"genes": [rec(s) for s in allg],
         "n_primary": len(sec) + len(org7) + len(myof24) + len(hk3) + len(panel7) + len(ref2)
         + N_NULL_PRIMARY,
         "n_extended": len(allg),
         "by_role": {r: sorted(g for g in allg if roles[g] == r) for r in set(roles.values())}},
        indent=1, ensure_ascii=False) + "\n")
    (E / "work/null_pool_003.json").write_text(json.dumps(
        {"order_seed": NULL_ORDER_SEED, "n_candidates": int(len(cand)),
         "reserved_per_judged": N_RESERVE_PER_JUDGED, "n_reserved": len(taken),
         "conditions": {"abs_bias_max": NULL_ABS_BIAS_MAX, "min_detect": NULL_MIN_DETECT,
                        "min_donors": NULL_MIN_DONORS, "min_cap20": MIN_CELLS},
         "versions": {"A_300": [rec(s) for s in null300], "B_500": [rec(s) for s in null500]},
         "nested": null500[:N_NULL_PRIMARY] == null300}, indent=1, ensure_ascii=False) + "\n")
    (E / "work/subsets_003.json").write_text(json.dumps(
        {"or_split": OR_SPLIT, "or_with_goal": OR_WITH_GOAL,
         "subsets": {k: [rec(s) for s in v] for k, v in subsets.items()},
         "sizes": {k: len(v) for k, v in subsets.items()}}, indent=1, ensure_ascii=False) + "\n")

    cov = {}
    for sym in allg:
        e = sym2ens[sym]
        t = tokd[e]
        idx = Ms.indices[Ms.indptr[t]:Ms.indptr[t + 1]].tolist()
        per = {}
        for s in SEEDS:
            keep = ri.balanced_subset(idx, list(donors[sidx]), MAX_CELLS, MAX_PER_DONOR,
                                      random.Random(f"{s}-{e}"))
            per[str(s)] = {"n_cells": len(keep), "n_donors": len({donors[sidx][i] for i in keep})}
        cov[sym] = {"ensembl": e, "per_seed": per,
                    "min_cells": min(v["n_cells"] for v in per.values())}
    (E / "work/precheck_coverage_003.json").write_text(
        json.dumps(cov, indent=1, ensure_ascii=False) + "\n")
    bad = [s for s, v in cov.items() if v["min_cells"] < MIN_CELLS]

    print(f"対照候補 {len(cand)} 個 → 対照 {N_NULL_PRIMARY}（主判定）/ {N_NULL_TOTAL}（入れ子の感度分析）")
    print(f"摂動する遺伝子: 分泌 {len(sec)} + 目標マーカー {len(org7)} + 筋線維芽細胞 {len(myof24)} "
          f"+ 恒常発現 {len(hk3)} + 中立パネル {len(panel7)} + 参照 {len(ref2)} + 対照 {N_NULL_PRIMARY} "
          f"= {len(sec)+len(org7)+len(myof24)+len(hk3)+len(panel7)+len(ref2)+N_NULL_PRIMARY}")
    print("部分セット:", json.dumps({k: len(v) for k, v in subsets.items()}, ensure_ascii=False))
    print("120 細胞を下回る遺伝子:", bad if bad else "なし")
    if bad:
        sys.exit("抽出細胞数が下限を下回る遺伝子がある。設計を見直すこと")


if __name__ == "__main__":
    main()
