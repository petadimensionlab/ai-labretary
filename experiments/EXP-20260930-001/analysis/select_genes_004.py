#!/usr/bin/env python3
"""EXP-20260930-001 — 摂動する遺伝子と対照遺伝子群を決める（**結果を見る前に固定する**）。

使うのは (i) 仮説・文献であらかじめ決まっている遺伝子の一覧と、(ii) 出発点の細胞だけから
計算できる実行可能性の量（トークン化済み入力での検出細胞数・提供者数・入力中の位置・入力長）、
(iii) 目標と出発点での検出率の比（発現の偏り）だけである。
**摂動の結果は 1 つも使わない**（この時点で摂動を 1 回も実行していない）。

## 対照遺伝子群を引き直す理由（HYP-20260930-01 §5.3 の申し送り）

EXP-20260929-001 の対照 500 個をそのまま流用すると **COL5A2 が陽性対照 39 個と重なる**。
同じ遺伝子を陽性対照と対照の両方に置くことはできない。そこで**除外集合を作り直してから引き直す**。
除外するのは次のすべてである（1 個でも重なれば停止する）。

  - 本実験で判定に使う遺伝子（陽性対照 39・感度分析の 71・摂動対象 T1/T2 の候補すべて・
    恒常発現・中立対照パネル）
  - Lake ら 2023 が線維芽細胞系列の**どれかの状態のマーカー**として挙げた遺伝子すべて
    （目標側だけでなく出発点側・退化側も。対照に状態マーカーを入れないため）
  - 盲検中の 5 遺伝子
  - EXP-20260929-001 で目標・対照の状態の定義に使った遺伝子（分泌 25・目標マーカー 7・
    筋線維芽細胞マーカー 24・細胞外基質の候補 15）

出力:
  work/gene_stats_004.csv          全トークンの統計（検出率・偏り・位置・長さ・細胞数）
  work/genes_004.json              役割ごとの遺伝子（摂動する全遺伝子）
  work/null_pool_004.json          対照遺伝子 500（先頭 300 が入れ子の感度分析）
  work/subsets_004.json            判定に使うセット（結果を見る前に固定）
  work/overlap_check_004.json      **重複 0 件の確認**（陽性対照 × 対照など総当たり）
  work/precheck_coverage_004.json  遺伝子ごとの抽出細胞数・提供者数（3 seed）
"""
from __future__ import annotations

import json
import random
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path("/Users/petadimensionlab/workspace/research/ai-labretary")
E = ROOT / "experiments/EXP-20260930-001"
PREV = ROOT / "experiments/EXP-20260929-001"
LIT = ROOT / "outputs/literature/2026-09-30-independent-goal-state/markers"
sys.path.insert(0, str(ROOT / "scripts/geneformer"))
import run_isp as ri  # noqa: E402

SEEDS = [81, 82, 83]
MIN_CELLS, MAX_CELLS, MAX_PER_DONOR = 120, 200, 20
GOAL, START = "aStr", "FIB_MFIB"

# 対照遺伝子の条件（EXP-20260929-001 と同じ形。軸が 1 本になったので偏りの条件も 1 本）
NULL_ABS_BIAS_MAX = 0.5          # |bias_goal| < 0.5
NULL_MIN_DETECT, NULL_MIN_DONORS = 250, 12
N_NULL_PRIMARY, N_NULL_TOTAL = 300, 500
N_RESERVE_PER_JUDGED = 2         # 判定する遺伝子 1 個につき、性質の近い対照を 2 個ずつ先に確保する
NULL_ORDER_SEED = 20260930       # **並べ替えの seed。結果を見る前に固定**
TOL_N = {"lncells": float(np.log(1.5)), "log2pos": 0.25, "log2len": 0.25}


def read_marker_csv(p: Path) -> list[str]:
    return [x.strip() for x in p.read_text().splitlines()[1:] if x.strip()]


def main() -> None:
    from datasets import load_from_disk

    G = json.loads((E / "analysis/gene_sets_004.json").read_text())
    sym2ens = json.loads((E / "work/symbol_to_ensembl_004.json").read_text())
    tokd = ri.gene_token_dict()
    ens2sym: dict[str, str] = {}
    for s, e in sym2ens.items():
        ens2sym.setdefault(e, s)

    # ---------------------------------------------------------------- 文献との照合
    lake = json.loads((LIT / "lake2023_fibroblast_state_markers.json").read_text())
    mini42 = read_marker_csv(LIT / "adaptive_stromal_state_aFIB_plus_MYOF_minimal.csv")
    full76 = read_marker_csv(LIT / "adaptive_stromal_state_aFIB_plus_MYOF_full_conserved.csv")
    if not set(mini42) <= set(full76):
        sys.exit("最小 42 が保守 76 の部分集合になっていない（入れ子の感度分析が成り立たない）")

    def usable(g: str) -> bool:
        return g in sym2ens and sym2ens[g] in tokd

    pc39 = [g for g in mini42 if usable(g)]
    pc71 = [g for g in full76 if usable(g)]
    if sorted(pc39) != sorted(G["positive_control_39"]):
        sys.exit(f"陽性対照 39 が仮説と一致しない: {sorted(set(pc39) ^ set(G['positive_control_39']))}")
    if len(pc71) != 71:
        sys.exit(f"保守版の陽性対照が 71 個にならない（{len(pc71)}）")
    if not set(pc39) <= set(pc71):
        sys.exit("39 が 71 の部分集合になっていない")
    pc_extra32 = [g for g in pc71 if g not in set(pc39)]
    # **2026-09-30 の PI 判断**: 主判定の陽性対照を保守 76 個版（使える 71 個）にする。
    # 39 個版は 71 個版の部分集合なので、入れ子の感度分析として同じ 1 回の実行で評価できる。
    PC_PRIMARY = pc71

    # ---------------------------------------------------------------- 入力から数える
    ds = load_from_disk(str(E / "work/tokenized.dataset"))
    st = np.array([str(x) for x in ds["fib_anno_state"]])
    donors = np.array([str(x) for x in ds["donor_id"]])
    ids = [np.asarray(x) for x in ds["input_ids"]]
    lens = np.array([len(x) for x in ids], float)

    indptr = np.cumsum([0] + [len(x) for x in ids])
    col = np.concatenate(ids)
    ntok = max(int(col.max()) + 1, max(tokd.values()) + 1)
    posdata = np.concatenate([np.arange(1, len(x) + 1) for x in ids]).astype(np.int32)
    Mb = sp.csr_matrix((np.ones(len(col), bool), col, indptr), shape=(len(ids), ntok)).tocsc()
    Mp = sp.csr_matrix((posdata, col, indptr), shape=(len(ids), ntok)).tocsc()

    sets = {"start": st == START, "goal": st == GOAL, "dfib": st == "dFIB", "dmfib": st == "dMFIB"}
    d = pd.DataFrame({"token": np.arange(ntok)})
    for k, v in sets.items():
        n = np.asarray(Mb[v].sum(0)).ravel()
        d[f"n_{k}"] = n
        d[f"det_{k}"] = n / v.sum()
    eps = 1.0 / int(sets["goal"].sum())
    for k in ("goal", "dfib", "dmfib"):
        d[f"bias_{k}"] = np.log2((d[f"det_{k}"] + eps) / (d.det_start + eps))

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
        cnt = np.bincount(dcodes[r])
        ndon[t] = int((cnt > 0).sum())
        ncap[t] = int(min(MAX_CELLS, np.minimum(cnt, MAX_PER_DONOR).sum()))
        p = Mps.data[Mps.indptr[t]:Mps.indptr[t + 1]].astype(float)
        lp_med[t] = float(np.median(np.log2(p)))
        ll_med[t] = float(np.median(np.log2(slen[r])))
    d["n_donor_start"], d["n_cap20"] = ndon, ncap
    d["lp_med"], d["ll_med"] = lp_med, ll_med
    tok2ens = {v: k for k, v in tokd.items()}
    d["ensembl"] = [tok2ens.get(int(t), "") for t in d.token]
    d["symbol"] = [ens2sym.get(e, "") for e in d.ensembl]
    d = d[d.ensembl != ""].reset_index(drop=True)
    d.to_csv(E / "work/gene_stats_004.csv", index=False)

    ref = d.set_index("symbol")
    ref = ref[~ref.index.duplicated()]

    # ---------------------------------------------------------------- 摂動できるかの確認
    def cap20(g: str) -> int:
        return int(ref.loc[g].n_cap20) if g in ref.index else -1

    for name in ("perturb_T1_11", "perturb_T2_7", "positive_control_39", "housekeeping_3",
                 "negctrl_panel_7"):
        bad = [g for g in G[name] if cap20(g) < MIN_CELLS]
        if bad:
            sys.exit(f"{name}: 提供者上限 20 の後に {MIN_CELLS} 細胞に届かない遺伝子がある {bad}")
    dropped_T1 = [g for g in G["perturb_T1_candidates_16"] if g not in set(G["perturb_T1_11"])]
    dropped_T2 = [g for g in G["perturb_T2_candidates_12"] if g not in set(G["perturb_T2_7"])]
    for g in dropped_T1 + dropped_T2:
        if cap20(g) >= MIN_CELLS:
            sys.exit(f"{g} は落としたはずなのに下限を満たしている（仮説 §4 と食い違う）")
    pc71_bad = [g for g in pc71 if cap20(g) < MIN_CELLS]
    if pc71_bad:
        sys.exit(f"保守 71 個版の陽性対照のうち {len(pc71_bad)} 個が下限 {MIN_CELLS} 細胞に届かない: "
                 f"{ {g: cap20(g) for g in pc71_bad} }。主判定に使えないので設計を見直すこと")

    # ---------------------------------------------------------------- 除外集合
    prev_groups = json.loads((PREV / "analysis/excluded_symbols_003.json").read_text())["groups"]
    prev_states = json.loads((PREV / "work/input_states.json").read_text())
    lake_all = sorted({g for v in lake["states"].values() for g in v["genes"]})
    judged = (G["perturb_T1_11"] + G["perturb_T2_7"] + pc71
              + G["housekeeping_3"] + G["housekeeping_extra_descriptive"] + G["negctrl_panel_7"])
    if len(judged) != len(set(judged)):
        dup = sorted({g for g in judged if judged.count(g) > 1})
        sys.exit(f"役割の間で遺伝子が重複している: {dup}")

    excluded = (set(judged) | set(lake_all) | set(mini42) | set(full76)
                | set(G["perturb_T1_candidates_16"]) | set(G["perturb_T2_candidates_12"])
                | set(G["blinded_5"]) | set(G["reference_dFIB_markers_6_descriptive"])
                | set(prev_groups["secretory25"]) | set(prev_groups["org_markers7"])
                | set(prev_groups["myof_markers24"]) | set(prev_groups["housekeeping3"])
                | set(prev_groups["negctrl_panel7"]) | set(prev_groups["reference2"])
                | set(prev_states["ecm_candidates_order"]) | set(prev_states["ecm_program"]))

    # ---------------------------------------------------------------- 対照遺伝子を引き直す
    cand = d[(~d.symbol.isin(excluded)) & (d.symbol != "")
             & (d.n_start >= NULL_MIN_DETECT) & (d.n_donor_start >= NULL_MIN_DONORS)
             & (d.n_cap20 >= MIN_CELLS)
             & (d.bias_goal.abs() < NULL_ABS_BIAS_MAX)].copy()
    cand = cand.set_index("symbol")
    cand = cand[~cand.index.duplicated()]

    order: list[str] = []
    taken: set[str] = set()
    for s in judged:                       # 判定する遺伝子に性質の近い対照を 2 個ずつ確保する
        r = ref.loc[s]
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

    # ---------------------------------------------------------------- 重複の総当たり確認
    named = {"陽性対照 71（主判定）": PC_PRIMARY, "（内訳）最小 42 版の 39": pc39,
             "（内訳）保守版で加わる 32": pc_extra32,
             "摂動対象 T1 11": G["perturb_T1_11"], "摂動対象 T2 7": G["perturb_T2_7"],
             "恒常発現 3": G["housekeeping_3"],
             "恒常発現の追加 1（記述のみ）": G["housekeeping_extra_descriptive"],
             "中立対照パネル 7": G["negctrl_panel_7"], "対照遺伝子 500": null500}
    NESTED = {("陽性対照 71（主判定）", "（内訳）最小 42 版の 39"),
              ("陽性対照 71（主判定）", "（内訳）保守版で加わる 32")}   # 入れ子なので重なって当然
    pairs = {f"{a} × {b}": sorted(set(named[a]) & set(named[b]))
             for a, b in combinations(named, 2) if (a, b) not in NESTED}
    blinded_hit = sorted(set(G["blinded_still_closed_3"])
                         & {g for v in named.values() for g in v})
    overlap = {"各セットの大きさ": {k: len(v) for k, v in named.items()},
               "総当たりの重複": pairs,
               "重複の合計": int(sum(len(v) for v in pairs.values())),
               "盲検中の 3 遺伝子の混入": blinded_hit,
               "注": "陽性対照 71 は 39 と 32 の和なので、その 3 者の間の重なりだけは検査から外す",
               "EXP-20260929-001 の対照 500 をそのまま使った場合の陽性対照との重複": sorted(
                   {x["symbol"] for x in json.loads(
                       (PREV / "work/null_pool_003.json").read_text())["versions"]["B_500"]}
                   & set(pc71)),
               "引き直した対照 500 と EXP-20260929-001 の対照 500 の共通": len(
                   {x["symbol"] for x in json.loads(
                       (PREV / "work/null_pool_003.json").read_text())["versions"]["B_500"]}
                   & set(null500))}
    (E / "work/overlap_check_004.json").write_text(
        json.dumps(overlap, indent=1, ensure_ascii=False) + "\n")
    if overlap["重複の合計"] or blinded_hit:
        sys.exit(f"重複がある: {overlap}")

    # ---------------------------------------------------------------- セットと役割
    subsets = {
        "T1": G["perturb_T1_11"],                    # **主判定 P-1**
        "T1_no_CXCL12": G["perturb_T1_10_no_CXCL12"],  # 入れ子の感度分析
        "T2": G["perturb_T2_7"],                     # P-2
        "T1_T2": G["perturb_T1_11"] + G["perturb_T2_7"],   # P-3
        # **P-4 健全性チェックの主判定は保守 71 個版**（PI 判断 2026-09-30）
        "positive_control71": pc71,
        # 以下はすべて 71 の部分集合＝入れ子の感度分析
        "positive_control39": pc39,                  # 最小 42 版（PI 判断の前の主判定）
        "positive_control65": [g for g in pc71
                               if g not in set(G["positive_control71_opposite_6"])],
        "positive_control70": [g for g in pc71 if g != "LTBP1"],
        "positive_control35": [g for g in pc39
                               if g not in set(G["positive_control_opposite_4"])],
        # LTBP1 は入力中の位置が上位すぎて比べられる対照が 1〜2 個しかない（事前チェックで判明）。
        # **セット平均には残す**が、これを除いた版での再判定を入れ子の感度分析として事前登録する。
        "positive_control38": [g for g in pc39 if g != "LTBP1"],
        "housekeeping": G["housekeeping_3"],         # P-5
        "negctrl_panel": G["negctrl_panel_7"],       # P-5
    }
    if len(subsets["T1_T2"]) != 18:
        sys.exit("T1 ∪ T2 が 18 個にならない")

    roles = ({g: "target_T1" for g in G["perturb_T1_11"]}
             | {g: "target_T2" for g in G["perturb_T2_7"]}
             | {g: "positive_control" for g in PC_PRIMARY}
             | {g: "housekeeping" for g in G["housekeeping_3"]}
             | {g: "housekeeping_extra" for g in G["housekeeping_extra_descriptive"]}
             | {g: "negctrl_panel" for g in G["negctrl_panel_7"]}
             | {g: "negative" for g in null500})
    allg = judged + null500
    if len(allg) != len(set(allg)):
        sys.exit("摂動する遺伝子が重複している")
    if set(allg) & set(G["blinded_still_closed_3"]):
        sys.exit("盲検中の遺伝子が入っている")

    def rec(s: str) -> dict:
        r = ref.loc[s]
        return {"symbol": s, "ensembl": sym2ens[s], "role": roles[s],
                "bias_goal": round(float(r.bias_goal), 3),
                "bias_dfib": round(float(r.bias_dfib), 3),
                "n_start": int(r.n_start), "n_goal": int(r.n_goal),
                "det_start": round(float(r.det_start), 4), "det_goal": round(float(r.det_goal), 4),
                "n_donor_start": int(r.n_donor_start), "n_cap20": int(r.n_cap20),
                "lp_med": round(float(r.lp_med), 4), "ll_med": round(float(r.ll_med), 4)}

    (E / "work/genes_004.json").write_text(json.dumps(
        {"genes": [rec(s) for s in allg], "n_total": len(allg),
         "by_role": {r: sorted(g for g in allg if roles[g] == r) for r in sorted(set(roles.values()))},
         "dropped": {"T1": {g: cap20(g) for g in dropped_T1},
                     "T2": {g: cap20(g) for g in dropped_T2},
                     "主判定の陽性対照 71 個で下限を割る": {g: cap20(g) for g in pc71_bad}}},
        indent=1, ensure_ascii=False) + "\n")
    (E / "work/null_pool_004.json").write_text(json.dumps(
        {"order_seed": NULL_ORDER_SEED, "n_candidates": int(len(cand)),
         "reserved_per_judged": N_RESERVE_PER_JUDGED, "n_reserved": len(taken),
         "n_excluded_symbols": len(excluded),
         "conditions": {"abs_bias_goal_max": NULL_ABS_BIAS_MAX, "min_detect_start": NULL_MIN_DETECT,
                        "min_donors_start": NULL_MIN_DONORS, "min_cap20": MIN_CELLS},
         "versions": {"A_300": [rec(s) for s in null300], "B_500": [rec(s) for s in null500]},
         "nested": null500[:N_NULL_PRIMARY] == null300}, indent=1, ensure_ascii=False) + "\n")
    (E / "work/subsets_004.json").write_text(json.dumps(
        {"subsets": {k: [rec(s) for s in v] for k, v in subsets.items()},
         "sizes": {k: len(v) for k, v in subsets.items()}}, indent=1, ensure_ascii=False) + "\n")

    # ---------------------------------------------------------------- seed ごとの抽出細胞数
    cov = {}
    for s in allg:
        e = sym2ens[s]
        t = tokd[e]
        idx = Ms.indices[Ms.indptr[t]:Ms.indptr[t + 1]].tolist()
        per = {}
        for sd in SEEDS:
            keep = ri.balanced_subset(idx, list(donors[sidx]), MAX_CELLS, MAX_PER_DONOR,
                                      random.Random(f"{sd}-{e}"))
            per[str(sd)] = {"n_cells": len(keep), "n_donors": len({donors[sidx][i] for i in keep})}
        cov[s] = {"ensembl": e, "per_seed": per,
                  "min_cells": min(v["n_cells"] for v in per.values())}
    (E / "work/precheck_coverage_004.json").write_text(
        json.dumps(cov, indent=1, ensure_ascii=False) + "\n")
    bad = {s: v["min_cells"] for s, v in cov.items() if v["min_cells"] < MIN_CELLS}

    print(f"対照候補 {len(cand)} 個 → 対照 {N_NULL_TOTAL}（先頭 {N_NULL_PRIMARY} が入れ子）")
    print("除外した記号:", len(excluded))
    print("摂動する遺伝子:", json.dumps(
        {r: len([g for g in allg if roles[g] == r]) for r in sorted(set(roles.values()))},
        ensure_ascii=False), "= 合計", len(allg))
    print("セットの大きさ:", json.dumps({k: len(v) for k, v in subsets.items()}, ensure_ascii=False))
    print("落とした遺伝子 T1:", {g: cap20(g) for g in dropped_T1})
    print("落とした遺伝子 T2:", {g: cap20(g) for g in dropped_T2})
    print("主判定の陽性対照 71 個で下限を割る:", {g: cap20(g) for g in pc71_bad} or "なし")
    print("71 個版で目標の細胞の平均が負のマーカー:", G["positive_control71_opposite_6"])
    print("重複の合計:", overlap["重複の合計"], "／盲検混入:", blinded_hit)
    print("EXP-20260929-001 の対照 500 × 陽性対照 71 の重複（引き直す理由）:",
          overlap["EXP-20260929-001 の対照 500 をそのまま使った場合の陽性対照との重複"])
    print("120 細胞を下回る遺伝子:", bad if bad else "なし")
    if bad:
        sys.exit("抽出細胞数が下限を下回る遺伝子がある。設計を見直すこと")


if __name__ == "__main__":
    main()
