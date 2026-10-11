#!/usr/bin/env python3
"""EXP-20261002-002 — 追加する対照遺伝子 500 個を決める（**結果を見る前に固定する**）。

## 原理（**新しい規則を 1 つも作らない**）

EXP-20261001-001 の `analysis/select_genes_005.py` は、対照の候補 9,266 個を
**ひとつの決まった並び** に並べ、その**先頭 500 個**を対照として採用した。並びは

  (1) 判定する遺伝子 1 個につき、性質（細胞数・元の順位・入力の長さ）の近い候補を
      2 個ずつ確保する（434 個）
  (2) 残りを固定 seed 20261001 で並べ替えて後ろにつなぐ

で決まる。**本スクリプトはこの並びをやり直し、先頭 500 個が事前登録済みの対照 500 個と
1 個の違いもなく一致することを確かめた上で、501〜1000 番目を追加の対照とする。**

> **「同じ規則で引く」とは、同じ並びの続きを取ることである。**
> 引き直しでも、別の seed での抽選でもない。既存の 500 個は 1 個も入れ替わらない。

一致の確認が取れなければ **その場で止まる**（規則が同じであることを示せないため）。

## 除外するもの（EXP-20261001-001 の除外集合をそのまま使う。550 記号）

介入セット 119 と経路の候補表すべて／陽性対照 71（最小 42 版・保守 76 版の全体）／
目標の定義に関わる遺伝子（Lake 2023 の状態マーカー全体）／中立対照の候補 11／
盲検 5（TP53・SPP1・LCN2 を含む）／前回の摂動対象 T1 11・T2 7／過去実験の各セット。
**既存の対照 500 個は、並びの先頭 500 個としてそのまま残るので重複しない。**

出力:
  work/extra_null_007.json        追加する対照 500 個（並びの 501〜1000 番目）
  work/overlap_check_007.json     **総当たりの重なりの確認**
  work/precheck_coverage_007.json 追加分の seed ごとの抽出細胞数・提供者数
"""
from __future__ import annotations

import json
import random
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "experiments/EXP-20261002-002"
P01 = ROOT / "experiments/EXP-20261001-001"
P29 = ROOT / "experiments/EXP-20260929-001"
DESIGN = ROOT / "outputs/design/2026-10-01-intervention-gene-set"
LIT = ROOT / "outputs/literature/2026-09-30-independent-goal-state/markers"
sys.path.insert(0, str(ROOT / "scripts/geneformer"))
import run_isp as ri  # noqa: E402

# ---- EXP-20261001-001 の select_genes_005.py と**同一の定数**（1 つも変えない） ----
SEEDS = [81, 82, 83]
MIN_CELLS, MAX_CELLS, MAX_PER_DONOR = 120, 200, 20
NULL_ABS_BIAS_MAX = 0.5
NULL_MIN_DETECT, NULL_MIN_DONORS = 250, 12
N_NULL_EXISTING = 500              # 既存の対照（**入れ替えない**）
N_RESERVE_PER_JUDGED = 2
NULL_ORDER_SEED = 20261001         # **並べ替えの seed も変えない**（同じ並びを再現するため）
TOL_N = {"lncells": float(np.log(1.5)), "log2pos": 0.25, "log2len": 0.25}

# ---- この実験で決める唯一の数（根拠は plan.md §3） ----
N_EXTRA = 500                      # 追加する対照の個数 → 合計 1,000


def read_marker_csv(p: Path) -> list[str]:
    return [x.strip() for x in p.read_text().splitlines()[1:] if x.strip()]


def main() -> None:
    G = json.loads((P01 / "analysis/gene_sets_005.json").read_text())

    # gene_stats_005.csv は EXP-20261001-001 の確定済みの出力（**読むだけ**）。
    # CSV では空文字が NaN として読まれるので "" に戻す。元のスクリプトは
    # メモリ上で "" を入れていた。この 1 行が無いと候補が 9,267 個になり並びがずれる。
    d = pd.read_csv(P01 / "work/gene_stats_005.csv")
    d = d[d.ensembl.notna() & (d.ensembl != "")].reset_index(drop=True)
    d["symbol"] = d.symbol.fillna("")
    ref = d.set_index("symbol")
    ref = ref[~ref.index.duplicated()]

    sym2ens = json.loads((P01 / "work/symbol_to_ensembl_005.json").read_text())
    tokd = ri.gene_token_dict()

    # ---------------------------------------------------------------- 除外集合（そのまま）
    lake = json.loads((LIT / "lake2023_fibroblast_state_markers.json").read_text())
    mini42 = read_marker_csv(LIT / "adaptive_stromal_state_aFIB_plus_MYOF_minimal.csv")
    full76 = read_marker_csv(LIT / "adaptive_stromal_state_aFIB_plus_MYOF_full_conserved.csv")

    def usable(g: str) -> bool:
        return g in sym2ens and sym2ens[g] in tokd

    pc39 = [g for g in mini42 if usable(g)]
    pc71 = [g for g in full76 if usable(g)]
    counts = json.loads((DESIGN / "counts.json").read_text())
    INT = [x["symbol"] for x in counts["候補"]["P1_経路セット120"]["採用した遺伝子"]]
    T1, T2 = G["perturb_T1_11"], G["perturb_T2_7"]
    pathway_all = ([x["symbol"] for x in counts["候補"]["P1_経路セット120"]["採用した遺伝子"]]
                   + [x["symbol"] for x in counts["候補"]["P2_候補表の全体"]["採用した遺伝子"]])
    lake_all = sorted({g for v in lake["states"].values() for g in v["genes"]})
    prev_groups = json.loads((P29 / "analysis/excluded_symbols_003.json").read_text())["groups"]
    prev_states = json.loads((P29 / "work/input_states.json").read_text())
    judged = list(dict.fromkeys(INT + pc71 + T1 + T2 + G["housekeeping_3"]
                                + G["housekeeping_extra_descriptive"] + G["negctrl_panel_7"]))
    excluded = (set(judged) | set(lake_all) | set(mini42) | set(full76) | set(pathway_all)
                | set(G["perturb_T1_candidates_16"]) | set(G["perturb_T2_candidates_12"])
                | set(G["blinded_5"]) | set(G["reference_dFIB_markers_6_descriptive"])
                | set(prev_groups["secretory25"]) | set(prev_groups["org_markers7"])
                | set(prev_groups["myof_markers24"]) | set(prev_groups["housekeeping3"])
                | set(prev_groups["negctrl_panel7"]) | set(prev_groups["reference2"])
                | set(prev_states["ecm_candidates_order"]) | set(prev_states["ecm_program"]))

    # ---------------------------------------------------------------- 並びをやり直す
    cand = d[(~d.symbol.isin(excluded)) & (d.symbol != "")
             & (d.n_start >= NULL_MIN_DETECT) & (d.n_donor_start >= NULL_MIN_DONORS)
             & (d.n_cap20 >= MIN_CELLS)
             & (d.bias_goal.abs() < NULL_ABS_BIAS_MAX)].copy()
    cand = cand.set_index("symbol")
    cand = cand[~cand.index.duplicated()]

    order: list[str] = []
    taken: set[str] = set()
    for s in judged:
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

    # ---------------------------------------------------------------- **一致の確認（止まりうる）**
    npool = json.loads((P01 / "work/null_pool_005.json").read_text())
    null500 = [x["symbol"] for x in npool["versions"]["B_500"]]
    repro = {
        "候補の数": {"やり直し": len(cand), "事前登録": npool["n_candidates"]},
        "確保した近傍": {"やり直し": len(taken), "事前登録": npool["n_reserved"]},
        "除外した記号": {"やり直し": len(excluded), "事前登録": npool["n_excluded_symbols"]},
        "並べ替えの seed": {"やり直し": NULL_ORDER_SEED, "事前登録": npool["order_seed"]},
        "先頭 500 が事前登録の対照 500 と完全一致": order[:N_NULL_EXISTING] == null500,
    }
    if not repro["先頭 500 が事前登録の対照 500 と完全一致"]:
        print(json.dumps(repro, ensure_ascii=False, indent=1))
        sys.exit("**止める**: 対照の並びを再現できなかった。『同じ規則で引いた』と言えない")
    if len(order) < N_NULL_EXISTING + N_EXTRA:
        sys.exit(f"候補が足りない（{len(order)} < {N_NULL_EXISTING + N_EXTRA}）")

    extra = order[N_NULL_EXISTING:N_NULL_EXISTING + N_EXTRA]
    if set(extra) & set(null500):
        sys.exit("追加分が既存の対照と重なっている（並びの取り方が誤っている）")

    # ---------------------------------------------------------------- 総当たりの重なりの確認
    pc_extra32 = [g for g in pc71 if g not in set(pc39)]
    adopted = [x["symbol"] for x in
               json.loads((P01 / "work/subsets_005.json").read_text())["subsets"]["neutral_adopted"]]
    #
    # 検査は 2 本立てにする。**ひとつの総当たりに混ぜてはいけない。**
    #
    #  (A) **判定に使うセットどうしの総当たり** — EXP-20261001-001 と同じ組に
    #      「追加する対照 500」を足したもの。重なりは許した 2 件以外 0 でなければならない。
    #  (B) **追加する対照 × 除外する側の一覧** — 目標の定義（Lake 2023 の状態マーカー全体）・
    #      経路の候補表 406・盲検 5 は、判定に使うセットではなく「対照から外す根拠の一覧」である。
    #      これらは互いに重なって当然で（陽性対照 71 はそもそも Lake 2023 から作った）、
    #      (A) に混ぜると意味のない違反が出る。**追加の対照と 1 個も重ならないこと**だけを見る。
    #
    judged_named = {
        "【介入】主判定 119": INT,
        "陽性対照 71": pc71,
        "（内訳）最小 42 版の 39": pc39,
        "（内訳）保守 76 版で加わる 32": pc_extra32,
        "前回の摂動対象 T1 11": T1,
        "前回の摂動対象 T2 7": T2,
        "恒常発現 3": G["housekeeping_3"],
        "恒常発現の追加 1": G["housekeeping_extra_descriptive"],
        "中立対照パネル 7": G["negctrl_panel_7"],
        "採用した中立対照": adopted,
        "既存の対照 500": null500,
        "**追加する対照 500**": extra,
    }
    excl_named = {
        "目標の定義（Lake 2023 の状態マーカー全体）": lake_all,
        "陽性対照の最小 42 版（全体）": mini42,
        "陽性対照の保守 76 版（全体）": full76,
        "経路の候補表すべて": sorted(set(pathway_all)),
        "中立対照の候補 11": (G["housekeeping_3"] + G["housekeeping_extra_descriptive"]
                              + G["negctrl_panel_7"]),
        "盲検 5": G["blinded_5"],
        "盲検のうち解除されていない 3": G["blinded_still_closed_3"],
        "除外集合すべて": sorted(excluded),
    }
    # 入れ子の族（重なって当然。包含関係で検算し、総当たりからは外す）
    FAMILIES = [["陽性対照 71", "（内訳）最小 42 版の 39", "（内訳）保守 76 版で加わる 32"],
                ["恒常発現 3", "恒常発現の追加 1", "中立対照パネル 7", "採用した中立対照"]]
    nested = {"39 ⊂ 71": set(pc39) <= set(pc71),
              "保守で加わる 32 ⊂ 71": set(pc_extra32) <= set(pc71),
              "39 ∪ 32 = 71": sorted(set(pc39) | set(pc_extra32)) == sorted(pc71),
              "採用した中立対照 ⊂ 候補 11": set(adopted) <= set(excl_named["中立対照の候補 11"]),
              "介入 119 ⊂ 経路の候補表": set(INT) <= set(pathway_all),
              "既存 500 ∩ 追加 500 = 空": not (set(null500) & set(extra))}
    # EXP-20261001-001 で許した重なり（**追加の対照には 1 件も許さない**）
    ALLOWED = {frozenset(("【介入】主判定 119", "前回の摂動対象 T1 11")): {"CCL2"},
               frozenset(("【介入】主判定 119", "前回の摂動対象 T2 7")): {"TIMP1"}}
    NESTED_FS = {frozenset((a, b)) for fam in FAMILIES for a, b in combinations(fam, 2)}
    pairs, violations = {}, {}
    for a, b in combinations(judged_named, 2):
        if frozenset((a, b)) in NESTED_FS:
            continue
        ov = sorted(set(judged_named[a]) & set(judged_named[b]))
        pairs[f"{a} × {b}"] = ov
        if set(ov) - ALLOWED.get(frozenset((a, b)), set()):
            violations[f"{a} × {b}"] = ov
    extra_pairs = {k: v for k, v in pairs.items() if "**追加する対照 500**" in k}
    extra_vs_excl = {k: sorted(set(extra) & set(v)) for k, v in excl_named.items()}
    blinded_hit = sorted(set(G["blinded_still_closed_3"]) & set(extra))

    overlap = {
        "(A) 判定に使うセットの大きさ": {k: len(v) for k, v in judged_named.items()},
        "(A) 包含関係の検算": nested,
        "(A) 総当たりの重なり": pairs,
        "(A) 重なりの合計": int(sum(len(v) for v in pairs.values())),
        "(A) 違反（許していない重なり）": violations,
        "(A) 許している重なり": {"介入 × T1": ["CCL2"], "介入 × T2": ["TIMP1"],
                                 "理由": "どちらも過去の介入対象であって目標の定義にも対照にも"
                                         "使われていない（HYP-20261001-01 §6.3）。"
                                         "**追加する対照には 1 件も許していない**"},
        "**(A) 追加する対照 500 と、判定に使う他のセットとの重なり**": extra_pairs,
        "**(A) 追加する対照の重なりの合計**": int(sum(len(v) for v in extra_pairs.values())),
        "(B) 除外する側の一覧の大きさ": {k: len(v) for k, v in excl_named.items()},
        "**(B) 追加する対照 500 と、除外する側の一覧との重なり**": extra_vs_excl,
        "**(B) 重なりの合計**": int(sum(len(v) for v in extra_vs_excl.values())),
        "盲検中の 3 遺伝子の混入（追加分）": blinded_hit,
    }
    (E / "work/overlap_check_007.json").write_text(
        json.dumps(overlap, indent=1, ensure_ascii=False) + "\n")
    if (violations or blinded_hit or not all(nested.values())
            or overlap["**(A) 追加する対照の重なりの合計**"] != 0
            or overlap["**(B) 重なりの合計**"] != 0):
        sys.exit(f"重なりがある: (A) {violations} / (B) "
                 f"{ {k: v for k, v in extra_vs_excl.items() if v} } / 盲検: {blinded_hit}")

    # ---------------------------------------------------------------- seed ごとの抽出細胞数
    from datasets import load_from_disk
    import scipy.sparse as sp
    ds = load_from_disk(str(P01 / "work/tokenized.dataset"))
    st = np.array([str(x) for x in ds["fib_anno_state"]])
    donors = np.array([str(x) for x in ds["donor_id"]])
    ids = [np.asarray(x) for x in ds["input_ids"]]
    indptr = np.cumsum([0] + [len(x) for x in ids])
    col = np.concatenate(ids)
    ntok = max(int(col.max()) + 1, max(tokd.values()) + 1)
    Mb = sp.csr_matrix((np.ones(len(col), bool), col, indptr), shape=(len(ids), ntok)).tocsc()
    sidx = np.flatnonzero(st == "FIB_MFIB")
    Ms = Mb[sidx].tocsc()
    cov = {}
    for s in extra:
        e = sym2ens[s]
        t = tokd[e]
        idx = Ms.indices[Ms.indptr[t]:Ms.indptr[t + 1]].tolist()
        per = {}
        for sd in SEEDS:
            keep = ri.balanced_subset(idx, list(donors[sidx]), MAX_CELLS, MAX_PER_DONOR,
                                      random.Random(f"{sd}-{e}"))
            per[str(sd)] = {"n_cells": len(keep),
                            "n_donors": len({donors[sidx][i] for i in keep})}
        cov[s] = {"ensembl": e, "role": "negative_extra", "per_seed": per,
                  "min_cells": min(v["n_cells"] for v in per.values())}
    (E / "work/precheck_coverage_007.json").write_text(
        json.dumps(cov, indent=1, ensure_ascii=False) + "\n")
    bad = {s: v["min_cells"] for s, v in cov.items() if v["min_cells"] < MIN_CELLS}

    def rec(s: str) -> dict:
        r = ref.loc[s]
        return {"symbol": s, "ensembl": sym2ens[s], "role": "negative_extra",
                "order_index": order.index(s),
                "bias_goal": round(float(r.bias_goal), 3),
                "n_start": int(r.n_start), "n_goal": int(r.n_goal),
                "det_start": round(float(r.det_start), 4),
                "det_goal": round(float(r.det_goal), 4),
                "n_donor_start": int(r.n_donor_start), "n_cap20": int(r.n_cap20),
                "lp_med": round(float(r.lp_med), 4), "ll_med": round(float(r.ll_med), 4)}

    (E / "work/extra_null_007.json").write_text(json.dumps({
        "引き方": "EXP-20261001-001 と同じ規則で作った並びの 501〜1000 番目。"
                  "**既存の 500 個は 1 個も入れ替えない**",
        "再現の確認": repro,
        "条件": {"abs_bias_goal_max": NULL_ABS_BIAS_MAX, "min_detect_start": NULL_MIN_DETECT,
                 "min_donors_start": NULL_MIN_DONORS, "min_cap20": MIN_CELLS},
        "n_candidates": int(len(cand)), "n_existing": N_NULL_EXISTING, "n_extra": N_EXTRA,
        "n_total": N_NULL_EXISTING + N_EXTRA,
        "追加する対照": [rec(s) for s in extra],
        "合計の対照（既存 500 → 追加 500 の順）": null500 + extra,
    }, indent=1, ensure_ascii=False) + "\n")

    print(f"並び全体 {len(order)} 個 ← 候補 {len(cand)}（事前登録 {npool['n_candidates']}）")
    print(f"**先頭 500 が事前登録の対照 500 と完全一致: {repro['先頭 500 が事前登録の対照 500 と完全一致']}**")
    print(f"追加する対照: {N_EXTRA} 個（並びの 501〜{N_NULL_EXISTING + N_EXTRA} 番目）"
          f"→ 合計 {N_NULL_EXISTING + N_EXTRA} 個")
    print("**追加する対照と判定に使う他セットの重なり:",
          overlap["**(A) 追加する対照の重なりの合計**"], "件／除外する側の一覧との重なり:",
          overlap["**(B) 重なりの合計**"], "件**")
    print("盲検中の遺伝子の混入:", blinded_hit or "なし")
    print("120 細胞を下回る遺伝子:", bad if bad else "なし")
    print("追加分の n_cap20: 最小", min(v["n_cap20"] for v in
                                        json.loads((E / "work/extra_null_007.json").read_text())["追加する対照"]))
    if bad:
        sys.exit("抽出細胞数が下限を下回る遺伝子がある。設計を見直すこと")


if __name__ == "__main__":
    main()
