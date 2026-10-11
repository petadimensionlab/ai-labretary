#!/usr/bin/env python3
"""候補ごとに、出発点・目標の細胞数と、出発点での摂動対象候補の検出細胞数を数える。

procedure.md の F-1〜F-5 を機械的に当てはめる。
入力: fib_A1_panel.csv.gz（extract_a1_fib.py の出力）
出力: candidate_counts.json, pert_feasibility.csv
"""
import json
import itertools
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score

OUT = Path(__file__).resolve().parent
D = pd.read_csv(OUT / "fib_A1_panel.csv.gz")
META = ("row", "donor", "disease", "subclass", "n_counts", "n_genes")
GENES = [c for c in D.columns if c not in META]

CAP = 20          # 提供者あたりの上限
MIN_CELLS = 120   # 1 遺伝子あたりの下限

CHEMO3 = ["CXCL13", "CCL19", "CCL21"]

# 摂動の対象にしうる遺伝子（分泌プログラム）。目標の定義に使うものは候補ごとに除く。
PERT_POOL = [
    "CXCL13", "CCL19", "CCL21", "LTB", "TNFSF14", "CXCL16", "CXCL1", "CXCL14",
    "CCL2", "CCL11", "CCL28", "CSF1", "IL7", "IL15", "IL33", "IL34", "IL11",
    "IL18", "IL32", "TNFSF10", "TNFSF15", "CXCL12", "VEGFA", "PGF", "ANGPT1",
]

# 目標の定義に使う候補（三次リンパ組織を形づくる細胞のマーカー。ケモカイン 3 遺伝子を含まない）
ORG_MARKERS_FULL = ["PDPN", "NGFR", "ICAM1", "VCAM1", "TNFSF13B", "MADCAM1",
                    "ENPP2", "CR2", "TNFSF13", "BST1"]

FIB_ORG = (D[CHEMO3] > 0).any(axis=1).values   # 従来の目標（716 細胞 / 64 提供者の想定）


def capped(mask, gene, cap=CAP):
    """mask の細胞のうち gene を検出する細胞を、提供者あたり cap で打ち切った数。"""
    sub = D.loc[mask & (D[gene] > 0)]
    if len(sub) == 0:
        return 0, 0, 0
    per = sub.groupby("donor").size()
    return int(len(sub)), int(per.size), int(np.minimum(per, cap).sum())


def depth_auc(goal, start):
    """深さ 2 変数だけで目標と出発点を分けたときの提供者分割 AUC。"""
    m = goal | start
    X = np.column_stack([np.log1p(D.loc[m, "n_counts"]), np.log1p(D.loc[m, "n_genes"])])
    y = goal[m].astype(int)
    grp = D.loc[m, "donor"].values
    if y.sum() < 20 or (1 - y).sum() < 20 or len(set(grp)) < 5:
        return None
    aucs = []
    for tr, te in GroupKFold(n_splits=5).split(X, y, grp):
        if len(set(y[te])) < 2:
            continue
        mu, sd = X[tr].mean(0), X[tr].std(0) + 1e-9
        clf = LogisticRegression(max_iter=2000).fit((X[tr] - mu) / sd, y[tr])
        aucs.append(roc_auc_score(y[te], clf.predict_proba((X[te] - mu) / sd)[:, 1]))
    return float(np.mean(aucs)) if aucs else None


def jac(a, b):
    u = (a | b).sum()
    return float((a & b).sum() / u) if u else 0.0


def evaluate(name, goal, start, goal_genes, pert_pool=None, note=""):
    pert_pool = pert_pool or [g for g in PERT_POOL if g not in goal_genes]
    r = {
        "name": name, "note": note,
        "goal_genes": goal_genes,
        "n_goal": int(goal.sum()), "n_goal_donors": int(D.loc[goal, "donor"].nunique()),
        "n_start": int(start.sum()), "n_start_donors": int(D.loc[start, "donor"].nunique()),
        "goal_median_ngenes": float(D.loc[goal, "n_genes"].median()) if goal.sum() else None,
        "start_median_ngenes": float(D.loc[start, "n_genes"].median()) if start.sum() else None,
        "overlap_with_goal_genes_and_pert": sorted(set(goal_genes) & set(pert_pool)),
    }
    r["ngenes_ratio"] = (round(r["goal_median_ngenes"] / r["start_median_ngenes"], 3)
                         if r["start_median_ngenes"] else None)
    r["depth_only_auc"] = depth_auc(goal, start)
    # 従来の FIB_org との関係
    r["FIB_org_n"] = int(FIB_ORG.sum())
    r["jaccard_goal_vs_FIB_org"] = round(jac(goal, FIB_ORG), 4)
    r["frac_goal_that_is_FIB_org"] = round(float((goal & FIB_ORG).sum() / max(goal.sum(), 1)), 4)
    r["frac_FIB_org_in_goal"] = round(float((goal & FIB_ORG).sum() / FIB_ORG.sum()), 4)
    r["FIB_org_cells_left_in_start"] = int((start & FIB_ORG).sum())
    # 摂動対象の実行可能性
    tab = []
    for g in pert_pool:
        nd, ndon, ncap = capped(start, g)
        tab.append(dict(candidate=name, gene=g, n_detected_start=nd,
                        n_donors=ndon, n_after_cap20=ncap, pass120=bool(ncap >= MIN_CELLS)))
    r["pert_table"] = tab
    ok = [t["gene"] for t in tab if t["pass120"]]
    r["pert_pass"] = ok
    r["n_pert_pass"] = len(ok)
    r["chemo3_pass"] = [g for g in CHEMO3 if g in ok]
    # F 判定
    r["F1_goal"] = bool(r["n_goal"] >= 300 and r["n_goal_donors"] >= 10)
    r["F2_start"] = bool(r["n_start"] >= 2000 and r["n_start_donors"] >= 20)
    r["F3_pert"] = bool(len(ok) >= 3 and len(r["chemo3_pass"]) >= 1)
    r["F4_circular"] = bool(len(r["overlap_with_goal_genes_and_pert"]) == 0)
    r["F5_depth"] = bool(r["ngenes_ratio"] is not None and r["ngenes_ratio"] <= 1.5
                         and (r["depth_only_auc"] is None or r["depth_only_auc"] < 0.75))
    r["PASS"] = all(r[k] for k in ("F1_goal", "F2_start", "F3_pert", "F4_circular", "F5_depth"))
    return r


def main():
    res = {"basis": {
        "n_fibroblasts": int(len(D)), "n_donors": int(D.donor.nunique()),
        "cap_per_donor": CAP, "min_cells_per_gene": MIN_CELLS,
        "FIB_org_n_cells": int(FIB_ORG.sum()),
        "FIB_org_n_donors": int(D.loc[FIB_ORG, "donor"].nunique()),
    }}

    # 上限の確認: 出発点を全線維芽細胞にしたときの摂動対象の実行可能性
    allmask = np.ones(len(D), bool)
    up = []
    for g in PERT_POOL:
        nd, ndon, ncap = capped(allmask, g)
        up.append(dict(candidate="上限（全線維芽細胞）", gene=g, n_detected_start=nd,
                       n_donors=ndon, n_after_cap20=ncap, pass120=bool(ncap >= MIN_CELLS)))
    res["upper_bound_all_fibroblasts"] = up

    cands = []

    # ---- (A) 三次リンパ組織を形づくる細胞のマーカーで定義する ----
    for markers in [ORG_MARKERS_FULL,
                    ["PDPN", "NGFR", "ICAM1", "VCAM1", "TNFSF13B", "MADCAM1", "ENPP2"],
                    ["PDPN", "ICAM1", "VCAM1", "TNFSF13B"],
                    ["PDPN", "NGFR", "MADCAM1", "ENPP2", "CR2", "TNFSF13B"]]:
        sc = (D[markers] > 0).sum(axis=1).values
        for k in (1, 2, 3):
            goal = sc >= k
            start = ~goal
            cands.append(evaluate(
                f"A: オーガナイザー様マーカー {len(markers)} 個中 {k} 個以上",
                goal, start, markers,
                note="/".join(markers)))

    # ---- (C) 提供者単位の三次リンパ組織スコア × サブタイプ ----
    dn = pd.read_csv(OUT.parent / "2026-09-25-tls-feasibility" / "donor_A1.csv")
    dn = dn.drop_duplicates("donor_id").set_index("donor_id")
    # 循環しないスコアだけを使う（細胞の組成のみ。CXCL13/CCL19/CCL21/LTB を使わない）
    for scorecol in ("frac_BPL_lymph", "frac_BPL_all"):
        s = D["donor"].map(dn[scorecol]).values.astype(float)
        valid = ~np.isnan(s)
        q = np.nanquantile(s[valid], [1 / 3, 2 / 3])
        hi, lo = (s >= q[1]) & valid, (s <= q[0]) & valid
        for sub in (["FIB"], ["FIB", "aFIB"], ["FIB", "M-FIB", "aFIB", "dFIB", "dM-FIB"]):
            insub = D["subclass"].isin(sub).values
            goal = hi & insub
            start = lo & insub
            cands.append(evaluate(
                f"C: {scorecol} 上位 1/3 の提供者 × {'+'.join(sub)}",
                goal, start, [],
                note=f"スコアは B 細胞・形質細胞の割合のみ（遺伝子を使わない）。"
                     f"出発点は同スコア下位 1/3 の同サブタイプ"))
            # 出発点を「目標以外の全線維芽細胞」に広げた版
            cands.append(evaluate(
                f"C': {scorecol} 上位 1/3 の提供者 × {'+'.join(sub)}（出発点=残り全部）",
                goal, ~goal, [],
                note="出発点を目標以外の全線維芽細胞に広げた版"))

    # ---- (A+C) 組み合わせ ----
    markers = ["PDPN", "NGFR", "ICAM1", "VCAM1", "TNFSF13B", "MADCAM1", "ENPP2"]
    sc = (D[markers] > 0).sum(axis=1).values
    s = D["donor"].map(dn["frac_BPL_lymph"]).values.astype(float)
    valid = ~np.isnan(s)
    q = np.nanquantile(s[valid], [1 / 3, 2 / 3])
    hi = (s >= q[1]) & valid
    for k in (1, 2):
        goal = (sc >= k) & hi
        cands.append(evaluate(
            f"A+C: マーカー {k} 個以上 かつ 三次リンパ組織スコア上位 1/3 の提供者",
            goal, ~goal, markers, note="/".join(markers)))

    res["candidates"] = cands
    (OUT / "candidate_counts.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False))
    rows = list(up)
    for c in cands:
        rows += c["pert_table"]
    pd.DataFrame(rows).to_csv(OUT / "pert_feasibility.csv", index=False)

    # 画面表示
    print(f"線維芽細胞 {len(D)} / 提供者 {D.donor.nunique()} / 従来 FIB_org "
          f"{FIB_ORG.sum()} 細胞 {D.loc[FIB_ORG,'donor'].nunique()} 提供者")
    print("\n== 上限（出発点=全線維芽細胞）で 120 を満たす遺伝子 ==")
    print(", ".join(f"{t['gene']}({t['n_after_cap20']})" for t in up if t["pass120"]))
    print("満たさない:", ", ".join(f"{t['gene']}({t['n_after_cap20']})"
                                   for t in up if not t["pass120"]))
    print("\n== 候補 ==")
    for c in cands:
        print(f"[{'PASS' if c['PASS'] else 'fail'}] {c['name']}")
        print(f"   目標 {c['n_goal']}細胞/{c['n_goal_donors']}提供者  "
              f"出発点 {c['n_start']}細胞/{c['n_start_donors']}提供者  "
              f"深さ比 {c['ngenes_ratio']} AUC {c['depth_only_auc']}")
        print(f"   F1{int(c['F1_goal'])} F2{int(c['F2_start'])} F3{int(c['F3_pert'])} "
              f"F4{int(c['F4_circular'])} F5{int(c['F5_depth'])} | "
              f"3遺伝子合格 {c['chemo3_pass']} | 合格遺伝子 {c['n_pert_pass']}個")
        print(f"   FIB_org と Jaccard {c['jaccard_goal_vs_FIB_org']}, "
              f"出発点に残る FIB_org {c['FIB_org_cells_left_in_start']}")


if __name__ == "__main__":
    main()
