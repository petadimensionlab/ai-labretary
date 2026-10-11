#!/usr/bin/env python3
"""採用候補（深さの十分位内でマーカー得点 上位 5%）の最終集計。
提供者上限の感度、循環性（目標の定義に使う遺伝子と摂動対象の重なり・共起）、従来 FIB_org との関係。"""
import json
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

OUT = Path(__file__).resolve().parent
D = pd.read_csv(OUT / "fib_A1_panel.csv.gz")
ORG7 = ["PDPN", "NGFR", "ICAM1", "VCAM1", "TNFSF13B", "MADCAM1", "ENPP2"]
CHEMO3 = ["CXCL13", "CCL19", "CCL21"]
PERT_POOL = ["CXCL13", "CCL19", "CCL21", "LTB", "TNFSF14", "CXCL16", "CXCL1",
             "CXCL14", "CCL2", "CCL11", "CCL28", "CSF1", "IL7", "IL15", "IL33",
             "IL34", "IL11", "IL18", "IL32", "TNFSF10", "TNFSF15", "CXCL12",
             "VEGFA", "PGF", "ANGPT1"]
FIB_ORG = (D[CHEMO3] > 0).any(axis=1).values
cp = D[ORG7].values / D["n_counts"].values[:, None] * 1e4
SCORE = np.log1p(cp).sum(1)


def strat_goal(tf, nbin=10):
    dq = pd.qcut(D["n_genes"].rank(method="first"), nbin, labels=False).values
    g = np.zeros(len(D), bool)
    for b in range(nbin):
        idx = np.where(dq == b)[0]
        k = max(1, int(round(tf * len(idx))))
        g[idx[np.lexsort((idx, -SCORE[idx]))[:k]]] = True
    return g


goal = strat_goal(0.05)
start = ~goal
res = {"definition": {
    "goal": "線維芽細胞（MYOF/cycMYOF を除く 19,049）を検出遺伝子数の十分位に分け、"
            "各十分位内でオーガナイザー様マーカー 7 個の得点 log1p(CP10K) の和の上位 5% を目標とする",
    "markers": ORG7, "top_frac": 0.05, "n_depth_bins": 10,
    "start": "残りの線維芽細胞すべて"}}
res["n_goal"] = int(goal.sum()); res["n_goal_donors"] = int(D.loc[goal, "donor"].nunique())
res["n_start"] = int(start.sum()); res["n_start_donors"] = int(D.loc[start, "donor"].nunique())
res["goal_by_subclass"] = D.loc[goal, "subclass"].value_counts().to_dict()
res["goal_by_disease"] = D.loc[goal, "disease"].value_counts().to_dict()
res["start_by_subclass"] = D.loc[start, "subclass"].value_counts().to_dict()
pg = D.loc[goal].groupby("donor").size()
res["goal_max_donor_share"] = round(float(pg.max() / pg.sum()), 3)
res["goal_med_ngenes"] = float(D.loc[goal, "n_genes"].median())
res["start_med_ngenes"] = float(D.loc[start, "n_genes"].median())
res["goal_med_ncounts"] = float(D.loc[goal, "n_counts"].median())
res["start_med_ncounts"] = float(D.loc[start, "n_counts"].median())

# 提供者上限の感度
tab = []
for g in PERT_POOL:
    sub = D.loc[start & (D[g] > 0)]
    per = sub.groupby("donor").size() if len(sub) else pd.Series(dtype=int)
    row = dict(gene=g, n_detected_start=int(len(sub)), n_donors=int(per.size),
               max_donor_share=round(float(per.max() / per.sum()), 3) if len(sub) else None)
    for c in (10, 20, 40, 100):
        row[f"cap{c}"] = int(np.minimum(per, c).sum()) if len(sub) else 0
    # 循環性: 目標との共起
    pos = (D[g] > 0).values
    a, b = int((pos & goal).sum()), int((pos & start).sum())
    cc, d = int((~pos & goal).sum()), int((~pos & start).sum())
    orr, p = stats.fisher_exact([[a, b], [cc, d]])
    row["n_in_goal"] = a
    row["pct_of_gene_cells_in_goal"] = round(100 * a / max(a + b, 1), 2)
    row["OR_with_goal"] = round(float(orr), 3)
    row["p_with_goal"] = float(p)
    tab.append(row)
T = pd.DataFrame(tab)
T.to_csv(OUT / "adopted_pert_table.csv", index=False)
res["pert_table"] = tab
res["pass120_cap20"] = [r["gene"] for r in tab if r["cap20"] >= 120]
res["pass120_cap40"] = [r["gene"] for r in tab if r["cap40"] >= 120]
res["overlap_goalgenes_pertgenes"] = sorted(set(ORG7) & set(PERT_POOL))

a = int((goal & FIB_ORG).sum())
orr, p = stats.fisher_exact([[a, int(goal.sum()) - a,],
                             [int(FIB_ORG.sum()) - a, int(start.sum()) - (int(FIB_ORG.sum()) - a)]])
res["vs_FIB_org"] = dict(
    FIB_org_n=int(FIB_ORG.sum()), FIB_org_donors=int(D.loc[FIB_ORG, "donor"].nunique()),
    goal_and_FIB_org=a, pct_of_goal=round(100 * a / goal.sum(), 2),
    pct_of_FIB_org=round(100 * a / FIB_ORG.sum(), 2),
    jaccard=round(float(a / ((goal | FIB_ORG).sum())), 4),
    OR=round(float(orr), 3), p=float(p),
    FIB_org_left_in_start=int((start & FIB_ORG).sum()))

# マーカーごとの目標への寄与
res["marker_detail"] = {m: dict(n_detected=int((D[m] > 0).sum()),
                                n_in_goal=int(((D[m] > 0) & goal).sum()),
                                pct_goal_positive=round(100 * float(((D[m] > 0) & goal).sum() / goal.sum()), 2))
                        for m in ORG7}
(OUT / "adopted.json").write_text(json.dumps(res, indent=1, ensure_ascii=False, default=str))

print(f"目標 {res['n_goal']} 細胞 / {res['n_goal_donors']} 提供者（最大占有 {res['goal_max_donor_share']}）")
print(f"出発点 {res['n_start']} 細胞 / {res['n_start_donors']} 提供者")
print(f"検出遺伝子数 中央値 目標 {res['goal_med_ngenes']} / 出発点 {res['start_med_ngenes']}")
print(f"サブタイプ 目標: {res['goal_by_subclass']}")
print(f"疾患 目標: {res['goal_by_disease']}")
print(f"従来 FIB_org: {res['vs_FIB_org']}")
print("\nマーカーごと:", json.dumps(res["marker_detail"], ensure_ascii=False))
print("\n", T[["gene", "n_detected_start", "n_donors", "cap10", "cap20", "cap40", "cap100",
               "n_in_goal", "pct_of_gene_cells_in_goal", "OR_with_goal"]].to_string(index=False))
print("\ncap20 で 120 以上:", len(res["pass120_cap20"]), res["pass120_cap20"])
print("cap40 で 120 以上:", len(res["pass120_cap40"]))
