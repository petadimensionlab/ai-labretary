#!/usr/bin/env python3
"""採用候補を確定するための最終の数え上げ。
- マーカー得点を連続値にして同点の任意性を除く
- 目標の割合を掃引して余裕（120 からの距離）を見る
- トークン化で切り捨てられる危険（検出遺伝子数 > 入力長）を数える
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score

OUT = Path(__file__).resolve().parent
D = pd.read_csv(OUT / "fib_A1_panel.csv.gz")
CHEMO3 = ["CXCL13", "CCL19", "CCL21"]
FIB_ORG = (D[CHEMO3] > 0).any(axis=1).values
ORG7 = ["PDPN", "NGFR", "ICAM1", "VCAM1", "TNFSF13B", "MADCAM1", "ENPP2"]
PERT_POOL = ["CXCL13", "CCL19", "CCL21", "LTB", "TNFSF14", "CXCL16", "CXCL1",
             "CXCL14", "CCL2", "CCL11", "CCL28", "CSF1", "IL7", "IL15", "IL33",
             "IL34", "IL11", "IL18", "IL32", "TNFSF10", "TNFSF15", "CXCL12",
             "VEGFA", "PGF", "ANGPT1"]
CAP = 20
MIN = 120

cp10k = D[ORG7].values / D["n_counts"].values[:, None] * 1e4
SCORE = np.log1p(cp10k).sum(1)          # 連続値（深さで正規化した発現の和）
NMARK = (D[ORG7] > 0).sum(1).values     # 検出マーカー数（参考）


def capped(mask, gene, cap=CAP):
    sub = D.loc[mask & (D[gene] > 0)]
    if not len(sub):
        return 0, 0
    per = sub.groupby("donor").size()
    return int(np.minimum(per, cap).sum()), int(per.size)


def strat_goal(score, top_frac, nbin=10):
    dq = pd.qcut(D["n_genes"].rank(method="first"), nbin, labels=False).values
    goal = np.zeros(len(D), bool)
    for b in range(nbin):
        idx = np.where(dq == b)[0]
        k = max(1, int(round(top_frac * len(idx))))
        order = np.lexsort((idx, -score[idx]))
        goal[idx[order[:k]]] = True
    return goal


def depth_auc(goal):
    start = ~goal
    X = np.column_stack([np.log1p(D["n_counts"]), np.log1p(D["n_genes"])])
    y = goal.astype(int)
    grp = D["donor"].values
    a = []
    for tr, te in GroupKFold(n_splits=5).split(X, y, grp):
        mu, sd = X[tr].mean(0), X[tr].std(0) + 1e-9
        c = LogisticRegression(max_iter=2000).fit((X[tr] - mu) / sd, y[tr])
        a.append(roc_auc_score(y[te], c.predict_proba((X[te] - mu) / sd)[:, 1]))
    return round(float(np.mean(a)), 4)


def main():
    res = {}
    # 同点の程度
    res["score_ties"] = {
        "n_score_zero": int((SCORE == 0).sum()),
        "pct_score_zero": round(100 * float((SCORE == 0).mean()), 2),
        "n_nmark_ge1": int((NMARK >= 1).sum()),
        "n_unique_scores": int(len(np.unique(SCORE))),
    }
    # トークン化の切り捨ての危険
    res["truncation_risk"] = {
        "median_n_genes_fib": float(D["n_genes"].median()),
        "n_cells_ngenes_gt_2048": int((D["n_genes"] > 2048).sum()),
        "n_cells_ngenes_gt_4096": int((D["n_genes"] > 4096).sum()),
    }
    for g in CHEMO3 + ["LTB", "TNFSF14"]:
        s = D.loc[D[g] > 0]
        res["truncation_risk"][g] = {
            "n_detected": int(len(s)),
            "median_n_genes": float(s["n_genes"].median()),
            "n_gt_2048": int((s["n_genes"] > 2048).sum()),
            "n_gt_4096": int((s["n_genes"] > 4096).sum()),
            "median_count": float(s[g].median()),
        }

    # 目標の割合の掃引
    sweep = []
    for tf in (0.05, 0.075, 0.10, 0.125, 0.15, 0.20, 0.25):
        goal = strat_goal(SCORE, tf)
        start = ~goal
        a = int((goal & FIB_ORG).sum()); b = int((goal & ~FIB_ORG).sum())
        c = int((start & FIB_ORG).sum()); d = int((start & ~FIB_ORG).sum())
        orr, p = stats.fisher_exact([[a, b], [c, d]])
        pg = D.loc[goal].groupby("donor").size()
        row = dict(
            top_frac=tf, n_goal=int(goal.sum()), n_goal_donors=int(D.loc[goal, "donor"].nunique()),
            goal_max_donor_share=round(float(pg.max() / pg.sum()), 3),
            n_start=int(start.sum()), n_start_donors=int(D.loc[start, "donor"].nunique()),
            goal_med_ngenes=float(D.loc[goal, "n_genes"].median()),
            start_med_ngenes=float(D.loc[start, "n_genes"].median()),
            depth_only_auc=depth_auc(goal),
            FIBorg_in_goal=a, pct_FIBorg_in_goal=round(100 * a / goal.sum(), 2),
            OR_FIBorg=round(float(orr), 3), p_FIBorg=float(p),
            FIBorg_left_in_start=c,
        )
        for g in PERT_POOL:
            n, nd = capped(start, g)
            row[f"cap20_{g}"] = n
        row["n_pass120"] = sum(1 for g in PERT_POOL if row[f"cap20_{g}"] >= MIN)
        row["chemo3_pass"] = [g for g in CHEMO3 if row[f"cap20_{g}"] >= MIN]
        row["margin_CXCL13"] = row["cap20_CXCL13"] - MIN
        sweep.append(row)
    res["sweep"] = sweep
    (OUT / "final_counts.json").write_text(json.dumps(res, indent=1, ensure_ascii=False, default=str))

    t = res["score_ties"]
    print(f"得点 0 の細胞 {t['n_score_zero']}（{t['pct_score_zero']}%）/ "
          f"マーカー 1 個以上 {t['n_nmark_ge1']} / 得点の異なる値 {t['n_unique_scores']}")
    r = res["truncation_risk"]
    print(f"\n検出遺伝子数 中央値 {r['median_n_genes_fib']} / >2048 {r['n_cells_ngenes_gt_2048']} 細胞 "
          f"/ >4096 {r['n_cells_ngenes_gt_4096']} 細胞")
    for g in CHEMO3 + ["LTB", "TNFSF14"]:
        v = r[g]
        print(f"  {g:<9} 検出 {v['n_detected']:<5} 中央 n_genes {v['median_n_genes']:<7} "
              f">2048 {v['n_gt_2048']:<4} >4096 {v['n_gt_4096']:<4} カウント中央 {v['median_count']}")

    print("\n== 目標の割合の掃引（深さの十分位内で連続得点の上位）==")
    hdr = ("割合 目標 提供者 最大占有 出発点 深さAUC 目標n_genes/出発点 "
           "FIBorg濃縮 CXCL13 CCL19 CCL21 LTB TNFSF14 合格数")
    print(hdr)
    for s in sweep:
        print(f"{s['top_frac']:<5} {s['n_goal']:<5} {s['n_goal_donors']:<4} "
              f"{s['goal_max_donor_share']:<6} {s['n_start']:<6} {s['depth_only_auc']:<7} "
              f"{s['goal_med_ngenes']:.0f}/{s['start_med_ngenes']:.0f}   "
              f"OR{s['OR_FIBorg']:<5} {s['cap20_CXCL13']:<7}{s['cap20_CCL19']:<6}"
              f"{s['cap20_CCL21']:<6}{s['cap20_LTB']:<5}{s['cap20_TNFSF14']:<8}{s['n_pass120']}")


if __name__ == "__main__":
    main()
