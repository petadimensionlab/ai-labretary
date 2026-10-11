#!/usr/bin/env python3
"""上位候補について、(1) 提供者上限の感度、(2) 深さをそろえた層別の定義、
(3) 従来 FIB_org との濃縮、(4) 提供者の交絡 を数える。"""
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
META = ("row", "donor", "disease", "subclass", "n_counts", "n_genes")
CHEMO3 = ["CXCL13", "CCL19", "CCL21"]
FIB_ORG = (D[CHEMO3] > 0).any(axis=1).values
PERT_POOL = ["CXCL13", "CCL19", "CCL21", "LTB", "TNFSF14", "CXCL16", "CXCL1",
             "CXCL14", "CCL2", "CCL11", "CCL28", "CSF1", "IL7", "IL15", "IL33",
             "IL34", "IL11", "IL18", "IL32", "TNFSF10", "TNFSF15", "CXCL12",
             "VEGFA", "PGF", "ANGPT1"]
ORG7 = ["PDPN", "NGFR", "ICAM1", "VCAM1", "TNFSF13B", "MADCAM1", "ENPP2"]


def cap_counts(mask, gene, caps=(10, 20, 40, 100, 10**9)):
    sub = D.loc[mask & (D[gene] > 0)]
    if not len(sub):
        return {f"cap{c}": 0 for c in caps} | {"n_donors": 0}
    per = sub.groupby("donor").size()
    out = {f"cap{c}": int(np.minimum(per, c).sum()) for c in caps}
    out["n_donors"] = int(per.size)
    out["max_donor_share"] = round(float(per.max() / per.sum()), 3)
    return out


def depth_auc(goal, start):
    m = goal | start
    X = np.column_stack([np.log1p(D.loc[m, "n_counts"]), np.log1p(D.loc[m, "n_genes"])])
    y = goal[m].astype(int)
    grp = D.loc[m, "donor"].values
    aucs = []
    for tr, te in GroupKFold(n_splits=5).split(X, y, grp):
        if len(set(y[te])) < 2:
            continue
        mu, sd = X[tr].mean(0), X[tr].std(0) + 1e-9
        clf = LogisticRegression(max_iter=2000).fit((X[tr] - mu) / sd, y[tr])
        aucs.append(roc_auc_score(y[te], clf.predict_proba((X[te] - mu) / sd)[:, 1]))
    return round(float(np.mean(aucs)), 4) if aucs else None


def enrich(goal):
    a = int((goal & FIB_ORG).sum()); b = int((goal & ~FIB_ORG).sum())
    c = int((~goal & FIB_ORG).sum()); d = int((~goal & ~FIB_ORG).sum())
    orr, p = stats.fisher_exact([[a, b], [c, d]])
    return dict(goal_and_FIBorg=a, OR=round(float(orr), 3), p=float(p),
                pct_FIBorg_in_goal=round(100 * a / max(goal.sum(), 1), 2),
                pct_FIBorg_overall=round(100 * FIB_ORG.mean(), 2))


def report(name, goal, start, goal_genes):
    r = dict(name=name, goal_genes=goal_genes,
             n_goal=int(goal.sum()), n_goal_donors=int(D.loc[goal, "donor"].nunique()),
             n_start=int(start.sum()), n_start_donors=int(D.loc[start, "donor"].nunique()),
             goal_med_ngenes=float(D.loc[goal, "n_genes"].median()),
             start_med_ngenes=float(D.loc[start, "n_genes"].median()),
             depth_only_auc=depth_auc(goal, start),
             enrichment_vs_FIB_org=enrich(goal))
    # 目標の提供者の偏り
    pg = D.loc[goal].groupby("donor").size()
    r["goal_max_donor_share"] = round(float(pg.max() / pg.sum()), 3)
    r["genes"] = {g: cap_counts(start, g) for g in PERT_POOL if g not in goal_genes}
    r["pass120_cap20"] = [g for g, v in r["genes"].items() if v["cap20"] >= 120]
    r["chemo3_pass_cap20"] = [g for g in CHEMO3 if g in r["pass120_cap20"]]
    return r


def depth_stratified_goal(score, nbin=10, top_frac=0.20):
    """深さの十分位ごとに、マーカー得点の上位 top_frac を目標にする（深さの交絡を構造的に断つ）。"""
    dq = pd.qcut(D["n_genes"].rank(method="first"), nbin, labels=False)
    goal = np.zeros(len(D), bool)
    for b in range(nbin):
        idx = np.where(dq.values == b)[0]
        s = score[idx]
        k = max(1, int(round(top_frac * len(idx))))
        # 得点の高い順。同点は乱数ではなく行番号順（再現性のため）
        order = np.lexsort((idx, -s))
        goal[idx[order[:k]]] = True
    return goal


def main():
    res = {}
    # --- 上限（出発点 = 全線維芽細胞）の提供者上限の感度 ---
    allm = np.ones(len(D), bool)
    res["cap_sensitivity_all_fibroblasts"] = {g: cap_counts(allm, g) for g in PERT_POOL}

    cands = {}
    # C' frac_BPL_all × FIB
    dn = pd.read_csv(OUT.parent / "2026-09-25-tls-feasibility" / "donor_A1.csv")
    dn = dn.drop_duplicates("donor_id").set_index("donor_id")
    sc7 = (D[ORG7] > 0).sum(axis=1).values
    for col in ("frac_BPL_all", "frac_BPL_lymph"):
        s = D["donor"].map(dn[col]).values.astype(float)
        valid = ~np.isnan(s)
        q = np.nanquantile(s[valid], [1 / 3, 2 / 3])
        hi = (s >= q[1]) & valid
        g = hi & D["subclass"].isin(["FIB"]).values
        cands[f"C': {col} 上位1/3 × FIB"] = (g, ~g, [])
        ga = (sc7 >= 1) & hi
        cands[f"A+C: マーカー7個中1個以上 かつ {col} 上位1/3"] = (ga, ~ga, ORG7)

    # 深さの十分位内で上位 20% / 10%（深さの交絡を構造的に断つ）
    for tf in (0.20, 0.10):
        g = depth_stratified_goal(sc7.astype(float), 10, tf)
        cands[f"A-depth: 深さ十分位内でマーカー得点 上位{int(tf*100)}%"] = (g, ~g, ORG7)
    # マーカー得点 × 提供者スコアの両方を深さ層内で
    s = D["donor"].map(dn["frac_BPL_all"]).values.astype(float)
    valid = ~np.isnan(s)
    q = np.nanquantile(s[valid], [1 / 3, 2 / 3])
    hi = (s >= q[1]) & valid
    for tf in (0.30, 0.20):
        g = depth_stratified_goal(sc7.astype(float), 10, tf) & hi
        cands[f"A+C-depth: 深さ十分位内 上位{int(tf*100)}% かつ 提供者スコア上位1/3"] = (g, ~g, ORG7)

    res["candidates"] = [report(k, *v) for k, v in cands.items()]
    (OUT / "refine.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))

    print("== 提供者上限の感度（出発点 = 全線維芽細胞 19,049。これが上限）==")
    print(f"{'遺伝子':<10}{'検出':>7}{'提供者':>7}{'最大占有':>8}"
          f"{'cap10':>7}{'cap20':>7}{'cap40':>7}{'cap100':>7}")
    for g in PERT_POOL:
        v = res["cap_sensitivity_all_fibroblasts"][g]
        tot = int((D[g] > 0).sum())
        print(f"{g:<10}{tot:>7}{v['n_donors']:>7}{v['max_donor_share']:>8}"
              f"{v['cap10']:>7}{v['cap20']:>7}{v['cap40']:>7}{v['cap100']:>7}")

    print("\n== 候補 ==")
    for c in res["candidates"]:
        e = c["enrichment_vs_FIB_org"]
        print(f"\n[{c['name']}]")
        print(f"  目標 {c['n_goal']}細胞/{c['n_goal_donors']}提供者（最大占有 {c['goal_max_donor_share']}）"
              f"  出発点 {c['n_start']}細胞/{c['n_start_donors']}提供者")
        print(f"  検出遺伝子数 中央値 目標 {c['goal_med_ngenes']} / 出発点 {c['start_med_ngenes']}"
              f"  深さだけの AUC {c['depth_only_auc']}")
        print(f"  従来 FIB_org の濃縮: 目標の {e['pct_FIBorg_in_goal']}% が FIB_org"
              f"（全体 {e['pct_FIBorg_overall']}%）OR={e['OR']} p={e['p']:.2e}")
        print(f"  3 遺伝子の合格(cap20): {c['chemo3_pass_cap20']}  "
              f"合格遺伝子 {len(c['pass120_cap20'])} 個")
        for g in CHEMO3 + ["LTB", "TNFSF14", "CXCL16"]:
            if g in c["genes"]:
                v = c["genes"][g]
                print(f"    {g:<9} cap20={v['cap20']:<5} cap40={v['cap40']:<5} "
                      f"cap100={v['cap100']:<5} 提供者={v['n_donors']}")


if __name__ == "__main__":
    main()
