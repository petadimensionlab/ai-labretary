#!/usr/bin/env python3
"""段階 3（第 2 版）の数え上げ — procedure-v2.md（2026-09-30T14:20:01+09:00 固定）のとおり。

目標 = 活性化した間質の状態（aFIB + MYOF、Lake 2023 の著者定義）
陽性対照 = Lake 2023 補足表 14 の最小 42（主）／保守 76（感度分析）

入力はすべて読み取りのみ。摂動は 1 回も行わない。
"""
from __future__ import annotations

import csv
import json
import pickle
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp

ROOT = Path("/Users/petadimensionlab/workspace/research/ai-labretary")
OUT = ROOT / "outputs/design/2026-09-30-annotation-goal-state"
LIT = ROOT / "outputs/literature/2026-09-30-independent-goal-state/markers"
H5 = ROOT / "experiments/EXP-20260929-001/work/input.h5ad"
TOKDICT = Path("/Users/petadimensionlab/workspace/research/Geneformer/geneformer_hf/"
               "geneformer/token_dictionary_gc104M.pkl")

# procedure-v2.md §4（順番を固定）
CANDIDATES = [
    ("A", ["aFIB", "MYOF"], ["FIB", "M-FIB"], ["dFIB"], ["dM-FIB"]),
    ("B", ["aFIB", "MYOF"], ["FIB", "M-FIB"], ["dFIB", "dM-FIB"], []),
    ("C", ["aFIB", "MYOF"], ["FIB"], ["M-FIB"], ["dFIB"]),
]
BLINDED = {"ENSG00000141510": "TP53", "ENSG00000118785": "SPP1", "ENSG00000148346": "LCN2"}
BLINDED_SYM = {"TP53", "SPP1", "LCN2"}
MAX_CELLS_PER_DONOR, MIN_CELLS_PER_GENE = 20, 120

T1 = ["CXCL13", "CCL19", "CCL21", "CXCL12", "LTB", "LTBR", "TNFSF13B", "TNFSF14",
      "ICAM1", "VCAM1", "IL7", "CXCL9", "CXCL10", "CXCL11", "CCL2", "MADCAM1"]
T2 = ["CDKN1A", "CDKN2A", "SERPINE1", "IGFBP3", "IGFBP7", "IL6", "IL1B", "TGFB1",
      "GDF15", "TIMP1", "MMP3", "CXCL8"]


def read_csv(name):
    return [r["gene_symbol"] for r in csv.DictReader(open(LIT / name))]


def depth_auc(F, y, donors, seed=20260930):
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import GroupKFold
    from sklearn.preprocessing import StandardScaler
    Z = np.column_stack([np.log1p(F[:, 0]), np.log1p(F[:, 1])])
    aucs, oof = [], np.full(len(y), np.nan)
    for tr, te in GroupKFold(n_splits=5).split(Z, y, groups=donors):
        if len(np.unique(y[tr])) < 2 or len(np.unique(y[te])) < 2:
            continue
        sc = StandardScaler().fit(Z[tr])
        m = LogisticRegression(max_iter=2000, random_state=seed).fit(sc.transform(Z[tr]), y[tr])
        oof[te] = m.predict_proba(sc.transform(Z[te]))[:, 1]
        aucs.append(roc_auc_score(y[te], oof[te]))
    ok = ~np.isnan(oof)
    return {"auc_per_fold": [float(a) for a in aucs],
            "auc_mean": float(np.mean(aucs)) if aucs else float("nan"),
            "auc_pooled_oof": float(roc_auc_score(y[ok], oof[ok])) if ok.sum() else float("nan")}


def main() -> None:
    minimal, conserved = read_csv("adaptive_stromal_state_aFIB_plus_MYOF_minimal.csv"), \
        read_csv("adaptive_stromal_state_aFIB_plus_MYOF_full_conserved.csv")
    gd_proxy = [g for g in conserved if g not in set(minimal)]

    A = ad.read_h5ad(str(H5))
    X = A.X.tocsr() if sp.issparse(A.X) else sp.csr_matrix(np.asarray(A.X))
    sub = A.obs["subclass_l2"].astype(str).to_numpy()
    don = A.obs["donor_id"].astype(str).to_numpy()
    dis = A.obs["disease"].astype(str).to_numpy()
    ng, nc = A.obs["n_genes"].to_numpy(float), A.obs["n_counts"].to_numpy(float)
    sym = A.var["gene_name"].astype(str).to_numpy()
    ens = A.var["ensembl_id"].astype(str).to_numpy()
    sym2i: dict[str, int] = {}
    for i, s in enumerate(sym):
        sym2i.setdefault(s, i)
    tok = pickle.loads(TOKDICT.read_bytes())
    in_tok = np.array([e in tok for e in ens])

    rep = {"generated_at": __import__("datetime").datetime.now().astimezone().isoformat(timespec="seconds"),
           "procedure": "outputs/design/2026-09-30-annotation-goal-state/procedure-v2.md",
           "procedure_fixed_at": "2026-09-30T14:20:01+09:00",
           "positive_control_source": "Lake BB et al. 2023 Nature 619:585-594 (PMID 37468583) Supp Table 14",
           "literature_sets": {"minimal": len(minimal), "conserved": len(conserved),
                               "minimal_is_subset_of_conserved": set(minimal) <= set(conserved),
                               "goal_defining_proxy_conserved_minus_minimal": len(gd_proxy)},
           "candidates": []}

    for name, gl, sl, rl, ol in CANDIDATES:
        g, s = np.isin(sub, gl), np.isin(sub, sl)
        r, o = np.isin(sub, rl), np.isin(sub, ol) if ol else np.zeros(len(sub), bool)
        c = {"candidate": name, "goal_labels": gl, "start_labels": sl,
             "reference_labels": rl, "other_labels": ol,
             "F1": {"n_cells": int(g.sum()), "n_donors": int(len(set(don[g]))),
                    "by_disease": {d: int(((dis == d) & g).sum()) for d in sorted(set(dis[g]))},
                    "pass": bool(g.sum() >= 300 and len(set(don[g])) >= 10)},
             "F2": {"n_cells": int(s.sum()), "n_donors": int(len(set(don[s]))),
                    "by_disease": {d: int(((dis == d) & s).sum()) for d in sorted(set(dis[s]))},
                    "pass": bool(s.sum() >= 2000 and len(set(don[s])) >= 20)},
             "reference": {"n_cells": int(r.sum()), "n_donors": int(len(set(don[r])))},
             "other": {"n_cells": int(o.sum()), "n_donors": int(len(set(don[o]))) if o.sum() else 0}}

        Xs, dsr = X[s], don[s]
        pop = g | s
        Xpop, gmask_pop = X[pop], g[pop]

        def stats(symbols, label):
            rows, kept = [], []
            for sy in symbols:
                rec = {"symbol": sy}
                if sy in BLINDED_SYM:
                    rec.update({"status": "盲検中のため除外"}); rows.append(rec); continue
                i = sym2i.get(sy)
                if i is None:
                    rec.update({"status": "var に無い"}); rows.append(rec); continue
                if not in_tok[i]:
                    rec.update({"status": "Geneformer の辞書に無い", "ensembl": str(ens[i])})
                    rows.append(rec); continue
                col_s = np.asarray(Xs[:, i].todense()).ravel() > 0
                dd = dsr[col_s]
                _, cnt = np.unique(dd, return_counts=True)
                cap = int(np.minimum(cnt, MAX_CELLS_PER_DONOR).sum()) if len(cnt) else 0
                col_g = np.asarray(X[g][:, i].todense()).ravel() > 0
                col_p = np.asarray(Xpop[:, i].todense()).ravel() > 0
                rec.update({"ensembl": str(ens[i]),
                            "detection_goal": float(col_g.mean()),
                            "detection_population": float(col_p.mean()),
                            "n_cells_start": int(col_s.sum()),
                            "n_donors_start": int(len(set(dd))),
                            "n_cells_after_cap": cap,
                            "status": "使える" if cap >= MIN_CELLS_PER_GENE else
                                      f"摂動できない（上限後 {cap} < {MIN_CELLS_PER_GENE}）"})
                if cap >= MIN_CELLS_PER_GENE:
                    kept.append(sy)
                rows.append(rec)
            det = [x["detection_goal"] for x in rows if "detection_goal" in x]
            return {"set": label, "n_requested": len(symbols), "K": len(kept), "genes_used": kept,
                    "detection_goal_median": float(np.median(det)) if det else float("nan"),
                    "detection_goal_min": float(min(det)) if det else float("nan"),
                    "detection_goal_max": float(max(det)) if det else float("nan"),
                    "per_gene": rows}

        pc_min = stats(minimal, "陽性対照（最小 42）")
        pc_con = stats(conserved, "陽性対照（保守 76・感度分析）")
        gdp = stats(gd_proxy, "目標の定義の代理（保守 76 − 最小 42）")
        c["K"] = {"primary_set": "minimal_42", "K": pc_min["K"], "pass": pc_min["K"] >= 20,
                  "minimal": pc_min, "conserved": pc_con, "goal_defining_proxy": gdp}

        t1, t2 = stats(T1, "T1"), stats(T2, "T2")
        pc_used = set(pc_min["genes_used"])
        for row in t1["per_gene"] + t2["per_gene"]:
            row["in_positive_control"] = row["symbol"] in pc_used
        t1_use = [x for x in t1["genes_used"] if x not in pc_used]
        t2_use = [x for x in t2["genes_used"] if x not in pc_used]

        # 全遺伝子で摂動できる数
        Dg = (Xs > 0)
        capall = np.zeros(len(ens), int)
        for j in range(len(ens)):
            col = Dg[:, j]
            if col.nnz:
                dd = dsr[np.asarray(col.todense()).ravel()]
                _, cnt = np.unique(dd, return_counts=True)
                capall[j] = int(np.minimum(cnt, MAX_CELLS_PER_DONOR).sum())
        n_pert = int(((capall >= MIN_CELLS_PER_GENE) & in_tok).sum())
        c["F3"] = {"n_genes_perturbable_whole_universe": n_pert,
                   "T1": t1, "T2": t2, "T1_usable": t1_use, "T2_usable": t2_use,
                   "T1_usable_without_CXCL12": [x for x in t1_use if x != "CXCL12"],
                   "pass": n_pert > 0}

        tgt = set(t1_use) | set(t2_use)
        c["F4"] = {"goal_defined_by": "著者注釈 subclass.l2（Lake 2023 補足表 10 の状態の種別）。遺伝子は使っていない",
                   "n_genes_defining_goal": 0,
                   "overlap_defining_x_positive_control": 0,
                   "overlap_defining_x_targets": 0,
                   "overlap_positive_control_x_targets": len(pc_used & tgt),
                   "overlap_triple": 0,
                   "targets_dropped_because_in_positive_control":
                       sorted((set(t1["genes_used"]) | set(t2["genes_used"])) & pc_used),
                   "goal_defining_proxy_x_positive_control":
                       len(set(gdp["genes_used"]) & pc_used),
                   "CXCL12_in_positive_control": "CXCL12" in set(minimal) | set(conserved),
                   "pass": len(pc_used & tgt) == 0}

        mg, ms = float(np.median(ng[g])), float(np.median(ng[s]))
        ratio = max(mg, ms) / min(mg, ms)
        auc = depth_auc(np.column_stack([ng[pop], nc[pop]]), g[pop].astype(int), don[pop])
        c["F5"] = {"median_n_genes_goal": mg, "median_n_genes_start": ms, "ratio": ratio,
                   "ratio_pass": bool(ratio <= 1.5),
                   "median_n_counts_goal": float(np.median(nc[g])),
                   "median_n_counts_start": float(np.median(nc[s])),
                   "depth_only_auc": auc, "auc_pass": bool(auc["auc_pooled_oof"] < 0.75),
                   "pass": bool(ratio <= 1.5 and auc["auc_pooled_oof"] < 0.75)}

        c["all_pass"] = bool(c["K"]["pass"] and c["F1"]["pass"] and c["F2"]["pass"]
                             and c["F3"]["pass"] and c["F4"]["pass"] and c["F5"]["pass"])
        rep["candidates"].append(c)
        print(f"[{name}] K={pc_min['K']}/42 (保守版 {pc_con['K']}/76) F1={c['F1']['pass']} "
              f"F2={c['F2']['pass']} F3={n_pert} F4={c['F4']['pass']} F5={c['F5']['pass']} "
              f"→ {c['all_pass']}", flush=True)
        if c["all_pass"]:
            rep["adopted"] = name
            break

    (OUT / "candidate_counts_v2.json").write_text(json.dumps(rep, ensure_ascii=False, indent=1) + "\n")
    print("→", OUT / "candidate_counts_v2.json")


if __name__ == "__main__":
    main()
