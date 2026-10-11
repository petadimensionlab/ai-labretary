#!/usr/bin/env python3
"""段階 3 の数え上げ — procedure.md（2026-09-30T13:59:44+09:00 固定）のとおりに数える。

入力はすべて読み取りのみ。書き出すのは outputs/design/2026-09-30-annotation-goal-state/ の下だけ。
摂動は 1 回も行わない。
"""
from __future__ import annotations

import json
import pickle
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp

ROOT = Path("/Users/petadimensionlab/workspace/research/ai-labretary")
OUT = ROOT / "outputs/design/2026-09-30-annotation-goal-state"
H5 = ROOT / "experiments/EXP-20260929-001/work/input.h5ad"
TOKDICT = Path("/Users/petadimensionlab/workspace/research/Geneformer/geneformer_hf/"
               "geneformer/token_dictionary_gc104M.pkl")

# --- procedure.md §4 の候補（順番を固定） ---
CANDIDATES = [
    ("I",   ["aFIB"],          ["FIB", "M-FIB"], ["MYOF"],  ["dFIB"]),
    ("II",  ["dFIB"],          ["FIB", "M-FIB"], ["MYOF"],  ["aFIB"]),
    ("III", ["aFIB", "MYOF"],  ["FIB", "M-FIB"], ["dFIB"],  ["M-FIB"]),
    ("IV",  ["M-FIB"],         ["FIB"],          ["MYOF"],  ["aFIB"]),
]

# --- procedure.md §5.1 で外す遺伝子 ---
BLINDED = {"ENSG00000141510": "TP53", "ENSG00000118785": "SPP1", "ENSG00000148346": "LCN2"}
MYOF_PC = [  # EXP-20260929-001 config.yaml controls.positive と同一
    "ENSG00000133392", "ENSG00000156113", "ENSG00000078098", "ENSG00000150471",
    "ENSG00000164176", "ENSG00000132334", "ENSG00000135842", "ENSG00000119681",
    "ENSG00000137941", "ENSG00000151320", "ENSG00000113657", "ENSG00000078018",
    "ENSG00000186340", "ENSG00000186479", "ENSG00000182667", "ENSG00000122641",
    "ENSG00000081923", "ENSG00000206190", "ENSG00000196923", "ENSG00000088756",
    "ENSG00000250305", "ENSG00000161791", "ENSG00000134243", "ENSG00000164483"]
HOUSEKEEPING = ["ENSG00000075624", "ENSG00000111640", "ENSG00000196262"]
NEUTRAL_PANEL = ["ENSG00000090661", "ENSG00000147687", "ENSG00000081870", "ENSG00000101558",
                 "ENSG00000054148", "ENSG00000166295", "ENSG00000104529"]

# --- procedure.md §5.3 の閾値（固定） ---
MIN_DET_GOAL = 0.20
MIN_CONSISTENCY = 0.70
MIN_DELTA = 0.05
MIN_DONOR_CELLS = 10
N_PC = 24
N_GD = 12          # goal_defining に回す個数（§5.5）

# --- procedure.md §2 ---
MAX_CELLS_PER_DONOR = 20
MIN_CELLS_PER_GENE = 120

# --- procedure.md §6 の摂動対象（記号） ---
T1 = ["CXCL13", "CCL19", "CCL21", "CXCL12", "LTB", "LTBR", "TNFSF13B", "TNFSF14",
      "ICAM1", "VCAM1", "IL7", "CXCL9", "CXCL10", "CXCL11", "CCL2", "MADCAM1"]
T2 = ["CDKN1A", "CDKN2A", "SERPINE1", "IGFBP3", "IGFBP7", "IL6", "IL1B", "TGFB1",
      "GDF15", "TIMP1", "MMP3", "CXCL8"]


def capped_cells(det_mask: np.ndarray, donors: np.ndarray, cap: int) -> int:
    """提供者上限 cap を当てた後の細胞数。"""
    if det_mask.sum() == 0:
        return 0
    d = donors[det_mask]
    _, cnt = np.unique(d, return_counts=True)
    return int(np.minimum(cnt, cap).sum())


def depth_auc(X_feat: np.ndarray, y: np.ndarray, donors: np.ndarray, seed: int = 20260930) -> dict:
    """深さ 2 変数だけのロジスティック回帰、提供者で分けた 5 分割の AUC。"""
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import GroupKFold
    from sklearn.preprocessing import StandardScaler
    Z = np.column_stack([np.log1p(X_feat[:, 0]), np.log1p(X_feat[:, 1])])
    gkf = GroupKFold(n_splits=5)
    aucs, oof = [], np.full(len(y), np.nan)
    for tr, te in gkf.split(Z, y, groups=donors):
        if len(np.unique(y[tr])) < 2 or len(np.unique(y[te])) < 2:
            continue
        sc = StandardScaler().fit(Z[tr])
        m = LogisticRegression(max_iter=2000, random_state=seed).fit(sc.transform(Z[tr]), y[tr])
        p = m.predict_proba(sc.transform(Z[te]))[:, 1]
        oof[te] = p
        aucs.append(roc_auc_score(y[te], p))
    ok = ~np.isnan(oof)
    return {"auc_per_fold": [float(a) for a in aucs],
            "auc_mean": float(np.mean(aucs)) if aucs else float("nan"),
            "auc_pooled_oof": float(roc_auc_score(y[ok], oof[ok])) if ok.sum() else float("nan"),
            "n_folds": len(aucs)}


def main() -> None:
    A = ad.read_h5ad(str(H5))
    X = A.X.tocsr() if sp.issparse(A.X) else sp.csr_matrix(np.asarray(A.X))
    D = (X > 0)                                   # 検出行列（生カウント ≥ 1）
    sub = A.obs["subclass_l2"].astype(str).to_numpy()
    donors = A.obs["donor_id"].astype(str).to_numpy()
    disease = A.obs["disease"].astype(str).to_numpy()
    n_genes = A.obs["n_genes"].to_numpy(dtype=float)
    n_counts = A.obs["n_counts"].to_numpy(dtype=float)
    ens = A.var["ensembl_id"].astype(str).to_numpy()
    sym = A.var["gene_name"].astype(str).to_numpy()
    ens2i = {e: i for i, e in enumerate(ens)}
    sym2i: dict[str, int] = {}
    for i, s in enumerate(sym):
        sym2i.setdefault(s, i)

    tok = pickle.loads(TOKDICT.read_bytes())
    in_tok = np.array([e in tok for e in ens])

    # §5.1 の除外
    drop = np.zeros(len(ens), bool)
    drop |= np.isin(ens, list(BLINDED))
    drop |= np.isin(ens, MYOF_PC)
    drop |= np.isin(ens, HOUSEKEEPING + NEUTRAL_PANEL)
    drop |= np.array([s.startswith("MT-") or s.startswith("RPL") or s.startswith("RPS")
                      for s in sym])
    eligible = in_tok & ~drop

    report: dict = {
        "generated_at": __import__("datetime").datetime.now().astimezone().isoformat(timespec="seconds"),
        "procedure": "outputs/design/2026-09-30-annotation-goal-state/procedure.md",
        "procedure_fixed_at": "2026-09-30T13:59:44+09:00",
        "input_h5ad": str(H5.relative_to(ROOT)),
        "n_cells_total": int(A.n_obs), "n_genes_total": int(A.n_vars),
        "gene_universe": {
            "in_var": int(len(ens)),
            "in_geneformer_v2_dict": int(in_tok.sum()),
            "excluded_blinded": int(np.isin(ens, list(BLINDED)).sum()),
            "excluded_myof_positive_control": int(np.isin(ens, MYOF_PC).sum()),
            "excluded_housekeeping_neutral": int(np.isin(ens, HOUSEKEEPING + NEUTRAL_PANEL).sum()),
            "excluded_mt_ribo": int(np.array([s.startswith("MT-") or s.startswith("RPL")
                                              or s.startswith("RPS") for s in sym]).sum()),
            "eligible_for_positive_control": int(eligible.sum()),
        },
        "state_counts": {},
        "candidates": [],
    }
    for lab in sorted(set(sub)):
        m = sub == lab
        report["state_counts"][lab] = {
            "n_cells": int(m.sum()), "n_donors": int(len(set(donors[m]))),
            "by_disease": {d: int(((disease == d) & m).sum()) for d in sorted(set(disease[m]))}}

    for name, goal_l, start_l, ref_l, other_l in CANDIDATES:
        g = np.isin(sub, goal_l)
        s = np.isin(sub, start_l)
        r = np.isin(sub, ref_l)
        o = np.isin(sub, other_l)
        c: dict = {"candidate": name, "goal_labels": goal_l, "start_labels": start_l,
                   "reference_labels": ref_l, "other_labels": other_l}

        c["F1"] = {"n_cells": int(g.sum()), "n_donors": int(len(set(donors[g]))),
                   "pass": bool(g.sum() >= 300 and len(set(donors[g])) >= 10)}
        c["F2"] = {"n_cells": int(s.sum()), "n_donors": int(len(set(donors[s]))),
                   "pass": bool(s.sum() >= 2000 and len(set(donors[s])) >= 20)}
        c["reference_n_cells"] = int(r.sum()); c["other_n_cells"] = int(o.sum())

        # ---- §5.2 提供者ごとの検出率 ----
        dg, ds = donors[g], donors[s]
        ug, cg = np.unique(dg, return_counts=True)
        us, cs = np.unique(ds, return_counts=True)
        okg = set(ug[cg >= MIN_DONOR_CELLS]); oks = set(us[cs >= MIN_DONOR_CELLS])
        paired = sorted(okg & oks)
        c["paired_donors"] = {"n": len(paired), "min_cells_each_side": MIN_DONOR_CELLS}

        Dg = D[g][:, eligible]; Ds = D[s][:, eligible]
        idx_el = np.flatnonzero(eligible)
        det_goal_all = np.asarray(Dg.mean(axis=0)).ravel()

        P_g = np.zeros((len(paired), Dg.shape[1])); P_s = np.zeros_like(P_g)
        for k, d in enumerate(paired):
            P_g[k] = np.asarray(Dg[dg == d].mean(axis=0)).ravel()
            P_s[k] = np.asarray(Ds[ds == d].mean(axis=0)).ravel()
        diff = P_g - P_s
        delta = np.median(diff, axis=0)
        cons = (diff > 0).mean(axis=0)

        passing = (det_goal_all >= MIN_DET_GOAL) & (cons >= MIN_CONSISTENCY) & (delta >= MIN_DELTA)
        n_pass = int(passing.sum())
        order = np.argsort(-delta)
        ranked = [int(i) for i in order if passing[i]]
        pc_local = ranked[:N_PC]
        gd_local = ranked[N_PC:N_PC + N_GD]

        def pack(loc):
            return [{"symbol": str(sym[idx_el[i]]), "ensembl": str(ens[idx_el[i]]),
                     "delta": float(delta[i]), "consistency": float(cons[i]),
                     "detection_goal": float(det_goal_all[i])} for i in loc]

        c["K"] = {"n_passing_filters": n_pass, "K": len(pc_local), "pass": len(pc_local) >= 20,
                  "filters": {"detection_goal_min": MIN_DET_GOAL,
                              "consistency_min": MIN_CONSISTENCY, "delta_min": MIN_DELTA},
                  "positive_control": pack(pc_local),
                  "goal_defining_proxy": pack(gd_local)}

        # ---- F-3 摂動できる遺伝子 ----
        Ds_full = D[s]
        cap_counts = np.zeros(len(ens), dtype=int)
        for j in range(len(ens)):
            col = Ds_full[:, j]
            if col.nnz == 0:
                continue
            cap_counts[j] = capped_cells(np.asarray(col.todense()).ravel() > 0, ds, MAX_CELLS_PER_DONOR)
        n_ok = int(((cap_counts >= MIN_CELLS_PER_GENE) & in_tok).sum())

        def gene_row(symbol):
            i = sym2i.get(symbol)
            if i is None:
                return {"symbol": symbol, "found": False}
            return {"symbol": symbol, "found": True, "ensembl": str(ens[i]),
                    "in_geneformer_dict": bool(in_tok[i]),
                    "n_cells_detected_start": int(np.asarray(Ds_full[:, i].todense()).ravel().__gt__(0).sum()),
                    "n_donors_start": int(len(set(ds[np.asarray(Ds_full[:, i].todense()).ravel() > 0]))),
                    "n_cells_after_cap": int(cap_counts[i]),
                    "pass_floor": bool(cap_counts[i] >= MIN_CELLS_PER_GENE and in_tok[i])}

        pc_ens = {p["ensembl"] for p in c["K"]["positive_control"]}
        t1 = [gene_row(x) for x in T1]; t2 = [gene_row(x) for x in T2]
        for row in t1 + t2:
            row["in_positive_control"] = bool(row.get("ensembl") in pc_ens)
            row["usable"] = bool(row.get("pass_floor") and not row["in_positive_control"]
                                 and row.get("ensembl") not in BLINDED)
        c["F3"] = {"min_cells_per_gene": MIN_CELLS_PER_GENE,
                   "max_cells_per_donor": MAX_CELLS_PER_DONOR,
                   "n_genes_perturbable_whole_universe": n_ok,
                   "T1": t1, "T2": t2,
                   "T1_usable": [r["symbol"] for r in t1 if r["usable"]],
                   "T2_usable": [r["symbol"] for r in t2 if r["usable"]],
                   "pass": bool(n_ok > 0)}

        # ---- F-4 循環性 ----
        gd_ens = {p["ensembl"] for p in c["K"]["goal_defining_proxy"]}
        tgt = {r["ensembl"] for r in t1 + t2 if r.get("usable")}
        dropped = [r["symbol"] for r in t1 + t2 if r.get("in_positive_control")]
        c["F4"] = {
            "goal_defined_by": "著者注釈 subclass.l2（Lake 2023）。遺伝子は使っていない",
            "n_genes_defining_goal": 0,
            "overlap_defining_x_positive_control": 0,
            "overlap_defining_x_targets": 0,
            "overlap_positive_control_x_targets": len(pc_ens & tgt),
            "overlap_triple": 0,
            "targets_dropped_because_in_positive_control": dropped,
            "goal_defining_proxy_x_positive_control": len(gd_ens & pc_ens),
            "pass": bool(len(pc_ens & tgt) == 0)}

        # ---- F-5 深さ ----
        med_g, med_s = float(np.median(n_genes[g])), float(np.median(n_genes[s]))
        ratio = max(med_g, med_s) / min(med_g, med_s)
        mask = g | s
        auc = depth_auc(np.column_stack([n_genes[mask], n_counts[mask]]),
                        g[mask].astype(int), donors[mask])
        c["F5"] = {"median_n_genes_goal": med_g, "median_n_genes_start": med_s,
                   "ratio": ratio, "ratio_pass": bool(ratio <= 1.5),
                   "median_n_counts_goal": float(np.median(n_counts[g])),
                   "median_n_counts_start": float(np.median(n_counts[s])),
                   "depth_only_auc": auc,
                   "auc_pass": bool(auc["auc_pooled_oof"] < 0.75),
                   "pass": bool(ratio <= 1.5 and auc["auc_pooled_oof"] < 0.75)}

        c["all_pass"] = bool(c["K"]["pass"] and c["F1"]["pass"] and c["F2"]["pass"]
                             and c["F3"]["pass"] and c["F4"]["pass"] and c["F5"]["pass"])
        report["candidates"].append(c)
        print(f"[{name}] K={c['K']['K']}(通過 {n_pass}) F1={c['F1']['pass']} F2={c['F2']['pass']} "
              f"F3={n_ok} F4={c['F4']['pass']} F5={c['F5']['pass']} → {c['all_pass']}", flush=True)
        if c["all_pass"]:
            report["adopted"] = name
            break

    (OUT / "candidate_counts.json").write_text(json.dumps(report, ensure_ascii=False, indent=1) + "\n")
    print("→", OUT / "candidate_counts.json")


if __name__ == "__main__":
    sys.exit(main())
