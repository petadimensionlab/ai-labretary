#!/usr/bin/env python3
"""EXP-20260929-001 事前チェック P-B の追加項目。

run_isp.py --precheck が出す提供者分割 AUC に加えて、次の 2 つを出す（HYP-20260929-01 §9-1）。
  (1) 並べ替え検定（ラベルを無作為に入れ替えたときの AUC の分布）。**p < 0.01 を要求する**
  (2) **深さ 2 変数（検出遺伝子数・総カウント）だけ**の提供者分割 AUC
      … 生カウントでの実測は 0.6025（HYP §2.2）。**深さをそろえた後の細胞**で測り直す

出力: work/precheck_auc_extra_003.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path("/Users/petadimensionlab/workspace/research/ai-labretary")
E = ROOT / "experiments/EXP-20260929-001"
SEED = 61
N_PERM = 200
FOLDS = 5


def cv_auc(X, y, g, rng=None, C=1.0, max_iter=5000):
    aucs = []
    for tr, te in GroupKFold(n_splits=FOLDS).split(X, y, groups=g):
        if len(set(y[te])) < 2:
            continue
        clf = make_pipeline(StandardScaler(), LogisticRegression(C=C, max_iter=max_iter))
        clf.fit(X[tr], y[tr])
        aucs.append(roc_auc_score(y[te], clf.predict_proba(X[te])[:, 1]))
    return float(np.mean(aucs)), aucs


def main() -> None:
    from datasets import load_from_disk
    emb = pd.read_csv(E / "results-precheck/precheck/cell_embs.csv")
    key, dkey = "fib_org_state", "donor_id"
    X = emb.drop(columns=[key, dkey, "row_id"]).to_numpy(dtype=np.float64)
    y = (emb[key] == "FIB_orgHigh").astype(int).to_numpy()
    g = emb[dkey].astype(str).to_numpy()
    rows = emb["row_id"].astype(int).to_numpy()

    # 埋め込みでの AUC（run_isp の出力と一致することを確かめる）
    auc_emb, folds_emb = cv_auc(X, y, g)

    # 深さ 2 変数だけ（**深さをそろえた後の同じ細胞**で測る）
    ds = load_from_disk(str(E / "work/tokenized.dataset"))
    spec = json.loads((E / f"work/state_rows-seed{SEED}.json").read_text())
    keep = sorted(spec["rows"]["FIB_orgHigh"] + spec["rows"]["FIB_orgLow"])
    sub = ds.select(keep)
    ng = np.array(sub["n_genes"], float)
    nc = np.array(sub["n_counts"], float)
    ln = np.array(sub["length"], float)
    st = np.array([str(x) for x in sub[key]])
    dn = np.array([str(x) for x in sub[dkey]])
    # cell_embs.csv の row_id は「状態重心用 dataset（= sorted(keep)）の中での行番号」である
    # （run_isp.state_emb_dataset_with_rowid が with_indices で付ける）。位置ではなく ID で結合する。
    order = rows
    assert (st[order] == emb[key].to_numpy()).all(), "行番号で結合できていない"
    Xd = np.column_stack([np.log(ng[order]), np.log(nc[order])])
    auc_depth, folds_depth = cv_auc(Xd, y, g)
    Xl = np.log(ln[order])[:, None]
    auc_len, _ = cv_auc(Xl, y, g)

    # 並べ替え検定（提供者の構造を保つため、提供者ごとにラベルを入れ替える）
    rng = np.random.default_rng(20260930)
    perm = []
    for _ in range(N_PERM):
        yp = y.copy()
        for d in np.unique(g):
            m = g == d
            yp[m] = rng.permutation(y[m])
        perm.append(cv_auc(X, yp, g, max_iter=1000)[0])
    perm = np.array(perm)
    p = float((perm >= auc_emb).sum() + 1) / (N_PERM + 1)

    res = {"seed": SEED, "n_cells": {"FIB_orgHigh": int(y.sum()), "FIB_orgLow": int((1 - y).sum())},
           "auc_embedding_group_by_donor": round(auc_emb, 4),
           "fold_aucs": [round(a, 4) for a in folds_emb],
           "threshold_F6": 0.65, "pass_F6": bool(auc_emb >= 0.65),
           "auc_depth_two_vars": round(auc_depth, 4),
           "auc_depth_fold_aucs": [round(a, 4) for a in folds_depth],
           "auc_token_length_only": round(auc_len, 4),
           "auc_depth_raw_counts_reference": 0.6025,
           "permutation": {"n": N_PERM, "mean": round(float(perm.mean()), 4),
                           "sd": round(float(perm.std(ddof=1)), 4),
                           "max": round(float(perm.max()), 4), "p_one_sided": round(p, 5),
                           "pass_p_lt_0.01": bool(p < 0.01)}}
    (E / "work/precheck_auc_extra_003.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False) + "\n")
    print(json.dumps(res, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
