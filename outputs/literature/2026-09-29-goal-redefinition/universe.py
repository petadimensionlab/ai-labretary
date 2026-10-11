#!/usr/bin/env python3
"""採用候補（深さの十分位内でマーカー得点 上位 5%）を固定し、
(1) 入力長 4096 の切り捨てを考慮した実行可能な遺伝子数、
(2) 出発点で 120 細胞（提供者上限 20）を満たす遺伝子の総数（対照の母集団）を数える。"""
import json
from pathlib import Path
import h5py, numpy as np, pandas as pd, scipy.sparse as sp
import pickle

ROOT = Path("/Users/petadimensionlab/workspace/research/ai-labretary")
RAW = ROOT / "data/raw/kidney_tls/A1_KPMP_v1.5_snRNA_Lake2023.h5ad"
OUT = ROOT / "outputs/literature/2026-09-29-goal-redefinition"
TOKDICT = Path("/Users/petadimensionlab/workspace/research/Geneformer/geneformer_hf/"
               "geneformer/token_dictionary_gc104M.pkl")
CAP, MIN, MAXLEN = 20, 120, 4096

D = pd.read_csv(OUT / "fib_A1_panel.csv.gz")
ORG7 = ["PDPN", "NGFR", "ICAM1", "VCAM1", "TNFSF13B", "MADCAM1", "ENPP2"]
CHEMO3 = ["CXCL13", "CCL19", "CCL21"]
cp = D[ORG7].values / D["n_counts"].values[:, None] * 1e4
SCORE = np.log1p(cp).sum(1)


def strat_goal(score, top_frac, nbin=10):
    dq = pd.qcut(D["n_genes"].rank(method="first"), nbin, labels=False).values
    g = np.zeros(len(D), bool)
    for b in range(nbin):
        idx = np.where(dq == b)[0]
        k = max(1, int(round(top_frac * len(idx))))
        g[idx[np.lexsort((idx, -score[idx]))[:k]]] = True
    return g


def main():
    goal = strat_goal(SCORE, 0.05)
    start = ~goal
    tok = pickle.loads(TOKDICT.read_bytes())
    res = {"goal_top_frac": 0.05, "n_goal": int(goal.sum()), "n_start": int(start.sum())}

    # 切り捨てを考慮した版（入力長の上限 4096 を超える細胞を出発点から外した場合）
    safe = start & (D["n_genes"].values <= MAXLEN)
    res["n_start_len_le_4096"] = int(safe.sum())
    donors = D["donor"].values

    def cap_n(mask, col):
        s = np.where(mask & (col > 0))[0]
        if not len(s):
            return 0
        per = pd.Series(donors[s]).value_counts()
        return int(np.minimum(per, CAP).sum())

    res["chemo_with_len_filter"] = {}
    for g in CHEMO3 + ["LTB", "TNFSF14", "CXCL16", "CXCL1"]:
        res["chemo_with_len_filter"][g] = {
            "cap20_start": cap_n(start, D[g].values),
            "cap20_start_len_le_4096": cap_n(safe, D[g].values),
        }

    # 全遺伝子での母集団（Geneformer の辞書にある遺伝子のみ）
    rows = D["row"].values
    with h5py.File(RAW, "r") as h:
        v = h["raw/var"]
        ens = np.array(v["_index"].asstr()[...])
        fn = v["feature_name"]
        sym = (np.array(fn["categories"].asstr()[...])[fn["codes"][...]]
               if isinstance(fn, h5py.Group) else np.array(fn.asstr()[...]))
        in_dict = np.array([str(e) in tok for e in ens])
        Xg = h["raw/X"]; indptr = Xg["indptr"][...]
        n = len(indptr) - 1
        # 出発点だけを行として集め、遺伝子ごとに提供者上限つきで数える
        keep_rows = set(rows[start].tolist())
        don_of = dict(zip(rows.tolist(), donors.tolist()))
        # 提供者 -> インデックス
        dlist = sorted(set(donors.tolist()))
        dix = {d: i for i, d in enumerate(dlist)}
        C = np.zeros((len(dlist), len(ens)), np.int32)
        step = 20000
        for s in range(0, n, step):
            e = min(n, s + step)
            sel = np.array([i in keep_rows for i in range(s, e)])
            if not sel.any():
                continue
            a, b = indptr[s], indptr[e]
            m = sp.csr_matrix((Xg["data"][a:b], Xg["indices"][a:b],
                               indptr[s:e + 1] - a), shape=(e - s, len(ens)))
            loc = np.where(sel)[0]
            m = (m[loc] > 0).astype(np.int32)
            dd = np.array([dix[don_of[i + s]] for i in loc])
            for k in range(len(dlist)):
                w = np.where(dd == k)[0]
                if len(w):
                    C[k] += np.asarray(m[w].sum(0)).ravel()
            print(f"  {e}/{n}", end="\r")
    capped_tot = np.minimum(C, CAP).sum(0)
    ok = (capped_tot >= MIN) & in_dict
    res["universe"] = {
        "n_genes_in_var": int(len(ens)),
        "n_genes_in_geneformer_dict": int(in_dict.sum()),
        "n_genes_pass120_cap20_in_start": int(ok.sum()),
        "n_genes_pass120_and_in_dict": int(ok.sum()),
    }
    pd.DataFrame({"ensembl": ens[ok], "symbol": sym[ok],
                  "n_cells_cap20": capped_tot[ok]}).sort_values(
        "n_cells_cap20", ascending=False).to_csv(OUT / "start_gene_universe.csv", index=False)
    (OUT / "universe.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print("\n", json.dumps(res, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
