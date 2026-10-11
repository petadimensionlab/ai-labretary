#!/usr/bin/env python3
"""線維芽細胞の細胞数と 12 遺伝子の検出率を data/raw/kidney_tls/*.h5ad の raw.X から数える（読み取り専用）。"""
import json, sys, warnings
from pathlib import Path
import h5py, numpy as np, pandas as pd, scipy.sparse as sp
warnings.filterwarnings("ignore")
ROOT = Path("/Users/petadimensionlab/workspace/research/ai-labretary")
RAW = ROOT / "data/raw/kidney_tls"
OUT = ROOT / "outputs/literature/2026-09-29-scrna-chemokine"
OUT.mkdir(parents=True, exist_ok=True)

PRIMARY = ["CXCL13", "CCL19", "CCL21", "CXCL8", "IL6", "LTB"]
EXTRA = ["LTBR", "RELB", "NFKB2", "TNFSF14", "CCL2", "CXCL12", "IL1B"]
GENES = PRIMARY + EXTRA

def rd(h, c):
    g = h["obs"][c]
    if isinstance(g, h5py.Group):
        cats = np.array(g["categories"].asstr()[...]); codes = g["codes"][...]
        return np.where(codes >= 0, cats[np.clip(codes, 0, None)], "nan")
    v = g[...]
    return np.array(g.asstr()[...]) if v.dtype.kind in "OS" else v

def var_names(h, grp):
    v = h[grp]
    ens = np.array(v["_index"].asstr()[...])
    fn = v["feature_name"]
    if isinstance(fn, h5py.Group):
        sym = np.array(fn["categories"].asstr()[...])[fn["codes"][...]]
    else:
        sym = np.array(fn.asstr()[...])
    return ens, sym

# データセット定義: file, 細胞種の列, 線維芽細胞のラベル, 疾患列, 提供者列
DS = {
 "A1": dict(file="A1_KPMP_v1.5_snRNA_Lake2023", susp="single nucleus", ctcol="subclass.l2",
            fib=["FIB","M-FIB","aFIB","dFIB","dM-FIB"], myof=["MYOF"], dis="disease", donor="donor_id",
            bcol="subclass.l2", b=["B"], pl=["PL"]),
 "A2": dict(file="A2_KPMP_v1.5_scRNA_Lake2023", susp="single cell", ctcol="author_cell_type",
            fib=["FIB","aFIB","dFIB"], myof=["MyoF"], dis="disease", donor="donor_id",
            bcol="author_cell_type", b=["B"], pl=["PL"]),
 "A3": dict(file="A3_RCC_LiR2022", susp="single cell", ctcol="annotation",
            fib=["CYGB+ Fibro","CNN1+ Fibro","Stress Fibro","MHC-II Fibro","MMP Fibro",
                 "NFKBIZ+ Fibro","SMC-like Fibro","Collagen Fibro","CDH19+ Fibro"], myof=[],
            dis="summaryDescription", donor="donor_id", bcol="broad_type", b=["B-cell"], pl=["Plasma"]),
 "A4": dict(file="A4_ccRCC_Krishna2021", susp="single cell", ctcol="author_cluster_name",
            fib=["CD45- Myofibroblast"], myof=[], dis="author_type", donor="donor_id",
            bcol="author_cluster_name", b=["B cell"], pl=[]),
 "A5a": dict(file="A5a_ccRCC_Zhang2021", susp="single cell", ctcol="author_cell_type",
            fib=[], myof=[], dis="disease", donor="donor_id",
            bcol="author_cell_type", b=["Bcell"], pl=["Plasma"]),
 "A5b": dict(file="A5b_normal_Zhang2021", susp="single cell", ctcol="author_cell_type",
            fib=[], myof=[], dis="disease", donor="donor_id",
            bcol="author_cell_type", b=["Bcell"], pl=[]),
}

def run(key):
    d = DS[key]
    res = {"key": key, "file": d["file"], "cell_type_column": d["ctcol"],
           "fibroblast_labels": d["fib"], "myofibroblast_labels": d["myof"]}
    with h5py.File(RAW / f"{d['file']}.h5ad", "r") as h:
        n_cells = h["obs"][list(h["obs"].attrs["column-order"])[0]]
        ct = rd(h, d["ctcol"]); dis = rd(h, d["dis"]); don = rd(h, d["donor"])
        bt = rd(h, d["bcol"]); susp = rd(h, "suspension_type")
        n = len(ct)
        res["n_cells"] = int(n)
        res["suspension_type"] = dict(pd.Series(susp).value_counts().astype(int))
        res["n_donors_total"] = int(len(set(don)))
        res["n_B"] = int(np.isin(bt, d["b"]).sum()); res["n_plasma"] = int(np.isin(bt, d["pl"]).sum())
        res["B_plasma_frac"] = round((res["n_B"] + res["n_plasma"]) / n, 5)
        fib = np.isin(ct, d["fib"]); myof = np.isin(ct, d["myof"]) if d["myof"] else np.zeros(n, bool)
        res["n_fibroblast"] = int(fib.sum()); res["n_myofibroblast"] = int(myof.sum())
        has_raw = "raw" in h
        res["has_raw"] = bool(has_raw)
        grp = "raw/var" if has_raw else "var"; Xg = h["raw/X"] if has_raw else h["X"]
        ens, sym = var_names(h, grp)
        s2i = {}
        for i, s in enumerate(sym): s2i.setdefault(s, i)
        present = [g for g in GENES if g in s2i]
        res["genes_missing_in_var"] = [g for g in GENES if g not in s2i]
        res["gene_ensembl"] = {g: str(ens[s2i[g]]) for g in present}
        cols = np.array([s2i[g] for g in present])
        sel = fib | myof
        if sel.sum() == 0:
            res["note"] = "線維芽細胞ラベルが存在しない"
            return res, None
        indptr = Xg["indptr"][...]
        # 生カウントか（先頭 200 細胞）
        a, b = indptr[0], indptr[min(200, n)]
        v = Xg["data"][a:b][:100000]
        res["counts_integer_like"] = bool(np.allclose(v, np.round(v)))
        blocks, tot, ng, rows = [], [], [], []
        step = 20000
        for s in range(0, n, step):
            e = min(n, s + step)
            if not sel[s:e].any(): continue
            a, b = indptr[s], indptr[e]
            m = sp.csr_matrix((Xg["data"][a:b], Xg["indices"][a:b], indptr[s:e+1]-a), shape=(e-s, len(ens)))
            loc = np.where(sel[s:e])[0]
            m = m[loc]
            tot.append(np.asarray(m.sum(1)).ravel()); ng.append(np.diff(m.indptr))
            blocks.append(m[:, cols].toarray().astype(np.float32))
            rows.append(loc + s)
            print(f"  {key} {e}/{n}", file=sys.stderr, end="\r")
        rows = np.concatenate(rows)
        G = pd.DataFrame(np.vstack(blocks), columns=present)
        G["_total"] = np.concatenate(tot); G["_ngenes"] = np.concatenate(ng)
        G["donor"] = don[rows]; G["disease"] = dis[rows]
        G["celltype"] = ct[rows]; G["is_myof"] = myof[rows]
    return res, G

def summarise(key, res, G):
    if G is None: return res, None, None
    fibG = G[~G["is_myof"]]
    res["fib_median_ngenes"] = float(np.median(fibG["_ngenes"])) if len(fibG) else None
    res["fib_median_counts"] = float(np.median(fibG["_total"])) if len(fibG) else None
    res["n_donors_with_fibroblast"] = int(fibG["donor"].nunique())
    rows = []
    for g in GENES:
        for lab, sub in [("fibroblast", fibG), ("fib+myof", G)]:
            if g not in G.columns:
                rows.append(dict(dataset=key, gene=g, cellset=lab, n_cells=len(sub), n_detected=None,
                                 pct=None, pass120=None, note="var に存在しない")); continue
            det = (sub[g] > 0).values
            nd = int(det.sum()); nc = len(sub)
            nd_donor = int(sub.loc[det, "donor"].nunique())
            rows.append(dict(dataset=key, gene=g, cellset=lab, n_cells=nc, n_detected=nd,
                             pct=round(100*nd/nc, 3) if nc else None, n_donors_detected=nd_donor,
                             pass120=bool(nd >= 120), mean_count=round(float(sub[g].mean()), 4)))
    genes_df = pd.DataFrame(rows)
    # 疾患・提供者ごとの内訳
    br = fibG.groupby(["disease"]).agg(n_fib=("_total", "size"), n_donor=("donor", "nunique"))
    for g in PRIMARY:
        if g in fibG.columns: br[g] = fibG.groupby("disease")[g].apply(lambda s: int((s > 0).sum()))
    br = br.reset_index(); br.insert(0, "dataset", key)
    donor_df = fibG.groupby(["donor", "disease"]).size().rename("n_fib").reset_index()
    donor_df.insert(0, "dataset", key)
    return res, genes_df, (br, donor_df)

if __name__ == "__main__":
    keys = sys.argv[1:] or list(DS)
    allres, allg, allbr, alldon = {}, [], [], []
    for k in keys:
        print(f"== {k}", file=sys.stderr)
        res, G = run(k)
        res, gdf, br = summarise(k, res, G)
        allres[k] = res
        if G is not None:
            G.to_csv(OUT / f"fib_cells_{k}.csv.gz", index=False, compression="gzip")
        if gdf is not None:
            allg.append(gdf); allbr.append(br[0]); alldon.append(br[1])
        (OUT / "counts_summary.json").write_text(json.dumps(allres, indent=1, ensure_ascii=False, default=str))
        if allg:
            pd.concat(allg).to_csv(OUT / "gene_detection.csv", index=False)
            pd.concat(allbr).to_csv(OUT / "fib_by_disease.csv", index=False)
            pd.concat(alldon).to_csv(OUT / "fib_by_donor.csv", index=False)
    print("done", file=sys.stderr)
