#!/usr/bin/env python3
"""h5ad を読み取り専用（backed）で開き、Geneformer 投入可否を確認する。
usage: inspect_h5ad.py FILE.h5ad [...] > report.json"""
import json, pickle, sys, warnings
import anndata as ad
import numpy as np
warnings.filterwarnings("ignore")
# scripts/geneformer/ が本物の geneformer パッケージを隠すので、スクリプト自身のディレクトリを外す
sys.path = [p for p in sys.path if not p.endswith("/scripts")]

def gf_vocab():
    try:
        from geneformer import TOKEN_DICTIONARY_FILE
        with open(TOKEN_DICTIONARY_FILE, "rb") as f:
            return {k for k in pickle.load(f) if str(k).startswith("ENS")}
    except Exception:
        return None

KEYS = ("disease", "condition", "cell_type", "celltype", "author_cell_type", "tissue", "region", "sample",
        "age", "development_stage", "sex", "donor", "status", "stage", "location", "compartment", "summary")
vocab = gf_vocab()
out = {}
for f in sys.argv[1:]:
    a = ad.read_h5ad(f, backed="r")
    r = {"shape": list(a.shape), "has_raw": a.raw is not None, "obs_columns": list(a.obs.columns),
         "var_names_head": list(a.var_names[:3])}
    ens = [g for g in a.var_names if str(g).startswith("ENSG")]
    r["ensembl_frac"] = round(len(ens) / a.n_vars, 3)
    if vocab is not None:
        r["gf_vocab_overlap"] = len(set(ens) & vocab)
    # 生カウントか: raw.X（なければ X）の先頭 200 細胞が整数か
    X = (a.raw.X if a.raw is not None else a.X)[:200]
    v = X.data if hasattr(X, "data") else np.asarray(X).ravel()
    v = v[:100000]
    r["counts_integer_like"] = bool(np.allclose(v, np.round(v))) if v.size else None
    r["labels"] = {}
    for c in a.obs.columns:
        if any(k in c.lower() for k in KEYS):
            vc = a.obs[c].astype(str).value_counts()
            r["labels"][c] = dict(vc.head(30).astype(int)) | ({"__n_levels": int(len(vc))} if len(vc) > 30 else {})
    out[f] = r
    a.file.close()
print(json.dumps(out, indent=1, ensure_ascii=False, default=int))
