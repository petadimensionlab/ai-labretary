#!/usr/bin/env python3
"""Q1: 目標の状態（aFIB+MYOF）で、文献が指標にしている線維化マーカーが実際に高いか。
読み取りのみ。摂動は回さない。"""
import json, sys
import numpy as np, h5py
from pathlib import Path
from scipy.sparse import csr_matrix

E = Path("/Users/petadimensionlab/workspace/research/ai-labretary")
H = E / "experiments/EXP-20260929-001/work/input.h5ad"
OUT = E / "experiments/EXP-20261001-001/posthoc-eng"

GENES = ["ACTA2","CCN2","CTGF","SERPINE1","ENG","COL1A1","COL1A2","COL3A1","FN1","POSTN",
         "TAGLN","PDGFRB","PDGFRA","DCN","LUM","MGP","ELN","NOX4","TGFB1","VIM","MYH11",
         "CA12","HNF1B","KCNJ1","SLC4A9","SLC13A3","RBMS3","PLAT","QSOX1","ABCB8","SEMA5A",
         "MADCAM1","CXCL13","CCL21","GSN","LTBP1","ITGA11","SULF1","ZEB1","CDKN1A"]

f = h5py.File(H, "r")
_gc = f["var/gene_name"]; _cats = np.array([x.decode() if isinstance(x,bytes) else x for x in _gc["categories"][:]]); gn = _cats[_gc["codes"][:]]
# subclass_l2 categorical
sc = f["obs/subclass_l2"]
cats = np.array([x.decode() if isinstance(x,bytes) else x for x in sc["categories"][:]])
codes = sc["codes"][:]
lab = cats[codes]
donor = f["obs/donor_id"]
dcats = np.array([x.decode() if isinstance(x,bytes) else x for x in donor["categories"][:]])
dlab = dcats[donor["codes"][:]]

START = {"FIB","M-FIB"}; GOAL = {"aFIB","MYOF"}
groups = {"start_FIB_MFIB": np.isin(lab, list(START)),
          "goal_aStr": np.isin(lab, list(GOAL)),
          "aFIB": lab=="aFIB", "MYOF": lab=="MYOF",
          "FIB": lab=="FIB", "M-FIB": lab=="M-FIB",
          "dFIB": lab=="dFIB", "dMFIB": lab=="dM-FIB"}
# find exact label spellings
print("labels:", sorted(set(lab.tolist())), file=sys.stderr)

X = f["X"]; indptr = X["indptr"][:]; n_obs = len(lab)
ncounts = f["obs/n_counts"][:] if "n_counts" in f["obs"] else None

# column-wise extraction: read by row blocks
gi = {g: int(np.flatnonzero(gn==g)[0]) for g in GENES if (gn==g).sum()==1}
missing = [g for g in GENES if g not in gi]
cols = np.array(sorted(gi.values()))
colpos = {c:i for i,c in enumerate(cols)}
dense = np.zeros((n_obs, len(cols)), dtype=np.float32)
tot = np.zeros(n_obs, dtype=np.float64)
B = 2000
for s in range(0, n_obs, B):
    e = min(s+B, n_obs)
    lo, hi = indptr[s], indptr[e]
    idx = X["indices"][lo:hi]; dat = X["data"][lo:hi]; ip = indptr[s:e+1]-lo
    m = csr_matrix((dat, idx, ip), shape=(e-s, len(gn)))
    tot[s:e] = np.asarray(m.sum(axis=1)).ravel()
    dense[s:e,:] = m[:, cols].toarray()

cp10k = dense / np.maximum(tot,1)[:,None] * 1e4
res = {"missing_genes": missing, "n_cells": int(n_obs), "groups": {}}
for gname, mask in groups.items():
    if mask.sum()==0: continue
    res["groups"][gname] = {"n_cells": int(mask.sum()), "n_donors": int(len(set(dlab[mask])))}
rows = []
for g, c in gi.items():
    j = colpos[c]
    r = {"gene": g}
    for gname, mask in groups.items():
        if mask.sum()==0: continue
        v = cp10k[mask, j]; raw = dense[mask, j]
        r[gname+"_mean_cp10k"] = round(float(v.mean()),4)
        r[gname+"_pct_pos"] = round(float((raw>0).mean()*100),2)
    rows.append(r)
res["genes"] = rows
json.dump(res, open(OUT/"q1_expression.json","w"), ensure_ascii=False, indent=1)
print(json.dumps(res["groups"], ensure_ascii=False))
import pandas as pd
d = pd.DataFrame(rows).set_index("gene")
keep=[c for c in d.columns if c.startswith(("start_","goal_","aFIB_","MYOF_","FIB_","M-FIB_"))]
pd.set_option("display.width",250); pd.set_option("display.max_columns",40)
print(d[["start_FIB_MFIB_mean_cp10k","goal_aStr_mean_cp10k","start_FIB_MFIB_pct_pos","goal_aStr_pct_pos"]]
      .assign(log2FC=lambda t: np.log2((t.goal_aStr_mean_cp10k+0.1)/(t.start_FIB_MFIB_mean_cp10k+0.1)))
      .sort_values("log2FC", ascending=False).to_string())
d.to_csv(OUT/"q1_expression.csv")
