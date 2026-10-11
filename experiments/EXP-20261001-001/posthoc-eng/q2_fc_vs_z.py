#!/usr/bin/env python3
"""全 1,217 遺伝子について「目標の状態で上がっているか（log2FC）」と、
スクリーンの z（目標の向きへの射影）の関係を調べる。摂動は回さない。"""
import numpy as np, pandas as pd, h5py, json
from pathlib import Path
from scipy.sparse import csr_matrix
from scipy.stats import spearmanr, pearsonr

E = Path("/Users/petadimensionlab/workspace/research/ai-labretary")
OUT = E/"experiments/EXP-20261001-001/posthoc-eng"
f = h5py.File(E/"experiments/EXP-20260929-001/work/input.h5ad","r")
_g=f["var/gene_name"]; _c=np.array([x.decode() if isinstance(x,bytes) else x for x in _g["categories"][:]])
gn=_c[_g["codes"][:]]
_e=f["var/ensembl_id"]; ens=np.array([x.decode() if isinstance(x,bytes) else x for x in _e["values"][:]])
sc=f["obs/subclass_l2"]; cats=np.array([x.decode() if isinstance(x,bytes) else x for x in sc["categories"][:]])
lab=cats[sc["codes"][:]]
dn=f["obs/donor_id"]; dcats=np.array([x.decode() if isinstance(x,bytes) else x for x in dn["categories"][:]])
dlab=dcats[dn["codes"][:]]

M = {"start": np.isin(lab,["FIB","M-FIB"]), "goal": np.isin(lab,["aFIB","MYOF"]),
     "dFIB": lab=="dFIB", "dMFIB": lab=="dM-FIB"}
n_obs=len(lab); G=len(gn)
X=f["X"]; indptr=X["indptr"][:]
sums={k:np.zeros(G) for k in M}; npos={k:np.zeros(G) for k in M}; tot={k:0.0 for k in M}; ncell={k:int(v.sum()) for k,v in M.items()}
B=2000
for s in range(0,n_obs,B):
    e=min(s+B,n_obs); lo,hi=indptr[s],indptr[e]
    m=csr_matrix((X["data"][lo:hi],X["indices"][lo:hi],indptr[s:e+1]-lo),shape=(e-s,G))
    rs=np.asarray(m.sum(axis=1)).ravel()
    cp=m.multiply(1e4/np.maximum(rs,1)[:,None]).tocsr()
    for k,mask in M.items():
        sub=mask[s:e]
        if sub.any():
            sums[k]+=np.asarray(cp[sub].sum(axis=0)).ravel()
            npos[k]+=np.asarray((m[sub]>0).sum(axis=0)).ravel()
mean={k:sums[k]/ncell[k] for k in M}
frac={k:npos[k]/ncell[k] for k in M}

ex=pd.DataFrame({"ensembl":ens,"symbol":gn,
                 "mean_start":mean["start"],"mean_goal":mean["goal"],
                 "frac_start":frac["start"],"frac_goal":frac["goal"],
                 "mean_dFIB":mean["dFIB"],"mean_dMFIB":mean["dMFIB"]})
ex=ex[~ex.ensembl.isin(["NaN","nan",""])].drop_duplicates("ensembl")
ex["log2FC"]=np.log2((ex.mean_goal+0.1)/(ex.mean_start+0.1))
ex["dfrac"]=ex.frac_goal-ex.frac_start

sr=pd.read_csv(E/"outputs/screen/2026-10-03-全遺伝子ランキング.csv")
d=sr.merge(ex,on="ensembl",how="left",suffixes=("","_h5"))
print("merged", d.shape, "missing expr", d.log2FC.isna().sum())
d.to_csv(OUT/"q2_fc_vs_z.csv",index=False)

res={}
for sub,name in [(d,"all_1217"),(d[d.role.isin(["negative","extra_control"])],"controls_1000"),
                 (d[d.role=="intervention"],"intervention_119"),(d[d.role=="positive_control"],"posctl_71")]:
    s=sub.dropna(subset=["log2FC","mean"])
    rho,p=spearmanr(s.log2FC,s["mean"]); r,pp=pearsonr(s.log2FC,s["mean"])
    rho2,p2=spearmanr(s.dfrac,s["mean"])
    res[name]={"n":int(len(s)),"spearman_log2FC_vs_z":round(float(rho),3),"p":float(p),
               "pearson":round(float(r),3),"spearman_dfrac_vs_z":round(float(rho2),3),"p_dfrac":float(p2)}
print(json.dumps(res,indent=1))
json.dump(res,open(OUT/"q2_corr.json","w"),indent=1)

key=["ENG","CCN2","SERPINE1","CA12","HNF1B","RBMS3","PLAT","QSOX1","ABCB8","SEMA5A",
     "MADCAM1","KCNJ1","SLC4A9","SLC13A3","COL1A1","COL1A2","FN1","ELN","GSN","NOX4","LTBP1","CXCL13","CCL21"]
pd.set_option("display.width",250)
k=d[d.symbol.isin(key)][["symbol","role","mean","z81","z82","z83","log2FC","dfrac","mean_start","mean_goal","frac_start","frac_goal"]]
k=k.copy(); k["z_pct"]=k["mean"].rank(pct=True)*0  # placeholder
k["log2FC_pct"]=(d.log2FC.rank(pct=True)[k.index]*100).round(1)
k["z_pct"]=(d["mean"].rank(pct=True)[k.index]*100).round(1)
print(k.sort_values("mean",ascending=False).round(3).to_string(index=False))
