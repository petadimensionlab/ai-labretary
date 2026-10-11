#!/usr/bin/env python3
"""Q6: 軸は「線維芽細胞らしさから離れること」を測っていないか。
全 1,217 遺伝子について、腎臓全体（細胞種 16 区分）での発現特異性を測り、z と比べる。
data/raw は読み取りのみ。摂動は回さない。"""
import numpy as np, pandas as pd, h5py, json, time
from pathlib import Path
from scipy.sparse import csr_matrix
from scipy.stats import spearmanr
E=Path("/Users/petadimensionlab/workspace/research/ai-labretary")
OUT=E/"experiments/EXP-20261001-001/posthoc-eng"
RAW=E/"data/raw/kidney_tls/A1_KPMP_v1.5_snRNA_Lake2023.h5ad"

sr=pd.read_csv(E/"outputs/screen/2026-10-03-全遺伝子ランキング.csv")
want=set(sr.ensembl)
f=h5py.File(RAW,"r")
vi=f["var"]
print("var keys",list(vi.keys()))
def cat(o):
    c=np.array([x.decode() if isinstance(x,bytes) else x for x in o["categories"][:]]); return c[o["codes"][:]]
# ensembl id の場所を探す
idx=None
for k in list(vi.keys()):
    o=vi[k]
    try:
        v=cat(o) if isinstance(o,h5py.Group) and "categories" in o else (
            np.array([x.decode() if isinstance(x,bytes) else x for x in o["values"][:]]) if isinstance(o,h5py.Group) else
            np.array([x.decode() if isinstance(x,bytes) else x for x in o[:]]))
    except Exception: continue
    if len(set(v)&want)>500: idx=v; print("ensembl 列:",k,len(set(v)&want)); break
assert idx is not None
G=len(idx)
sub=cat(f["obs"]["subclass.l1"]); cls=cat(f["obs"]["class"])
n_obs=len(sub); types=sorted(set(sub))
tmap={t:i for i,t in enumerate(types)}
tcode=np.array([tmap[x] for x in sub])
cols=np.array(sorted({i for i,e in enumerate(idx) if e in want}))
print("対象遺伝子列",len(cols))
S=np.zeros((len(types),len(cols))); N=np.zeros(len(types))
X=f["X"]; ip=X["indptr"][:]; B=4000; t0=time.time()
for s in range(0,n_obs,B):
    e=min(s+B,n_obs); lo,hi=ip[s],ip[e]
    m=csr_matrix((X["data"][lo:hi],X["indices"][lo:hi],ip[s:e+1]-lo),shape=(e-s,G))
    rs=np.asarray(m.sum(axis=1)).ravel()
    cp=m.multiply(1e4/np.maximum(rs,1)[:,None]).tocsr()[:,cols].toarray()
    tc=tcode[s:e]
    for ti in range(len(types)):
        k=tc==ti
        if k.any(): S[ti]+=cp[k].sum(axis=0); N[ti]+=k.sum()
    if s% 40000==0: print(s,round(time.time()-t0,1),"秒")
mean=S/N[:,None]
ex=pd.DataFrame(mean.T,columns=types); ex["ensembl"]=idx[cols]
ex=ex.groupby("ensembl").mean().reset_index()
ex.to_csv(OUT/"q6_celltype_mean.csv",index=False)

TUB=["PT","TAL","DCT","CNT","PC","IC","DTL","ATL","PEC","PapE"]
d=sr.merge(ex,on="ensembl",how="left")
d["tub"]=d[TUB].max(axis=1); d["fib"]=d["FIB"]; d["vsm"]=d["VSM/P"]
d["尿細管らしさ"]=np.log2((d.tub+0.05)/(d.fib+0.05))
d["非線維芽細胞らしさ"]=np.log2((d[[c for c in types if c!="FIB"]].max(axis=1)+0.05)/(d.fib+0.05))
d["線維芽細胞での発現"]=np.log2(d.fib+0.05)
d.to_csv(OUT/"q6_merged.csv",index=False)
ctl=d[d.role.isin(["negative","extra_control"])]
out={}
for v in ["尿細管らしさ","非線維芽細胞らしさ","線維芽細胞での発現"]:
    a=spearmanr(d[v],d["mean"],nan_policy="omit"); b=spearmanr(ctl[v],ctl["mean"],nan_policy="omit")
    out[v]={"rho_全1217":round(float(a.statistic),3),"p_全":float(a.pvalue),
            "rho_対照1000":round(float(b.statistic),3),"p_対照":float(b.pvalue)}
print(json.dumps(out,ensure_ascii=False,indent=1))
json.dump({"corr":out,"types":types,"n_per_type":dict(zip(types,N.astype(int).tolist()))},
          open(OUT/"q6_corr.json","w"),ensure_ascii=False,indent=1)
pd.set_option("display.width",320); pd.set_option("display.max_columns",40)
key=["ENG","CCN2","SERPINE1","CA12","HNF1B","RBMS3","PLAT","QSOX1","ABCB8","SEMA5A","MADCAM1",
     "KCNJ1","SLC4A9","SLC13A3","COL1A1","ELN","GSN","NOX4","LTBP1","CXCL13","CCL21","ZFHX3","ASRGL1","PIK3C2G","DPYSL3","COL5A3","PDZRN4","ABI3BP","ITGA7","COL6A1"]
print(d[d.symbol.isin(key)][["symbol","role","mean","尿細管らしさ","非線維芽細胞らしさ","線維芽細胞での発現","FIB","PT","TAL","EC","VSM/P"]]
      .sort_values("mean",ascending=False).round(3).to_string(index=False))
# 上位 21 / 下位 6 のまとめ
top=d.nlargest(21,"mean"); bot=d.nsmallest(6,"mean")
for nm,g_ in [("上位21",top),("下位6",bot),("対照1000",ctl),("全1217",d)]:
    print(f"{nm}: 尿細管らしさ中央値={g_['尿細管らしさ'].median():+.3f}  非線維芽細胞らしさ={g_['非線維芽細胞らしさ'].median():+.3f}  線維芽細胞発現={g_['線維芽細胞での発現'].median():+.3f}")
