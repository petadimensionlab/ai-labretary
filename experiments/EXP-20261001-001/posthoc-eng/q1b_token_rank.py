#!/usr/bin/env python3
"""Geneformer が実際に見ている「入力の中での順位」を、出発点の細胞と目標の細胞で比べる。
過剰発現は『順位を 1 位へ上げる』操作なので、これが摂動の意味を直接表す。摂動は回さない。"""
import numpy as np, pandas as pd, pickle, json
from pathlib import Path
from datasets import load_from_disk
from geneformer import TOKEN_DICTIONARY_FILE
E=Path("/Users/petadimensionlab/workspace/research/ai-labretary")
OUT=E/"experiments/EXP-20261001-001/posthoc-eng"
tokd=pickle.load(open(TOKEN_DICTIONARY_FILE,"rb"))
ds=load_from_disk(str(E/"experiments/EXP-20261001-001/work/tokenized.dataset"))
sr=pd.read_csv(E/"outputs/screen/2026-10-03-全遺伝子ランキング.csv")
s2e=dict(zip(sr.symbol,sr.ensembl))
SYMS=["ENG","CCN2","SERPINE1","CA12","HNF1B","RBMS3","PLAT","QSOX1","SEMA5A","ABCB8",
      "MADCAM1","CXCL13","CCL21","KCNJ1","SLC4A9","SLC13A3","COL1A1","COL1A2","ELN","GSN","NOX4","LTBP1","FN1","ZFHX3","PIK3C2G","ASRGL1"]
tk={}
for s in SYMS:
    e=s2e.get(s)
    if e and e in tokd: tk[s]=tokd[e]
print("token が引けた遺伝子",len(tk))
state=np.array(ds["fib_anno_state"]); sub=np.array(ds["subclass_l2"])
ids=ds["input_ids"]
rank={s:{} for s in tk}
groups={"start":state=="FIB_MFIB","goal":state=="aStr","aFIB":sub=="aFIB","MYOF":sub=="MYOF"}
recs=[]
pos={s:np.full(len(ids),np.nan) for s in tk}
for i,row in enumerate(ids):
    a=np.asarray(row)
    for s,t in tk.items():
        w=np.flatnonzero(a==t)
        if len(w): pos[s][i]=w[0]/len(a)      # 相対順位（0 = 最上位）
for s in tk:
    r={"symbol":s}
    for gname,m in groups.items():
        v=pos[s][m]; r[gname+"_検出率%"]=round(float(np.isfinite(v).mean()*100),2)
        r[gname+"_相対順位中央値"]=round(float(np.nanmedian(v)),4) if np.isfinite(v).any() else np.nan
        r[gname+"_上位10%入り率%"]=round(float(np.nanmean(v<0.10)*100),2)
    r["Δ相対順位(目標-出発点)"]=round(r["goal_相対順位中央値"]-r["start_相対順位中央値"],4)
    r["Δ上位10%入り率"]=round(r["goal_上位10%入り率%"]-r["start_上位10%入り率%"],2)
    recs.append(r)
t=pd.DataFrame(recs).merge(sr[["symbol","mean","role"]],on="symbol")
pd.set_option("display.width",320); pd.set_option("display.max_columns",40)
print(t[["symbol","role","mean","start_検出率%","goal_検出率%","start_相対順位中央値","goal_相対順位中央値",
         "Δ相対順位(目標-出発点)","start_上位10%入り率%","goal_上位10%入り率%","Δ上位10%入り率"]]
      .sort_values("mean",ascending=False).to_string(index=False))
t.to_csv(OUT/"q1b_token_rank.csv",index=False)
from scipy.stats import spearmanr
print("\nΔ上位10%入り率 と z の順位相関:",spearmanr(t["Δ上位10%入り率"],t["mean"]))
print("Δ相対順位 と z の順位相関:",spearmanr(t["Δ相対順位(目標-出発点)"],t["mean"]))
