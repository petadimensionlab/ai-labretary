#!/usr/bin/env python3
"""Q3: ENG の変位はどこへ向かうのか。全 717 遺伝子（3 seed）の変位を
   いくつかの向きに分解する。保存済みの埋め込みだけを読む。摂動は回さない。"""
import numpy as np, pandas as pd, json
from pathlib import Path
E=Path("/Users/petadimensionlab/workspace/research/ai-labretary")
R=E/"experiments/EXP-20261001-001/results"
OUT=E/"experiments/EXP-20261001-001/posthoc-eng"
sr=pd.read_csv(E/"outputs/screen/2026-10-03-全遺伝子ランキング.csv").set_index("ensembl")

def unit(v): return v/np.linalg.norm(v)
rows=[]; meand={}
for seed in (81,82,83):
    run=R/f"run-{seed}-L-1"
    sz=np.load(run/"state_emb_vectors.npz")
    g,s,df,dm=[sz[k].astype(np.float64) for k in ("aStr","FIB_MFIB","dFIB","dMFIB")]
    dirs={"目標(start→goal)":unit(unit(g)-unit(s)),
          "dFIB(start→dFIB)":unit(unit(df)-unit(s)),
          "dMFIB(start→dMFIB)":unit(unit(dm)-unit(s)),
          "出発点の重心そのもの":unit(s)}
    genes=sorted(p.name for p in (run/"genes").iterdir() if p.is_dir())
    D={}
    for gene in genes:
        o=np.load(run/"genes"/gene/"emb/orig.npy").astype(np.float64)
        p=np.load(run/"genes"/gene/"emb/pert.npy").astype(np.float64)
        d=p-o; n=np.linalg.norm(d,axis=1); ok=n>0
        D[gene]=d[ok].mean(axis=0)
        r={"ensembl":gene,"seed":seed,"n":int(ok.sum()),"dnorm_med":float(np.median(n[ok]))}
        for k,u in dirs.items(): r["cos_"+k]=float(np.median((d[ok]@u)/n[ok]))
        # 重心からの距離の変化（コサインでなく実距離）
        for k,v in [("goal",g),("start",s),("dFIB",df),("dMFIB",dm)]:
            r["dist_"+k]=float(np.median(np.linalg.norm(p[ok]-v,axis=1)-np.linalg.norm(o[ok]-v,axis=1)))
        rows.append(r)
    meand[seed]=D
    print("seed",seed,"done",len(genes))

t=pd.DataFrame(rows); t["symbol"]=t.ensembl.map(sr.symbol); t["role"]=t.ensembl.map(sr.role)

# 全遺伝子に共通の変位方向 w（摂動がどの遺伝子でも生む共通の動き）
for seed in (81,82,83):
    M=np.array([meand[seed][k] for k in sorted(meand[seed])])
    w=unit(M.mean(axis=0))
    gl=sorted(meand[seed])
    cw={gl[i]:float(np.dot(unit(M[i]),w)) for i in range(len(gl))}
    t.loc[t.seed==seed,"cos_共通方向w"]=t.loc[t.seed==seed,"ensembl"].map(cw)
    # w 自身が目標の向きとどれだけ重なるか
    run=R/f"run-{seed}-L-1"; sz=np.load(run/"state_emb_vectors.npz")
    ug=unit(unit(sz["aStr"].astype(np.float64))-unit(sz["FIB_MFIB"].astype(np.float64)))
    print(f"seed{seed}: 共通方向 w と目標の向きのコサイン = {float(np.dot(w,ug)):+.4f}")

t.to_csv(OUT/"q3_geometry.csv",index=False)
agg=t.groupby(["ensembl","symbol","role"]).median(numeric_only=True).reset_index().drop(columns=["seed"])
agg["z"]=agg.ensembl.map(sr["mean"])
agg.to_csv(OUT/"q3_geometry_median.csv",index=False)
pd.set_option("display.width",300); pd.set_option("display.max_columns",40)
key=["ENG","PLAT","CCL21","CXCL13","MADCAM1","COL1A1","ELN","NOX4","GSN","LTBP1","SERPINE1","TGFB1","VIM","CDKN1A","LTBR"]
print("\n=== 代表遺伝子（3 seed の中央値）===")
print(agg[agg.symbol.isin(key)].sort_values("z",ascending=False).round(4).to_string(index=False))
print("\n=== 全 717 の要約 ===")
print(agg[[c for c in agg.columns if c.startswith(("cos_","dist_"))]+["z"]].describe().round(4).to_string())
from scipy.stats import spearmanr
print("\n=== z との順位相関（717 個）===")
for c in [c for c in agg.columns if c.startswith(("cos_","dist_"))]:
    r=spearmanr(agg[c],agg.z,nan_policy="omit"); print(f"  {c:28s} rho={r.statistic:+.3f}  p={r.pvalue:.2e}")
