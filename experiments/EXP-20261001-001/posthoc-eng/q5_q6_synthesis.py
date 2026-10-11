#!/usr/bin/env python3
"""Q2 条件づけ・Q5 逆向き 5 件の共通説明・Q6 尿細管遺伝子の濃縮検定。摂動は回さない。"""
import numpy as np, pandas as pd, json
from pathlib import Path
from scipy.stats import spearmanr, fisher_exact, mannwhitneyu
E=Path("/Users/petadimensionlab/workspace/research/ai-labretary")
OUT=E/"experiments/EXP-20261001-001/posthoc-eng"
d4=pd.read_csv(OUT/"q4_merged.csv")
d6=pd.read_csv(OUT/"q6_merged.csv")[["ensembl","尿細管らしさ","非線維芽細胞らしさ","線維芽細胞での発現"]]
g3=pd.read_csv(OUT/"q3_geometry_median.csv").rename(columns={"dnorm_med":"dnorm_med_g3","n":"n_cells_g3"})
d=d4.merge(d6,on="ensembl").merge(g3.drop(columns=["symbol","role","z"]),on="ensembl",how="left")
d.to_csv(OUT/"all_merged.csv",index=False)
ctl=d[d.role.isin(["negative","extra_control"])]
R={}

# --- Q2 条件づけ: log2FC が同じくらいの対照の中で ENG はどれだけ外れているか
lo=ctl[ctl.log2FC<=ctl.log2FC.quantile(.05)]
R["log2FC下位5%の対照"]={"n":int(len(lo)),"z_平均":round(float(lo["mean"].mean()),3),
    "z_標準偏差":round(float(lo["mean"].std()),3),"z_最小":round(float(lo["mean"].min()),3),
    "z_5%":round(float(lo["mean"].quantile(.05)),3)}
X=np.column_stack([np.ones(len(ctl)),ctl.log2FC]); b=np.linalg.lstsq(X,ctl["mean"],rcond=None)[0]
resid=ctl["mean"]-X@b; sd=float(resid.std(ddof=2))
for s in ["ENG","RBMS3","QSOX1","PLAT","CCN2","SERPINE1","CA12","HNF1B","MADCAM1","CCL21","SLC4A9","KCNJ1"]:
    r=d[d.symbol==s]
    if len(r)==0: continue
    r=r.iloc[0]; pred=b[0]+b[1]*r.log2FC
    R.setdefault("log2FC で条件づけた残差",{})[s]={
        "z":round(float(r["mean"]),2),"log2FC":round(float(r.log2FC),3),
        "log2FC から予測される z":round(float(pred),3),"残差(標準偏差単位)":round(float((r["mean"]-pred)/sd),2)}

# --- Q3 まとめ: 目標の重心から実距離で遠ざかった遺伝子
R["目標の重心から実距離で遠ざかった遺伝子"]={
  "n_717中":int((g3.dist_goal>0).sum()),
  "遺伝子":g3[g3.dist_goal>0].sort_values("dist_goal",ascending=False)[["symbol","dist_goal","z"]].round(5).to_dict("records")}
R["ENG_幾何"]={k:round(float(g3[g3.symbol=="ENG"].iloc[0][k]),5) for k in
   ["cos_目標(start→goal)","cos_dFIB(start→dFIB)","cos_dMFIB(start→dMFIB)","cos_出発点の重心そのもの",
    "dist_goal","dist_start","dist_dFIB","dist_dMFIB","cos_共通方向w","dnorm_med_g3"]}
for k in list(R["ENG_幾何"]):
    v=g3[k]; R["ENG_幾何"][k]={"値":R["ENG_幾何"][k],"717個中の順位":int((v<g3[g3.symbol=="ENG"].iloc[0][k]).sum())+1,
                              "717個の中央値":round(float(v.median()),5)}

# --- Q5 逆向き 5 件（＋参考 2 件）の横並び
cols=["symbol","role","mean","log2FC","dfrac","mean_start","mean_goal","尿細管らしさ",
      "線維芽細胞での発現","lp_med","ll_med","dnorm_med","uniform","cos_目標(start→goal)","dist_goal"]
rev=["ENG","CA12","HNF1B","RBMS3","PLAT","QSOX1","CCN2","SERPINE1"]
tb=d[d.symbol.isin(rev)][cols].copy()
for c in ["log2FC","mean_start","尿細管らしさ","線維芽細胞での発現","lp_med","ll_med","dnorm_med","uniform","cos_目標(start→goal)"]:
    tb[c+"_対照内分位%"]=[round(float((ctl[c]<v).mean()*100),1) for v in tb[c]]
pd.set_option("display.width",400); pd.set_option("display.max_columns",60)
print("=== Q5 逆向きの遺伝子の横並び ===")
print(tb[["symbol","role","mean","log2FC","log2FC_対照内分位%","mean_start","mean_start_対照内分位%",
          "尿細管らしさ","尿細管らしさ_対照内分位%","lp_med_対照内分位%","dnorm_med_対照内分位%","uniform_対照内分位%","dist_goal"]]
      .sort_values("mean").round(3).to_string(index=False))
tb.round(4).to_csv(OUT/"q5_reversed_table.csv",index=False)

# --- Q6 濃縮検定: 上位 21 個に尿細管特異的遺伝子が多いか
d=d.sort_values("mean",ascending=False)
top21=d.head(21); 
for thr in (2.0,1.5):
    a=int((top21["尿細管らしさ"]>thr).sum()); bb=21-a
    cc=int((ctl["尿細管らしさ"]>thr).sum()); dd=len(ctl)-cc
    odds,p=fisher_exact([[a,bb],[cc,dd]],alternative="greater")
    R.setdefault("Q6_上位21の尿細管濃縮",{})[f"閾値{thr}"]={
        "上位21中":a,"上位21の割合%":round(a/21*100,1),"対照1000中":cc,"対照の割合%":round(cc/len(ctl)*100,1),
        "オッズ比":round(float(odds),2),"片側p":round(float(p),4)}
u=mannwhitneyu(top21["尿細管らしさ"],ctl["尿細管らしさ"],alternative="greater")
R["Q6_上位21の尿細管らしさ(順位和)"]={"p":round(float(u.pvalue),4),
   "上位21中央値":round(float(top21["尿細管らしさ"].median()),3),"対照中央値":round(float(ctl["尿細管らしさ"].median()),3)}
# 線維芽細胞での発現が低いことの濃縮
u2=mannwhitneyu(top21["線維芽細胞での発現"],ctl["線維芽細胞での発現"],alternative="less")
R["Q6_上位21は線維芽細胞で低発現か"]={"p":round(float(u2.pvalue),4),
   "上位21中央値":round(float(top21["線維芽細胞での発現"].median()),3),"対照中央値":round(float(ctl["線維芽細胞での発現"].median()),3)}

# --- 多変量: z を何が説明するか（対照 1000 個で当てる）
import itertools
V=["log2FC","線維芽細胞での発現","尿細管らしさ","lp_med","ll_med","relpos_med","dnorm_med","uniform"]
cc=ctl.dropna(subset=V+["mean"])
Xq=np.column_stack([np.ones(len(cc))]+[cc[v].rank(pct=True).to_numpy() for v in V])
bb2,res2,_,_=np.linalg.lstsq(Xq,cc["mean"].to_numpy(),rcond=None)
yhat=Xq@bb2; r2=1-((cc["mean"]-yhat)**2).sum()/((cc["mean"]-cc["mean"].mean())**2).sum()
R["対照1000で z を説明する重回帰(分位)"]={"R2":round(float(r2),4),
   **{v:round(float(bb2[i+1]),3) for i,v in enumerate(V)}}
print("\n"+json.dumps(R,ensure_ascii=False,indent=1))
json.dump(R,open(OUT/"q5_q6_synthesis.json","w"),ensure_ascii=False,indent=1)
