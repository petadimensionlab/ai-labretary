#!/usr/bin/env python3
"""Q4: 交絡の確認。発現量・元の順位・入力の長さ・移動量・一様成分と z の関係。
ENG がそれぞれの対照分布のどこにいるかを出す。摂動は回さない。"""
import numpy as np, pandas as pd, json
from pathlib import Path
from scipy.stats import spearmanr

E=Path("/Users/petadimensionlab/workspace/research/ai-labretary")
OUT=E/"experiments/EXP-20261001-001/posthoc-eng"

# --- 遺伝子×seed の交絡量（1001-001: 717 個 / 1002-002: 500 個）
a=pd.read_csv(E/"experiments/EXP-20261001-001/results/eval/gene_seed_005.csv")
a=a[a.col=="m_goal"][["gene","seed","D_raw","lp_med","ll_med","relpos_med","dnorm_med","n_cells"]]
au=pd.read_csv(E/"experiments/EXP-20261001-001/results/eval/gene_seed_005.csv")
au=au[au.col=="m_uniform"][["gene","seed","D_raw"]].rename(columns={"D_raw":"uniform"})
a=a.merge(au,on=["gene","seed"])
bs=[]
for s in (81,82,83):
    b=pd.read_csv(E/f"experiments/EXP-20261002-002/results/gene_seed_007-{s}.csv")
    bs.append(b)
b=pd.concat(bs,ignore_index=True)
print("1002-002 columns:",list(b.columns))
b=b[["gene","seed","dnorm_med","lp_med","ll_med","relpos_med","uniform","n_cells"]]
common=set(a.gene)&set(b.gene); print("重複遺伝子",len(common))
cov=pd.concat([a[["gene","seed","lp_med","ll_med","relpos_med","dnorm_med","uniform","n_cells"]],
               b[~b.gene.isin(common)]],ignore_index=True)
cov=cov.groupby("gene").median(numeric_only=True).reset_index()

d=pd.read_csv(OUT/"q2_fc_vs_z.csv").merge(cov,left_on="ensembl",right_on="gene",how="left")
print("交絡量が付いた遺伝子", d.lp_med.notna().sum(),"/",len(d))

d["expr_start"]=d.mean_start; d["logexpr"]=np.log2(d.mean_start+0.01)
ctl=d[d.role.isin(["negative","extra_control"])].dropna(subset=["lp_med"])
rows=[]
for v in ["logexpr","frac_start","lp_med","ll_med","relpos_med","dnorm_med","uniform","log2FC","dfrac"]:
    r_all=spearmanr(d[v],d["mean"],nan_policy="omit")
    r_c=spearmanr(ctl[v],ctl["mean"],nan_policy="omit")
    rows.append({"変数":v,"rho_全1217":round(float(r_all.statistic),3),"p_全":float(r_all.pvalue),
                 "rho_対照1000":round(float(r_c.statistic),3),"p_対照":float(r_c.pvalue)})
t=pd.DataFrame(rows); print(t.to_string(index=False))
t.to_csv(OUT/"q4_corr_table.csv",index=False)

# ENG の位置（対照分布の分位）
def pct(col,val): 
    x=ctl[col].dropna().to_numpy(); return round(float((x<val).mean()*100),1)
eng=d[d.symbol=="ENG"].iloc[0]
pos={}
for v in ["logexpr","frac_start","lp_med","ll_med","relpos_med","dnorm_med","uniform","log2FC"]:
    pos[v]={"ENG":round(float(eng[v]),4),"対照内分位_%":pct(v,eng[v]),
            "対照_中央値":round(float(ctl[v].median()),4),
            "対照_5%":round(float(ctl[v].quantile(.05)),4),"対照_95%":round(float(ctl[v].quantile(.95)),4)}
print(json.dumps(pos,ensure_ascii=False,indent=1))
json.dump(pos,open(OUT/"q4_eng_position.json","w"),ensure_ascii=False,indent=1)
d.to_csv(OUT/"q4_merged.csv",index=False)
