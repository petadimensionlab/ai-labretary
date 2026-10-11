#!/usr/bin/env python3
"""Q4: 歪みの予算を ENG 単独（および逆向き 5 遺伝子）に当てる。摂動は回さない。"""
import sys, json
import numpy as np, pandas as pd
from pathlib import Path
E=Path("/Users/petadimensionlab/workspace/research/ai-labretary")
sys.path.insert(0,str(E/"scripts/analysis"))
from bias_budget import bias_budget, CANONICAL_COVARIATES
OUT=E/"experiments/EXP-20261001-001/posthoc-eng"

ens2sym=pd.read_csv(E/"outputs/screen/2026-10-03-全遺伝子ランキング.csv").set_index("ensembl")
fc=pd.read_csv(OUT/"q2_fc_vs_z.csv").set_index("ensembl")[["log2FC","mean_start"]]
fc["logexpr"]=np.log2(fc.mean_start+0.01)

tables={}
for s in (81,82,83):
    t=pd.read_csv(E/f"experiments/EXP-20261002-002/results/gene_seed_007-{s}.csv")
    t=t.join(fc,on="gene")
    t["symbol"]=t.gene.map(ens2sym.symbol); t["role"]=t.gene.map(ens2sym.role)
    tables[s]=t
ctl=[g for g,r in ens2sym.role.items() if r in ("negative","extra_control")]
print("対照",len(ctl))

TARGETS={"ENG":["ENSG00000106991"],"RBMS3":["ENSG00000144642"],"PLAT":["ENSG00000104368"],
         "CA12":["ENSG00000074410"],"HNF1B":["ENSG00000275410"],"QSOX1":None,"CCN2":None,
         "逆向き5件(ENG,RBMS3,PLAT,CA12,HNF1B)":None}
sym2ens={v:k for k,v in ens2sym.symbol.items()}
TARGETS={k:(v if v else [sym2ens[x] for x in (k.split("(")[1].rstrip(")").split(",") if "(" in k else [k])]) for k,v in TARGETS.items()}

res={}
for name,genes in TARGETS.items():
    per={}
    for s,t in tables.items():
        r=bias_budget(t,genes,ctl)                        # 事前登録の 5 候補
        r2=bias_budget(t,genes,ctl,covariates=list(CANONICAL_COVARIATES)+["log2FC","logexpr"])
        per[s]={"事前登録5候補":{k:r[k] for k in("B_joint","se_B","upper_joint","pass","sum_abs_marginal")},
                "＋発現(探索的)":{k:r2[k] for k in("B_joint","se_B","upper_joint","pass","sum_abs_marginal")},
                "内訳_事前登録":{v:r["covariates"][v]["歪み b_v（同時）"] for v in CANONICAL_COVARIATES},
                "内訳_探索的":{v:r2["covariates"][v]["歪み b_v（同時）"] for v in r2["covariates"]}}
    z=[float(tables[s].set_index("gene").loc[genes,"z"].mean()) for s in (81,82,83)]
    res[name]={"遺伝子数":len(genes),"z_seed81_82_83":[round(x,3) for x in z],"per_seed":per}
    print(f"\n### {name}  z={[round(x,2) for x in z]}")
    for s in (81,82,83):
        a=per[s]["事前登録5候補"]; b=per[s]["＋発現(探索的)"]
        print(f"  seed{s}: 事前登録 B={a['B_joint']:+.3f} 上側={a['upper_joint']:.3f} 合格={a['pass']}"
              f" | ＋発現 B={b['B_joint']:+.3f} 上側={b['upper_joint']:.3f}")
    print("  内訳(seed81, 事前登録):",per[81]["内訳_事前登録"])
    print("  内訳(seed81, 探索的):",per[81]["内訳_探索的"])
json.dump(res,open(OUT/"q4_bias_budget.json","w"),ensure_ascii=False,indent=1)
