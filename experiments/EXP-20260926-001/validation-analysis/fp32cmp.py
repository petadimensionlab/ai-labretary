import pandas as pd, numpy as np, glob
E='experiments/EXP-20260926-001'
SP='experiments/EXP-20260926-001/validation-analysis'
SYM={'ENSG00000124762':'CDKN1A','ENSG00000118785':'SPP1','ENSG00000148346':'LCN2','ENSG00000162692':'VCAM1','ENSG00000113249':'HAVCR1'}
g=pd.read_csv(f'{SP}/gene_seed_tls_high.csv')
rows=[]
for f in sorted(glob.glob(f'{E}/results-recheck-fp32/tls_high/run-*-L-1/genes/*/cells.csv')):
    seed=int(f.split('/run-')[1][0]); gene=f.split('/genes/')[1].split('/')[0]
    c=pd.read_csv(f,dtype={'donor_id':str}); b=pd.read_csv(f'{E}/results/tls_high/run-{seed}-L-1/genes/{gene}/cells.csv',dtype={'donor_id':str})
    same=set(c.cell_uid)==set(b.cell_uid)
    m=c.merge(b,on='cell_uid',suffixes=('_fp32','_bf16'))
    s=g[g.run==f'run-{seed}-L-1']; null=s[s.role=='negative']
    r=dict(sym=SYM.get(gene,gene),seed=seed,n=len(c),same_cells=same)
    for met,fx in [('goal',lambda d,t: d[f'shift_goal_{t}']),('gms',lambda d,t: d[f'shift_goal_{t}']-d[f'shift_start_{t}'])]:
        xf,xb=fx(m,'fp32').mean(),fx(m,'bf16').mean(); nl=null[met].to_numpy()
        r|={f'{met}_bf16':xb,f'{met}_fp32':xf,f'{met}_absdiff_rel':abs(xf-xb)/abs(xb),
            f'{met}_pct_bf16':(nl<xb).mean(),f'{met}_pct_fp32':(nl<xf).mean(),
            f'{met}_q90_bf16':xb>np.quantile(nl,.9),f'{met}_q90_fp32':xf>np.quantile(nl,.9),
            f'{met}_q95_bf16':xb>np.quantile(nl,.95),f'{met}_q95_fp32':xf>np.quantile(nl,.95),
            f'{met}_cell_r':np.corrcoef(fx(m,'fp32'),fx(m,'bf16'))[0,1]}
    rows.append(r)
t=pd.DataFrame(rows); t.to_csv('experiments/EXP-20260926-001/results-recheck-fp32/fp32_vs_bf16.csv',index=False)
pd.set_option('display.width',250)
for met in ['goal','gms']:
    print(t[['sym','seed','n','same_cells']+[c for c in t.columns if c.startswith(met)]].to_string())
