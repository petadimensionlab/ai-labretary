import pandas as pd, numpy as np
from scipy import stats
E='experiments/EXP-20260926-001'
SYM={'ENSG00000141510':'TP53','ENSG00000124762':'CDKN1A','ENSG00000111321':'LTBR','ENSG00000118785':'SPP1','ENSG00000148346':'LCN2','ENSG00000162692':'VCAM1','ENSG00000113249':'HAVCR1','ENSG00000075624':'ACTB','ENSG00000111640':'GAPDH','ENSG00000196262':'PPIA'}
TG=list(SYM)[:5]; PR=list(SYM)[:3]; POS=list(SYM)[5:7]; HK=list(SYM)[7:]
def bh(p):
    p=np.asarray(p);n=len(p);o=np.argsort(p);q=np.empty(n);q[o]=np.minimum.accumulate((p[o]*n/np.arange(1,n+1))[::-1])[::-1];return np.clip(q,0,1)
for r in ['tls_high','tls_low']:
    c=pd.read_csv(f'{E}/results/{r}/cells.csv',dtype={'donor_id':str}); c['gms']=c.shift_goal-c.shift_start
    print('#####',r)
    for m in ['shift_goal','gms']:
        res={g:[] for g in TG+POS+HK}
        donors=sorted(c.donor_id.unique())
        for run,sub in c.groupby('run'):
            for d in [None]+donors:
                s=sub if d is None else sub[sub.donor_id!=d]
                gm=s.groupby(['gene','role'])[m].mean().reset_index()
                null=gm[gm.role=='negative'][m].to_numpy()
                mu,sd=null.mean(),null.std(ddof=1)
                rows={}
                for g in TG+POS+HK:
                    if g not in set(gm.gene): continue
                    x=gm[gm.gene==g][m].iloc[0]
                    rows[g]=dict(pct=(null<x).mean(),q95=x>np.quantile(null,.95),q90=x>np.quantile(null,.9),le80=x<=np.quantile(null,.8),neg=x<0,pn=stats.norm.sf((x-mu)/sd))
                q=bh([rows[g]['pn'] for g in TG]);
                for i,g in enumerate(TG): rows[g]['fdr']=q[i]
                for g,v in rows.items(): res[g].append(dict(run=run,drop=d,**v))
        for g in PR+POS+HK:
            t=pd.DataFrame(res[g]); full=t[t['drop'].isna()]; lo=t[t['drop'].notna()]
            # per-drop: count seeds passing
            agg=lo.groupby('drop').agg(n_pass=('q95',lambda x: 0), )
            npass=lo.assign(pp=lo.q95&(lo.get('fdr',pd.Series(1,index=lo.index))<0.05)).groupby('drop').pp.sum()
            npos=lo.groupby('drop').q90.sum()
            fal_q=lo.groupby('drop').le80.all(); fal_any=lo.assign(f=lo.le80|lo.neg).groupby('drop').f.all()
            print(f"  {m:10s} {SYM[g]:7s} full pct={list(full.pct.round(2))} | LODO pct range per seed: " +
                  ' '.join(f"[{s.pct.min():.2f},{s.pct.max():.2f}]" for _,s in lo.groupby('run')) +
                  f" | primary-pass seeds max={npass.max() if g in PR else '-'} | q90 seeds min..max={npos.min()}..{npos.max()} | falsify(q80 all) {fal_q.mean():.2f} of drops, falsify(q80|neg) {fal_any.mean():.2f}")
