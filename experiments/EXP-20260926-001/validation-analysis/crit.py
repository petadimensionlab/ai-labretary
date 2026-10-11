import pandas as pd, numpy as np, json
from scipy import stats
E='experiments/EXP-20260926-001'
SYM={'ENSG00000141510':'TP53','ENSG00000124762':'CDKN1A','ENSG00000111321':'LTBR','ENSG00000118785':'SPP1','ENSG00000148346':'LCN2','ENSG00000162692':'VCAM1','ENSG00000113249':'HAVCR1','ENSG00000075624':'ACTB','ENSG00000111640':'GAPDH','ENSG00000196262':'PPIA'}
out={}
for r in ['tls_high','tls_low']:
    c=pd.read_csv(f'{E}/results/{r}/cells.csv',dtype={'donor_id':str})
    alt=[x for x in c.columns if x.startswith('shift_alt')][0]
    c['alt']=c[alt]; c['gms']=c.shift_goal-c.shift_start
    g=c.groupby(['run','gene','role']).agg(goal=('shift_goal','mean'),start=('shift_start','mean'),alt=('alt','mean'),gms=('gms','mean'),n=('cell_uid','size'),nd=('donor_id','nunique'),L=('length','mean')).reset_index()
    print(f'\n##### {r}')
    # (a) drift
    print('(a) per role: mean goal/start/alt (x1e-4), frac gene-seed with all three <0')
    for role,s in g.groupby('role'):
        print(f"  {role:12s} n={len(s):3d} goal={s.goal.mean()*1e4:+.2f}±{s.goal.std()*1e4:.2f} start={s.start.mean()*1e4:+.2f} alt={s.alt.mean()*1e4:+.2f} gms={s.gms.mean()*1e4:+.3f}±{s.gms.std()*1e4:.3f} all3neg={((s.goal<0)&(s.start<0)&(s.alt<0)).mean():.2f}")
    print('  overall gene-seed all3neg', ((g.goal<0)&(g.start<0)&(g.alt<0)).mean(), 'n',len(g), '| cells goal<0', (c.shift_goal<0).mean().round(3),'start<0',(c.shift_start<0).mean().round(3),'alt<0',(c.alt<0).mean().round(3))
    print('  gene-level spearman goal vs start', stats.spearmanr(g.goal,g.start).statistic.round(3), 'goal vs alt',stats.spearmanr(g.goal,g.alt).statistic.round(3), '| cell-level pearson goal vs start', np.corrcoef(c.shift_goal,c.shift_start)[0,1].round(3),'start vs alt',np.corrcoef(c.shift_start,c.alt)[0,1].round(3))
    print('  mean |start-alt| cell', np.abs(c.shift_start-c.alt).mean(), 'mean |goal| cell', np.abs(c.shift_goal).mean(), ' gene-level mean |start-alt|', np.abs(g.start-g.alt).mean(), 'vs |goal-start|',np.abs(g.gms).mean())
    # (b) confounds
    det=json.load(open(f'{E}/work/detection-'+('aPT_CKD_TLShigh' if r=='tls_high' else 'aPT_CKD_TLSlow')+'.json'))
    g['det']=g.gene.map(det)
    for m in ['goal','gms']:
        print(f'(b) {m}: gene-seed spearman vs det {stats.spearmanr(g.det,g[m]).statistic:.3f} vs meanLen {stats.spearmanr(g.L,g[m]).statistic:.3f} vs n_cells {stats.spearmanr(g.n,g[m]).statistic:.3f}')
    c['v_gms']=c.gms
    for m in ['shift_goal','v_gms']:
        rr=[stats.spearmanr(s.length,s[m]).statistic for _,s in c.groupby(['run','gene']) if s.length.nunique()>5]
        print(f'   cell-level within gene-seed spearman({m}, length): median {np.nanmedian(rr):.3f} IQR {np.nanpercentile(rr,25):.3f}..{np.nanpercentile(rr,75):.3f}; frac len==4096 {(c.length==4096).mean():.3f}')
    # (d) small n genes
    print('(d) n cells / donors per seed for named genes')
    for gene in SYM:
        s=g[g.gene==gene]
        if len(s): print(f'   {SYM[gene]:7s} n={list(s.n)} donors={list(s.nd)} det={det.get(gene)} goal={s.goal.mean()*1e4:+.2f}±{s.goal.std()*1e4:.2f} gms={s.gms.mean()*1e4:+.3f}±{s.gms.std()*1e4:.3f}')
    # null stats per seed
    for m in ['goal','gms']:
        ns=g[g.role=='negative'].groupby('run')[m].agg(['mean','std','median',lambda x: np.quantile(x,.95), lambda x:(x<0).mean()])
        ns.columns=['mean','sd','median','q95','frac_neg']
        print(f'   null {m} per seed (x1e-4):'); print((ns*[1e4,1e4,1e4,1e4,1]).round(3).to_string())
    g.to_csv(f'experiments/EXP-20260926-001/validation-analysis/gene_seed_{r}.csv',index=False)
