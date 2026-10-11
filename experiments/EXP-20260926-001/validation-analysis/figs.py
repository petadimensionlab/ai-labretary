import pandas as pd, numpy as np, matplotlib, os
matplotlib.use('Agg'); import matplotlib.pyplot as plt
E='experiments/EXP-20260926-001'; F=f'{E}/figures'
SP='experiments/EXP-20260926-001/validation-analysis'
SYM={'ENSG00000141510':'TP53','ENSG00000124762':'CDKN1A','ENSG00000111321':'LTBR','ENSG00000118785':'SPP1','ENSG00000148346':'LCN2','ENSG00000162692':'VCAM1','ENSG00000113249':'HAVCR1','ENSG00000075624':'ACTB','ENSG00000111640':'GAPDH','ENSG00000196262':'PPIA'}
COL={'target':'#2a78d6','positive':'#eb6834','housekeeping':'#1baf7a','negative':'#9a9993'}
plt.rcParams.update({'font.family':['Hiragino Sans','DejaVu Sans'],'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
LAB={'goal':'shift_goal (×1e-4)','gms':'shift_goal − shift_start (×1e-4)'}
for r in ['tls_high','tls_low']:
    g=pd.read_csv(f'{SP}/gene_seed_{r}.csv')
    fig,axes=plt.subplots(2,1,figsize=(8,7.5))
    for ax,m in zip(axes,['goal','gms']):
        order=list(SYM); 
        for k,(run,s) in enumerate(sorted(g.groupby('run'))):
            null=s[s.role=='negative'][m]*1e4
            x0=k*(len(order)+2)
            ax.scatter(np.full(len(null),x0)+np.random.default_rng(k).uniform(-.3,.3,len(null)),null,s=10,color=COL['negative'],label='null（ランダム 20）' if k==0 else None)
            for q,ls in [(.95,'-'),(.80,':')]:
                ax.hlines(np.quantile(null,q),x0-.5,x0+len(order)+.5,color='#52514e',lw=.8,ls=ls,label=(f'null q{int(q*100)}' if k==0 else None))
            for j,gene in enumerate(order):
                v=s[s.gene==gene][m]
                if len(v):
                    role=s[s.gene==gene].role.iloc[0]
                    ax.scatter(x0+1+j,v.iloc[0]*1e4,s=28,color=COL[role],zorder=3,label=({'target':'対象','positive':'陽性対照','housekeeping':'HK'}[role] if k==0 and gene in ('ENSG00000141510','ENSG00000162692','ENSG00000075624') else None))
            ax.text(x0+len(order)/2,ax.get_ylim()[1] if False else None or 0,'',)
        ticks=[];labels=[]
        for k in range(3):
            x0=k*(len(order)+2); ticks+= [x0]+[x0+1+j for j in range(len(order))]; labels+=[f'null s{k}']+[SYM[x] for x in order]
        ax.set_xticks(ticks); ax.set_xticklabels(labels,rotation=90,fontsize=7)
        ax.axhline(0,color='k',lw=.5); ax.set_ylabel(LAB[m])
        ax.set_title(f'{r}: {LAB[m]} — 遺伝子ごとの細胞平均（1 点 = 1 seed, n≈100 細胞）',fontsize=9)
    axes[0].legend(fontsize=7,ncol=6,loc='lower left',frameon=False)
    fig.tight_layout(); fig.savefig(f'{F}/null_vs_genes_{r}.png',dpi=180); plt.close(fig)
# drift figure
fig,axes=plt.subplots(1,2,figsize=(8,3.4),sharey=True)
for ax,r in zip(axes,['tls_high','tls_low']):
    g=pd.read_csv(f'{SP}/gene_seed_{r}.csv')
    roles=['negative','housekeeping','target','positive']
    for i,role in enumerate(roles):
        s=g[g.role==role]
        for j,(m,lab) in enumerate([('goal','goal'),('start','start'),('alt','alt')]):
            y=s[m]*1e4; x=i*4+j
            ax.bar(x,y.mean(),yerr=y.std(),color=COL[role],alpha=[1,.65,.4][j],width=.85,capsize=2)
    ax.set_xticks([i*4+j for i in range(4) for j in range(3)]); ax.set_xticklabels(['goal','start','alt']*4,rotation=90,fontsize=7)
    for i,role in enumerate(roles): ax.text(i*4+1,-3.75,{'negative':'null','housekeeping':'HK','target':'対象','positive':'陽性'}[role],ha='center',fontsize=8)
    ax.set_ylim(-3.9,0.1); ax.axhline(0,color='k',lw=.5); ax.set_title(f'{r}: 3 重心へのシフト（平均±SD, 遺伝子×seed）',fontsize=9)
axes[0].set_ylabel('cos 変化 ×1e-4')
fig.tight_layout(); fig.savefig(f'{F}/deletion_drift.png',dpi=180); plt.close(fig)
# fp32 recheck
p='experiments/EXP-20260926-001/results-recheck-fp32/fp32_vs_bf16.csv'
if os.path.exists(p):
    t=pd.read_csv(p)
    fig,axes=plt.subplots(1,2,figsize=(7,3.3))
    for ax,m in zip(axes,['goal','gms']):
        ax.scatter(t[f'{m}_bf16']*1e4,t[f'{m}_fp32']*1e4,s=20,color='#2a78d6')
        for _,rw in t.iterrows(): ax.annotate(f"{rw.sym}/s{rw.seed}",(rw[f'{m}_bf16']*1e4,rw[f'{m}_fp32']*1e4),fontsize=6,xytext=(2,2),textcoords='offset points')
        lo=min(t[f'{m}_bf16'].min(),t[f'{m}_fp32'].min())*1e4; hi=max(t[f'{m}_bf16'].max(),t[f'{m}_fp32'].max())*1e4
        ax.plot([lo,hi],[lo,hi],color='#52514e',lw=.7,ls=':'); ax.set_xlabel(f'bf16 {LAB[m]}'); ax.set_ylabel('fp32')
    fig.suptitle('fp32 再計算（run 1 の境界遺伝子, 同一細胞・同一 seed）',fontsize=9)
    fig.tight_layout(); fig.savefig(f'{F}/fp32_recheck.png',dpi=180); plt.close(fig)
print(sorted(os.listdir(F)))
