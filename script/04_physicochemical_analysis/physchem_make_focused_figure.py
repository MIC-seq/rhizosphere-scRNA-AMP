#!/usr/bin/env python3
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

ROOT=Path('/public3/labmember/qianqh/soil_amp_analysis/AMP_four_model_consensus_out/physchem_comparison')
T=ROOT/'tables'; F=ROOT/'figures'
files={
 'Conservative consensus':T/'physchem_features_conservative_consensus.tsv.gz',
 'Expanded consensus (>=2/4)':T/'physchem_features_expanded_consensus_ge2of4.tsv.gz',
 'APD6':T/'physchem_features_apd6.tsv.gz','AMPSphere':T/'physchem_features_ampsphere.tsv.gz'}
order=list(files); data=pd.concat([pd.read_csv(p,sep='\t') for p in files.values()],ignore_index=True)
colors={'Conservative consensus':'#007C70','Expanded consensus (>=2/4)':'#4C78A8','APD6':'#F28E2B','AMPSphere':'#C23B75'}
features=['length','basic_residues_pct','isoelectric_point','charge_ph7','instability_index','boman_index','hydrophobic_moment','gravy']
titles=['Length (residues)','Basic residues (%)','Isoelectric point','Charge at pH 7.0','Instability index','Boman index','Hydrophobic moment','GRAVY']
rng=np.random.default_rng(20260718)
fig,axes=plt.subplots(2,4,figsize=(15,6.8)); axes=axes.ravel()
for ax,f,title in zip(axes,features,titles):
    vals={}; lim=[]
    for label in order:
        x=data.loc[data.dataset==label,f].dropna().to_numpy(); lo,hi=np.quantile(x,[.01,.99]); x=x[(x>=lo)&(x<=hi)]
        if len(x)>100000:x=rng.choice(x,100000,replace=False)
        vals[label]=x; lim += [x.min(),x.max()]
    pad=(max(lim)-min(lim))*.03; grid=np.linspace(min(lim)-pad,max(lim)+pad,300)
    for level,label in enumerate(order[::-1]):
        x=vals[label]; den=gaussian_kde(x)(grid); den=den/den.max()*.78
        ax.fill_between(grid,level,level+den,color=colors[label],alpha=.22,lw=0); ax.plot(grid,level+den,color=colors[label],lw=1.6)
        ax.axhline(level,color=colors[label],lw=.8); med=np.median(x); ax.plot([med,med],[level,level+.22],color=colors[label],lw=2)
    ax.set_title(title,loc='left',fontsize=12); ax.set_yticks([]); ax.spines[['top','right','left']].set_visible(False); ax.set_ylim(-.08,3.9)
handles=[plt.Line2D([0],[0],color=colors[x],lw=5,label=x) for x in order]
fig.legend(handles=handles,loc='lower center',ncol=2,frameon=False,bbox_to_anchor=(.5,-.02))
fig.suptitle('Key physicochemical properties of predicted and reference AMPs',fontsize=16)
fig.tight_layout(rect=(0,.09,1,.95))
for ext,dpi in [('png',600),('pdf',None),('svg',None)]: fig.savefig(F/f'physchem_figure3_focused_8_features.{ext}',dpi=dpi,bbox_inches='tight')
