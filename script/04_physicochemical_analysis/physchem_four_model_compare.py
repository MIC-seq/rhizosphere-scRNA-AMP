#!/usr/bin/env python3
import gzip
import json
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from modlamp.descriptors import GlobalDescriptor, PeptideDescriptor
from scipy.stats import gaussian_kde, mannwhitneyu

BASE = Path('/public3/labmember/qianqh/soil_amp_analysis/AMP_four_model_consensus_out')
OUT = BASE/'physchem_comparison'
FIG, TAB, SCRIPT = OUT/'figures', OUT/'tables', OUT/'scripts'
for p in [OUT, FIG, TAB, SCRIPT]: p.mkdir(parents=True, exist_ok=True)

SOURCES = {
    'Conservative consensus': BASE/'fasta/conservative_consensus_esm_or_macrel_plus.faa',
    'Expanded consensus (>=2/4)': BASE/'fasta/candidates_at_least_2_of_4.faa',
    'APD6': Path('/public2/labmember/wyc_rawdata/amp/ref_apd6/apd6_clean.faa'),
    'AMPSphere': Path('/public2/labmember/wyc_rawdata/amp/ref_ampsphere/AMP_fullinfo.fa'),
}
ORDER = list(SOURCES)
COLORS = {'Conservative consensus':'#007C70','Expanded consensus (>=2/4)':'#4C78A8',
          'APD6':'#F28E2B','AMPSphere':'#C23B75'}
VALID = set('ACDEFGHIKLMNPQRSTVWY')
SMALL=set('ACDGNPSTV'); BASIC=set('HKR'); POLAR=set('DEHKNQRST')
AROMATIC=set('FHWY'); CHARGED=set('DEHKR')

def read_fasta(path):
    rows=[]
    for r in SeqIO.parse(str(path),'fasta'):
        s=str(r.seq).upper().rstrip('*')
        if s and set(s)<=VALID: rows.append((r.id,s))
    return rows

def descriptor_block(ids, seqs):
    gd=GlobalDescriptor(seqs); gd.calculate_all(ph=7.0,amide=False)
    g=pd.DataFrame(gd.descriptor,columns=gd.featurenames)
    def pep(scale, moment=False):
        d=PeptideDescriptor(seqs,scale)
        if moment: d.calculate_moment(window=1000,angle=100,modality='max')
        else: d.calculate_global(window=1000,modality='max')
        return d.descriptor[:,0]
    n=np.array([len(s) for s in seqs],dtype=float)
    frac=lambda chars: np.array([sum(a in chars for a in s)/len(s)*100 for s in seqs])
    return pd.DataFrame({
        'sequence_id':ids,'sequence':seqs,'length':g.Length,
        'small_residues_pct':frac(SMALL),'basic_residues_pct':frac(BASIC),
        'polar_residues_pct':frac(POLAR),'aromatic_residues_pct':frac(AROMATIC),
        'charged_residues_pct':frac(CHARGED),'isoelectric_point':g.pI,
        'charge_ph7':g.Charge,'aliphatic_index':g.AliphaticInd,
        'instability_index':g.InstabilityInd,'boman_index':g.BomanInd,
        'hydrophobic_moment':pep('eisenberg',True),'hydrophobic_ratio':g.HydrophRatio,
        'hydrophobicity':pep('eisenberg'),'flexibility':pep('flexibility'),'gravy':pep('gravy')})

def compute_dataset(label,path,batch=10000):
    cache=TAB/f"physchem_features_{label.lower().replace(' ','_').replace('(','').replace(')','').replace('>=','ge').replace('/','of')}.tsv.gz"
    if cache.exists(): return pd.read_csv(cache,sep='\t')
    rec=read_fasta(path); chunks=[]
    for i in range(0,len(rec),batch):
        b=rec[i:i+batch]
        chunks.append(descriptor_block([x[0] for x in b],[x[1] for x in b]))
        print(f'{label}: {min(i+batch,len(rec))}/{len(rec)}',flush=True)
    d=pd.concat(chunks,ignore_index=True); d.insert(0,'dataset',label)
    d.to_csv(cache,sep='\t',index=False,compression='gzip',float_format='%.7g')
    return d

frames=[compute_dataset(k,v) for k,v in SOURCES.items()]
data=pd.concat(frames,ignore_index=True)
features=['length','small_residues_pct','basic_residues_pct','polar_residues_pct',
          'aromatic_residues_pct','charged_residues_pct','isoelectric_point','charge_ph7',
          'aliphatic_index','instability_index','boman_index','hydrophobic_moment',
          'hydrophobic_ratio','hydrophobicity','flexibility','gravy']
titles=['Length (residues)','Small residues (%)','Basic residues (%)','Polar residues (%)',
        'Aromatic residues (%)','Charged residues (%)','Isoelectric point','Charge at pH 7.0',
        'Aliphatic index','Instability index','Boman index','Hydrophobic moment',
        'Hydrophobic ratio','Hydrophobicity','Flexibility','GRAVY']

summary=data.groupby('dataset')[features].agg(['count','mean','std','median',lambda x:x.quantile(.25),lambda x:x.quantile(.75)])
summary.columns=['__'.join([a, {'<lambda_0>':'q25','<lambda_1>':'q75'}.get(b,b)]) for a,b in summary.columns]
summary.to_csv(TAB/'physchem_summary_statistics.tsv',sep='\t')

tests=[]
for f in features:
    for a,b in combinations(ORDER,2):
        x=data.loc[data.dataset==a,f].dropna(); y=data.loc[data.dataset==b,f].dropna()
        stat,p=mannwhitneyu(x,y,alternative='two-sided')
        tests.append({'feature':f,'dataset_a':a,'dataset_b':b,'n_a':len(x),'n_b':len(y),'U':stat,'p_raw':p,
                      'median_a':x.median(),'median_b':y.median()})
tests=pd.DataFrame(tests)
order=np.argsort(tests.p_raw.to_numpy()); m=len(tests); raw=tests.p_raw.to_numpy()[order]
step=np.array([1-(1-p)**(m-i) for i,p in enumerate(raw)])
adj_sorted=np.maximum.accumulate(step); adj=np.empty(m); adj[order]=np.minimum(adj_sorted,1.0)
tests['p_holm_sidak']=adj
tests['significant_0.05']=tests.p_holm_sidak<.05
tests.to_csv(TAB/'physchem_mannwhitney_holm_sidak.tsv',sep='\t',index=False)

# AMPSphere Fig. SI1-style ridgeline densities. KDE is deterministically capped at 100k/group.
rng=np.random.default_rng(20260718)
fig,axes=plt.subplots(4,4,figsize=(15,12)); axes=axes.ravel()
for ax,f,title in zip(axes,features,titles):
    pooled=[]
    clipped={}
    for label in ORDER:
        x=data.loc[data.dataset==label,f].dropna().to_numpy()
        lo,hi=np.quantile(x,[.01,.99]); x=x[(x>=lo)&(x<=hi)]
        if len(x)>100000: x=rng.choice(x,100000,replace=False)
        clipped[label]=x; pooled.extend([x.min(),x.max()])
    xmin,xmax=min(pooled),max(pooled); pad=(xmax-xmin)*.03 if xmax>xmin else 1
    grid=np.linspace(xmin-pad,xmax+pad,300)
    for level,label in enumerate(ORDER[::-1]):
        x=clipped[label]
        try: den=gaussian_kde(x)(grid) if np.std(x)>0 else np.zeros_like(grid)
        except Exception: den=np.zeros_like(grid)
        den=den/den.max()*.78 if den.max()>0 else den
        y=level+den
        ax.fill_between(grid,level,y,color=COLORS[label],alpha=.22,lw=0)
        ax.plot(grid,y,color=COLORS[label],lw=1.5)
        ax.axhline(level,color=COLORS[label],lw=.8,alpha=.8)
    ax.set_title(title,loc='left',fontsize=11); ax.set_yticks([])
    ax.spines[['top','right','left']].set_visible(False); ax.tick_params(axis='x',labelsize=8)
    ax.set_ylim(-.08,3.9)
handles=[plt.Line2D([0],[0],color=COLORS[x],lw=5,alpha=.8,label=x) for x in ORDER]
fig.legend(handles=handles,loc='lower center',ncol=2,frameon=False,bbox_to_anchor=(.5,-.015),fontsize=10)
fig.suptitle('Physicochemical profiles of four-model AMP candidates and reference AMP datasets',fontsize=16,y=.995)
fig.tight_layout(rect=(0,.055,1,.975))
for ext,dpi in [('png',600),('pdf',None),('svg',None)]:
    fig.savefig(FIG/f'physchem_figure1_ridgeline_16_features.{ext}',dpi=dpi,bbox_inches='tight')
plt.close(fig)

# Median/IQR standardized heatmap aids interpretation across unlike units.
med=data.groupby('dataset')[features].median().reindex(ORDER)
z=(med-med.mean(axis=0))/med.std(axis=0,ddof=0)
fig,ax=plt.subplots(figsize=(13,4.2)); im=ax.imshow(z,cmap='RdBu_r',vmin=-2,vmax=2,aspect='auto')
ax.set_yticks(range(len(ORDER)),ORDER); ax.set_xticks(range(len(features)),titles,rotation=45,ha='right')
for i in range(len(ORDER)):
    for j in range(len(features)): ax.text(j,i,f'{z.iloc[i,j]:.1f}',ha='center',va='center',fontsize=7)
fig.colorbar(im,ax=ax,label='Standardized dataset median'); ax.set_title('Relative physicochemical profile by dataset')
fig.tight_layout()
for ext,dpi in [('png',600),('pdf',None),('svg',None)]: fig.savefig(FIG/f'physchem_figure2_standardized_median_heatmap.{ext}',dpi=dpi,bbox_inches='tight')
plt.close(fig)

meta={'datasets':{k:int((data.dataset==k).sum()) for k in ORDER},'features':features,
      'outlier_display_rule':'Within each dataset and feature, values below p1 or above p99 are omitted from density plots only.',
      'kde_max_sample_per_dataset':100000,'random_seed':20260718,'modlamp_version':'4.3.2'}
(OUT/'physchem_analysis_metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
print(json.dumps(meta,indent=2))
