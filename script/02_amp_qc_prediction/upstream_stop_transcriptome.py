#!/usr/bin/env python3
"""Transcript-aware upstream in-frame stop analysis for soil smORF AMP sets.

Important: absence of an upstream stop on an assembled transcript is censored,
not biological evidence that no stop exists in the genome.
"""
from pathlib import Path
from collections import defaultdict
import re, json, math

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from scipy.stats import fisher_exact, chi2
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = Path('/public3/labmember/qianqh/soil_amp_analysis')
BASE = ROOT / 'AMP_four_model_consensus_out'
FASTA_DIR = BASE / 'fasta'
TD_ROOT = ROOT / '3_td_orf'
ASM_ROOT = ROOT / '2_assembly_rnaspades'
OUT = BASE / 'upstream_stop_transcriptome'
TABLES = OUT / 'tables'; FIGURES = OUT / 'figures'; SCRIPTS = OUT / 'scripts'
for p in (OUT, TABLES, FIGURES, SCRIPTS): p.mkdir(parents=True, exist_ok=True)

SETS = {
    'Conservative': FASTA_DIR / 'conservative_consensus_esm_or_macrel_plus.faa',
    'Extended': FASTA_DIR / 'candidates_at_least_2_of_4.faa',
}
STOP_CODONS = {'TAA', 'TAG', 'TGA'}

def fasta_ids(path):
    return {r.id for r in SeqIO.parse(path, 'fasta')}

set_ids = {k: fasta_ids(v) for k,v in SETS.items()}
all_ids = set.union(*set_ids.values())

# Build robust prefixes from actual crop/sample directory names. Longest first
# handles sample names that contain underscores or hyphens.
prefix_map = []
for crop_dir in sorted(p for p in TD_ROOT.iterdir() if p.is_dir()):
    for sample_dir in sorted(p for p in crop_dir.iterdir() if p.is_dir()):
        prefix_map.append((f'{crop_dir.name}_{sample_dir.name}_', crop_dir.name, sample_dir.name))
prefix_map.sort(key=lambda x: len(x[0]), reverse=True)

def split_tagged_id(tagged):
    for prefix,crop,sample in prefix_map:
        if tagged.startswith(prefix):
            return crop,sample,tagged[len(prefix):]
    return None,None,None

needed = defaultdict(set)
id_meta = {}
for tagged in sorted(all_ids):
    crop,sample,orf_id = split_tagged_id(tagged)
    if orf_id is None:
        id_meta[tagged] = {'mapping_status':'prefix_not_mapped'}
    else:
        needed[(crop,sample)].add(orf_id)
        id_meta[tagged] = {'mapping_status':'prefix_mapped','crop':crop,'sample':sample,'orf_id':orf_id}

def attrs_dict(raw):
    out={}
    for item in raw.split(';'):
        if '=' in item:
            k,v=item.split('=',1); out[k]=v.strip('"')
    return out

orf_records = {}
for (crop,sample), target_orfs in sorted(needed.items()):
    sample_dir = TD_ROOT/crop/sample
    cds_gff = sample_dir/'nonoverlap.cds.gff'
    td_gff = sample_dir/'transcripts.fasta.transdecoder.gff3'
    transcript_fa = ASM_ROOT/crop/f'rnaspades_{sample}'/'transcripts.fasta'
    target_tx = {x.rsplit('.p',1)[0] for x in target_orfs}

    coords={}
    with open(cds_gff) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip(): continue
            c=line.rstrip().split('\t')
            if len(c)<9 or c[2]!='CDS': continue
            a=attrs_dict(c[8]); oid=a.get('Parent') or a.get('ID','').removeprefix('cds.')
            if oid in target_orfs:
                coords[oid]={'transcript_id':c[0], 'cds_start':int(c[3]), 'cds_end':int(c[4]),
                             'strand':c[6], 'phase':c[7]}

    orf_type={}
    with open(td_gff) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip(): continue
            c=line.rstrip().split('\t')
            if len(c)<9 or c[2] != 'mRNA': continue
            a=attrs_dict(c[8]); oid=a.get('ID','')
            if oid not in target_orfs: continue
            name=a.get('Name','')
            m=re.search(r'ORF type:([^\s(,"]+)',name)
            orf_type[oid]=m.group(1) if m else 'unknown'

    txseq={}
    for rec in SeqIO.parse(transcript_fa,'fasta'):
        if rec.id in target_tx:
            txseq[rec.id]=str(rec.seq).upper()

    for oid in target_orfs:
        rec={'crop':crop,'sample':sample,'orf_id':oid,'source_cds_gff':str(cds_gff),
             'source_transcript_fasta':str(transcript_fa)}
        if oid not in coords:
            rec['analysis_status']='cds_not_found'; orf_records[(crop,sample,oid)]=rec; continue
        rec.update(coords[oid]); rec['orf_type']=orf_type.get(oid,'unknown')
        tid=rec['transcript_id']; seq=txseq.get(tid)
        if seq is None:
            rec['analysis_status']='transcript_not_found'; orf_records[(crop,sample,oid)]=rec; continue
        start,end,strand=rec['cds_start'],rec['cds_end'],rec['strand']
        if strand=='+':
            upstream=seq[:start-1]
            oriented=upstream
            start_codon=seq[start-1:start+2]
        else:
            upstream=seq[end:]
            oriented=str(Seq(upstream).reverse_complement())
            start_codon=str(Seq(seq[end-3:end]).reverse_complement())
        nearest=None
        for i in range(0,len(oriented)-2,3):
            codon=oriented[i:i+3]
            if codon in STOP_CODONS:
                nearest=i+3; break
        rec.update({'analysis_status':'ok','transcript_length':len(seq),
                    'upstream_available_bp':len(oriented),'has_upstream_stop':int(nearest is not None),
                    'nearest_stop_bp':nearest,'start_codon':start_codon,
                    'five_prime_complete':int(rec['orf_type'] in {'complete','3prime_partial'}),
                    'five_prime_partial':int(rec['orf_type'] in {'5prime_partial','internal'})})
        orf_records[(crop,sample,oid)]=rec

rows=[]
for tagged in sorted(all_ids):
    meta=id_meta[tagged].copy(); rec={'candidate_id':tagged}
    rec.update(meta)
    if meta.get('mapping_status')=='prefix_mapped':
        rec.update(orf_records.get((meta['crop'],meta['sample'],meta['orf_id']),{'analysis_status':'record_not_found'}))
    rec['in_conservative']=tagged in set_ids['Conservative']
    rec['in_extended']=tagged in set_ids['Extended']
    if rec.get('analysis_status')=='ok':
        d=rec.get('nearest_stop_bp')
        if pd.isna(d): rec['stop_evidence']='Not observed (right-censored)'
        elif d<=30: rec['stop_evidence']='Stop ≤30 bp'
        elif d<=90: rec['stop_evidence']='Stop 31–90 bp'
        else: rec['stop_evidence']='Stop >90 bp'
    rows.append(rec)
df=pd.DataFrame(rows)
df.to_csv(TABLES/'upstream_stop_candidate_results.tsv.gz',sep='\t',index=False,compression='gzip')

ok=df[df.analysis_status.eq('ok')].copy()
groups={
    'Conservative': ok[ok.in_conservative],
    'Extended': ok[ok.in_extended],
    'Extended-only': ok[ok.in_extended & ~ok.in_conservative],
    'Conservative-only': ok[ok.in_conservative & ~ok.in_extended],
}

def wilson(k,n,z=1.95996398454):
    if n==0: return (np.nan,np.nan)
    p=k/n; den=1+z*z/n; cen=(p+z*z/(2*n))/den
    half=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return cen-half,cen+half

summary=[]
group_input_n={'Conservative':len(set_ids['Conservative']),'Extended':len(set_ids['Extended']),
               'Extended-only':len(set_ids['Extended']-set_ids['Conservative']),
               'Conservative-only':len(set_ids['Conservative']-set_ids['Extended'])}
for name,g in groups.items():
    n=len(g); k=int(g.has_upstream_stop.sum()); lo,hi=wilson(k,n)
    summary.append({'set':name,'input_n':group_input_n[name], 'analyzable_n':n,
                    'stop_observed_n':k,'stop_observed_fraction':k/n if n else np.nan,
                    'stop_observed_ci_low':lo,'stop_observed_ci_high':hi,
                    'five_prime_complete_n':int(g.five_prime_complete.sum()),
                    'five_prime_complete_fraction':g.five_prime_complete.mean() if n else np.nan,
                    'median_upstream_available_bp':g.upstream_available_bp.median() if n else np.nan,
                    'median_observed_stop_bp':g.loc[g.has_upstream_stop.eq(1),'nearest_stop_bp'].median()})
summary=pd.DataFrame(summary)
summary.to_csv(TABLES/'upstream_stop_set_summary.tsv',sep='\t',index=False,float_format='%.6f')

evidence_order=['Stop ≤30 bp','Stop 31–90 bp','Stop >90 bp','Not observed (right-censored)']
comp=[]
for name in ['Conservative','Extended','Extended-only']:
    g=groups[name]
    counts=g.stop_evidence.value_counts()
    for ev in evidence_order:
        comp.append({'set':name,'stop_evidence':ev,'n':int(counts.get(ev,0)),
                     'fraction':float(counts.get(ev,0)/len(g)) if len(g) else np.nan})
pd.DataFrame(comp).to_csv(TABLES/'upstream_stop_evidence_composition.tsv',sep='\t',index=False)

# Fisher test only uses disjoint groups to avoid pseudo-replication.
a=groups['Conservative']; b=groups['Extended-only']
cont=np.array([[int(a.has_upstream_stop.sum()),int((1-a.has_upstream_stop).sum())],
               [int(b.has_upstream_stop.sum()),int((1-b.has_upstream_stop).sum())]])
odds,p_fisher=fisher_exact(cont) if cont.min()>=0 else (np.nan,np.nan)

def logrank(g1,g2):
    # Time is bp scanned; event is first in-frame stop. No-stop observations
    # are right-censored at available transcript length rounded to full codons.
    t1=np.where(g1.has_upstream_stop.eq(1),g1.nearest_stop_bp,(g1.upstream_available_bp//3)*3).astype(float)
    t2=np.where(g2.has_upstream_stop.eq(1),g2.nearest_stop_bp,(g2.upstream_available_bp//3)*3).astype(float)
    e1=g1.has_upstream_stop.to_numpy(int); e2=g2.has_upstream_stop.to_numpy(int)
    times=np.unique(np.r_[t1[e1==1],t2[e2==1]])
    oe=var=0.0
    for t in times:
        n1=(t1>=t).sum(); n2=(t2>=t).sum(); d1=((t1==t)&(e1==1)).sum(); d2=((t2==t)&(e2==1)).sum()
        n=n1+n2; d=d1+d2
        if n<=1: continue
        oe += d1-d*n1/n
        var += n1*n2*d*(n-d)/(n*n*(n-1))
    stat=oe*oe/var if var>0 else np.nan
    return stat,chi2.sf(stat,1) if np.isfinite(stat) else np.nan

lr_stat,lr_p=logrank(a,b)
pd.DataFrame([{'comparison':'Conservative vs Extended-only','fisher_odds_ratio':odds,
               'fisher_p':p_fisher,'logrank_chi2':lr_stat,'logrank_p':lr_p,
               'note':'Disjoint groups; no-stop distances treated as right-censored.'}]).to_csv(
                   TABLES/'upstream_stop_statistical_tests.tsv',sep='\t',index=False)

def km_cdf(g):
    t=np.where(g.has_upstream_stop.eq(1),g.nearest_stop_bp,(g.upstream_available_bp//3)*3).astype(float)
    e=g.has_upstream_stop.to_numpy(int); event_times=np.unique(t[e==1]); surv=1.; xs=[0.]; ys=[0.]
    for et in event_times:
        nr=(t>=et).sum(); d=((t==et)&(e==1)).sum()
        if nr: surv*=1-d/nr
        xs.append(et); ys.append(1-surv)
    return np.array(xs),np.array(ys)

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.linewidth':.8,
                      'pdf.fonttype':42,'ps.fonttype':42})
palette={'Conservative':'#0072B2','Extended':'#D55E00','Extended-only':'#CC79A7'}
fig,axs=plt.subplots(1,3,figsize=(10.2,3.45),constrained_layout=True)

# a: evidence composition
bottom=np.zeros(2); x=np.arange(2); colors=['#009E73','#56B4E9','#E69F00','#BDBDBD']
for ev,col in zip(evidence_order,colors):
    vals=[]
    for name in ['Conservative','Extended']:
        g=groups[name]; vals.append((g.stop_evidence==ev).mean())
    axs[0].bar(x,vals,bottom=bottom,color=col,width=.68,label=ev,edgecolor='white',linewidth=.5)
    bottom+=vals
axs[0].set_xticks(x,['Conservative','Extended']); axs[0].set_ylim(0,1); axs[0].set_ylabel('Fraction of candidates')
axs[0].set_title('Upstream-stop evidence'); axs[0].legend(frameon=False,fontsize=6.8,loc='lower center')

# b: censoring-aware cumulative incidence
for name in ['Conservative','Extended']:
    xs,ys=km_cdf(groups[name]); mask=xs<=300
    axs[1].step(xs[mask],ys[mask],where='post',lw=2,color=palette[name],label=f'{name} (n={len(groups[name]):,})')
axs[1].set(xlim=(0,300),ylim=(0,1),xlabel='Distance scanned upstream (bp)',ylabel='Cumulative probability of in-frame stop',title='Censoring-aware stop distance')
axs[1].legend(frameon=False,fontsize=7,loc='lower right')

# c: independent completeness diagnostic
names=['Conservative','Extended']; vals=[groups[n].five_prime_complete.mean() for n in names]
bars=axs[2].bar(np.arange(2),vals,color=[palette[n] for n in names],width=.62)
axs[2].set_xticks(np.arange(2),names); axs[2].set_ylim(0,1); axs[2].set_ylabel('Fraction of candidates')
axs[2].set_title("TransDecoder 5′ completeness")
for bar,v,n in zip(bars,vals,names):
    axs[2].text(bar.get_x()+bar.get_width()/2,v+.025,f'{v:.1%}\n(n={int(groups[n].five_prime_complete.sum()):,})',ha='center',va='bottom',fontsize=8)

for tag,ax in zip('abc',axs):
    ax.spines[['top','right']].set_visible(False); ax.text(-.16,1.04,tag,transform=ax.transAxes,fontweight='bold',fontsize=13)
for ext in ('png','pdf','svg'):
    fig.savefig(FIGURES/f'upstream_stop_conservative_extended.{ext}',dpi=600 if ext=='png' else None,bbox_inches='tight')
plt.close(fig)

# Recommended main figure: feasibility/QC rather than a genomic-validity claim.
# Conservative and Extended-only are disjoint, preventing nested-set comparison.
qgroups={'Conservative':groups['Conservative'],'Extended-only':groups['Extended-only']}
context_order=['Insufficient (<3 bp)','3–29 bp','30–89 bp','90–299 bp','≥300 bp']
context_bins=[-1,2,29,89,299,np.inf]
context_colors=['#BDBDBD','#A6CEE3','#1F78B4','#FDBF6F','#FF7F00']
orf_order=['complete','3prime_partial','5prime_partial','internal']
orf_labels=['Complete','3′ partial','5′ partial','Internal']
orf_colors=['#0072B2','#56B4E9','#E69F00','#BDBDBD']
qc_rows=[]
for name,g in qgroups.items():
    cats=pd.cut(g.upstream_available_bp,context_bins,labels=context_order)
    for cat in context_order:
        qc_rows.append({'set':name,'measure':'upstream_context','category':cat,
                        'n':int((cats==cat).sum()),'fraction':float((cats==cat).mean())})
    for typ,label in zip(orf_order,orf_labels):
        qc_rows.append({'set':name,'measure':'orf_type','category':label,
                        'n':int((g.orf_type==typ).sum()),'fraction':float((g.orf_type==typ).mean())})
pd.DataFrame(qc_rows).to_csv(TABLES/'upstream_stop_feasibility_qc.tsv',sep='\t',index=False)

fig,axs=plt.subplots(1,3,figsize=(12.2,3.75),constrained_layout=False)
# Reserve a dedicated legend rail at the right so legends never cover data.
fig.subplots_adjust(left=.065,right=.79,bottom=.23,top=.88,wspace=.42)
x=np.arange(2); names=list(qgroups); bottom=np.zeros(2)
for cat,col in zip(context_order,context_colors):
    vals=[]
    for name in names:
        g=qgroups[name]; c=pd.cut(g.upstream_available_bp,context_bins,labels=context_order)
        vals.append(float((c==cat).mean()))
    axs[0].bar(x,vals,bottom=bottom,color=col,width=.44,label=cat,edgecolor='white',linewidth=.5)
    bottom+=vals
axs[0].set_xticks(x,['Conservative','Extended-only']); axs[0].set_ylim(0,1)
axs[0].set_ylabel('Fraction of candidates'); axs[0].set_title('Available upstream sequence')

bottom=np.zeros(2)
for typ,label,col in zip(orf_order,orf_labels,orf_colors):
    vals=[float((qgroups[n].orf_type==typ).mean()) for n in names]
    axs[1].bar(x,vals,bottom=bottom,color=col,width=.44,label=label,edgecolor='white',linewidth=.5)
    bottom+=vals
axs[1].set_xticks(x,['Conservative','Extended-only']); axs[1].set_ylim(0,1)
axs[1].set_ylabel('Fraction of candidates'); axs[1].set_title('TransDecoder ORF status')

metric_labels=['≥3 bp context','5′ complete','Upstream stop observed']
vals=np.array([[float((qgroups[n].upstream_available_bp>=3).mean()),
                float(qgroups[n].five_prime_complete.mean()),
                float(qgroups[n].has_upstream_stop.mean())] for n in names])
xx=np.arange(len(metric_labels)); width=.34
for i,name in enumerate(names):
    bars=axs[2].bar(xx+(i-.5)*width,vals[i],width,color=palette[name],label=f'{name} (n={len(qgroups[name]):,})')
    for bar,v in zip(bars,vals[i]):
        label=f'{v:.1%}' if v>=.005 else f'{v:.2%}'
        axs[2].text(bar.get_x()+bar.get_width()/2,v+.006,label,ha='center',va='bottom',fontsize=7)
axs[2].set_xticks(xx,metric_labels,rotation=18,ha='right'); axs[2].set_ylim(0,.22)
axs[2].set_ylabel('Fraction of candidates'); axs[2].set_title('Evidence recoverability')

# Three aligned legend blocks in the shared right margin.
context_handles=[Patch(facecolor=c,edgecolor='white',label=l) for l,c in zip(context_order,context_colors)]
orf_handles=[Patch(facecolor=c,edgecolor='white',label=l) for l,c in zip(orf_labels,orf_colors)]
set_handles=[Patch(facecolor=palette[n],label=f'{n} (n={len(qgroups[n]):,})') for n in names]
leg1=fig.legend(handles=context_handles,title='Available upstream context',loc='upper left',
                bbox_to_anchor=(.815,.88),frameon=False,fontsize=7.4,title_fontsize=8)
leg2=fig.legend(handles=orf_handles,title='TransDecoder ORF status',loc='upper left',
                bbox_to_anchor=(.815,.55),frameon=False,fontsize=7.4,title_fontsize=8)
leg3=fig.legend(handles=set_handles,title='Candidate set',loc='upper left',
                bbox_to_anchor=(.815,.27),frameon=False,fontsize=7.4,title_fontsize=8)

for tag,ax in zip('abc',axs):
    ax.spines[['top','right']].set_visible(False); ax.text(-.16,1.04,tag,transform=ax.transAxes,fontweight='bold',fontsize=13)
for ext in ('png','pdf','svg'):
    fig.savefig(FIGURES/f'upstream_stop_feasibility_qc.{ext}',dpi=600 if ext=='png' else None,bbox_inches='tight')
plt.close(fig)

# Sensitivity figure with disjoint groups.
fig,axs=plt.subplots(1,2,figsize=(7.0,3.25),constrained_layout=True)
for name in ['Conservative','Extended-only']:
    g=groups[name]; k=int(g.has_upstream_stop.sum()); lo,hi=wilson(k,len(g)); v=k/len(g)
    axs[0].errorbar([name],[v],yerr=[[v-lo],[hi-v]],fmt='o',ms=7,capsize=4,color=palette[name],label=name)
    xs,ys=km_cdf(g); mask=xs<=300
    axs[1].step(xs[mask],ys[mask],where='post',lw=2,color=palette[name],label=f'{name} (n={len(g):,})')
axs[0].set(ylim=(0,1),ylabel='Fraction with observed upstream stop',title='Observed stop proportion')
axs[1].set(xlim=(0,300),ylim=(0,1),xlabel='Distance scanned upstream (bp)',ylabel='Cumulative probability',title='Stop distance (right-censoring retained)')
axs[1].legend(frameon=False,fontsize=7,loc='lower right')
for tag,ax in zip('ab',axs): ax.spines[['top','right']].set_visible(False); ax.text(-.15,1.04,tag,transform=ax.transAxes,fontweight='bold',fontsize=13)
for ext in ('png','pdf','svg'):
    fig.savefig(FIGURES/f'upstream_stop_disjoint_sensitivity.{ext}',dpi=600 if ext=='png' else None,bbox_inches='tight')
plt.close(fig)

audit={
    'conservative_input_n':len(set_ids['Conservative']), 'extended_input_n':len(set_ids['Extended']),
    'overlap_n':len(set_ids['Conservative']&set_ids['Extended']),
    'conservative_only_n':len(set_ids['Conservative']-set_ids['Extended']),
    'extended_only_n':len(set_ids['Extended']-set_ids['Conservative']),
    'union_n':len(all_ids), 'analyzable_union_n':len(ok),
    'status_counts':df.analysis_status.fillna(df.mapping_status).value_counts().to_dict(),
    'interpretation':'Transcript-assembly boundary evidence; absence of a stop is right-censored and not genomic absence.'
}
(OUT/'analysis_metadata.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False))
print(json.dumps(audit,indent=2,ensure_ascii=False))
print(summary.to_string(index=False))
