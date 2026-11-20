#!/usr/bin/env python3
import json
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from upsetplot import UpSet, from_indicators

BASE = Path('/public3/labmember/qianqh/soil_amp_analysis')
OLD = BASE/'AMP_consensus_out'/'tables'/'amp_predictions_integrated.tsv.gz'
ASCAN = BASE/'AMPScannerV2_out'/'AMPScannerV2_predictions.csv'
OUT = BASE/'AMP_four_model_consensus_out'
FIG, TAB, FAA, SCRIPTS = OUT/'figures', OUT/'tables', OUT/'fasta', OUT/'scripts'
for p in [OUT, FIG, TAB, FAA, SCRIPTS]: p.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(OLD, sep='\t')
asc = pd.read_csv(ASCAN).rename(columns={
    'SeqID':'sequence_id', 'Sequence':'sequence',
    'Prediction_Class':'ampscanner_prediction',
    'Prediction_Probability':'ampscanner_probability'})
asc['ampscanner_pred'] = asc.ampscanner_prediction.eq('AMP')
if len(asc) != 18454 or asc.sequence_id.nunique()!=len(asc) or asc.sequence.nunique()!=len(asc):
    raise RuntimeError('AMP Scanner output is incomplete or duplicated')
df = df.merge(asc[['sequence_id','sequence','ampscanner_pred','ampscanner_probability']],
              on=['sequence_id','sequence'], how='left', validate='one_to_one')

pred_cols = ['esm2_pred','macrel_pred','amplify_pred','ampscanner_pred']
names = ['ESM2-LGBM','Macrel','AMPlify','AMP Scanner v2']
for c in pred_cols:
    if df[c].isna().any(): raise RuntimeError(f'Missing predictions in {c}')
    if df[c].dtype != bool:
        df[c] = df[c].astype(str).str.lower().isin(['true','1','amp'])
df['support_count_4'] = df[pred_cols].sum(axis=1).astype(int)
df['support_tier'] = df.support_count_4.map({4:'Tier 1: 4/4',3:'Tier 2: 3/4',2:'Tier 3: 2/4',1:'Tier 4: 1/4',0:'Not predicted'})

def pattern(row):
    active = [n for n,c in zip(names,pred_cols) if row[c]]
    return ' + '.join(active) if active else 'None'
df['exact_pattern'] = df.apply(pattern, axis=1)

df.to_csv(TAB/'amp_predictions_four_models.tsv.gz', sep='\t', index=False, compression='gzip')
for k in [1,2,3,4]:
    df[df.support_count_4>=k].to_csv(TAB/f'candidates_at_least_{k}_of_4.tsv', sep='\t', index=False)
conservative = df[(df.support_count_4>=2) & (df.esm2_pred | df.macrel_pred)]
conservative.to_csv(TAB/'conservative_consensus_esm_or_macrel_plus.tsv', sep='\t', index=False)
tiers = df.support_tier.value_counts().reindex(['Tier 1: 4/4','Tier 2: 3/4','Tier 3: 2/4','Tier 4: 1/4','Not predicted'], fill_value=0)
tiers.rename_axis('tier').rename('n_sequences').to_csv(TAB/'candidate_tier_counts.tsv', sep='\t')
patterns = df.exact_pattern.value_counts().rename_axis('exact_pattern').rename('n_sequences').reset_index()
patterns.to_csv(TAB/'exact_intersection_counts.tsv', sep='\t', index=False)

pair_rows=[]
for (na,ca),(nb,cb) in combinations(zip(names,pred_cols),2):
    inter=int((df[ca]&df[cb]).sum()); union=int((df[ca]|df[cb]).sum())
    pair_rows.append({'tool_a':na,'tool_b':nb,'intersection':inter,'union':union,
                      'jaccard':inter/union,'overlap_coefficient':inter/min(df[ca].sum(),df[cb].sum())})
pd.DataFrame(pair_rows).to_csv(TAB/'pairwise_overlap_metrics.tsv',sep='\t',index=False)

def write_fasta(sub, path):
    with path.open('w') as h:
        for r in sub.itertuples(): h.write(f'>{r.sequence_id}\n{r.sequence}\n')
for k in [1,2,3,4]: write_fasta(df[df.support_count_4>=k], FAA/f'candidates_at_least_{k}_of_4.faa')
write_fasta(conservative, FAA/'conservative_consensus_esm_or_macrel_plus.faa')
for tier in range(1,5): write_fasta(df[df.support_count_4==tier], FAA/f'exactly_{tier}_of_4.faa')

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':14,'axes.labelsize':11})
colors=['#4C78A8','#59A14F','#E15759','#B279A2']
def save(fig,name):
    fig.savefig(FIG/f'{name}.png',dpi=600,bbox_inches='tight')
    fig.savefig(FIG/f'{name}.pdf',bbox_inches='tight')
    fig.savefig(FIG/f'{name}.svg',bbox_inches='tight')
    plt.close(fig)

# Main intersection figure: union only, avoiding the dominating all-negative set.
ind = df.loc[df.support_count_4>0,pred_cols].rename(columns=dict(zip(pred_cols,names)))
updata = from_indicators(names, ind)
fig=plt.figure(figsize=(11,6.8)); UpSet(updata,subset_size='count',show_counts=True,
    sort_by='cardinality',sort_categories_by='input',element_size=33).plot(fig=fig)
fig.suptitle('Four-model intersection structure of predicted AMPs',y=1.01,fontsize=16)
save(fig,'figure1_upset_four_models')

# Positive calls.
counts=[int(df[c].sum()) for c in pred_cols]
fig,ax=plt.subplots(figsize=(7.5,4.8)); bars=ax.bar(names,counts,color=colors)
ax.set_ylabel('Predicted AMP sequences'); ax.set_title('Positive calls by AMP predictor')
ax.spines[['top','right']].set_visible(False); ax.grid(axis='y',alpha=.25); ax.set_axisbelow(True)
for b,v in zip(bars,counts): ax.text(b.get_x()+b.get_width()/2,v+max(counts)*.015,f'{v:,}\n({v/len(df):.1%})',ha='center',va='bottom')
ax.set_ylim(0,max(counts)*1.16); save(fig,'figure2_positive_counts_four_models')

# Jaccard matrix.
jac=np.eye(4)
for i,j in combinations(range(4),2):
    inter=(df[pred_cols[i]]&df[pred_cols[j]]).sum(); union=(df[pred_cols[i]]|df[pred_cols[j]]).sum()
    jac[i,j]=jac[j,i]=inter/union
fig,ax=plt.subplots(figsize=(6.8,5.6)); im=ax.imshow(jac,cmap='Blues',vmin=0,vmax=1)
ax.set_xticks(range(4),names,rotation=30,ha='right'); ax.set_yticks(range(4),names)
for i in range(4):
    for j in range(4): ax.text(j,i,f'{jac[i,j]:.3f}',ha='center',va='center',color='white' if jac[i,j]>.55 else 'black')
fig.colorbar(im,ax=ax,label='Jaccard index'); ax.set_title('Pairwise agreement of positive sets'); save(fig,'figure3_jaccard_four_models')

# Candidate confidence tiers.
order=['Tier 1: 4/4','Tier 2: 3/4','Tier 3: 2/4','Tier 4: 1/4']
vals=tiers[order].values
fig,ax=plt.subplots(figsize=(7.6,4.8)); bars=ax.bar(order,vals,color=['#2F6B4F','#5B8E6E','#D19A3E','#9E9E9E'])
ax.set_ylabel('Number of sequences'); ax.set_title('Cross-model support tiers'); ax.spines[['top','right']].set_visible(False)
ax.grid(axis='y',alpha=.25); ax.set_axisbelow(True)
for b,v in zip(bars,vals): ax.text(b.get_x()+b.get_width()/2,v+max(vals)*.015,f'{v:,}',ha='center')
ax.set_ylim(0,max(vals)*1.13); save(fig,'figure4_candidate_support_tiers')

# Exact positive patterns, excluding all-negative; log scale preserves rare consensus groups.
p=patterns[patterns.exact_pattern!='None'].sort_values('n_sequences')
fig,ax=plt.subplots(figsize=(10,7.2)); bars=ax.barh(p.exact_pattern,p.n_sequences,color='#4C78A8')
ax.set_xscale('log'); ax.set_xlabel('Number of sequences (log scale)'); ax.set_title('Exact four-model positive-call combinations')
ax.spines[['top','right']].set_visible(False); ax.grid(axis='x',which='both',alpha=.25); ax.set_axisbelow(True)
for b,v in zip(bars,p.n_sequences): ax.text(v*1.08,b.get_y()+b.get_height()/2,f'{v:,}',va='center')
save(fig,'figure5_exact_combinations_four_models')

summary={'input_sequences':len(df),'positive_counts':dict(zip(names,counts)),
         'union_any_model':int((df.support_count_4>=1).sum()),
         'at_least_2_of_4':int((df.support_count_4>=2).sum()),
         'at_least_3_of_4':int((df.support_count_4>=3).sum()),
         'all_4_models':int((df.support_count_4==4).sum()),
         'conservative_esm_or_macrel_plus_another':len(conservative),
         'ampscanner_model':'021820_FULL_MODEL.h5','threshold':0.5}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2)); print('\n',patterns.to_string(index=False))
