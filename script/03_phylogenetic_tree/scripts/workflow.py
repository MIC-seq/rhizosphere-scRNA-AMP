#!/usr/bin/env python3
"""Pooled environmental transcript-read composition and reference phylogeny."""
import argparse, bz2, csv, gzip, hashlib, json, os, shutil, subprocess, sys, time
from pathlib import Path
from collections import Counter, defaultdict

ROOT = Path(__file__).resolve().parents[1]
BASE = Path('/public3/labmember/qianqh/soil_amp_analysis')
REF = Path('/public3/labmember/qianqh/ncbi_refseq_new')
ENV = Path('/public/home/qianqh/softwares/micromamba/envs/microbe_analysis')
DB = Path('/public2/labmember/qianqh/sc_microbio/kraken_ncbi_refseq_db_202309')
PY = str(ENV/'bin/python')
THREADS = int(os.environ.get('PHYLO_THREADS', '32'))
R = '/public/home/qianqh/softwares/micromamba/envs/use_R/bin/Rscript'
os.environ['PATH'] = str(ENV/'bin')+':/public/home/qianqh/softpath:'+os.environ['PATH']
os.environ['MAFFT_BINARIES'] = str(ENV/'libexec/mafft')

def dump_tsv(path, rows, fields=None):
    rows=list(rows); fields=fields or list(rows[0])
    with open(path,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)

def run(cmd, log):
    print(time.strftime('%F %T'), 'RUN', ' '.join(map(str,cmd)),flush=True)
    with open(log,'a') as f:
        f.write('\nCOMMAND '+json.dumps(list(map(str,cmd)))+'\n');f.flush()
        env=os.environ.copy()
        if str(cmd[0]).endswith('/kraken2'):
            env['LD_LIBRARY_PATH']=str(ENV/'lib')+':'+env.get('LD_LIBRARY_PATH','')
        subprocess.run(list(map(str,cmd)),stdout=f,stderr=subprocess.STDOUT,check=True,env=env)

def status(stage, **extra):
    (ROOT/'status.json').write_text(json.dumps(dict(stage=stage,time=time.strftime('%F %T'),**extra),indent=2))

def init():
    os.chdir(ROOT)
    for d in ['config','logs','reports','source_data','inputs/genomes','databases','tree','figures']:
        Path(d).mkdir(parents=True,exist_ok=True)
    rows=[]
    for crop in ['soybean','rice','wheat']:
        for raw in sorted(Path('/public4/home/qianqh/sc_data',crop).glob('*_2.fq.gz')):
            sample=raw.name[:-8]; key=crop+'__'+sample
            fq=BASE/'1_star_rrna_out'/crop/'non_rrna'/(sample+'_non_rRNA.R2.fq')
            asm=BASE/'2_assembly_rnaspades'/crop/('rnaspades_'+sample)/'transcripts.fasta'
            tax=BASE/'9_taxonomy_anno'/crop/(sample+'_taxonomy.tsv')
            for p in [fq,asm,tax]:
                if not p.is_file() or p.stat().st_size==0: raise RuntimeError('Missing input '+str(p))
            rows.append(dict(key=key,sample=sample,crop=crop,reads=str(fq),reads_bytes=fq.stat().st_size,assembly=str(asm),annotation=str(tax)))
    assert Counter(r['crop'] for r in rows)=={'soybean':37,'rice':10,'wheat':2}
    dump_tsv('config/samples.tsv',rows)
    (Path('config/analysis.json')).write_text(json.dumps(dict(n_samples=49,kraken_database=str(DB),kraken_confidence=0.1,min_reads=100,min_pooled_fraction=0.0001,domains=['Bacteria','Archaea'],denominator='all reads classified to Bacteria or Archaea, including assignments above species',aggregation='sum read counts across 49 datasets; sequencing-depth weighted',genome_policy='latest non-excluded local RefSeq assembly, reference/representative status then completeness then size',tree='PhyloPhlAn400, DIAMOND, MAFFT, trimAl, IQ-TREE; minimum 100 markers',scope='relative classified transcript-read abundance, not cellular abundance'),indent=2))
    return rows

def classify(rows):
    for i,s in enumerate(rows,1):
        out=Path('reports')/(s['key']+'.report.tsv'); done=out.with_suffix('.done.json')
        fingerprint=dict(input=s['reads'],bytes=s['reads_bytes'],mtime=Path(s['reads']).stat().st_mtime,db=str(DB),db_bytes=(DB/'hash.k2d').stat().st_size,confidence=0.1)
        if out.exists() and done.exists() and json.loads(done.read_text())==fingerprint:continue
        status('classifying',sample=s['key'],index=i,total=49)
        tmp=out.with_suffix('.partial')
        run(['/public/home/qianqh/softpath/kraken2','--db',str(DB),'--threads',THREADS,'--confidence','0.1','--report',str(tmp),'--output','/dev/null',s['reads']], 'logs/'+s['key']+'.kraken.log')
        if not tmp.exists() or tmp.stat().st_size==0:raise RuntimeError('Missing Kraken report')
        tmp.replace(out);done.write_text(json.dumps(fingerprint))

def aggregate(rows):
    total=Counter(); meta={}; long=[]; qc=[]
    for s in rows:
        stack=[]; microbial=0; unclassified=0; classified=0
        with open('reports/'+s['key']+'.report.tsv') as f:
            for line in f:
                parts=line.rstrip('\n').split('\t')
                if len(parts)!=6: raise RuntimeError('Unexpected Kraken report columns')
                pct,clade,direct,rank,taxid,label=parts; n=int(clade); depth=len(label)-len(label.lstrip());name=label.strip()
                while stack and stack[-1][0]>=depth: stack.pop()
                ancestry={r:(t,nm) for _,r,t,nm in stack}
                if taxid in ['2','2157']:microbial+=n
                if taxid=='0':unclassified=n
                if taxid=='1':classified=n
                micro = taxid in ['2','2157'] or any(t in ['2','2157'] for _,_,t,_ in stack)
                if rank=='S' and micro:
                    total[taxid]+=n
                    m=dict(taxid=taxid,species=name,phylum=ancestry.get('P',('','Unclassified'))[1],family=ancestry.get('F',('','Unclassified'))[1],genus=ancestry.get('G',('','Unclassified'))[1],domain='Bacteria' if any(t=='2' for _,_,t,_ in stack) else 'Archaea')
                    if taxid in meta and meta[taxid]!=m:raise RuntimeError('Inconsistent taxonomy '+taxid)
                    meta[taxid]=m;long.append(dict(key=s['key'],crop=s['crop'],taxid=taxid,reads=n))
                stack.append((depth,rank,taxid,name))
        qc.append(dict(key=s['key'],classified=classified,unclassified=unclassified,prokaryotic_reads=microbial))
    denominator=sum(x['prokaryotic_reads'] for x in qc)
    if not denominator:raise RuntimeError('No prokaryotic reads classified')
    if sum(total.values())>denominator:raise RuntimeError('Species counts exceed prokaryotic denominator; check report hierarchy')
    presence=Counter(x['taxid'] for x in long if x['reads']>0)
    pooled=[]
    for taxid,n in sorted(total.items(),key=lambda x:-x[1]):
        reasons=[]
        if n<100:reasons.append('reads_below_100')
        if n/denominator<0.0001:reasons.append('pooled_fraction_below_0.0001')
        pooled.append(dict(**meta[taxid],reads=n,pooled_fraction=n/denominator,pooled_percent=100*n/denominator,samples_detected=presence[taxid],eligible=int(not reasons),exclusion=';'.join(reasons)))
    dump_tsv('source_data/sample_read_qc.tsv',qc);dump_tsv('source_data/sample_species_reads.tsv',long)
    dump_tsv('source_data/all_species.tsv',pooled)
    Path('source_data/denominator.json').write_text(json.dumps(dict(prokaryotic_reads=denominator,total_input_reads=sum(x['classified']+x['unclassified'] for x in qc),species_assigned_reads=sum(total.values())),indent=2))
    return pooled

def select_genomes(pooled):
    wanted={x['taxid']:x for x in pooled if x['eligible']};candidates=defaultdict(list)
    review_file=Path('config/reference_supplement/supplement_audit.json')
    reviewed={r['taxid']:r for r in json.loads(review_file.read_text())} if review_file.exists() else {}
    with open(REF/'assembly_summary.txt') as f:
        header=None
        for line in f:
            if line.startswith('#assembly_accession'):header=line.lstrip('#').rstrip('\n').split('\t');break
        if not header:raise RuntimeError('Missing RefSeq header')
        for v in csv.reader(f,delimiter='\t'):
            a=dict(zip(header,v)); tid=a.get('species_taxid')
            if tid not in wanted:continue
            if a.get('version_status')!='latest' or a.get('excluded_from_refseq','na') not in ['na','']:continue
            fname=a.get('ftp_path','').rsplit('/',1)[-1]+'_genomic.fna.gz'
            fp=REF/'raw'/fname
            if not fp.is_file() or fp.stat().st_size==0:continue
            quality=({'reference genome':0,'representative genome':1}.get(a.get('refseq_category'),2),{'Complete Genome':0,'Chromosome':1,'Scaffold':2,'Contig':3}.get(a.get('assembly_level'),4),-int(a.get('genome_size','0') or 0),a['assembly_accession'])
            candidates[tid].append((quality,a,fp))
    selected=[];audit=[]
    for x in pooled:
        y=dict(x)
        if not x['eligible']:audit.append(y);continue
        reviewed_ref=reviewed.get(x['taxid'])
        if reviewed_ref and reviewed_ref['status']=='excluded_after_reference_review':
            y.update(eligible=0,exclusion=reviewed_ref['reason']);audit.append(y);continue
        if reviewed_ref and reviewed_ref['status']=='verified_reference':
            fp=Path(reviewed_ref['reference_file'])
            if not fp.exists():raise RuntimeError('Verified supplemental genome missing: '+str(fp))
            if hashlib.sha256(fp.read_bytes()).hexdigest()!=reviewed_ref['sha256']:raise RuntimeError('Supplemental genome checksum mismatch')
            a={'assembly_accession':reviewed_ref['assembly_accession'],'assembly_level':reviewed_ref['assembly_level']}
            candidates[x['taxid']]=[((-1,0,0,a['assembly_accession']),a,fp)]
        cand=candidates[x['taxid']]
        if not cand:
            y.update(eligible=0,exclusion='no_eligible_local_reference_genome');audit.append(y);continue
        _,a,fp=sorted(cand,key=lambda v:v[0])[0];tip='TAXID'+x['taxid'];dest=Path('inputs/genomes')/(tip+'.fna')
        if not dest.exists():
            opener=gzip.open if str(fp).endswith('.gz') else open
            with opener(fp,'rb') as f,open(str(dest)+'.partial','wb') as g:shutil.copyfileobj(f,g)
            Path(str(dest)+'.partial').replace(dest)
        y.update(tip=tip,assembly_accession=a['assembly_accession'],assembly_level=a['assembly_level'],reference_file=str(fp));selected.append(y);audit.append(y)
    fields=list(pooled[0])+['tip','assembly_accession','assembly_level','reference_file']
    dump_tsv('source_data/species_selection_audit.tsv',audit,fields)
    dump_tsv('source_data/tree_species.tsv',selected,fields)
    if len(selected)<4:raise RuntimeError('Fewer than four eligible genomes')
    status('references_ready',selected=len(selected),species_detected=len(pooled))

def infer_tree():
    run([PY,str(ENV/'bin/phylophlan_write_config_file'),'-o','config/supermatrix.cfg','-d','a','--db_aa','diamond','--map_dna','diamond','--map_aa','diamond','--msa','mafft','--trim','trimal','--tree1','iqtree','--overwrite'],'logs/config.log')
    import configparser
    c=configparser.ConfigParser();c.read('config/supermatrix.cfg')
    # Keep PhyloPhlAn's filename placeholders while specifying reproducible ML settings.
    c['tree1']['params']='-quiet -nt '+str(THREADS)+' -m LG+G4 -bb 1000 -alrt 1000 -seed 42'
    with open('config/supermatrix.cfg','w') as f:c.write(f)
    status('phylogeny')
    pkg=ENV/'lib/python3.8/site-packages/phylophlan'
    tree_folder=os.environ.get('PHYLO_TREE_FOLDER','tree')
    run([PY,str(ENV/'bin/phylophlan'),'-i','genomes','--input_folder','inputs/genomes','-d','phylophlan','--databases_folder','databases','--submat_folder',str(pkg/'phylophlan_substitution_matrices'),'--submod_folder',str(pkg/'phylophlan_substitution_models'),'--diversity','high','--accurate','--remove_fragmentary_entries','--min_num_markers','100','-f','config/supermatrix.cfg','--nproc',THREADS,'--output_folder',tree_folder,'--verbose'],'logs/phylophlan.log')
    trees=list(Path(tree_folder).rglob('*.treefile'))
    if not trees:trees=[p for p in Path(tree_folder).rglob('*.tre') if p.stat().st_size>0]
    if len(trees)!=1:raise RuntimeError('Expected one final tree, found '+repr(trees))
    shutil.copyfile(trees[0],'tree/species_phylogeny.nwk')
    status('plotting')
    run([R,'scripts/plot_tree.R'],'logs/plot.log')
    status('complete_pending_visual_review')

def main():
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['init','classify','finish','all']);a=p.parse_args()
    rows=init()
    if a.stage=='init':return
    if a.stage in ['all','classify']:classify(rows)
    if a.stage in ['all','finish']:select_genomes(aggregate(rows));infer_tree()

if __name__=='__main__':
    try:main()
    except Exception as e:
        if ROOT.exists():status('failed',error=str(e))
        raise
