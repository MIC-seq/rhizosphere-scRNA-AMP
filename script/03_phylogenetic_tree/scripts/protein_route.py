"""Reference-proteome phylogeny with frozen transcript abundance metadata."""
import argparse, concurrent.futures, csv, hashlib, json, os, shutil, subprocess, time
from pathlib import Path
from Bio import SeqIO

ROOT=Path(__file__).resolve().parents[1]
ORIGINAL=Path('/public3/labmember/qianqh/soil_amp_analysis/10_phylogenetic tree')
SOURCE=ORIGINAL/'comparison/supplemented_8refs'
ENV=Path('/public/home/qianqh/softwares/micromamba/envs/microbe_analysis')
PRODIGAL='/public/home/qianqh/softwares/micromamba/envs/metagenome/bin/prodigal'
R='/public/home/qianqh/softwares/micromamba/envs/use_R/bin/Rscript'
os.environ['PATH']=str(ENV/'bin')+':'+os.environ['PATH']
os.environ['MAFFT_BINARIES']=str(ENV/'libexec/mafft')
os.environ['OPENBLAS_NUM_THREADS']='1'
THREADS=int(os.environ.get('PHYLO_THREADS','32'))

def load(path):
    with open(path) as f:return list(csv.DictReader(f,delimiter='\t'))

def write_table(path,rows):
    with open(path,'w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)

def state(stage,**kw):
    p=ROOT/'status.tmp';p.write_text(json.dumps(dict(stage=stage,time=time.strftime('%F %T'),**kw),indent=2));p.replace(ROOT/'status.json')

def run(command,log,cwd,timeout=None):
    print(time.strftime('%F %T'),json.dumps(list(map(str,command))),flush=True)
    with open(log,'a') as f:
        f.write('\n'+json.dumps(list(map(str,command)))+'\n');f.flush()
        subprocess.run(list(map(str,command)),stdout=f,stderr=subprocess.STDOUT,cwd=str(cwd),check=True,timeout=timeout)

def predict(row):
    tip=row['tip'];source=SOURCE/'inputs/genomes'/(tip+'.fna')
    out=ROOT/'proteomes'/(tip+'.faa');done=ROOT/'proteomes'/(tip+'.json')
    digest=hashlib.sha256(source.read_bytes()).hexdigest()
    if done.exists() and out.exists():
        record=json.loads(done.read_text())
        if record['input_sha256']==digest and record['output_sha256']==hashlib.sha256(out.read_bytes()).hexdigest():return record
        raise RuntimeError('Protein cache differs from frozen reference: '+tip)
    started=time.monotonic()
    run([PRODIGAL,'-i',source,'-a',str(out)+'.partial','-o',ROOT/'proteomes'/(tip+'.genes.gff'),'-f','gff','-p','single','-g','11','-q'],ROOT/'logs'/(tip+'.prodigal.log'),ROOT,timeout=900)
    records=list(SeqIO.parse(str(out)+'.partial','fasta'))
    assert records and len({r.id for r in records})==len(records)
    assert all(len(r.seq)>0 for r in records)
    Path(str(out)+'.partial').replace(out)
    record=dict(tip=tip,proteins=len(records),seconds=time.monotonic()-started,input_sha256=digest,output_sha256=hashlib.sha256(out.read_bytes()).hexdigest(),predictor='Prodigal2.6.3 single; genetic code11')
    done.write_text(json.dumps(record,indent=2))
    print('PROTEOME_COMPLETED',tip,len(records),flush=True)
    return record

def prepare_branch(name,rows):
    branch=ROOT/name
    for d in ['inputs/proteomes','logs','tree','figures','source_data','config']:(branch/d).mkdir(parents=True,exist_ok=True)
    for d,target in [('databases',ORIGINAL/'databases'),('software',ORIGINAL/'software')]:
        if not (branch/d).exists():(branch/d).symlink_to(target,target_is_directory=True)
    for row in rows:
        p=branch/'inputs/proteomes'/(row['tip']+'.faa')
        if not p.exists():p.symlink_to(ROOT/'proteomes'/p.name)
    assert {p.stem for p in (branch/'inputs/proteomes').glob('*.faa')}=={r['tip'] for r in rows}
    write_table(branch/'source_data/tree_species.tsv',rows)
    shutil.copy2(SOURCE/'source_data/denominator.json',branch/'source_data/denominator.json')
    shutil.copy2(ROOT/'config/proteins.cfg',branch/'config/proteins.cfg')
    return branch

def infer(branch,timeout=None):
    pkg=ENV/'lib/python3.8/site-packages/phylophlan'
    run([ENV/'bin/python',ROOT/'scripts/phylophlan_proteins.py','-i','proteomes','--input_folder','inputs/proteomes','-d','phylophlan','--databases_folder','databases','--submat_folder',pkg/'phylophlan_substitution_matrices','--submod_folder',pkg/'phylophlan_substitution_models','--diversity','high','--accurate','--remove_fragmentary_entries','--min_num_markers','100','-f','config/proteins.cfg','--nproc',str(THREADS),'--output_folder','tree','--verbose'],branch/'logs/phylophlan.log',branch,timeout)
    trees=list((branch/'tree').rglob('*.treefile'));assert len(trees)==1 and trees[0].stat().st_size
    shutil.copy2(trees[0],branch/'tree/species_phylogeny.nwk')
    run([R,ROOT/'scripts/plot_tree.R'],branch/'logs/plot.log',branch,timeout=300)
    assert all((branch/('figures/Figure1b_phylogeny.'+ext)).stat().st_size>1000 for ext in ['pdf','svg','png','tiff'])
    # Count markers represented in the extracted marker records for every reference.
    import bz2
    counts=[]
    for path in (branch/'tree/proteomes_phylophlan/tmp/markers_aa').glob('*.bz2'):
        with bz2.open(path,'rt') as f:count=sum(1 for r in SeqIO.parse(f,'fasta'))
        counts.append(dict(tip=path.name.replace('.faa.bz2',''),markers=count))
    assert counts and all(r['markers']>=100 for r in counts)
    write_table(branch/'source_data/marker_counts.tsv',counts)
    return counts

def main(mode):
    os.chdir(ROOT)
    all_rows=load(SOURCE/'source_data/tree_species.tsv');assert len(all_rows)==207
    if mode=='pilot':
        targets={'2962065','477697','312279','669359','118562','454171'}
        rows=[r for r in all_rows if r['taxid'] in targets];assert len(rows)==6
        state('real_genome_pilot_prediction')
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(predict,rows))
        branch=prepare_branch('pilot6',rows)
        state('real_genome_pilot_phylogeny')
        counts=infer(branch,timeout=1200)
        assert len(counts)==6
        passed=dict(status='passed',marker_counts=counts,time=time.strftime('%F %T'),input='six real reference genomes including the previous stalled genome',pipeline='Prodigal -> DIAMOND blastp ->400 markers -> MAFFT -> trimming/subsampling -> IQ-TREE1000UFBoot1000SHaLRT -> four plot formats')
        (ROOT/'PILOT_PASSED.json').write_text(json.dumps(passed,indent=2));state('pilot_passed')
    else:
        assert (ROOT/'PILOT_PASSED.json').exists()
        state('predicting_207_proteomes')
        with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:prediction=list(pool.map(predict,all_rows))
        write_table(ROOT/'source_data/proteome_prediction_audit.tsv',prediction)
        baseline=load(ORIGINAL/'source_data/tree_species.tsv');assert len(baseline)==199
        by_tip={r['tip']:r for r in all_rows}
        assert all(by_tip[r['tip']]['reads']==r['reads'] and by_tip[r['tip']]['assembly_accession']==r['assembly_accession'] for r in baseline)
        supp=prepare_branch('supplemented207',all_rows);base=prepare_branch('baseline199',baseline)
        # The six fully validated pilot proteomes are unchanged: reuse their completed mapping.
        for folder in ['map_aa','markers_aa']:
            src=ROOT/'pilot6/tree/proteomes_phylophlan/tmp'/folder
            dst=supp/'tree/proteomes_phylophlan/tmp'/folder;dst.mkdir(parents=True,exist_ok=True)
            for file in src.iterdir():
                if file.is_file() and file.stat().st_size and not (dst/file.name).exists():
                    (dst/file.name).symlink_to(file)
        state('supplemented207_phylogeny');infer(supp)
        # Reuse only finished per-proteome marker mapping; redo all joint alignments for the199 set.
        common={r['tip'] for r in baseline}
        for folder in ['map_aa','markers_aa']:
            src=supp/'tree/proteomes_phylophlan/tmp'/folder
            dst=base/'tree/proteomes_phylophlan/tmp'/folder;dst.mkdir(parents=True,exist_ok=True)
            for file in src.iterdir():
                tip=file.name.split('.')[0]
                if tip in common and file.is_file() and file.stat().st_size and not (dst/file.name).exists():
                    (dst/file.name).symlink_to(file)
        state('baseline199_phylogeny');infer(base)
        run([R,ROOT/'scripts/compare_trees.R'],ROOT/'logs/comparison.log',ROOT,timeout=300)
        state('both_complete_pending_visual_review')

if __name__=='__main__':
    mode=argparse.ArgumentParser();mode.add_argument('mode',choices=['pilot','full']);args=mode.parse_args()
    try:main(args.mode)
    except Exception as e:state('failed',error=str(e));raise
