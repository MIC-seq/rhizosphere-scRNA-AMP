#!/usr/bin/env python3
"""Convert featureCounts CORE assignments into collision-safe gene x cell UMI matrices."""
import argparse,csv,gzip,re,sys
from collections import Counter,defaultdict

QNAME_RE=re.compile(r'_([ACGTN]{20})_([ACGTN]{8})$')

def load_whitelist(path):
    with open(path,encoding="utf-8") as fh: return {line.split("\t",1)[0].strip() for line in fh if line.strip()}

def extract(args):
    whitelist=load_whitelist(args.whitelist); qc=Counter()
    with open(args.core,encoding="utf-8",errors="replace") as fh, open(args.output,"w",encoding="utf-8",newline="\n") as out:
        for line in fh:
            qc["core_records"]+=1; x=line.rstrip("\n").split("\t")
            if len(x)<4 or x[1]!="Assigned": qc["not_assigned"]+=1; continue
            m=QNAME_RE.search(x[0])
            if not m: qc["bad_read_id"]+=1; continue
            cell,umi=m.groups()
            if cell not in whitelist: qc["not_whitelisted"]+=1; continue
            genes=[g for g in x[3].split(",") if g and g!="NA"]
            if len(genes)!=1: qc["ambiguous_read_gene"]+=1; continue
            out.write("%s\t%s\t%s\n"%(cell,umi,genes[0])); qc["assigned_triples"]+=1
    with open(args.qc,"w",encoding="utf-8",newline="") as fh:
        w=csv.writer(fh,delimiter="\t",lineterminator="\n"); w.writerow(["metric","value"]); w.writerows(sorted(qc.items()))

def resolve(args):
    qc=Counter(); last=None; genes=[]
    with open(args.input,encoding="utf-8") as fh, open(args.output,"w",encoding="utf-8",newline="\n") as out:
        for line in fh:
            cell,umi,gene=line.rstrip("\n").split("\t"); key=(cell,umi)
            if last is not None and key!=last:
                qc["molecules_seen"]+=1
                if len(genes)==1: out.write("%s\t%s\t1\n"%(genes[0],last[0])); qc["molecules_assigned"]+=1
                else: qc["molecules_ambiguous_gene"]+=1
                genes=[]
            if not genes or gene!=genes[-1]: genes.append(gene)
            last=key
        if last is not None:
            qc["molecules_seen"]+=1
            if len(genes)==1: out.write("%s\t%s\t1\n"%(genes[0],last[0])); qc["molecules_assigned"]+=1
            else: qc["molecules_ambiguous_gene"]+=1
    with open(args.qc,"a",encoding="utf-8",newline="") as fh:
        w=csv.writer(fh,delimiter="\t",lineterminator="\n"); w.writerows(sorted(qc.items()))

def annotation(path):
    genes=[]; meta={}
    for r in csv.DictReader(open(path,encoding="utf-8"),delimiter="\t"):
        gid=r["gene_id"]; genes.append(gid)
        text=(r.get("gene_biotype","")+" "+r.get("product","")).lower()
        excluded=("rrna" in text or "ribosomal rna" in text or "trna" in text or "transfer rna" in text)
        meta[gid]=(excluded,r["species"],r["assembly_accession"])
    return genes,meta

def build(args):
    genes,meta=annotation(args.annotation); gene_index={g:i+1 for i,g in enumerate(genes)}
    cells=[line.split("\t",1)[0].strip() for line in open(args.whitelist,encoding="utf-8") if line.strip()]
    cell_index={c:i+1 for i,c in enumerate(cells)}; raw_nnz=non_nnz=0
    with open(args.counts,encoding="utf-8") as fh:
        for line in fh:
            g,c,n=line.rstrip("\n").split("\t"); raw_nnz+=1
            if g in meta and not meta[g][0]: non_nnz+=1
    def headers(path,nnz):
        f=gzip.open(path,"wt",encoding="utf-8"); f.write("%%%%MatrixMarket matrix coordinate integer general\n%%\n%d %d %d\n"%(len(genes),len(cells),nnz)); return f
    raw=headers(args.raw_matrix,raw_nnz); non=headers(args.non_rrna_matrix,non_nnz)
    cell_total=Counter(); cell_non=Counter(); cell_genes=Counter(); cell_non_genes=Counter(); species=Counter(); missing=Counter()
    with open(args.counts,encoding="utf-8") as fh:
        for line in fh:
            g,c,n=line.rstrip("\n").split("\t"); n=int(n)
            if g not in gene_index or c not in cell_index: missing["unknown_gene_or_cell"]+=1; continue
            i=gene_index[g]; j=cell_index[c]; raw.write("%d %d %d\n"%(i,j,n)); cell_total[c]+=n; cell_genes[c]+=1
            excluded,sp,acc=meta[g]; species[(c,sp,acc)]+=n
            if not excluded: non.write("%d %d %d\n"%(i,j,n)); cell_non[c]+=n; cell_non_genes[c]+=1
    raw.close(); non.close()
    with gzip.open(args.features,"wt",encoding="utf-8") as f:
        for g in genes: f.write(g+"\n")
    with gzip.open(args.barcodes,"wt",encoding="utf-8") as f:
        for c in cells: f.write(args.sample+"__"+c+"\n")
    with open(args.cell_qc,"w",encoding="utf-8",newline="") as fh:
        w=csv.writer(fh,delimiter="\t",lineterminator="\n"); w.writerow(["sample","cell_barcode","cell_id","total_gene_umi","detected_genes","total_non_rrna_umi","detected_non_rrna_genes","rrna_trna_umi","rrna_trna_fraction"])
        for c in cells:
            rr=cell_total[c]-cell_non[c]; frac=rr/cell_total[c] if cell_total[c] else 0
            w.writerow([args.sample,c,args.sample+"__"+c,cell_total[c],cell_genes[c],cell_non[c],cell_non_genes[c],rr,frac])
    with gzip.open(args.species_umi,"wt",encoding="utf-8",newline="") as fh:
        w=csv.writer(fh,delimiter="\t",lineterminator="\n"); w.writerow(["sample","cell_barcode","cell_id","species","assembly_accession","umi_count"])
        for (c,sp,acc),n in sorted(species.items()): w.writerow([args.sample,c,args.sample+"__"+c,sp,acc,n])

def main():
    p=argparse.ArgumentParser(); sub=p.add_subparsers(dest="cmd")
    x=sub.add_parser("extract"); x.add_argument("--core",required=True); x.add_argument("--whitelist",required=True); x.add_argument("--output",required=True); x.add_argument("--qc",required=True)
    x=sub.add_parser("resolve"); x.add_argument("--input",required=True); x.add_argument("--output",required=True); x.add_argument("--qc",required=True)
    x=sub.add_parser("build"); x.add_argument("--counts",required=True); x.add_argument("--annotation",required=True); x.add_argument("--whitelist",required=True); x.add_argument("--sample",required=True); x.add_argument("--raw-matrix",required=True); x.add_argument("--non-rrna-matrix",required=True); x.add_argument("--features",required=True); x.add_argument("--barcodes",required=True); x.add_argument("--cell-qc",required=True); x.add_argument("--species-umi",required=True)
    a=p.parse_args(); {"extract":extract,"resolve":resolve,"build":build}[a.cmd](a)
if __name__=="__main__": main()
