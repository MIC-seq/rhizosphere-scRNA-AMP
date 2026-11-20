#!/usr/bin/env python3
"""Translate CDS features from prefixed directed-reference GTF and genomic FASTA."""
import argparse, csv, re
from collections import defaultdict
from pathlib import Path

ATTR=re.compile(r'(\S+) "([^"]*)"')
COMP=str.maketrans("ACGTNacgtn","TGCANtgcan")
CODE={
"TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S","TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"*","TGG":"W",
"CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
"ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T","AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
"GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A","GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G"}

def fasta(path):
    d={}; name=None; seq=[]
    for line in open(path,encoding="utf-8"):
        if line.startswith(">"):
            if name: d[name]="".join(seq).upper()
            name=line[1:].split()[0]; seq=[]
        else: seq.append(line.strip())
    if name: d[name]="".join(seq).upper()
    return d

def wrap(seq,n=80): return "\n".join(seq[i:i+n] for i in range(0,len(seq),n))

def main():
    p=argparse.ArgumentParser(); p.add_argument("--targets",required=True); p.add_argument("--output",required=True); p.add_argument("--summary",required=True); a=p.parse_args()
    targets=list(csv.DictReader(open(a.targets,encoding="utf-8"),delimiter="\t")); counts=[]; seen=set(); Path(a.output).parent.mkdir(parents=True,exist_ok=True)
    with open(a.output,"w",encoding="utf-8") as out:
        for t in targets:
            root=Path(t["target_dir"])/"reference"; genome=fasta(root/"combined_genome.fa"); cds=defaultdict(list); info={}
            for line in open(root/"combined_annotation.gtf",encoding="utf-8"):
                x=line.rstrip("\n").split("\t")
                if len(x)!=9 or x[2]!="CDS": continue
                at=dict(ATTR.findall(x[8])); gene=at.get("gene_id")
                if not gene: continue
                cds[gene].append((x[0],int(x[3]),int(x[4]),x[6],int(x[7]) if x[7].isdigit() else 0)); info[gene]=at
            written=0
            for gene,parts in cds.items():
                strand=parts[0][3]; parts=sorted(parts,key=lambda z:z[1],reverse=strand=="-"); nts=[]
                for contig,start,end,st,phase in parts:
                    s=genome[contig][start-1:end]
                    if st=="-": s=s.translate(COMP)[::-1]
                    nts.append(s[phase:])
                nt="".join(nts); aa="".join(CODE.get(nt[i:i+3],"X") for i in range(0,len(nt)-2,3)).rstrip("*")
                if not aa: continue
                if gene in seen: raise SystemExit(f"duplicate protein gene {gene}")
                seen.add(gene); at=info[gene]
                out.write(f">{gene} protein_id={at.get('protein_id','')} species={t['species']} assembly={t['assembly_accession']}\n{wrap(aa)}\n"); written+=1
            counts.append({"species":t["species"],"assembly_accession":t["assembly_accession"],"cds_genes":len(cds),"proteins_written":written})
    with open(a.summary,"w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(counts[0]),delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(counts)
    print("proteins",len(seen),"targets",len(targets))
if __name__=="__main__": main()
