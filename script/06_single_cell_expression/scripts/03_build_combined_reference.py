#!/usr/bin/env python3
"""Stream selected RefSeq FNA/GTF files into one collision-safe reference."""
import argparse,csv,gzip,json,re
from collections import Counter
from pathlib import Path

ATTR_KEYS=("gene_id","transcript_id","locus_tag","old_locus_tag","protein_id")
ATTR_RE=re.compile(r'(?P<key>gene_id|transcript_id|locus_tag|old_locus_tag|protein_id) "(?P<val>[^"]*)"')

def op(path):
    return gzip.open(path,"rt",encoding="utf-8",errors="replace") if str(path).endswith(".gz") else open(path,encoding="utf-8",errors="replace")

def prefix_attrs(text,prefix):
    def f(m):
        val=m.group("val")
        return '%s "%s"'%(m.group("key"),prefix+"__"+val) if val else m.group(0)
    return ATTR_RE.sub(f,text)

def attrs(text):
    return {m.group(1):m.group(2) for m in re.finditer(r'([A-Za-z_][A-Za-z0-9_]*) "([^"]*)"',text)}

def main():
    p=argparse.ArgumentParser(); p.add_argument("--assemblies",required=True); p.add_argument("--outdir",required=True); a=p.parse_args()
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    assemblies=[x for x in csv.DictReader(open(a.assemblies,encoding="utf-8"),delimiter="\t") if x["selected"]=="1"]
    feature_counts=Counter(); bases=contigs=0; gene_rows={}
    with open(out/"combined_genome.fa","w",encoding="utf-8",newline="\n") as fa, open(out/"combined_annotation.gtf","w",encoding="utf-8",newline="\n") as gt, open(out/"contig_map.tsv","w",encoding="utf-8",newline="") as cm:
        mw=csv.writer(cm,delimiter="\t",lineterminator="\n"); mw.writerow(["assembly_accession","species","original_contig","combined_contig"])
        for r in assemblies:
            acc=r["assembly_accession"]; prefix=acc
            with op(r["fasta_path"]) as fh:
                for line in fh:
                    if line.startswith(">"):
                        orig=line[1:].split()[0]; new=prefix+"__"+orig; fa.write(">"+new+"\n"); mw.writerow([acc,r["species"],orig,new]); contigs+=1
                    else:
                        seq=line.strip(); fa.write(seq+"\n"); bases+=len(seq)
            with op(r["gff_path"]) as fh:
                for line in fh:
                    if line.startswith("#") or not line.strip(): continue
                    x=line.rstrip("\n").split("\t")
                    if len(x)!=9: raise SystemExit("Malformed GTF in %s"%r["gff_path"])
                    x[0]=prefix+"__"+x[0]; x[8]=prefix_attrs(x[8],prefix); gt.write("\t".join(x)+"\n"); feature_counts[x[2]]+=1
                    if x[2] in ("gene","CDS"):
                        d=attrs(x[8]); gid=d.get("gene_id","")
                        if gid:
                            z=gene_rows.setdefault(gid,{"gene_id":gid,"locus_tag":d.get("locus_tag",""),"protein_id":d.get("protein_id",""),"species":r["species"],"taxid":r["taxid"],"assembly_accession":acc,"contig":x[0],"start":x[3],"end":x[4],"strand":x[6],"gene_biotype":d.get("gene_biotype",""),"product":d.get("product","")})
                            for k in ("locus_tag","protein_id","gene_biotype","product"):
                                if not z.get(k) and d.get(k): z[k]=d[k]
    fields=["gene_id","locus_tag","protein_id","species","taxid","assembly_accession","contig","start","end","strand","gene_biotype","product"]
    with open(out/"gene_annotation_refseq.tsv","w",encoding="utf-8",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(gene_rows[k] for k in sorted(gene_rows))
    summary={"assemblies":len(assemblies),"species":len(set(x["taxid"] for x in assemblies)),"bases":bases,"contigs":contigs,"genes":len(gene_rows),"features":dict(feature_counts)}
    with open(out/"reference_summary.json","w",encoding="utf-8") as fh: json.dump(summary,fh,indent=2,sort_keys=True)
    print(json.dumps(summary,sort_keys=True))
if __name__=="__main__": main()
