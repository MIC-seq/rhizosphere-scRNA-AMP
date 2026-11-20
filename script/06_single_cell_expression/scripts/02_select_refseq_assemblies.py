#!/usr/bin/env python3
"""Rank local RefSeq assemblies and retain explicit inclusion/exclusion evidence."""
import argparse,csv,json,re
from datetime import datetime
from pathlib import Path

LEVEL={"Complete Genome":0,"Chromosome":1,"Scaffold":2,"Contig":3}
CAT={"reference genome":0,"representative genome":1,"na":2,"":2}
def clean(x): return x.strip().lstrip("# ")
def main():
    p=argparse.ArgumentParser(); p.add_argument("--config",required=True); p.add_argument("--species",required=True); p.add_argument("--assembly-summary",required=True); p.add_argument("--refseq-root",required=True); p.add_argument("--outdir",required=True); a=p.parse_args()
    cfg=json.load(open(a.config,encoding="utf-8"))["refseq"]; out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    selection_column=cfg.get("selection_column","selected")
    wanted={(r["taxid"],r["species"]) for r in csv.DictReader(open(a.species,encoding="utf-8"),delimiter="\t") if r[selection_column] in ("1","True","true")}
    wanted_by_taxid={x[0]:x for x in wanted}; audit=[]; header=None
    with open(a.assembly_summary,encoding="utf-8") as fh:
      for line in fh:
        if line.startswith("# assembly_accession") or line.startswith("#assembly_accession"):
            header=[clean(x) for x in line.rstrip().split("\t")]; continue
        if line.startswith("#") or not line.strip(): continue
        if header is None: raise SystemExit("assembly_summary header not found")
        r=dict(zip(header,line.rstrip("\n").split("\t")))
        taxid=r.get("species_taxid",r.get("taxid","")); name=r.get("organism_name","")
        if taxid not in wanted_by_taxid: continue
        matches=[wanted_by_taxid[taxid]]
        wtax,wname=matches[0]; reasons=[]; level=r.get("assembly_level",""); cat=r.get("refseq_category","").lower()
        if cfg["require_exact_species_taxid"] and taxid!=wtax: reasons.append("taxid_mismatch")
        if level not in cfg["allowed_assembly_levels"]: reasons.append("assembly_level_not_allowed")
        txt=" ".join(r.values()).lower()
        if cfg["exclude_anomalous"] and "anomal" in txt: reasons.append("anomalous")
        if cfg["exclude_contaminated"] and "contamin" in txt: reasons.append("contaminated")
        acc=r["assembly_accession"]; ftp=r.get("ftp_path",""); base=ftp.rstrip("/").split("/")[-1]
        local=Path(a.refseq_root)
        fna=local/(base+"_genomic.fna.gz")
        gff=local/(base+cfg.get("annotation_extension","_genomic.gff.gz"))
        faa=local/(base+"_protein.faa.gz")
        if not fna.exists(): reasons.append("missing_fasta")
        if cfg["require_gff"] and not gff.exists(): reasons.append("missing_gff")
        if cfg["require_faa"] and not faa.exists(): reasons.append("missing_faa")
        date=r.get("seq_rel_date",""); date_score=-(int(re.sub(r"\D","",date) or 0))
        audit.append({"species":wname,"taxid":wtax,"organism_name":name,"assembly_accession":acc,"assembly_level":level,"refseq_category":r.get("refseq_category",""),"genome_size":r.get("genome_size",""),"contig_count":r.get("contig_count",""),"seq_rel_date":date,"genome_rep":r.get("genome_rep",""),"fasta_path":str(fna),"gff_path":str(gff),"faa_path":str(faa) if faa.exists() else "","protein_strategy":cfg.get("protein_strategy","existing_FAA"),"eligible":int(not reasons),"exclusion_reasons":"PASS" if not reasons else ";".join(reasons),"_score":(CAT.get(cat,3),LEVEL.get(level,9),date_score,acc)})
    chosen=set()
    for sp in {x["species"] for x in audit}:
        eligible=sorted((x for x in audit if x["species"]==sp and x["eligible"]),key=lambda x:x["_score"])
        chosen.update(x["assembly_accession"] for x in eligible[:cfg["assemblies_per_species_initial"]])
    for x in audit: x["selected"]=int(x["assembly_accession"] in chosen); x.pop("_score")
    fields=list(audit[0]) if audit else ["species","taxid","assembly_accession","eligible","selected","exclusion_reasons"]
    with open(out/"refseq_assembly_audit.tsv","w",newline="",encoding="utf-8") as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(audit)
    with open(out/"selected_refseq_assemblies.tsv","w",newline="",encoding="utf-8") as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(x for x in audit if x["selected"])
if __name__=="__main__": main()
