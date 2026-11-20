#!/usr/bin/env python3
"""Create one-row assembly manifests for strict-analysis species-directed mapping."""
import argparse,csv,json,re
from pathlib import Path
def main():
    p=argparse.ArgumentParser(); p.add_argument("--config",required=True); p.add_argument("--species",required=True); p.add_argument("--assemblies",required=True); p.add_argument("--outdir",required=True); a=p.parse_args()
    cfg=json.load(open(a.config,encoding="utf-8")); col=cfg["directed_mapping"]["target_selection_column"]; excluded={str(x) for x in cfg["species_filter"].get("excluded_taxids",[])}; out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    wanted={r["taxid"]:r for r in csv.DictReader(open(a.species,encoding="utf-8"),delimiter="\t") if r[col]=="1" and r["taxid"] not in excluded}
    assemblies=[r for r in csv.DictReader(open(a.assemblies,encoding="utf-8"),delimiter="\t") if r["selected"]=="1" and r["taxid"] in wanted]
    fields=list(assemblies[0]) if assemblies else []
    target=[]
    for r in assemblies:
        d=out/r["assembly_accession"]; d.mkdir(parents=True,exist_ok=True); mf=d/"assembly.tsv"
        with open(mf,"w",encoding="utf-8",newline="") as fh:
            w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerow(r)
        target.append({"taxid":r["taxid"],"species":r["species"],"assembly_accession":r["assembly_accession"],"target_dir":str(d),"assembly_manifest":str(mf),"total_cells":wanted[r["taxid"]]["total_cells"],"sample_cell_counts":wanted[r["taxid"]]["sample_cell_counts"]})
    with open(out/"targets.tsv","w",encoding="utf-8",newline="") as fh:
        f=list(target[0]) if target else ["taxid","species","assembly_accession","target_dir","assembly_manifest","total_cells","sample_cell_counts"]; w=csv.DictWriter(fh,fieldnames=f,delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(target)
    missing=[r for t,r in wanted.items() if t not in {x["taxid"] for x in assemblies}]
    with open(out/"missing_targets.tsv","w",encoding="utf-8",newline="") as fh:
        f=list(missing[0]) if missing else ["taxid","species"]; w=csv.DictWriter(fh,fieldnames=f,delimiter="\t",lineterminator="\n",extrasaction="ignore"); w.writeheader(); w.writerows(missing)
    print("targets",len(target),"missing",len(missing))
if __name__=="__main__": main()
