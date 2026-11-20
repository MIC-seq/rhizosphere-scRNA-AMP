#!/usr/bin/env python3
"""Assess mapping-fraction thresholds and sample/species pseudobulk eligibility."""
import argparse, csv, gzip, json
from collections import Counter


def quantile(x, q):
    if not x: return ""
    x = sorted(x); pos = (len(x) - 1) * q; lo = int(pos); hi = min(lo + 1, len(x) - 1); f = pos - lo
    return x[lo] * (1 - f) + x[hi] * f


def main():
    p=argparse.ArgumentParser(); p.add_argument("--config",required=True); p.add_argument("--input",required=True); p.add_argument("--out-prefix",required=True); a=p.parse_args()
    cfg=json.load(open(a.config,encoding="utf-8")); cq=cfg["cell_qc"]; aa=cfg["amp_analysis"]
    with gzip.open(a.input,"rt",encoding="utf-8") as fh: rows=list(csv.DictReader(fh,delimiter="\t"))
    base=[]
    for r in rows:
        r["_base"]=(int(r["total_non_rrna_umi"])>=cq["min_non_rrna_umi"] and int(r["detected_non_rrna_genes"])>=cq["min_detected_non_rrna_genes"] and int(r["amp_valid"])==1)
        if r["_base"]: base.append(r)
    samples=sorted({r["sample"] for r in rows}); thresholds=[0,0.001,0.002,0.005,0.01,0.02,0.05]
    with open(a.out_prefix+".thresholds.tsv","w",encoding="utf-8",newline="") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n"); w.writerow(["sample","base_cells","p10","p25","median","p75","p90"]+[f"cells_ge_{t:g}" for t in thresholds])
        for s in samples:
            sr=[r for r in base if r["sample"]==s]; x=[float(r["unique_mapping_fraction"]) for r in sr]
            w.writerow([s,len(sr)]+[quantile(x,q) for q in (.1,.25,.5,.75,.9)]+[sum(v>=t for v in x) for t in thresholds])
    selected=cq["min_unique_mapping_fraction"]; counts=Counter()
    for r in base:
        if float(r["unique_mapping_fraction"])>=selected: counts[(r["species"],r["sample"])]+=1
    species=sorted({k[0] for k in counts})
    min_group=aa["min_cells_per_state_per_sample"]; need_total=5*min_group
    with open(a.out_prefix+".species_sample.tsv","w",encoding="utf-8",newline="") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n"); w.writerow(["species","sample","qc_cells","estimated_tail_cells","sample_eligible"])
        for sp in species:
            for s in samples:
                n=counts[(sp,s)]
                if n: w.writerow([sp,s,n,n//5,int(n>=need_total)])
    with open(a.out_prefix+".species_summary.tsv","w",encoding="utf-8",newline="") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n"); w.writerow(["species","qc_cells","samples_with_cells","samples_ge_required_total","pseudobulk_eligible"])
        for sp in species:
            ns=[counts[(sp,s)] for s in samples]; elig=sum(n>=need_total for n in ns)
            w.writerow([sp,sum(ns),sum(n>0 for n in ns),elig,int(elig>=aa["min_samples"])])
    with open(a.out_prefix+".threshold_species_eligibility.tsv","w",encoding="utf-8",newline="") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n"); w.writerow(["min_unique_mapping_fraction","retained_cells","species_with_cells","pseudobulk_eligible_species","eligible_species"])
        for threshold in thresholds:
            tc=Counter()
            for r in base:
                if float(r["unique_mapping_fraction"])>=threshold: tc[(r["species"],r["sample"])]+=1
            spp=sorted({k[0] for k in tc}); eligible=[]
            for sp in spp:
                if sum(tc[(sp,s)]>=need_total for s in samples)>=aa["min_samples"]: eligible.append(sp)
            w.writerow([threshold,sum(tc.values()),len(spp),len(eligible),";".join(eligible)])
    print(f"base_cells={len(base)} selected_mapping_fraction={selected}")
if __name__=="__main__": main()
