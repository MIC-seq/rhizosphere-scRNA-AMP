#!/usr/bin/env python3
"""Normalize per-barcode taxonomy and select major species with audited reasons."""
import argparse, csv, gzip, json
from collections import Counter, defaultdict
from pathlib import Path

def op(path, mode="rt"):
    return gzip.open(path, mode, encoding="utf-8") if str(path).endswith(".gz") else open(path, mode, encoding="utf-8")

def rows(path):
    with op(path) as fh:
        first = fh.readline(); fh.seek(0)
        if "\t" in first and any(x in first.lower() for x in ("barcode", "species", "taxid")):
            yield from csv.DictReader(fh, delimiter="\t")
        else:
            for x in csv.reader(fh, delimiter="\t"):
                if x: yield {"cell_barcode": x[0], "label": x[1] if len(x)>1 else ""}

def first(d, *names, default=""):
    for n in names:
        if n in d and str(d[n]).strip(): return str(d[n]).strip()
    return default

def read_amp_cells(path):
    if not path or not Path(path).exists(): return set()
    ans=set()
    for r in rows(path):
        cb=first(r,"cell_barcode","barcode","Cell")
        if cb: ans.add(cb)
    return ans

def read_barcode_taxonomy(path):
    """MICtools barcode_count: barcode,total,classified,genus_taxid,genus_support,
    genus_name,species_taxid,species_support,species_name."""
    out={}
    with op(path) as fh:
        for x in csv.reader(fh,delimiter="\t"):
            if len(x)<9: continue
            cb=x[0].strip(); species=x[8].strip()
            if species.startswith("s_"): species=species[2:]
            genus=x[5].strip()
            if genus.startswith("g_"): genus=genus[2:]
            taxid=x[6].strip()
            identifiable=taxid not in ("","0") and species.lower() not in ("unindentifiable","unidentifiable")
            out[cb]={"taxid":taxid,"rank":"species" if identifiable else "unclassified","genus":genus,"species":species if identifiable else "","support_reads":x[1].strip(),"classification_support":x[7].strip()}
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--config",required=True); ap.add_argument("--samples",required=True); ap.add_argument("--outdir",required=True)
    a=ap.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    cfg=json.load(open(a.config,encoding="utf-8")); th=cfg["species_filter"]; mr=cfg["mapping_reference_filter"]
    manifest=list(csv.DictReader(open(a.samples,encoding="utf-8"),delimiter="\t"))
    cells=[]
    for s in manifest:
        authoritative=read_barcode_taxonomy(s["barcode_support"])
        amp=read_amp_cells(s.get("amp_matrix",""))
        for cb,t in authoritative.items():
            cells.append({"cell_id":s["sample_id"]+"__"+cb,"sample_id":s["sample_id"],"crop":s["crop"],"batch":s["batch"],"cell_barcode":cb,**t,"amp_valid":int(cb in amp)})
    keys=["cell_id","sample_id","crop","batch","cell_barcode","taxid","rank","genus","species","support_reads","classification_support","amp_valid"]
    with open(out/"cell_taxonomy.tsv","w",newline="",encoding="utf-8") as fh:
        w=csv.DictWriter(fh,fieldnames=keys,delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(cells)
    by=defaultdict(list)
    for c in cells:
        if c["species"]: by[(c["taxid"],c["species"])].append(c)
    result=[]; sample_species_totals=Counter(x["sample_id"] for x in cells if x["rank"]=="species"); all_species_cells=sum(sample_species_totals.values())
    for (taxid,species), z in sorted(by.items(),key=lambda kv:(-len(kv[1]),kv[0][1])):
        per=Counter(x["sample_id"] for x in z); crops=sorted(set(x["crop"] for x in z)); ranks=Counter(x["rank"].lower() for x in z)
        checks={"species_rank":ranks["species"]==len(z),"not_excluded_taxid":taxid not in {str(x) for x in th.get("excluded_taxids",[])},"sample_coverage":len(per)>=th["min_samples"],"total_cells":len(z)>=th["min_total_cells"],"samples_with_min_cells":sum(n>=th["min_cells_per_sample"] for n in per.values())>=th["min_samples_at_cell_threshold"]}
        failed=[k for k,v in checks.items() if not v]
        vals=[]; read_support=[]
        for x in z:
            try: vals.append(float(x["classification_support"])); read_support.append(float(x["support_reads"]))
            except ValueError: pass
        global_fraction=len(z)/all_species_cells
        max_sample_fraction=max(per[s]/sample_species_totals[s] for s in per)
        prevalence=len(per)/len(manifest)
        map_checks={"global_fraction":global_fraction>=mr["min_global_cell_fraction"],"max_sample_fraction":max_sample_fraction>=mr["min_max_sample_cell_fraction"],"sample_prevalence":prevalence>=mr["min_sample_prevalence_fraction"]}
        map_reasons=[k for k,v in map_checks.items() if v]
        map_selected=taxid not in {str(x) for x in mr.get("exclude_taxids",[])} and any(map_checks.values())
        result.append({"taxid":taxid,"species":species,"total_cells":len(z),"global_cell_fraction":global_fraction,"n_samples":len(per),"sample_prevalence_fraction":prevalence,"max_sample_cell_fraction":max_sample_fraction,"n_crops":len(crops),"crops":";".join(crops),"samples_ge_threshold":sum(n>=th["min_cells_per_sample"] for n in per.values()),"amp_valid_cells":sum(x["amp_valid"] for x in z),"classification_support_mean":sum(vals)/len(vals) if vals else "","support_reads_mean":sum(read_support)/len(read_support) if read_support else "","sample_cell_counts":";".join(f"{k}:{v}" for k,v in sorted(per.items())),"mapping_selected":int(map_selected),"mapping_selection_reasons":";".join(map_reasons) if map_selected else "BELOW_ALL_UNION_THRESHOLDS_OR_EXCLUDED","analysis_selected":int(not failed),"analysis_failure_reasons":"PASS" if not failed else ";".join(failed),"selected":int(not failed),"failure_reasons":"PASS" if not failed else ";".join(failed)})
    if sum(x["mapping_selected"] for x in result)>mr["max_species"]:
        raise SystemExit("mapping reference exceeds max_species; increase configured prevalence/abundance thresholds")
    fields=list(result[0]) if result else ["taxid","species","selected","failure_reasons"]
    with open(out/"species_selection.tsv","w",newline="",encoding="utf-8") as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(result)

if __name__=="__main__": main()
