#!/usr/bin/env python3
"""Single-pass split of corrected R2 reads by Kraken barcode for directed validation."""
import argparse,csv,gzip,json,re
from collections import Counter
from pathlib import Path
QNAME_RE=re.compile(r'_([ACGTN]{20})_([ACGTN]{8})$')
def main():
    p=argparse.ArgumentParser(); p.add_argument("--config",required=True); p.add_argument("--sample",required=True); p.add_argument("--fastq",required=True); p.add_argument("--taxonomy",required=True); p.add_argument("--targets",required=True); p.add_argument("--outdir",required=True); a=p.parse_args()
    cfg=json.load(open(a.config,encoding="utf-8")); level=cfg["directed_mapping"]["gzip_compresslevel"]; out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    targets={r["taxid"]:r for r in csv.DictReader(open(a.targets,encoding="utf-8"),delimiter="\t")}; bc={}
    with open(a.taxonomy,encoding="utf-8") as fh:
        for r in csv.DictReader(fh,delimiter="\t"):
            if r["sample_id"]==a.sample and r["taxid"] in targets: bc[r["cell_barcode"]]=targets[r["taxid"]]["assembly_accession"]
    handles={}; counts=Counter(); cells=Counter()
    try:
        with gzip.open(a.fastq,"rt",encoding="utf-8",errors="replace") as fh:
            while True:
                rec=[fh.readline() for _ in range(4)]
                if not rec[0]: break
                if any(x=="" for x in rec): raise SystemExit("truncated FASTQ")
                counts["input_reads"]+=1; m=QNAME_RE.search(rec[0].rstrip())
                if not m: counts["bad_read_id"]+=1; continue
                barcode=m.group(1); acc=bc.get(barcode)
                if not acc: counts["non_target_reads"]+=1; continue
                if acc not in handles: handles[acc]=gzip.open(str(out/(acc+".fq.gz")),"wt",encoding="utf-8",compresslevel=level)
                handles[acc].writelines(rec); counts["target_reads"]+=1; cells[(acc,barcode)]+=1
    finally:
        for h in handles.values(): h.close()
    with open(out/"split_qc.tsv","w",encoding="utf-8",newline="") as fh:
        w=csv.writer(fh,delimiter="\t",lineterminator="\n"); w.writerow(["metric","value"]); w.writerows(sorted(counts.items()))
    with open(out/"target_barcodes.tsv","w",encoding="utf-8",newline="") as fh:
        w=csv.writer(fh,delimiter="\t",lineterminator="\n"); w.writerow(["assembly_accession","cell_barcode","reads"]); w.writerows((acc,cb,n) for (acc,cb),n in sorted(cells.items()))
if __name__=="__main__": main()
