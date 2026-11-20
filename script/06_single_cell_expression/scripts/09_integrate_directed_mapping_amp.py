#!/usr/bin/env python3
"""Integrate directed mapping QC, gene UMI summaries, Kraken metadata and AMP counts."""
import argparse
import csv
import gzip
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

QNAME_RE = re.compile(r'_([ACGTN]{20})_([ACGTN]{8})$')


def table(path):
    return csv.DictReader(open(path, encoding="utf-8"), delimiter="\t")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--config", required=True)
    p.add_argument("--samples", required=True)
    p.add_argument("--taxonomy", required=True)
    p.add_argument("--split-root", required=True)
    p.add_argument("--mapping-root", required=True)
    p.add_argument("--outdir", required=True)
    args = p.parse_args()
    cfg = json.load(open(args.config, encoding="utf-8"))
    scale = cfg["amp_analysis"]["scale_factor"]
    stable_min = cfg["amp_analysis"]["stable_amp_min_umi"]
    qcfg = cfg["cell_qc"]
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    sample_rows = list(table(args.samples))
    tax = {}
    for r in table(args.taxonomy):
        tax[(r["sample_id"], r["cell_barcode"])] = r

    records = []
    summary = []
    for sm in sample_rows:
        sample = sm["sample_id"]
        split = Path(args.split_root) / sample
        mapped = Path(args.mapping_root) / sample
        input_reads = Counter()
        accession = {}
        for r in table(split / "target_barcodes.tsv"):
            input_reads[r["cell_barcode"]] += int(r["reads"])
            accession[r["cell_barcode"]] = r["assembly_accession"]

        unique_reads = Counter()
        assigned_reads = Counter()
        for core in mapped.glob("tmp/*/unique.bam.featureCounts.gz"):
            with gzip.open(core, "rt", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    x = line.rstrip("\n").split("\t")
                    m = QNAME_RE.search(x[0]) if x else None
                    if not m:
                        continue
                    bc = m.group(1)
                    unique_reads[bc] += 1
                    if len(x) >= 4 and x[1] == "Assigned" and len([g for g in x[3].split(",") if g and g != "NA"]) == 1:
                        assigned_reads[bc] += 1

        gene_qc = {r["cell_barcode"]: r for r in table(mapped / f"{sample}_cell_qc.tsv")}
        amp_total = Counter()
        stable = Counter()
        amp_detected = Counter()
        with gzip.open(sm["amp_matrix"], "rt", encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                n = int(r["umi_count"])
                amp_total[r["cell_barcode"]] += n
                amp_detected[r["cell_barcode"]] += 1
                if n >= stable_min:
                    stable[r["cell_barcode"]] += 1

        counts = Counter()
        for bc, g in gene_qc.items():
            t = tax.get((sample, bc), {})
            inp = input_reads[bc]
            uniq = unique_reads[bc]
            assigned = assigned_reads[bc]
            non = int(g["total_non_rrna_umi"])
            genes = int(g["detected_non_rrna_genes"])
            uf = uniq / inp if inp else 0.0
            af = assigned / inp if inp else 0.0
            amp = amp_total[bc]
            amp_valid = int(t.get("amp_valid", "0") or 0)
            qc_umi = non >= qcfg["min_non_rrna_umi"]
            qc_gene = genes >= qcfg["min_detected_non_rrna_genes"]
            qc_map = uf >= qcfg["min_unique_mapping_fraction"]
            qc_pass = qc_umi and qc_gene and qc_map and amp_valid == 1
            reasons = []
            if not qc_umi: reasons.append("low_non_rrna_umi")
            if not qc_gene: reasons.append("low_detected_genes")
            if not qc_map: reasons.append("low_unique_mapping_fraction")
            if not amp_valid: reasons.append("no_valid_amp_profile")
            rec = {
                "sample": sample, "crop": sm["crop"], "batch": sm["batch"],
                "cell_barcode": bc, "cell_id": f"{sample}__{bc}",
                "kraken_taxid": t.get("taxid", ""), "genus": t.get("genus", ""), "species": t.get("species", ""),
                "assembly_accession": accession.get(bc, ""), "kraken_support_reads": t.get("support_reads", ""),
                "kraken_classification_support": t.get("classification_support", ""),
                "target_input_reads": inp, "unique_mapped_reads": uniq, "assigned_gene_reads": assigned,
                "unique_mapping_fraction": uf, "assigned_gene_read_fraction": af,
                "total_gene_umi": int(g["total_gene_umi"]), "detected_genes": int(g["detected_genes"]),
                "total_non_rrna_umi": non, "detected_non_rrna_genes": genes,
                "rrna_trna_umi": int(g["rrna_trna_umi"]), "rrna_trna_fraction": float(g["rrna_trna_fraction"]),
                "amp_umi": amp, "detected_amp_count": amp_detected[bc], "stable_amp_count": stable[bc],
                "amp_load": math.log1p(scale * amp / non) if non else "",
                "amp_positive": int(amp > 0), "amp_valid": amp_valid,
                "qc_pass": int(qc_pass), "qc_failure_reasons": "PASS" if qc_pass else ";".join(reasons),
            }
            records.append(rec)
            counts["directed_cells"] += 1
            counts["qc_pass"] += int(qc_pass)
            counts["amp_positive"] += int(amp > 0)
            counts["input_reads"] += inp
            counts["unique_reads"] += uniq
            counts["assigned_reads"] += assigned
            counts["non_rrna_umi"] += non
        summary.append({"sample": sample, **counts})

    fields = list(records[0])
    with gzip.open(out / "cell_directed_amp_integrated.tsv.gz", "wt", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(records)
    sfields = list(summary[0])
    with open(out / "sample_directed_qc_summary.tsv", "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sfields, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(summary)
    print(f"cells={len(records)} qc_pass={sum(int(r['qc_pass']) for r in records)}")


if __name__ == "__main__":
    main()
