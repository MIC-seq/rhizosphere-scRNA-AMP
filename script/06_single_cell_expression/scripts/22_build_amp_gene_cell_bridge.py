#!/usr/bin/env python3
"""Join the established Figure 4 AMP UMAP to directed gene-expression metadata."""
import argparse
import csv
import gzip
from collections import Counter
from pathlib import Path


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else open(path, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--amp-umap", required=True)
    parser.add_argument("--gene-cells", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    with open_text(args.gene_cells) as handle:
        gene_rows = {row["cell_id"]: row for row in csv.DictReader(handle, delimiter="\t")}

    gene_fields = [
        "assembly_accession", "kraken_taxid", "genus", "species",
        "target_input_reads", "unique_mapped_reads", "assigned_gene_reads",
        "unique_mapping_fraction", "total_gene_umi", "detected_genes",
        "total_non_rrna_umi", "detected_non_rrna_genes", "rrna_trna_fraction",
        "amp_load", "amp_valid", "qc_pass", "qc_failure_reasons",
    ]
    output_fields = [
        "cell_id", "sample", "crop", "cell_barcode", "AMP_UMAP_1", "AMP_UMAP_2",
        "AMP_cluster", "AMP_UMI_Fig4", "stable_AMP_count_Fig4",
        "AMP_positive_Fig4", "AMP_species_Fig4", "AMP_genus_Fig4",
        "gene_expression_available",
    ] + ["gene_" + field for field in gene_fields]

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    counts = Counter()
    targets = {"Chloracidobacterium validum", "Halomonas sp. JS92-SW72", "Pseudomonas aeruginosa"}
    with open_text(args.amp_umap) as source, gzip.open(args.output, "wt", encoding="utf-8", newline="") as destination:
        reader = csv.DictReader(source, delimiter="\t")
        writer = csv.DictWriter(destination, fieldnames=output_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for amp in reader:
            gene = gene_rows.get(amp["cell_id"])
            row = {
                "cell_id": amp["cell_id"], "sample": amp["sample"], "crop": amp["crop"],
                "cell_barcode": amp["cell_barcode"], "AMP_UMAP_1": amp["UMAP_1"],
                "AMP_UMAP_2": amp["UMAP_2"], "AMP_cluster": amp["cluster"],
                "AMP_UMI_Fig4": amp["amp_umi"], "stable_AMP_count_Fig4": amp["n_amp"],
                "AMP_positive_Fig4": amp["amp_positive"], "AMP_species_Fig4": amp["species"],
                "AMP_genus_Fig4": amp["genus"], "gene_expression_available": int(gene is not None),
            }
            for field in gene_fields:
                row["gene_" + field] = gene.get(field, "") if gene else ""
            writer.writerow(row)
            counts["amp_umap_cells"] += 1
            if gene:
                counts["matched_gene_cells"] += 1
                if gene.get("qc_pass") == "1":
                    counts["matched_gene_qc_cells"] += 1
                    if gene.get("species") in targets:
                        counts["matched_target_species_qc_cells"] += 1

    counts["directed_gene_cells"] = len(gene_rows)
    counts["directed_gene_qc_cells"] = sum(row.get("qc_pass") == "1" for row in gene_rows.values())
    with open(args.summary, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        for key in sorted(counts):
            writer.writerow([key, counts[key]])
    print(dict(counts))


if __name__ == "__main__":
    main()
