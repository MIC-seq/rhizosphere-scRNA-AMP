#!/usr/bin/env python3
"""Merge RefSeq gene metadata, translated proteins, and eggNOG COG annotations."""
import argparse
import csv
from pathlib import Path


def read_fasta(path):
    seqs, name, parts = {}, None, []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(parts)
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip())
    if name is not None:
        seqs[name] = "".join(parts)
    return seqs


def read_eggnog(path):
    rows, header = {}, None
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##") or not line.strip():
                continue
            if line.startswith("#"):
                header = line[1:].rstrip("\n").split("\t")
                continue
            if header is None:
                raise SystemExit("eggNOG header line was not found")
            values = line.rstrip("\n").split("\t")
            row = dict(zip(header, values))
            query = row.get("query") or row.get("#query")
            if query:
                rows[query] = row
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene-annotation", required=True)
    parser.add_argument("--proteins", required=True)
    parser.add_argument("--eggnog", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    proteins = read_fasta(args.proteins)
    eggnog = read_eggnog(args.eggnog)
    with open(args.gene_annotation, encoding="utf-8") as handle:
        genes = list(csv.DictReader(handle, delimiter="\t"))

    core = list(genes[0]) if genes else ["gene_id"]
    extra = ["protein_sequence", "COG_category", "COG_description",
             "eggNOG_ortholog", "eggNOG_evalue", "eggNOG_score"]
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    annotated = 0
    with open(args.output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=core + extra,
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for gene in genes:
            gene_id = gene["gene_id"]
            hit = eggnog.get(gene_id, {})
            category = hit.get("COG_category", "-")
            if category not in ("", "-"):
                annotated += 1
            gene.update({
                "protein_sequence": proteins.get(gene_id, ""),
                "COG_category": category,
                "COG_description": hit.get("Description", "-"),
                "eggNOG_ortholog": hit.get("seed_ortholog", "-"),
                "eggNOG_evalue": hit.get("evalue", "-"),
                "eggNOG_score": hit.get("score", "-"),
            })
            writer.writerow(gene)

    summary = [
        ("total_genes", len(genes)),
        ("proteins_present", sum(g["gene_id"] in proteins for g in genes)),
        ("eggnog_hits", sum(g["gene_id"] in eggnog for g in genes)),
        ("genes_with_cog_category", annotated),
    ]
    with open(args.summary, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerows(summary)
    print(dict(summary))


if __name__ == "__main__":
    main()
