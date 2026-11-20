#!/usr/bin/env python3
"""Combine gene annotations produced for all species-directed references."""
import argparse
import csv
from pathlib import Path


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--targets", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    rows = list(csv.DictReader(open(args.targets, encoding="utf-8"), delimiter="\t"))
    header = None
    genes = set()
    n = 0
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8", newline="") as dst:
        writer = None
        for target in rows:
            source = Path(target["target_dir"]) / "reference" / "gene_annotation_refseq.tsv"
            with open(source, encoding="utf-8") as src:
                reader = csv.DictReader(src, delimiter="\t")
                if header is None:
                    header = reader.fieldnames
                    writer = csv.DictWriter(dst, fieldnames=header, delimiter="\t", lineterminator="\n")
                    writer.writeheader()
                elif reader.fieldnames != header:
                    raise SystemExit(f"annotation columns differ: {source}")
                for row in reader:
                    gene = row["gene_id"]
                    if gene in genes:
                        raise SystemExit(f"duplicate gene_id: {gene}")
                    genes.add(gene)
                    writer.writerow(row)
                    n += 1
    print(f"targets={len(rows)} genes={n} output={out}")


if __name__ == "__main__":
    main()
