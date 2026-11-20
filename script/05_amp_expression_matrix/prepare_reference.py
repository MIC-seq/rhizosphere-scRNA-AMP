#!/usr/bin/env python3
"""Merge AMP FASTAs with source-safe IDs before MMseqs2 clustering."""
from argparse import ArgumentParser
from collections import Counter
from pathlib import Path
import csv
import re

AA = set("ACDEFGHIKLMNPQRSTVWYXBZUOJ")

def records(path):
    header, seq = None, []
    with open(path, encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq).upper().rstrip("*")
                header, seq = line[1:], []
            else:
                seq.append(line)
    if header is not None:
        yield header, "".join(seq).upper().rstrip("*")

def safe_id(text):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text)

def main():
    p = ArgumentParser()
    p.add_argument("--apd6", required=True)
    p.add_argument("--ampsphere", required=True)
    p.add_argument("--novel", required=True)
    p.add_argument("--out-fasta", required=True)
    p.add_argument("--out-records", required=True)
    args = p.parse_args()
    inputs = [("APD6", args.apd6), ("AMPSphere", args.ampsphere), ("Novel", args.novel)]
    Path(args.out_fasta).parent.mkdir(parents=True, exist_ok=True)
    seen_tag = Counter()
    with open(args.out_fasta, "w", encoding="utf-8") as fasta_out, open(args.out_records, "w", newline="", encoding="utf-8") as tab_out:
        w = csv.writer(tab_out, delimiter="\t", lineterminator="\n")
        w.writerow(["tagged_id", "source", "source_id", "original_header", "length", "eligible_10_100"])
        for source, path in inputs:
            for header, seq in records(path):
                if not seq or any(c not in AA for c in seq):
                    raise ValueError(f"Invalid sequence in {path}: {header}")
                source_id = header.split()[0]
                base = f"{source}__{safe_id(source_id)}"
                seen_tag[base] += 1
                tag = base if seen_tag[base] == 1 else f"{base}__dup{seen_tag[base]}"
                fasta_out.write(f">{tag}\n{seq}\n")
                w.writerow([tag, source, source_id, header, len(seq), int(10 <= len(seq) <= 100)])

if __name__ == "__main__":
    main()
