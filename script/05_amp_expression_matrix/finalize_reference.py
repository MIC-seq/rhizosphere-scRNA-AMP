#!/usr/bin/env python3
"""Create unified representative IDs and source membership after MMseqs2 clustering."""
from argparse import ArgumentParser
from collections import defaultdict
from pathlib import Path
import csv

PRIORITY = {"APD6": 0, "AMPSphere": 1, "Novel": 2}

def fasta(path):
    h, seq = None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line: continue
            if line.startswith(">"):
                if h is not None: yield h.split()[0], "".join(seq)
                h, seq = line[1:], []
            else: seq.append(line)
    if h is not None: yield h.split()[0], "".join(seq)

def read_tsv(path):
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def main():
    p = ArgumentParser()
    p.add_argument("--records", required=True)
    p.add_argument("--clusters", required=True, help="MMseqs createtsv representative<TAB>member")
    p.add_argument("--tagged-fasta", required=True)
    p.add_argument("--outdir", required=True)
    a = p.parse_args()
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    recs = {r["tagged_id"]: r for r in read_tsv(a.records)}
    seqs = dict(fasta(a.tagged_fasta))
    members = defaultdict(list)
    with open(a.clusters, encoding="utf-8") as fh:
        for line in fh:
            rep, member = line.rstrip("\n").split("\t")[:2]
            members[rep].append(member)
    reps = sorted(members)
    with open(out/"amp_reference_nr.faa", "w", encoding="utf-8") as fo, \
         open(out/"amp_reference_metadata.tsv", "w", newline="", encoding="utf-8") as mo, \
         open(out/"amp_reference_membership.tsv", "w", newline="", encoding="utf-8") as mm:
        mw = csv.writer(mo, delimiter="\t", lineterminator="\n")
        xw = csv.writer(mm, delimiter="\t", lineterminator="\n")
        mw.writerow(["amp_id","representative_tag","length","primary_source","source_pattern","n_members","n_APD6","n_AMPSphere","n_Novel","eligible_10_100"])
        xw.writerow(["amp_id","representative_tag","member_tag","source","source_id","original_header","member_length"])
        for i, rep in enumerate(reps, 1):
            amp_id = f"AMPCLUST_{i:07d}"
            cluster = members[rep]
            sources = sorted({recs[m]["source"] for m in cluster}, key=PRIORITY.get)
            primary = sources[0]
            counts = {s: sum(recs[m]["source"] == s for m in cluster) for s in PRIORITY}
            seq = seqs[rep]
            fo.write(f">{amp_id} primary_source={primary} source_pattern={'+'.join(sources)} representative={rep}\n{seq}\n")
            mw.writerow([amp_id,rep,len(seq),primary,"+".join(sources),len(cluster),counts["APD6"],counts["AMPSphere"],counts["Novel"],int(10 <= len(seq) <= 100)])
            for m in cluster:
                r = recs[m]
                xw.writerow([amp_id,rep,m,r["source"],r["source_id"],r["original_header"],r["length"]])

if __name__ == "__main__":
    main()
