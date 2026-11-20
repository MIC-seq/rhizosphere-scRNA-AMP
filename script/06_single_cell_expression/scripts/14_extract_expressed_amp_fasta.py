#!/usr/bin/env python3
"""Extract AMPCLUST representative sequences observed in the sample AMP matrices."""
import argparse, csv, gzip

def main():
    p=argparse.ArgumentParser(); p.add_argument("--samples",required=True); p.add_argument("--reference",required=True); p.add_argument("--output",required=True); p.add_argument("--ids",required=True); a=p.parse_args()
    wanted=set()
    for s in csv.DictReader(open(a.samples,encoding="utf-8"),delimiter="\t"):
        with gzip.open(s["amp_matrix"],"rt",encoding="utf-8") as f:
            for r in csv.DictReader(f,delimiter="\t"): wanted.add(r["amp_id"])
    found=set(); write=False
    with open(a.reference,encoding="utf-8") as src, open(a.output,"w",encoding="utf-8") as out:
        for line in src:
            if line.startswith(">"):
                amp=line[1:].split()[0]; write=amp in wanted
                if write: found.add(amp)
            if write: out.write(line)
    with open(a.ids,"w",encoding="utf-8") as f:
        f.write("amp_id\n"); f.writelines(x+"\n" for x in sorted(found))
    missing=wanted-found
    if missing: raise SystemExit(f"missing AMP sequences: {len(missing)}")
    print("expressed_amp",len(found))
if __name__=="__main__": main()
