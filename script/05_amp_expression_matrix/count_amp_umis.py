#!/usr/bin/env python3
"""Filter DIAMOND hits and assign each Cell-UMI molecule to one unambiguous AMP."""
from argparse import ArgumentParser
from collections import Counter, defaultdict
from pathlib import Path
import csv, gzip, re, sqlite3

BC_RE = re.compile(r"^[ACGTN]{20}$")
UMI_RE = re.compile(r"^[ACGTN]{8}$")

def op(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)

def load_cells(path):
    cells = set()
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            x = line.split("\t", 1)[0].strip()
            if x: cells.add(x)
    return cells

def parse_molecule(qid):
    bits = qid.split("_")
    if len(bits) < 3: return None
    cb, umi = bits[-2], bits[-1].split()[0]
    return (cb, umi) if BC_RE.match(cb) and UMI_RE.match(umi) else None

def main():
    p = ArgumentParser()
    p.add_argument("--hits", required=True)
    p.add_argument("--whitelist", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--mode", choices=["primary","sensitivity"], default="primary")
    p.add_argument("--min-score-gap", type=float, default=5.0)
    a = p.parse_args()
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    cells = load_cells(a.whitelist)
    if a.mode == "primary":
        min_ident, min_cov, min_aln, max_eval, min_bits = 90.0, 0.80, 10, 1e-5, 30.0
    else:
        min_ident, min_cov, min_aln, max_eval, min_bits = 70.0, 0.70, 8, 1e-5, 25.0
    db = sqlite3.connect(str(out/f"{a.sample}.{a.mode}.sqlite"))
    db.execute("PRAGMA journal_mode=WAL")
    db.execute("PRAGMA synchronous=OFF")
    db.execute("CREATE TABLE IF NOT EXISTS mol(cell TEXT, umi TEXT, best_amp TEXT, best REAL, second REAL DEFAULT -1, PRIMARY KEY(cell,umi))")
    qc = Counter(); batch=[]
    with op(a.hits) as fh:
        for line in fh:
            qc["hit_lines"] += 1
            if qc["hit_lines"] % 1000000 == 0:
                print("[%s:%s] processed %d hit lines; %d passed" % (a.sample,a.mode,qc["hit_lines"],qc["passing_hits"]), flush=True)
            x = line.rstrip("\n").split("\t")
            if len(x) < 14: qc["malformed"] += 1; continue
            mol = parse_molecule(x[0])
            if mol is None: qc["bad_read_id"] += 1; continue
            cb, umi = mol
            if cb not in cells: qc["not_whitelisted"] += 1; continue
            ident, aln, evalue, bits, slen = float(x[2]), int(x[3]), float(x[10]), float(x[11]), int(x[13])
            scov = aln / slen if slen else 0
            if not (ident >= min_ident and scov >= min_cov and aln >= min_aln and evalue <= max_eval and bits >= min_bits):
                qc["failed_filter"] += 1; continue
            batch.append((cb,umi,x[1],bits))
            qc["passing_hits"] += 1
            if len(batch) >= 100000:
                update(db,batch); batch=[]
        if batch: update(db,batch)
    counts=defaultdict(int); assigned=ambiguous=0
    for cell,umi,amp,best,second in db.execute("SELECT cell,umi,best_amp,best,second FROM mol"):
        if second >= 0 and best-second < a.min_score_gap:
            ambiguous += 1; continue
        counts[(cell,amp)] += 1; assigned += 1
    with gzip.open(out/f"{a.sample}.{a.mode}.cell_amp_umi.tsv.gz","wt",newline="") as fh:
        w=csv.writer(fh,delimiter="\t",lineterminator="\n"); w.writerow(["sample","cell_barcode","amp_id","umi_count"])
        for (cell,amp),n in sorted(counts.items()): w.writerow([a.sample,cell,amp,n])
    qc.update({"whitelist_cells":len(cells),"assigned_molecules":assigned,"ambiguous_molecules":ambiguous,"expressed_cell_amp_pairs":len(counts)})
    with open(out/f"{a.sample}.{a.mode}.qc.tsv","w",newline="") as fh:
        w=csv.writer(fh,delimiter="\t",lineterminator="\n"); w.writerow(["metric","value"]); w.writerows(sorted(qc.items()))
    db.close()

def update(db, rows):
    cur=db.cursor()
    for cell,umi,amp,score in rows:
        old=cur.execute("SELECT best_amp,best,second FROM mol WHERE cell=? AND umi=?",(cell,umi)).fetchone()
        if old is None:
            cur.execute("INSERT INTO mol(cell,umi,best_amp,best,second) VALUES(?,?,?,?,?)",(cell,umi,amp,score,-1.0)); continue
        best_amp,best,second=old
        if amp == best_amp:
            if score > best: cur.execute("UPDATE mol SET best=? WHERE cell=? AND umi=?",(score,cell,umi))
        elif score > best:
            cur.execute("UPDATE mol SET best_amp=?,best=?,second=? WHERE cell=? AND umi=?",(amp,score,max(best,second),cell,umi))
        elif score > second:
            cur.execute("UPDATE mol SET second=? WHERE cell=? AND umi=?",(score,cell,umi))
    db.commit()

if __name__ == "__main__": main()
