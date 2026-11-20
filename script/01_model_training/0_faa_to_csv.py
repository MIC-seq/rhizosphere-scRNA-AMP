from Bio import SeqIO
import pandas as pd

fasta_path="dataset/dataset.faa"
csv_out="dataset/dataset.csv"

def fasta_to_csv(fasta_path, csv_out):
    seqs, lens = [], []
    for r in SeqIO.parse(fasta_path, "fasta"):
        seqs.append(str(r.seq))
        lens.append(len(r.seq))
    pd.DataFrame({"Sequence": seqs, "Length": lens, "label": [None]*len(seqs)}).to_csv(csv_out, index=False)

fasta_to_csv(fasta_path, csv_out)

