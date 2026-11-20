#!/usr/bin/env python3
import gzip
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import numpy as np
import pandas as pd
from Bio import SeqIO
from upsetplot import UpSet, from_indicators


BASE = Path("/public3/labmember/qianqh/soil_amp_analysis")
OUT = BASE / "AMP_consensus_out"
FIG = OUT / "figures"
TAB = OUT / "tables"
FAA = OUT / "fasta"
for directory in (OUT, FIG, TAB, FAA):
    directory.mkdir(parents=True, exist_ok=True)

INPUT_FASTA = BASE / "6_smORF_filter/blastp/smORF_prefiltered.faa"
MACREL_FILE = BASE / "macrel_out/macrel.out.prediction.gz"
ESM_FILE = BASE / "8_esm2lgbm10/smorf_pred_lgbm_esm2_mean_max.csv"
amplify_files = sorted((BASE / "AMPlify_out").glob("AMPlify_imbalanced_results_*.tsv"))
if not amplify_files:
    raise FileNotFoundError("AMPlify result TSV not found")
AMPLIFY_FILE = amplify_files[-1]


def parse_bool(value):
    return str(value).strip().lower() in {"true", "1", "yes", "amp"}


records = list(SeqIO.parse(str(INPUT_FASTA), "fasta"))
master = pd.DataFrame({
    "sequence_id": [r.id for r in records],
    "sequence": [str(r.seq).upper() for r in records],
})
master["length"] = master["sequence"].str.len()
if master["sequence_id"].duplicated().any() or master["sequence"].duplicated().any():
    raise ValueError("Input FASTA IDs and sequences must be unique")

macrel = pd.read_csv(MACREL_FILE, sep="\t", comment="#", compression="gzip")
macrel = macrel.rename(columns={
    "Access": "sequence_id", "Sequence": "macrel_sequence", "is_AMP": "macrel_pred",
    "AMP_probability": "macrel_probability", "AMP_family": "macrel_family",
    "Hemolytic": "macrel_hemolytic", "Hemolytic_probability": "macrel_hemolytic_probability",
})
macrel["macrel_pred"] = macrel["macrel_pred"].map(parse_bool)

amplify = pd.read_csv(AMPLIFY_FILE, sep="\t").rename(columns={
    "Sequence_ID": "sequence_id", "Sequence": "sequence",
    "Probability_score": "amplify_probability", "Prediction": "amplify_prediction",
    "Charge": "amplify_charge", "AMPlify_log_scaled_score": "amplify_log_score",
})
amplify["amplify_pred"] = amplify["amplify_prediction"].eq("AMP")

esm = pd.read_csv(ESM_FILE).rename(columns={
    "Sequence": "sequence", "prob_mean_max": "esm2_probability", "pred_mean_max": "esm2_pred"
})
esm["sequence"] = esm["sequence"].str.upper()
esm["esm2_pred"] = esm["esm2_pred"].astype(int).astype(bool)

df = master.merge(
    macrel[["sequence_id", "macrel_sequence", "macrel_pred", "macrel_probability", "macrel_family",
            "macrel_hemolytic", "macrel_hemolytic_probability"]],
    on="sequence_id", how="left", validate="one_to_one"
).merge(
    amplify[["sequence_id", "sequence", "amplify_pred", "amplify_probability",
             "amplify_log_score", "amplify_charge"]],
    on=["sequence_id", "sequence"], how="left", validate="one_to_one"
).merge(
    esm[["sequence", "esm2_pred", "esm2_probability"]],
    on="sequence", how="left", validate="one_to_one"
)

pred_cols = ["esm2_pred", "macrel_pred", "amplify_pred"]
if len(df) != len(master) or df[pred_cols].isna().any().any():
    missing = df[pred_cols].isna().sum().to_dict()
    raise RuntimeError(f"At least one tool did not return exactly one prediction per input sequence: {missing}")
macrel_sequence_ok = (df["sequence"] == df["macrel_sequence"]) | (df["sequence"].str[1:] == df["macrel_sequence"])
if not macrel_sequence_ok.all():
    raise RuntimeError(f"Macrel ID/sequence mismatch for {(~macrel_sequence_ok).sum()} records")
for col in pred_cols:
    df[col] = df[col].astype(bool)
df["support_count"] = df[pred_cols].sum(axis=1).astype(int)
df["consensus_pattern"] = np.select([
    df.esm2_pred & df.macrel_pred & df.amplify_pred,
    df.esm2_pred & df.macrel_pred & ~df.amplify_pred,
    df.esm2_pred & ~df.macrel_pred & df.amplify_pred,
    ~df.esm2_pred & df.macrel_pred & df.amplify_pred,
    df.esm2_pred & ~df.macrel_pred & ~df.amplify_pred,
    ~df.esm2_pred & df.macrel_pred & ~df.amplify_pred,
    ~df.esm2_pred & ~df.macrel_pred & df.amplify_pred,
], [
    "All three", "ESM2 + Macrel", "ESM2 + AMPlify", "Macrel + AMPlify",
    "ESM2 only", "Macrel only", "AMPlify only",
], default="None")

ordered = ["All three", "ESM2 + Macrel", "ESM2 + AMPlify", "Macrel + AMPlify",
           "ESM2 only", "Macrel only", "AMPlify only", "None"]
counts = df["consensus_pattern"].value_counts().reindex(ordered, fill_value=0)
counts.rename_axis("consensus_pattern").reset_index(name="n_sequences").to_csv(
    TAB / "intersection_counts.tsv", sep="\t", index=False)

df.to_csv(TAB / "amp_predictions_integrated.tsv.gz", sep="\t", index=False, compression="gzip")
df[df.support_count > 0].to_csv(TAB / "predicted_by_any_tool.tsv", sep="\t", index=False)
df[df.support_count >= 2].to_csv(TAB / "consensus_at_least_two.tsv", sep="\t", index=False)
df[df.support_count == 3].to_csv(TAB / "consensus_all_three.tsv", sep="\t", index=False)


def write_fasta(subset, path):
    with open(path, "w") as handle:
        for row in subset.itertuples(index=False):
            handle.write(f">{row.sequence_id}\n{row.sequence}\n")


write_fasta(df[df.support_count > 0], FAA / "predicted_by_any_tool.faa")
write_fasta(df[df.support_count >= 2], FAA / "consensus_at_least_two.faa")
write_fasta(df[df.support_count == 3], FAA / "consensus_all_three.faa")
for pattern in ordered[:-1]:
    safe = pattern.lower().replace(" + ", "_").replace(" ", "_")
    write_fasta(df[df.consensus_pattern == pattern], FAA / f"{safe}.faa")

names = ["ESM2-LGBM", "Macrel", "AMPlify"]
cols = ["esm2_pred", "macrel_pred", "amplify_pred"]
pair_rows = []
jaccard = np.eye(3)
for i, (name_a, col_a) in enumerate(zip(names, cols)):
    for j, (name_b, col_b) in enumerate(zip(names, cols)):
        a, b = df[col_a], df[col_b]
        inter = int((a & b).sum())
        union = int((a | b).sum())
        jac = inter / union if union else 1.0
        jaccard[i, j] = jac
        if i < j:
            pair_rows.append({"tool_a": name_a, "tool_b": name_b, "intersection": inter,
                              "union": union, "jaccard": jac,
                              "overlap_coefficient": inter / min(int(a.sum()), int(b.sum()))})
pd.DataFrame(pair_rows).to_csv(TAB / "pairwise_overlap_metrics.tsv", sep="\t", index=False)

summary = {
    "input_sequences": int(len(df)),
    "esm2_lgbm_positive": int(df.esm2_pred.sum()),
    "macrel_positive": int(df.macrel_pred.sum()),
    "amplify_positive": int(df.amplify_pred.sum()),
    "union_any_tool": int((df.support_count > 0).sum()),
    "at_least_two_tools": int((df.support_count >= 2).sum()),
    "intersection_all_three": int((df.support_count == 3).sum()),
    "source_files": {"input": str(INPUT_FASTA), "macrel": str(MACREL_FILE),
                     "esm2_lgbm": str(ESM_FILE), "amplify": str(AMPLIFY_FILE)},
}
with open(OUT / "summary.json", "w") as handle:
    json.dump(summary, handle, indent=2)

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                     "axes.spines.top": False, "axes.spines.right": False})
colors = ["#4477AA", "#228833", "#CC6677"]

# Figure 1: Venn diagram
fig, ax = plt.subplots(figsize=(7.2, 6.4))
v = venn3(subsets=[set(df.index[df[c]]) for c in cols], set_labels=names, ax=ax,
          set_colors=colors, alpha=0.65)
ax.set_title("Agreement among AMP prediction tools")
fig.tight_layout()
for ext in ("png", "pdf", "svg"):
    fig.savefig(FIG / f"figure1_venn_three_predictors.{ext}", dpi=600, bbox_inches="tight")
plt.close(fig)

# Figure 2: UpSet plot
indicator = df[cols].rename(columns=dict(zip(cols, names)))
upset_data = from_indicators(names, indicator)
fig = plt.figure(figsize=(9.0, 5.8))
up = UpSet(upset_data, subset_size="count", show_counts=True, sort_by="cardinality",
           facecolor="#4477AA")
up.plot(fig=fig)
fig.suptitle("Intersection structure of predicted AMPs", y=1.01)
for ext in ("png", "pdf", "svg"):
    fig.savefig(FIG / f"figure2_upset_intersections.{ext}", dpi=600, bbox_inches="tight")
plt.close(fig)

# Figure 3: positive counts
positive_counts = [int(df[c].sum()) for c in cols]
fig, ax = plt.subplots(figsize=(6.8, 4.8))
bars = ax.bar(names, positive_counts, color=colors, width=0.62)
ax.set_ylabel("Predicted AMP sequences")
ax.set_title("Positive calls by prediction tool")
ax.bar_label(bars, labels=[f"{n:,}\n({n/len(df):.2%})" for n in positive_counts], padding=4)
ax.set_ylim(0, max(positive_counts) * 1.22)
fig.tight_layout()
for ext in ("png", "pdf", "svg"):
    fig.savefig(FIG / f"figure3_positive_counts.{ext}", dpi=600, bbox_inches="tight")
plt.close(fig)

# Figure 4: Jaccard heatmap
fig, ax = plt.subplots(figsize=(5.6, 4.8))
im = ax.imshow(jaccard, cmap="Blues", vmin=0, vmax=1)
ax.set_xticks(range(3), names, rotation=25, ha="right")
ax.set_yticks(range(3), names)
for i in range(3):
    for j in range(3):
        ax.text(j, i, f"{jaccard[i, j]:.3f}", ha="center", va="center",
                color="white" if jaccard[i, j] > 0.55 else "black")
fig.colorbar(im, ax=ax, label="Jaccard index", fraction=0.046, pad=0.04)
ax.set_title("Pairwise agreement of positive sets")
fig.tight_layout()
for ext in ("png", "pdf", "svg"):
    fig.savefig(FIG / f"figure4_jaccard_heatmap.{ext}", dpi=600, bbox_inches="tight")
plt.close(fig)

# Figure 5: exact consensus patterns excluding common negatives
plot_counts = counts.drop("None")
fig, ax = plt.subplots(figsize=(7.6, 4.8))
bars = ax.barh(plot_counts.index[::-1], plot_counts.values[::-1], color="#4477AA")
ax.set_xlabel("Number of sequences")
ax.set_title("Exact positive-call combinations")
ax.bar_label(bars, padding=3, fmt="%d")
ax.set_xlim(0, max(plot_counts.max() * 1.18, 1))
fig.tight_layout()
for ext in ("png", "pdf", "svg"):
    fig.savefig(FIG / f"figure5_exact_combinations.{ext}", dpi=600, bbox_inches="tight")
plt.close(fig)

print(json.dumps(summary, indent=2))
