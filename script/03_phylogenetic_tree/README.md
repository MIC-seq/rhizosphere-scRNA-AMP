# Reference phylogeny and Figure 1b

The tree relates reference taxa identified across 49 rhizosphere datasets. The abundance bars summarize pooled transcript reads.

## Analysis steps

1. Summarize Kraken2 assignments and select reference genomes. [workflow.py](scripts/workflow.py)

2. Infer the protein-marker phylogeny from the selected references. [protein_route.py](scripts/protein_route.py)

3. Render the annotated 207-tip tree from the supplied Newick tree and species table. [01_render_figure1b.R](scripts/01_render_figure1b.R), [plot_tree_reviewed.R](scripts/plot_tree_reviewed.R)

## Render Figure 1b

From this directory, run

```bash
Rscript scripts/01_render_figure1b.R
```

The renderer writes the annotated tree to `supplemented207/figures/`. It requires ape, svglite, ragg and an R build with Cairo support.
