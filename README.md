# rhizosphere-scRNA-AMP

Analysis code for the study of antimicrobial peptide (AMP) expression in individual bacterial cells from soybean, rice and wheat rhizospheres. The study combines peptide discovery across 49 datasets with single-cell expression analysis to examine how candidate AMPs vary across taxa and bacterial functional states.

## Data

The Figure 1b inputs include the [207-tip reference tree](script/03_phylogenetic_tree/supplemented207/tree/species_phylogeny.nwk) and [taxon abundance table](script/03_phylogenetic_tree/supplemented207/source_data/tree_species.tsv). The [example sample manifest](config/species_expression/samples.example.tsv) lists the eight study samples. Copy it to `config/samples.tsv` in your species-expression analysis directory and replace the example paths with local input paths. Set the read, database, software and output paths for your installation before running the cluster jobs.

## Analysis scripts

| Module | Analysis |
|---|---|
| [AMP classifier training](script/01_model_training/README.md) | The training scripts take labeled peptide sequences and produce a fine-tuned ESM2 encoder and a ten-fold LightGBM ensemble. |
| [Candidate AMP discovery](script/02_amp_qc_prediction/README.md) | This module processes rhizosphere transcript assemblies and combines predictions from ESM2–LightGBM, Macrel, AMPlify and AMP Scanner v2. |
| [Reference phylogeny and Figure 1b](script/03_phylogenetic_tree/README.md) | The tree relates reference taxa identified across 49 rhizosphere datasets. The abundance bars summarize pooled transcript reads. |
| [Peptide physicochemical properties](script/04_physicochemical_analysis/README.md) | These scripts compare sequence-derived properties of candidate peptides with reference AMPs. |
| [Single-cell AMP expression matrix](script/05_amp_expression_matrix/README.md) | The workflow maps reads to a nonredundant peptide reference and produces cell-by-feature UMI counts. |
| [Species-resolved gene expression and functional states](script/06_single_cell_expression/README.md) | Species-specific gene matrices and AMP counts are linked by cell barcode. COG profiles provide the functional-state labels used downstream. |
| [AMP expression and COG functional states](script/07_cog_amp_centroid_association/README.md) | Three analyses use matched cells to relate AMP expression to independently defined functional states within species. |
| [Pathway and metabolic enrichment](script/08_pathway_enrichment/README.md) | Gene counts, state labels and KEGG pathway membership are used to compare metabolic activity across cell states. |

## Data visualization

To redraw Figure 1b from the repository root, run

```bash
Rscript script/03_phylogenetic_tree/scripts/01_render_figure1b.R
```

The renderer writes PDF, SVG, PNG and TIFF files under `script/03_phylogenetic_tree/supplemented207/figures/`.

### Requirements

- **Tree rendering** uses R with Cairo support and the ape, svglite and ragg packages.
- **Expression and pathway analyses** use Seurat, Matrix, data.table and GSVA. The species expression and association parameters are in [config/species_expression](config/species_expression/).
- **Peptide modeling and characterization** use Python with PyTorch, fair-esm, LightGBM, scikit-learn, Biopython and modlAMP. AMP Scanner v2 uses a separate TensorFlow/Keras environment compatible with its trained model.
- **Read processing and annotation** use Linux command-line tools including STAR, rnaSPAdes, TransDecoder, MMseqs2, DIAMOND and eggNOG-mapper. The supplied job scripts use PBS.
