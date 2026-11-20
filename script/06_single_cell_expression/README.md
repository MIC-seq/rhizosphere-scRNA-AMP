# Species-resolved gene expression and functional states

Species-specific gene matrices and AMP counts are linked by cell barcode. COG profiles provide the functional-state labels used downstream.

## Analysis steps

1. Select taxa with adequate cell coverage and identify their reference assemblies. [01_select_major_species.py](scripts/01_select_major_species.py), [02_select_refseq_assemblies.py](scripts/02_select_refseq_assemblies.py)

2. Build species references, split reads by taxonomic assignment and count mapped molecules. [06_build_directed_references.pbs](pbs/06_build_directed_references.pbs), [07_split_sample.pbs](pbs/07_split_sample.pbs), [08_directed_map_count_sample.pbs](pbs/08_directed_map_count_sample.pbs)

3. Integrate expression and taxonomy and assess cell eligibility. [09_integrate_directed_mapping_amp.py](scripts/09_integrate_directed_mapping_amp.py), [10_assess_directed_qc_eligibility.py](scripts/10_assess_directed_qc_eligibility.py)

4. Annotate reference proteins and identify cells shared by the gene and AMP matrices. [13_eggnog_cog_annotation.pbs](pbs/13_eggnog_cog_annotation.pbs), [22_build_amp_gene_cell_bridge.py](scripts/22_build_amp_gene_cell_bridge.py)

5. Build functional matrices, cluster cells by COG profile and export state labels for association analysis. [23_prepare_species_function_matrices.py](scripts/23_prepare_species_function_matrices.py), [24_species_function_umap.R](scripts/24_species_function_umap.R), [28_prepare_optimized_figure5_data.R](scripts/28_prepare_optimized_figure5_data.R)
