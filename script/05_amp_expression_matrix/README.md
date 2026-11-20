# Single-cell AMP expression matrix

The workflow maps reads to a nonredundant peptide reference and produces cell-by-feature UMI counts.

## Analysis steps

1. Build the reference from candidate peptides and the APD6 and AMPSphere collections. [01_build_reference.pbs](01_build_reference.pbs)

2. Align each sample, assign reads to features and collapse UMIs within cell barcodes. [02_align_sample.pbs](02_align_sample.pbs), [03_count_sample.pbs](03_count_sample.pbs), [count_amp_umis.py](count_amp_umis.py)

3. Schedule the per-sample alignment and counting jobs after the reference job. [submit_pipeline.sh](submit_pipeline.sh)

4. Merge the sample matrices, filter recurrent features and cluster cells by AMP expression. [analyze_and_plot.R](analyze_and_plot.R)
