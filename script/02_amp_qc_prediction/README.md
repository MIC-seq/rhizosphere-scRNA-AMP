# Candidate AMP discovery

This module processes rhizosphere transcript assemblies and combines predictions from ESM2–LightGBM, Macrel, AMPlify and AMP Scanner v2.

## Analysis steps

1. Remove rRNA reads, assemble transcripts and extract open reading frames. [1_star_rrna.pbs](1_star_rrna.pbs), [2_assembly.pbs](2_assembly.pbs), [3_transdecoder_orf.pbs](3_transdecoder_orf.pbs)

2. Remove redundant sequences and filter short ORFs using reference annotations. The Macrel call in the first job belongs to the early screening stage. [4_nr_macrel.pbs](4_nr_macrel.pbs), [5_smORF_filter.pbs](5_smORF_filter.pbs)

3. Score peptide sequences with AMP Scanner v2 using its TensorFlow/Keras model. [run_ampscanner_v2.pbs](run_ampscanner_v2.pbs), [amp_scanner_v2_predict_tf1.py](amp_scanner_v2_predict_tf1.py)

4. Merge the three predictor outputs, then add AMP Scanner votes to obtain the candidate sets. [amp_consensus_analysis.py](amp_consensus_analysis.py), [amp_four_model_analysis.py](amp_four_model_analysis.py)

5. Assess upstream in-frame stop codons in transcript context. [upstream_stop_transcriptome.py](upstream_stop_transcriptome.py)
