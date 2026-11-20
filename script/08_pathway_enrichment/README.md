# Pathway and metabolic enrichment

Gene counts, state labels and KEGG pathway membership are used to compare metabolic activity across cell states.

## Analysis steps

1. Prepare gene sets, align state labels with cells and build the analysis objects. [01_prepare_inputs.py](01_prepare_inputs.py), [02_analyse_metabolism.R](02_analyse_metabolism.R)

2. Compute pathway ssGSEA scores and metabolic enrichment scores (MES), then combine results across species and states. [run_author.R](run_author.R), [aggregate_author.R](aggregate_author.R)

3. Score six broad metabolic classes. [score_six_classes.R](score_six_classes.R)

4. Summarize paired contrasts across soils and states and adjust for multiple testing. [analyse_fg.R](analyse_fg.R), [statistics_v6.R](statistics_v6.R)

5. Render the main and supplementary pathway panels from the prepared score and contrast tables. [main_supp_v8.R](main_supp_v8.R)
