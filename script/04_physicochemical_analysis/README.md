# Peptide physicochemical properties

These scripts compare sequence-derived properties of candidate peptides with reference AMPs.

## Analysis steps

1. Calculate per-peptide descriptors, including length, isoelectric point, charge and hydrophobicity, and summarize group differences. [physchem_four_model_compare.py](physchem_four_model_compare.py)

2. Plot the selected descriptors from the feature table. [physchem_make_focused_figure.py](physchem_make_focused_figure.py)

3. Prepare the comparison panels and assemble Figure 2. [prepare_panels.R](prepare_panels.R), [build_Fig2_revision.R](build_Fig2_revision.R)
