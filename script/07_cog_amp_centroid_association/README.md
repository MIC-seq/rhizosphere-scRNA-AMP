# AMP expression and COG functional states

Three analyses use matched cells to relate AMP expression to independently defined functional states within species.

## Analysis steps

1. Test associations between AMP expression clusters and COG states. [25_amp_function_association.R](scripts/25_amp_function_association.R)

2. Predict the state of cells in a held-out sample from AMP profiles using a nearest-centroid classifier. [30_amp_expression_predicts_functional_state.R](scripts/30_amp_expression_predicts_functional_state.R)

3. Test associations between individual candidate AMPs and functional modules. [31_candidate_amp_function_association.R](scripts/31_candidate_amp_function_association.R)
