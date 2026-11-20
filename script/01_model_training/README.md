# AMP classifier training

The training scripts take labeled peptide sequences and produce a fine-tuned ESM2 encoder and a ten-fold LightGBM ensemble.

## Analysis steps

1. Split the training pool into encoder fine-tuning and classifier-training subsets. [split_traindata.py](split_traindata.py)

2. Fine-tune the encoder and extract peptide embeddings. [1_finetune_esm2.py](1_finetune_esm2.py), [2_finetuned_esm2_embedding.py](2_finetuned_esm2_embedding.py)

3. Train the fold models and evaluate their ensemble on the test set. [3_train_lgbm_10fold.py](3_train_lgbm_10fold.py), [4_ensemble_predict.py](4_ensemble_predict.py)

4. Convert candidate FASTA sequences to a table, embed them and score them with the discovery models. [0_faa_to_csv.py](0_faa_to_csv.py), [discovery_embedding.py](discovery_embedding.py), [discovery_predict.py](discovery_predict.py)
