#!/usr/bin/env python3
"""
ensemble_predict.py

Usage:
 python ensemble_predict.py --model_dirs results/cls_10fold results/mean_10fold --emb emb/test_cls.npy --meta emb/test_meta.csv --out results/ensemble

- model_dirs: one or more directories containing model_fold_*.joblib
- emb: embeddings that match the model type (if you pass multiple model_dirs with different embedding shapes,
       you should pass multiple --emb and --meta matching order; for simplicity we assume same embedding file.)
Outputs:
 - out/ensemble_preds.csv  (ID, Label, Prob_ensemble)
 - out/ensemble_metrics.txt  (AUC/AP)
"""
import argparse
import os
import glob
import numpy as np
import pandas as pd
import joblib
from sklearn.metrics import roc_auc_score, average_precision_score

def load_models_from_dir(d):
    files = sorted(glob.glob(os.path.join(d, "model_fold_*.joblib")))
    models = [joblib.load(f) for f in files]
    return models

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--model_dirs", nargs="+", required=True, help="one or more model directories with model_fold_*.joblib")
    p.add_argument("--emb", required=True, help="embedding npy for test set (matching model input)")
    p.add_argument("--meta", required=True, help="meta csv with ID and Label")
    p.add_argument("--out", required=True)
    args = p.parse_args()

    X = np.load(args.emb)
    meta = pd.read_csv(args.meta)
    if 'Label' not in meta.columns or 'ID' not in meta.columns:
        raise ValueError("meta csv must contain ID and Label columns")
    y = meta['Label'].astype(int).values
    ids = meta['ID'].astype(str).values

    all_model_probs = []  # shape (n_models, N)
    for d in args.model_dirs:
        models = load_models_from_dir(d)
        print(f"Loaded {len(models)} models from {d}")
        for m in models:
            probs = m.predict_proba(X)[:, 1]
            all_model_probs.append(probs)

    all_model_probs = np.array(all_model_probs)  # (M, N)
    # average across models
    avg_probs = all_model_probs.mean(0)

    out_dir = args.out
    os.makedirs(out_dir, exist_ok=True)
    out_df = pd.DataFrame({
        'ID': ids,
        'Label': y,
        'Prob': avg_probs
    })
    out_df.to_csv(os.path.join(out_dir, "ensemble_preds.csv"), index=False)

    auc = roc_auc_score(y, avg_probs)
    ap = average_precision_score(y, avg_probs)
    with open(os.path.join(out_dir, "ensemble_metrics.txt"), "w") as f:
        f.write(f"AUC\t{auc:.6f}\nAP\t{ap:.6f}\n")
    print(f"✔ Ensemble results saved to {out_dir}. AUC={auc:.4f} AP={ap:.4f}")

if __name__ == "__main__":
    main()

