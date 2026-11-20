#!/usr/bin/env python3
"""
Compatible with ALL LightGBM versions (3.x / 4.x)
"""

import argparse
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, average_precision_score
import joblib
import os
import lightgbm as lgb


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--emb", required=True)
    p.add_argument("--meta", required=True)
    p.add_argument("--out_dir", required=True)
    p.add_argument("--n_splits", type=int, default=10)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--num_leaves", type=int, default=31)
    p.add_argument("--n_estimators", type=int, default=1000)
    args = p.parse_args()

    X = np.load(args.emb)
    meta = pd.read_csv(args.meta)

    # Label column
    if "Label" in meta.columns:
        label_col = "Label"
    elif "label" in meta.columns:
        label_col = "label"
    else:
        raise ValueError("meta csv must contain Label or label column")

    # ID exists or generate
    if "ID" not in meta.columns:
        meta = meta.copy()
        meta["ID"] = [f"seq_{i:06d}" for i in range(len(meta))]
        print("⚠️ 'ID' column not found — auto-generated IDs seq_000000 ...")
    meta["ID"] = meta["ID"].astype(str)

    y = meta[label_col].astype(int).values
    ids = meta["ID"].astype(str).values

    if len(X) != len(y):
        raise AssertionError("embedding rows != meta rows")

    os.makedirs(args.out_dir, exist_ok=True)

    skf = StratifiedKFold(
        n_splits=args.n_splits, shuffle=True, random_state=args.seed
    )

    aucs, aps = [], []
    fold = 0

    for train_idx, val_idx in skf.split(X, y):
        print(f"→ Fold {fold}")

        X_train, X_val = X[train_idx], X[val_idx]
        y_train, y_val = y[train_idx], y[val_idx]

        clf = lgb.LGBMClassifier(
            objective="binary",
            boosting_type="gbdt",
            num_leaves=args.num_leaves,
            n_estimators=args.n_estimators,
            learning_rate=0.05,
            random_state=args.seed + fold,
            n_jobs=-1,
        )

        # ⭐ Works on ALL versions of LightGBM
        clf.fit(
            X_train, y_train,
            eval_set=[(X_val, y_val)],
            eval_metric="auc",
            callbacks=[lgb.early_stopping(50)]
        )

        # save model
        joblib.dump(clf, os.path.join(args.out_dir, f"model_fold_{fold}.joblib"))

        # validate
        probs = clf.predict_proba(X_val)[:, 1]
        auc = roc_auc_score(y_val, probs)
        ap = average_precision_score(y_val, probs)
        aucs.append(auc)
        aps.append(ap)
        print(f"  Fold {fold}  AUC={auc:.4f}  AP={ap:.4f}")

        # save predictions
        pd.DataFrame({
            "ID": ids[val_idx],
            "Label": y_val,
            "Prob": probs
        }).to_csv(os.path.join(args.out_dir, f"preds_fold_{fold}.csv"), index=False)

        fold += 1

    pd.DataFrame({
        "metric": ["AUC", "AP"],
        "mean": [np.mean(aucs), np.mean(aps)],
        "std": [np.std(aucs), np.std(aps)],
    }).to_csv(os.path.join(args.out_dir, "metrics_summary.csv"), index=False)

    print("✔ Training complete.")


if __name__ == "__main__":
    main()

