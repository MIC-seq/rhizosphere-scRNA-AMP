#!/usr/bin/env python3
"""
FP16-accelerated ESM2 embedding generator
"""

import argparse
import torch
from torch.amp import autocast
import esm
import pandas as pd
import numpy as np
from tqdm import tqdm
import re
import os

def sanitize_seq(s):
    return re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', 'X', s.upper())

def load_finetuned_esm(ckpt_path):
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model.eval()
    if ckpt_path is None:
        return model, alphabet

    sd = torch.load(ckpt_path, map_location='cpu')
    esm_sd = {}

    for k, v in sd.items():
        if k.startswith('esm.'):
            esm_sd[k[len('esm.'):]] = v
        else:
            idx = k.find('esm.')
            if idx >= 0:
                newk = k[idx+len('esm.'):]
                esm_sd[newk] = v

    if len(esm_sd) > 0:
        model_sd = model.state_dict()
        for k_in, v in list(esm_sd.items()):
            if k_in in model_sd:
                model_sd[k_in].copy_(v)
        model.load_state_dict(model_sd)
        print(f"✔ Loaded ESM weights from {ckpt_path}")
    else:
        print("⚠️ No 'esm.' keys found — using pretrained model.")

    return model, alphabet

def main():
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--input_csv", required=True)
    p.add_argument("--ckpt", required=False, default="/public2/labmember/wyc_rawdata/amp/REF_DBs/5_esm2lgbm_modified/finetuned_esm2.pth")
    p.add_argument("-o", "--out_prefix", required=True)
    p.add_argument("--batch_size", type=int, default=16)
    p.add_argument("--max_len", type=int, default=1022)
    args = p.parse_args()

    df = pd.read_csv(args.input_csv)
    if "Sequence" not in df.columns:
        raise ValueError("Input CSV must contain 'Sequence' column")

    ids = df['ID'].astype(str).tolist() if 'ID' in df.columns else [f"seq_{i:06d}" for i in range(len(df))]
    labels = df['Label'].tolist() if 'Label' in df.columns else [None] * len(df)
    seqs = [sanitize_seq(str(s))[:args.max_len] for s in df["Sequence"].astype(str)]

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model, alphabet = load_finetuned_esm(args.ckpt)
    model = model.to(device)
    model.eval()
    batch_converter = alphabet.get_batch_converter()

    N = len(seqs)
    D = 1280
    emb_cls = np.zeros((N, D), dtype=np.float32)
    emb_mean = np.zeros((N, D), dtype=np.float32)
    emb_mean_max = np.zeros((N, D * 2), dtype=np.float32)

    ptr = 0

    for i in tqdm(range(0, N, args.batch_size), desc="Embedding (FP16)"):
        batch_ids = ids[i:i+args.batch_size]
        batch_seqs = seqs[i:i+args.batch_size]
        batch = [(bid, seq) for bid, seq in zip(batch_ids, batch_seqs)]

        _, _, tokens = batch_converter(batch)
        tokens = tokens.to(device)

        with torch.no_grad():
            # ------------ FP16 Autocast 加速 -----------------
            with autocast("cuda",dtype=torch.float16):
                out = model(tokens, repr_layers=[33])
            # -------------------------------------------------

        reps = out["representations"][33]   # (B, L, 1280)

        cls_emb = reps[:, 0, :].cpu().numpy()
        seq_reps = reps[:, 1:, :]

        mean_emb = seq_reps.mean(1).cpu().numpy()
        max_emb = seq_reps.amax(1).cpu().numpy()

        mean_max = np.concatenate([mean_emb, max_emb], axis=1)

        B = cls_emb.shape[0]
        emb_cls[ptr:ptr+B] = cls_emb
        emb_mean[ptr:ptr+B] = mean_emb
        emb_mean_max[ptr:ptr+B] = mean_max
        ptr += B

    out_dir = os.path.dirname(args.out_prefix)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    np.save(f"{args.out_prefix}_cls.npy", emb_cls)
    np.save(f"{args.out_prefix}_mean.npy", emb_mean)
    np.save(f"{args.out_prefix}_mean_max.npy", emb_mean_max)
    np.save(f"{args.out_prefix}_ids.npy", np.array(ids))

    meta = pd.DataFrame({"ID": ids, "Label": labels})
    meta.to_csv(f"{args.out_prefix}_meta.csv", index=False)

    print("✔ Saved embeddings.")
    print(f"  {args.out_prefix}_cls.npy")
    print(f"  {args.out_prefix}_mean.npy")
    print(f"  {args.out_prefix}_mean_max.npy")
    print(f"  {args.out_prefix}_meta.csv")

if __name__ == "__main__":
    main()

