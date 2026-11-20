#!/usr/bin/env python3
import torch
import torch.nn as nn
import torch.optim as optim
import esm
import pandas as pd
from torch.utils.data import Dataset, DataLoader
from tqdm import tqdm
import argparse
import numpy as np

# ================================
# 数据集（已修改以适配 Sequence, Length, label）
# ================================
class AMP_Dataset(Dataset):
    def __init__(self, df, alphabet, max_len=1022):

        # --- 新增：如果没有 ID 列，则自动创建 ---
        if "ID" not in df.columns:
            df["ID"] = range(1, len(df)+1)

        # --- 新增：兼容 label / Label 两种形式 ---
        if "Label" in df.columns:
            self.labels = df["Label"].astype(int).tolist()
        elif "label" in df.columns:
            self.labels = df["label"].astype(int).tolist()
        else:
            raise ValueError("CSV 中必须包含 label 或 Label 列。")

        # 序列列名固定为 Sequence
        if "Sequence" not in df.columns:
            raise ValueError("CSV 中必须包含 Sequence 列。")

        self.ids = df["ID"].tolist()
        seqs = df["Sequence"].astype(str)

        # 清洗非法氨基酸字符
        seqs = seqs.str.replace(r"[^ACDEFGHIKLMNPQRSTVWY]", "X", regex=True)
        self.seqs = seqs.tolist()

        self.batch_converter = alphabet.get_batch_converter()
        self.max_len = max_len

    def __len__(self):
        return len(self.ids)

    def __getitem__(self, idx):
        seq = self.seqs[idx][:self.max_len]
        return self.ids[idx], seq, self.labels[idx]

    def collate_fn(self, batch):
        labels = [x[2] for x in batch]
        batch = [(x[0], x[1]) for x in batch]  # (ID, seq)
        _, _, tokens = self.batch_converter(batch)
        return tokens, torch.tensor(labels).long()


# ================================
# Finetune Model
# ================================
class ESM2ForAMP(nn.Module):
    def __init__(self, esm_model, freeze_until=29):
        super().__init__()
        self.esm = esm_model

        # 冻结前29层
        for name, param in self.esm.named_parameters():
            if any(f"layers.{i}." in name for i in range(freeze_until)):
                param.requires_grad = False

        # 分类头
        self.classifier = nn.Sequential(
            nn.Linear(1280, 512),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(512, 2)
        )

    def forward(self, tokens):
        out = self.esm(tokens, repr_layers=[33])
        reps = out["representations"][33]
        cls = reps[:, 0, :]
        logits = self.classifier(cls)
        return logits


# ================================
# Train Loop
# ================================
def train(model, loader, device, epochs, lr):
    model = model.to(device)
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.AdamW(filter(lambda p: p.requires_grad, model.parameters()),
                            lr=lr, weight_decay=0.01)
    scaler = torch.cuda.amp.GradScaler()

    for epoch in range(1, epochs+1):
        model.train()
        total_loss = 0

        loop = tqdm(loader, desc=f"Epoch {epoch}/{epochs}", total=len(loader))
        for tokens, labels in loop:
            tokens = tokens.to(device)
            labels = labels.to(device)

            optimizer.zero_grad()
            with torch.cuda.amp.autocast():
                logits = model(tokens)
                loss = criterion(logits, labels)

            scaler.scale(loss).backward()
            scaler.unscale_(optimizer)
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            scaler.step(optimizer)
            scaler.update()

            total_loss += loss.item()
            loop.set_postfix(loss=loss.item())

        print(f"Epoch {epoch} Loss = {total_loss/len(loader):.4f}")

    print("✔ Finetune complete.")
    return model


# ================================
# Main
# ================================
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--train_csv", required=True)
    parser.add_argument("--out_ckpt", default="finetuned_esm2.pth")
    parser.add_argument("--epochs", type=int, default=2)
    parser.add_argument("--batch_size", type=int, default=8)
    parser.add_argument("--lr", type=float, default=1e-5)
    args = parser.parse_args()

    print("▶ Loading ESM2-650M model...")
    esm_model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    esm_model.eval()

    # --- 修改：读取后统一列名，让 Dataset 能处理 ---
    df = pd.read_csv(args.train_csv)

    # 如果用户用的是原始格式 Sequence,Length,label → 自动重命名
    if "Sequence" not in df.columns and "sequence" in df.columns:
        df.rename(columns={"sequence": "Sequence"}, inplace=True)
    if "Label" not in df.columns and "label" in df.columns:
        df.rename(columns={"label": "Label"}, inplace=True)

    dataset = AMP_Dataset(df, alphabet)
    loader = DataLoader(dataset, batch_size=args.batch_size,
                        shuffle=True, collate_fn=dataset.collate_fn)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print("▶ Using device:", device)

    model = ESM2ForAMP(esm_model)

    model = train(
        model,
        loader,
        device,
        epochs=args.epochs,
        lr=args.lr
    )

    torch.save(model.state_dict(), args.out_ckpt)
    print("✔ Saved finetuned checkpoint:", args.out_ckpt)


if __name__ == "__main__":
    main()

