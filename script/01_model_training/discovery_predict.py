#!/usr/bin/env python3
import numpy as np
import pandas as pd
import joblib, json, os, glob

# ===========================
#        配置部分
# ===========================

# 三种 embedding 输入文件
EMB_FILES = {
    "mean": "smorf_mean.npy",
    "mean_max": "smorf_mean_max.npy",
    "cls": "smorf_cls.npy"
}

SEQ_CSV = "smorf.csv"       # 序列文件（必须包含 Sequence 列）

# 每种模型的 LightGBM 模型目录
MODEL_ROOT = "/public2/labmember/wyc_rawdata/amp/model_test/5_esm2lgbm_modified"
MODEL_DIRS = {
    "mean": os.path.join(MODEL_ROOT, "train_mean"),
    "mean_max": os.path.join(MODEL_ROOT, "train_mean_max"),
    "cls": os.path.join(MODEL_ROOT, "train_cls"),
}

OUT_PREFIX = "smorf_pred_lgbm_esm2"   # 输出文件前缀


# ===========================
#    加载 LightGBM 模型
# ===========================
def load_models(path):
    if os.path.isdir(path):
        model_files = sorted(glob.glob(os.path.join(path, "model_fold_*.joblib")))
        if len(model_files) > 0:
            print(f"✔ {path} 检测到 {len(model_files)} 个 fold 模型 -> ensemble")
            return [joblib.load(m) for m in model_files]

        model_files = glob.glob(os.path.join(path, "*.joblib"))
        if len(model_files) == 1:
            print(f"✔ {path} 检测到单模型")
            return [joblib.load(model_files[0])]
        elif len(model_files) > 1:
            print(f"✔ {path} 检测到 {len(model_files)} 个模型 -> ensemble")
            return [joblib.load(m) for m in model_files]

        raise ValueError(f"❌ {path} 内未找到 joblib 模型文件")

    else:
        print(f"✔ 加载单模型: {path}")
        return [joblib.load(path)]


# ===========================
#     ensemble 预测
# ===========================
def predict_ensemble(models, X):
    all_probs = []
    for m in models:
        if hasattr(m, "predict_proba"):
            p = m.predict_proba(X)[:, 1]
        else:
            raw = m.predict(X)
            if raw.ndim == 1:
                raw = raw.astype(float)
            p = raw
        all_probs.append(p)
    return np.vstack(all_probs).mean(axis=0)


# ===========================
#   主程序：三种 embedding
# ===========================
seq_df = pd.read_csv(SEQ_CSV)
if "Sequence" not in seq_df.columns:
    raise ValueError("❌ smorf.csv 必须包含 Sequence 列")

print("====================================")
print("   🚀 开始运行 mean / mean_max / cls 三种模型")
print("====================================\n")

for mode in ["mean", "mean_max", "cls"]:
    print(f"\n========== 🧩 处理 {mode} 模型 ==========")

    emb_file = EMB_FILES[mode]
    model_dir = MODEL_DIRS[mode]

    if not os.path.exists(emb_file):
        print(f"⚠ 跳过 {mode}: 找不到 embedding 文件 {emb_file}")
        continue

    # 加载 embedding
    X = np.load(emb_file)
    if len(X) != len(seq_df):
        raise ValueError(f"❌ {mode} embedding 数量与 CSV 行数不匹配")

    # 加载模型
    models = load_models(model_dir)

    # 预测
    probs = predict_ensemble(models, X)

    # 加载 best threshold
    metrics_path = os.path.join(model_dir, "metrics.json")
    best_thresh = 0.5
    if os.path.exists(metrics_path):
        try:
            with open(metrics_path) as f:
                md = json.load(f)
            best_thresh = float(md.get("best_threshold", 0.5))
            print(f"✔ 使用最佳阈值 best_threshold={best_thresh}")
        except:
            print("⚠ metrics.json 读取失败，使用默认 0.5")

    # 输出文件名
    out_csv = f"{OUT_PREFIX}_{mode}.csv"
    out_pos = f"predicted_50_{mode}.faa"
    out_high = f"predicted_80_{mode}.faa"

    # 生成输出 CSV
    out = seq_df.copy()
    out[f"prob_{mode}"] = probs
    out[f"pred_{mode}"] = (probs >= best_thresh).astype(int)

    out.to_csv(out_csv, index=False)
    print("✔ 输出 CSV:", out_csv)

    # 导出 positive faa
    pos = out[out[f"pred_{mode}"] == 1]
    with open(out_pos, "w") as f:
        for idx, row in pos.iterrows():
            f.write(f">seq_{idx}\n{row['Sequence']}\n")
    print("✔ Positive FASTA:", out_pos)

    # 导出 high confidence faa
    high = out[out[f"prob_{mode}"] >= 0.8]
    with open(out_high, "w") as f:
        for idx, row in high.iterrows():
            f.write(f">highconf_{idx}\n{row['Sequence']}\n")
    print("✔ High-confidence FASTA:", out_high)

print("\n🎉 全部模型运行完成！")

