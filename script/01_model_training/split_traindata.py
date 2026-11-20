#!/usr/bin/env python3
import pandas as pd
from sklearn.model_selection import train_test_split

def main():
    input_file = "train.csv"     # 原始训练集
    output_10 = "finetune_10.csv"   # 10% 数据（给 ESM2 finetune）
    output_90 = "train_90.csv"   # 90% 剩余训练集

    # 读取原始数据
    df = pd.read_csv(input_file)

    # 按 10% 切分（stratify 保持 label 比例一致）
    df_90, df_10 = train_test_split(
        df, 
        test_size=0.10, 
        random_state=42, 
        stratify=df["Label"]
    )

    # 保存结果
    df_10.to_csv(output_10, index=False)
    df_90.to_csv(output_90, index=False)

    print(f"Done! 总样本数: {len(df)}")
    print(f"finetune_10.csv: {len(df_10)}")
    print(f"train_90.csv: {len(df_90)}")

if __name__ == "__main__":
    main()

