#!/usr/bin/env python3
"""
Step 1b: Universal script to merge per-sample HTSeq count files into a gene-by-sample expression matrix.
- Automatically handles union of features (missing values filled with 0).
- Excludes HTSeq special summary rows (__no_feature, __ambiguous, etc.) by default.
"""

import os
import glob
import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge any set of HTSeq count files into a universal expression matrix."
    )
    parser.add_argument(
        "--counts-dir",
        required=True,
        help="Directory containing per-sample *.counts.txt files."
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output file path for the merged TSV matrix (e.g., expression_matrix.tsv)."
    )
    parser.add_argument(
        "--keep-special-rows",
        action="store_true",
        help="Keep HTSeq special rows like __no_feature, __ambiguous, etc."
    )
    return parser.parse_args()


def main():
    args = parse_args()

    counts_dir = os.path.abspath(args.counts_dir)
    out_path = os.path.abspath(args.out)

    if not os.path.isdir(counts_dir):
        raise FileNotFoundError(f"❌ Counts directory not found: {counts_dir}")

    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    # 收集資料夾下所有 counts 檔案
    count_files = sorted(glob.glob(os.path.join(counts_dir, "*.counts.txt")))

    if not count_files:
        raise RuntimeError(f"❌ No .counts.txt files found in: {counts_dir}")

    print(f"[INFO] Counts directory: {counts_dir}")
    print(f"[INFO] Found {len(count_files)} count files.")

    matrix = None

    for path in count_files:
        fname = os.path.basename(path)                  # e.g. sample01.counts.txt
        sample_id = fname.replace(".counts.txt", "")   # -> sample01

        print(f"[INFO] Reading {fname} -> sample_id = {sample_id}")

        df = pd.read_csv(
            path,
            sep=r"\s+",
            header=None,
            names=["feature_id", sample_id],
            dtype={"feature_id": str}
        )

        if not args.keep_special_rows:
            df = df[~df["feature_id"].str.startswith("__")].copy()

        df = df.set_index("feature_id")

        if matrix is None:
            matrix = df
        else:
            matrix = matrix.join(df, how="outer")

    # 缺值補 0，轉成整數
    matrix = matrix.fillna(0).astype(int)
    print(f"[SUCCESS] Merged matrix shape (Genes x Samples): {matrix.shape}")

    # 輸出最終總表
    matrix.to_csv(out_path, sep="\t")
    print(f"[DONE] Universal expression matrix successfully written to: {out_path}")


if __name__ == "__main__":
    main()