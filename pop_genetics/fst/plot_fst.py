#!/usr/bin/env python3
"""Summarize and plot pairwise windowed FST files."""
import argparse
import re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def natural_key(value):
    return [int(x) if x.isdigit() else x.lower() for x in re.split(r"(\d+)", str(value))]

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input-dir", required=True, type=Path)
    p.add_argument("--pattern", default="*.windowed.weir.fst")
    p.add_argument("--fst-column", choices=["WEIGHTED_FST", "MEAN_FST"], default="WEIGHTED_FST")
    p.add_argument("--minimum-variants", type=int, default=1)
    p.add_argument("--contigs-file", type=Path)
    p.add_argument("--output-prefix", required=True, type=Path)
    args = p.parse_args()
    files = sorted(args.input_dir.glob(args.pattern))
    if not files or args.minimum_variants < 1:
        raise ValueError("No input files found or invalid minimum variants")
    frames = []
    for path in files:
        df = pd.read_csv(path, sep="\t")
        needed = {"CHROM", "BIN_START", "BIN_END", "N_VARIANTS", args.fst_column}
        if not needed.issubset(df.columns):
            raise ValueError(f"Missing columns in {path}: {sorted(needed - set(df.columns))}")
        df["FST"] = pd.to_numeric(df[args.fst_column], errors="coerce")
        df["N_VARIANTS"] = pd.to_numeric(df["N_VARIANTS"], errors="coerce")
        df = df.dropna(subset=["FST", "N_VARIANTS", "BIN_START", "BIN_END"])
        df = df[df.N_VARIANTS >= args.minimum_variants].copy()
        df["comparison"] = path.name.removesuffix(".windowed.weir.fst")
        frames.append(df)
    data = pd.concat(frames, ignore_index=True)
    if data.empty:
        raise ValueError("No windows remain after filtering")
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    summary = data.groupby("comparison").agg(n_windows=("FST", "size"), mean_fst=("FST", "mean"), median_fst=("FST", "median"), q95_fst=("FST", lambda x: x.quantile(.95)), min_fst=("FST", "min"), max_fst=("FST", "max")).reset_index()
    summary.to_csv(f"{args.output_prefix}.summary.tsv", sep="\t", index=False)

    order = summary.sort_values("median_fst").comparison
    fig, ax = plt.subplots(figsize=(max(8, len(order)*1.3), 6))
    sns.violinplot(data=data, x="comparison", y="FST", order=order, inner=None, cut=0, color="#B8C4D9", ax=ax)
    sns.boxplot(data=data, x="comparison", y="FST", order=order, width=.18, showfliers=False, color="white", ax=ax)
    ax.tick_params(axis="x", rotation=35); ax.set_xlabel(""); ax.set_ylabel(args.fst_column.replace("_", " "))
    fig.tight_layout()
    for ext in ("png", "pdf"): fig.savefig(f"{args.output_prefix}.distribution.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)

    if args.contigs_file:
        contigs = [x.strip() for x in args.contigs_file.read_text().splitlines() if x.strip()]
        data = data[data.CHROM.astype(str).isin(contigs)].copy()
    else:
        contigs = sorted(data.CHROM.astype(str).unique(), key=natural_key)
    comparisons = sorted(data.comparison.unique())
    ncols = min(3, len(comparisons)); nrows = int(np.ceil(len(comparisons)/ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.5*ncols, 4*nrows), squeeze=False)
    for ax, comparison in zip(axes.flat, comparisons):
        sub = data[data.comparison == comparison].copy(); offset=0; ticks=[]; labels=[]
        for ci, chrom in enumerate(contigs):
            c = sub[sub.CHROM.astype(str) == chrom].copy()
            if c.empty: continue
            c["POS"] = (pd.to_numeric(c.BIN_START)+pd.to_numeric(c.BIN_END))/2 + offset
            ax.scatter(c.POS, c.FST, s=6, color=["#5C66A8", "#AAB5D5"][ci%2], linewidths=0)
            ticks.append((c.POS.min()+c.POS.max())/2); labels.append(chrom); offset=c.POS.max()+1_000_000
        ax.axhline(sub.FST.quantile(.95), color="#C97C7C", ls="--", lw=1)
        ax.set_title(comparison); ax.set_xticks(ticks, labels, rotation=45, ha="right", fontsize=7); ax.set_ylabel("FST")
    for ax in axes.flat[len(comparisons):]: ax.set_visible(False)
    fig.tight_layout()
    for ext in ("png", "pdf"): fig.savefig(f"{args.output_prefix}.genome_scan.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)

if __name__ == "__main__": main()
