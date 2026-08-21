#!/usr/bin/env python3
"""Summarize and plot PLINK --het inbreeding-coefficient estimates."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize and plot the F column from a PLINK .het file."
    )
    parser.add_argument("--het", required=True, type=Path)
    parser.add_argument("--output-prefix", required=True, type=Path)
    parser.add_argument("--bins", type=int, default=30)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument(
        "--label-extremes",
        type=int,
        default=5,
        help="Label this many samples at each end of ranked F plot (default: 5)",
    )
    return parser.parse_args()


def read_het(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"PLINK .het file not found: {path}")

    frame = pd.read_csv(path, sep=r"\s+")
    required = {"FID", "IID", "O(HOM)", "E(HOM)", "N(NM)", "F"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Missing required .het columns: {sorted(missing)}")

    if frame.empty:
        raise ValueError("PLINK .het file has no sample rows")
    if frame["IID"].astype(str).duplicated().any():
        raise ValueError("IID values are not unique")

    numeric_columns = ["O(HOM)", "E(HOM)", "N(NM)", "F"]
    for column in numeric_columns:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")

    invalid = frame[numeric_columns].isna().any(axis=1)
    if invalid.any():
        samples = frame.loc[invalid, "IID"].astype(str).tolist()
        raise ValueError(
            "Non-numeric or missing PLINK estimates for samples: "
            + ", ".join(samples[:10])
        )

    if (frame["N(NM)"] <= 0).any():
        raise ValueError("At least one sample has no nonmissing genotype observations")

    frame["IID"] = frame["IID"].astype(str)
    frame["observed_heterozygous"] = frame["N(NM)"] - frame["O(HOM)"]
    frame["expected_heterozygous"] = frame["N(NM)"] - frame["E(HOM)"]
    frame["observed_heterozygosity_rate"] = (
        frame["observed_heterozygous"] / frame["N(NM)"]
    )
    frame["expected_heterozygosity_rate"] = (
        frame["expected_heterozygous"] / frame["N(NM)"]
    )
    return frame


def save_tables(frame: pd.DataFrame, output_prefix: Path) -> None:
    samples_path = Path(f"{output_prefix}.samples.tsv")
    frame.sort_values("F").to_csv(samples_path, sep="\t", index=False)

    summary = pd.DataFrame(
        {
            "metric": [
                "n_samples",
                "mean_F",
                "median_F",
                "sd_F",
                "min_F",
                "q25_F",
                "q75_F",
                "max_F",
                "n_F_below_zero",
                "n_F_above_zero",
            ],
            "value": [
                len(frame),
                frame["F"].mean(),
                frame["F"].median(),
                frame["F"].std(ddof=1),
                frame["F"].min(),
                frame["F"].quantile(0.25),
                frame["F"].quantile(0.75),
                frame["F"].max(),
                int((frame["F"] < 0).sum()),
                int((frame["F"] > 0).sum()),
            ],
        }
    )
    summary_path = Path(f"{output_prefix}.summary.tsv")
    summary.to_csv(summary_path, sep="\t", index=False)
    print(f"Saved: {samples_path}")
    print(f"Saved: {summary_path}")


def plot_distribution(frame: pd.DataFrame, args: argparse.Namespace) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.histplot(
        data=frame,
        x="F",
        bins=args.bins,
        color="#5C66A8",
        edgecolor="white",
        ax=ax,
    )
    ax.axvline(0, color="#C97C7C", linestyle="--", linewidth=1.5, label="F = 0")
    ax.axvline(
        frame["F"].median(),
        color="#6F9969",
        linestyle="-",
        linewidth=1.5,
        label=f"Median = {frame['F'].median():.3f}",
    )
    ax.set_xlabel("PLINK method-of-moments inbreeding coefficient (F)")
    ax.set_ylabel("Number of samples")
    ax.set_title("Distribution of individual inbreeding coefficients")
    ax.legend(frameon=False)
    fig.tight_layout()

    for extension in ("png", "pdf"):
        output = Path(f"{args.output_prefix}.distribution.{extension}")
        fig.savefig(output, dpi=args.dpi, bbox_inches="tight")
        print(f"Saved: {output}")
    plt.close(fig)


def plot_ranked(frame: pd.DataFrame, args: argparse.Namespace) -> None:
    ordered = frame.sort_values("F").reset_index(drop=True)
    positions = np.arange(1, len(ordered) + 1)

    fig_width = max(9.0, min(20.0, len(ordered) * 0.06))
    fig, ax = plt.subplots(figsize=(fig_width, 6))
    colors = np.where(ordered["F"] >= 0, "#5C66A8", "#C97C7C")
    ax.scatter(positions, ordered["F"], c=colors, s=24, alpha=0.9)
    ax.axhline(0, color="#444444", linestyle="--", linewidth=1)
    ax.set_xlabel("Samples ranked from lowest to highest F")
    ax.set_ylabel("PLINK method-of-moments F")
    ax.set_title("Ranked individual inbreeding coefficients")
    ax.spines[["top", "right"]].set_visible(False)

    label_count = min(args.label_extremes, len(ordered) // 2)
    label_indices = list(range(label_count)) + list(
        range(len(ordered) - label_count, len(ordered))
    )
    for index in sorted(set(label_indices)):
        row = ordered.iloc[index]
        ax.annotate(
            row["IID"],
            (positions[index], row["F"]),
            xytext=(0, 5),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=7,
            rotation=45,
        )

    fig.tight_layout()
    for extension in ("png", "pdf"):
        output = Path(f"{args.output_prefix}.ranked.{extension}")
        fig.savefig(output, dpi=args.dpi, bbox_inches="tight")
        print(f"Saved: {output}")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    if args.bins < 1:
        raise ValueError("--bins must be at least 1")
    if args.label_extremes < 0:
        raise ValueError("--label-extremes must be zero or greater")

    args.het = args.het.expanduser().resolve()
    args.output_prefix = args.output_prefix.expanduser().resolve()
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)

    frame = read_het(args.het)
    save_tables(frame, args.output_prefix)
    plot_distribution(frame, args)
    plot_ranked(frame, args)


if __name__ == "__main__":
    main()
