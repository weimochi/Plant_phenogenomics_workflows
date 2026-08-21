#!/usr/bin/env python3
"""Plot ADMIXTURE ancestry proportions and cross-validation errors."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


CLUSTER_COLORS = [
    "#EFC86E",
    "#97C684",
    "#6F9969",
    "#AAB5D5",
    "#5C66A8",
    "#454A74",
    "#E6A4B4",
    "#C97C7C",
    "#8FB8B8",
    "#D8C3A5",
]

CV_PATTERN = re.compile(
    r"CV error\s*\(K=(\d+)\):\s*([0-9.eE+-]+)"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create one ancestry bar plot per ADMIXTURE K and summarize "
            "cross-validation errors."
        )
    )
    parser.add_argument("--base-dir", required=True, type=Path)
    parser.add_argument(
        "--prefix",
        required=True,
        help="Shared prefix, e.g. brassica_cohort.admixture",
    )
    parser.add_argument("--k-min", required=True, type=int)
    parser.add_argument("--k-max", required=True, type=int)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--fam-col",
        type=int,
        default=1,
        help="Zero-based .fam column containing sample IDs (default: 1/IID)",
    )
    parser.add_argument("--fig-width", type=float, default=18.0)
    parser.add_argument("--fig-height", type=float, default=4.8)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument(
        "--show-sample-labels",
        action="store_true",
        help="Display sample IDs below bars; best for small cohorts",
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        choices=("png", "pdf", "svg"),
        default=("png",),
    )
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if args.k_min < 1 or args.k_max < 1:
        raise ValueError("K values must be positive integers")
    if args.k_min > args.k_max:
        raise ValueError("--k-min must be <= --k-max")
    if args.k_max > len(CLUSTER_COLORS):
        raise ValueError(
            f"K={args.k_max} exceeds the {len(CLUSTER_COLORS)} fixed colors. "
            "Extend CLUSTER_COLORS deliberately before plotting larger K."
        )
    if args.fam_col < 0:
        raise ValueError("--fam-col must be zero or greater")


def read_sample_ids(fam_path: Path, fam_col: int) -> list[str]:
    if not fam_path.is_file():
        raise FileNotFoundError(f"FAM file not found: {fam_path}")

    fam = pd.read_csv(fam_path, sep=r"\s+", header=None, dtype=str)
    if fam_col >= fam.shape[1]:
        raise ValueError(
            f"--fam-col {fam_col} is outside the {fam.shape[1]}-column FAM file"
        )

    sample_ids = fam.iloc[:, fam_col].astype(str).tolist()
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("Sample IDs selected from the FAM file are not unique")
    return sample_ids


def read_q_matrix(q_path: Path, sample_ids: list[str], k: int) -> pd.DataFrame:
    if not q_path.is_file():
        raise FileNotFoundError(f"Q file not found: {q_path}")

    q_matrix = pd.read_csv(q_path, sep=r"\s+", header=None)
    if q_matrix.shape != (len(sample_ids), k):
        raise ValueError(
            f"Unexpected Q dimensions for K={k}: {q_matrix.shape}; "
            f"expected ({len(sample_ids)}, {k})"
        )

    q_matrix = q_matrix.apply(pd.to_numeric, errors="raise")
    row_sums = q_matrix.sum(axis=1).to_numpy()
    if not np.allclose(row_sums, 1.0, atol=1e-3):
        raise ValueError(f"Ancestry proportions do not sum to 1 in {q_path}")

    q_matrix.index = sample_ids
    q_matrix.columns = list(range(k))
    return q_matrix


def sort_by_major_cluster(q_matrix: pd.DataFrame) -> pd.DataFrame:
    cluster_columns = list(q_matrix.columns)
    sortable = q_matrix.copy()
    sortable["major_cluster"] = sortable.idxmax(axis=1)
    sortable["major_proportion"] = sortable[cluster_columns].max(axis=1)

    # Place ancestry components from Cluster 1 to Cluster K. Within each
    # component, place the least-admixed individuals first by sorting their
    # membership in the assigned major component from high to low.
    return sortable.sort_values(
        ["major_cluster", "major_proportion"],
        ascending=[True, False],
        kind="stable",
    )[cluster_columns]


def parse_cv_error(log_path: Path, expected_k: int) -> float | None:
    if not log_path.is_file():
        return None

    for line in log_path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = CV_PATTERN.search(line)
        if match and int(match.group(1)) == expected_k:
            return float(match.group(2))
    return None


def plot_q_matrix(
    q_matrix: pd.DataFrame,
    k: int,
    cv_error: float | None,
    args: argparse.Namespace,
) -> None:
    sorted_q = sort_by_major_cluster(q_matrix)
    sample_ids = sorted_q.index.tolist()
    positions = np.arange(len(sample_ids))
    bottom = np.zeros(len(sample_ids))

    fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height))
    for cluster in range(k):
        values = sorted_q[cluster].to_numpy()
        ax.bar(
            positions,
            values,
            bottom=bottom,
            width=1.0,
            color=CLUSTER_COLORS[cluster],
            edgecolor="none",
            label=f"Cluster {cluster + 1}",
        )
        bottom += values

    title = f"ADMIXTURE K={k}"
    if cv_error is not None:
        title += f" | CV error={cv_error:.6g}"

    ax.set_title(title)
    ax.set_ylabel("Ancestry proportion")
    ax.set_ylim(0, 1)
    ax.set_xlim(-0.5, len(sample_ids) - 0.5)
    ax.set_yticks(np.linspace(0, 1, 6))
    ax.spines[["top", "right"]].set_visible(False)

    if args.show_sample_labels:
        ax.set_xticks(positions)
        ax.set_xticklabels(sample_ids, rotation=90, fontsize=6)
        ax.set_xlabel("Sample")
    else:
        ax.set_xticks([])
        ax.set_xlabel("Individuals sorted by major ancestry component")

    ax.legend(
        loc="upper left",
        bbox_to_anchor=(1.005, 1),
        frameon=False,
        fontsize=8,
    )
    fig.tight_layout()

    for output_format in args.formats:
        output_path = args.output_dir / f"{args.prefix}.K{k}.ancestry.{output_format}"
        fig.savefig(output_path, dpi=args.dpi, bbox_inches="tight")
        print(f"Saved: {output_path}")
    plt.close(fig)


def save_cv_summary(cv_errors: dict[int, float], args: argparse.Namespace) -> None:
    if not cv_errors:
        print("No CV errors found; skipping CV summary plot")
        return

    summary = pd.DataFrame(
        sorted(cv_errors.items()), columns=["K", "CV_error"]
    )
    summary_path = args.output_dir / f"{args.prefix}.cv_errors.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    print(f"Saved: {summary_path}")

    best_row = summary.loc[summary["CV_error"].idxmin()]
    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    ax.plot(
        summary["K"],
        summary["CV_error"],
        marker="o",
        color="#454A74",
        linewidth=1.8,
    )
    ax.scatter(
        [best_row["K"]],
        [best_row["CV_error"]],
        color="#C97C7C",
        zorder=3,
        label=f"Minimum: K={int(best_row['K'])}",
    )
    ax.set_xlabel("Number of ancestral clusters (K)")
    ax.set_ylabel("Cross-validation error")
    ax.set_xticks(summary["K"])
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False)
    fig.tight_layout()

    for output_format in args.formats:
        output_path = args.output_dir / f"{args.prefix}.cv_error.{output_format}"
        fig.savefig(output_path, dpi=args.dpi, bbox_inches="tight")
        print(f"Saved: {output_path}")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    validate_args(args)
    args.base_dir = args.base_dir.expanduser().resolve()
    args.output_dir = args.output_dir.expanduser().resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    fam_path = args.base_dir / f"{args.prefix}.fam"
    sample_ids = read_sample_ids(fam_path, args.fam_col)
    cv_errors: dict[int, float] = {}

    for k in range(args.k_min, args.k_max + 1):
        q_path = args.base_dir / f"{args.prefix}.{k}.Q"
        log_path = args.base_dir / f"{args.prefix}.K{k}.log"
        q_matrix = read_q_matrix(q_path, sample_ids, k)
        cv_error = parse_cv_error(log_path, k)
        if cv_error is not None:
            cv_errors[k] = cv_error
        plot_q_matrix(q_matrix, k, cv_error, args)

    save_cv_summary(cv_errors, args)


if __name__ == "__main__":
    main()
