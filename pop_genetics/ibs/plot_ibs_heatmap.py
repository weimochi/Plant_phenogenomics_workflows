#!/usr/bin/env python3
"""Cluster a PLINK IBS matrix and draw a consistently reordered heatmap."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform


HEATMAP_COLORS = [
    "#FFFFFF",
    "#E8ECF5",
    "#B8C4D9",
    "#6F7FA2",
    "#384157",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Hierarchically cluster PLINK IBS similarity or distance and "
            "draw a reordered similarity heatmap."
        )
    )
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--id-file", required=True, type=Path)
    parser.add_argument(
        "--matrix-kind",
        required=True,
        choices=("similarity", "distance"),
        help="Use similarity for .mibs and distance for .mdist",
    )
    parser.add_argument(
        "--id-col",
        type=int,
        default=1,
        help="Zero-based ID-file column containing sample IDs (default: 1/IID)",
    )
    parser.add_argument(
        "--linkage-method",
        choices=("average", "complete", "single"),
        default="average",
    )
    parser.add_argument("--output-prefix", required=True, type=Path)
    parser.add_argument("--fig-size", type=float, default=11.0)
    parser.add_argument("--dpi", type=int, default=350)
    parser.add_argument(
        "--hide-labels",
        action="store_true",
        help="Hide sample labels for very large cohorts",
    )
    return parser.parse_args()


def read_sample_ids(id_path: Path, id_col: int) -> list[str]:
    if not id_path.is_file():
        raise FileNotFoundError(f"ID file not found: {id_path}")
    ids = pd.read_csv(id_path, sep=r"\s+", header=None, dtype=str)
    if id_col < 0 or id_col >= ids.shape[1]:
        raise ValueError(
            f"--id-col {id_col} is outside the {ids.shape[1]}-column ID file"
        )
    sample_ids = ids.iloc[:, id_col].astype(str).tolist()
    if len(sample_ids) < 2:
        raise ValueError("At least two samples are required")
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("Selected sample IDs are not unique")
    return sample_ids


def read_matrix(matrix_path: Path, sample_ids: list[str]) -> np.ndarray:
    if not matrix_path.is_file():
        raise FileNotFoundError(f"Matrix file not found: {matrix_path}")
    matrix = pd.read_csv(matrix_path, sep=r"\s+", header=None).to_numpy(dtype=float)
    expected_shape = (len(sample_ids), len(sample_ids))
    if matrix.shape != expected_shape:
        raise ValueError(
            f"Matrix shape is {matrix.shape}; expected {expected_shape} from ID file"
        )
    if not np.isfinite(matrix).all():
        raise ValueError("Matrix contains NaN or infinite values")
    if not np.allclose(matrix, matrix.T, atol=1e-8):
        max_difference = np.max(np.abs(matrix - matrix.T))
        raise ValueError(
            f"Matrix is not symmetric; maximum difference={max_difference:.6g}"
        )
    return matrix


def prepare_matrices(
    matrix: np.ndarray, matrix_kind: str
) -> tuple[np.ndarray, np.ndarray]:
    tolerance = 1e-6
    if matrix_kind == "similarity":
        if matrix.min() < -tolerance or matrix.max() > 1 + tolerance:
            raise ValueError("IBS similarity values must lie between 0 and 1")
        similarity = np.clip(matrix, 0, 1)
        distance = 1.0 - similarity
    else:
        if matrix.min() < -tolerance or matrix.max() > 1 + tolerance:
            raise ValueError("IBS distance values must lie between 0 and 1")
        distance = np.clip(matrix, 0, 1)
        similarity = 1.0 - distance

    expected_diagonal = 1.0 if matrix_kind == "similarity" else 0.0
    if not np.allclose(np.diag(matrix), expected_diagonal, atol=1e-5):
        raise ValueError(
            f"Unexpected matrix diagonal for {matrix_kind}; expected "
            f"{expected_diagonal:g}"
        )

    # Remove negligible floating-point asymmetry and enforce an exact zero
    # diagonal before converting to SciPy condensed-distance form.
    distance = (distance + distance.T) / 2.0
    np.fill_diagonal(distance, 0.0)
    similarity = (similarity + similarity.T) / 2.0
    np.fill_diagonal(similarity, 1.0)
    return similarity, distance


def build_linkage(distance: np.ndarray, method: str) -> np.ndarray:
    condensed = squareform(distance, checks=True)
    return linkage(condensed, method=method, optimal_ordering=True)


def save_tables(
    similarity: np.ndarray,
    sample_ids: list[str],
    order: np.ndarray,
    output_prefix: Path,
) -> pd.DataFrame:
    ordered_ids = [sample_ids[index] for index in order]
    ordered_similarity = similarity[np.ix_(order, order)]
    ordered_df = pd.DataFrame(
        ordered_similarity, index=ordered_ids, columns=ordered_ids
    )

    order_path = Path(f"{output_prefix}.sample_order.tsv")
    pd.DataFrame(
        {
            "plot_position": np.arange(1, len(ordered_ids) + 1),
            "sample": ordered_ids,
        }
    ).to_csv(order_path, sep="\t", index=False)

    matrix_path = Path(f"{output_prefix}.similarity.reordered.tsv")
    ordered_df.to_csv(matrix_path, sep="\t", index=True, index_label="sample")
    print(f"Saved: {order_path}")
    print(f"Saved: {matrix_path}")
    return ordered_df


def draw_heatmap(
    ordered_df: pd.DataFrame,
    args: argparse.Namespace,
) -> None:
    cmap = LinearSegmentedColormap.from_list("ibs_blue", HEATMAP_COLORS)
    labels = not args.hide_labels
    label_size = max(3.0, min(8.0, 180.0 / len(ordered_df)))

    # Linkage already has optimal leaf ordering. Disable clustermap's internal
    # reordering and display the explicitly ordered matrix on both axes.
    grid = sns.clustermap(
        ordered_df,
        row_cluster=False,
        col_cluster=False,
        cmap=cmap,
        vmin=0,
        vmax=1,
        square=True,
        linewidths=0,
        xticklabels=labels,
        yticklabels=labels,
        figsize=(args.fig_size, args.fig_size),
        cbar_kws={"label": "IBS similarity"},
        cbar_pos=(0.02, 0.80, 0.03, 0.16),
    )

    if labels:
        grid.ax_heatmap.tick_params(axis="x", labelrotation=90, labelsize=label_size)
        grid.ax_heatmap.tick_params(axis="y", labelrotation=0, labelsize=label_size)
    grid.ax_heatmap.set_xlabel("Sample")
    grid.ax_heatmap.set_ylabel("Sample")
    grid.fig.suptitle(
        f"IBS similarity clustered by {args.linkage_method}-linkage distance",
        y=1.01,
    )

    for extension in ("png", "pdf"):
        output_path = Path(f"{args.output_prefix}.heatmap.{extension}")
        grid.fig.savefig(output_path, dpi=args.dpi, bbox_inches="tight")
        print(f"Saved: {output_path}")
    plt.close(grid.fig)


def main() -> None:
    args = parse_args()
    args.matrix = args.matrix.expanduser().resolve()
    args.id_file = args.id_file.expanduser().resolve()
    args.output_prefix = args.output_prefix.expanduser().resolve()
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)

    sample_ids = read_sample_ids(args.id_file, args.id_col)
    raw_matrix = read_matrix(args.matrix, sample_ids)
    similarity, distance = prepare_matrices(raw_matrix, args.matrix_kind)
    linkage_matrix = build_linkage(distance, args.linkage_method)
    order = leaves_list(linkage_matrix)
    ordered_df = save_tables(similarity, sample_ids, order, args.output_prefix)
    draw_heatmap(ordered_df, args)


if __name__ == "__main__":
    main()
