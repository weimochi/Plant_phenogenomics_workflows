#!/usr/bin/env python3
"""Convert a sample-by-feature table into a labeled distance matrix."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import pairwise_distances
from sklearn.preprocessing import normalize


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--sample-axis", choices=("rows", "columns"), required=True)
    parser.add_argument(
        "--sample-id-column",
        default="sample",
        help="Sample-ID column when --sample-axis rows (default: sample)",
    )
    parser.add_argument(
        "--exclude-columns",
        default="",
        help="Comma-separated metadata columns to exclude in row-oriented input",
    )
    parser.add_argument(
        "--feature-prefix",
        default=None,
        help="Optional prefix selecting feature columns, e.g. feat_",
    )
    parser.add_argument(
        "--metric",
        choices=("cosine", "euclidean", "correlation", "manhattan"),
        required=True,
    )
    parser.add_argument(
        "--transform",
        choices=("none", "zscore", "l2"),
        default="none",
    )
    parser.add_argument(
        "--missing",
        choices=("error", "feature-mean"),
        default="error",
    )
    parser.add_argument("--output-prefix", required=True, type=Path)
    return parser.parse_args()


def read_features(args: argparse.Namespace) -> tuple[list[str], pd.DataFrame]:
    if not args.input.is_file():
        raise FileNotFoundError(f"Input table not found: {args.input}")

    table = pd.read_csv(args.input, sep=None, engine="python")
    if args.sample_axis == "rows":
        if args.sample_id_column not in table.columns:
            raise ValueError(f"Missing sample-ID column: {args.sample_id_column}")
        samples = table[args.sample_id_column].astype(str).str.strip().tolist()
        excluded = {
            value.strip()
            for value in args.exclude_columns.split(",")
            if value.strip()
        }
        excluded.add(args.sample_id_column)
        columns = [column for column in table.columns if column not in excluded]
        if args.feature_prefix is not None:
            columns = [
                column for column in columns if str(column).startswith(args.feature_prefix)
            ]
        features = table[columns].copy()
    else:
        if table.shape[1] < 3:
            raise ValueError("Column-oriented input needs one feature-ID column and >=2 samples")
        samples = [str(column).strip() for column in table.columns[1:]]
        features = table.iloc[:, 1:].T.copy()
        features.columns = table.iloc[:, 0].astype(str).tolist()

    if len(samples) < 3:
        raise ValueError("At least three samples are required")
    if len(samples) != len(set(samples)):
        raise ValueError("Sample IDs are not unique")
    if not len(features.columns):
        raise ValueError("No feature columns were selected")

    features = features.apply(pd.to_numeric, errors="coerce")
    all_missing = features.columns[features.isna().all()].tolist()
    if all_missing:
        raise ValueError(f"Features contain no numeric values: {all_missing[:10]}")

    if features.isna().any().any():
        if args.missing == "error":
            raise ValueError("Feature matrix contains missing/non-numeric values")
        features = features.fillna(features.mean(axis=0))

    features.index = samples
    return samples, features


def transform_features(
    features: pd.DataFrame, transform: str
) -> tuple[np.ndarray, list[str], list[str]]:
    variances = features.var(axis=0, ddof=0)
    zero_variance = variances.index[np.isclose(variances, 0)].astype(str).tolist()
    retained = features.loc[:, ~np.isclose(variances, 0)].copy()
    if retained.shape[1] == 0:
        raise ValueError("All selected features have zero variance")

    values = retained.to_numpy(dtype=float)
    if transform == "zscore":
        values = (values - values.mean(axis=0)) / values.std(axis=0, ddof=0)
    elif transform == "l2":
        zero_norm_samples = np.isclose(np.linalg.norm(values, axis=1), 0)
        if zero_norm_samples.any():
            bad = retained.index[zero_norm_samples].tolist()
            raise ValueError(f"Cannot L2-normalize zero vectors: {bad[:10]}")
        values = normalize(values, norm="l2", axis=1)

    return values, retained.columns.astype(str).tolist(), zero_variance


def main() -> None:
    args = parse_args()
    args.input = args.input.expanduser().resolve()
    args.output_prefix = args.output_prefix.expanduser().resolve()
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)

    samples, features = read_features(args)
    values, retained_features, removed_features = transform_features(
        features, args.transform
    )

    if args.metric in {"cosine", "correlation"}:
        zero_norm = np.isclose(np.linalg.norm(values, axis=1), 0)
        if zero_norm.any():
            bad = np.asarray(samples)[zero_norm].tolist()
            raise ValueError(f"Metric is undefined for zero vectors: {bad[:10]}")

    distance = pairwise_distances(values, metric=args.metric)
    distance = (distance + distance.T) / 2.0
    np.fill_diagonal(distance, 0.0)
    if not np.isfinite(distance).all() or (distance < -1e-10).any():
        raise ValueError("Calculated distance matrix contains invalid values")

    matrix_path = Path(f"{args.output_prefix}.distance.tsv")
    pd.DataFrame(distance, index=samples, columns=samples).to_csv(
        matrix_path, sep="\t", index_label="sample"
    )

    metadata = {
        "input": str(args.input),
        "sample_axis": args.sample_axis,
        "metric": args.metric,
        "transform": args.transform,
        "missing_policy": args.missing,
        "n_samples": len(samples),
        "n_features_input": features.shape[1],
        "n_features_retained": len(retained_features),
        "zero_variance_features_removed": removed_features,
        "output_matrix": str(matrix_path),
    }
    metadata_path = Path(f"{args.output_prefix}.metadata.json")
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Saved: {matrix_path}")
    print(f"Saved: {metadata_path}")


if __name__ == "__main__":
    main()
