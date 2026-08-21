#!/usr/bin/env python3
"""Align and compare two or more sample distance matrices."""

from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial import procrustes
from scipy.stats import pearsonr, rankdata, spearmanr


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--permutations", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--pcoa-dimensions", type=int, default=2)
    parser.add_argument("--minimum-samples", type=int, default=4)
    return parser.parse_args()


def read_labeled(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep=None, engine="python", index_col=0)
    frame.index = frame.index.astype(str).str.strip()
    frame.columns = frame.columns.astype(str).str.strip()
    return frame


def read_plink(path: Path, id_file: Path, id_column: int) -> pd.DataFrame:
    ids = pd.read_csv(id_file, sep=r"\s+", header=None, dtype=str)
    if id_column < 0 or id_column >= ids.shape[1]:
        raise ValueError(f"id_column={id_column} is outside {id_file}")
    samples = ids.iloc[:, id_column].astype(str).str.strip().tolist()
    matrix = pd.read_csv(path, sep=r"\s+", header=None)
    if matrix.shape != (len(samples), len(samples)):
        raise ValueError(f"PLINK matrix dimensions do not match {id_file}")
    matrix.index = samples
    matrix.columns = samples
    return matrix


def validate_distance(name: str, frame: pd.DataFrame) -> pd.DataFrame:
    if frame.index.duplicated().any() or frame.columns.duplicated().any():
        raise ValueError(f"{name}: duplicate sample labels")
    if set(frame.index) != set(frame.columns):
        raise ValueError(f"{name}: row and column sample sets differ")
    frame = frame.loc[frame.index, frame.index].apply(pd.to_numeric, errors="coerce")
    values = frame.to_numpy(dtype=float)
    if values.shape[0] < 2 or not np.isfinite(values).all():
        raise ValueError(f"{name}: invalid or nonfinite matrix")
    if not np.allclose(values, values.T, atol=1e-8):
        raise ValueError(f"{name}: matrix is not symmetric")
    if not np.allclose(np.diag(values), 0, atol=1e-8):
        raise ValueError(f"{name}: distance diagonal is not zero")
    if (values < -1e-10).any():
        raise ValueError(f"{name}: distance contains negative values")
    values = np.maximum((values + values.T) / 2.0, 0)
    np.fill_diagonal(values, 0)
    return pd.DataFrame(values, index=frame.index, columns=frame.index)


def load_manifest(path: Path) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    manifest = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    required = {"name", "path", "format", "kind"}
    if not required.issubset(manifest.columns):
        raise ValueError(f"Manifest needs columns: {sorted(required)}")
    if manifest.name.duplicated().any() or len(manifest) < 2:
        raise ValueError("Manifest needs at least two uniquely named matrices")

    matrices = {}
    for row in manifest.itertuples(index=False):
        matrix_path = Path(row.path).expanduser().resolve()
        if row.format == "labeled":
            frame = read_labeled(matrix_path)
        elif row.format == "plink":
            if not getattr(row, "id_file", ""):
                raise ValueError(f"{row.name}: PLINK format requires id_file")
            id_column = int(getattr(row, "id_column", "1") or 1)
            frame = read_plink(
                matrix_path, Path(row.id_file).expanduser().resolve(), id_column
            )
        else:
            raise ValueError(f"{row.name}: format must be labeled or plink")

        frame = frame.loc[frame.index, frame.index]
        values = frame.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        if row.kind == "similarity":
            if not np.isfinite(values).all() or values.min() < -1e-8 or values.max() > 1 + 1e-8:
                raise ValueError(f"{row.name}: similarity must lie in [0,1]")
            values = 1.0 - values
            np.fill_diagonal(values, 0)
            frame = pd.DataFrame(values, index=frame.index, columns=frame.columns)
        elif row.kind != "distance":
            raise ValueError(f"{row.name}: kind must be distance or similarity")
        matrices[row.name] = validate_distance(row.name, frame)
    return manifest, matrices


def upper(matrix: np.ndarray) -> np.ndarray:
    return matrix[np.triu_indices_from(matrix, k=1)]


def mantel(
    first: np.ndarray, second: np.ndarray, permutations: int, rng: np.random.Generator
) -> tuple[float, float]:
    index = np.triu_indices_from(first, k=1)
    x = first[index]
    observed = pearsonr(x, second[index]).statistic
    exceed = 0
    for _ in range(permutations):
        order = rng.permutation(len(first))
        permuted = second[np.ix_(order, order)][index]
        if abs(pearsonr(x, permuted).statistic) >= abs(observed):
            exceed += 1
    return float(observed), (exceed + 1) / (permutations + 1)


def pcoa(distance: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = len(distance)
    centering = np.eye(n) - np.ones((n, n)) / n
    gram = -0.5 * centering @ (distance ** 2) @ centering
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    positive = eigenvalues > max(1e-12, eigenvalues[0] * 1e-10)
    coordinates = eigenvectors[:, positive] * np.sqrt(eigenvalues[positive])
    return coordinates, eigenvalues, eigenvalues[positive]


def procrustes_test(
    first: np.ndarray,
    second: np.ndarray,
    dimensions: int,
    permutations: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, float, float, float, int]:
    dims = min(dimensions, first.shape[1], second.shape[1])
    if dims < 1:
        raise ValueError("No positive PCoA axes available for Procrustes analysis")
    aligned_first, aligned_second, disparity = procrustes(
        first[:, :dims], second[:, :dims]
    )
    statistic = np.sqrt(max(0.0, 1.0 - disparity))
    exceed = 0
    for _ in range(permutations):
        order = rng.permutation(len(second))
        _, _, permuted_disparity = procrustes(
            first[:, :dims], second[order, :dims]
        )
        permuted_statistic = np.sqrt(max(0.0, 1.0 - permuted_disparity))
        exceed += permuted_statistic >= statistic
    p_value = (exceed + 1) / (permutations + 1)
    return aligned_first, aligned_second, disparity, statistic, p_value, dims


def bh_adjust(values: pd.Series) -> np.ndarray:
    p = values.to_numpy(dtype=float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty_like(ranked)
    adjusted[order] = np.minimum(ranked, 1.0)
    return adjusted


def main() -> None:
    args = parse_args()
    if args.permutations < 1 or args.pcoa_dimensions < 1:
        raise ValueError("Permutations and PCoA dimensions must be positive")
    args.output_dir = args.output_dir.expanduser().resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest, matrices = load_manifest(args.manifest.expanduser().resolve())

    first_name = next(iter(matrices))
    common = [
        sample
        for sample in matrices[first_name].index
        if all(sample in matrix.index for matrix in matrices.values())
    ]
    if len(common) < args.minimum_samples:
        raise ValueError(f"Only {len(common)} common samples; minimum is {args.minimum_samples}")

    audit_rows = []
    aligned = {}
    coordinates = {}
    for name, matrix in matrices.items():
        excluded = [sample for sample in matrix.index if sample not in common]
        audit_rows.append(
            {
                "matrix": name,
                "n_input": len(matrix),
                "n_common": len(common),
                "n_excluded": len(excluded),
                "excluded_samples": ",".join(excluded),
            }
        )
        aligned[name] = matrix.loc[common, common]
        aligned[name].to_csv(
            args.output_dir / f"{name}.aligned.distance.tsv",
            sep="\t",
            index_label="sample",
        )
        coords, eigvals, positive = pcoa(aligned[name].to_numpy())
        coordinates[name] = coords
        pd.DataFrame(
            {
                "axis": np.arange(1, len(eigvals) + 1),
                "eigenvalue": eigvals,
                "is_positive": eigvals > 0,
            }
        ).to_csv(args.output_dir / f"{name}.pcoa_eigenvalues.tsv", sep="\t", index=False)

        sns.heatmap(aligned[name], cmap="viridis", square=True, xticklabels=False, yticklabels=False)
        plt.title(f"{name} distance | n={len(common)}")
        plt.tight_layout()
        plt.savefig(args.output_dir / f"{name}.distance_heatmap.png", dpi=300)
        plt.close()

    pd.DataFrame(audit_rows).to_csv(
        args.output_dir / "sample_alignment_audit.tsv", sep="\t", index=False
    )
    pd.DataFrame({"sample": common}).to_csv(
        args.output_dir / "common_samples.tsv", sep="\t", index=False
    )

    rng = np.random.default_rng(args.seed)
    results = []
    for first_name, second_name in itertools.combinations(aligned, 2):
        first = aligned[first_name].to_numpy()
        second = aligned[second_name].to_numpy()
        x, y = upper(first), upper(second)
        pearson_r = pearsonr(x, y).statistic
        spearman_r = spearmanr(x, y).statistic
        mantel_r, mantel_p = mantel(first, second, args.permutations, rng)
        a1, a2, disparity, proc_r, proc_p, dims = procrustes_test(
            coordinates[first_name], coordinates[second_name],
            args.pcoa_dimensions, args.permutations, rng
        )
        results.append({
            "matrix_1": first_name, "matrix_2": second_name,
            "n_samples": len(common), "n_distances": len(x),
            "pearson_r": pearson_r, "spearman_r": spearman_r,
            "mantel_r": mantel_r, "mantel_p": mantel_p,
            "procrustes_dimensions": dims, "procrustes_disparity": disparity,
            "procrustes_r": proc_r, "procrustes_p": proc_p,
        })

        fig, ax = plt.subplots(figsize=(6.5, 5.5))
        ax.scatter(x, y, s=18, alpha=0.55, color="#5C66A8")
        ax.set_xlabel(f"{first_name} distance"); ax.set_ylabel(f"{second_name} distance")
        ax.set_title(f"{first_name} vs {second_name}\nMantel r={mantel_r:.3f}, p={mantel_p:.4g}")
        fig.tight_layout(); fig.savefig(args.output_dir / f"{first_name}_vs_{second_name}.distance_scatter.png", dpi=300); plt.close(fig)

        fig, ax = plt.subplots(figsize=(7, 6))
        ax.scatter(a1[:, 0], a1[:, 1] if dims > 1 else np.zeros(len(a1)), label=first_name, color="#5C66A8")
        ax.scatter(a2[:, 0], a2[:, 1] if dims > 1 else np.zeros(len(a2)), label=second_name, color="#C97C7C")
        for i in range(len(common)):
            ax.plot([a1[i,0],a2[i,0]], [a1[i,1] if dims>1 else 0,a2[i,1] if dims>1 else 0], color="#999999", lw=.6, alpha=.6)
        ax.legend(frameon=False); ax.set_title(f"Procrustes: {first_name} vs {second_name}\nr={proc_r:.3f}, p={proc_p:.4g}")
        fig.tight_layout(); fig.savefig(args.output_dir / f"{first_name}_vs_{second_name}.procrustes.png", dpi=300); plt.close(fig)

    results = pd.DataFrame(results)
    results["mantel_q_bh"] = bh_adjust(results["mantel_p"])
    results["procrustes_q_bh"] = bh_adjust(results["procrustes_p"])
    results.to_csv(args.output_dir / "distance_concordance_summary.tsv", sep="\t", index=False)
    print(f"Compared {len(matrices)} matrices across {len(common)} common samples")
    print(f"Results: {args.output_dir}")


if __name__ == "__main__":
    main()
