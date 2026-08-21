#!/usr/bin/env python3
"""Create population lists from an ADMIXTURE Q matrix."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--q-file", required=True, type=Path)
    p.add_argument("--fam-file", required=True, type=Path)
    p.add_argument("--k", required=True, type=int)
    p.add_argument("--threshold", type=float, default=0.8)
    p.add_argument("--minimum-size", type=int, default=2)
    p.add_argument("--fam-col", type=int, default=1)
    p.add_argument("--output-dir", required=True, type=Path)
    args = p.parse_args()
    if args.k < 1 or not 0 <= args.threshold <= 1 or args.minimum_size < 2:
        raise ValueError("Require K >= 1, threshold in [0,1], and minimum size >= 2")
    q = pd.read_csv(args.q_file, sep=r"\s+", header=None)
    fam = pd.read_csv(args.fam_file, sep=r"\s+", header=None, dtype=str)
    if q.shape != (len(fam), args.k):
        raise ValueError(f"Q shape {q.shape} does not match {len(fam)} samples and K={args.k}")
    if args.fam_col >= fam.shape[1]:
        raise ValueError("--fam-col is outside the FAM file")
    ids = fam.iloc[:, args.fam_col].astype(str)
    if ids.duplicated().any():
        raise ValueError("Sample IDs are not unique")
    q = q.apply(pd.to_numeric, errors="raise")
    if not np.allclose(q.sum(axis=1), 1, atol=1e-3):
        raise ValueError("Q rows do not sum to 1")
    major = q.to_numpy().argmax(axis=1) + 1
    prop = q.max(axis=1).to_numpy()
    assignment = q.copy()
    assignment.columns = [f"cluster_{i}" for i in range(1, args.k + 1)]
    assignment.insert(0, "sample", ids.to_numpy())
    assignment["major_cluster"] = major
    assignment["major_proportion"] = prop
    assignment["status"] = np.where(prop >= args.threshold, "assigned", "admixed")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    assignment.to_csv(args.output_dir / "admixture_assignments.tsv", sep="\t", index=False)
    for cluster in range(1, args.k + 1):
        members = assignment.loc[(assignment.major_cluster == cluster) & (assignment.status == "assigned"), "sample"]
        if len(members) < args.minimum_size:
            raise ValueError(f"cluster_{cluster} has {len(members)} samples; minimum is {args.minimum_size}")
        members.to_csv(args.output_dir / f"cluster_{cluster}.samples.txt", index=False, header=False)
    assignment.loc[assignment.status == "admixed", "sample"].to_csv(
        args.output_dir / "admixed.samples.txt", index=False, header=False
    )
    print(assignment.groupby(["status", "major_cluster"]).size().to_string())

if __name__ == "__main__":
    main()
