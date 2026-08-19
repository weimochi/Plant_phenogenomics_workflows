import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import pdist

def parse_args():
    parser = argparse.ArgumentParser(description="RNA-seq Expression Matrix Preprocessing, QC, and Filtering.")
    parser.add_argument("--count_matrix", type=str, required=True, help="Path to raw count matrix (Gene x Sample).")
    parser.add_argument("--annotation", type=str, required=True, help="Path to gene annotation file (with gene length).")
    parser.add_argument("--meta", type=str, required=True, help="Path to sample metadata file.")
    parser.add_argument("--out_dir", type=str, default="processed_output", help="Output directory.")
    parser.add_argument("--min_library_size", type=int, default=1000000, help="Minimum total assigned counts per sample for library size QC.")
    parser.add_argument("--length_col", type=str, default="length", help="Annotation column containing gene/transcript length for TPM calculation.")
    parser.add_argument("--min_cpm", type=float, default=1.0, help="Min CPM threshold for gene filtering.")
    parser.add_argument("--min_samples_cpm", type=int, default=3, help="Min samples required to pass CPM threshold.")
    parser.add_argument("--mad_threshold", type=float, default=1.0, help="MAD threshold for high-variability gene filtering.")
    parser.add_argument("--hue_col", type=str, default=None, help="Column name in metadata for PCA coloring.")
    return parser.parse_args()

def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    
    print(">>> Step 1: Loading data...")
    counts = pd.read_csv(args.count_matrix, index_col=0)
    annot = pd.read_csv(args.annotation, index_col=0)
    meta = pd.read_csv(args.meta, index_col=0)
    
    # Record original dimensions
    original_n_genes = counts.shape[0]
    original_n_samples = counts.shape[1]
    print(f"Original count matrix: {original_n_genes:,} genes x {original_n_samples:,} samples")
    
    # =========================================================
    # Step 2: Sample Filtering (Library Size QC)
    # =========================================================
    print(">>> Step 2: Sample QC (Filtering by total assigned counts)...")
    sample_depth = counts.sum(axis=0)
    valid_samples = sample_depth[sample_depth >= args.min_library_size].index
    
    counts = counts.loc[:, valid_samples]
    
    # Safe metadata alignment with validation
    meta = meta.reindex(valid_samples)
    missing_meta_samples = meta.isna().all(axis=1).sum()
    if missing_meta_samples > 0:
        print(f"Warning: {missing_meta_samples} retained samples have no matching metadata entry (resulted in NaN).")
    
    # Keep original library sizes for retained samples (crucial for accurate CPM denominator)
    library_size = sample_depth.loc[valid_samples]
    
    after_qc_n_samples = counts.shape[1]
    removed_samples_step2 = original_n_samples - after_qc_n_samples
    print(f"Samples: {original_n_samples:,} -> {after_qc_n_samples:,} (removed {removed_samples_step2:,})")
    
    # =========================================================
    # Step 3: Annotation Matching & Gene-Length QC
    # =========================================================
    print(">>> Step 3: Annotation Matching & Gene-Length QC...")
    common_genes = counts.index.intersection(annot.index)
    counts = counts.loc[common_genes]
    annot = annot.loc[common_genes]
    annotation_n_genes = counts.shape[0]
    removed_genes_step3 = original_n_genes - annotation_n_genes
    
    if args.length_col not in annot.columns:
        raise ValueError(f"Length column '{args.length_col}' not found in annotation. Available columns: {annot.columns.tolist()}")
    gene_lens = annot[args.length_col]
    
    # Check valid gene lengths (remove missing, zero, or negative values)
    valid_length = gene_lens.notna() & (gene_lens > 0)
    n_invalid_length = (~valid_length).sum()
    if n_invalid_length > 0:
        print(f"Warning: removing {n_invalid_length:,} genes with missing or non-positive gene length.")
        counts = counts.loc[valid_length]
        annot = annot.loc[valid_length]
        gene_lens = gene_lens.loc[valid_length]
        
    valid_length_n_genes = counts.shape[0]
    removed_genes_step4 = annotation_n_genes - valid_length_n_genes
    
    print(f"  Original genes:             {original_n_genes:,}")
    print(f"  Annotation matched:         {annotation_n_genes:,}")
    print(f"  Valid gene length:          {valid_length_n_genes:,}")
    
    # =========================================================
    # Step 4: Normalization & Transformation
    # =========================================================
    print(">>> Step 4: Computing CPM, TPM, and log2(TPM+1)...")
    # CPM uses original library size of retained samples
    cpm = counts.div(library_size, axis=1) * 1e6
    
    # TPM is calculated from annotation-matched genes using specified length column
    rate = counts.div(gene_lens, axis=0)
    tpm = rate.div(rate.sum(axis=0), axis=1) * 1e6
    log_tpm = np.log2(tpm + 1)
    
    cpm.to_csv(os.path.join(args.out_dir, "matrix_cpm.csv"))
    tpm.to_csv(os.path.join(args.out_dir, "matrix_tpm.csv"))
    log_tpm.to_csv(os.path.join(args.out_dir, "matrix_log_tpm.csv"))
    
    # =========================================================
    # Step 5: Low-expression Gene Filtering
    # =========================================================
    print(">>> Step 5: Low-expression Gene Filtering...")
    pass_cpm = (cpm >= args.min_cpm).sum(axis=1) >= args.min_samples_cpm
    filtered_log_tpm = log_tpm.loc[pass_cpm]
    after_cpm_n_genes = filtered_log_tpm.shape[0]
    removed_genes_step5 = valid_length_n_genes - after_cpm_n_genes
    print(f"Genes after CPM filter: {valid_length_n_genes:,} -> {after_cpm_n_genes:,} (removed {removed_genes_step5:,})")
    
    # Robustness checks before downstream QC/PCA
    if filtered_log_tpm.shape[0] == 0:
        raise ValueError("No genes remain after CPM filtering. Consider lowering --min_cpm or --min_samples_cpm.")
    if filtered_log_tpm.shape[1] < 2:
        raise ValueError("Fewer than 2 samples remain after sample QC. PCA cannot be performed.")

    # =========================================================
    # Step 6: Expression QC plots (on filtered expression matrix)
    # =========================================================
    print(">>> Step 6: Generating Expression QC plots...")
    
    # PCA
    pca = PCA(n_components=2)
    pca_res = pca.fit_transform(filtered_log_tpm.T)
    pca_df = pd.DataFrame(pca_res, index=filtered_log_tpm.columns, columns=['PC1', 'PC2'])
    pca_df = pca_df.join(meta)
    
    # Determine hue column safely
    if args.hue_col and args.hue_col in meta.columns:
        hue_to_use = args.hue_col
    elif len(meta.columns) > 0:
        hue_to_use = meta.columns[0]
        print(f"Notice: --hue_col not specified or invalid. Falling back to first metadata column: '{hue_to_use}'")
    else:
        hue_to_use = None

    plt.figure(figsize=(6, 6))
    sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue=hue_to_use)
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    plt.title("PCA of Log2(TPM+1)")
    plt.savefig(os.path.join(args.out_dir, "qc_pca.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Sample Correlation Heatmap
    corr_matrix = filtered_log_tpm.corr(method='pearson')
    plt.figure(figsize=(8, 8))
    sns.heatmap(corr_matrix, cmap='coolwarm', vmin=None, vmax=1.0, annot=False)
    plt.title("Sample-to-Sample Pearson Correlation Heatmap")
    plt.savefig(os.path.join(args.out_dir, "qc_sample_correlation_heatmap.png"), dpi=300, bbox_inches='tight')
    plt.close()

    # Sample Clustering Dendrogram
    dist_matrix = pdist(filtered_log_tpm.T, metric='euclidean')
    linkage_matrix = sch.linkage(dist_matrix, method='average')

    plt.figure(figsize=(10, 5))
    sch.dendrogram(linkage_matrix, labels=filtered_log_tpm.columns, leaf_rotation=90)
    plt.title("Sample Clustering Dendrogram (Outlier Check)")
    plt.xlabel("Samples")
    plt.ylabel("Euclidean Distance")
    plt.savefig(os.path.join(args.out_dir, "qc_sample_dendrogram.png"), dpi=300, bbox_inches='tight')
    plt.close()
    
    # =========================================================
    # Step 7: High-Variability Gene Selection (Feature Selection)
    # =========================================================
    print(">>> Step 7: High-Variability Gene Selection (MAD filtering)...")
    gene_median = filtered_log_tpm.median(axis=1)
    mad = filtered_log_tpm.sub(gene_median, axis=0).abs().median(axis=1)
    variable_genes = mad >= args.mad_threshold
    final_log_tpm = filtered_log_tpm.loc[variable_genes]
    final_n_genes = final_log_tpm.shape[0]
    removed_genes_step7 = after_cpm_n_genes - final_n_genes
    
    if final_n_genes == 0:
        print("Warning: No genes passed the MAD threshold. Consider lowering --mad_threshold.")
        
    print(f"Genes after MAD filter: {after_cpm_n_genes:,} -> {final_n_genes:,} (removed {removed_genes_step7:,})")

    # =========================================================
    # Print Terminal Summary
    # =========================================================
    print("\n========== PREPROCESSING SUMMARY ==========")
    print(f"Samples: {original_n_samples:,} -> {after_qc_n_samples:,} (removed {removed_samples_step2:,})")
    print(f"Genes: {original_n_genes:,} -> {final_n_genes:,} (removed {original_n_genes - final_n_genes:,})")
    print(f"  Original genes:             {original_n_genes:,}")
    print(f"  Annotation matched:         {annotation_n_genes:,}")
    print(f"  Valid gene length:          {valid_length_n_genes:,}")
    print(f"  After CPM filtering:        {after_cpm_n_genes:,}")
    print(f"  After MAD filtering:        {final_n_genes:,}")
    print("===========================================\n")

    # =========================================================
    # Save Outputs with Clear, Self-Documenting Filenames
    # =========================================================
    filtered_log_tpm.to_csv(os.path.join(args.out_dir, "matrix_cpm_filtered_log_tpm.csv"))
    final_log_tpm.to_csv(os.path.join(args.out_dir, "matrix_high_variability_log_tpm.csv"))

    summary = pd.DataFrame({
        "step": [
            "Original",
            "Sample library size QC",
            "Annotation matching",
            "Valid gene length",
            "CPM filtering",
            "MAD filtering"
        ],
        "n_samples": [
            original_n_samples,
            after_qc_n_samples,
            after_qc_n_samples,
            after_qc_n_samples,
            after_qc_n_samples,
            after_qc_n_samples
        ],
        "n_genes": [
            original_n_genes,
            original_n_genes,
            annotation_n_genes,
            valid_length_n_genes,
            after_cpm_n_genes,
            final_n_genes
        ],
        "removed_samples": [
            0,
            removed_samples_step2,
            0,
            0,
            0,
            0
        ],
        "removed_genes": [
            0,
            0,
            removed_genes_step3,
            removed_genes_step4,
            removed_genes_step5,
            removed_genes_step7
        ]
    })
    summary.to_csv(os.path.join(args.out_dir, "preprocessing_summary.csv"), index=False)

    print(">>> Preprocessing & QC completed successfully!")

if __name__ == "__main__":
    main()