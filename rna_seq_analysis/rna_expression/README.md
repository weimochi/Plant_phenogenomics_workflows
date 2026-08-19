# RNA Expression Quantification & Matrix Integration

This module handles gene-level read quantification from sorted BAM files and integrates individual sample counts into a comprehensive raw expression matrix.

---

## Step 01a: Read Quantification (`01_gene_counts/01a_rna_htseq.slurm`)

This script processes sorted BAM files in parallel using Slurm array jobs to calculate gene-level read counts.

### Prerequisite: Inspecting Your GTF / GFF Annotation File

Before running `htseq-count`, **you must inspect your annotation file** to ensure the feature types and attribute tags match your data. Different genomes may use different conventions (e.g., `CDS` vs. `exon`, `gene_id` vs. `geneID` or `Parent`).

Run these quick inspection commands in your terminal to check the structure:

```bash
# 1. Check available feature types (column 3) in your GTF/GFF
zcat /path/to/annotation.gff3.gz | grep -v "^#" | awk '{print $3}' | sort | uniq -c

# 2. Check attribute formats (column 9) to see how gene IDs are named
zcat /path/to/annotation.gff3.gz | grep -v "^#" | head -n 5 | awk '{print $9}'

```

* *Common Feature Types (`--feature-type`)*: `CDS` (protein-coding regions) or `exon`.
* *Common ID Attributes (`--id-attr`)*: `gene_id`, `geneID`, `gene`, or `Parent`.

### Adjustable Parameters (Slurm `--export` Variables)

* **`BAM_DIR`** *(Required)*: Path to the directory containing sorted BAM files.
* **`BAM_LIST`** *(Required)*: Text file listing one BAM filename per line.
* **`GTF`** *(Required)*: Path to the GTF/GFF annotation file.
* **`OUT_DIR`** *(Required)*: Output directory for per-sample `.counts.txt` files.
* **`STRAND`** *(Optional, default: `no`)*: Library strandedness (`no`, `yes`, or `reverse`).
* **`FEATURE_TYPE`** *(Optional, default: `CDS`)*: Feature type to count from the GTF.
* **`ID_ATTR`** *(Optional, default: `gene_id`)*: Attribute tag in the GTF to use as the gene identifier.

### Usage Example

```bash
export BAM_DIR="/path/to/sorted_bams"
export BAM_LIST="/path/to/sorted_bams/bam_list.txt"
export GTF="/path/to/annotation.gff3.gz"
export OUT_DIR="/path/to/rna_expression/01_gene_counts/results"

# Calculate total number of samples from the list
N=$(wc -l < "$BAM_LIST")

# Submit array job
sbatch --array=0-$((N-1)) \
       --export=ALL,BAM_DIR,BAM_LIST,GTF,OUT_DIR,STRAND=no,FEATURE_TYPE=CDS,ID_ATTR=gene_id \
       01_gene_counts/01a_rna_htseq.slurm

```

---

## Step 01b: Expression Matrix Integration (`01_gene_counts/01b_make_expression_matrix.py`)

After all parallel tasks in Step 1a are successfully completed, this universal script merges all individual `.counts.txt` files into a single `Genes × Samples` raw count matrix using an outer join (missing values are filled with 0).

### Adjustable Command-Line Arguments

* **`--counts-dir`** *(Required)*: Directory containing the per-sample `*.counts.txt` files.
* **`--out`** *(Required)*: File path for the output merged TSV matrix.
* **`--keep-special-rows`** *(Optional)*: Flag to retain HTSeq special summary rows like `__no_feature` or `__ambiguous` (by default, these are filtered out).

### Usage Example

```bash
python 01_gene_counts/01b_make_expression_matrix.py \
    --counts-dir /path/to/rna_expression/01_gene_counts/results \
    --out /path/to/rna_expression/expression_matrix_all_samples.tsv
```
## Step 02: Expression Matrix Preprocessing, QC, and Feature Selection (`02_matrix_preprocessing`)

This module takes the raw gene-level count matrix and performs rigorous quality control, annotation matching, normalization, transformation, expression QC, and high-variability gene selection. It provides two separate script options depending on how you want to filter high-variability genes for downstream co-expression analyses (e.g., WGCNA or Pearson correlation networks).

---

## Complete Preprocessing Pipeline

```text
Raw count matrix (Genes × Samples)
        ↓
1. Sample library-size QC (Filter out low assigned count samples)
        ↓
2. Annotation matching & Gene-length QC (Ensure valid structural lengths)
        ↓
3. Normalization & Transformation (CPM, TPM, and log2(TPM + 1))
        ↓
4. Low-expression gene filtering (CPM ≥ 1 in ≥ 3 samples)
        ↓
5. Expression-profile QC (PCA, Correlation Heatmap, and Hierarchical Clustering)
        ↓
6. High-variability gene selection (Absolute MAD threshold OR Top-fraction ranking)
        ↓
Final Output Matrices (Preprocessed vs. Feature-Selected)

```

---

## Detailed Step-by-Step Rationale

* **Sample Library-Size QC**:
* **Why**: Removes unreliable samples with extremely low total assigned counts to prevent technical noise from dominating downstream correlations.
* **Note**: Uses pre-filtered library sizes as the denominator for CPM to avoid artificial inflation.


* **Annotation Matching & Gene-Length QC**:
* **Why**: Intersects gene IDs with reference annotation and drops genes with missing, zero, or negative lengths to prevent computational errors (NaN/inf) during TPM calculation.


* **Normalization & Transformation (log2(TPM + 1))**:
* TPM (Transcripts Per Million) corrects for both gene length and sequencing depth.
* Adding +1 provides mathematical safety against log(0) = −∞.
* log2 transformation compresses extreme outliers (e.g., housekeeping genes) and establishes linear symmetry for fold-change and correlation calculations.


* **Low-Expression Gene Filtering**:
* **Why**: Removes uninformative background noise (genes failing to meet CPM thresholds across a minimum number of samples).
* **Critical Order**: Performed before expression QC and PCA so that near-zero noise genes do not dilute biological structures.


* **Expression Profile QC**:
* Generates PCA scatter plots, sample-to-sample Pearson correlation heatmaps, and hierarchical clustering dendrograms to detect batch effects or outlier samples.


* **High-Variability Gene Selection**:
* **Option A (`02a`)**: Uses an **absolute MAD threshold** (`--mad_threshold`, e.g., MAD ≥ 1.0) to filter genes based on an explicit variability cutoff.
* **Option B (`02b`)**: Uses a **rank-based top-fraction filter** (`--mad_top_fraction`, e.g., top 25%) to retain a controllable, scale-invariant set of the most variable genes.



---

## Adjustable Command-Line Arguments

* **`--count_matrix`** (Required): Path to the raw gene count matrix.
* **`--annotation`** (Required): Path to the gene annotation file containing length information.
* **`--meta`** (Required): Path to the sample metadata file.
* **`--out_dir`** (Optional, default: `processed_output`): Directory to save output matrices, QC plots, and statistics.
* **`--min_library_size`** (Optional, default: `1000000`): Minimum total assigned counts per sample for library size QC.
* **`--length_col`** (Optional, default: `length`): Annotation column name for gene/transcript length.
* **`--min_cpm`** (Optional, default: `1.0`): Minimum CPM threshold for low-expression gene filtering.
* **`--min_samples_cpm`** (Optional, default: `3`): Minimum number of samples required to pass the CPM threshold.
* **`--mad_threshold`** (Specific to `02a`, default: `1.0`): Absolute MAD threshold for high-variability gene filtering.
* **`--mad_top_fraction`** (Specific to `02b`, default: `0.25`): Fraction of genes with the highest MAD to retain (e.g., 0.25 for top 25%).
* **`--hue_col`** (Optional): Metadata column name used for coloring PCA scatter plots.

---

## How to Explain the Feature Selection Choice to Supervisors / Reviewers

* **For Option A (`02a` - Absolute Threshold)**:
> "We applied an absolute median absolute deviation threshold after low-expression filtering to isolate genes with robust, predefined variability across biological samples."


* **For Option B (`02b` - Top Fraction)**:
> "We ranked genes by median absolute deviation after low-expression filtering and retained the top 25% most variable genes to reduce low-information features for downstream exploratory co-expression analyses."



---

## Usage Examples

### Option A: Absolute MAD Filtering (`02a`)

```bash
python 02a_preprocess_and_qc.py \
    --count_matrix ../01_gene_counts/expression_matrix_all_samples.tsv \
    --annotation /path/to/reference/gene_annotations.csv \
    --meta /path/to/metadata/sample_meta.csv \
    --out_dir ./processed_output_02a \
    --min_library_size 1000000 \
    --length_col length \
    --min_cpm 1.0 \
    --min_samples_cpm 3 \
    --mad_threshold 1.0 \
    --hue_col treatment

```

### Option B: Top-Fraction MAD Filtering (`02b`)

```bash
python 02b_preprocess_and_qc.py \
    --count_matrix ../01_gene_counts/expression_matrix_all_samples.tsv \
    --annotation /path/to/reference/gene_annotations.csv \
    --meta /path/to/metadata/sample_meta.csv \
    --out_dir ./processed_output_02b \
    --min_library_size 1000000 \
    --length_col length \
    --min_cpm 1.0 \
    --min_samples_cpm 3 \
    --mad_top_fraction 0.25 \
    --hue_col treatment

```