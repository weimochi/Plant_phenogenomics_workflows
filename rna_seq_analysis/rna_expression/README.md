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
