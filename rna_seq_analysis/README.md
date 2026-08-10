# RNA-seq Analysis

This directory contains a streamlined, scalable pipeline for processing RNA-seq data, specifically optimized for competitive mapping across complex genomes (e.g., *Brassica* species).

## Step 0: Data Preprocessing & Quality Control

A rigorous bioinformatics pipeline requires assessing data quality *before* and *after* preprocessing to ensure the reliability of downstream alignments.

### a. Initial Quality Control (Raw Reads)
First, evaluate the raw sequencing data using `FastQC` to identify adapter contamination, k-mer biases, and per-base sequence quality drops.

```bash
export INPUT_DIR="/path/to/raw_reads"
export OUTPUT_DIR="/path/to/fastqc_raw"

# Count raw fastq files
N=$(find "$INPUT_DIR" -maxdepth 1 -type f -name "*.fastq.gz" | wc -l)

# Submit the FastQC array job
sbatch --array=0-$((N-1)) \
  --export=ALL,INPUT_DIR,OUTPUT_DIR \
  00_fastqc.slurm
```

### b. Adapter Trimming & Quality Filtering
Based on the initial QC results, use `fastp` to trim adapters, filter low-quality bases, and optionally remove 5' biases (e.g., trimming 15bp to remove hexamer priming bias).

```bash
export INPUT_DIR="/path/to/raw_reads"
export OUTPUT_DIR="/path/to/clean_reads"

# [Optional] Override default Illumina adapters if using a different platform (e.g., BGI)
# export AD1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
# export AD2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

# [Optional] Override default fastp filtering options
# export FASTP_OPTS="-5 -3 -W 16 -M 20 -q 20 -u 40 -l 80 --trim_poly_g --trim_poly_x"

# Count the number of R1 files
N=$(find "$INPUT_DIR" -maxdepth 1 -type f -name "*_R1_*.fastq.gz" | wc -l)

# Submit the trimming job (default TRIM_FRONT=15)
sbatch --array=0-$((N-1)) \
  --export=ALL,INPUT_DIR,OUTPUT_DIR,TRIM_FRONT=15 \
  00a_fastp_trim.slurm
```

### c. Post-trimming Quality Control (Clean Reads)
Finally, verify the trimming efficacy. While `fastp` generates its own HTML report, running `FastQC` again on the clean reads provides a standardized metric for comparison.

```bash
export INPUT_DIR="/path/to/clean_reads"
export OUTPUT_DIR="/path/to/fastqc_clean"

# Count cleaned fastq files
N=$(find "$INPUT_DIR" -maxdepth 1 -type f -name "*.fq.gz" | wc -l)

# Submit the FastQC array job
sbatch --array=0-$((N-1)) \
  --export=ALL,INPUT_DIR,OUTPUT_DIR \
  00_fastqc.slurm
```
## Step 1: Sequence Alignment (Mapping)

This pipeline supports two state-of-the-art splice-aware aligners. Choose either GSNAP or STAR based on your experimental requirements (e.g., handling complex plant genomes vs. ultrafast processing).

### a. Splice-aware Mapping via GSNAP (`01a_align_gsnap.slurm`)
GSNAP is particularly well-suited for handling complex genomes with novel splicing events and homoeolog-specific alignments. 

**Advanced Parameter Tuning:**
This script is highly flexible. You can override the default alignment behavior by exporting these variables before submission depending on your target species and analysis goals:

*   **`NOVEL`** (default: `1`): Enable (`1`) or disable (`0`) novel splicing discovery. Keep this on if you are looking for unannotated transcripts, or turn it off to strictly rely on known splice sites for faster processing.
*   **`NPATHS`** (default: `10`): Maximum number of alignments to report per read. Decrease this (e.g., `1`) if you only want unique alignments, or increase it when dealing with highly repetitive sequences or recent polyploids (like *Brassica*).
*   **`ONLY_CONCORDANT`** (default: `1`): Restrict output to properly paired reads. 
*   **`MAX_INTRON_MID` / `MAX_INTRON_END`** (default: `20000` / `5000`): Maximum allowed intron lengths. You may want to adjust these limits based on the specific intron size distribution of your target plant genome.
*   **`PAIRMAX_RNA`** (default: `1000`): Maximum expected genomic distance between paired-end reads.

```bash
export READS="/path/to/clean_reads"
export BASE="/path/to/output_directory"
export DB_DIR="/path/to/gmap_db_directory"
export DB_NAME="genome_db_name"

# [Optional] Example: Override defaults for a strict, unique-mapping run without novel splicing
# export NOVEL=0
# export NPATHS=1

# Count the number of R1 files
N=$(find "$READS" -maxdepth 1 -type f -name "*_R1_*.fq.gz" | wc -l)

# Submit the GSNAP array job
sbatch --array=0-$((N-1)) \
  --export=ALL,READS,BASE,DB_DIR,DB_NAME \
  01a_align_gsnap.slurm
```

### b. Ultrafast Alignment via STAR (`01b_align_star.slurm`)
STAR provides ultrafast two-pass mapping to reference genomes, generating sorted BAM files ready for downstream quantification.

**Advanced Parameter Tuning:**
STAR is extremely fast but highly sensitive to parameter tuning, especially when dealing with complex or polyploid plant genomes. You can override the defaults by exporting these variables before submission:

*   **`TWOPASS_MODE`** (default: `Basic`): Enables two-pass mapping to discover novel splice junctions during the first pass and use them to inform the second pass. Essential for accurate novel isoform detection.
*   **`OUT_FILTER_MULTIMAP_NMAX`** (default: `10`): Maximum number of multiple alignments allowed for a read. Increase this if mapping to highly homologous sub-genomes or repetitive gene families; reads exceeding this limit are considered unmapped.
*   **`ALIGN_INTRON_MAX` / `ALIGN_MATES_GAP_MAX`** (default: `20000`): Maximum intron size and maximum genomic gap between paired-end reads. Plant introns are generally shorter than mammalian ones, so tuning this down (e.g., to `10000` or `5000`) can reduce false-positive alignments spanning across distant genomic regions.
*   **`STAR_EXTRA`** (default: empty): A flexible catch-all for any additional STAR parameters. For example, if you want STAR to calculate read counts per gene alongside mapping, you can pass `--quantMode GeneCounts`.

```bash
export READS="/path/to/clean_reads"
export BASE="/path/to/output_directory"
export GENOME_DIR="/path/to/star_index_dir"

# [Optional] Example: Tuning for highly homologous genomes and generating gene counts
# export OUT_FILTER_MULTIMAP_NMAX=20
# export ALIGN_INTRON_MAX=10000
# export ALIGN_MATES_GAP_MAX=10000
# export STAR_EXTRA="--quantMode GeneCounts"

# Count the number of R1 files
N=$(find "$READS" -maxdepth 1 -type f -name "*_R1_*.fq.gz" | wc -l)

# Submit the STAR array job
sbatch --array=0-$((N-1)) \
  --export=ALL,READS,BASE,GENOME_DIR \
  01b_align_star.slurm
```
## Step 2: Mapping Statistics & Quality Control

This module contains scripts to evaluate alignment quality, compute coverage, and resolve competitive mapping between sub-genomes.

### a. Competitive NM & Coverage Analysis (`02a_compare_nm.slurm`)
This script evaluates sample alignments across multiple reference genomes. It extracts mismatch rates (NM) and calculates full-contig weighted coverage metrics to determine which sub-genome (e.g., AA vs. CC) best represents the sample. 

**Usage:**
Pass the target genomes as a space-separated string containing `NAME:DIRECTORY` pairs. 

```bash
# Define your genomes and output directory
export GENOMES="AA:/path/to/01_alignment/AA/bam CC:/path/to/01_alignment/CC/bam"
export OUTDIR="/path/to/02_stats_and_qc/results"

# [Optional] Adjust mapping quality thresholds or parallelization
# export MAPQ=20
# export JOBS=4
# export THREADS=2
# export UNCERT_DELTA=0.50

# Submit the comparison job
sbatch --export=ALL,GENOMES,OUTDIR 02a_compare_nm.slurm
```

### b. Depth & Gene-Body Coverage via Mosdepth (`02b_calc_depth.slurm`)
To ensure that reads are appropriately enriched within target regions rather than distributed across intergenic background noise, this script utilizes `mosdepth` to calculate global mapping depth, $1\text{x}$ genome coverage, and length-weighted gene-body mean depth.

**Usage:**
```bash
export BAMDIR="/path/to/01_alignment/AA/bam"
export OUT="/path/to/output_directory/AA_depth.tsv"
export GENEBED="/path/to/reference/gene_body.bed"

# [Optional] Adjust mapping quality threshold (default: Q=10)
# export Q=20

# Submit the job
sbatch --export=ALL,BAMDIR,OUT,GENEBED 02b_calc_depth.slurm
```
### c. Final Metadata & Report Synthesis (`02c_merge_metadata.py`)

This script acts as the final reporting step in the pipeline. It automatically detects the latest alignment statistics or competitive NM comparison report, merges the metrics with your laboratory sample metadata (e.g., sample tracking IDs, cultivar names, and batch info), and outputs a unified, ready-to-use tabular report.

**Prerequisites:**
*   Python 3 with the `pandas` library installed.
*   A completed sample list file (TSV format) containing at least a sample identifier column (`sample_id` or `RNAseq`).

**Usage:**

You can configure the input paths via environment variables and execute the script directly:

```bash
# Set environment variables for your paths
export SAMPLE_LIST="/path/to/SampleList.txt"
export MAP_DIR="/path/to/02_stats_and_qc/results"
export META_OUT="/path/to/merged_metadata"

# Run the merger script
python3 02c_merge_metadata.py
```
Alternative (Direct File Specification):
If you want to target a specific summary TSV file directly instead of letting the script auto-detect the latest one, you can pass it as a command-line argument:
```
python3 02c_merge_metadata.py /path/to/results/custom_compare_report.tsv
```
