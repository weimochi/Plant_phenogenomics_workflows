# RNA-seq Variant Calling Pipeline

This directory contains a scalable, high-performance pipeline for calling variants (SNPs and Indels) from RNA-seq data. The workflow adheres to modified GATK best practices and leverages Google DeepVariant (GPU) for highly accurate variant calling in complex plant genomes (e.g., *Brassica*).

### SLURM Resource Adjustments (Optional)

The script comes with default SLURM resource allocations (`8 CPUs`, `32GB RAM`, `12 hours`). You can easily override these defaults directly from the command line without editing the script. 

The number of threads used by `samtools` dynamically adapts to the CPUs you allocate via `--cpus-per-task`.

*   **`--cpus-per-task`**: Increase if you want faster sorting and BAM compression (e.g., `--cpus-per-task=16`).
*   **`--mem`**: Increase if you encounter Out-Of-Memory (OOM) errors from GATK or SLURM (e.g., `--mem=64G`).
*   **`--time`**: Adjust the time limit if you have massive BAM files (e.g., `--time=24:00:00`).

*(Note: The GATK Java heap size is hardcoded to `-Xmx16g -Xms8g` inside the script. If you allocate significantly more memory via SLURM, you may also want to manually open the script and increase the `-Xmx` value for GATK).*

### Advanced Usage Example

If you are dealing with very large BAM files and want to allocate 16 CPUs and 64GB of memory per sample, you can submit the job like this:

```bash
export BAMDIR="/path/to/01_alignment/bam"
export OUTDIR="/path/to/02_variant_calling/01_markdup"

N=$(find "$BAMDIR" -maxdepth 1 -type f -name "*.bam" ! -name "*.md.bam" | wc -l)

# Override default SLURM resources on the fly
sbatch --array=0-$((N-1)) \
       --cpus-per-task=16 \
       --mem=64G \
       --export=ALL,BAMDIR,OUTDIR \
       01_markdup.slurm
```
## Pipeline Structure
1. `01_mark_duplicates/`: PCR duplicate marking and coordinate sorting.
2. `02_deepvariant/`: Variant calling per sample using DeepVariant (GPU).
3. `03_joint_calling/`: Cohort-level joint genotyping using GLnexus.
4. `04_variant_filtering/`: Quality control and variant filtering via bcftools.

---

## Step 1: Mark Duplicates (`01_mark_duplicates/01_markdup.slurm`)

This step processes aligned BAM files to identify and tag read pairs that are likely PCR duplicates. If the input BAM files are not coordinate-sorted, the script will automatically sort them before running GATK `MarkDuplicates`.

### Adjustable Parameters (Environment Variables)

You can customize the script's behavior by exporting the following variables before submission:

*   **`BAMDIR`** *(Required)*: The absolute path to the directory containing your input `.bam` files (usually the output from the alignment step).
*   **`OUTDIR`** *(Required)*: The directory where the deduplicated BAMs (`*.md.bam`) and metrics files will be saved.
*   **`TAG`** *(Optional, default: `md`)*: The suffix appended to the output files. For example, leaving it as default will produce `sample.md.bam`.

### Usage Example

This script uses SLURM Job Arrays to process multiple samples in parallel. First, count the number of input BAM files to determine the array size, then submit the job.

```bash
# 1. Define your input and output directories
export BAMDIR="/path/to/01_alignment/bam"
export OUTDIR="/path/to/02_variant_calling/01_markdup"

# 2. Count the number of samples (BAM files)
N=$(find "$BAMDIR" -maxdepth 1 -type f -name "*.bam" ! -name "*.md.bam" | wc -l)

# 3. Submit the array job
sbatch --array=0-$((N-1)) --export=ALL,BAMDIR,OUTDIR 01_markdup.slurm
```
## Step 2: DeepVariant (GPU) (`02_deepvariant/02_deepvariant_gpu.slurm`)

This step performs high-accuracy variant calling using Google DeepVariant. It leverages GPU acceleration and Singularity containers. To prevent network file system (NFS) I/O bottlenecks during execution, intermediate files are aggressively written to a configurable fast local storage path (e.g., a local SSD or NVMe drive).

### Adjustable Parameters (Environment Variables)

*   **`BAMDIR`** *(Required)*: Directory containing the `.md.bam` files from Step 1.
*   **`REF`** *(Required)*: Absolute path to the reference genome FASTA file.
*   **`OUTDIR`** *(Required)*: Output directory for the resulting `.vcf.gz` and `.g.vcf.gz` files.
*   **`TAG`** *(Optional)*: If provided, it will be appended to the output filenames (e.g., `sample.TAG.vcf.gz`). If left blank, filenames remain clean (`sample.vcf.gz`).
*   **`MODEL`** *(Optional, default: `WES`)*: DeepVariant model type (`WES` or `WGS`).
*   **`FAST_STORAGE`** *(Optional, default: `/tmp`)*: Path to a fast local scratch disk (SSD/NVMe). Setting this to a high-speed disk dramatically reduces compute time.
*   **`EXTRA_BINDS`** *(Optional)*: Additional Singularity bind mounts required to access your input data on the cluster (e.g., `-B /nas -B /data`).

### Usage Example (Array Job)

Since *Brassica* samples are often mapped competitively to multiple sub-genomes (e.g., AA and CC), you can submit separate array jobs for each reference without modifying the script.

```bash
export BAMDIR="/path/to/01_markdup"
export REF="/path/to/reference_genome.fa"
export OUTDIR_AA="/path/to/02_deepvariant/AA"
export TAG="group1"

# 1. Provide necessary cluster mount points (Replace with your actual data drives)
export EXTRA_BINDS="-B /path/to/nas -B /path/to/workdir"

# 2. Set this to your cluster's local scratch/SSD path for maximum performance
export FAST_STORAGE="/path/to/local_scratch" 

# 3. Count the number of MD BAMs
N=$(find "$BAMDIR" -maxdepth 1 -type f -name "*.md.bam" | wc -l)

# 4. Submit DeepVariant array job
sbatch --array=0-$((N-1)) \
       --export=ALL,BAMDIR,REF=$REF_AA,OUTDIR=$OUTDIR_AA,TAG,EXTRA_BINDS,FAST_STORAGE \
       02_deepvariant_gpu.slurm
```
## Step 3: Joint Calling (`03_joint_calling/03_glnexus.slurm`)

This step merges individual DeepVariant `.g.vcf.gz` files into a single, comprehensive cohort-level `.vcf.gz` matrix using GLnexus. Joint calling is critical for distinguishing true homozygous reference sites from sites with insufficient coverage across the population.

*Note: This is a single, heavily parallelized task (not an array job). It utilizes local SSD space for its internal database to prevent NFS I/O bottlenecks.*

### Adjustable Parameters (Environment Variables)

*   **`IN_DIR`** *(Required)*: Directory containing the `.g.vcf.gz` files from Step 2.
*   **`OUT_VCF`** *(Required)*: Absolute path and filename for the final joint VCF (e.g., `/path/to/AA_joint.vcf.gz`).
*   **`GLNEXUS_CONFIG`** *(Optional, default: `DeepVariantWES`)*: GLnexus configuration preset.
*   **`EXCLUDE_SAMPLES`** *(Optional)*: A comma-separated list of sample IDs to ignore (e.g., `sample1,sample2`).
*   **`EXCLUDE_FILE`** *(Optional)*: Path to a text file containing one sample ID per line to exclude.

### Usage Example

```bash
export IN_DIR="/path/to/02_deepvariant/AA"
export OUT_VCF="/path/to/03_joint/AA_joint.vcf.gz"
export EXCLUDE_SAMPLES="bad_sample_01,bad_sample_02"

sbatch --export=ALL,IN_DIR,OUT_VCF,EXCLUDE_SAMPLES 03_glnexus.slurm
```
## Step 4: Variant Filtering (`04_filtering/04_filter_vcf.slurm`)

This step rigorously filters the joint VCF file using `bcftools`. It performs masking of low-quality genotypes (based on Depth and Genotype Quality), recalculates missing rates and allele frequencies, and filters out monomorphic sites, high-missingness sites, and rare variants.

### Adjustable Parameters (Environment Variables)

*   **`IN_VCF`** *(Required)*: The input joint `.vcf.gz` file from Step 3.
*   **`OUT_VCF`** *(Required)*: The desired output filtered `.vcf.gz` file.
*   **`MIN_DP`** *(Optional, default: 10)*: Minimum Depth. Genotypes below this are set to missing (`./.`).
*   **`MIN_GQ`** *(Optional, default: 20)*: Minimum Genotype Quality. Genotypes below this are set to missing (`./.`).
*   **`MIN_MAF`** *(Optional, default: 0.05)*: Minimum Minor Allele Frequency.
*   **`MAX_MISSING`** *(Optional, default: 0.2)*: Maximum allowed missing rate (e.g., 0.2 means max 20% missing data allowed).
*   **`FAST_STORAGE`** *(Optional, default: `/tmp`)*: Fast local scratch disk for sorting temporary files.

---

### A. Production Mode (Automated via SLURM)

Once your parameters are decided, use the SLURM script to process large files quickly in a single piped stream without generating intermediate files.

```bash
export IN_VCF="/path/to/03_joint/cohort_joint.vcf.gz"
export OUT_VCF="/path/to/04_filtered/cohort_DP10_GQ20_MAF05.vcf.gz"
export MIN_DP="10"
export MIN_GQ="20"
export MIN_MAF="0.05"
export MAX_MISSING="0.2"

sbatch --export=ALL,IN_VCF,OUT_VCF,MIN_DP,MIN_GQ,MIN_MAF,MAX_MISSING \
       04_filter_vcf.slurm
```
### B. Parameter Tuning (Interactive Mode)
If you are exploring a new dataset and need to determine how many SNPs are lost at each filtering step, run these commands interactively in your terminal instead of using the SLURM script. This step-by-step method writes intermediate .bcf files so you can check the retained SNP counts.

```bash
INPUT="/path/to/raw_joint.vcf.gz"
OUTDIR="./filter_tuning_test"
mkdir -p "$OUTDIR"

# Check original SNP count
echo "Original SNPs: $(bcftools view -H -v snps "$INPUT" | wc -l)"

# Step 1: Keep only biallelic SNPs
bcftools view -m2 -M2 -v snps -Ob -o "${OUTDIR}/step1_snps.bcf" "$INPUT"
bcftools index -f "${OUTDIR}/step1_snps.bcf"
echo "After Step 1: $(bcftools view -H${OUTDIR}/step1_snps.bcf | wc -l)"

# Step 2: Mask low DP and low GQ as Missing (./.)
# (Does not remove rows, just alters genotypes)
bcftools filter -S . -e 'FMT/DP<10 | FMT/GQ<20' -Ob -o "${OUTDIR}/step2_setmiss.bcf" "${OUTDIR}/step1_snps.bcf"
bcftools index -f "${OUTDIR}/step2_setmiss.bcf"

# Step 3: Recalculate allele frequencies and missing rates
bcftools +fill-tags "${OUTDIR}/step2_setmiss.bcf" -Ob -o "${OUTDIR}/step3_tags.bcf" -- -t AC,AN,AF,MAF,F_MISSING
bcftools index -f "${OUTDIR}/step3_tags.bcf"

# Step 4: Strict filtering based on calculated tags (Adjust parameters here)
bcftools view -i 'AC>0 && F_MISSING<=0.2 && MAF>=0.05' -Oz -o "${OUTDIR}/final_filtered.vcf.gz" "${OUTDIR}/step3_tags.bcf"
tabix -f -p vcf "${OUTDIR}/final_filtered.vcf.gz"

# Check final retained SNPs
echo "Final SNPs: $(bcftools view -H${OUTDIR}/final_filtered.vcf.gz | wc -l)"
```