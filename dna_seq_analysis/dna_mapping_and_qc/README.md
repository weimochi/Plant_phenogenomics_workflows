# DNA-seq quality control, trimming, mapping, and mapping QC

This directory contains a modular Slurm workflow for paired-end DNA sequencing
data. It covers raw-read quality assessment, optional read trimming,
BWA-MEM2 alignment, duplicate marking, and per-sample/cohort-level mapping
quality control.

The workflow is designed for array jobs and uses explicit file lists or TSV
manifests to keep sample-to-file relationships reproducible. Each processing
stage validates its inputs and writes sample-specific outputs. Expensive
reference preparation and cross-sample summary generation are kept outside
array tasks to avoid concurrent writes.

## Workflow overview

```text
Raw paired-end FASTQ
        |
        +--> 00a FastQC on raw reads
        |
        +--> Reads acceptable -------------------------------+
        |                                                     |
        +--> Reads require trimming                           |
                |                                             |
                +--> 00b fastp                                |
                        |                                     |
                        +--> 00c FastQC on trimmed reads       |
                                                              |
Final selected reads (raw or trimmed) <-----------------------+
        |
        +--> 01 BWA-MEM2 mapping and coordinate sorting
                |
                +--> 01b GATK MarkDuplicates
                        |
                        +--> 01c samtools/mosdepth mapping QC
                                |
                                +--> 01d cohort QC summary
```

Trimming is optional. Samples that pass raw-read QC can proceed directly to
mapping, while samples requiring cleanup can be processed with fastp first.
The final mapping manifest may therefore contain a mixture of raw and trimmed
FASTQ paths.

## Scripts

| Script | Purpose | Execution mode |
|---|---|---|
| `00a_fastqc_raw.slurm` | FastQC on raw FASTQ files | Slurm array; one FASTQ per task |
| `00b_fastp_trim.slurm` | Adapter/quality trimming of paired reads | Slurm array; one sample per task |
| `00c_fastqc_trimmed.slurm` | FastQC on trimmed FASTQ files | Slurm array; one FASTQ per task |
| `01_bwa_mem2_mapping.slurm` | BWA-MEM2 alignment and coordinate sorting | Slurm array; one sample per task |
| `01b_gatk_markduplicates.slurm` | Mark duplicate reads and produce metrics | Slurm array; one BAM per task |
| `01c_mapping_qc.slurm` | Mapping, depth, and coverage QC | Slurm array; one BAM per task |
| `01d_mapping_qc_summary.slurm` | Merge per-sample QC metrics | One non-array job |

## Software requirements

The following commands must be available in the batch-job environment:

- `fastqc`
- `fastp`
- `bwa-mem2`
- `samtools`
- `gatk`
- `mosdepth`
- standard Unix tools including Bash, `awk`, `find`, `sort`, and `sed`

The scripts do not assume universal module names because HPC installations
vary. Load the required modules before submission or add the appropriate
cluster-specific `module load` commands to the scripts.

Example only:

```bash
module load fastqc
module load fastp
module load bwa-mem2
module load samtools
module load gatk
module load mosdepth
```

Record the exact versions used for a production analysis. The mapping,
MarkDuplicates, and mapping-QC scripts also write tool-version information to
their output directories.

## General input conventions

Use absolute paths wherever possible. Input lists and manifests must:

- contain no header unless explicitly stated otherwise;
- contain no blank lines or comments;
- use Unix line endings;
- contain unique sample identifiers;
- remain unchanged while the corresponding array job is running.

Sample identifiers should contain only letters, numbers, `.`, `_`, and `-`.
The same identifier must be used consistently across mapping, duplicate
marking, QC, and downstream variant calling.

Slurm stdout and stderr are written to the directory from which `sbatch` is
submitted because the scripts use `%x_%A_%a.out` or `%x_%j.out`. Submit jobs
from a dedicated run/log directory if these files should not appear in the
current working directory.

---

## Step 0a: Initial quality control of raw reads

### Rationale

FastQC provides a standardized overview of raw sequencing quality. Important
modules include per-base quality, adapter content, sequence duplication,
overrepresented sequences, GC distribution, and k-mer enrichment. These
results should guide whether trimming is necessary; FastQC warnings alone do
not automatically imply that reads must be trimmed.

`00a_fastqc_raw.slurm` processes one FASTQ file per array task. R1 and R2 are
therefore separate entries in the file list.

### Prepare the raw FASTQ list

```bash
RAW_READ_DIR=/absolute/path/to/raw_reads
RAW_FASTQ_LIST=/absolute/path/to/manifests/raw_fastq.list

find "${RAW_READ_DIR}" \
    -maxdepth 1 \
    -type f \
    \( -iname '*.fastq.gz' -o -iname '*.fq.gz' \) \
    | LC_ALL=C sort \
    > "${RAW_FASTQ_LIST}"

N_RAW_FASTQ=$(wc -l < "${RAW_FASTQ_LIST}")
(( N_RAW_FASTQ > 0 )) || echo "ERROR: no raw FASTQ files found"
```

Example file list:

```text
/data/raw/sample01_R1.fastq.gz
/data/raw/sample01_R2.fastq.gz
/data/raw/sample02_R1.fastq.gz
/data/raw/sample02_R2.fastq.gz
```

### Submit FastQC

```bash
RAW_FASTQC_OUT=/absolute/path/to/results/00a_fastqc_raw

sbatch \
    --array="0-$((N_RAW_FASTQ - 1))" \
    --export="ALL,FILELIST=${RAW_FASTQ_LIST},OUTPUT_DIR=${RAW_FASTQC_OUT}" \
    00a_fastqc_raw.slurm
```

### Main configurable variables

- `FILELIST`: required; one FASTQ path per line.
- `OUTPUT_DIR`: required; destination for FastQC HTML/ZIP files and logs.
- `FASTQC_BIN` (default: `fastqc`): executable name or absolute path.

The script uses one CPU because each array task processes only one FASTQ.
FastQC's thread option mainly controls concurrent files rather than making a
single file scale efficiently across many CPUs.

Expected output:

```text
00a_fastqc_raw/
├── fastqc/
│   ├── sample01_R1_fastqc.html
│   ├── sample01_R1_fastqc.zip
│   └── ...
├── logs/
└── tmp/
```

---

## Step 0b: Optional adapter and quality trimming

### Rationale

fastp can remove adapter sequence, trim low-quality read ends, remove a fixed
number of 5-prime bases, filter reads by length/quality, and produce HTML/JSON
reports. Trimming should be driven by library design and raw-read QC. Excessive
trimming can shorten reads, reduce mappability, and introduce coverage bias.

Only samples that require trimming need to be included in the fastp manifest.

### Prepare the paired-read manifest

The manifest has three tab-separated columns and no header:

```text
sample<TAB>R1<TAB>R2
```

Example:

```text
sample01	/data/raw/sample01_R1.fastq.gz	/data/raw/sample01_R2.fastq.gz
sample03	/data/raw/sample03_R1.fastq.gz	/data/raw/sample03_R2.fastq.gz
```

```bash
FASTP_MANIFEST=/absolute/path/to/manifests/fastp_samples.tsv
N_FASTP=$(wc -l < "${FASTP_MANIFEST}")
(( N_FASTP > 0 )) || echo "ERROR: fastp manifest is empty"
```

### Submit fastp

```bash
FASTP_OUT=/absolute/path/to/results/00b_fastp_trim

sbatch \
    --array="0-$((N_FASTP - 1))" \
    --export="ALL,MANIFEST=${FASTP_MANIFEST},OUTPUT_DIR=${FASTP_OUT}" \
    00b_fastp_trim.slurm
```

### Main configurable variables

- `TRIM_FRONT` (default: `0`): fixed number of bases removed from the 5-prime
  end of both reads. Use only when library design or QC supports fixed
  trimming.
- `QUALIFIED_QUALITY` (default: `20`): Phred score used to define a qualified
  base.
- `UNQUALIFIED_PERCENT` (default: `40`): maximum percentage of low-quality
  bases allowed in a retained read.
- `MIN_LENGTH` (default: `80`): minimum read length after trimming.
- `CUT_FRONT` / `CUT_TAIL` (default: `1` / `1`): enable sliding-window quality
  trimming at the 5-prime and 3-prime ends.
- `CUT_WINDOW_SIZE` (default: `16`): window size for end-quality trimming.
- `CUT_MEAN_QUALITY` (default: `20`): minimum mean quality within the trimming
  window.
- `TRIM_POLY_G` (default: `1`): remove poly-G tails. This is especially
  relevant to two-colour Illumina systems; verify that it is appropriate for
  the sequencing platform.
- `TRIM_POLY_X` (default: `0`): remove general homopolymer tails. This is off by
  default to avoid unnecessary trimming of biological homopolymers.
- `ADAPTER_R1` / `ADAPTER_R2` (default: empty): explicit paired-end adapter
  sequences. Both must be provided together. When both are empty, fastp uses
  paired-end adapter detection.
- `OUT_SUFFIX` (default: `trimmed`): label used in output filenames.

Example parameter override:

```bash
sbatch \
    --array="0-$((N_FASTP - 1))" \
    --export="ALL,MANIFEST=${FASTP_MANIFEST},OUTPUT_DIR=${FASTP_OUT},TRIM_FRONT=10,MIN_LENGTH=70,TRIM_POLY_G=0" \
    00b_fastp_trim.slurm
```

Expected output:

```text
00b_fastp_trim/
├── reads/
│   ├── sample01_R1.trimmed.fastq.gz
│   └── sample01_R2.trimmed.fastq.gz
├── reports/
│   ├── sample01.trimmed.fastp.html
│   └── sample01.trimmed.fastp.json
└── logs/
```

---

## Step 0c: Quality control of trimmed reads

### Rationale

fastp already provides before/after metrics, so this step is optional. Running
FastQC again is useful when a standardized direct comparison between raw and
trimmed reads is required or when FastQC results are aggregated with MultiQC.

### Prepare the trimmed FASTQ list

```bash
TRIMMED_FASTQ_LIST=/absolute/path/to/manifests/trimmed_fastq.list

find "${FASTP_OUT}/reads" \
    -maxdepth 1 \
    -type f \
    -name '*.trimmed.fastq.gz' \
    | LC_ALL=C sort \
    > "${TRIMMED_FASTQ_LIST}"

N_TRIMMED_FASTQ=$(wc -l < "${TRIMMED_FASTQ_LIST}")
```

### Submit post-trimming FastQC

```bash
TRIMMED_FASTQC_OUT=/absolute/path/to/results/00c_fastqc_trimmed

sbatch \
    --array="0-$((N_TRIMMED_FASTQ - 1))" \
    --export="ALL,FILELIST=${TRIMMED_FASTQ_LIST},OUTPUT_DIR=${TRIMMED_FASTQC_OUT}" \
    00c_fastqc_trimmed.slurm
```

The configurable variables and output structure are equivalent to Step 0a.

---

## Step 1: BWA-MEM2 alignment and coordinate sorting

### Rationale

BWA-MEM2 aligns paired-end DNA reads to a reference genome. Its SAM output is
piped directly into `samtools sort`, avoiding a large intermediate SAM file.
The resulting coordinate-sorted BAM can be indexed and processed by duplicate
marking, coverage tools, and variant callers.

The script adds a read group containing `ID`, `SM`, `LB`, and `PL`. Correct and
unique `SM` values are essential because downstream tools use them to identify
samples.

### Prepare the reference once

Reference indexes must be built before submitting the mapping array. Do not
allow multiple array tasks to build the same indexes concurrently.

```bash
REFERENCE=/absolute/path/to/reference.fa

bwa-mem2 index "${REFERENCE}"
samtools faidx "${REFERENCE}"
```

The mapping script expects the modern BWA-MEM2 index files:

```text
reference.fa.0123
reference.fa.amb
reference.fa.ann
reference.fa.bwt.2bit.64
reference.fa.pac
```

### Prepare the final mapping manifest

The mapping manifest has three tab-separated columns and no header:

```text
sample<TAB>R1<TAB>R2
```

It may mix untrimmed and trimmed inputs:

```text
sample01	/results/00b_fastp_trim/reads/sample01_R1.trimmed.fastq.gz	/results/00b_fastp_trim/reads/sample01_R2.trimmed.fastq.gz
sample02	/data/raw/sample02_R1.fastq.gz	/data/raw/sample02_R2.fastq.gz
sample03	/results/00b_fastp_trim/reads/sample03_R1.trimmed.fastq.gz	/results/00b_fastp_trim/reads/sample03_R2.trimmed.fastq.gz
```

Each biological sample must appear exactly once. If a sample contains multiple
sequencing lanes, either merge the lanes deliberately or use unique lane IDs
and perform an explicit BAM merge later. Do not allow lane files to overwrite
one another under the same sample name.

```bash
MAPPING_MANIFEST=/absolute/path/to/manifests/mapping_samples.tsv
N_SAMPLES=$(wc -l < "${MAPPING_MANIFEST}")
(( N_SAMPLES > 0 )) || echo "ERROR: mapping manifest is empty"
```

### Submit mapping

```bash
MAPPING_OUT=/absolute/path/to/results/01_bwa_mem2_mapping

sbatch \
    --array="0-$((N_SAMPLES - 1))" \
    --export="ALL,MANIFEST=${MAPPING_MANIFEST},REFERENCE=${REFERENCE},OUTPUT_DIR=${MAPPING_OUT}" \
    01_bwa_mem2_mapping.slurm
```

### Main configurable variables

- `BWA_THREADS` (default: `12`): BWA-MEM2 alignment threads.
- `SORT_THREADS` (default: `3`): additional samtools sorting/compression
  threads.
- `SORT_MEMORY` (default: `2G`): approximate memory per samtools sort thread,
  not total sort memory.
- `RG_PLATFORM` (default: `ILLUMINA`): read-group platform (`PL`).
- `RG_LIBRARY_PREFIX` (default: `lib`): prefix used for the read-group library
  value (`LB`).
- `RG_CENTER` (default: empty): optional sequencing centre (`CN`).
- `BWA_MEM2_BIN` / `SAMTOOLS_BIN`: executable names or absolute paths.

The BWA and samtools processes run concurrently in a pipe. Their combined
thread allocation must not exceed `SLURM_CPUS_PER_TASK` after allowing for the
samtools main thread.

Expected output:

```text
01_bwa_mem2_mapping/
├── bam/
│   ├── sample01.sorted.bam
│   └── sample01.sorted.bam.bai
└── logs/
    ├── sample01.bwa-mem2.log
    ├── sample01.samtools-sort.log
    └── sample01.tool-versions.log
```

The BAM is first written with a partial filename, checked with
`samtools quickcheck`, and only then renamed to its final name. If a valid BAM
exists but its BAI is missing, the script rebuilds only the index.

---

## Step 1b: Duplicate marking with GATK

### Rationale

PCR or optical duplicates can create multiple read pairs that appear to
support the same molecule. GATK MarkDuplicates identifies these records and
sets the duplicate SAM flag. By default, reads are retained rather than
physically removed, allowing downstream tools to decide how duplicates should
be handled.

### Prepare the MarkDuplicates manifest

The manifest has two tab-separated columns and no header:

```text
sample<TAB>coordinate_sorted_bam
```

It can be generated from the mapping manifest because output names are
deterministic:

```bash
MARKDUP_MANIFEST=/absolute/path/to/manifests/markduplicates_samples.tsv

awk -F '\t' -v bam_dir="${MAPPING_OUT}/bam" '
    BEGIN { OFS="\t" }
    { print $1, bam_dir "/" $1 ".sorted.bam" }
' "${MAPPING_MANIFEST}" > "${MARKDUP_MANIFEST}"
```

### Submit MarkDuplicates

```bash
MARKDUP_OUT=/absolute/path/to/results/01b_gatk_markduplicates

sbatch \
    --array="0-$((N_SAMPLES - 1))" \
    --export="ALL,MANIFEST=${MARKDUP_MANIFEST},OUTPUT_DIR=${MARKDUP_OUT}" \
    01b_gatk_markduplicates.slurm
```

### Main configurable variables

- `JAVA_HEAP` (default: `40g`): Java maximum heap. Keep this below Slurm
  `--mem` to leave room for JVM/native overhead.
- `REMOVE_DUPLICATES` (default: `false`): mark and retain duplicates when
  `false`; physically remove duplicate records when `true`. Retaining marked
  reads is the safer general default for germline variant calling.
- `OPTICAL_DUPLICATE_PIXEL_DISTANCE` (default: `100`): pixel-distance
  threshold used for optical duplicate detection. Appropriate values depend on
  the flow-cell/instrument generation.
- `GATK_BIN` / `SAMTOOLS_BIN`: executable names or paths.

Expected output:

```text
01b_gatk_markduplicates/
├── bam/
│   ├── sample01.markdup.bam
│   └── sample01.markdup.bam.bai
├── metrics/
│   └── sample01.duplicate_metrics.txt
└── logs/
```

The script requires a BAM declared as coordinate-sorted. It writes BAM and
metrics to partial paths, validates the BAM, publishes the outputs, and builds
the index. A missing index can be rebuilt without repeating duplicate marking.

---

## Step 1c: Per-sample mapping and coverage QC

### Rationale

This stage combines complementary mapping metrics:

- `samtools flagstat`: mapped, paired, properly paired, duplicate, and related
  flag-based counts;
- `samtools idxstats`: mapped/unmapped records per reference sequence;
- `samtools stats`: detailed alignment and insert-size statistics;
- `mosdepth`: mean depth and cumulative genome coverage.

The script derives coverage percentages such as 1x, 5x, 10x, 20x, and 30x
from mosdepth's cumulative global distribution. It does not incorrectly treat
columns in `mosdepth.summary.txt` as threshold coverage.

### Prepare the QC manifest

```bash
QC_MANIFEST=/absolute/path/to/manifests/mapping_qc_samples.tsv

awk -F '\t' -v bam_dir="${MARKDUP_OUT}/bam" '
    BEGIN { OFS="\t" }
    { print $1, bam_dir "/" $1 ".markdup.bam" }
' "${MAPPING_MANIFEST}" > "${QC_MANIFEST}"
```

### Submit mapping QC

```bash
MAPPING_QC_OUT=/absolute/path/to/results/01c_mapping_qc

sbatch \
    --array="0-$((N_SAMPLES - 1))" \
    --export="ALL,MANIFEST=${QC_MANIFEST},OUTPUT_DIR=${MAPPING_QC_OUT}" \
    01c_mapping_qc.slurm
```

### Main configurable variables

- `MOSDEPTH_EXCLUDE_FLAG` (default: `1796`): excludes reads carrying any of
  the unmapped, secondary, QC-fail, or duplicate flags. Change only when the
  intended coverage definition is clear.
- `MOSDEPTH_MAPQ` (default: `0`): minimum mapping quality included in coverage.
  Increasing this removes low-confidence alignments from depth estimates.
- `MOSDEPTH_FAST_MODE` (default: `1`): faster calculation that does not inspect
  internal CIGAR operations or correct mate overlaps. Disable for a more
  detailed coverage interpretation when overlapping paired reads matter.
- `MOSDEPTH_PRECISION` (default: `5`): decimal precision in mosdepth
  distribution output.
- `COVERAGE_LEVELS` (default: `1,5,10,20,30`): depths reported in the compact
  per-sample metrics table.
- `SAMTOOLS_BIN` / `MOSDEPTH_BIN`: executable names or paths.

Expected output:

```text
01c_mapping_qc/
└── samples/
    └── sample01/
        ├── sample01.flagstat.txt
        ├── sample01.idxstats.txt
        ├── sample01.samtools.stats.txt
        ├── sample01.mosdepth.mosdepth.summary.txt
        ├── sample01.mosdepth.mosdepth.global.dist.txt
        ├── sample01.mapping_qc_metrics.tsv
        ├── sample01.samtools.log
        ├── sample01.mosdepth.log
        └── sample01.tool_versions.txt
```

Each sample is written to a staging directory. The entire directory is renamed
to its final path only after all required outputs pass validation.

---

## Step 1d: Cohort mapping-QC summary

### Rationale

Array tasks must not append concurrently to a shared summary file. This final
non-array job runs only after all per-sample QC tasks finish. It validates that
all per-sample tables have identical headers, rejects duplicate sample IDs,
sorts rows by sample, and publishes one cohort TSV.

### Submit the summary job

```bash
QC_ROOT="${MAPPING_QC_OUT}/samples"
MAPPING_SUMMARY="${MAPPING_QC_OUT}/mapping_qc_summary.tsv"

sbatch \
    --export="ALL,QC_ROOT=${QC_ROOT},OUTPUT_TSV=${MAPPING_SUMMARY}" \
    01d_mapping_qc_summary.slurm
```

Expected summary columns:

```text
sample
mapped_pct
mean_depth
cov_1x_pct
cov_5x_pct
cov_10x_pct
cov_20x_pct
cov_30x_pct
```

`COVERAGE_LEVELS` must be kept consistent across all `01c` array tasks so that
their headers can be merged.

---

## Submitting the mapping stages with Slurm dependencies

Raw-read QC and trimming include a human decision point, so they are commonly
run and reviewed before the mapping chain is submitted. Once the final mapping
manifest has been prepared, Steps 1 through 1d can be chained automatically.

Prepare all three manifests before submission:

```bash
N_SAMPLES=$(wc -l < "${MAPPING_MANIFEST}")

awk -F '\t' -v bam_dir="${MAPPING_OUT}/bam" '
    BEGIN { OFS="\t" }
    { print $1, bam_dir "/" $1 ".sorted.bam" }
' "${MAPPING_MANIFEST}" > "${MARKDUP_MANIFEST}"

awk -F '\t' -v bam_dir="${MARKDUP_OUT}/bam" '
    BEGIN { OFS="\t" }
    { print $1, bam_dir "/" $1 ".markdup.bam" }
' "${MAPPING_MANIFEST}" > "${QC_MANIFEST}"
```

Submit the chain:

```bash
MAP_JOB=$(
    sbatch \
        --parsable \
        --array="0-$((N_SAMPLES - 1))" \
        --export="ALL,MANIFEST=${MAPPING_MANIFEST},REFERENCE=${REFERENCE},OUTPUT_DIR=${MAPPING_OUT}" \
        01_bwa_mem2_mapping.slurm
)

MARKDUP_JOB=$(
    sbatch \
        --parsable \
        --dependency="afterok:${MAP_JOB}" \
        --array="0-$((N_SAMPLES - 1))" \
        --export="ALL,MANIFEST=${MARKDUP_MANIFEST},OUTPUT_DIR=${MARKDUP_OUT}" \
        01b_gatk_markduplicates.slurm
)

QC_JOB=$(
    sbatch \
        --parsable \
        --dependency="afterok:${MARKDUP_JOB}" \
        --array="0-$((N_SAMPLES - 1))" \
        --export="ALL,MANIFEST=${QC_MANIFEST},OUTPUT_DIR=${MAPPING_QC_OUT}" \
        01c_mapping_qc.slurm
)

SUMMARY_JOB=$(
    sbatch \
        --parsable \
        --dependency="afterok:${QC_JOB}" \
        --export="ALL,QC_ROOT=${MAPPING_QC_OUT}/samples,OUTPUT_TSV=${MAPPING_QC_OUT}/mapping_qc_summary.tsv" \
        01d_mapping_qc_summary.slurm
)

printf 'mapping_job\t%s\n' "${MAP_JOB}"
printf 'markduplicates_job\t%s\n' "${MARKDUP_JOB}"
printf 'mapping_qc_job\t%s\n' "${QC_JOB}"
printf 'summary_job\t%s\n' "${SUMMARY_JOB}"
```

`afterok` ensures that a downstream stage is submitted only after the previous
job or array completes successfully.

## Resume and incomplete-output behaviour

The scripts distinguish complete outputs from partial or inconsistent outputs:

- complete validated outputs are skipped;
- a valid BAM without its BAI is re-indexed where supported;
- partial BAMs/staging directories are removed by cleanup traps;
- incomplete final output directories cause an error instead of silent
  overwriting;
- summary files are generated by one non-array job to avoid write races.

If a script reports an incomplete final output, inspect the corresponding log
and output directory before deciding whether to rename, archive, or remove the
incomplete result.

## Interpretation notes for complex plant genomes

- High duplication or multi-mapping may reflect real repeats, recent whole
  genome duplication, or homoeologous regions rather than library failure.
- A high mapped percentage does not guarantee uniquely or correctly mapped
  reads.
- Mean depth can hide strong differences among chromosomes, organellar
  contigs, repeats, and low-mappability regions; inspect per-contig outputs.
- High-depth regions may indicate collapsed repeats, organellar sequence, or
  homoeologous mapping.
- Duplicate percentage depends on library complexity, sequencing depth, insert
  size, and optical duplicate settings.
- Coverage and mapping-quality thresholds should be chosen for the biological
  question and reference quality, not copied blindly from human WGS workflows.

## Recommended handoff to variant calling

The analysis-ready inputs for downstream variant calling are:

```text
<MARKDUP_OUT>/bam/<sample>.markdup.bam
<MARKDUP_OUT>/bam/<sample>.markdup.bam.bai
```

Review `mapping_qc_summary.tsv` and sample-level QC before passing these BAMs to
DeepVariant or another caller. DeepVariant's standard germline models assume
diploid genotype states, so species and ploidy must be considered explicitly
for Brassica datasets.
