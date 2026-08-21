# DNA variant calling, joint calling, filtering, and QC

This directory contains a modular workflow for germline small-variant analysis
from analysis-ready BAM files. It runs DeepVariant per sample, performs cohort
joint calling with GLnexus, generates baseline cohort QC, provides an
interactive filtering protocol, and compares the cohort before and after
filtering.

The workflow is separated into automated and decision-dependent stages:

```text
MarkDuplicates BAM files
        |
        +--> 01 DeepVariant (one sample per Slurm array task)
                |
                +--> per-sample VCF + gVCF
                        |
                        +--> 02 GLnexus joint calling (one cohort job)
                                |
                                +--> unfiltered cohort VCF/BCF
                                        |
                                        +--> 03 baseline variant QC
                                        |
                                        +--> 04 interactive filtering protocol
                                                |
                                                +--> filtered cohort VCF
                                                        |
                                                        +--> 05 pre/post-filter QC
```

`01`, `02`, and `03` can be chained with Slurm dependencies. Step `04` is a
Markdown protocol rather than a fixed batch script because genotype and site
filters should be selected after examining the actual cohort. Once filtering
thresholds have been accepted, Step `05` checks what changed.

## Scripts and documents

| File | Purpose | Execution mode |
|---|---|---|
| `01_deepvariant.slurm` | Per-sample small-variant calling | GPU Slurm array; one BAM per task |
| `02_glnexus_joint_calling.slurm` | Joint-call all completed DeepVariant gVCFs | One CPU cohort job |
| `03_variant_qc.slurm` | Generate filtering-before baseline QC | One CPU cohort job |
| `04_filter_variants.md` | Select and test genotype/site filters interactively | Interactive Slurm session |
| `05_filtered_variant_qc.slurm` | Compare raw and filtered cohort VCFs | One CPU cohort job |

## Important biological scope

DeepVariant's standard germline models emit diploid genotype states and are
primarily trained and validated using human data. The workflow can technically
run on non-human references, but technical completion does not guarantee that
the genotype model is biologically appropriate.

Before running a Brassica cohort, document:

- species and accession;
- expected ploidy;
- reference assembly and subgenome representation;
- whether homoeologous regions can be distinguished reliably;
- sequencing platform and expected coverage;
- whether samples are inbred, outbred, related, or population structured.

Diploid species such as *B. rapa* and *B. oleracea* are closer to the standard
DeepVariant genotype model. Allotetraploid material such as *B. napus* requires
additional consideration. A diploidized analysis, subgenome-specific calling,
or another polyploid-aware strategy may be necessary depending on the research
question.

DeepVariant documentation:

- <https://github.com/google/deepvariant>
- <https://github.com/google/deepvariant/blob/r1.10/docs/deepvariant-quick-start.md>

GLnexus documentation:

- <https://github.com/dnanexus-rnd/GLnexus>

## Software requirements

The following host commands must be available in the batch-job environment:

- `samtools`
- `bcftools`
- `apptainer` or `singularity`
- Bash and standard Unix tools

DeepVariant and GLnexus run in containers by default. The scripts can use a
remote `docker://` URI, but production array jobs should use pre-pulled local
SIF files whenever possible. This prevents many tasks from downloading or
writing to a shared container cache simultaneously.

Example image preparation:

```bash
CONTAINER_DIR=/absolute/path/to/containers
mkdir -p "${CONTAINER_DIR}"

apptainer pull \
    "${CONTAINER_DIR}/deepvariant_1.10.0-gpu.sif" \
    docker://google/deepvariant:1.10.0-gpu

apptainer pull \
    "${CONTAINER_DIR}/glnexus_1.4.1.sif" \
    docker://ghcr.io/dnanexus-rnd/glnexus:v1.4.1
```

Image downloads require network access and are best performed once from an
appropriate cluster node, not independently from every array task.

## Input from the mapping workflow

The expected DeepVariant inputs are the coordinate-sorted, duplicate-marked,
and indexed BAM files produced by `dna_mapping_and_qc`:

```text
<MARKDUP_OUTPUT>/bam/
├── sample01.markdup.bam
├── sample01.markdup.bam.bai
├── sample02.markdup.bam
└── sample02.markdup.bam.bai
```

The BAM filenames must have unique sample prefixes. The default DeepVariant
input pattern is `*.markdup.bam`, and the sample name is derived by removing
the `.markdup.bam` suffix.

The reference used for variant calling must be the same assembly used for
alignment. It must have a samtools FASTA index:

```bash
REFERENCE=/absolute/path/to/reference.fa
samtools faidx "${REFERENCE}"
```

---

## Step 1: Per-sample variant calling with DeepVariant

### Rationale

DeepVariant converts aligned reads into candidate examples, evaluates those
examples with a neural-network model, and converts predictions into VCF/gVCF
records. Its `run_deepvariant` wrapper coordinates three main stages:

1. `make_examples`: reads BAM/CRAM alignments and builds candidate examples;
2. `call_variants`: applies the trained model;
3. `postprocess_variants`: creates final VCF and gVCF outputs.

In GPU mode, `call_variants` uses one GPU. `make_examples` and much of the
remaining workflow still use CPU resources, so `NUM_SHARDS` remains important.

The VCF contains called variant sites for a single sample. The gVCF additionally
contains reference-confidence blocks and is the input required for GLnexus
joint calling.

### Count BAM inputs

The script discovers BAMs directly and does not require a separate manifest.
The input directory must remain unchanged while the array is running.

```bash
INPUT_BAM_DIR=/absolute/path/to/markduplicates/bam

N_BAMS=$(
    find "${INPUT_BAM_DIR}" \
        -maxdepth 1 \
        -type f \
        -name '*.markdup.bam' \
        | wc -l
)

(( N_BAMS > 0 )) || echo "ERROR: no markdup BAM files found"
```

### Submit DeepVariant

```bash
DEEPVARIANT_OUT=/absolute/path/to/results/01_deepvariant
DEEPVARIANT_IMAGE=/absolute/path/to/containers/deepvariant_1.10.0-gpu.sif

sbatch \
    --array="0-$((N_BAMS - 1))" \
    --export="ALL,INPUT_BAM_DIR=${INPUT_BAM_DIR},REFERENCE=${REFERENCE},OUTPUT_DIR=${DEEPVARIANT_OUT},DEEPVARIANT_IMAGE=${DEEPVARIANT_IMAGE}" \
    01_deepvariant.slurm
```

The script currently requests:

```text
GPU partition
1 GPU
8 CPUs
48 GB memory
```

Clusters use different GPU resource syntax. Confirm whether the site expects
`--gres=gpu:1`, `--gpus=1`, or a GPU type-specific request. A GPU partition
alone may not allocate a GPU.

### Main configurable variables

- `INPUT_BAM_DIR`: required directory containing analysis-ready BAM files.
- `REFERENCE`: required FASTA used for alignment and calling.
- `OUTPUT_DIR`: required DeepVariant output root.
- `BAM_SUFFIX` (default: `.markdup.bam`): suffix removed to derive the sample
  name.
- `BAM_PATTERN` (default: `*${BAM_SUFFIX}`): pattern used to discover BAMs.
- `MODEL_TYPE` (default: `WGS`): allowed values include `WGS`, `WES`, `PACBIO`,
  `ONT_R104`, and `HYBRID_PACBIO_ILLUMINA`. Select the model matching the data,
  not merely the reference genome.
- `NUM_SHARDS` (default: allocated CPUs): parallelism for `make_examples`. It
  may not exceed `SLURM_CPUS_PER_TASK`.
- `USE_GPU` (default: `1`): adds container GPU support with `--nv`. Setting this
  to `0` selects the CPU image, although the Slurm resource directives should
  also be changed to avoid reserving an unused GPU.
- `DEEPVARIANT_VERSION` (default: `1.10.0`): container tag used when an explicit
  image is not supplied.
- `DEEPVARIANT_IMAGE`: local SIF path or a `docker://` URI.
- `CONTAINER_RUNTIME`: optional explicit `apptainer` or `singularity` command.
- `EXTRA_BINDS`: optional additional container bind specification.
- `REGIONS` (default: whole genome): region string or region file. This is
  useful for pilot runs but produces a region-restricted call set.
- `CUSTOM_MODEL`: optional customized DeepVariant model/checkpoint.
- `VCF_STATS_REPORT` (default: `true`): request DeepVariant's visual stats
  report.
- `MAKE_EXAMPLES_EXTRA_ARGS`, `CALL_VARIANTS_EXTRA_ARGS`, and
  `POSTPROCESS_VARIANTS_EXTRA_ARGS`: advanced wrapper arguments. Record and
  validate these carefully because they change calling semantics.
- `SAMTOOLS_BIN` / `BCFTOOLS_BIN`: host executable names or paths.

### Pilot run

Before running an entire cohort, test a small region and a small number of
samples:

```bash
sbatch \
    --array="0-1" \
    --export="ALL,INPUT_BAM_DIR=${INPUT_BAM_DIR},REFERENCE=${REFERENCE},OUTPUT_DIR=${DEEPVARIANT_OUT},DEEPVARIANT_IMAGE=${DEEPVARIANT_IMAGE},REGIONS=chrA01" \
    01_deepvariant.slurm
```

Do not mix region-restricted and whole-genome gVCFs in the same joint-calling
cohort.

### Expected output

```text
01_deepvariant/
├── samples/
│   ├── sample01/
│   │   ├── sample01.deepvariant.vcf.gz
│   │   ├── sample01.deepvariant.vcf.gz.tbi
│   │   ├── sample01.deepvariant.g.vcf.gz
│   │   ├── sample01.deepvariant.g.vcf.gz.tbi
│   │   └── stage_logs/
│   └── sample02/
└── logs/
    ├── sample01.deepvariant.log
    └── sample01.deepvariant.config.tsv
```

Each sample is first written to a staging directory. VCF, gVCF, and indexes are
validated with bcftools before the sample directory is renamed to its final
path. Complete outputs are skipped on rerun; incomplete final directories
cause an error rather than being overwritten.

---

## Step 2: Cohort joint calling with GLnexus

### Rationale

Single-sample gVCFs contain both variant candidates and reference-confidence
information. GLnexus combines all gVCFs, unifies candidate alleles, and
revises/genotypes sites across the cohort. This produces one multisample cohort
call set suitable for filtering and population-level analysis.

`02_glnexus_joint_calling.slurm` is one non-array job. It automatically finds:

```text
<DEEPVARIANT_OUTPUT_DIR>/samples/*/*.deepvariant.g.vcf.gz
```

Every gVCF and index is validated. Each input must contain exactly one sample,
and sample names must be unique.

### RocksDB scratch requirements

GLnexus uses RocksDB and performs substantial random I/O. The database should
be placed on sufficiently large node-local scratch or another storage system
approved for this workload.

The script uses:

```bash
GLNEXUS_DB_ROOT="${GLNEXUS_DB_ROOT:-${SLURM_TMPDIR:-}}"
```

If the cluster does not define `SLURM_TMPDIR`, `GLNEXUS_DB_ROOT` is required.
The script deliberately does not fall back to the output directory.

Example:

```bash
GLNEXUS_DB_ROOT=/absolute/path/to/large/local_scratch
df -h "${GLNEXUS_DB_ROOT}"
```

Required scratch depends on sample count, genome size, variant density, and
gVCF structure. Confirm both capacity and I/O suitability before submitting a
large plant cohort.

### Submit GLnexus

```bash
GLNEXUS_OUT=/absolute/path/to/results/02_glnexus_joint_calling
GLNEXUS_IMAGE=/absolute/path/to/containers/glnexus_1.4.1.sif
COHORT_NAME=cohort

sbatch \
    --export="ALL,DEEPVARIANT_OUTPUT_DIR=${DEEPVARIANT_OUT},OUTPUT_DIR=${GLNEXUS_OUT},COHORT_NAME=${COHORT_NAME},GLNEXUS_IMAGE=${GLNEXUS_IMAGE},GLNEXUS_DB_ROOT=${GLNEXUS_DB_ROOT}" \
    02_glnexus_joint_calling.slurm
```

If the cluster provides `SLURM_TMPDIR` automatically inside the job, omit
`GLNEXUS_DB_ROOT` from `--export` and allow the script to use it.

### Main configurable variables

- `DEEPVARIANT_OUTPUT_DIR`: required root produced by Step 1.
- `OUTPUT_DIR`: required GLnexus output root.
- `COHORT_NAME` (default: `cohort`): output directory and filename prefix.
- `GVCF_PATTERN` (default: `*.deepvariant.g.vcf.gz`): discovered input pattern.
- `MIN_SAMPLES` (default: `2`): minimum number of gVCFs required.
- `GLNEXUS_CONFIG` (default: `DeepVariantWGS`): preset used for cohort calling.
- `GLNEXUS_VERSION` (default: `1.4.1`): default image tag.
- `GLNEXUS_IMAGE`: local SIF path or `docker://` URI.
- `GLNEXUS_BIN_IN_CONTAINER` (default: `/usr/local/bin/glnexus_cli`): executable
  path inside the image.
- `GLNEXUS_DB_ROOT`: required scratch root when `SLURM_TMPDIR` is unavailable.
- `CONTAINER_RUNTIME` / `EXTRA_BINDS`: container runtime customization.
- `BCFTOOLS_BIN`: host bcftools executable.

### Choosing a GLnexus preset

- `DeepVariantWGS`: official starting preset for DeepVariant WGS gVCFs.
- `DeepVariantWES`: preset for DeepVariant WES gVCFs.
- `DeepVariant_unfiltered`: retains more information and applies fewer preset
  decisions, but may greatly increase low-quality/complex sites and downstream
  filtering burden.
- `DeepVariant`: another available preset accepted by the script.

For non-human or complex plant cohorts, compare presets on a pilot cohort
rather than assuming that `DeepVariant_unfiltered` is automatically superior.
Evaluate record counts, missingness, heterozygosity, depth, allele balance,
Ti/Tv, and repetitive/homoeologous regions.

### Expected output

```text
02_glnexus_joint_calling/
└── cohort/
    ├── cohort.gvcf.list
    ├── cohort.samples.txt
    ├── cohort.deepvariant.glnexus.bcf
    ├── cohort.deepvariant.glnexus.bcf.csi
    ├── cohort.deepvariant.glnexus.vcf.gz
    ├── cohort.deepvariant.glnexus.vcf.gz.tbi
    ├── cohort.deepvariant.glnexus.bcftools.stats.txt
    ├── cohort.glnexus.config.tsv
    └── cohort.glnexus.log
```

The output sample count is checked against the number of input gVCFs. The
RocksDB directory is removed after successful completion, and the cohort
staging directory is renamed only after BCF, VCF, indexes, and stats pass
validation.

---

## Step 3: Baseline cohort variant QC

### Rationale

Filtering thresholds should be chosen from the unfiltered cohort rather than
copied from another organism. `03_variant_qc.slurm` creates a baseline without
changing any variants.

It reports:

- overall records, SNPs, INDELs, multiallelic sites, and related bcftools
  summary metrics;
- per-sample SNP reference-homozygous, alternate-homozygous, and heterozygous
  genotype counts;
- per-sample SNP heterozygosity fraction;
- transitions, transversions, and Ti/Tv;
- INDEL and singleton counts;
- average genotype depth when `FORMAT/DP` is available;
- missing genotype count and percentage;
- FILTER-label distribution;
- record counts per contig.

The `PSC` ref/het/hom fields produced by `bcftools stats` refer to SNPs; the
script names those derived columns explicitly to avoid treating them as
all-variant genotype counts.

### Submit baseline QC

```bash
JOINT_VCF="${GLNEXUS_OUT}/${COHORT_NAME}/${COHORT_NAME}.deepvariant.glnexus.vcf.gz"
VARIANT_QC_OUT=/absolute/path/to/results/03_variant_qc

sbatch \
    --export="ALL,INPUT_VCF=${JOINT_VCF},OUTPUT_DIR=${VARIANT_QC_OUT},COHORT_NAME=${COHORT_NAME}" \
    03_variant_qc.slurm
```

### Main configurable variables

- `INPUT_VCF`: required indexed cohort VCF or BCF.
- `OUTPUT_DIR`: required QC output root.
- `COHORT_NAME` (default: `cohort`): output subdirectory and filename prefix.
- `BCFTOOLS_BIN`: executable name or path.
- `THREADS`: inherited from `SLURM_CPUS_PER_TASK`.

### Expected output

```text
03_variant_qc/
└── cohort/
    ├── cohort.bcftools.stats.txt
    ├── cohort.site_summary.tsv
    ├── cohort.sample_qc.tsv
    ├── cohort.filter_summary.tsv
    ├── cohort.contig_summary.tsv
    ├── cohort.samples.txt
    ├── cohort.variant_qc.config.tsv
    └── cohort.bcftools.log
```

Review sample-level outliers before filtering sites. A poor sample can inflate
missingness, distort allele-frequency filters, and remove otherwise usable
variants across the entire cohort.

---

## Step 4: Interactive variant filtering

`04_filter_variants.md` is the filtering protocol and should be followed one
section at a time in an interactive Slurm allocation.

The protocol covers:

1. input and annotation validation;
2. reference-aware normalization;
3. decomposition of multiallelic records;
4. separate SNP and INDEL branches;
5. genotype masking with `bcftools +setGT`;
6. recalculation of `AC`, `AN`, `AF`, `MAF`, `NS`, and `F_MISSING`;
7. site-level filtering;
8. immediate QC after each trial;
9. optional recombination of filtered SNPs and INDELs;
10. recording final thresholds and biological assumptions.

Open the protocol:

```bash
less 04_filter_variants.md
```

### Why filtering is not initially a fixed Slurm script

The appropriate values for `MIN_GQ`, `MIN_DP`, `MAX_DP`, `MAX_MISSING`,
`MIN_MAF`, and `MIN_QUAL` depend on the cohort. They may also differ between
SNPs and INDELs. Inbred lines, wild populations, related samples, polyploid
material, and collapsed homoeologous regions can produce very different
distributions.

The Markdown protocol creates separately indexed intermediate files so each
decision can be inspected and compared. Once thresholds and filter order are
accepted, the commands can be frozen into a production
`04_filter_variants.slurm` if repeated execution is required.

Important cautions:

- mask low-quality genotypes before calculating cohort missingness/allele
  frequency;
- recalculate cohort tags after genotype masking;
- do not apply Hardy-Weinberg filtering routinely without considering
  population structure, selfing, breeding design, and ploidy;
- retain the separate filtered SNP and INDEL files even if a combined VCF is
  also generated;
- never overwrite the unfiltered GLnexus cohort VCF.

---

## Step 5: QC after filtering

### Rationale

Filtering should improve data quality without unexpectedly changing sample
identity or eliminating biologically important variation. Step 5 calculates
the same statistics for raw and filtered VCFs with the same bcftools version,
then produces direct pre/post comparisons.

The script requires sample names and sample order to remain identical. It
compares:

- site metrics and retained percentages;
- called SNP genotype retention per sample;
- SNP heterozygosity fraction;
- Ti/Tv;
- average depth;
- missingness;
- filtered-VFC FILTER distribution;
- filtered record counts per contig.

### Submit filtered QC

```bash
FILTERED_VCF=/absolute/path/to/variant_filtering/final_trial/cohort.filtered.vcf.gz
FILTERED_QC_OUT=/absolute/path/to/results/05_filtered_variant_qc

sbatch \
    --export="ALL,RAW_VCF=${JOINT_VCF},FILTERED_VCF=${FILTERED_VCF},OUTPUT_DIR=${FILTERED_QC_OUT},COHORT_NAME=${COHORT_NAME}" \
    05_filtered_variant_qc.slurm
```

### Main configurable variables

- `RAW_VCF`: required indexed unfiltered cohort VCF/BCF.
- `FILTERED_VCF`: required indexed filtered cohort VCF/BCF.
- `OUTPUT_DIR`: required comparison output root.
- `COHORT_NAME` (default: `cohort`): output subdirectory and prefix.
- `BCFTOOLS_BIN`: executable name or path.
- `THREADS`: inherited from the Slurm allocation.

### Expected output

```text
05_filtered_variant_qc/
└── cohort/
    ├── cohort.raw.bcftools.stats.txt
    ├── cohort.filtered.bcftools.stats.txt
    ├── cohort.raw.site_summary.tsv
    ├── cohort.filtered.site_summary.tsv
    ├── cohort.raw.sample_qc.tsv
    ├── cohort.filtered.sample_qc.tsv
    ├── cohort.pre_post_site_comparison.tsv
    ├── cohort.pre_post_sample_comparison.tsv
    ├── cohort.filtered.filter_summary.tsv
    ├── cohort.filtered.contig_summary.tsv
    ├── cohort.samples.txt
    └── cohort.filtered_variant_qc.config.tsv
```

The two primary comparison tables are:

```text
cohort.pre_post_site_comparison.tsv
cohort.pre_post_sample_comparison.tsv
```

Investigate unexpectedly low retention, increased missingness, strongly
shifted Ti/Tv, loss of particular contigs, or sample-specific outliers before
declaring the filtered VCF final.

---

## Automated submission through baseline QC

Steps 1 through 3 can be chained. Filtering remains an intentional manual
checkpoint.

```bash
N_BAMS=$(
    find "${INPUT_BAM_DIR}" \
        -maxdepth 1 \
        -type f \
        -name '*.markdup.bam' \
        | wc -l
)

DEEPVARIANT_JOB=$(
    sbatch \
        --parsable \
        --array="0-$((N_BAMS - 1))" \
        --export="ALL,INPUT_BAM_DIR=${INPUT_BAM_DIR},REFERENCE=${REFERENCE},OUTPUT_DIR=${DEEPVARIANT_OUT},DEEPVARIANT_IMAGE=${DEEPVARIANT_IMAGE}" \
        01_deepvariant.slurm
)

GLNEXUS_JOB=$(
    sbatch \
        --parsable \
        --dependency="afterok:${DEEPVARIANT_JOB}" \
        --export="ALL,DEEPVARIANT_OUTPUT_DIR=${DEEPVARIANT_OUT},OUTPUT_DIR=${GLNEXUS_OUT},COHORT_NAME=${COHORT_NAME},GLNEXUS_IMAGE=${GLNEXUS_IMAGE},GLNEXUS_DB_ROOT=${GLNEXUS_DB_ROOT}" \
        02_glnexus_joint_calling.slurm
)

VARIANT_QC_JOB=$(
    sbatch \
        --parsable \
        --dependency="afterok:${GLNEXUS_JOB}" \
        --export="ALL,INPUT_VCF=${GLNEXUS_OUT}/${COHORT_NAME}/${COHORT_NAME}.deepvariant.glnexus.vcf.gz,OUTPUT_DIR=${VARIANT_QC_OUT},COHORT_NAME=${COHORT_NAME}" \
        03_variant_qc.slurm
)

printf 'deepvariant_job\t%s\n' "${DEEPVARIANT_JOB}"
printf 'glnexus_job\t%s\n' "${GLNEXUS_JOB}"
printf 'baseline_qc_job\t%s\n' "${VARIANT_QC_JOB}"
```

If `SLURM_TMPDIR` is created only inside jobs, remove
`,GLNEXUS_DB_ROOT=${GLNEXUS_DB_ROOT}` from the GLnexus export line and allow
the script to resolve `SLURM_TMPDIR` at runtime.

## Resume and failure behaviour

The automated scripts use staging directories and explicit validation:

- complete validated outputs are skipped;
- incomplete final directories cause an error instead of being overwritten;
- partial output directories are removed by cleanup traps;
- VCF/BCF headers and indexes are checked with bcftools;
- DeepVariant publishes a sample only after VCF/gVCF validation;
- GLnexus publishes the cohort only after BCF, VCF, indexes, stats, and sample
  counts pass validation;
- QC directories are published only after all required tables exist.

Inspect retained log files before retrying a failed job. Do not delete an
incomplete final directory until the failure and recovery plan are understood.

## Recommended final deliverables

Preserve at least:

- the unfiltered GLnexus cohort BCF/VCF and indexes;
- the final filtered SNP and INDEL VCFs and indexes;
- the optional combined filtered VCF and index;
- DeepVariant and GLnexus configuration logs;
- baseline and filtered-QC tables;
- the exact filtering thresholds and commands;
- tool/container versions;
- reference assembly metadata;
- sample inclusion/exclusion records;
- species, ploidy, and population-design notes.

These records are required to reproduce the call set and to interpret whether
observed genotype patterns are technical, biological, or a consequence of the
chosen filtering strategy.
