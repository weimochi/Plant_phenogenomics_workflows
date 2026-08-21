# 04 — Cohort variant filtering

This document is an interactive filtering protocol for the joint-called cohort
VCF produced by `02_glnexus_joint_calling.slurm`. It is intended to be run one
section at a time after reviewing the baseline results from
`03_variant_qc.slurm`.

Do not run large cohort operations on a login node. Request an interactive
Slurm allocation or convert the finalized commands into a batch script.

## Filtering principles

- Keep the original GLnexus VCF unchanged.
- Normalize variants before applying allele-based filters.
- Mask low-quality **genotypes** before calculating cohort missingness and
  allele frequency.
- Apply different site-level decisions to SNPs and INDELs when appropriate.
- Derive thresholds from the actual cohort instead of copying human defaults.
- Run the same QC before and after filtering.
- Record every threshold and tool version used for the final dataset.

DeepVariant and the standard GLnexus presets use diploid genotypes. Confirm
that this interpretation is appropriate for the Brassica material before
using genotype-based filters.

---

## 1. Start an interactive Slurm session

Adjust the resource syntax to match the cluster configuration.

```bash
srun \
    --partition=CPU \
    --cpus-per-task=8 \
    --mem=32G \
    --time=08:00:00 \
    --pty bash
```

Inside the allocated session:

```bash
set -u
set -o pipefail
export LC_ALL=C
```

Load the cluster's bcftools module if needed:

```bash
module load bcftools
```

Confirm the executable and version:

```bash
command -v bcftools
bcftools --version
```

---

## 2. Configure one filtering run

Use absolute paths. Create a new `RUN_ID` whenever thresholds change so that
an earlier result is not overwritten.

```bash
INPUT_VCF=/absolute/path/to/cohort.deepvariant.glnexus.vcf.gz
REFERENCE=/absolute/path/to/reference.fa
FILTER_ROOT=/absolute/path/to/variant_filtering

COHORT_NAME=cohort
THREADS="${SLURM_CPUS_PER_TASK:-8}"
RUN_ID="${RUN_ID:-trial_01}"
WORK_DIR="${FILTER_ROOT}/${RUN_ID}"

mkdir -p "${WORK_DIR}"
```

Validate the inputs:

```bash
[[ -s "${INPUT_VCF}" ]] || echo "ERROR: missing INPUT_VCF"
[[ -s "${INPUT_VCF}.tbi" || -s "${INPUT_VCF}.csi" ]] || echo "ERROR: missing VCF index"
[[ -s "${REFERENCE}" ]] || echo "ERROR: missing reference FASTA"
[[ -s "${REFERENCE}.fai" ]] || echo "ERROR: missing reference FASTA index"

bcftools view --header-only "${INPUT_VCF}" >/dev/null
bcftools index --stats "${INPUT_VCF}" | head
bcftools query -l "${INPUT_VCF}" | wc -l
```

Record the starting configuration:

```bash
{
    printf 'run_id\t%s\n' "${RUN_ID}"
    printf 'input_vcf\t%s\n' "${INPUT_VCF}"
    printf 'reference\t%s\n' "${REFERENCE}"
    printf 'threads\t%s\n' "${THREADS}"
    printf '\n'
    bcftools --version
} > "${WORK_DIR}/run_config.txt"
```

---

## 3. Inspect available VCF annotations

Filtering expressions must only use fields that are actually present.

```bash
bcftools view --header-only "${INPUT_VCF}" \
    > "${WORK_DIR}/input.header.txt"

grep '^##FORMAT=' "${WORK_DIR}/input.header.txt"
grep '^##INFO=' "${WORK_DIR}/input.header.txt"
grep '^##FILTER=' "${WORK_DIR}/input.header.txt"
```

Confirm that genotype fields such as `GT`, `DP`, `GQ`, and `AD` are present
before using them below.

Inspect the first few records without truncating the upstream process through
`head`:

```bash
bcftools view --max-alleles 10 --types snps,indels "${INPUT_VCF}" \
    | awk 'BEGIN { n=0 } /^#/ { next } n<5 { print; n++ }'
```

---

## 4. Normalize and split multiallelic records

`--check-ref e` stops on a REF mismatch. Do not use automatic REF swapping as
a substitute for diagnosing reference incompatibility.

```bash
NORMALIZED_VCF="${WORK_DIR}/${COHORT_NAME}.norm.split.vcf.gz"

bcftools norm \
    --threads "${THREADS}" \
    --fasta-ref "${REFERENCE}" \
    --check-ref e \
    --multiallelics -any \
    --output-type z \
    --output "${NORMALIZED_VCF}" \
    "${INPUT_VCF}" \
    2> "${WORK_DIR}/normalize.log"

bcftools index --threads "${THREADS}" --tbi "${NORMALIZED_VCF}" \
    || bcftools index --threads "${THREADS}" --csi "${NORMALIZED_VCF}"
```

Immediate checks:

```bash
bcftools view --header-only "${NORMALIZED_VCF}" >/dev/null
bcftools index --stats "${NORMALIZED_VCF}" | head
bcftools stats --threads "${THREADS}" "${NORMALIZED_VCF}" \
    > "${WORK_DIR}/${COHORT_NAME}.norm.split.stats.txt"
```

If normalization reports REF errors, stop and verify that DeepVariant,
GLnexus, and this step use the same reference assembly.

---

## 5. Split SNPs and INDELs

```bash
SNP_RAW="${WORK_DIR}/${COHORT_NAME}.snps.raw.vcf.gz"
INDEL_RAW="${WORK_DIR}/${COHORT_NAME}.indels.raw.vcf.gz"

bcftools view \
    --threads "${THREADS}" \
    --types snps \
    --min-alleles 2 \
    --max-alleles 2 \
    --output-type z \
    --output "${SNP_RAW}" \
    "${NORMALIZED_VCF}"

bcftools view \
    --threads "${THREADS}" \
    --types indels \
    --min-alleles 2 \
    --max-alleles 2 \
    --output-type z \
    --output "${INDEL_RAW}" \
    "${NORMALIZED_VCF}"

bcftools index --threads "${THREADS}" --tbi "${SNP_RAW}" \
    || bcftools index --threads "${THREADS}" --csi "${SNP_RAW}"

bcftools index --threads "${THREADS}" --tbi "${INDEL_RAW}" \
    || bcftools index --threads "${THREADS}" --csi "${INDEL_RAW}"
```

Check record counts:

```bash
printf 'SNP records\t'
bcftools index --stats "${SNP_RAW}" | awk '{ n += $3 } END { print n+0 }'

printf 'INDEL records\t'
bcftools index --stats "${INDEL_RAW}" | awk '{ n += $3 } END { print n+0 }'
```

This protocol intentionally retains only biallelic SNPs and INDELs at this
stage. If multiallelic sites are scientifically important, create a separate
analysis branch rather than silently discarding them.

---

## 6. Choose genotype-level thresholds

Do not copy these example values blindly. Determine them from sequencing
depth, `03_variant_qc.slurm`, library design, sample quality, and ploidy.

Example trial values:

```bash
MIN_GQ=20
MIN_DP=5
MAX_DP=60
```

Questions to answer before continuing:

- Does every sample have `FORMAT/DP` and `FORMAT/GQ`?
- Is `MAX_DP` appropriate for all samples, or should high-depth samples be
  reviewed separately?
- Are low-depth samples expected from the experimental design?
- Could high depth indicate collapsed repeats or homoeologous mapping?

Inspect per-sample values from the unfiltered baseline:

```bash
column -t -s $'\t' \
    /absolute/path/to/03_variant_qc/cohort/cohort.sample_qc.tsv \
    | less -S
```

---

## 7. Mask low-quality genotypes

This step changes failing genotypes to missing (`./.`) instead of deleting an
entire site because one sample has low DP or GQ.

Check that the `setGT` plugin is available:

```bash
bcftools plugin -l | grep '^setGT$'
```

Mask SNP genotypes:

```bash
SNP_GT_MASKED="${WORK_DIR}/${COHORT_NAME}.snps.gt_masked.vcf.gz"

bcftools +setGT \
    "${SNP_RAW}" \
    --threads "${THREADS}" \
    --output-type z \
    --output "${SNP_GT_MASKED}" \
    -- \
    --target-genotypes q \
    --new-gt . \
    --include "FMT/GQ<${MIN_GQ} || FMT/DP<${MIN_DP} || FMT/DP>${MAX_DP}"

bcftools index --threads "${THREADS}" --tbi "${SNP_GT_MASKED}" \
    || bcftools index --threads "${THREADS}" --csi "${SNP_GT_MASKED}"
```

Mask INDEL genotypes:

```bash
INDEL_GT_MASKED="${WORK_DIR}/${COHORT_NAME}.indels.gt_masked.vcf.gz"

bcftools +setGT \
    "${INDEL_RAW}" \
    --threads "${THREADS}" \
    --output-type z \
    --output "${INDEL_GT_MASKED}" \
    -- \
    --target-genotypes q \
    --new-gt . \
    --include "FMT/GQ<${MIN_GQ} || FMT/DP<${MIN_DP} || FMT/DP>${MAX_DP}"

bcftools index --threads "${THREADS}" --tbi "${INDEL_GT_MASKED}" \
    || bcftools index --threads "${THREADS}" --csi "${INDEL_GT_MASKED}"
```

If `DP` or `GQ` is missing from a genotype, inspect how bcftools evaluates the
expression before finalizing this command. Do not assume missing FORMAT values
will behave like numeric zero.

---

## 8. Recalculate cohort tags after genotype masking

Allele counts and missingness must be recalculated after changing genotypes.

```bash
SNP_TAGGED="${WORK_DIR}/${COHORT_NAME}.snps.gt_masked.tags.vcf.gz"
INDEL_TAGGED="${WORK_DIR}/${COHORT_NAME}.indels.gt_masked.tags.vcf.gz"

bcftools +fill-tags \
    "${SNP_GT_MASKED}" \
    --threads "${THREADS}" \
    --output-type z \
    --output "${SNP_TAGGED}" \
    -- \
    --tags AC,AN,AF,MAF,NS,F_MISSING

bcftools +fill-tags \
    "${INDEL_GT_MASKED}" \
    --threads "${THREADS}" \
    --output-type z \
    --output "${INDEL_TAGGED}" \
    -- \
    --tags AC,AN,AF,MAF,NS,F_MISSING

bcftools index --threads "${THREADS}" --tbi "${SNP_TAGGED}" \
    || bcftools index --threads "${THREADS}" --csi "${SNP_TAGGED}"

bcftools index --threads "${THREADS}" --tbi "${INDEL_TAGGED}" \
    || bcftools index --threads "${THREADS}" --csi "${INDEL_TAGGED}"
```

Inspect the recalculated tags:

```bash
bcftools query \
    -f '%CHROM\t%POS\t%AC\t%AN\t%AF\t%MAF\t%NS\t%F_MISSING\n' \
    "${SNP_TAGGED}" \
    | awk 'NR<=10 { print }'
```

---

## 9. Choose site-level thresholds

Example trial values only:

```bash
MAX_MISSING=0.20
MIN_MAF=0.01
MIN_QUAL=20
```

Before selecting final values, compare multiple trials. For highly inbred
lines, a very low heterozygosity rate may be expected. For structured wild
populations, aggressive MAF or Hardy–Weinberg filters may remove real biology.

Avoid applying a routine HWE filter before population structure, breeding
design, selfing rate, and ploidy are understood.

---

## 10. Apply site-level filters

The following example retains `PASS` or unlabelled (`.`) records and applies
the trial QUAL, missingness, and MAF thresholds.

Filter SNPs:

```bash
SNP_FILTERED="${WORK_DIR}/${COHORT_NAME}.snps.filtered.vcf.gz"

bcftools view \
    --threads "${THREADS}" \
    --include "(FILTER='PASS' || FILTER='.') && QUAL>=${MIN_QUAL} && F_MISSING<=${MAX_MISSING} && MAF>=${MIN_MAF}" \
    --output-type z \
    --output "${SNP_FILTERED}" \
    "${SNP_TAGGED}"

bcftools index --threads "${THREADS}" --tbi "${SNP_FILTERED}" \
    || bcftools index --threads "${THREADS}" --csi "${SNP_FILTERED}"
```

Filter INDELs:

```bash
INDEL_FILTERED="${WORK_DIR}/${COHORT_NAME}.indels.filtered.vcf.gz"

bcftools view \
    --threads "${THREADS}" \
    --include "(FILTER='PASS' || FILTER='.') && QUAL>=${MIN_QUAL} && F_MISSING<=${MAX_MISSING} && MAF>=${MIN_MAF}" \
    --output-type z \
    --output "${INDEL_FILTERED}" \
    "${INDEL_TAGGED}"

bcftools index --threads "${THREADS}" --tbi "${INDEL_FILTERED}" \
    || bcftools index --threads "${THREADS}" --csi "${INDEL_FILTERED}"
```

SNPs and INDELs do not have to use identical `MIN_QUAL`, missingness, or other
site-level thresholds. Use separate variables once the cohort distributions
have been reviewed.

---

## 11. Check each filtered branch immediately

```bash
bcftools stats \
    --threads "${THREADS}" \
    --samples - \
    "${SNP_FILTERED}" \
    > "${WORK_DIR}/${COHORT_NAME}.snps.filtered.stats.txt"

bcftools stats \
    --threads "${THREADS}" \
    --samples - \
    "${INDEL_FILTERED}" \
    > "${WORK_DIR}/${COHORT_NAME}.indels.filtered.stats.txt"
```

Quick comparison:

```bash
for stats_file in \
    "${WORK_DIR}/${COHORT_NAME}.norm.split.stats.txt" \
    "${WORK_DIR}/${COHORT_NAME}.snps.filtered.stats.txt" \
    "${WORK_DIR}/${COHORT_NAME}.indels.filtered.stats.txt"
do
    printf '\nFILE\t%s\n' "${stats_file}"
    awk -F '\t' '$1=="SN" { print $3 "\t" $4 }' "${stats_file}"
done
```

Review at least:

- number of records, SNPs, and INDELs;
- Ti/Tv;
- per-sample missingness;
- per-sample heterozygosity;
- mean depth;
- singleton count;
- samples that behave as outliers.

---

## 12. Optionally combine filtered SNPs and INDELs

Keep the separate files even if a combined file is needed downstream.

```bash
FINAL_VCF="${WORK_DIR}/${COHORT_NAME}.filtered.vcf.gz"
SORT_TMP="${SLURM_TMPDIR:-${WORK_DIR}/sort_tmp}"

mkdir -p "${SORT_TMP}"

bcftools concat \
    --allow-overlaps \
    --output-type u \
    "${SNP_FILTERED}" \
    "${INDEL_FILTERED}" \
    | bcftools sort \
        --temp-dir "${SORT_TMP}" \
        --output-type z \
        --output "${FINAL_VCF}"

bcftools index --threads "${THREADS}" --tbi "${FINAL_VCF}" \
    || bcftools index --threads "${THREADS}" --csi "${FINAL_VCF}"
```

Validate the combined file:

```bash
bcftools view --header-only "${FINAL_VCF}" >/dev/null
bcftools index --stats "${FINAL_VCF}" | head

bcftools stats \
    --threads "${THREADS}" \
    --samples - \
    "${FINAL_VCF}" \
    > "${WORK_DIR}/${COHORT_NAME}.filtered.stats.txt"
```

---

## 13. Record the chosen thresholds

```bash
{
    printf 'MIN_GQ\t%s\n' "${MIN_GQ}"
    printf 'MIN_DP\t%s\n' "${MIN_DP}"
    printf 'MAX_DP\t%s\n' "${MAX_DP}"
    printf 'MAX_MISSING\t%s\n' "${MAX_MISSING}"
    printf 'MIN_MAF\t%s\n' "${MIN_MAF}"
    printf 'MIN_QUAL\t%s\n' "${MIN_QUAL}"
    printf 'final_vcf\t%s\n' "${FINAL_VCF}"
} > "${WORK_DIR}/filter_thresholds.tsv"
```

Also record:

- Brassica species and ploidy;
- reference assembly version;
- sequencing platform and expected depth;
- whether samples are inbred, outbred, related, or population structured;
- excluded samples and the reason for exclusion;
- whether repetitive or low-mappability regions were excluded;
- the comparison used to choose `DeepVariantWGS` or
  `DeepVariant_unfiltered`.

---

## 14. When to convert this protocol into a Slurm script

Create `04_filter_variants.slurm` only after the thresholds and filtering
order have been accepted. The production script should then:

1. accept the finalized thresholds as explicit variables;
2. refuse to overwrite existing outputs;
3. preserve SNP and INDEL branches;
4. validate every output and index;
5. record all tool versions and parameters;
6. publish the final directory only after all stages succeed.

The next automated QC stage should run the equivalent of
`03_variant_qc.slurm` on the final filtered VCF and can be named
`05_filtered_variant_qc.slurm`.
