# Inbreeding coefficient analysis with PLINK

This document describes a reusable workflow for estimating individual inbreeding coefficients from autosomal diploid SNP genotypes with PLINK 1.9 and visualizing their distribution.

## Biological principle: Hardy–Weinberg equilibrium

For a biallelic locus with allele frequencies `p` and `q = 1 - p`, Hardy–Weinberg equilibrium (HWE) predicts genotype frequencies:

```text
AA = p²
Aa = 2pq
aa = q²
```

Therefore, the expected heterozygosity and homozygosity are:

```text
H_expected = 2pq
HOM_expected = p² + q² = 1 - 2pq
```

Inbreeding increases the probability that the two alleles carried by an individual are identical and generally produces a heterozygote deficit relative to HWE expectations. A conceptual form of the inbreeding coefficient is:

```text
F = 1 - H_observed / H_expected
```

Interpretation:

```text
F ≈ 0    observed heterozygosity is close to HWE expectation
F > 0    heterozygote deficit / excess homozygosity
F < 0    heterozygote excess relative to the estimated expectation
F → 1    very high homozygosity relative to expectation
```

A negative estimate is possible and does not automatically indicate a software error. It can reflect excess heterozygosity, sampling variation, hybrid ancestry, incorrect allele-frequency reference, genotype error, or biological processes that violate the simple HWE model.

## What PLINK `--het` calculates

PLINK counts observed and expected autosomal homozygous genotypes for each sample and reports a method-of-moments estimate:

```text
F = [O(HOM) - E(HOM)] / [N(NM) - E(HOM)]
```

where:

- `O(HOM)` is the observed number of homozygous genotype calls;
- `E(HOM)` is the expected homozygous count based on estimated or loaded allele frequencies;
- `N(NM)` is the number of nonmissing genotypes used for that sample;
- `F` is the method-of-moments inbreeding coefficient.

Since `N(NM) - O(HOM)` is the observed heterozygous count and `N(NM) - E(HOM)` is the expected heterozygous count, this is closely related to `1 - H_observed/H_expected`.

This is a genome-wide marker-based estimate. It is not a pedigree inbreeding coefficient and should not be interpreted as a direct probability without considering the marker set and allele-frequency reference.

## Important assumptions and limitations

The estimate depends on several assumptions:

- markers are autosomal and diploid;
- genotypes are reliably called;
- allele frequencies represent an appropriate reference population;
- missingness is not systematically related to genotype or population;
- enough informative polymorphic markers are available;
- strong LD does not cause a few genomic regions to receive excessive weight.

Population structure is especially important. If genetically distinct populations are pooled to estimate allele frequencies, the Wahlund effect can generate an apparent heterozygote deficit and inflate `F` even when individuals are not recently inbred. For strongly differentiated AA, CC, species, or ancestry groups, estimate allele frequencies and inbreeding within biologically meaningful groups, or explicitly use an appropriate external frequency file.

PLINK 1.9 `--het` is designed for autosomal diploid data. Results from true polyploids, mixed-ploidy datasets, aneuploid regions, or genotype calls that do not represent diploid dosages require a ploidy-aware method instead. Brassica datasets must therefore be checked carefully to confirm that the analyzed genotype representation is diploid for the selected samples and reference genome.

## Workflow overview

```text
Filtered cohort VCF
        │
        ▼
PLINK conversion and genotype QC
        │
        ▼
LD pruning
        │
        ▼
Define an appropriate allele-frequency population
        │
        ▼
PLINK --het
        │
        ▼
Per-sample F estimates and diagnostic plots
```

## Required software

```bash
module load plink
```

The plotting program requires Python 3 with NumPy, pandas, Matplotlib, and seaborn:

```bash
python -c 'import numpy, pandas, matplotlib, seaborn; print("Plotting dependencies are available")'
```

## Set paths and parameters

```bash
VCF="/path/to/final.filtered.vcf.gz"
OUT_DIR="/path/to/pop_genetics/inbreeding/results"
PREFIX="brassica_cohort"

GENO_MAX=0.10
MIND_MAX=0.10
MAF_MIN=0.05
LD_WINDOW=50
LD_STEP=10
LD_R2=0.20
```

```bash
mkdir -p "${OUT_DIR}"
```

The thresholds above are examples. Record the final values actually used.

## Step 1: Prepare a QC-filtered PLINK dataset

If a suitable PLINK dataset already exists, set `BFILE` to its prefix and continue to LD pruning.

Otherwise convert the VCF:

```bash
plink \
    --vcf "${VCF}" \
    --double-id \
    --allow-extra-chr \
    --vcf-half-call missing \
    --biallelic-only strict \
    --snps-only just-acgt \
    --make-bed \
    --out "${OUT_DIR}/${PREFIX}.raw"
```

Measure missingness:

```bash
plink \
    --bfile "${OUT_DIR}/${PREFIX}.raw" \
    --allow-extra-chr \
    --missing \
    --out "${OUT_DIR}/${PREFIX}.raw_missingness"
```

Apply sample and variant filters:

```bash
plink \
    --bfile "${OUT_DIR}/${PREFIX}.raw" \
    --allow-extra-chr \
    --mind "${MIND_MAX}" \
    --geno "${GENO_MAX}" \
    --maf "${MAF_MIN}" \
    --make-bed \
    --out "${OUT_DIR}/${PREFIX}.qc"
```

```bash
BFILE="${OUT_DIR}/${PREFIX}.qc"
```

## Step 2: LD pruning

PLINK's inbreeding estimators do not model LD. Use an approximately linkage-independent marker set so long linked regions do not dominate the estimate:

```bash
plink \
    --bfile "${BFILE}" \
    --allow-extra-chr \
    --indep-pairwise "${LD_WINDOW}" "${LD_STEP}" "${LD_R2}" \
    --out "${OUT_DIR}/${PREFIX}.ld"
```

```bash
plink \
    --bfile "${BFILE}" \
    --allow-extra-chr \
    --extract "${OUT_DIR}/${PREFIX}.ld.prune.in" \
    --make-bed \
    --out "${OUT_DIR}/${PREFIX}.inbreeding_input"
```

Check the retained sample and marker counts:

```bash
wc -l "${OUT_DIR}/${PREFIX}.inbreeding_input.fam"
wc -l "${OUT_DIR}/${PREFIX}.inbreeding_input.bim"
```

## Step 3A: Estimate F within the current cohort

When the current samples form an appropriate allele-frequency reference population:

```bash
plink \
    --bfile "${OUT_DIR}/${PREFIX}.inbreeding_input" \
    --allow-extra-chr \
    --het \
    --out "${OUT_DIR}/${PREFIX}.inbreeding"
```

Primary output:

```text
${PREFIX}.inbreeding.het
```

Expected columns:

```text
FID  IID  O(HOM)  E(HOM)  N(NM)  F
```

## Step 3B: Use explicitly defined allele frequencies

For small groups or projections onto a reference population, estimate frequencies in an appropriate reference set and reuse them with `--read-freq`.

The keep file contains FID and IID:

```text
family1 sample1
family2 sample2
```

Estimate reference frequencies:

```bash
REFERENCE_KEEP="/path/to/reference_samples.keep"

plink \
    --bfile "${OUT_DIR}/${PREFIX}.inbreeding_input" \
    --allow-extra-chr \
    --keep "${REFERENCE_KEEP}" \
    --freqx \
    --out "${OUT_DIR}/${PREFIX}.reference_frequency"
```

Calculate F using that frequency reference:

```bash
plink \
    --bfile "${OUT_DIR}/${PREFIX}.inbreeding_input" \
    --allow-extra-chr \
    --read-freq "${OUT_DIR}/${PREFIX}.reference_frequency.frqx" \
    --het \
    --out "${OUT_DIR}/${PREFIX}.inbreeding_reference_freq"
```

When the immediate dataset contains very few samples, an appropriate `--read-freq` file is important because frequency estimates from the small target set can be unstable. The reference frequency file must be derived from compatible variants and allele coding.

## Optional: small-sample correction

PLINK 1.9 omits the `n/(n-1)` multiplier in Nei's expected homozygosity formula by default. The `small-sample` modifier includes this correction while forcing allele frequencies to be estimated from founders in the current dataset:

```bash
plink \
    --bfile "${OUT_DIR}/${PREFIX}.inbreeding_input" \
    --allow-extra-chr \
    --het small-sample \
    --out "${OUT_DIR}/${PREFIX}.inbreeding_small_sample"
```

Do not add this modifier mechanically. Record whether it was used, and do not combine its interpretation with an external-frequency analysis as if they were identical estimators.

## Step 4: Plot the estimates

```bash
python plot_inbreeding.py \
    --het "${OUT_DIR}/${PREFIX}.inbreeding.het" \
    --output-prefix "${OUT_DIR}/${PREFIX}.inbreeding"
```

The plotting program exports:

```text
*.summary.tsv
*.samples.tsv
*.distribution.png
*.distribution.pdf
*.ranked.png
*.ranked.pdf
```

The distribution plot shows the overall F distribution and marks `F = 0`. The ranked plot orders samples from lowest to highest F so samples with unusually high homozygosity or unusually negative estimates can be inspected directly.

## Interpretation and diagnostic checks

Do not assign a universal cutoff such as `F > 0.1` without a biological and methodological justification. Instead, inspect the empirical distribution and compare it across relevant biological groups.

Check whether F is associated with:

- sample missingness;
- sequencing depth and genotype quality;
- species or ancestry group;
- sequencing batch;
- heterozygous genotype calling bias;
- reference-genome divergence;
- the selected MAF and LD-pruning thresholds.

High positive F may indicate inbreeding or selfing, but can also be inflated by population mixture, allele dropout, mapping bias, or under-called heterozygotes. Negative F may reflect hybrid ancestry or heterozygote excess, but can also arise from frequency misspecification or technical artifacts.

For selfing plant species, high homozygosity may be biologically expected. The result should still be interpreted relative to an appropriate population-specific allele-frequency baseline.

## `--het` versus `--hardy`

These commands answer different questions:

- `--het` summarizes observed versus expected homozygosity for each individual and reports individual F.
- `--hardy` tests Hardy–Weinberg equilibrium separately at each variant.

Running `--hardy` is useful as a site-level diagnostic, but it is not a replacement for individual inbreeding estimates:

```bash
plink \
    --bfile "${OUT_DIR}/${PREFIX}.inbreeding_input" \
    --allow-extra-chr \
    --hardy midp \
    --out "${OUT_DIR}/${PREFIX}.hwe_diagnostic"
```

HWE tests should generally be performed within populations rather than on a strongly structured pooled cohort.

## Suggested directory structure

```text
inbreeding/
├── run_inbreeding.md
├── plot_inbreeding.py
└── results/
    ├── cohort.raw.*
    ├── cohort.qc.*
    ├── cohort.ld.prune.in
    ├── cohort.inbreeding_input.*
    ├── cohort.inbreeding.het
    ├── cohort.inbreeding.summary.tsv
    └── cohort.inbreeding.*.png
```

## References

- [PLINK 1.9 inbreeding coefficients and Hardy–Weinberg statistics](https://www.cog-genomics.org/plink/1.9/basic_stats)
- [PLINK 1.9 LD pruning](https://www.cog-genomics.org/plink/1.9/ld)
- [PLINK 1.9 input filtering](https://www.cog-genomics.org/plink/1.9/filter)
