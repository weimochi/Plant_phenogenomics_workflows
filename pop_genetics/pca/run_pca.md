# PCA with PLINK

This document describes a reproducible, step-by-step workflow for principal component analysis (PCA) of a filtered DNA variant dataset with PLINK 1.9. Commands are intentionally presented separately so that every intermediate result can be inspected before continuing.

## What PCA represents

PCA summarizes major axes of genetic variation among samples. Samples with similar genotypes tend to occur near one another in PC space, while separation can reflect population structure, species differences, ancestry, relatedness, batch effects, or uneven missingness.

PLINK 1.9 `--pca` first constructs a variance-standardized genomic relationship matrix and then eigendecomposes that matrix. In simplified diploid notation, the genotype dosage at variant `j` for sample `i` is:

```text
g_ij = 0, 1, or 2 copies of the counted allele
```

Each variant is centered by its expected dosage and scaled by its expected variance:

```text
z_ij = (g_ij - 2p_j) / sqrt(2p_j(1-p_j))
```

where `p_j` is the counted-allele frequency at variant `j`. The relationship matrix is derived from these standardized genotype values. Its eigenvectors give the sample PC coordinates, and its eigenvalues describe the amount of variation represented by each PC.

PCA should normally be calculated from common, high-quality, approximately independent biallelic SNPs. It should not be run directly on an unfiltered cohort VCF.

## Required software

```bash
module load plink
```

Confirm the executable and version:

```bash
command -v plink
plink --version
```

This workflow is written for PLINK 1.9. Record the exact version in the analysis log.

## Set paths and parameters

Edit these values before running the remaining commands:

```bash
VCF="/path/to/final.filtered.vcf.gz"
OUT_DIR="/path/to/population_structure/pca/results"
PREFIX="brassica_cohort"

GENO_MAX=0.10
MIND_MAX=0.10
MAF_MIN=0.05
LD_WINDOW=50
LD_STEP=10
LD_R2=0.20
N_PCS=20
```

Create the output directory:

```bash
mkdir -p "${OUT_DIR}"
```

`GENO_MAX`, `MIND_MAX`, and `MAF_MIN` are example starting values rather than universal biological cutoffs. Use thresholds appropriate for the cohort, sequencing depth, ploidy, and downstream question.

## Step 1: Convert the filtered VCF to PLINK format

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

The options mean:

- `--double-id`: copies the VCF sample ID into both the FID and IID fields.
- `--allow-extra-chr`: permits Brassica chromosome or contig names that do not follow the human chromosome naming scheme.
- `--vcf-half-call missing`: converts half-called genotypes such as `0/.` and `./1` to missing instead of interpreting them as complete diploid genotypes.
- `--biallelic-only strict`: retains strictly biallelic variants.
- `--snps-only just-acgt`: retains SNPs whose alleles are single A/C/G/T bases.

Expected files:

```text
${PREFIX}.raw.bed
${PREFIX}.raw.bim
${PREFIX}.raw.fam
${PREFIX}.raw.log
```

Inspect the PLINK log before continuing:

```bash
less "${OUT_DIR}/${PREFIX}.raw.log"
```

Pay attention to sample count, variant count, half-call warnings, duplicate variant IDs, chromosome warnings, and whether any samples contain no valid genotype calls.

## Step 2: Measure missingness before filtering

```bash
plink \
    --bfile "${OUT_DIR}/${PREFIX}.raw" \
    --allow-extra-chr \
    --missing \
    --out "${OUT_DIR}/${PREFIX}.raw_missingness"
```

Outputs:

```text
${PREFIX}.raw_missingness.imiss
${PREFIX}.raw_missingness.lmiss
```

- `.imiss` reports missingness per sample.
- `.lmiss` reports missingness per variant.
- `F_MISS = 0` means no missing genotypes for that sample or variant.
- `F_MISS = 1` means all relevant genotype calls are missing.

Review the distributions before choosing final thresholds. A sample with unusually high missingness can distort PCA and may indicate poor sequencing, mapping, or genotype calling.

## Step 3: Apply sample and variant QC

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

Parameter interpretation:

- `--mind 0.10`: removes samples with more than 10% missing genotypes.
- `--geno 0.10`: removes variants with more than 10% missing genotypes.
- `--maf 0.05`: removes variants with minor allele frequency below 5%.

Because sample and variant missingness depend on the current dataset, strict projects may perform sample and variant filtering in separate passes and recheck missingness afterward.

Recalculate missingness in the retained dataset:

```bash
plink \
    --bfile "${OUT_DIR}/${PREFIX}.qc" \
    --allow-extra-chr \
    --missing \
    --out "${OUT_DIR}/${PREFIX}.qc_missingness"
```

Check the retained counts:

```bash
wc -l "${OUT_DIR}/${PREFIX}.qc.fam"
wc -l "${OUT_DIR}/${PREFIX}.qc.bim"
```

## Step 4: LD pruning

Closely linked variants contain redundant information and can cause a few genomic regions to dominate the PCs. Identify a relatively independent marker set:

```bash
plink \
    --bfile "${OUT_DIR}/${PREFIX}.qc" \
    --allow-extra-chr \
    --indep-pairwise "${LD_WINDOW}" "${LD_STEP}" "${LD_R2}" \
    --out "${OUT_DIR}/${PREFIX}.ld"
```

With `50 10 0.20`, PLINK examines windows of 50 variants, advances by 10 variants, and removes variants until no remaining pair in the window exceeds the specified `r²` threshold.

Extract the retained variants:

```bash
plink \
    --bfile "${OUT_DIR}/${PREFIX}.qc" \
    --allow-extra-chr \
    --extract "${OUT_DIR}/${PREFIX}.ld.prune.in" \
    --make-bed \
    --out "${OUT_DIR}/${PREFIX}.pca_input"
```

Confirm how many variants remain:

```bash
wc -l "${OUT_DIR}/${PREFIX}.ld.prune.in"
wc -l "${OUT_DIR}/${PREFIX}.ld.prune.out"
wc -l "${OUT_DIR}/${PREFIX}.pca_input.bim"
```

## Step 5: Run PCA

```bash
plink \
    --bfile "${OUT_DIR}/${PREFIX}.pca_input" \
    --allow-extra-chr \
    --pca "${N_PCS}" header tabs \
    --out "${OUT_DIR}/${PREFIX}.pca"
```

Primary outputs:

```text
${PREFIX}.pca.eigenvec
${PREFIX}.pca.eigenval
${PREFIX}.pca.log
```

- `.eigenvec` contains FID, IID, and sample coordinates for each requested PC.
- `.eigenval` contains the corresponding eigenvalues.
- `.log` records the number of samples and variants actually used and any exclusions or warnings.

## How missing genotypes are handled

The missing-data behavior needs to be distinguished from explicit data filtering:

1. `--mind` removes an entire sample only when its missing rate is above the chosen threshold.
2. `--geno` removes an entire variant only when its missing rate is above the chosen threshold.
3. Missing calls that remain after filtering are not changed in the `.bed` file by this workflow.
4. During construction of the variance-standardized relationship matrix, a remaining missing genotype contributes the variant mean. In centered form, this has a contribution of zero and is equivalent to mean-dosage imputation for that PCA calculation.

For example, if the counted-allele frequency is `p = 0.30`, the expected diploid dosage is:

```text
2p = 0.60
```

A missing call is therefore treated as the mean dosage for the PCA calculation, not as genotype `0` and not as homozygous reference.

Mean imputation prevents a small amount of random missingness from forcing sample deletion. It does not repair systematically poor samples or variants. This is why missingness filtering and inspection must occur before PCA.

Do not add PLINK's `--fill-missing-a2` here. That option replaces missing calls with homozygous A2 genotypes, which is different from PCA mean imputation and can introduce artificial structure.

## What does genotype zero mean?

The values `0`, `1`, and `2` are valid allele dosages:

```text
0 = zero copies of the counted allele
1 = one copy of the counted allele
2 = two copies of the counted allele
```

PLINK does not remove a sample merely because it has many genotype dosages equal to `0`. A zero dosage is not a missing genotype.

The following cases are different:

- A variant with minor allele frequency `0` is monomorphic. It has zero genotype variance and cannot contribute information to variance-standardized PCA. The explicit `--maf 0.05` filter removes it before PCA.
- A sample with no usable genotype data cannot be analyzed. Current PLINK 1.9 builds report an error instead of returning an all-zero PCA result for such a sample.
- A PC eigenvalue can be zero or numerically near zero when more components are requested than the effective rank of the genotype matrix supports. PLINK does not interpret this as a reason to remove samples. Those PCs simply contain no additional informative variation and should not be plotted or interpreted.

The maximum number of informative PCs cannot exceed the effective rank of the centered genotype matrix, which is bounded by both the number of samples minus one and the number of informative variants. In practice, inspect `.eigenval` and retain only PCs with meaningful positive eigenvalues.

## Variance explained

For the returned PCs, a commonly reported descriptive percentage is:

```text
variance explained by PCk = eigenvalue_k / sum(returned eigenvalues) × 100
```

If only the first 20 PCs were requested, the denominator is the sum of those 20 returned eigenvalues, not necessarily the total variance of every possible PC. Figure labels should state this limitation, or PCA should be rerun with enough PCs to represent the denominator intended for reporting.

## Recommended checks before plotting

```bash
head "${OUT_DIR}/${PREFIX}.pca.eigenvec"
head "${OUT_DIR}/${PREFIX}.pca.eigenval"
grep -Ei 'warning|error|people|variants|remaining' "${OUT_DIR}/${PREFIX}.pca.log"
```

Confirm that:

- sample IDs match the annotation table;
- no expected samples disappeared unexpectedly;
- the PCA input contains enough LD-pruned variants;
- early eigenvalues are positive and decrease sensibly;
- PCs are not primarily associated with missingness, sequencing batch, or coverage;
- AA, CC, or other biologically distinct datasets are combined only when that comparison is intentional.

## Suggested output organization

```text
pca/
├── run_pca.md
└── results/
    ├── cohort.raw.*
    ├── cohort.raw_missingness.*
    ├── cohort.qc.*
    ├── cohort.qc_missingness.*
    ├── cohort.ld.prune.in
    ├── cohort.ld.prune.out
    ├── cohort.pca_input.*
    ├── cohort.pca.eigenvec
    ├── cohort.pca.eigenval
    └── cohort.pca.log
```

Plotting programs can be added later as independent, reusable modules. They are intentionally not required by this command-line PCA workflow.

## References

- [PLINK 1.9 population stratification and PCA](https://www.cog-genomics.org/plink/1.9/strat)
- [PLINK 1.9 input filtering](https://www.cog-genomics.org/plink/1.9/filter)
- [PLINK 1.9 relationship matrix calculation](https://www.cog-genomics.org/plink/1.9/distance)
- [PLINK 1.9 standard data input](https://www.cog-genomics.org/plink/1.9/input)
