# ADMIXTURE analysis

This document describes a reusable workflow for model-based ancestry analysis with PLINK 1.9 and ADMIXTURE. Commands are separated so that every intermediate output can be inspected before continuing.

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
ADMIXTURE-ready PLINK BED/BIM/FAM
        │
        ▼
ADMIXTURE for user-defined K range
        │
        ├── ancestry proportions (*.K.Q)
        ├── cluster allele frequencies (*.K.P)
        ├── cross-validation logs
        └── one ancestry plot per K
```

ADMIXTURE performs the ancestry estimation. PLINK does not run ADMIXTURE itself; it prepares the biallelic, quality-filtered, LD-pruned binary genotype files that ADMIXTURE reads.

## Required software

```bash
module load plink
module load admixture
```

```bash
plink --version
admixture --help | head
```

Record the exact versions in the project log.

The plotting program requires Python 3 with NumPy, pandas, and Matplotlib:

```bash
python -c 'import numpy, pandas, matplotlib; print("Plotting dependencies are available")'
```

## Set paths and adjustable parameters

```bash
VCF="/path/to/final.filtered.vcf.gz"
OUT_DIR="/path/to/pop_genetics/admixture/results"
PREFIX="brassica_cohort"

GENO_MAX=0.10
MIND_MAX=0.10
MAF_MIN=0.05
LD_WINDOW=50
LD_STEP=10
LD_R2=0.20

K_MIN=1
K_MAX=10
CV_FOLDS=5
THREADS=8
SEED=20260821
```

The tested number of ancestral clusters is fully configurable. For example:

```bash
K_MIN=2
K_MAX=6
```

will test K=2, 3, 4, 5, and 6. K must be a positive integer, `K_MIN` must not exceed `K_MAX`, and K must remain smaller than the number of usable samples.

Create the output directory:

```bash
mkdir -p "${OUT_DIR}"
```

## Step 1: Convert VCF to PLINK binary format

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

Inspect the conversion log:

```bash
less "${OUT_DIR}/${PREFIX}.raw.log"
```

Expected files are `${PREFIX}.raw.bed`, `${PREFIX}.raw.bim`, and `${PREFIX}.raw.fam`.

## Step 2: Inspect missingness

```bash
plink \
    --bfile "${OUT_DIR}/${PREFIX}.raw" \
    --allow-extra-chr \
    --missing \
    --out "${OUT_DIR}/${PREFIX}.raw_missingness"
```

Review both files before choosing final thresholds:

```text
*.imiss  missingness per sample
*.lmiss  missingness per variant
```

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

- `--mind 0.10` removes samples with more than 10% missing genotypes.
- `--geno 0.10` removes variants with more than 10% missing genotypes.
- `--maf 0.05` removes rare and monomorphic variants below 5% MAF.

These are starting values, not universal thresholds. Population composition and sample size affect MAF, while sequencing depth affects missingness.

## Step 4: LD pruning

ADMIXTURE assumes markers provide approximately independent ancestry information. Dense linked variants can overweight long haplotypes or genomic regions, so use an LD-pruned marker set.

```bash
plink \
    --bfile "${OUT_DIR}/${PREFIX}.qc" \
    --allow-extra-chr \
    --indep-pairwise "${LD_WINDOW}" "${LD_STEP}" "${LD_R2}" \
    --out "${OUT_DIR}/${PREFIX}.ld"
```

```bash
plink \
    --bfile "${OUT_DIR}/${PREFIX}.qc" \
    --allow-extra-chr \
    --extract "${OUT_DIR}/${PREFIX}.ld.prune.in" \
    --make-bed \
    --out "${OUT_DIR}/${PREFIX}.admixture_input"
```

Inspect retained sample and marker counts:

```bash
wc -l "${OUT_DIR}/${PREFIX}.admixture_input.fam"
wc -l "${OUT_DIR}/${PREFIX}.admixture_input.bim"
wc -l "${OUT_DIR}/${PREFIX}.ld.prune.in"
```

## Step 5: Make Brassica chromosome codes readable by ADMIXTURE

ADMIXTURE reads PLINK BED/BIM/FAM files, but many versions reject non-human chromosome strings such as Brassica chromosome or assembly-contig names. Preserve the scientifically meaningful PLINK files and create a separate ADMIXTURE copy whose BIM chromosome column is set to `0`:

```bash
ADMIX_PREFIX="${OUT_DIR}/${PREFIX}.admixture"

cp "${OUT_DIR}/${PREFIX}.admixture_input.bed" "${ADMIX_PREFIX}.bed"
cp "${OUT_DIR}/${PREFIX}.admixture_input.fam" "${ADMIX_PREFIX}.fam"

awk 'BEGIN {OFS="\t"} {$1=0; print}' \
    "${OUT_DIR}/${PREFIX}.admixture_input.bim" \
    > "${ADMIX_PREFIX}.bim.partial"

mv "${ADMIX_PREFIX}.bim.partial" "${ADMIX_PREFIX}.bim"
```

This changes only the chromosome-label field in the ADMIXTURE-specific BIM copy. It does not change genotypes, alleles, variant order, or sample order. Perform LD pruning before this conversion because chromosome boundaries are needed for sensible LD calculations.

Confirm that the three files have matching prefixes:

```bash
ls -lh "${ADMIX_PREFIX}.bed" "${ADMIX_PREFIX}.bim" "${ADMIX_PREFIX}.fam"
head "${ADMIX_PREFIX}.bim"
```

## Step 6: Run a configurable range of K values

Move into the output directory because ADMIXTURE writes `.Q` and `.P` files into the current working directory:

```bash
cd "${OUT_DIR}"
ADMIX_BASENAME="${PREFIX}.admixture"
```

Validate the requested range:

```bash
[[ "${K_MIN}" =~ ^[1-9][0-9]*$ ]] || { echo "K_MIN must be a positive integer" >&2; exit 2; }
[[ "${K_MAX}" =~ ^[1-9][0-9]*$ ]] || { echo "K_MAX must be a positive integer" >&2; exit 2; }
(( K_MIN <= K_MAX )) || { echo "K_MIN must be <= K_MAX" >&2; exit 2; }

N_SAMPLES=$(wc -l < "${ADMIX_BASENAME}.fam")
(( K_MAX < N_SAMPLES )) || { echo "K_MAX must be smaller than sample count (${N_SAMPLES})" >&2; exit 2; }
```

Run ADMIXTURE with cross-validation:

```bash
for K in $(seq "${K_MIN}" "${K_MAX}"); do
    echo "[$(date)] Running ADMIXTURE K=${K}"

    admixture \
        --cv="${CV_FOLDS}" \
        -s "${SEED}" \
        "${ADMIX_BASENAME}.bed" \
        "${K}" \
        -j"${THREADS}" \
        2>&1 | tee "${ADMIX_BASENAME}.K${K}.log"

    [[ -s "${ADMIX_BASENAME}.${K}.Q" ]] || { echo "Missing Q output for K=${K}" >&2; exit 3; }
    [[ -s "${ADMIX_BASENAME}.${K}.P" ]] || { echo "Missing P output for K=${K}" >&2; exit 3; }
done
```

Important adjustable parameters:

- `K_MIN` and `K_MAX`: the range of ancestral-cluster models to test.
- `CV_FOLDS`: number of cross-validation folds; ADMIXTURE uses 5-fold CV by default.
- `THREADS`: CPU threads passed with `-j`.
- `SEED`: random seed used to make the run reproducible.

K=1 is a useful baseline but does not show population subdivision. If the aim is only to compare structured models, begin with K=2.

For publication-grade inference, consider multiple independent seeds per K because ADMIXTURE optimization may converge to different local solutions. The single-seed loop above is a clear baseline workflow; replicate management and cluster-label alignment should be added before interpreting instability across runs.

## Step 7: Extract CV errors

```bash
grep -h 'CV error' "${ADMIX_BASENAME}".K*.log \
    | tee "${ADMIX_BASENAME}.cv_errors.txt"
```

A lower CV error indicates better prediction of masked genotypes, but the minimum should not be treated as an automatic biological truth. Also consider:

- whether the improvement over neighboring K values is substantial;
- stability across random seeds;
- consistency with PCA, IBS, geography, taxonomy, and sampling design;
- whether small clusters are supported by enough individuals;
- whether uneven relatedness or missingness is creating artificial structure.

Cluster numbers and colors are arbitrary labels. Cluster 1 at K=3 is not guaranteed to correspond to Cluster 1 at K=4.

## Step 8: Plot all requested K values

The accompanying generic plotting program creates one ancestry bar plot for every available K in the requested range and a CV-error plot:

```bash
python plot_admixture.py \
    --base-dir "${OUT_DIR}" \
    --prefix "${PREFIX}.admixture" \
    --k-min "${K_MIN}" \
    --k-max "${K_MAX}" \
    --output-dir "${OUT_DIR}/plots"
```

Samples are ordered first by their assigned major ancestry component (Cluster 1 through Cluster K). Within each component, samples are sorted by their membership in that major component from highest to lowest. Consequently, the least-admixed representatives appear first within each cluster block and increasingly admixed individuals appear later.

For K=1 through K=10, this produces ten ancestry bar plots plus one CV-error summary plot. The fixed cluster colors begin with:

```text
#EFC86E  soft yellow
#97C684  light green
#6F9969  muted green
#AAB5D5  pale blue-lilac
#5C66A8  medium blue
#454A74  deep blue
#E6A4B4  soft pink
#C97C7C  muted brick red
#8FB8B8  desaturated teal
#D8C3A5  warm beige
```

These reproduce the intended palette from the earlier Brassica figures while keeping the new plotting program independent of cohort names and filesystem paths.

## Output interpretation

For each K:

```text
*.K.Q       one row per sample; K ancestry proportions per row
*.K.P       one row per variant; allele-frequency estimates for K clusters
*.K.log     optimization progress and cross-validation error
```

The row order of `.Q` matches the row order of `.fam`. Never sort one file independently before joining sample IDs to ancestry proportions.

ADMIXTURE ancestry components are statistical clusters, not automatically species, subspecies, or historical populations. Biological names should be assigned only after comparison with metadata and independent evidence.

## Suggested directory structure

```text
admixture/
├── run_admixture.md
├── plot_admixture.py
└── results/
    ├── cohort.raw.*
    ├── cohort.qc.*
    ├── cohort.ld.prune.in
    ├── cohort.ld.prune.out
    ├── cohort.admixture.bed
    ├── cohort.admixture.bim
    ├── cohort.admixture.fam
    ├── cohort.admixture.K*.log
    ├── cohort.admixture.*.Q
    ├── cohort.admixture.*.P
    └── plots/
```

## References

- [ADMIXTURE software manual](https://dalexander.github.io/admixture/admixture-manual.pdf)
- [PLINK 1.9 input formats and chromosome handling](https://www.cog-genomics.org/plink/1.9/input)
- [PLINK 1.9 LD pruning](https://www.cog-genomics.org/plink/1.9/ld)
- [PLINK 1.9 input filtering](https://www.cog-genomics.org/plink/1.9/filter)
