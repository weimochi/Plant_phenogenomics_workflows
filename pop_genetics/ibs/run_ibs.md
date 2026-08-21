# Identity-by-state analysis with PLINK

This document describes a reusable workflow for calculating pairwise identity by state (IBS) with PLINK 1.9 and plotting a hierarchically clustered similarity heatmap.

## What IBS measures

IBS measures allele sharing between two samples without requiring the shared alleles to have been inherited from a known common ancestor. At a diploid biallelic locus, a pair of samples can share:

```text
IBS2  two alleles identical by state
IBS1  one allele identical by state
IBS0  no alleles identical by state
```

A simple per-locus similarity contribution can be represented as:

```text
IBS similarity = 1.0 for IBS2
                 0.5 for IBS1
                 0.0 for IBS0
```

PLINK aggregates information across usable variants for every pair of samples. Its `ibs` output is a similarity matrix, while `1-ibs` is a genetic-distance matrix:

```text
IBS distance = 1 - IBS similarity
```

Consequently:

- larger IBS similarity means genetically more similar samples;
- smaller IBS distance means genetically more similar samples;
- the diagonal of a similarity matrix should be 1;
- the diagonal of a distance matrix should be 0.

IBS is not the same as Pearson genotype correlation, kinship, or identity by descent (IBD). The clustered heatmap in this workflow displays allele-sharing similarity and uses IBS distance to determine sample order.

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
Pairwise IBS similarity and distance
        │
        ▼
Hierarchical clustering with 1 - IBS
        │
        ▼
Heatmap with genetically similar samples adjacent
```

## Required software

```bash
module load plink
```

The plotting program requires Python 3 with NumPy, pandas, SciPy, Matplotlib, and seaborn:

```bash
python -c 'import numpy, pandas, scipy, matplotlib, seaborn; print("IBS plotting dependencies are available")'
```

## Set paths and parameters

```bash
VCF="/path/to/final.filtered.vcf.gz"
OUT_DIR="/path/to/pop_genetics/ibs/results"
PREFIX="brassica_cohort"

GENO_MAX=0.10
MIND_MAX=0.10
MAF_MIN=0.05
LD_WINDOW=50
LD_STEP=10
LD_R2=0.20
THREADS=8
```

Create the output directory:

```bash
mkdir -p "${OUT_DIR}"
```

The example QC thresholds are starting values, not universal requirements. Use the same well-documented sample and variant set when comparing IBS with PCA or ADMIXTURE.

## Step 1: Prepare a PLINK dataset

If an appropriate QC-filtered PLINK dataset already exists, set `BFILE` to its prefix and continue to LD pruning.

Otherwise convert the filtered VCF:

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

Inspect missingness:

```bash
plink \
    --bfile "${OUT_DIR}/${PREFIX}.raw" \
    --allow-extra-chr \
    --missing \
    --out "${OUT_DIR}/${PREFIX}.raw_missingness"
```

Apply sample and variant QC:

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

Dense linked markers can cause particular genomic regions to dominate genome-wide similarity. For a population-structure overview, calculate IBS from an LD-pruned marker set:

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
    --out "${OUT_DIR}/${PREFIX}.ibs_input"
```

Check sample and marker counts:

```bash
wc -l "${OUT_DIR}/${PREFIX}.ibs_input.fam"
wc -l "${OUT_DIR}/${PREFIX}.ibs_input.bim"
```

Whether LD pruning is appropriate depends on the question. It is recommended here for genome-wide sample relationships and visualization. Analyses focused on local haplotype sharing require a different design.

## Step 3: Calculate IBS similarity and distance

PLINK can write both matrices in the same run:

```bash
plink \
    --bfile "${OUT_DIR}/${PREFIX}.ibs_input" \
    --allow-extra-chr \
    --threads "${THREADS}" \
    --distance square ibs 1-ibs \
    --out "${OUT_DIR}/${PREFIX}.ibs"
```

Expected outputs:

```text
${PREFIX}.ibs.mibs       IBS similarity matrix
${PREFIX}.ibs.mibs.id    row/column sample IDs for .mibs
${PREFIX}.ibs.mdist      IBS distance matrix (1 - IBS)
${PREFIX}.ibs.mdist.id   row/column sample IDs for .mdist
${PREFIX}.ibs.log
```

The matrix files do not contain sample names in their rows or columns. Never interpret or reorder a matrix without its corresponding `.id` file.

Inspect dimensions:

```bash
N_SAMPLES=$(wc -l < "${OUT_DIR}/${PREFIX}.ibs.mibs.id")
N_ROWS=$(wc -l < "${OUT_DIR}/${PREFIX}.ibs.mibs")

echo "samples=${N_SAMPLES} matrix_rows=${N_ROWS}"
[[ "${N_SAMPLES}" -eq "${N_ROWS}" ]] || { echo "IBS matrix dimension mismatch" >&2; exit 2; }
```

## Missing-genotype handling

Only jointly usable genotype information contributes directly to a sample pair. PLINK 1.9 applies a missingness correction to observed genomic distances. By default, the correction accounts for the expected average distance contribution of missing variants; this assumes allele frequency is approximately independent of missingness.

The optional `flat-missing` modifier uses a simpler correction based on the missing-call fraction:

```bash
--distance square ibs 1-ibs flat-missing
```

Do not add `flat-missing` automatically. Use it only when its assumption is intended and record the decision. More importantly, remove samples and variants with excessive missingness before calculating IBS, then check whether clustering is associated with missing rate or sequencing batch.

## Step 4: Cluster genetically similar samples and draw the heatmap

```bash
python plot_ibs_heatmap.py \
    --matrix "${OUT_DIR}/${PREFIX}.ibs.mibs" \
    --id-file "${OUT_DIR}/${PREFIX}.ibs.mibs.id" \
    --matrix-kind similarity \
    --linkage-method average \
    --output-prefix "${OUT_DIR}/${PREFIX}.ibs_clustered"
```

The plotting program:

1. validates that the matrix is numeric, square, finite, and symmetric;
2. converts similarity to distance with `1 - IBS`;
3. converts the distance matrix to condensed form;
4. performs agglomerative hierarchical clustering;
5. uses optimal leaf ordering to place nearby samples next to one another where possible;
6. applies exactly the same sample order to heatmap rows and columns;
7. exports the reordered similarity matrix and sample order.

The default is average-linkage clustering. It defines the distance between two clusters as the mean of all cross-cluster pairwise IBS distances. This is a transparent choice for visualization, but it is not the only possible clustering definition. `complete` and `single` linkage are also supported and can produce different dendrograms.

Alternatively, plot a PLINK distance matrix directly:

```bash
python plot_ibs_heatmap.py \
    --matrix "${OUT_DIR}/${PREFIX}.ibs.mdist" \
    --id-file "${OUT_DIR}/${PREFIX}.ibs.mdist.id" \
    --matrix-kind distance \
    --linkage-method average \
    --output-prefix "${OUT_DIR}/${PREFIX}.ibs_clustered"
```

Do not calculate `1 - matrix` when the input is already `.mdist`; `.mdist` is already an IBS-distance matrix.

## Heatmap interpretation

The displayed cell value is IBS similarity even when `.mdist` is supplied as input. Higher values and darker colors represent greater allele sharing. The hierarchical sample ordering is calculated from IBS distance.

Nearby samples in the heatmap are those grouped by the selected hierarchical-clustering algorithm. Adjacency is a visualization order, not proof of a discrete biological population. Interpret the pattern together with PCA, ADMIXTURE, taxonomy, geography, relatedness, missingness, and sequencing metadata.

Pairs with unexpectedly high similarity may represent:

- duplicated or renamed samples;
- clones or near-isogenic accessions;
- close relatives;
- sample contamination or swaps;
- genuinely low diversity within a group.

## Outputs from the plotting program

```text
*.heatmap.png
*.heatmap.pdf
*.sample_order.tsv
*.similarity.reordered.tsv
```

The sample-order table makes the clustering order explicit and reproducible instead of leaving it embedded only in the figure.

## Suggested directory structure

```text
ibs/
├── run_ibs.md
├── plot_ibs_heatmap.py
└── results/
    ├── cohort.raw.*
    ├── cohort.qc.*
    ├── cohort.ld.prune.in
    ├── cohort.ibs_input.*
    ├── cohort.ibs.mibs
    ├── cohort.ibs.mibs.id
    ├── cohort.ibs.mdist
    ├── cohort.ibs.mdist.id
    └── cohort.ibs_clustered.*
```

## Reference

- [PLINK 1.9 distance matrices and IBS calculation](https://www.cog-genomics.org/plink/1.9/distance)
- [PLINK 1.9 input filtering](https://www.cog-genomics.org/plink/1.9/filter)
- [PLINK 1.9 LD pruning](https://www.cog-genomics.org/plink/1.9/ld)
