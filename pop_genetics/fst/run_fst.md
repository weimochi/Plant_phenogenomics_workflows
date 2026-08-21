# Pairwise FST analysis

This module calculates pairwise Weir–Cockerham FST between predefined populations and visualizes both the genome-wide distribution and windowed genome scan.

## Principle

FST measures genetic differentiation among populations by comparing allele-frequency variation within populations with variation among populations. Conceptually:

```text
FST ≈ (total genetic variation - within-population variation)
      / total genetic variation
```

`FST ≈ 0` indicates similar allele frequencies, while larger positive values indicate stronger differentiation. VCFtools implements the Weir and Cockerham (1984) estimator, which accounts for finite sample sizes and unequal population sizes.

Negative estimates can occur through sampling variation, especially with small populations or weak differentiation. They do not imply biologically negative differentiation. Retain raw values for transparency; if a downstream summary requires a bounded descriptive scale, state explicitly when negative values are truncated to zero.

## Workflow

```text
Filtered cohort VCF
        │
        ▼
Population sample lists
        │
        ├── supplied independently, or
        └── derived from ADMIXTURE with a declared threshold
        │
        ▼
All pairwise population comparisons
        │
        ▼
Windowed Weir–Cockerham FST
        │
        ├── distribution plot
        └── chromosome genome scan
```

## Required software

```bash
module load vcftools
module load bcftools
```

Plotting requires Python 3 with NumPy, pandas, Matplotlib, and seaborn.

## Step 1: Define populations

Each population file must contain exactly one VCF sample ID per line:

```text
sample_001
sample_002
sample_003
```

Use biological or independently defined populations whenever possible. Place lists in one directory and name them with the suffix `.samples.txt`:

```text
groups/
├── population_A.samples.txt
├── population_B.samples.txt
└── population_C.samples.txt
```

### Optional: derive groups from ADMIXTURE

```bash
python make_fst_groups_from_admixture.py \
    --q-file /path/to/cohort.admixture.4.Q \
    --fam-file /path/to/cohort.admixture.fam \
    --k 4 \
    --threshold 0.80 \
    --minimum-size 3 \
    --output-dir /path/to/fst/groups
```

Each sample is assigned to its largest ancestry component. Samples whose maximum ancestry proportion is below the threshold are written to `admixed.samples.txt` and excluded from population lists.

This definition must be reported, for example:

> Populations were operationally defined from the K=4 ADMIXTURE solution using a minimum major-ancestry proportion of 0.80.

Using the same variants to infer ADMIXTURE groups and then demonstrate differentiation between those groups is partly circular. Such FST is useful for describing the clusters, but it is not independent evidence that the clusters exist. Compare the result with taxonomy, geography, PCA, IBS, and externally defined populations.

## Step 2: Set calculation parameters

```bash
VCF="/path/to/final.filtered.vcf.gz"
GROUP_DIR="/path/to/fst/groups"
OUT_DIR="/path/to/fst/results"
WINDOW_SIZE=50000
WINDOW_STEP=10000
```

The example uses 50-kb windows advancing by 10 kb. Overlapping windows produce a smoother scan but are not independent observations. Appropriate values depend on marker density, LD decay, chromosome size, recombination, and the biological scale of interest.

## Step 3: Run all pairwise comparisons

```bash
bash run_pairwise_fst.sh \
    "${VCF}" \
    "${GROUP_DIR}" \
    "${OUT_DIR}" \
    "${WINDOW_SIZE}" \
    "${WINDOW_STEP}"
```

The script validates population lists against VCF samples, rejects overlapping groups, and runs every unique pair once.

For populations A and B, the underlying command is:

```bash
vcftools \
    --gzvcf "${VCF}" \
    --weir-fst-pop "${GROUP_DIR}/population_A.samples.txt" \
    --weir-fst-pop "${GROUP_DIR}/population_B.samples.txt" \
    --fst-window-size "${WINDOW_SIZE}" \
    --fst-window-step "${WINDOW_STEP}" \
    --out "${OUT_DIR}/population_A_vs_population_B"
```

VCFtools produces `*.windowed.weir.fst` with columns including:

```text
CHROM  BIN_START  BIN_END  N_VARIANTS  WEIGHTED_FST  MEAN_FST
```

- `WEIGHTED_FST` combines Weir–Cockerham numerator and denominator contributions across sites and is generally preferred for window summaries.
- `MEAN_FST` is the arithmetic mean of per-site estimates and can be unstable when sites have small denominators.
- `N_VARIANTS` records how many variants informed the window; low-density windows should be interpreted cautiously.

## Step 4: Plot results

```bash
python plot_fst.py \
    --input-dir "${OUT_DIR}" \
    --pattern '*.windowed.weir.fst' \
    --fst-column WEIGHTED_FST \
    --minimum-variants 5 \
    --output-prefix "${OUT_DIR}/pairwise_fst"
```

Outputs:

```text
pairwise_fst.summary.tsv
pairwise_fst.distribution.png
pairwise_fst.distribution.pdf
pairwise_fst.genome_scan.png
pairwise_fst.genome_scan.pdf
```

The distribution plot compares all pairwise population contrasts. The genome scan lays chromosomes end to end and marks the empirical 95th percentile separately for each comparison. This percentile is an exploratory outlier guide, not a formal significance threshold.

Use `--contigs-file` to restrict and order chromosomes. The file contains one chromosome name per line in desired order:

```bash
python plot_fst.py \
    --input-dir "${OUT_DIR}" \
    --contigs-file /path/to/main_chromosomes.txt \
    --output-prefix "${OUT_DIR}/pairwise_fst"
```

## Major interpretation cautions

- FST depends on the population definition and allele-frequency spectrum.
- Very small or unequal groups produce noisier estimates.
- Genotype missingness and depth differences between populations can mimic differentiation.
- Reference bias can be severe when divergent Brassica groups are mapped to one reference.
- MAF filtering can change FST and must be reported.
- Overlapping windows must not be treated as independent replicates in statistical tests.
- High FST identifies differentiation, not necessarily selection; low diversity, linked selection, inversions, assembly problems, and low recombination can also generate peaks.
- Compare FST peaks with π, coverage, missingness, recombination, gene annotation, and independent biological evidence.

## References

- [VCFtools FST documentation](https://vcftools.github.io/man_latest.html)
- Weir BS, Cockerham CC. 1984. Estimating F-statistics for the analysis of population structure. Evolution 38:1358–1370.
