# Nucleotide diversity (π)

This module calculates nucleotide diversity for a whole cohort or predefined populations with VCFtools and visualizes windowed estimates.

## Important tool distinction

Nucleotide diversity is calculated here with VCFtools, not PLINK.

PLINK is useful for genotype QC, PCA, ADMIXTURE input preparation, IBS, and individual inbreeding estimates, but PLINK 1.9 does not provide the standard `--window-pi` calculation used in this workflow. Population sample lists created for PLINK or ADMIXTURE can still be reused as VCFtools `--keep` files when they contain one VCF sample ID per line.

## Biological principle

Nucleotide diversity, written as π, is the expected proportion of nucleotide positions that differ between two chromosomes sampled at random from a population.

For a biallelic site with allele frequencies `p` and `q = 1 - p`, expected pairwise diversity is closely related to expected heterozygosity:

```text
π_site = 2pq
```

More generally, for multiple alleles with frequencies `p₁ ... pₖ`:

```text
π_site = 1 - Σ pᵢ²
```

A sample-size correction may be applied by the estimator. Windowed π summarizes pairwise differences across genomic intervals and expresses diversity per nucleotide position.

Interpretation:

```text
higher π  → greater within-population genetic diversity
lower π   → lower within-population diversity
```

Low π may reflect selfing, bottlenecks, selective sweeps, background selection, low recombination, or technical loss of variants. High π may reflect large effective population size, population mixture, introgression, hybrid ancestry, paralogous mapping, or genotype error.

π measures diversity within a population. FST measures differentiation between populations. They answer different questions and are most informative when interpreted together.

## Critical input choice

Do not calculate π from the LD-pruned dataset used for PCA or ADMIXTURE. LD pruning deliberately removes real variants and invalidates a nucleotide-diversity estimate.

Avoid applying a strong MAF cutoff before π calculation. Rare variants are genuine contributors to diversity; removing them biases π downward. Use a high-quality SNP dataset with genotype-level and site-level QC, but retain the allele-frequency spectrum appropriate to the study.

Absolute π also depends on the callable denominator. Variant-only files do not by themselves distinguish confidently invariant sequence from regions where variants could not be called. For rigorous cross-population or cross-species comparisons, use the same callable genomic regions, coverage criteria, reference representation, and accessibility mask for every group. Otherwise, apparent diversity differences can be produced by unequal callability.

## Workflow

```text
Quality-filtered, non-LD-pruned cohort VCF
        │
        ├── whole cohort
        └── population sample lists
                 │
                 ▼
        VCFtools windowed π
                 │
                 ├── distribution comparison
                 └── chromosome diversity scan
```

## Required software

```bash
module load vcftools
module load bcftools
```

Plotting requires Python 3 with NumPy, pandas, Matplotlib, and seaborn.

## Prepare inputs

```bash
VCF="/path/to/high_quality.non_ld_pruned.vcf.gz"
OUT_DIR="/path/to/pop_genetics/nucleotide_diversity/results"
GROUP_DIR="/path/to/population_sample_lists"
WINDOW_SIZE=50000
WINDOW_STEP=10000
```

Population files contain one VCF sample ID per line and end in `.samples.txt`:

```text
groups/
├── population_A.samples.txt
├── population_B.samples.txt
└── population_C.samples.txt
```

Validate VCF samples:

```bash
bcftools query -l "${VCF}" | head
bcftools index --stats "${VCF}" | head
```

## Calculate π for all samples

```bash
vcftools \
    --gzvcf "${VCF}" \
    --window-pi "${WINDOW_SIZE}" \
    --window-pi-step "${WINDOW_STEP}" \
    --out "${OUT_DIR}/whole_cohort"
```

## Calculate π for one population

```bash
vcftools \
    --gzvcf "${VCF}" \
    --keep "${GROUP_DIR}/population_A.samples.txt" \
    --window-pi "${WINDOW_SIZE}" \
    --window-pi-step "${WINDOW_STEP}" \
    --out "${OUT_DIR}/population_A"
```

Expected output:

```text
population_A.windowed.pi
```

Typical columns are:

```text
CHROM  BIN_START  BIN_END  N_VARIANTS  PI
```

`N_VARIANTS` should always be retained as a window-level diagnostic. A π estimate based on very few variants is less stable and may reflect low callability.

## Run every population automatically

```bash
bash run_nucleotide_diversity.sh \
    "${VCF}" \
    "${GROUP_DIR}" \
    "${OUT_DIR}" \
    "${WINDOW_SIZE}" \
    "${WINDOW_STEP}"
```

To calculate only the whole cohort, use `-` instead of a group directory:

```bash
bash run_nucleotide_diversity.sh \
    "${VCF}" - "${OUT_DIR}" "${WINDOW_SIZE}" "${WINDOW_STEP}"
```

The script validates group members against the VCF and calculates one output per population. Unlike pairwise FST, population lists do not have to be mutually exclusive for the software to run, but overlapping groups must be biologically intentional and disclosed.

## Window size and step

Example:

```text
WINDOW_SIZE = 50,000 bp
WINDOW_STEP = 10,000 bp
```

This creates overlapping 50-kb windows beginning every 10 kb.

- Smaller windows provide greater spatial resolution but contain fewer informative sites and are noisier.
- Larger windows are more stable but can hide local variation.
- `WINDOW_STEP < WINDOW_SIZE` creates overlapping windows, which are not statistically independent.
- Window choices should consider SNP density, LD decay, assembly quality, chromosome size, and the expected scale of diversity changes.

## Plot the results

```bash
python plot_nucleotide_diversity.py \
    --input-dir "${OUT_DIR}" \
    --pattern '*.windowed.pi' \
    --minimum-variants 5 \
    --output-prefix "${OUT_DIR}/nucleotide_diversity"
```

Optional chromosome ordering:

```bash
python plot_nucleotide_diversity.py \
    --input-dir "${OUT_DIR}" \
    --contigs-file /path/to/main_chromosomes.txt \
    --output-prefix "${OUT_DIR}/nucleotide_diversity"
```

Outputs:

```text
nucleotide_diversity.summary.tsv
nucleotide_diversity.distribution.png
nucleotide_diversity.distribution.pdf
nucleotide_diversity.genome_scan.png
nucleotide_diversity.genome_scan.pdf
```

The distribution plot compares windowed π among populations. The genome scan shows spatial variation along chromosomes. These windows are genomic observations, not independent biological replicates; formal population comparisons require methods that account for linkage and genomic block structure.

## Important interpretation cautions

- Use comparable sample sizes or evaluate sample-size sensitivity.
- Estimate each population from an appropriate sample list.
- Do not use LD-pruned variants.
- Avoid strong MAF filtering.
- Use consistent callable regions and filtering across groups.
- Inspect coverage, missingness, mapping quality, and reference bias.
- Exclude or separately examine duplicated and repetitive regions.
- Population mixture can inflate π.
- A low-diversity reference-biased group may simply map less effectively to the chosen genome.
- True polyploid genotype likelihoods require methods consistent with the organism's ploidy; diploid VCF genotypes do not automatically provide a valid polyploid π estimate.

## Suggested structure

```text
nucleotide_diversity/
├── run_nucleotide_diversity.md
├── run_nucleotide_diversity.sh
├── plot_nucleotide_diversity.py
└── results/
```

## Reference

- [VCFtools nucleotide-diversity documentation](https://vcftools.github.io/man_latest.html)
