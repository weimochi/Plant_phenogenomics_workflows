# Linkage disequilibrium decay with PLINK

This module calculates pairwise linkage disequilibrium (LD) with PLINK 1.9 and plots the decline of mean `r²` with increasing physical distance.

## Principle

LD is the non-random association of alleles at different loci. If alleles at two variants are inherited independently, their association is low. Nearby variants often remain associated because recombination has had fewer opportunities to separate them.

PLINK's unphased `--r2` calculation uses the squared correlation between diploid genotype allele counts at two variants:

```text
r² = correlation(genotype dosage at variant A,
                 genotype dosage at variant B)²
```

Interpretation:

```text
r² ≈ 1  strong association
r² ≈ 0  weak association
```

LD decay is visualized by grouping SNP pairs according to physical distance and calculating mean `r²` in every distance bin. A rapidly declining curve indicates shorter-range LD; a slowly declining curve indicates longer-range LD.

LD patterns reflect recombination, effective population size, selfing, selection, demographic history, population mixture, marker ascertainment, and genome structure. LD decay is not a direct recombination-rate estimate.

## Critical input rule

Do not use an LD-pruned dataset. LD pruning removes correlated variants—the signal this analysis is intended to measure.

Use a quality-filtered, non-LD-pruned, biallelic SNP dataset. Missingness and genotype quality filters are appropriate. A moderate MAF filter is commonly applied because `r²` estimates for extremely rare variants are unstable, but the selected cutoff changes the result and must be reported.

## Workflow

```text
QC-filtered, non-LD-pruned PLINK dataset
        │
        ├── whole cohort, or
        └── population sample lists
                 │
                 ▼
        Pairwise SNP r² within maximum distance
                 │
                 ▼
        Bin pairs by physical distance
                 │
                 ▼
        Mean r² decay curves
```

## Required software

```bash
module load plink
```

Plotting requires Python 3 with pandas, NumPy, and Matplotlib.

## Prepare a non-pruned PLINK dataset

```bash
VCF="/path/to/final.filtered.vcf.gz"
OUT_DIR="/path/to/pop_genetics/ld_decay/results"
PREFIX="brassica_cohort"

GENO_MAX=0.10
MIND_MAX=0.10
MAF_MIN=0.05
```

```bash
mkdir -p "${OUT_DIR}"
```

```bash
plink \
    --vcf "${VCF}" \
    --double-id \
    --allow-extra-chr \
    --vcf-half-call missing \
    --biallelic-only strict \
    --snps-only just-acgt \
    --mind "${MIND_MAX}" \
    --geno "${GENO_MAX}" \
    --maf "${MAF_MIN}" \
    --make-bed \
    --out "${OUT_DIR}/${PREFIX}.ld_input"
```

Inspect the log, sample count, and variant count before continuing.

PLINK LD calculations consider founders. If meaningful pedigree relationships are encoded in `.fam`, retain that interpretation deliberately. If all samples should be treated as unrelated founders and pedigree fields are merely placeholders, prepare the `.fam` consistently rather than adding `--make-founders` without understanding its effect.

## Set LD parameters

```bash
BFILE="${OUT_DIR}/${PREFIX}.ld_input"
GROUP_DIR="/path/to/population_sample_lists"
MAX_DISTANCE_KB=1000
MAX_VARIANTS=999999
MIN_R2=0
THREADS=8
```

- `MAX_DISTANCE_KB`: maximum physical distance between SNPs included in the report.
- `MAX_VARIANTS`: prevents the variant-count window from truncating pairs before the physical-distance limit.
- `MIN_R2=0`: retains low-r² pairs. PLINK otherwise omits table entries below its default threshold of 0.2, which would severely bias an LD-decay curve upward.
- `THREADS`: CPUs used by PLINK where supported.

Pair counts grow rapidly with variant density and maximum distance. A request covering millions of dense variants can create extremely large files. Pilot one chromosome or a smaller maximum distance first and inspect output size.

## Calculate the whole cohort

```bash
plink \
    --bfile "${BFILE}" \
    --allow-extra-chr \
    --r2 gz \
    --ld-window "${MAX_VARIANTS}" \
    --ld-window-kb "${MAX_DISTANCE_KB}" \
    --ld-window-r2 "${MIN_R2}" \
    --threads "${THREADS}" \
    --out "${OUT_DIR}/whole_cohort.ld"
```

Expected output:

```text
whole_cohort.ld.ld.gz
```

The repeated `.ld.ld.gz` comes from the chosen output prefix plus PLINK's `.ld.gz` suffix. The automated script uses a cleaner population prefix.

## Calculate each population automatically

Population files contain one FID and IID pair per line because PLINK `--keep` expects two ID columns:

```text
sample_001 sample_001
sample_002 sample_002
```

```bash
bash run_ld_decay.sh \
    "${BFILE}" \
    "${GROUP_DIR}" \
    "${OUT_DIR}" \
    "${MAX_DISTANCE_KB}" \
    "${THREADS}"
```

Use `-` for `GROUP_DIR` to calculate the entire cohort only:

```bash
bash run_ld_decay.sh \
    "${BFILE}" - "${OUT_DIR}" "${MAX_DISTANCE_KB}" "${THREADS}"
```

The wrapper fixes `--ld-window-r2 0` and uses a very large variant-count window so the user-defined physical distance is the primary limit.

## Plot LD decay

```bash
python plot_ld_decay.py \
    --input-dir "${OUT_DIR}" \
    --pattern '*.ld.gz' \
    --bin-size-bp 10000 \
    --maximum-distance-bp 1000000 \
    --output-prefix "${OUT_DIR}/ld_decay"
```

The plotting program reads large PLINK reports in chunks, calculates distance as `abs(BP_B - BP_A)`, and accumulates the mean `r²` and pair count per distance bin.

Outputs:

```text
ld_decay.binned.tsv
ld_decay.png
ld_decay.pdf
```

Always inspect `n_pairs` in the TSV. The farthest bins or small populations may contain few SNP pairs and yield unstable means.

## Choosing distance bins

- Smaller bins give finer resolution but noisier estimates.
- Larger bins smooth the curve and are more stable.
- Use the same bin size and maximum distance for populations being compared.
- Report whether x coordinates represent the bin start, center, or end; this module plots bin centers.

A visual threshold such as the distance where mean `r²` falls below 0.2 can be reported descriptively, but it is sensitive to MAF, sample size, marker density, and smoothing. It should not be treated as a universal biological boundary.

## Interpretation cautions

- Do not use LD-pruned input.
- Use comparable MAF and missingness filters across groups.
- Population mixture can create long-range LD.
- Small sample sizes inflate variance and may bias `r²`.
- Selfing populations commonly show slower LD decay.
- Structural variants, inversions, low recombination, and assembly errors can extend LD locally.
- Reference bias and paralogous mapping are important in divergent Brassica genomes.
- PLINK's ordinary `--r2` is an unphased genotype-correlation measure; haplotype-based D′ or phased LD answers a related but different question.
- True polyploid genotypes require a ploidy-aware LD method.

## Suggested structure

```text
ld_decay/
├── run_ld_decay.md
├── run_ld_decay.sh
├── plot_ld_decay.py
└── results/
```

## Reference

- [PLINK 1.9 linkage-disequilibrium documentation](https://www.cog-genomics.org/plink/1.9/ld)
