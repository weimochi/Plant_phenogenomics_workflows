# Population genetics workflow

This directory contains reusable workflows for population structure, genetic similarity, differentiation, diversity, inbreeding, and linkage disequilibrium analyses of quality-filtered DNA variant datasets.

```text
Quality-filtered cohort VCF
│
├── Common SNPs + LD pruning
│   │
│   ├── pca
│   │   └── Major axes of genetic variation
│   │
│   ├── admixture
│   │   └── Individual ancestry proportions across K models
│   │
│   ├── ibs
│   │   └── Pairwise genetic similarity and clustered heatmap
│   │
│   └── inbreeding
│       └── Individual observed versus expected homozygosity
│
├── Quality-filtered, non-LD-pruned variants
│   │
│   ├── fst
│   │   └── Pairwise population differentiation and genome scans
│   │
│   └── nucleotide_diversity
│       └── Within-population nucleotide diversity (π)
│
└── Quality-filtered, non-LD-pruned SNPs
    │
    └── ld_decay
        └── Decline of genotype r² with physical distance
```

## Modules

| Module | Main question | Core tool | Entry document |
|---|---|---|---|
| PCA | What are the major axes of genetic variation among samples? | PLINK | [`pca/run_pca.md`](pca/run_pca.md) |
| ADMIXTURE | What ancestry components are inferred for each individual across different K values? | PLINK + ADMIXTURE | [`admixture/run_admixture.md`](admixture/run_admixture.md) |
| IBS | Which samples share the most alleles and cluster genetically? | PLINK | [`ibs/run_ibs.md`](ibs/run_ibs.md) |
| Inbreeding | Does an individual show excess homozygosity relative to Hardy–Weinberg expectations? | PLINK | [`inbreeding/run_inbreeding.md`](inbreeding/run_inbreeding.md) |
| FST | How strongly are predefined populations differentiated? | VCFtools | [`fst/run_fst.md`](fst/run_fst.md) |
| Nucleotide diversity | How much genetic diversity exists within each population? | VCFtools | [`nucleotide_diversity/run_nucleotide_diversity.md`](nucleotide_diversity/run_nucleotide_diversity.md) |
| LD decay | How quickly does linkage disequilibrium decline with physical distance? | PLINK | [`ld_decay/run_ld_decay.md`](ld_decay/run_ld_decay.md) |

## Important input distinction

The same filtered VCF can be the common starting point, but the final marker set should not be reused blindly across every analysis.

- PCA, ADMIXTURE, IBS, and marker-based inbreeding estimates generally benefit from common, approximately independent SNPs produced by LD pruning.
- FST uses population-defined, non-LD-pruned variants. Windowed results retain physical linkage by design.
- Nucleotide diversity must not use an LD-pruned dataset, and strong MAF filtering should be avoided because removing rare variants biases π downward.
- LD decay must not use an LD-pruned dataset because correlated markers are the signal being measured. A documented MAF filter is often used to stabilize `r²` estimates.

All population comparisons should use consistent genotype quality, missingness, callable regions, sample definitions, chromosome representation, and reference-genome treatment.

## Population definitions

Population sample lists may be based on taxonomy, geography, breeding history, or another independently justified classification. ADMIXTURE-derived groups are also supported for descriptive analyses.

When ADMIXTURE is used to define populations for FST or diversity calculations, report the selected K, ancestry threshold, treatment of admixed individuals, and minimum group size. Results derived from those groups describe the ADMIXTURE-defined clusters and should not be treated as independent validation of the same clustering solution.

## Recommended analysis order

```text
Variant and sample QC
        │
        ├── PCA / IBS
        │       └── identify broad structure, relatedness, and outliers
        │
        ├── ADMIXTURE across a user-defined K range
        │       └── examine CV error and ancestry stability
        │
        ├── define biologically meaningful populations
        │       ├── FST
        │       ├── nucleotide diversity
        │       └── inbreeding within suitable frequency groups
        │
        └── LD decay within suitable populations
```

The analyses are complementary rather than a single mandatory linear pipeline. Their order should follow the biological question and population definitions.

## General cautions for Brassica data

- Confirm that genotype calls and analysis methods match the biological ploidy.
- PLINK 1.9 estimators used here generally assume diploid genotype coding.
- Analyze AA, CC, allopolyploid, or mixed-ploidy material separately when their genotype representations are not directly comparable.
- Reference bias can affect missingness, heterozygosity, π, FST, IBS, and LD when divergent accessions are aligned to one reference genome.
- Population mixture can inflate apparent inbreeding, diversity, and long-range LD.
- High FST or low π does not by itself demonstrate selection.

Each module contains its own principles, adjustable parameters, input validation, commands, outputs, and interpretation notes.
