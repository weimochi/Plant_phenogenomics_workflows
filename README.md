# Plant Phenogenomics Workflows

Reusable workflows and analysis utilities for integrating plant genomic, transcriptomic, population-genetic, and image-derived phenotypic data.

This repository collects practical workflows developed for plant phenogenomics projects, with particular attention to complex plant genomes such as *Brassica*. It includes HPC-oriented sequencing pipelines, population-genetic analyses, image preprocessing and feature interpretation tools, and cross-modality distance-matrix comparison.

> **Project status:** This is an evolving research workflow collection rather than a single turnkey software package. Modules differ in maturity and should be reviewed and adapted before use on a new dataset or computing environment.

## Repository overview

```text
Plant_phenogenomics_workflows
│
├── dna_seq_analysis
│   ├── DNA read QC, trimming, and BWA-MEM2 alignment
│   └── DeepVariant calling, GLnexus joint calling, filtering, and QC
│
├── rna_seq_analysis
│   ├── RNA read QC and splice-aware alignment
│   ├── Gene-expression quantification and matrix analysis
│   └── RNA-seq variant-calling workflow
│
├── pop_genetics
│   ├── PCA and ADMIXTURE
│   ├── IBS and inbreeding
│   ├── FST and nucleotide diversity
│   └── LD decay
│
├── phenotyping_vision
│   ├── Image color calibration and plant segmentation
│   └── ResNet50 feature interpretation and exploration
│
├── multiomics_integration
│   └── Distance-matrix concordance across data modalities
│
└── docs
    └── Recorded software versions
```

## Modules

| Module | Purpose | Main inputs | Main outputs | Current contents |
|---|---|---|---|---|
| [`dna_seq_analysis`](dna_seq_analysis/) | DNA read processing and germline small-variant analysis | Paired-end FASTQ, reference genome | Analysis-ready BAM, per-sample gVCF, cohort VCF, QC summaries | Slurm workflows and detailed documentation |
| [`rna_seq_analysis`](rna_seq_analysis/) | RNA alignment, expression analysis, and RNA-derived variant calling | RNA-seq FASTQ, genome, annotation | BAM, raw count matrix, expression QC, co-expression results, RNA-derived VCF | Workflow documentation; implementation coverage varies by submodule |
| [`pop_genetics`](pop_genetics/) | Population structure, genetic diversity, differentiation, and LD analyses | Quality-filtered cohort VCF | PCA, ancestry estimates, IBS, FST, π, inbreeding, and LD-decay results | Analysis protocols, shell/Python tools, and interpretation notes |
| [`phenotyping_vision`](phenotyping_vision/) | Standardize plant images and explore pretrained visual representations | Plant images and optional sample metadata | Calibrated/segmented images and feature-interpretation figures | Interactive preprocessing and ResNet50 visualization scripts |
| [`multiomics_integration`](multiomics_integration/distance_matrix_concordance/) | Compare sample relationships across omics and phenotype modalities | Feature tables or labeled distance/similarity matrices | Mantel and Procrustes results, PCoA diagnostics, heatmaps, and audit tables | Executable Python workflow and documentation |

## End-to-end concept

The modules can be combined into a broader phenogenomics workflow:

```text
DNA-seq FASTQ
    │
    ├── QC and alignment
    ├── Variant calling and joint genotyping
    └── Population-genetic analyses
                 │
                 └── Genotype distance ──────────────┐

RNA-seq FASTQ                                      │
    │                                               │
    ├── QC and splice-aware alignment               │
    ├── Gene-level expression matrix                │
    └── Expression distance ────────────────────────┤
                                                    │
Plant images                                       ├── Distance-matrix
    │                                               │   concordance
    ├── Color calibration and segmentation          │
    ├── Visual feature exploration                  │
    └── Image-feature distance ─────────────────────┘
```

Each modality should be quality-controlled independently before cross-modality comparison. Sample identifiers and biological observation units must be harmonized explicitly.

## DNA-seq workflow

The DNA-seq module provides the most complete HPC workflow in this repository:

```text
Paired-end FASTQ
    ↓
FastQC
    ↓
Optional fastp trimming
    ↓
BWA-MEM2 alignment
    ↓
GATK MarkDuplicates
    ↓
samtools / mosdepth mapping QC
    ↓
DeepVariant per-sample calling
    ↓
Per-sample VCF + gVCF
    ↓
GLnexus joint calling
    ↓
Cohort VCF
    ↓
Dataset-specific filtering and pre/post-filter QC
```

Start with:

- [`dna_seq_analysis/README.md`](dna_seq_analysis/README.md)
- [`dna_seq_analysis/dna_mapping_and_qc/README.md`](dna_seq_analysis/dna_mapping_and_qc/README.md)
- [`dna_seq_analysis/dna_variant_calling/README.md`](dna_seq_analysis/dna_variant_calling/README.md)

The supplied batch scripts are designed for Slurm-based HPC systems. Resource requests, module names, container paths, reference files, and scratch-storage locations must be adapted to the target cluster.

## RNA-seq workflow

The RNA-seq documentation describes three related branches:

```text
RNA-seq FASTQ
    ↓
FastQC / fastp
    ↓
STAR or GSNAP splice-aware alignment
    ↓
Mapping and coverage QC
    ├── Gene-level counting
    │       ↓
    │   Raw count matrix
    │       ↓
    │   Normalization, PCA, correlation, and WGCNA
    │
    └── Duplicate marking
            ↓
        DeepVariant
            ↓
        GLnexus joint calling
            ↓
        Filtered RNA-derived variant set
```

Module documentation:

- [`rna_seq_analysis/rna_mapping_and_qc/`](rna_seq_analysis/rna_mapping_and_qc/)
- [`rna_seq_analysis/rna_expression/`](rna_seq_analysis/rna_expression/)
- [`rna_seq_analysis/rna_seq_variant_calling/`](rna_seq_analysis/rna_seq_variant_calling/)

RNA-derived variant calls require especially careful interpretation. Coverage depends on gene expression, transcript structure, tissue, condition, allele-specific expression, mapping bias, and RNA processing. They should not automatically be treated as equivalent to variants obtained from genomic DNA.

## Population genetics

The [`pop_genetics`](pop_genetics/) module begins with a quality-filtered cohort VCF and supports:

| Analysis | Main question |
|---|---|
| PCA | What are the major axes of genetic variation? |
| ADMIXTURE | What ancestry components are inferred across different values of K? |
| IBS | Which samples are most genetically similar? |
| Inbreeding | Which individuals show excess homozygosity? |
| FST | How strongly are predefined populations differentiated? |
| Nucleotide diversity | How much genetic diversity exists within each population? |
| LD decay | How quickly does genotype correlation decline with physical distance? |

Marker preparation is analysis-dependent. In particular, LD-pruned SNPs may be useful for PCA or ADMIXTURE but must not be reused blindly for nucleotide-diversity or LD-decay analyses.

See [`pop_genetics/README.md`](pop_genetics/README.md) for recommended analysis order, input distinctions, and cautions for *Brassica* datasets.

## Image phenotyping

The [`phenotyping_vision`](phenotyping_vision/) module currently contains two complementary components.

### Image preprocessing

[`phenotyping_vision/image_preprocessing`](phenotyping_vision/image_preprocessing/) provides:

1. reference-patch-relative white balance using a Kodak Q-13 color card;
2. interactive plant segmentation using Segment Anything Model priors and HSV constraints;
3. standardized cropped and padded plant-image export;
4. calibration logging for later batch-effect auditing.

### Feature interpretation

[`phenotyping_vision/feature_extraction`](phenotyping_vision/feature_extraction/) uses an ImageNet-pretrained ResNet50 to generate exploratory visualizations:

- Eigen-CAM activation maps;
- multi-scale spatial feature clustering;
- feature-space occlusion sensitivity.

These tools help identify which image regions influence pretrained representations. Their outputs are exploratory and should not be interpreted automatically as validated biological traits.

## Multi-omics integration

The [`distance_matrix_concordance`](multiomics_integration/distance_matrix_concordance/) module compares sample geometry across genotype, transcriptome, image, phenotype, metabolite, or other high-dimensional data.

```text
Modality-specific feature tables
    ↓
Validated sample × sample distance matrices
    ↓
Sample intersection and matrix alignment
    ↓
Pearson / Spearman effect sizes
    ↓
Mantel permutation tests
    ↓
PCoA and Procrustes comparison
    ↓
Multiple-testing correction and diagnostic plots
```

The workflow includes explicit sample-alignment auditing because mislabeled samples, unmatched biological units, and pseudoreplication can produce misleading cross-modality associations.

Start with [`multiomics_integration/distance_matrix_concordance/README.md`](multiomics_integration/distance_matrix_concordance/README.md).

## Getting started

Clone the repository:

```bash
git clone https://github.com/weimochi/Plant_phenogenomics_workflows.git
cd Plant_phenogenomics_workflows
```

Then:

1. select the module matching your input data and research question;
2. read that module's README completely;
3. inspect script parameters, expected file naming, and input manifests;
4. confirm the biological assumptions, including ploidy and reference-genome treatment;
5. adapt software environments and Slurm resources to the local system;
6. run a small pilot dataset before launching a full cohort;
7. record all software versions, parameters, references, and sample exclusions.

Software versions recorded for previous analyses are listed in [`docs/software_versions.md`](docs/software_versions.md). They document prior environments and are not guaranteed to be universal requirements for every module.

## Important considerations

### Complex and polyploid genomes

Many standard variant and population-genetic tools assume diploid genotype representations. Before applying these workflows to polyploid or mixed-ploidy material, evaluate:

- expected ploidy;
- subgenome representation;
- homoeologous read mapping;
- reference bias;
- model compatibility;
- whether subgenome-specific analyses are required.

A workflow completing without an error does not guarantee biologically valid genotype calls.

### Reproducibility

For each analysis, preserve:

- reference assembly and annotation versions;
- sample and file manifests;
- software and container versions;
- command-line parameters;
- Slurm resource settings;
- random seeds;
- QC exclusions and their reasons;
- filtering thresholds and retained-variant counts.

### Validation

The scripts and default parameters should be treated as starting points. Thresholds must be evaluated against the sequencing platform, coverage, organism, tissue, experimental design, and biological question.

## Citation and reuse

If these workflows contribute to a publication, describe the exact module, version or commit, parameters, reference resources, and modifications used.

A formal citation file and software license have not yet been added. Check the repository status or contact the maintainer before redistributing or incorporating substantial portions into another project.

## Contributing

Issues and pull requests that improve documentation, portability, input validation, or support for additional plant-genome configurations are welcome.

When proposing a change, please include:

- the affected module;
- biological and computational assumptions;
- example input structure;
- expected output;
- required software versions;
- any HPC-specific requirements.

## Disclaimer

These workflows are provided for research use. Users are responsible for validating methods, parameters, genotype models, and biological interpretations for their own organisms and datasets.