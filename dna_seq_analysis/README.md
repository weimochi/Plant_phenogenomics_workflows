# DNA-seq analysis workflow

```text
dna_seq_analysis
│
├── dna_mapping_and_qc
│   │
│   └── Paired-end FASTQ
│        │
│        ▼
│   [00a_raw_read_qc]
│   FastQC
│        │
│        │  Assess raw-read quality, base composition,
│        │  adapter content, and sequence duplication
│        ▼
│   [00b_read_trimming]  (optional)
│   fastp
│        │
│        │  Remove adapters and low-quality bases,
│        │  and discard reads that fail quality criteria
│        ▼
│   [00c_trimmed_read_qc]  (after trimming)
│   FastQC
│        │
│        │  Confirm that trimming improved read quality
│        │  without introducing unexpected bias
│        ▼
│   [01_alignment]
│   BWA-MEM2
│        │
│        │  Align DNA reads to the reference genome
│        │  and produce coordinate-sorted BAM files
│        ▼
│   [01b_duplicate_marking]
│   GATK MarkDuplicates
│        │
│        │  Identify and mark PCR or optical duplicates
│        ▼
│   [01c_mapping_qc]
│   samtools
│        │
│        │  Calculate per-sample alignment, coverage,
│        │  insert-size, and duplication statistics
│        ▼
│   [01d_mapping_qc_summary]
│        │
│        │  Combine per-sample metrics into a
│        │  cohort-level mapping QC summary
│        ▼
│   Analysis-ready BAM + BAI
│        │
│        ▼
├── dna_variant_calling
│   │
│   ├── [01_per_sample_variant_calling]
│   │   DeepVariant
│   │        │
│   │        └── Per-sample VCF + gVCF
│   │                    │
│   │                    ▼
│   ├── [02_joint_calling]
│   │   GLnexus
│   │        │
│   │        └── Cohort VCF
│   │                    │
│   │                    ▼
│   ├── [03_baseline_variant_qc]
│   │   bcftools
│   │        │
│   │        │  Summarize sites, samples, filters,
│   │        │  contigs, and variant-level metrics
│   │        ▼
│   ├── [04_variant_filtering]
│   │   bcftools
│   │        │
│   │        │  Apply project-specific site and
│   │        │  genotype filtering criteria
│   │        ▼
│   └── [05_filtered_variant_qc]
│       bcftools
│            │
│            │  Recalculate QC metrics and compare
│            │  the dataset before and after filtering
│            ▼
│       Final filtered cohort VCF
│            │
│            ▼
└── Downstream analysis
    ├── Population structure and diversity
    ├── Phylogenetic analysis
    ├── GWAS
    ├── Genotype–phenotype association
    └── Candidate variant validation
```
