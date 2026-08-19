```bash

rna_seq_analysis
│
├── rna_mapping_and_qc
│   │
│   └── FASTQ
│        │
│        ▼
│   [00_preprocessing]
│   fastp / FastQC / MultiQC
│        │
│        │  Read trimming, adapter removal,
│        │  and quality control
│        ▼
│   [01_alignment]
│   STAR / GSNAP
│        │
│        │  Splice-aware alignment against
│        │  AA / CC reference genomes
│        ▼
│   [02_stats_and_qc]
│        │
│        │  Mapping statistics, coverage,
│        │  and alignment QC
│        │
│        ├───────────────────────────────┐
│        │                               │
│        ▼                               ▼
│   Expression Analysis            Variant Calling
│        │                               │
│        ▼                               ▼
│   Name-Sorted BAM               Coordinate-Sorted BAM
│        │                               │
│        ▼                               ▼
│
├── rna_expression                ├── rna_seq_variant_calling
│   │                             │
│   ├── Read Counting             ├── Mark Duplicates
│   │   htseq-count /             │   Picard
│   │   featureCounts             │
│   │                             ├── DeepVariant
│   ▼                             │   └── per-sample gVCF
│   Raw Count Matrix              │
│   Gene × Sample                 ├── GLnexus
│   │                             │   └── Cohort VCF
│   ▼                             │
│   Normalization                 ├── bcftools filtering
│   ├── CPM                       │
│   ├── TPM                       ▼
│   └── log2(TPM + 1)         Final SNP VCF
│   │                             │
│   ▼                             ▼
│   Downstream Analysis       Downstream Analysis
│   ├── PCA / Correlation         ├── Validation
│   ├── WGCNA                     ├── Population Analysis
│   ├── Differential Expression   └── GWAS
│   └── GO / KEGG
```