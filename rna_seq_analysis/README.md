```bash

rna_seq_analysis
│
├── rna_mapping_and_qc (Fastq Alignment & Quality Control)
│   └── FASTQ
│        │
│        ▼
│   [00_preprocessing] (fastp, FastQC, MultiQC)
│        │  * Raw reads trimming, adaptor removal, and quality filtering
│        ▼
│   [01_alignment] (STAR / GSNAP)
│        │  * Splice-aware alignment against reference genome (AA / CC reference)
│        ▼
│   [02_stats_and_qc] (mosdepth)
│        │  * Coverage distribution, uniquely mapped rates, and degradation checks
│        │
│        ├───────────────────────────────────────┐
│        ▼                                       ▼
│   (For Expression Analysis)           (For Variant Calling)
│        │                                       │
│        ▼                                       ▼
│  Name-Sorted BAM                       Coordinate-Sorted BAM
│        │                                       │
│        ▼                                       ▼
│  rna_expression                        rna_seq_variant_calling
│  ├── Read Counting                     ├── [01_mark_duplicates] (Picard)
│  │   └── htseq-count / featureCounts   ├── [02_deepvariant] (DeepVariant gVCF)
│  ├── Expression Matrix                 ├── [03_joint_calling] (GLnexus)
│  │   └── Gene × Sample (CPM/TPM)       └── [04_filtering] (bcftools filter)
│  └── Downstream Analysis                       │
│      ├── WGCNA Co-Expression Networks          ▼
│      ├── Normalization & PCA               Final SNP VCF
│      └── GO / KEGG Enrichment                  │
│                                                ▼
│                                        Validation / GWAS
```