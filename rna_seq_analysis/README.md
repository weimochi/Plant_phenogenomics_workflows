```bash

rna_seq_analysis/
├── 01. QC / Trimming
│   ├── fastp & FastQC / MultiQC
│   ├── 💡 Watch out: 
│   │   ├── Check 5'/3' biases (common in plant hexamer priming)
│   │   └── Enable --trim_poly_g for NovaSeq platforms
│   └── Output: cleaned FASTQ and read-quality reports
│
├── 02. Mapping
│   ├── STAR / GSNAP
│   ├── 💡 Watch out: 
│   │   ├── Tune plant-specific intron limits (shorter than mammals)
│   │   └── Handle multihitting reads carefully for polyploids/homeologs
│   └── Output: sorted BAM & mapping summary
│
├── 03. Mapping QC
│   ├── mosdepth / Alignment metrics
│   ├── 💡 Watch out: 
│   │   ├── Focus on Uniquely Mapped Rate, not just Overall Mapping Rate
│   │   └── Check gene-body coverage for RNA degradation (3'/5' bias)
│   └── Output: Flagged low-quality samples
│
├── 04. Read Counting
│   ├── htseq-count / featureCounts
│   ├── 💡 Watch out: 
│   │   ├── Match strandedness precisely (dUTP/directional vs unstranded)
│   │   └── Ensure GTF feature type and ID attribute alignment
│   └── Output: Per-sample gene count files
│
├── 05. Expression Matrix & Preprocessing
│   ├── Merge per-sample count files (Gene × Sample raw count matrix)
│   ├── 1. Sample QC (Remove low depth & inconsistent samples)
│   ├── 2. Annotation Matching (Verify valid CDS length)
│   ├── 3. Normalization & Transformation (CPM, TPM, log2(TPM+1))
│   ├── 4. Expression QC (PCA & Correlation Heatmap)
│   ├── 5. Gene Filtering (Low-expression & MAD variability filter)
│   └── 6. Matrix Transposition (Sample × Gene for WGCNA)
│
└── 06. Downstream Analysis
    ├── WGCNA
    ├── Module-trait correlation
    ├── Hub gene identification
    └── Functional analysis
        ├── GO enrichment
        └── KEGG enrichment
```