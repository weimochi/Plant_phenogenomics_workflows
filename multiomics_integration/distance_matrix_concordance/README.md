# Distance-matrix concordance

This module provides a modality-agnostic method for comparing sample relationships across genotype, transcriptome, image, phenotype, metabolite, or other high-dimensional data.

The central design rule is:

> Convert each modality to a validated sample-by-sample distance matrix, align the same samples in the same order, and only then compare their geometric structures.

## Workflow

```text
Modality-specific feature tables
│
├── Image embeddings (sample × DINOv2 features)
├── RNA expression (gene × sample or sample × gene)
├── Quantitative phenotypes (sample × traits)
└── Other high-dimensional representations
        │
        ▼
build_distance_matrix.py
        │
        └── labeled sample × sample distance matrices

Existing distance or similarity matrices
│
├── PLINK IBS distance / similarity
├── ecological distance
└── any validated labeled matrix
        │
        ▼
matrix manifest
        │
        ▼
compare_distance_matrices.py
        │
        ├── sample-intersection audit
        ├── aligned distance matrices
        ├── Pearson / Spearman effect sizes
        ├── Mantel permutation tests
        ├── PCoA diagnostics
        ├── Procrustes permutation tests
        ├── BH multiple-testing correction
        └── heatmaps and pairwise comparison plots
```

Tree construction is intentionally excluded from the core method. A neighbor-joining tree derived from image or expression distance is a visualization of similarity, not direct evidence of evolutionary ancestry.

## Files

```text
distance_matrix_concordance/
├── README.md
├── build_distance_matrix.py
├── compare_distance_matrices.py
└── matrix_manifest.example.tsv
```

## Software requirements

```bash
python -c 'import numpy, pandas, scipy, sklearn, matplotlib, seaborn; print("Dependencies are available")'
```

## Step 1: Define the biological observation unit

Before building any matrix, decide what one row represents:

```text
one accession
one biological replicate
one tissue from one accession
one accession × condition combination
```

All compared modalities must describe the same unit. Do not silently match an accession-level genotype to several image or RNA technical replicates and treat those replicates as independent individuals.

Create and preserve a sample manifest outside this module with at least:

```text
canonical_sample_id
original_image_id
original_rna_id
original_genotype_id
biological_group
replicate
included
exclusion_reason
```

Sample-ID cleaning must be explicit. A regular expression such as `s\d{4}` can accidentally merge biologically distinct IDs such as `s0001_1` and `s0001_2`.

## Step 2: Build image-feature distance

Assume an input table with one sample per row:

```text
sample  group  feat_0  feat_1  ...  feat_383
s0001   AA     ...
s0002   CC     ...
```

For DINOv2-style embeddings, cosine distance after L2 normalization is a reasonable documented starting point:

```bash
python build_distance_matrix.py \
    --input image_features.tsv \
    --sample-axis rows \
    --sample-id-column sample \
    --exclude-columns group \
    --feature-prefix feat_ \
    --transform l2 \
    --metric cosine \
    --missing error \
    --output-prefix results/image_dinov2_cosine
```

Outputs:

```text
image_dinov2_cosine.distance.tsv
image_dinov2_cosine.metadata.json
```

Do not automatically apply feature-wise z-scoring before cosine distance. Centering every embedding dimension changes vector directions and therefore changes cosine geometry. If z-scoring is scientifically intended, specify `--transform zscore` and report it.

When several images represent one accession, avoid data leakage and pseudoreplication. Common strategies include:

- extract one embedding per image and average embeddings within accession;
- model images as repeated measurements;
- select one image using a predefined QC rule.

Pixel-wise mean or median composite images can be useful for visualization, but they may create blurred synthetic morphology. They should not automatically replace accession-level embedding aggregation.

## Step 3: Build RNA-expression distance

For a gene-by-sample expression matrix whose first column contains gene IDs:

```text
gene_id  s0001  s0002  s0003
gene_1   ...
gene_2   ...
```

Correlation distance is `1 - Pearson correlation` between sample expression profiles:

```bash
python build_distance_matrix.py \
    --input gene_log2tpm_plus1.tsv \
    --sample-axis columns \
    --transform none \
    --metric correlation \
    --missing error \
    --output-prefix results/rna_expression_correlation
```

Expression preprocessing occurs before this module and must be documented:

- gene filtering;
- normalization and transformation;
- batch correction;
- tissue and developmental-stage matching;
- treatment or environmental covariates;
- handling of technical and biological replicates.

Distance concordance can otherwise reflect batch, tissue, or growth conditions rather than shared biology.

## Step 4: Use genotype IBS distance

PLINK `.mdist` is already an IBS-distance matrix and must not be transformed with `1 - matrix` again.

PLINK `.mibs` is a similarity matrix. Declare it as `similarity` in the manifest so the comparison program converts it once to `1 - IBS`.

Always provide the matching `.mdist.id` or `.mibs.id`; the numeric matrix has no embedded sample labels.

## Step 5: Create the matrix manifest

Copy and edit `matrix_manifest.example.tsv`:

```text
name             path                                  format   kind        id_file                      id_column
genotype_ibs     /path/cohort.mdist                    plink    distance    /path/cohort.mdist.id        1
image_dinov2     /path/image_cosine.distance.tsv       labeled  distance
rna_expression   /path/rna_correlation.distance.tsv    labeled  distance
```

Fields:

- `name`: unique, filesystem-safe modality label.
- `path`: numeric matrix path.
- `format`: `labeled` or `plink`.
- `kind`: `distance` or `similarity`.
- `id_file`: required for PLINK numeric matrices.
- `id_column`: zero-based sample-ID column in the PLINK ID file; normally `1` for IID.

All labeled matrices must have sample IDs in the first row and first column.

## Step 6: Compare all matrices

```bash
python compare_distance_matrices.py \
    --manifest matrix_manifest.tsv \
    --output-dir results/concordance \
    --permutations 9999 \
    --seed 42 \
    --pcoa-dimensions 2 \
    --minimum-samples 4
```

The program uses the global intersection of samples across all matrices. The order of the first manifest matrix is retained, then applied identically to every matrix.

## Statistical outputs

### Distance-vector correlations

The upper triangles of two aligned matrices are compared using Pearson and Spearman correlations. These coefficients are descriptive effect sizes.

Ordinary correlation p-values are not reported as primary inference because pairwise distances are not independent observations: every sample contributes to many matrix cells.

### Mantel permutation test

The Mantel statistic is the Pearson correlation between corresponding upper-triangle distances. Significance is estimated by permuting sample labels of one entire matrix while preserving its internal distance structure.

```text
p = (number of |permuted r| ≥ |observed r| + 1)
    / (number of permutations + 1)
```

Mantel analysis tests matrix-level association but has limitations, including sensitivity to spatial structure and limited ability to distinguish direct from confounded relationships. Report effect size, permutation design, sample count, and number of permutations.

### PCoA and Procrustes

Each distance matrix is separately converted to principal-coordinate space. Procrustes analysis then translates, scales, and rotates two configurations to minimize their disparity.

- lower disparity indicates more similar configurations;
- larger Procrustes `r` indicates greater concordance;
- sample connector length shows sample-level disagreement between modalities;
- significance is estimated by permuting sample correspondence.

The program exports every PCoA eigenvalue, including negative eigenvalues. Negative values indicate that a distance matrix is not perfectly Euclidean. They must be inspected rather than silently ignored. Procrustes uses the requested number of available positive axes.

### Multiple comparisons

With three modalities, there are three pairwise tests. With more modalities, the number grows quickly. The summary reports Benjamini–Hochberg-adjusted q-values separately for Mantel and Procrustes permutation p-values.

## Outputs

```text
sample_alignment_audit.tsv
common_samples.tsv
distance_concordance_summary.tsv

<modality>.aligned.distance.tsv
<modality>.pcoa_eigenvalues.tsv
<modality>.distance_heatmap.png

<modality1>_vs_<modality2>.distance_scatter.png
<modality1>_vs_<modality2>.procrustes.png
```

`sample_alignment_audit.tsv` is a required result, not a disposable log. It records the input sample count, common sample count, and excluded IDs for every modality.

## Interpretation

A positive image–genotype concordance means genetically similar samples tend to have similar image embeddings. It does not demonstrate that genomic differences directly caused visible morphology.

Possible confounders include:

- taxonomic group;
- population structure;
- growth environment;
- developmental stage;
- camera, background, scale, and viewpoint;
- sequencing or expression batch;
- sample-label mismatch;
- multiple images or RNA replicates per accession.

RNA–image concordance may reflect biological morphology–expression coupling, but it can also reflect shared species labels, tissue differences, or experimental batches. Stratified analyses or covariate-aware methods may be needed after the initial distance comparison.

## Recommended analysis sequence

```text
1. Freeze canonical sample manifest
2. QC each modality independently
3. Aggregate replicates at the declared biological level
4. Build one validated distance matrix per modality
5. Inspect distance distributions and PCoA eigenvalues
6. Align the global sample intersection
7. Run Mantel and Procrustes comparisons
8. Inspect sample-level discordance and confounders
9. Only then create optional trees or publication layouts
```
