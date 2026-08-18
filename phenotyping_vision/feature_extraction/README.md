# Plant Phenotyping Vision: Pretrained Feature Interpretation & Exploration Toolkit

A set of unsupervised feature interpretation and exploration tools designed for plant phenotyping images. Utilizing an ImageNet-pretrained **ResNet50** (without fine-tuning), this toolkit provides multi-angle perspectives to inspect spatial structures and feature dependencies captured by convolutional neural networks on plant imagery.

This toolkit includes three complementary command-line scripts with unified CLI conventions:
1. **`visualize_resnet_eigencam.py`**: Visualizes dominant activation patterns from the final convolutional layer.
2. **`visualize_resnet_spatial_clustering.py`**: Explores multi-scale spatial partitions of deep feature representations across ResNet50 layers using K-Means.
3. **`visualize_resnet_occlusion_sensitivity.py`**: Evaluates region-based feature importance through a perturbation-based occlusion method.

---

## Tool Overview & Core Principles

### 1. Eigen-CAM (`visualize_resnet_eigencam.py`)
* **Principle**: Applies Principal Component Analysis (PCA / extracting the first principal component) directly to the activation maps of the final ResNet50 convolutional layer (`layer4`) without requiring class logits.
* **Scientific Meaning**: Examines which spatial regions are most strongly associated with the dominant component of the final convolutional feature representation.
* **Use Case**: Inspecting whether dominant deep feature activations are spatially associated with plant structures rather than background elements, pots, or shadows.

### 2. Multi-Scale Spatial Clustering (`visualize_resnet_spatial_clustering.py`)
* **Principle**: Chains all four primary ResNet50 convolutional stages (Layers 1 to 4) within a single forward pass. Performs independent K-Means spatial feature clustering on spatial locations at each layer (e.g., $14 \times 14 = 196$ locations for Layer 3).
* **Scientific Meaning**: Explores multi-scale spatial partitions of deep feature representations across ResNet50 layers using K-Means, illustrating how feature granularity transitions from shallow textures to deep regional contours.
* **Use Case**: Exploring spatial feature groupings across layers to understand how unsupervised representations segment plant images at different receptive field scales.

### 3. Feature-Space Occlusion Sensitivity (`visualize_resnet_occlusion_sensitivity.py`)
* **Principle**: A perturbation-based method that systematically occludes local regions using the ImageNet RGB mean and measures the resulting L2-distance change in the global 2048-dimensional feature representation.
* **Scientific Meaning**: Measures how strongly the global deep feature representation changes when a specific local image region is occluded.
* **Use Case**: Identifying which local image regions have the highest impact on the global feature vector when perturbed.

---

## Important Caveats & Interpretational Boundaries

* **Domain Mismatch**: Models use ImageNet pre-trained weights without fine-tuning. Extracted features are based on generic visual priors and do not necessarily align directly with specific biological phenotypes.
* **Fixed-K Clustering**: K-Means always partitions spatial features into the user-specified number of clusters. The resulting cluster maps should therefore be interpreted as exploratory feature-space partitions, not as evidence that an image inherently contains exactly $K$ biological structures.
* **Cluster Semantics**: Cluster labels do not automatically correspond to biological structures such as veins, lamina, edges, or background. Such interpretations require independent validation.
* **Cluster Identity**: K-Means is fitted independently for each image and layer. Cluster IDs therefore have no shared semantic meaning across images or layers.
* **Occlusion Map Scaling**: Occlusion sensitivity maps are min-max normalized independently for visualization. Their absolute intensity values should therefore not be quantitatively compared across images.

---

## Installation & Requirements

Ensure your Python environment has the following core packages installed:

```bash
pip install torch torchvision numpy pandas opencv-python pillow scikit-learn matplotlib
```
## Usage
All scripts share a standardized set of command-line arguments (CLI) and utilize a display-only annotation isolation design (metadata is only used for chart title annotations and never interferes with unsupervised clustering or feature calculations).

### 1. Run Eigen-CAM Visualization
```bash
python visualize_resnet_eigencam.py \
  --image_folder output \
  --meta_path metadata.tsv \
  --id_col sample_id \
  --annotation_col call \
  --id_regex '^(kale\d+)' \
  --num_samples 3 \
  --sampling random \
  --seed 42 \
  --output results/eigencam_result.png
```
### 2. Run Multi-Scale Spatial Clustering (Layers 1–4)
```bash
python visualize_resnet_spatial_clustering.py \
  --image_folder output \
  --meta_path metadata.tsv \
  --id_col sample_id \
  --annotation_col call \
  --id_regex '^(kale\d+)' \
  --n_clusters 4 \
  --num_samples 3 \
  --sampling random \
  --seed 42 \
  --output results/multiscale_clustering.png
  ```
### 3. Run Feature-Space Occlusion Sensitivity Analysis
```bash
python visualize_resnet_occlusion_sensitivity.py \
  --image_folder output \
  --meta_path metadata.tsv \
  --id_col sample_id \
  --annotation_col call \
  --id_regex '^(kale\d+)' \
  --patch_size 16 \
  --stride 8 \
  --batch_size 32 \
  --num_samples 3 \
  --output results/occlusion_sensitivity.png
```
## Common Arguments

* `--image_folder` : Path to the directory containing plant images (default: `output`).
* `--meta_path` : Path to an optional annotation file (`.csv` or `.tsv`). Skipped if omitted.
* `--id_col` : Column name for sample IDs in the metadata (default: `sample_id`).
* `--annotation_col` : Column name for display annotations such as genotype or cultivar (default: `call`).
* `--id_regex` : Regular expression to extract the sample ID from filenames (e.g., `'^(kale\d+)'`).
* `--num_samples` : Number of samples to visualize in the figure layout (default: `3`).
* `--sampling` : Method to select samples (`first` or `random`, supported in clustering and Eigen-CAM).
* `--seed` : Random seed for reproducibility and sampling/K-Means initialization (default: `42`).
* `--batch_size` : Batch size for streaming inference in occlusion sensitivity (default: `32`).
* `--device` : Computation device (`auto`, `cpu`, `cuda`, `mps`).
* `--output` : Path to save the generated figure. Displays interactively in a window if omitted.