import os
import re
import argparse
import pandas as pd
import numpy as np
import cv2
from PIL import Image
import torch
import torchvision.models as models
import torchvision.transforms as transforms
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

# ==========================================
# 1. CLI Arguments & Setup
# ==========================================
def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Perform multi-scale unsupervised spatial feature clustering (K-Means) across "
            "ResNet50 layers (Layer 1 to Layer 4) for plant images. "
            "(Note: Cluster identities are image- and layer-specific and not comparable across panels.)"
        )
    )
    parser.add_argument("--image_folder", type=str, default="output", help="Path to the directory containing plant images.")
    parser.add_argument("--meta_path", type=str, default=None, help="Path to annotation file (CSV or TSV). Optional.")
    parser.add_argument("--id_col", type=str, default="sample_id", help="Column name for sample ID in annotation file.")
    parser.add_argument("--annotation_col", type=str, default="call", help="Column name for display annotation (e.g., genotype, species).")
    parser.add_argument("--id_regex", type=str, default=None, help="Regex to extract sample ID from filename.")
    parser.add_argument("--n_clusters", type=int, default=4, help="Number of clusters for spatial feature clustering.")
    parser.add_argument("--num_samples", type=int, default=3, help="Number of samples to visualize.")
    parser.add_argument("--sampling", type=str, choices=["first", "random"], default="first", help="Method to select samples.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility.")
    parser.add_argument("--device", type=str, choices=["auto", "cpu", "cuda", "mps"], default="auto", help="Computing device.")
    parser.add_argument("--output", type=str, default=None, help="Path to save the generated figure.")
    return parser.parse_args()


# ==========================================
# 2. Helper Functions
# ==========================================
def select_device(device_arg):
    if device_arg == "cpu":
        return torch.device("cpu")
    elif device_arg == "cuda":
        if not torch.cuda.is_available():
            raise RuntimeError("CUDA is requested but not available.")
        return torch.device("cuda")
    elif device_arg == "mps":
        if not torch.backends.mps.is_available():
            raise RuntimeError("MPS (Apple Silicon) is requested but not available.")
        return torch.device("mps")
    
    if torch.backends.mps.is_available():
        return torch.device("mps")
    elif torch.cuda.is_available():
        return torch.device("cuda")
    else:
        return torch.device("cpu")


def extract_sample_id(filename, id_regex=None):
    stem = os.path.splitext(os.path.basename(filename))[0]
    if id_regex is None:
        return stem
    try:
        match = re.search(id_regex, stem)
    except re.error as e:
        raise ValueError(f"Invalid --id_regex '{id_regex}': {e}") from e
    if not match:
        return None
    return match.group(1) if match.groups() else match.group(0)


def load_annotations(meta_path, id_col, annotation_col):
    """Load optional display annotations. Does NOT influence clustering."""
    if meta_path is None:
        return {}
    if not os.path.exists(meta_path):
        raise FileNotFoundError(f"Annotation file not found: {meta_path}")
    sep = '\t' if meta_path.endswith('.tsv') else ','
    meta_df = pd.read_csv(meta_path, sep=sep)
    
    required_cols = {id_col, annotation_col}
    missing_cols = required_cols - set(meta_df.columns)
    if missing_cols:
        raise ValueError(f"Missing annotation columns: {', '.join(sorted(missing_cols))}")
        
    meta_df[id_col] = meta_df[id_col].astype(str)
    return dict(zip(meta_df[id_col], meta_df[annotation_col]))


def get_preprocess():
    return transforms.Compose([
        transforms.Resize((224, 224)),
        transforms.ToTensor(),
        transforms.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225]),
    ])


# ==========================================
# 3. Model & Multi-Layer Extraction Core
# ==========================================
def load_resnet50(device):
    print(f"Loading ImageNet-pretrained ResNet50 to {device}...")
    model = models.resnet50(weights=models.ResNet50_Weights.DEFAULT)
    model.eval()
    model.to(device)
    return model


def extract_multilayer_features(model, x):
    """Single forward pass to collect feature maps from all 4 ResNet stages."""
    features = {}
    
    x = model.conv1(x)
    x = model.bn1(x)
    x = model.relu(x)
    x = model.maxpool(x)

    x = model.layer1(x)
    features["Layer 1 (56x56)"] = x

    x = model.layer2(x)
    features["Layer 2 (28x28)"] = x

    x = model.layer3(x)
    features["Layer 3 (14x14)"] = x

    x = model.layer4(x)
    features["Layer 4 (7x7)"] = x

    return features


def compute_clustering_map(features_tensor, n_clusters=4):
    features = features_tensor.squeeze(0).cpu()
    c, h, w = features.shape[0], features.shape[1], features.shape[2]
    
    features_flat = features.view(c, -1).permute(1, 0).numpy()
    
    # Ensure n_clusters does not exceed spatial points for shallower layers
    actual_k = min(n_clusters, h * w)
    
    kmeans = KMeans(n_clusters=actual_k, random_state=42, n_init=10)
    labels = kmeans.fit_predict(features_flat)

    cluster_map = labels.reshape(h, w)
    cluster_map_resized = cv2.resize(
        cluster_map.astype(np.int32), 
        (224, 224), 
        interpolation=cv2.INTER_NEAREST
    )
    return cluster_map_resized


# ==========================================
# 4. Main Execution Flow
# ==========================================
def main():
    args = parse_args()
    if args.n_clusters < 2:
        raise ValueError("--n_clusters must be at least 2.")
    if args.num_samples < 1:
        raise ValueError("--num_samples must be at least 1.")
        
    device = select_device(args.device)
    id_to_annotation = load_annotations(args.meta_path, args.id_col, args.annotation_col)
    
    if not os.path.exists(args.image_folder):
        raise FileNotFoundError(f"Image directory not found: {args.image_folder}")
        
    valid_extensions = ('.png', '.jpg', '.jpeg', '.bmp', '.tiff')
    files = sorted([f for f in os.listdir(args.image_folder) if f.lower().endswith(valid_extensions)])
    
    if not files:
        print(f"No valid images found in {args.image_folder}.")
        return
        
    # Sampling strategy
    if args.sampling == "random":
        import random
        random.seed(args.seed)
        random.shuffle(files)
        
    num_samples = min(args.num_samples, len(files))
    selected_files = files[:num_samples]

    model = load_resnet50(device)
    preprocess = get_preprocess()

    # Layout: 1 column for original image + 4 columns for layers (Layer 1 to 4)
    rows = len(selected_files)
    fig, axes = plt.subplots(rows, 5, figsize=(18, 4 * rows))
    if rows == 1:
        axes = np.array([axes])
        
    fig.suptitle(f"Multi-Scale ResNet50 Spatial Feature Clustering (K={args.n_clusters})", fontsize=16)
    cmap_name = "tab20" if args.n_clusters <= 20 else "nipy_spectral"

    for i, filename in enumerate(selected_files):
        path = os.path.join(args.image_folder, filename)
        sample_id = extract_sample_id(filename, args.id_regex)
        
        annotation = "N/A"
        if sample_id and sample_id in id_to_annotation:
            annotation = str(id_to_annotation[sample_id])
            
        img_pil = Image.open(path).convert('RGB')
        img_tensor = preprocess(img_pil).unsqueeze(0).to(device)

        # 1. Original Image with optional annotation display (Column 0)
        orig_resized = np.array(img_pil.resize((224, 224)))
        axes[i, 0].imshow(orig_resized)
        axes[i, 0].set_title(f"ID: {sample_id}\nAnnotation: {annotation}\nFile: {filename}", fontsize=9, fontweight='bold', loc='left')
        axes[i, 0].axis('off')

        # 2. Extract multi-layer features in a SINGLE forward pass
        with torch.no_grad():
            feature_maps = extract_multilayer_features(model, img_tensor)

        # 3. Compute and plot cluster maps for Layer 1 ~ Layer 4 (Columns 1-4)
        for j, (layer_name, feat_map) in enumerate(feature_maps.items(), start=1):
            cluster_map_resized = compute_clustering_map(feat_map, n_clusters=args.n_clusters)
            
            axes[i, j].imshow(cluster_map_resized, cmap=cmap_name, vmin=0, vmax=args.n_clusters - 1)
            axes[i, j].set_title(layer_name, fontsize=11)
            axes[i, j].axis('off')

    plt.tight_layout()
    
    if args.output:
        os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
        plt.savefig(args.output, dpi=300, bbox_inches="tight")
        print(f"Multi-scale clustering visualization successfully saved to: {args.output}")
    else:
        plt.show()
        
    plt.close(fig)

if __name__ == "__main__":
    main()