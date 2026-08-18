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
import matplotlib.pyplot as plt

# ==========================================
# 1. CLI Arguments & Setup
# ==========================================
def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Visualize dominant spatial feature activations from an "
            "ImageNet-pretrained ResNet50 using Eigen-CAM for plant images. "
            "(Note: Eigen-CAM visualizes dominant deep feature patterns rather than class-specific predictions.)"
        )
    )
    parser.add_argument("--image_folder", type=str, default="output", help="Path to the directory containing plant images.")
    parser.add_argument("--meta_path", type=str, default=None, help="Path to metadata file (CSV or TSV). Optional.")
    parser.add_argument("--id_col", type=str, default="sample_id", help="Column name for sample ID in metadata.")
    parser.add_argument("--group_col", type=str, default="call", help="Column name for grouping/phenotype call.")
    parser.add_argument("--id_regex", type=str, default=None, help="Regex to extract sample ID from filename (e.g., '^(kale\\d+)').")
    parser.add_argument("--num_samples", type=int, default=5, help="Number of samples to visualize.")
    parser.add_argument("--sampling", type=str, choices=["first", "random"], default="first", help="Method to select samples.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility.")
    parser.add_argument("--device", type=str, choices=["auto", "cpu", "cuda", "mps"], default="auto", help="Computing device.")
    parser.add_argument("--output", type=str, default=None, help="Path to save the generated figure instead of displaying it.")
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


def load_metadata(meta_path, id_col, group_col):
    if meta_path is None:
        return {}
    
    if not os.path.exists(meta_path):
        raise FileNotFoundError(f"Metadata file not found: {meta_path}")
    
    sep = '\t' if meta_path.endswith('.tsv') else ','
    meta_df = pd.read_csv(meta_path, sep=sep)
    
    required_cols = {id_col, group_col}
    missing_cols = required_cols - set(meta_df.columns)
    if missing_cols:
        raise ValueError(f"Missing metadata columns: {', '.join(sorted(missing_cols))}")
    
    # Check for duplicate sample IDs to prevent silent bugs
    if meta_df[id_col].duplicated().any():
        duplicated = meta_df.loc[meta_df[id_col].duplicated(), id_col].unique()
        raise ValueError(
            f"Duplicate sample IDs found in metadata: {', '.join(map(str, duplicated[:10]))}"
        )
    
    meta_df[id_col] = meta_df[id_col].astype(str)
    return dict(zip(meta_df[id_col], meta_df[group_col]))


def get_preprocess():
    return transforms.Compose([
        transforms.Resize((224, 224)),
        transforms.ToTensor(),
        transforms.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225]),
    ])


# ==========================================
# 3. Model & Eigen-CAM Core
# ==========================================
def load_resnet50_model(device):
    print(f"Loading ImageNet-pretrained ResNet50 model to {device}...")
    model = models.resnet50(weights=models.ResNet50_Weights.DEFAULT)
    model.eval()
    model.to(device)
    
    activations = []
    def hook_fn(module, input, output):
        activations.append(output)
        
    target_layer = model.layer4[-1]
    hook_handle = target_layer.register_forward_hook(hook_fn)
    return model, activations, hook_handle


def generate_eigencam(img_path, model, activations, device, preprocess):
    img_pil = Image.open(img_path).convert('RGB')
    img_tensor = preprocess(img_pil).unsqueeze(0).to(device)
    
    activations.clear()
    with torch.no_grad():
        model(img_tensor)
        
    feature_map = activations[0][0].detach().cpu()  # (C, H, W) -> (2048, 7, 7)
    c, h, w = feature_map.shape
    
    # Reshape feature activations to (H*W, C)
    spatial_features = feature_map.permute(1, 2, 0).reshape(-1, c)
    
    try:
        # Eigen-CAM: project activations onto the first principal component
        _, _, Vh = torch.linalg.svd(spatial_features, full_matrices=False)
        principal_component = Vh[0]  # (C,)
        cam = torch.matmul(spatial_features, principal_component).reshape(h, w)
    except RuntimeError as e:
        print(f"Warning: SVD failed for {img_path}: {e}. Falling back to channel mean.")
        cam = torch.mean(feature_map, dim=0)
        
    # Retain positive projections and normalize
    cam = torch.relu(cam)
    cam_max = cam.max()
    if cam_max > 0:
        cam = cam / cam_max
    cam = cam.numpy()
    
    cam_resized = cv2.resize(cam, (224, 224))
    heatmap = cv2.applyColorMap(np.uint8(255 * cam_resized), cv2.COLORMAP_JET)
    heatmap = cv2.cvtColor(heatmap, cv2.COLOR_BGR2RGB)
    
    original_resized = np.array(img_pil.resize((224, 224)))
    overlay = cv2.addWeighted(original_resized, 0.6, heatmap, 0.4, 0)
    
    return original_resized, overlay


# ==========================================
# 4. Main Execution Flow
# ==========================================
def main():
    args = parse_args()
    
    if args.num_samples < 1:
        raise ValueError("--num_samples must be at least 1.")
        
    device = select_device(args.device)
    
    # 1. Load metadata & mapping
    id_to_group = load_metadata(args.meta_path, args.id_col, args.group_col)
    
    if not os.path.exists(args.image_folder):
        raise FileNotFoundError(f"Image directory not found: {args.image_folder}")
        
    valid_extensions = ('.png', '.jpg', '.jpeg', '.bmp', '.tiff')
    files = sorted([f for f in os.listdir(args.image_folder) if f.lower().endswith(valid_extensions)])
    
    if not files:
        print(f"No valid images found in {args.image_folder}.")
        return
        
    # Validate regex matching before sampling
    if args.id_regex is not None:
        matched_count = sum(extract_sample_id(f, args.id_regex) is not None for f in files)
        if matched_count == 0:
            raise ValueError(f"--id_regex '{args.id_regex}' did not match any image filenames in {args.image_folder}.")
        print(f"Regex validation passed: {matched_count}/{len(files)} files matched.")
        
    # 2. Sampling strategy with bound check
    if args.sampling == "random":
        import random
        random.seed(args.seed)
        random.shuffle(files)
    
    num_samples = min(args.num_samples, len(files))
    if args.num_samples > len(files):
        print(
            f"Requested {args.num_samples} images, "
            f"but only {len(files)} are available. "
            f"Using all available images."
        )
    selected_files = files[:num_samples]
    
    # 3. Load Model, Hook, & Preprocess
    model, activations, hook_handle = load_resnet50_model(device)
    preprocess = get_preprocess()
    
    try:
        # 4. Process images
        rows = len(selected_files)
        fig, axes = plt.subplots(rows, 2, figsize=(8, 4 * rows))
        if rows == 1:
            axes = np.array([axes])
            
        fig.suptitle("ResNet50 Eigen-CAM Feature Activation", fontsize=16)
        
        for i, filename in enumerate(selected_files):
            path = os.path.join(args.image_folder, filename)
            sample_id = extract_sample_id(filename, args.id_regex)
            
            group = "Unknown"
            if sample_id and sample_id in id_to_group:
                group = str(id_to_group[sample_id])
                
            orig, attention = generate_eigencam(path, model, activations, device, preprocess)
            
            # Left: Original Image
            axes[i, 0].imshow(orig)
            axes[i, 0].set_title(f"ID: {sample_id} | Group: {group}\nFile: {filename}", fontsize=10, fontweight='bold', loc='left')
            axes[i, 0].axis('off')
            
            # Right: Eigen-CAM
            axes[i, 1].imshow(attention)
            axes[i, 1].set_title("ResNet50 Eigen-CAM", fontsize=11)
            axes[i, 1].axis('off')
            
        plt.tight_layout()
        
        if args.output:
            os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
            plt.savefig(args.output, dpi=300, bbox_inches="tight")
            print(f"Visualization successfully saved to: {args.output}")
        else:
            plt.show()
            
        plt.close(fig)
        
    finally:
        # Clean up hook handle properly
        hook_handle.remove()

if __name__ == "__main__":
    main()