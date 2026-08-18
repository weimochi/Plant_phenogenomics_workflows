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
            "Perform feature-space occlusion sensitivity analysis on plant images "
            "using an ImageNet-pretrained ResNet50. Spatial importance is quantified "
            "by the change in the global deep feature representation after local occlusion."
        )
    )
    parser.add_argument("--image_folder", type=str, required=True, help="Path to the directory containing plant images.")
    parser.add_argument("--meta_path", type=str, default=None, help="Path to annotation file (CSV or TSV). Optional.")
    parser.add_argument("--id_col", type=str, default="sample_id", help="Column name for sample ID.")
    parser.add_argument("--annotation_col", type=str, default="call", help="Column name for display annotation.")
    parser.add_argument("--id_regex", type=str, default=None, help="Regex to extract sample ID from filename.")
    parser.add_argument("--patch_size", type=int, default=16, help="Size of the sliding occlusion patch (default: 16).")
    parser.add_argument("--stride", type=int, default=8, help="Stride for the sliding occlusion patch (default: 8).")
    parser.add_argument("--batch_size", type=int, default=32, help="Batch size for streaming inference.")
    parser.add_argument("--num_samples", type=int, default=3, help="Number of samples to visualize.")
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
    if meta_path is None:
        return {}
    if not os.path.exists(meta_path):
        raise FileNotFoundError(f"Annotation file not found: {meta_path}")
    sep = '\t' if meta_path.endswith('.tsv') else ','
    meta_df = pd.read_csv(meta_path, sep=sep)
    meta_df[id_col] = meta_df[id_col].astype(str)
    return dict(zip(meta_df[id_col], meta_df[annotation_col]))


def get_preprocess():
    return transforms.Compose([
        transforms.Resize((224, 224)),
        transforms.ToTensor(),
        transforms.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225]),
    ])


def get_window_starts(length, patch_size, stride):
    """Ensure full spatial coverage including the exact image boundaries."""
    starts = list(range(0, length - patch_size + 1, stride))
    last_start = length - patch_size
    if starts[-1] != last_start:
        starts.append(last_start)
    return starts


# ==========================================
# 3. Feature-Space Occlusion Sensitivity Core
# ==========================================
def compute_occlusion_sensitivity(img_pil, model, device, preprocess, patch_size=16, stride=8, batch_size=32):
    model.eval()
    img_resized = img_pil.resize((224, 224))
    orig_np = np.array(img_resized).astype(np.float32) / 255.0
    
    # 1. Baseline feature representation without occlusion
    orig_tensor = preprocess(img_resized).unsqueeze(0).to(device)
    with torch.no_grad():
        baseline_out = model(orig_tensor)
        baseline_feature = torch.flatten(baseline_out, 1)  # (1, 2048)
        
    h, w = 224, 224
    sensitivity_map = np.zeros((h, w), dtype=np.float32)
    count_map = np.zeros((h, w), dtype=np.float32)

    # Fill occluded regions with the ImageNet RGB mean.
    # After normalization, the occluded region is approximately zero-valued.
    patch_color = [0.485, 0.456, 0.406]

    y_positions = get_window_starts(h, patch_size, stride)
    x_positions = get_window_starts(w, patch_size, stride)

    # 2. Streaming batch inference to prevent RAM exhaustion
    occluded_batch = []
    positions = []

    for y in y_positions:
        for x in x_positions:
            occluded_np = orig_np.copy()
            for c_idx, c_val in enumerate(patch_color):
                occluded_np[y:y+patch_size, x:x+patch_size, c_idx] = c_val
                
            occluded_pil = Image.fromarray((occluded_np * 255).astype(np.uint8))
            occluded_batch.append(preprocess(occluded_pil))
            positions.append((y, x))

            # When batch size is reached, process and flush
            if len(occluded_batch) >= batch_size:
                tensor_batch = torch.stack(occluded_batch).to(device)
                with torch.no_grad():
                    outputs = model(tensor_batch)
                    output_features = torch.flatten(outputs, 1)  # (B, 2048)
                    
                    # Compute L2 vector norm of feature difference
                    diffs = torch.linalg.vector_norm(
                        output_features - baseline_feature,
                        ord=2,
                        dim=1
                    ).cpu().numpy()
                
                for idx, diff in enumerate(diffs):
                    py, px = positions[idx]
                    sensitivity_map[py:py+patch_size, px:px+patch_size] += diff
                    count_map[py:py+patch_size, px:px+patch_size] += 1.0

                occluded_batch = []
                positions = []

    # Process remaining items in the last batch
    if occluded_batch:
        tensor_batch = torch.stack(occluded_batch).to(device)
        with torch.no_grad():
            outputs = model(tensor_batch)
            output_features = torch.flatten(outputs, 1)
            diffs = torch.linalg.vector_norm(
                output_features - baseline_feature,
                ord=2,
                dim=1
            ).cpu().numpy()
        
        for idx, diff in enumerate(diffs):
            py, px = positions[idx]
            sensitivity_map[py:py+patch_size, px:px+patch_size] += diff
            count_map[py:py+patch_size, px:px+patch_size] += 1.0

    # Avoid division by zero
    count_map[count_map == 0] = 1.0
    sensitivity_map /= count_map

    # Normalize map to [0, 1]
    if sensitivity_map.max() > sensitivity_map.min():
        sensitivity_map = (sensitivity_map - sensitivity_map.min()) / (sensitivity_map.max() - sensitivity_map.min())

    return orig_np, sensitivity_map


# ==========================================
# 4. Main Execution Flow
# ==========================================
def main():
    args = parse_args()
    
    # Parameter validations
    if args.patch_size < 1 or args.patch_size > 224:
        raise ValueError("--patch_size must be between 1 and 224.")
    if args.stride < 1:
        raise ValueError("--stride must be at least 1.")
    if args.num_samples < 1:
        raise ValueError("--num_samples must be at least 1.")
    if args.batch_size < 1:
        raise ValueError("--batch_size must be at least 1.")

    device = select_device(args.device)
    id_to_annotation = load_annotations(args.meta_path, args.id_col, args.annotation_col)
    
    if not os.path.exists(args.image_folder):
        raise FileNotFoundError(f"Image directory not found: {args.image_folder}")
        
    valid_extensions = ('.png', '.jpg', '.jpeg', '.bmp', '.tiff')
    files = sorted([f for f in os.listdir(args.image_folder) if f.lower().endswith(valid_extensions)])
    
    if not files:
        print(f"No valid images found in {args.image_folder}.")
        return
        
    num_samples = min(args.num_samples, len(files))
    selected_files = files[:num_samples]

    print(f"Loading ResNet50 model to {device}...")
    model = models.resnet50(weights=models.ResNet50_Weights.DEFAULT)
    # Stop at avgpool to get global feature representation (2048-D)
    model = torch.nn.Sequential(*list(model.children())[:-1])
    model.to(device)
    preprocess = get_preprocess()

    rows = len(selected_files)
    fig, axes = plt.subplots(rows, 2, figsize=(8, 4 * rows))
    if rows == 1:
        axes = np.array([axes])
        
    fig.suptitle("Feature-Space Occlusion Sensitivity Analysis", fontsize=16)

    for i, filename in enumerate(selected_files):
        path = os.path.join(args.image_folder, filename)
        sample_id = extract_sample_id(filename, args.id_regex)
        
        annotation = id_to_annotation.get(str(sample_id), "N/A") if sample_id else "N/A"
        
        img_pil = Image.open(path).convert('RGB')
        orig_np, sensitivity_map = compute_occlusion_sensitivity(
            img_pil, model, device, preprocess, 
            patch_size=args.patch_size, stride=args.stride, batch_size=args.batch_size
        )

        # Left: Original Image
        axes[i, 0].imshow(orig_np)
        axes[i, 0].set_title(f"ID: {sample_id} | Ann: {annotation}\nFile: {filename}", fontsize=9, fontweight='bold', loc='left')
        axes[i, 0].axis('off')

        # Right: Feature-Space Occlusion Sensitivity Heatmap
        axes[i, 1].imshow(orig_np)
        axes[i, 1].imshow(sensitivity_map, cmap='jet', alpha=0.5, vmin=0, vmax=1)
        axes[i, 1].set_title("Feature-Space Occlusion Sensitivity", fontsize=11)
        axes[i, 1].axis('off')

    plt.tight_layout()
    
    if args.output:
        os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
        plt.savefig(args.output, dpi=300, bbox_inches="tight")
        print(f"Feature-space occlusion sensitivity visualization saved to: {args.output}")
    else:
        plt.show()
        
    plt.close(fig)

if __name__ == "__main__":
    main()