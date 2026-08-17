import gradio as gr
import cv2
import numpy as np
import torch
from segment_anything import sam_model_registry, SamPredictor
import os
import time

# ==========================================
# 1. Configuration & Model Loading
# ==========================================
CHECKPOINT_PATH = "sam_vit_h_4b8939.pth"
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
if torch.backends.mps.is_available():
    DEVICE = "mps"
MODEL_TYPE = "vit_h"

print(f"Loading SAM model ({DEVICE})... Please wait")
if not os.path.exists(CHECKPOINT_PATH):
    raise FileNotFoundError(f"Model checkpoint not found: {CHECKPOINT_PATH}. Please download it first!")

sam = sam_model_registry[MODEL_TYPE](checkpoint=CHECKPOINT_PATH)
sam.to(device=DEVICE)
predictor = SamPredictor(sam)
print("Model loaded successfully!")

# Global State (Maintained for single-user research session)
state = {
    "image": None,
    "original_stem": None,
    "points": [],
    "labels": [],
    "sam_raw_mask": None,
    "final_mask": None
}


# ==========================================
# 2. Core Logic & Processing
# ==========================================

def reset_state():
    state["points"] = []
    state["labels"] = []
    state["sam_raw_mask"] = None
    state["final_mask"] = None


def process_image_upload(filepath):
    if filepath is None: 
        return None

    filename = os.path.basename(filepath)
    file_stem = os.path.splitext(filename)[0]
    state["original_stem"] = file_stem

    image = cv2.imread(filepath)
    image = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
    state["image"] = image

    reset_state()
    predictor.set_image(image)

    # 🛠️ Fixed: Exactly 2 outputs to match Gradio component binding
    return image, f"Loaded: {filename}"


def get_preset_hsv(preset_name):
    """
    Return HSV slider values based on the selected plant type preset (H: 0-180).
    """
    if "Green with Yellow" in preset_name:
        return True, 10, 90, 30, 255, 40, 255
    elif "Strict Green" in preset_name:
        return True, 30, 90, 40, 255, 40, 255
    elif "Dark / Pigmented Plant" in preset_name:
        # Broad HSV range for dark/saturated plants (relying on S & V instead of a tight hue wrap)
        return True, 0, 180, 20, 255, 20, 255
    else:  # "General / No Filter"
        return False, 0, 180, 0, 255, 0, 255


def apply_filters(sam_mask, use_hsv, h_min, h_max, s_min, s_max, v_min, v_max, glue_strength, opening_strength):
    """Filtering and Post-processing with Conservative Morphology"""
    if sam_mask is None or state["image"] is None:
        return None

    image = state["image"]
    final_mask = sam_mask.astype(np.uint8) * 255

    # 1. HSV Filtering
    if use_hsv:
        hsv = cv2.cvtColor(image, cv2.COLOR_RGB2HSV)

        if h_min > h_max:
            lower1 = np.array([h_min, s_min, v_min])
            upper1 = np.array([180, s_max, v_max])
            lower2 = np.array([0, s_min, v_min])
            upper2 = np.array([h_max, s_max, v_max])
            mask1 = cv2.inRange(hsv, lower1, upper1)
            mask2 = cv2.inRange(hsv, lower2, upper2)
            hsv_mask = cv2.bitwise_or(mask1, mask2)
        else:
            lower = np.array([h_min, s_min, v_min])
            upper = np.array([h_max, s_max, v_max])
            hsv_mask = cv2.inRange(hsv, lower, upper)

        final_mask = np.logical_and(sam_mask, hsv_mask > 0).astype(np.uint8) * 255

    # 2. Morphological Closing (Optional repair tool, default=0 to prevent unwanted merging)
    if glue_strength > 0:
        k_size = int(glue_strength * 2 + 1)
        kernel = np.ones((k_size, k_size), np.uint8)
        final_mask = cv2.morphologyEx(final_mask, cv2.MORPH_CLOSE, kernel)

    # 3. Morphological Opening (Optional / Conservative Denoising)
    if opening_strength > 0:
        k_open = int(opening_strength * 2 + 1)
        open_kernel = np.ones((k_open, k_open), np.uint8)
        final_mask = cv2.morphologyEx(final_mask, cv2.MORPH_OPEN, open_kernel)

    # 4. Component Selection Strategy
    # 🛠️ Fixed: Retain ALL components that contain positive clicks to support disconnected plant parts
    num_labels, labels, stats, _ = cv2.connectedComponentsWithStats(final_mask, connectivity=8)
    if num_labels > 1:
        positive_points = [p for p, l in zip(state["points"], state["labels"]) if l == 1]
        
        if positive_points:
            target_labels = set()
            for px, py in positive_points:
                ix, iy = int(px), int(py)
                if 0 <= iy < labels.shape[0] and 0 <= ix < labels.shape[1]:
                    lbl = labels[iy, ix]
                    if lbl > 0:
                        target_labels.add(lbl)
            
            if target_labels:
                selected_mask = np.zeros_like(final_mask)
                for lbl in target_labels:
                    selected_mask[labels == lbl] = 255
                final_mask = selected_mask
            else:
                # Fallback to largest if clicks don't hit any labeled component
                largest_label = 1 + np.argmax(stats[1:, cv2.CC_STAT_AREA])
                final_mask = np.zeros_like(final_mask)
                final_mask[labels == largest_label] = 255
        else:
            # Fallback to largest if no positive points exist yet
            largest_label = 1 + np.argmax(stats[1:, cv2.CC_STAT_AREA])
            final_mask = np.zeros_like(final_mask)
            final_mask[labels == largest_label] = 255

    state["final_mask"] = final_mask
    return final_mask


def update_view(mask):
    if state["image"] is None: 
        return None, None
    mask_view = mask if mask is not None else np.zeros_like(state["image"][:, :, 0])

    vis_img = state["image"].copy()
    if mask is not None:
        overlay = state["image"].copy()
        overlay[mask == 255] = [255, 0, 0]  # Red mask overlay
        vis_img = cv2.addWeighted(state["image"], 0.7, overlay, 0.3, 0)

    for p, l in zip(state["points"], state["labels"]):
        color = (0, 255, 0) if l == 1 else (0, 0, 255)
        cv2.circle(vis_img, (int(p[0]), int(p[1])), 8, color, -1)
        cv2.circle(vis_img, (int(p[0]), int(p[1])), 8, (255, 255, 255), 2)

    return mask_view, vis_img


def handle_click(evt: gr.SelectData, mode, use_hsv, h_min, h_max, s_min, s_max, v_min, v_max, glue, opening):
    if state["image"] is None: 
        return None, None

    label = 1 if "Add Point" in mode else 0
    state["points"].append([evt.index[0], evt.index[1]])
    state["labels"].append(label)

    masks, scores, _ = predictor.predict(
        point_coords=np.array(state["points"]),
        point_labels=np.array(state["labels"]),
        multimask_output=True
    )

    best_idx = np.argmax(scores)
    state["sam_raw_mask"] = masks[best_idx]

    processed_mask = apply_filters(state["sam_raw_mask"], use_hsv, h_min, h_max, s_min, s_max, v_min, v_max, glue, opening)
    return update_view(processed_mask)


def handle_slider_change(use_hsv, h_min, h_max, s_min, s_max, v_min, v_max, glue, opening):
    if state["sam_raw_mask"] is None: 
        return None, None
    processed_mask = apply_filters(state["sam_raw_mask"], use_hsv, h_min, h_max, s_min, s_max, v_min, v_max, glue, opening)
    return update_view(processed_mask)


def undo_last_point(use_hsv, h_min, h_max, s_min, s_max, v_min, v_max, glue, opening):
    if not state["points"]: 
        return None, None
    state["points"].pop()
    state["labels"].pop()

    if not state["points"]:
        state["sam_raw_mask"] = None
        state["final_mask"] = None
        return update_view(None)

    masks, scores, _ = predictor.predict(
        point_coords=np.array(state["points"]),
        point_labels=np.array(state["labels"]),
        multimask_output=True
    )
    state["sam_raw_mask"] = masks[np.argmax(scores)]
    processed_mask = apply_filters(state["sam_raw_mask"], use_hsv, h_min, h_max, s_min, s_max, v_min, v_max, glue, opening)
    return update_view(processed_mask)


def save_file(target_folder):
    if state["final_mask"] is None: 
        return None, "No image to save"
    mask = state["final_mask"]
    image = state["image"]
    y_indices, x_indices = np.where(mask > 0)
    if len(y_indices) == 0: 
        return None, "Mask is empty"

    padding = 20
    h, w = image.shape[:2]
    
    # 🛠️ Fixed Boundary Slicing: inclusive upper bound (+1)
    x_min, x_max = max(0, np.min(x_indices) - padding), min(w, np.max(x_indices) + padding + 1)
    y_min, y_max = max(0, np.min(y_indices) - padding), min(h, np.max(y_indices) + padding + 1)

    crop_img = image[y_min:y_max, x_min:x_max]
    crop_mask = mask[y_min:y_max, x_min:x_max]
    result = cv2.bitwise_and(crop_img, crop_img, mask=crop_mask)

    ch, cw = result.shape[:2]
    max_dim = max(ch, cw)
    square_img = np.zeros((max_dim, max_dim, 3), dtype=np.uint8)
    x_off, y_off = (max_dim - cw) // 2, (max_dim - ch) // 2
    square_img[y_off:y_off + ch, x_off:x_off + cw] = result

    if not target_folder.strip(): 
        target_folder = "output"
    os.makedirs(target_folder, exist_ok=True)
    stem = state["original_stem"] if state["original_stem"] else f"brassica_{int(time.time())}"
    final_filename = f"{stem}_segmented.png" # Updated suffix for clarity
    save_path = os.path.join(target_folder, final_filename)
    
    ext = final_filename.split('.')[-1]
    success, encoded_img = cv2.imencode(f'.{ext}', cv2.cvtColor(square_img, cv2.COLOR_RGB2BGR))
    if success:
        with open(save_path, mode='wb') as f:
            f.write(encoded_img)
        return square_img, f"✅ Saved: {final_filename}"
    else:
        return None, "❌ Failed to save image."


# ==========================================
# 3. Gradio Interface
# ==========================================

with gr.Blocks(title="Brassica Interactive Segmentation V6.1") as app:
    gr.Markdown("# 🌱 Plant Interactive Segmentation Tool (SAM-Powered)")
    gr.Markdown("Interactive pipeline supporting SAM, HSV constraint, multi-component tracking, and square export.")

    with gr.Row():
        with gr.Column(scale=5):
            input_img = gr.Image(label="1. Upload Image", type="filepath")

            with gr.Row():
                mode_radio = gr.Radio(["➕ Add Point", "➖ Remove Point"], value="➕ Add Point", label="Interaction Mode")
                undo_btn = gr.Button("↩️ Undo")

            gr.Markdown("### 2. Color Strategy Presets")
            color_preset = gr.Radio(
                [
                    "🟢+🟡 Green with Yellow (Natural Growth)",
                    "🟢 Strict Green (Exclude Yellow)",
                    "🟣 Dark / Pigmented Plant (Broad HSV)",
                    "⚪️ General (No HSV Filter)"
                ],
                value="🟢+🟡 Green with Yellow (Natural Growth)",
                label="🍃 Select Plant Phenotype"
            )

            with gr.Group():
                use_hsv_box = gr.Checkbox(value=True, label="Enable HSV Filtering")
                
                # 🛠️ Fixed: Default glue set to 0 (conservative repair tool)
                glue_slider = gr.Slider(0, 40, value=0, label="Morphological Closing (Glue / Default=0)")
                opening_slider = gr.Slider(0, 10, value=0, step=1, label="Morphological Opening (Denoise / Default=0)")

                with gr.Accordion("🎨 HSV Detailed Parameters", open=False):
                    with gr.Row():
                        h_min = gr.Slider(0, 180, value=10, label="H Min")
                        h_max = gr.Slider(0, 180, value=90, label="H Max")
                    with gr.Row():
                        s_min = gr.Slider(0, 255, value=30, label="S Min")
                        s_max = gr.Slider(0, 255, value=255, label="S Max")
                    with gr.Row():
                        v_min = gr.Slider(0, 255, value=40, label="V Min")
                        v_max = gr.Slider(0, 255, value=255, label="V Max")

            folder_input = gr.Textbox(value="output", label="Output Directory")
            save_btn = gr.Button("💾 Save Cropped Object", variant="primary")
            status_text = gr.Textbox(label="System Status", value="Ready")

        with gr.Column(scale=4):
            mask_output = gr.Image(label="Mask Preview", type="numpy", image_mode="L")
            vis_output = gr.Image(label="Interactive Overlay", type="numpy")
            final_output = gr.Image(label="Saved Result Preview", type="numpy")

    filter_inputs = [use_hsv_box, h_min, h_max, s_min, s_max, v_min, v_max, glue_slider, opening_slider]

    def on_preset_change(preset):
        return get_preset_hsv(preset)

    color_preset.change(
        on_preset_change,
        inputs=[color_preset],
        outputs=[use_hsv_box, h_min, h_max, s_min, s_max, v_min, v_max]
    )

    input_img.upload(process_image_upload, inputs=[input_img], outputs=[input_img, status_text])

    input_img.select(handle_click, inputs=[mode_radio] + filter_inputs, outputs=[mask_output, vis_output])

    for inp in filter_inputs:
        inp.change(handle_slider_change, inputs=filter_inputs, outputs=[mask_output, vis_output])

    undo_btn.click(undo_last_point, inputs=filter_inputs, outputs=[mask_output, vis_output])
    save_btn.click(save_file, inputs=[folder_input], outputs=[final_output, status_text])

if __name__ == "__main__":
    app.launch()