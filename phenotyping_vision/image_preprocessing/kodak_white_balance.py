import tkinter as tk
from tkinter import filedialog, messagebox
import cv2
import numpy as np
from PIL import Image, ImageTk
from pathlib import Path
import csv

class WhiteBalanceApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Reference-Patch Relative White Balance Tool")
        self.root.geometry("1000x800")

        # --- Variables Initialization ---
        self.original_image = None   # Original OpenCV image (BGR)
        self.processed_image = None  # Processed OpenCV image (BGR)
        self.display_image = None    # Image displayed on Canvas (PIL Object)
        self.scale_factor = 1.0      # Scale factor for coordinate mapping
        self._resize_timer = None    # Timer for resize debounce
        
        # File and QC tracking variables
        self.current_file_path = None
        self.last_calibration_meta = None
        self.ROI_RADIUS = 10         # 21x21 sampling window (radius=10)

        # --- UI Layout ---
        control_frame = tk.Frame(root, pady=10)
        control_frame.pack(side=tk.TOP, fill=tk.X)

        btn_load = tk.Button(control_frame, text="1. Load Image", command=self.load_image, font=("Arial", 12), bg="#e1f5fe")
        btn_load.pack(side=tk.LEFT, padx=10)

        self.lbl_instruction = tk.Label(control_frame, text="Please load an image, then click the reference white patch (Clicking again resets to original).", font=("Arial", 11), fg="gray")
        self.lbl_instruction.pack(side=tk.LEFT, padx=10)

        btn_save = tk.Button(control_frame, text="3. Save Result & Log", command=self.save_image, font=("Arial", 12), bg="#e8f5e9")
        btn_save.pack(side=tk.RIGHT, padx=10)

        btn_reset = tk.Button(control_frame, text="Reset", command=self.reset_image, font=("Arial", 12))
        btn_reset.pack(side=tk.RIGHT, padx=10)

        # Image Display Canvas
        self.canvas_frame = tk.Frame(root)
        self.canvas_frame.pack(fill=tk.BOTH, expand=True)

        self.canvas = tk.Canvas(self.canvas_frame, bg="#333333", cursor="cross")
        self.canvas.pack(fill=tk.BOTH, expand=True)

        self.canvas.bind("<Button-1>", self.on_image_click)

    def load_image(self):
        file_path = filedialog.askopenfilename(filetypes=[("Images", "*.jpg;*.jpeg;*.png;*.webp;*.bmp;*.tif")])
        if not file_path:
            return

        img = cv2.imread(file_path)
        if img is None:
            messagebox.showerror("Error", "Failed to load image. Please check if the path or filename contains special characters.")
            return

        self.current_file_path = file_path
        self.original_image = img
        self.processed_image = img.copy()
        self.last_calibration_meta = None
        self.show_image(self.processed_image)
        self.lbl_instruction.config(text=f"Loaded: {Path(file_path).name}. Click the reference white patch.", fg="blue")

    def show_image(self, cv_img):
        if cv_img is None:
            return

        img_rgb = cv2.cvtColor(cv_img, cv2.COLOR_BGR2RGB)
        h, w = img_rgb.shape[:2]

        canvas_w = self.canvas.winfo_width()
        canvas_h = self.canvas.winfo_height()

        if canvas_w < 100 or canvas_h < 100:
            canvas_w, canvas_h = 980, 700

        scale_w = canvas_w / w
        scale_h = canvas_h / h
        self.scale_factor = min(scale_w, scale_h, 1.0)

        new_w = max(1, int(w * self.scale_factor))
        new_h = max(1, int(h * self.scale_factor))

        img_resized = cv2.resize(img_rgb, (new_w, new_h), interpolation=cv2.INTER_AREA)

        pil_img = Image.fromarray(img_resized)
        self.display_image = ImageTk.PhotoImage(pil_img)

        self.canvas.delete("all")
        x_center = canvas_w // 2
        y_center = canvas_h // 2
        self.canvas.create_image(x_center, y_center, anchor=tk.CENTER, image=self.display_image)

    def on_image_click(self, event):
        if self.original_image is None:
            return

        canvas_w = self.canvas.winfo_width()
        canvas_h = self.canvas.winfo_height()

        img_w = int(self.original_image.shape[1] * self.scale_factor)
        img_h = int(self.original_image.shape[0] * self.scale_factor)
        offset_x = (canvas_w - img_w) // 2
        offset_y = (canvas_h - img_h) // 2

        click_x = event.x - offset_x
        click_y = event.y - offset_y

        if click_x < 0 or click_x >= img_w or click_y < 0 or click_y >= img_h:
            return

        real_x = int(click_x / self.scale_factor)
        real_y = int(click_y / self.scale_factor)

        self.perform_relative_white_balance(real_x, real_y)

    def perform_relative_white_balance(self, x, y):
        h_img, w_img = self.original_image.shape[:2]

        # 🛡️ Boundary protection using adjustable ROI_RADIUS
        x_min = max(0, x - self.ROI_RADIUS)
        x_max = min(w_img, x + self.ROI_RADIUS + 1)
        y_min = max(0, y - self.ROI_RADIUS)
        y_max = min(h_img, y + self.ROI_RADIUS + 1)

        roi = self.original_image[y_min:y_max, x_min:x_max]
        if roi.size == 0:
            return

        # 🔍 Pixel-level Clipping QC Check (Reject if > 5% pixels are saturated >= 250)
        clipped_fraction = np.mean(roi >= 250, axis=(0, 1))
        b_clip, g_clip, r_clip = clipped_fraction

        if max(r_clip, g_clip, b_clip) > 0.05:
            messagebox.showerror(
                "Calibration Unreliable (Clipping)",
                f"The selected reference patch contains clipped/saturated pixels.\n\n"
                f"Clipped fraction:\n"
                f"R: {r_clip*100:.1f}%\n"
                f"G: {g_clip*100:.1f}%\n"
                f"B: {b_clip*100:.1f}%\n\n"
                "Please mark this image as unreliable or use another image captured without saturation."
            )
            return

        # Use median to protect against dust, glare, or local blemishes
        raw_b, raw_g, raw_r = np.median(roi, axis=(0, 1))

        print(f"Sampled Patch Median (BGR) -> B: {raw_b:.1f}, G: {raw_g:.1f}, R: {raw_r:.1f}")

        # 🎯 Relative White Balance Logic (Exposure Preserved via Neutral Target Mean)
        target = (raw_r + raw_g + raw_b) / 3.0

        safe_r = max(raw_r, 1e-5)
        safe_g = max(raw_g, 1e-5)
        safe_b = max(raw_b, 1e-5)

        gain_r = target / safe_r
        gain_g = target / safe_g
        gain_b = target / safe_b

        # ⚠️ Gain QC Warning
        gains = np.array([gain_r, gain_g, gain_b])
        large_gain_warning = bool(np.any(gains < 0.7) or np.any(gains > 1.4))
        
        if large_gain_warning:
            resp = messagebox.askyesno(
                "Large White-Balance Correction Warning",
                f"Unusually large correction detected:\n\n"
                f"R gain = {gain_r:.3f}\n"
                f"G gain = {gain_g:.3f}\n"
                f"B gain = {gain_b:.3f}\n\n"
                "Please verify that the correct reference white patch was selected.\n"
                "Continue anyway?"
            )
            if not resp:
                return

        # Always calculate from the original image (prevents cumulative errors)
        img_float = self.original_image.astype(np.float32)
        img_b, img_g, img_r = cv2.split(img_float)

        img_r = cv2.multiply(img_r, gain_r)
        img_g = cv2.multiply(img_g, gain_g)
        img_b = cv2.multiply(img_b, gain_b)

        img_merged = cv2.merge([img_b, img_g, img_r])
        self.processed_image = np.clip(img_merged, 0, 255).astype(np.uint8)

        # Store comprehensive audit metadata
        self.last_calibration_meta = {
            "source_file": self.current_file_path,
            "x": x, "y": y,
            "x_min": x_min, "x_max": x_max, "y_min": y_min, "y_max": y_max,
            "roi_n_pixels": roi.shape[0] * roi.shape[1],
            "ref_R": raw_r, "ref_G": raw_g, "ref_B": raw_b,
            "target": target,
            "gain_R": gain_r, "gain_G": gain_g, "gain_B": gain_b,
            "r_clip": r_clip, "g_clip": g_clip, "b_clip": b_clip,
            "large_gain_warning": large_gain_warning
        }

        self.show_image(self.processed_image)
        self.lbl_instruction.config(
            text=f"Calibrated! Gains: R={gain_r:.3f}, G={gain_g:.3f}, B={gain_b:.3f} (Clicking again recalculates from original)", fg="green"
        )

    def save_image(self):
        # 🛡️ Protection: Block saving if no calibration has been applied yet
        if self.last_calibration_meta is None or self.processed_image is None:
            messagebox.showwarning(
                "Warning",
                "No white-balance calibration has been applied yet.\n"
                "Please select the reference white patch before saving."
            )
            return

        file_path = filedialog.asksaveasfilename(
            defaultextension=".jpg",
            filetypes=[("JPEG", "*.jpg"), ("PNG", "*.png"), ("TIFF", "*.tif")]
        )
        if file_path:
            ext = file_path.split('.')[-1]
            success, encoded_img = cv2.imencode(f'.{ext}', self.processed_image)
            if success:
                with open(file_path, mode='wb') as f:
                    f.write(encoded_img)
                
                # Automatically log calibration metadata to a CSV file for batch effect auditing
                log_csv = Path(file_path).parent / "calibration_log.csv"
                file_exists = log_csv.exists()
                
                # Use a copy to prevent mutating the original metadata dictionary
                log_row = self.last_calibration_meta.copy()
                log_row["output_file"] = file_path
                
                fieldnames = [
                    "source_file", "output_file", "x", "y", 
                    "x_min", "x_max", "y_min", "y_max", "roi_n_pixels",
                    "ref_R", "ref_G", "ref_B", "target", 
                    "gain_R", "gain_G", "gain_B", 
                    "r_clip", "g_clip", "b_clip", "large_gain_warning"
                ]
                
                with open(log_csv, mode="a", newline="", encoding="utf-8") as csv_file:
                    writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
                    if not file_exists:
                        writer.writeheader()
                    writer.writerow(log_row)
                    
                print(f"[LOG] Metadata successfully appended to {log_csv}")

                messagebox.showinfo("Success", f"Image successfully saved and logged to:\n{file_path}")
            else:
                messagebox.showerror("Error", "Failed to save image. Please check the file path.")

    def reset_image(self):
        if self.original_image is not None:
            self.processed_image = self.original_image.copy()
            self.last_calibration_meta = None
            self.show_image(self.processed_image)
            self.lbl_instruction.config(text="Image reset. Please select the reference white patch again.", fg="gray")

if __name__ == "__main__":
    root = tk.Tk()
    app = WhiteBalanceApp(root)

    def on_resize(event):
        if event.widget == root and app.processed_image is not None:
            if app._resize_timer:
                root.after_cancel(app._resize_timer)
            app._resize_timer = root.after(150, lambda: app.show_image(app.processed_image))

    root.bind("<Configure>", on_resize)
    root.mainloop()