# Image Preprocessing Pipeline for Plant Phenomics

This module contains production-ready preprocessing tools designed for high-throughput plant phenotyping images, optimized for Brassica species. The workflow ensures consistent color baselines and accurate object isolation prior to downstream feature extraction and machine learning.

---

## Step-by-Step Preprocessing Workflow

### Step 1: Color Calibration (`kodak_white_balance.py`)

Raw images captured across different batches or lighting conditions often suffer from chromatic cast. This tool corrects color without altering the overall exposure.

* **Core Methodology (Reference-Patch Relative White Balance):**
* **Fixed Reference:** Users click on the neutral white patch of the Kodak Q-13 color card.
* **Median Sampling:** Extracts median BGR values using a robust $21 \times 21$ ROI to avoid dust, glare, or blemishes.
* **Clipping Quality Control:** Rejects samples if over $5\%$ of pixels in the patch are saturated ($\ge 250$).
* **Gain Calculation:** Computes channel-wise multiplicative gains using the mean intensity of the patch as a neutral target (`target = (R + G + B) / 3`). This avoids explicit exposure normalization, preserving natural brightness variations.
* **Audit Logging:** Automatically appends calibration metadata to `calibration_log.csv` for batch effect auditing.


* **Usage:**
```bash
python kodak_white_balance.py

```



---

### Step 2: Object Isolation & Segmentation (`interactive_segmentation.py`)

After color calibration, images are processed using a hybrid segmentation pipeline combining deep learning object priors with color constraints.

* **Core Methodology:**
* **SAM (Segment Anything Model):** Leverages deep learning object prompts to locate plant canopy boundaries.
* **HSV Constraints:** Restricts the SAM mask with preset color ranges to filter out background soil, pots, and color charts.
* **Multi-Component Tracking:** Retains all connected components containing positive user clicks, handling segmented or multi-part plant canopies.
* **Conservative Morphology:** Defaults morphological closing and opening to zero to prevent unwanted destruction of delicate leaf tips and margins.
* **Standardized Export:** Crops the isolated plant object and pads it into a centered square RGB image (`*_segmented.png`).


* **Usage:**
```bash
python interactive_segmentation.py

```

---

## Best Practices for Batch Processing

1. Always click the **same neutral white patch** on the Kodak Q-13 card across all images to maintain a consistent baseline.
2. Review the `calibration_log.csv` generated during white balance to screen for outliers or large correction warnings (`large_gain_warning == True`).