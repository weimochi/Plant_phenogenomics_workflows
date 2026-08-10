#!/usr/bin/env python3
import pandas as pd
import glob, os, sys
from datetime import datetime

# ============================================================
# Final Metadata & Statistics Merger
# ============================================================
# Paths can be overridden via environment variables or arguments
BASE_DIR = os.environ.get("BASE_DIR", ".")
SAMPLE_LIST = os.environ.get("SAMPLE_LIST", "")
MAP_DIR = os.environ.get("MAP_DIR", ".")
META_OUT = os.environ.get("META_OUT", os.environ.get("OUTDIR", "./merged_tables"))

os.makedirs(META_OUT, exist_ok=True)

if not SAMPLE_LIST or not os.path.exists(SAMPLE_LIST):
    sys.exit("[ERR] Please provide a valid SAMPLE_LIST path via environment variable.")

# ---------- Find latest mapping/NM comparison TSV ----------
mapping_tsv = sys.argv[1] if len(sys.argv) > 1 else None
if mapping_tsv is None:
    files = sorted(glob.glob(os.path.join(MAP_DIR, "*_compare_*.tsv")) + glob.glob(os.path.join(MAP_DIR, "mapping_rates_*.tsv")))
    if not files:
        sys.exit(f"[ERR] No summary/comparison TSV found in {MAP_DIR}. Run 02a first.")
    mapping_tsv = files[-1]

print(f"[INFO] Sample list : {SAMPLE_LIST}")
print(f"[INFO] Stats report: {mapping_tsv}")
print(f"[INFO] Output dir : {META_OUT}")

# ---------- Read Sample List ----------
sl = pd.read_csv(SAMPLE_LIST, sep="\t", dtype=str).fillna("")
# Standardize sample column name if needed (e.g., RNAseq -> sample_id)
if "RNAseq" in sl.columns and "sample_id" not in sl.columns:
    sl = sl.rename(columns={"RNAseq": "sample_id"})

# ---------- Read Statistics Report ----------
mapdf = pd.read_csv(mapping_tsv, sep="\t", dtype=str).fillna("")

# ---------- Merge ----------
md = sl.merge(mapdf, how="left", on="sample_id")

# ---------- Reorder Columns (Put identifiers upfront) ----------
front_cols = [c for c in ["sample_id", "OurID", "Company", "CommonName", "ProductName"] if c in md.columns]
other_cols = [c for c in md.columns if c not in front_cols]
md = md[front_cols + other_cols]

# ---------- Export ----------
ts = datetime.now().strftime("%Y%m%d_%H%M")
out_tsv = os.path.join(META_OUT, f"metadata_merged_{ts}.tsv")
md.to_csv(out_tsv, sep="\t", index=False)

print(f"\n[OK] Metadata and statistics merged successfully.")
print(f" Output TSV: {out_tsv}")