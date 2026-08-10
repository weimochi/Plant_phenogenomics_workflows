#!/usr/bin/env bash
#SBATCH --job-name=filter_vcf
#SBATCH --partition=CPU
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -Eeuo pipefail

# =====================================================
# 1. User Configurable Variables (Passed via --export)
# =====================================================
IN_VCF="${IN_VCF:-}"
OUT_VCF="${OUT_VCF:-}"
MIN_DP="${MIN_DP:-10}"             # 預設最低深度 10
MIN_GQ="${MIN_GQ:-20}"             # 預設最低基因型品質 20
MIN_MAF="${MIN_MAF:-0.05}"         # 預設最低次要等位基因頻率 5%
MAX_MISSING="${MAX_MISSING:-0.2}"  # 預設最高容許缺失率 20%
THREADS="${SLURM_CPUS_PER_TASK:-8}"
FAST_STORAGE="${FAST_STORAGE:-/tmp}"

if [[ -z "$IN_VCF" || -z "$OUT_VCF" ]]; then
  echo "[ERR] Missing required vars: IN_VCF, OUT_VCF" >&2
  exit 1
fi

OUTDIR="$(dirname "$OUT_VCF")"
mkdir -p "$OUTDIR" "logs"

# =====================================================
# 2. Modules & Environment
# =====================================================
module purge || true
module load bcftools 2>/dev/null || true
module load htslib 2>/dev/null || true

command -v bcftools >/dev/null 2>&1 || { echo "[ERR] bcftools not found"; exit 3; }

# =====================================================
# 3. SSD Temp Directories Setup
# =====================================================
OUTBASE="$(basename "$OUT_VCF")"
TMPDIR="${FAST_STORAGE}/${USER}_bcftools_tmp_${OUTBASE}_${SLURM_JOB_ID:-$$}"
mkdir -p "$TMPDIR"
echo "[INFO] Using SSD TMPDIR: $TMPDIR"

trap 'echo "[INFO] Cleaning up TMPDIR..."; rm -rf "$TMPDIR" 2>/dev/null || true' EXIT

# =====================================================
# 4. VCF Filtering Execution (High-speed Pipeline)
# =====================================================
echo "[RUN] Filtering VCF: $IN_VCF -> $OUT_VCF"
echo "[CONF] DP < $MIN_DP | GQ < $MIN_GQ | MAF >= $MIN_MAF | F_MISSING <= $MAX_MISSING"

bcftools view -v snps -m2 -M2 -O u "$IN_VCF" \
  | bcftools +setGT -O u -- -t q -n . -i "FMT/DP<${MIN_DP} || FMT/GQ<${MIN_GQ}" \
  | bcftools +fill-tags -O u -- -t AC,AN,AF,MAF,F_MISSING \
  | bcftools view -i "AC>0 && MAF>=${MIN_MAF} && F_MISSING<=${MAX_MISSING}" -O u \
  | bcftools sort -T "$TMPDIR" -O z -o "$OUT_VCF"

# =====================================================
# 5. Indexing & Validation
# =====================================================
tabix -f -p vcf "$OUT_VCF"

echo "=== Filtering Summary ==="
echo "Number of taxa:  $(bcftools query -l "$OUT_VCF" | wc -l)"
echo "Number of sites: $(bcftools view -H "$OUT_VCF" | wc -l)"
echo "Filtering Finished Successfully!"