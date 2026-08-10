#!/usr/bin/env bash
#SBATCH --job-name=glnexus
#SBATCH --partition=CPU
#SBATCH --cpus-per-task=40
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -Eeuo pipefail

# =====================================================
# 1. User Configurable Variables (Passed via --export)
# =====================================================
IN_DIR="${IN_DIR:-}"
OUT_VCF="${OUT_VCF:-}"
GLNEXUS_CONFIG="${GLNEXUS_CONFIG:-DeepVariantWES}"
EXCLUDE_SAMPLES="${EXCLUDE_SAMPLES:-}"
EXCLUDE_FILE="${EXCLUDE_FILE:-}"
THREADS="${SLURM_CPUS_PER_TASK:-40}"
FAST_STORAGE="${FAST_STORAGE:-/tmp}"         # 高速暫存路徑，預設使用系統 /tmp

if [[ -z "$IN_DIR" || -z "$OUT_VCF" ]]; then
  echo "[ERR] Missing required vars: IN_DIR, OUT_VCF" >&2
  exit 1
fi

OUTDIR="$(dirname "$OUT_VCF")"
mkdir -p "$OUTDIR" "logs"

# =====================================================
# 2. Modules & Environment
# =====================================================
module purge || true
module load GLnexus/1.4.1 2>/dev/null || true
module load bcftools 2>/dev/null || true
module load htslib 2>/dev/null || true

command -v glnexus_cli >/dev/null 2>&1 || { echo "[ERR] glnexus_cli not found"; exit 3; }

# =====================================================
# 3. SSD Temp Directories Setup & Safe Cleanup
# =====================================================
WORKROOT="${FAST_STORAGE}/${USER}_glnexus_work"
WORKDIR="${WORKROOT}/$(basename "$OUT_VCF" .vcf.gz)_${SLURM_JOB_ID:-$$}"
mkdir -p "$WORKDIR"
echo "[INFO] Using SSD WORKDIR: $WORKDIR"

# 【修正 2】宣告 Trap：無論腳本如何結束，必定連根拔起清除整個 WORKDIR
trap 'echo "[INFO] Cleaning up SSD Temp..."; rm -rf "$WORKDIR" 2>/dev/null || true' EXIT

# =====================================================
# 4. Exclude List Processing
# =====================================================
EX_RE=""
EX_LIST_FILE="${WORKDIR}/exclude.samples.txt"
: > "$EX_LIST_FILE"

if [[ -n "$EXCLUDE_SAMPLES" ]]; then
  echo "$EXCLUDE_SAMPLES" | tr ',' '\n' | awk 'NF{print $1}' >> "$EX_LIST_FILE"
fi

if [[ -n "$EXCLUDE_FILE" && -s "$EXCLUDE_FILE" ]]; then
  awk 'NF{print $1}' "$EXCLUDE_FILE" >> "$EX_LIST_FILE"
fi

sort -u "$EX_LIST_FILE" -o "$EX_LIST_FILE"

if [[ -s "$EX_LIST_FILE" ]]; then
  EX_RE="$(paste -sd'|' "$EX_LIST_FILE")"
  echo "[CONF] Excluding $(wc -l < "$EX_LIST_FILE") samples"
fi

# =====================================================
# 5. File Discovery & Index Check
# =====================================================
LIST_FILE="${WORKDIR}/gvcf.list"
: > "$LIST_FILE"

if [[ -n "$EX_RE" ]]; then
  # 【修正 1】精準正則比對：只比對路徑最尾端的 "斜線 + 排除樣本名稱 + 任意字元 + g.vcf.gz"
  # 這樣絕對不會誤殺名稱相似的父資料夾！
  find -L "$IN_DIR" -type f -name "*.g.vcf.gz" -size +0c | grep -Ev "/(${EX_RE})[^/]*g\.vcf\.gz$" | sort > "$LIST_FILE"
else
  find -L "$IN_DIR" -type f -name "*.g.vcf.gz" -size +0c | sort > "$LIST_FILE"
fi

N_GVCF=$(wc -l < "$LIST_FILE" 2>/dev/null || echo 0)
[[ "$N_GVCF" -gt 0 ]] || { echo "[ERR] No gVCF files found in $IN_DIR"; exit 4; }
echo "[INFO] Total gVCFs to joint call: $N_GVCF"

echo "[INFO] Checking .tbi indices..."
while read -r V; do
  [[ -f "${V}.tbi" ]] || tabix -p vcf "$V"
done < "$LIST_FILE"

# =====================================================
# 6. GLnexus Execution
# =====================================================
TMPDIR="$(mktemp -d "${WORKDIR}/.tmp_joint.XXXXXX")"

echo "[RUN] GLnexus (Config: $GLNEXUS_CONFIG)"
(
  cd "$WORKDIR"
  rm -rf GLnexus.DB 2>/dev/null || true

  glnexus_cli --config "$GLNEXUS_CONFIG" \
              --threads "$THREADS" \
              --mem-gbytes 80 \
              --list "$LIST_FILE" \
    | bcftools sort -m 2G -T "$TMPDIR" \
    | bcftools view -Oz -o "$OUT_VCF"
)

tabix -p vcf "$OUT_VCF"
echo "Joint Calling Finished -> $OUT_VCF"