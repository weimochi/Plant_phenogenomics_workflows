#!/usr/bin/env bash
#SBATCH --job-name=markdup
#SBATCH --partition=CPU
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail
: "${SLURM_ARRAY_TASK_ID:?Submit with --array=0-(N-1)}"

# =====================================================
# 1. User Configurable Variables (Passed via --export)
# =====================================================
BAMDIR="${BAMDIR:-}"
OUTDIR="${OUTDIR:-}"
TAG="${TAG:-md}"          # 預設加上 md 標籤，用來跟原始 bam 區分
THREADS="${SLURM_CPUS_PER_TASK:-8}"

if [[ -z "$BAMDIR" || -z "$OUTDIR" ]]; then
  echo "[ERR] Required env vars missing: BAMDIR, OUTDIR. Set via --export" >&2
  exit 1
fi

mkdir -p "$OUTDIR" "logs"

# =====================================================
# 2. Modules & Tools
# =====================================================
module purge || true
module add gatk/4.3.0.0 2>/dev/null || true
module add samtools/1.16.1 2>/dev/null || true

# =====================================================
# 3. Dynamic File List & Array Assignment
# =====================================================
# TAG 動態處理：有值則為 ".TAG"，無值則不加點 (防止出現 sample..bam)
TAG_SUFFIX="${TAG:+.${TAG}}"

LIST="${OUTDIR}/bam_list.txt"
if [[ ! -s "$LIST" ]]; then
  # 防競爭危害 (Race Condition) 的安全寫入法
  TMP_LIST="${LIST}.${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.tmp"
  find -L "$BAMDIR" -maxdepth 1 -type f -name "*.bam" ! -name "*.md.bam" | sort > "$TMP_LIST"
  # mv -n 代表「不覆蓋」。只有第一個衝到的任務能建立成功，其他任務會優雅地失敗並刪除自己的暫存檔
  mv -n "$TMP_LIST" "$LIST" 2>/dev/null || rm -f "$TMP_LIST"
  sleep 2 # 給網路硬碟 (NFS) 一點時間同步，確保大家都能讀到剛建好的檔案
fi

# 依據號碼牌 (TASK_ID) 去菜單 (LIST) 裡面點菜
BAM="$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$LIST")"
[[ -n "$BAM" ]] || { echo "[ERR] Empty BAM at index $SLURM_ARRAY_TASK_ID" >&2; exit 2; }

SAMPLE="$(basename "$BAM" .bam)"
MD_BAM="${OUTDIR}/${SAMPLE}${TAG_SUFFIX}.bam"
MD_MET="${OUTDIR}/${SAMPLE}${TAG_SUFFIX}.metrics.txt"

if [[ -s "$MD_BAM" && -s "${MD_BAM}.bai" ]]; then
  echo "Already done: $SAMPLE"
  exit 0
fi

# =====================================================
# 4. Coordinate Sort Check & Execution
# =====================================================
USE_BAM="$BAM"
SORTED_BAM="${OUTDIR}/${SAMPLE}.sorted.bam"

# Check if BAM is coordinate-sorted
if ! samtools view -H "$BAM" | grep -q '^@HD.*SO:coordinate'; then
  echo "[SORT] $SAMPLE needs coordinate sorting -> $SORTED_BAM"
  samtools sort -@ "$THREADS" -o "$SORTED_BAM" "$BAM"
  USE_BAM="$SORTED_BAM"
fi

echo "[RUN] MarkDuplicates: $SAMPLE -> $MD_BAM"
gatk --java-options "-Xmx16g -Xms8g" MarkDuplicates \
  -I "$USE_BAM" \
  -O "$MD_BAM" \
  -M "$MD_MET" \
  --REMOVE_DUPLICATES false \
  --VALIDATION_STRINGENCY LENIENT

samtools index -@ "$THREADS" "$MD_BAM"

# Clean up intermediate sorted bam to save space
[[ "$USE_BAM" == "$SORTED_BAM" ]] && rm -f "$SORTED_BAM"

echo "Finished: $SAMPLE"