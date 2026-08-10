#!/bin/bash
#SBATCH --job-name=rna_htseq
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err

set -Eeuo pipefail

# =====================================================
# 1. User Configurable Variables (Passed via --export)
# =====================================================
BAM_DIR="${BAM_DIR:-}"
BAM_LIST="${BAM_LIST:-}"
GTF="${GTF:-}"
OUT_DIR="${OUT_DIR:-}"
STRAND="${STRAND:-no}"
FEATURE_TYPE="${FEATURE_TYPE:-CDS}"
ID_ATTR="${ID_ATTR:-gene_id}"

if [[ -z "$BAM_DIR" || -z "$BAM_LIST" || -z "$GTF" || -z "$OUT_DIR" ]]; then
    echo "[ERR] Missing required arguments. Please provide via --export:" >&2
    echo "      BAM_DIR, BAM_LIST, GTF, OUT_DIR" >&2
    exit 1
fi

mkdir -p "$OUT_DIR" "logs"

# =====================================================
# 2. Environment & Tool Check
# =====================================================
module purge || true
module load HTSeq/2.0.2 2>/dev/null || true

command -v htseq-count >/dev/null 2>&1 || { 
    echo "[ERR] htseq-count command not found. Check module load." >&2
    exit 3 
}

# =====================================================
# 3. Array Task Assignment & Execution
# =====================================================
: "${SLURM_ARRAY_TASK_ID:?Submit with --array=0-(N-1)}"

BAM_NAME=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$BAM_LIST")
[[ -n "${BAM_NAME:-}" ]] || { 
    echo "[ERR] No BAM entry found for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID" >&2
    exit 2 
}

BAM_PATH="$BAM_DIR/$BAM_NAME"
SAMPLE="${BAM_NAME%.bam}"
# 如果你的檔名帶有其他結尾（例如 .name_sorted.bam），可以根據實際情況調整：
# SAMPLE="${SAMPLE%.name_sorted}"

OUT_COUNTS="$OUT_DIR/${SAMPLE}.counts.txt"

# 防呆：檢查檔案是否存在
[[ -f "$BAM_PATH" ]] || { echo "[ERR] BAM not found: $BAM_PATH" >&2; exit 4; }
[[ -f "$GTF" ]]      || { echo "[ERR] GTF not found: $GTF" >&2; exit 5; }

# 防重複執行 (Resume 機制)
if [[ -s "$OUT_COUNTS" ]]; then
    echo "[INFO] Count file already exists and is non-empty. Skipping: $SAMPLE"
    exit 0
fi

echo "========================================="
echo "Task ID     : $SLURM_ARRAY_TASK_ID"
echo "Sample      : $SAMPLE"
echo "BAM Path    : $BAM_PATH"
echo "GTF Path    : $GTF"
echo "Stranded    : $STRAND"
echo "========================================="

htseq-count \
  -f bam \
  -r pos \
  -s "$STRAND" \
  -a 10 \
  -t "$FEATURE_TYPE" \
  -i "$ID_ATTR" \
  "$BAM_PATH" \
  "$GTF" \
  > "$OUT_COUNTS"

echo "[SUCCESS] Finished htseq-count for $SAMPLE -> $OUT_COUNTS"