#!/usr/bin/env bash
#SBATCH --job-name=RNA-fastqc
#SBATCH --partition=CPU
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2; }
error_exit() { log "ERROR: $*"; exit 1; }

# ========= Configurable Parameters =========
INPUT_DIR="${INPUT_DIR:-}"
OUTPUT_DIR="${OUTPUT_DIR:-}"
PATTERN="${PATTERN:-*.[fF]*[qQ]*.gz}" # Matches .fq.gz, .fastq.gz
FILELIST="${FILELIST:-}"
FASTQC_BIN="${FASTQC_BIN:-fastqc}"

log "=========================================="
log "FastQC Array Job Started"
log "=========================================="

if [[ -z "${OUTPUT_DIR}" || -z "${INPUT_DIR}" ]]; then
    error_exit "INPUT_DIR and OUTPUT_DIR must be exported."
fi

mkdir -p "logs" "${OUTPUT_DIR}/tmp" "${OUTPUT_DIR}/logs" || error_exit "Failed to create directories"

# Generate or use existing file list
LIST="${FILELIST:-${OUTPUT_DIR}/filelist.txt}"
if [[ ! -s "$LIST" ]]; then
    log "Generating file list from: ${INPUT_DIR}"
    find -L "${INPUT_DIR}" -maxdepth 1 -type f -name "${PATTERN}" | sort > "${LIST}"
fi

: "${SLURM_ARRAY_TASK_ID:?Submit with --array=0-(N-1)}"

F="$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "${LIST}" 2>/dev/null || true)"
[[ -n "${F}" && -r "${F}" ]] || error_exit "File not found or not readable at index ${SLURM_ARRAY_TASK_ID}"

BASENAME=$(basename "${F}")
START_TIME=$(date +%s)

log "Processing: ${BASENAME}"

# Execute FastQC
if ${FASTQC_BIN} \
    --extract \
    -t "${SLURM_CPUS_PER_TASK}" \
    -o "${OUTPUT_DIR}" \
    --dir "${OUTPUT_DIR}/tmp" \
    "${F}" 2>&1 | tee "${OUTPUT_DIR}/logs/${BASENAME}.fastqc.log"; then

    # Smart extension stripping for FastQC verification
    PREFIX="${BASENAME%.gz}"
    PREFIX="${PREFIX%.fastq}"
    PREFIX="${PREFIX%.fq}"
    
    EXPECTED_ZIP="${OUTPUT_DIR}/${PREFIX}_fastqc.zip"
    EXPECTED_DIR="${OUTPUT_DIR}/${PREFIX}_fastqc"

    if [[ ! -f "${EXPECTED_ZIP}" ]] && [[ ! -d "${EXPECTED_DIR}" ]]; then
        log "WARNING: Expected output files not found (${EXPECTED_ZIP})"
        exit 1
    else
        log "FastQC completed successfully in $(( $(date +%s) - START_TIME ))s"
    fi
else
    error_exit "FastQC execution failed for ${BASENAME}"
fi