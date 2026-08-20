#!/usr/bin/env bash
#SBATCH --job-name=fastp-trim
#SBATCH --partition=CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -Eeuo pipefail

# ==============================================================================
# 00b_fastp_trim.slurm
# Run fastp on one paired-end sample per Slurm array task.
#
# MANIFEST format, without header:
# sample<TAB>R1<TAB>R2
# ==============================================================================

# ---------- Required inputs ----------

MANIFEST="${MANIFEST:-}"
OUTPUT_DIR="${OUTPUT_DIR:-}"

# ---------- Executable and resources ----------

FASTP_BIN="${FASTP_BIN:-fastp}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# ---------- Trimming parameters ----------

# Fixed number of bases removed from the 5' end.
# Keep at 0 unless FastQC or library design supports fixed trimming.
TRIM_FRONT="${TRIM_FRONT:-0}"

QUALIFIED_QUALITY="${QUALIFIED_QUALITY:-20}"
UNQUALIFIED_PERCENT="${UNQUALIFIED_PERCENT:-40}"
MIN_LENGTH="${MIN_LENGTH:-80}"

CUT_WINDOW_SIZE="${CUT_WINDOW_SIZE:-16}"
CUT_MEAN_QUALITY="${CUT_MEAN_QUALITY:-20}"

# 1: enable; 0: disable
CUT_FRONT="${CUT_FRONT:-1}"
CUT_TAIL="${CUT_TAIL:-1}"
TRIM_POLY_G="${TRIM_POLY_G:-1}"
TRIM_POLY_X="${TRIM_POLY_X:-0}"

# Leave both empty to use fastp paired-end adapter detection.
# If adapters are supplied, both R1 and R2 sequences are required.
ADAPTER_R1="${ADAPTER_R1:-}"
ADAPTER_R2="${ADAPTER_R2:-}"

OUT_SUFFIX="${OUT_SUFFIX:-trimmed}"

# ---------- Logging ----------

log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2
}

die() {
    log "ERROR: $*"
    exit 1
}

# ---------- Helper functions ----------

require_nonnegative_integer() {
    local name="$1"
    local value="$2"

    [[ "${value}" =~ ^[0-9]+$ ]] ||
        die "${name} must be a non-negative integer: ${value}"
}

require_boolean() {
    local name="$1"
    local value="$2"

    [[ "${value}" == "0" || "${value}" == "1" ]] ||
        die "${name} must be 0 or 1: ${value}"
}

# ---------- Validation ----------

[[ -n "${MANIFEST}" ]] ||
    die "MANIFEST is required."

[[ -s "${MANIFEST}" ]] ||
    die "MANIFEST does not exist or is empty: ${MANIFEST}"

[[ -n "${OUTPUT_DIR}" ]] ||
    die "OUTPUT_DIR is required."

[[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] ||
    die "This script must be submitted as a Slurm array job."

[[ "${SLURM_ARRAY_TASK_ID}" =~ ^[0-9]+$ ]] ||
    die "Invalid SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"

command -v "${FASTP_BIN}" >/dev/null 2>&1 ||
    die "fastp executable not found: ${FASTP_BIN}"

require_nonnegative_integer "TRIM_FRONT" "${TRIM_FRONT}"
require_nonnegative_integer "QUALIFIED_QUALITY" "${QUALIFIED_QUALITY}"
require_nonnegative_integer "UNQUALIFIED_PERCENT" "${UNQUALIFIED_PERCENT}"
require_nonnegative_integer "MIN_LENGTH" "${MIN_LENGTH}"
require_nonnegative_integer "CUT_WINDOW_SIZE" "${CUT_WINDOW_SIZE}"
require_nonnegative_integer "CUT_MEAN_QUALITY" "${CUT_MEAN_QUALITY}"

require_boolean "CUT_FRONT" "${CUT_FRONT}"
require_boolean "CUT_TAIL" "${CUT_TAIL}"
require_boolean "TRIM_POLY_G" "${TRIM_POLY_G}"
require_boolean "TRIM_POLY_X" "${TRIM_POLY_X}"

if [[ -n "${ADAPTER_R1}" || -n "${ADAPTER_R2}" ]]; then
    [[ -n "${ADAPTER_R1}" && -n "${ADAPTER_R2}" ]] ||
        die "ADAPTER_R1 and ADAPTER_R2 must be supplied together."
fi

# ---------- Read one manifest row ----------

LINE_NUMBER=$((SLURM_ARRAY_TASK_ID + 1))
LINE="$(sed -n "${LINE_NUMBER}p" "${MANIFEST}")"
LINE="${LINE%$'\r'}"

[[ -n "${LINE}" ]] ||
    die "No manifest entry at line ${LINE_NUMBER}: ${MANIFEST}"

IFS=$'\t' read -r SAMPLE R1 R2 EXTRA <<< "${LINE}"

[[ -n "${SAMPLE}" ]] ||
    die "Missing sample name at manifest line ${LINE_NUMBER}"

[[ -n "${R1}" ]] ||
    die "Missing R1 at manifest line ${LINE_NUMBER}"

[[ -n "${R2}" ]] ||
    die "Missing R2 at manifest line ${LINE_NUMBER}"

[[ -z "${EXTRA:-}" ]] ||
    die "Manifest line ${LINE_NUMBER} contains more than three columns"

# Prevent sample names from accidentally creating paths.
[[ "${SAMPLE}" != */* && "${SAMPLE}" != "." && "${SAMPLE}" != ".." ]] ||
    die "Invalid sample name: ${SAMPLE}"

[[ -f "${R1}" && -r "${R1}" ]] ||
    die "R1 does not exist or is not readable: ${R1}"

[[ -f "${R2}" && -r "${R2}" ]] ||
    die "R2 does not exist or is not readable: ${R2}"

[[ "${R1}" != "${R2}" ]] ||
    die "R1 and R2 point to the same file: ${R1}"

# ---------- Output paths ----------

READS_DIR="${OUTPUT_DIR}/reads"
REPORT_DIR="${OUTPUT_DIR}/reports"
LOG_DIR="${OUTPUT_DIR}/logs"

mkdir -p "${READS_DIR}" "${REPORT_DIR}" "${LOG_DIR}"

OUT_R1="${READS_DIR}/${SAMPLE}_R1.${OUT_SUFFIX}.fastq.gz"
OUT_R2="${READS_DIR}/${SAMPLE}_R2.${OUT_SUFFIX}.fastq.gz"

JSON="${REPORT_DIR}/${SAMPLE}.${OUT_SUFFIX}.fastp.json"
HTML="${REPORT_DIR}/${SAMPLE}.${OUT_SUFFIX}.fastp.html"
FASTP_LOG="${LOG_DIR}/${SAMPLE}.${OUT_SUFFIX}.fastp.log"

# ---------- Resume protection ----------

if [[ -s "${OUT_R1}" &&
      -s "${OUT_R2}" &&
      -s "${JSON}" &&
      -s "${HTML}" ]]
then
    log "fastp outputs already exist; skipping sample: ${SAMPLE}"
    exit 0
fi

# ---------- Construct fastp arguments ----------

FASTP_ARGS=(
    --in1 "${R1}"
    --in2 "${R2}"
    --out1 "${OUT_R1}"
    --out2 "${OUT_R2}"
    --thread "${THREADS}"
    --qualified_quality_phred "${QUALIFIED_QUALITY}"
    --unqualified_percent_limit "${UNQUALIFIED_PERCENT}"
    --length_required "${MIN_LENGTH}"
    --json "${JSON}"
    --html "${HTML}"
)

# Use explicit adapters when supplied; otherwise let fastp detect PE adapters.
if [[ -n "${ADAPTER_R1}" ]]; then
    FASTP_ARGS+=(
        --adapter_sequence "${ADAPTER_R1}"
        --adapter_sequence_r2 "${ADAPTER_R2}"
    )
else
    FASTP_ARGS+=(--detect_adapter_for_pe)
fi

if (( TRIM_FRONT > 0 )); then
    FASTP_ARGS+=(
        --trim_front1 "${TRIM_FRONT}"
        --trim_front2 "${TRIM_FRONT}"
    )
fi

if (( CUT_FRONT == 1 )); then
    FASTP_ARGS+=(
        --cut_front
        --cut_front_window_size "${CUT_WINDOW_SIZE}"
        --cut_front_mean_quality "${CUT_MEAN_QUALITY}"
    )
fi

if (( CUT_TAIL == 1 )); then
    FASTP_ARGS+=(
        --cut_tail
        --cut_tail_window_size "${CUT_WINDOW_SIZE}"
        --cut_tail_mean_quality "${CUT_MEAN_QUALITY}"
    )
fi

if (( TRIM_POLY_G == 1 )); then
    FASTP_ARGS+=(--trim_poly_g)
fi

if (( TRIM_POLY_X == 1 )); then
    FASTP_ARGS+=(--trim_poly_x)
fi

# ---------- Run fastp ----------

log "Starting fastp"
log "Sample     : ${SAMPLE}"
log "R1         : ${R1}"
log "R2         : ${R2}"
log "Output R1  : ${OUT_R1}"
log "Output R2  : ${OUT_R2}"
log "Threads    : ${THREADS}"
log "Trim front: ${TRIM_FRONT}"

START_TIME="${SECONDS}"

if "${FASTP_BIN}" "${FASTP_ARGS[@]}" >"${FASTP_LOG}" 2>&1; then
    :
else
    status=$?
    die "fastp failed with exit status ${status}; see ${FASTP_LOG}"
fi

# ---------- Verify outputs ----------

[[ -s "${OUT_R1}" ]] ||
    die "Trimmed R1 was not created: ${OUT_R1}"

[[ -s "${OUT_R2}" ]] ||
    die "Trimmed R2 was not created: ${OUT_R2}"

[[ -s "${JSON}" ]] ||
    die "fastp JSON report was not created: ${JSON}"

[[ -s "${HTML}" ]] ||
    die "fastp HTML report was not created: ${HTML}"

log "fastp completed: ${SAMPLE}"
log "Elapsed time: $((SECONDS - START_TIME)) seconds"