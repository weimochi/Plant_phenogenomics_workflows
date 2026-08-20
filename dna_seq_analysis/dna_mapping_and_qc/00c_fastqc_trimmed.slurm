#!/usr/bin/env bash
#SBATCH --job-name=fastqc-trimmed
#SBATCH --partition=CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -Eeuo pipefail

# ==============================================================================
# 00c_fastqc_trimmed.slurm
# Run FastQC on one trimmed FASTQ file per Slurm array task.
#
# FILELIST format, without header:
# /path/to/sample01_R1.trimmed.fastq.gz
# /path/to/sample01_R2.trimmed.fastq.gz
# ==============================================================================

# ---------- Required inputs ----------

FILELIST="${FILELIST:-}"
OUTPUT_DIR="${OUTPUT_DIR:-}"

# ---------- Executable and resources ----------

FASTQC_BIN="${FASTQC_BIN:-fastqc}"
THREADS="${SLURM_CPUS_PER_TASK:-1}"

# ---------- Logging ----------

log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2
}

die() {
    log "ERROR: $*"
    exit 1
}

# ---------- Validation ----------

[[ -n "${FILELIST}" ]] ||
    die "FILELIST is required."

[[ -s "${FILELIST}" ]] ||
    die "FILELIST does not exist or is empty: ${FILELIST}"

[[ -n "${OUTPUT_DIR}" ]] ||
    die "OUTPUT_DIR is required."

[[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] ||
    die "This script must be submitted as a Slurm array job."

[[ "${SLURM_ARRAY_TASK_ID}" =~ ^[0-9]+$ ]] ||
    die "Invalid SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"

command -v "${FASTQC_BIN}" >/dev/null 2>&1 ||
    die "FastQC executable not found: ${FASTQC_BIN}"

# ---------- Select one trimmed FASTQ ----------

LINE_NUMBER=$((SLURM_ARRAY_TASK_ID + 1))
FASTQ="$(sed -n "${LINE_NUMBER}p" "${FILELIST}")"

# Remove a possible carriage return from Windows-formatted file lists.
FASTQ="${FASTQ%$'\r'}"

[[ -n "${FASTQ}" ]] ||
    die "No FASTQ found at line ${LINE_NUMBER} of ${FILELIST}"

[[ -f "${FASTQ}" ]] ||
    die "FASTQ does not exist: ${FASTQ}"

[[ -r "${FASTQ}" ]] ||
    die "FASTQ is not readable: ${FASTQ}"

# ---------- Determine filename prefix ----------

FASTQ_BASENAME="$(basename "${FASTQ}")"

case "${FASTQ_BASENAME}" in
    *.fastq.gz)
        PREFIX="${FASTQ_BASENAME%.fastq.gz}"
        ;;
    *.fq.gz)
        PREFIX="${FASTQ_BASENAME%.fq.gz}"
        ;;
    *.fastq)
        PREFIX="${FASTQ_BASENAME%.fastq}"
        ;;
    *.fq)
        PREFIX="${FASTQ_BASENAME%.fq}"
        ;;
    *)
        die "Unsupported FASTQ filename: ${FASTQ_BASENAME}"
        ;;
esac

# ---------- Output directories ----------

FASTQC_DIR="${OUTPUT_DIR}/fastqc"
LOG_DIR="${OUTPUT_DIR}/logs"
TMP_ROOT="${SLURM_TMPDIR:-${OUTPUT_DIR}/tmp}"
TASK_TMP="${TMP_ROOT}/fastqc_${SLURM_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID}"

mkdir -p \
    "${FASTQC_DIR}" \
    "${LOG_DIR}" \
    "${TASK_TMP}"

EXPECTED_HTML="${FASTQC_DIR}/${PREFIX}_fastqc.html"
EXPECTED_ZIP="${FASTQC_DIR}/${PREFIX}_fastqc.zip"
FASTQC_LOG="${LOG_DIR}/${PREFIX}.fastqc.log"

# ---------- Resume protection ----------

if [[ -s "${EXPECTED_HTML}" && -s "${EXPECTED_ZIP}" ]]; then
    log "FastQC outputs already exist; skipping: ${FASTQ_BASENAME}"
    exit 0
fi

# ---------- Run FastQC ----------

log "Starting trimmed-read FastQC"
log "Array task: ${SLURM_ARRAY_TASK_ID}"
log "Input     : ${FASTQ}"
log "Output    : ${FASTQC_DIR}"
log "Threads   : ${THREADS}"

START_TIME="${SECONDS}"

if "${FASTQC_BIN}" \
    --threads "${THREADS}" \
    --noextract \
    --dir "${TASK_TMP}" \
    --outdir "${FASTQC_DIR}" \
    "${FASTQ}" \
    2>&1 | tee "${FASTQC_LOG}"
then
    :
else
    status=$?
    die "FastQC failed with exit status ${status}: ${FASTQ_BASENAME}"
fi

# ---------- Verify outputs ----------

[[ -s "${EXPECTED_HTML}" ]] ||
    die "Expected FastQC HTML was not created: ${EXPECTED_HTML}"

[[ -s "${EXPECTED_ZIP}" ]] ||
    die "Expected FastQC ZIP was not created: ${EXPECTED_ZIP}"

log "Trimmed-read FastQC completed: ${FASTQ_BASENAME}"
log "Elapsed time: $((SECONDS - START_TIME)) seconds"