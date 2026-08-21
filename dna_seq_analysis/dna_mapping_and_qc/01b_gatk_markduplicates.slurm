#!/usr/bin/env bash
#SBATCH --job-name=gatk-markdup
#SBATCH --partition=CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=1-00:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -Eeuo pipefail

# ==============================================================================
# 01b_gatk_markduplicates.slurm
#
# One BAM per Slurm array task.
#
# MANIFEST format, without header:
# sample<TAB>coordinate_sorted_bam
# ==============================================================================

MANIFEST="${MANIFEST:-}"
OUTPUT_DIR="${OUTPUT_DIR:-}"

GATK_BIN="${GATK_BIN:-gatk}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}"

THREADS="${SLURM_CPUS_PER_TASK:-4}"
JAVA_HEAP="${JAVA_HEAP:-40g}"

# false: mark duplicates but retain reads
# true : physically remove duplicate reads
REMOVE_DUPLICATES="${REMOVE_DUPLICATES:-false}"

# Configure this according to the sequencing platform if optical duplicate
# detection is important. Picard/GATK traditionally defaults to 100.
OPTICAL_DUPLICATE_PIXEL_DISTANCE="${OPTICAL_DUPLICATE_PIXEL_DISTANCE:-100}"

log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2
}

die() {
    log "ERROR: $*"
    exit 1
}

PARTIAL_BAM=""
PARTIAL_METRICS=""

cleanup_partial_outputs() {
    if [[ -n "${PARTIAL_BAM}" && -f "${PARTIAL_BAM}" ]]; then
        rm -f -- "${PARTIAL_BAM}"
    fi

    if [[ -n "${PARTIAL_METRICS}" && -f "${PARTIAL_METRICS}" ]]; then
        rm -f -- "${PARTIAL_METRICS}"
    fi
}

trap cleanup_partial_outputs EXIT

# ---------- Validate configuration ----------

[[ -n "${MANIFEST}" ]] ||
    die "MANIFEST is required"

[[ -s "${MANIFEST}" ]] ||
    die "MANIFEST does not exist or is empty: ${MANIFEST}"

[[ -n "${OUTPUT_DIR}" ]] ||
    die "OUTPUT_DIR is required"

[[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] ||
    die "Submit this script as a Slurm array job"

[[ "${SLURM_ARRAY_TASK_ID}" =~ ^[0-9]+$ ]] ||
    die "Invalid SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"

[[ "${THREADS}" =~ ^[1-9][0-9]*$ ]] ||
    die "THREADS must be a positive integer: ${THREADS}"

[[ "${OPTICAL_DUPLICATE_PIXEL_DISTANCE}" =~ ^[0-9]+$ ]] ||
    die "OPTICAL_DUPLICATE_PIXEL_DISTANCE must be a non-negative integer"

case "${REMOVE_DUPLICATES}" in
    true|false)
        ;;
    *)
        die "REMOVE_DUPLICATES must be true or false"
        ;;
esac

command -v "${GATK_BIN}" >/dev/null 2>&1 ||
    die "GATK executable not found: ${GATK_BIN}"

command -v "${SAMTOOLS_BIN}" >/dev/null 2>&1 ||
    die "samtools executable not found: ${SAMTOOLS_BIN}"

# ---------- Read one manifest row ----------

LINE_NUMBER=$((SLURM_ARRAY_TASK_ID + 1))
LINE="$(sed -n "${LINE_NUMBER}p" "${MANIFEST}")"
LINE="${LINE%$'\r'}"

[[ -n "${LINE}" ]] ||
    die "No manifest entry at line ${LINE_NUMBER}: ${MANIFEST}"

IFS=$'\t' read -r SAMPLE INPUT_BAM EXTRA <<< "${LINE}"

[[ -n "${SAMPLE}" ]] ||
    die "Missing sample name at manifest line ${LINE_NUMBER}"

[[ -n "${INPUT_BAM}" ]] ||
    die "Missing BAM path at manifest line ${LINE_NUMBER}"

[[ -z "${EXTRA:-}" ]] ||
    die "Manifest line ${LINE_NUMBER} has more than two columns"

[[ "${SAMPLE}" =~ ^[A-Za-z0-9._-]+$ ]] ||
    die "Sample name contains unsupported characters: ${SAMPLE}"

[[ -s "${INPUT_BAM}" && -r "${INPUT_BAM}" ]] ||
    die "Input BAM does not exist or is not readable: ${INPUT_BAM}"

"${SAMTOOLS_BIN}" quickcheck -v "${INPUT_BAM}" ||
    die "Input BAM failed quickcheck: ${INPUT_BAM}"

# Require coordinate-sorted input.
SORT_ORDER="$(
    "${SAMTOOLS_BIN}" view -H "${INPUT_BAM}" |
        awk '
            /^@HD/ {
                for (i = 1; i <= NF; i++) {
                    if ($i ~ /^SO:/) {
                        sub(/^SO:/, "", $i)
                        print $i
                        exit
                    }
                }
            }
        '
)"

[[ "${SORT_ORDER}" == "coordinate" ]] ||
    die "Input BAM is not declared coordinate-sorted: ${INPUT_BAM}"

# ---------- Output paths ----------

BAM_DIR="${OUTPUT_DIR}/bam"
METRICS_DIR="${OUTPUT_DIR}/metrics"
LOG_DIR="${OUTPUT_DIR}/logs"

mkdir -p "${BAM_DIR}" "${METRICS_DIR}" "${LOG_DIR}"

FINAL_BAM="${BAM_DIR}/${SAMPLE}.markdup.bam"
FINAL_BAI="${FINAL_BAM}.bai"
FINAL_METRICS="${METRICS_DIR}/${SAMPLE}.duplicate_metrics.txt"

# ---------- Resume and recovery ----------

if [[ -e "${FINAL_BAM}" ]]; then
    [[ -s "${FINAL_BAM}" ]] ||
        die "Existing output BAM is empty: ${FINAL_BAM}"

    "${SAMTOOLS_BIN}" quickcheck -v "${FINAL_BAM}" ||
        die "Existing output BAM failed quickcheck: ${FINAL_BAM}"

    [[ -s "${FINAL_METRICS}" ]] ||
        die "Valid markdup BAM exists, but metrics are missing: ${FINAL_METRICS}"

    if [[ ! -s "${FINAL_BAI}" ]]; then
        log "Valid markdup BAM exists but its index is missing; rebuilding index"

        "${SAMTOOLS_BIN}" index \
            -@ "${THREADS}" \
            "${FINAL_BAM}"

        [[ -s "${FINAL_BAI}" ]] ||
            die "Failed to rebuild BAM index: ${FINAL_BAI}"
    fi

    log "Valid MarkDuplicates outputs already exist; skipping: ${SAMPLE}"
    exit 0
fi

[[ ! -e "${FINAL_BAI}" ]] ||
    die "BAM index exists without its BAM: ${FINAL_BAI}"

[[ ! -e "${FINAL_METRICS}" ]] ||
    die "Duplicate metrics exist without their BAM: ${FINAL_METRICS}"

# ---------- Temporary paths ----------

JOB_TAG="${SLURM_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID}"

PARTIAL_BAM="${BAM_DIR}/.${SAMPLE}.${JOB_TAG}.partial.bam"
PARTIAL_METRICS="${METRICS_DIR}/.${SAMPLE}.${JOB_TAG}.partial.duplicate_metrics.txt"

TMP_ROOT="${SLURM_TMPDIR:-${OUTPUT_DIR}/tmp}"
TASK_TMP="${TMP_ROOT}/markduplicates_${JOB_TAG}"

mkdir -p "${TASK_TMP}"

GATK_LOG="${LOG_DIR}/${SAMPLE}.gatk-markduplicates.log"
VERSION_LOG="${LOG_DIR}/${SAMPLE}.tool-versions.log"

# ---------- Record versions and configuration ----------

{
    printf 'sample\t%s\n' "${SAMPLE}"
    printf 'input_bam\t%s\n' "${INPUT_BAM}"
    printf 'remove_duplicates\t%s\n' "${REMOVE_DUPLICATES}"
    printf 'java_heap\t%s\n' "${JAVA_HEAP}"
    printf 'optical_duplicate_pixel_distance\t%s\n' \
        "${OPTICAL_DUPLICATE_PIXEL_DISTANCE}"

    printf '\n[gatk]\n'
    "${GATK_BIN}" --version 2>&1 || true

    printf '\n[samtools]\n'
    "${SAMTOOLS_BIN}" --version 2>&1 || true
} > "${VERSION_LOG}"

# ---------- Run GATK MarkDuplicates ----------

log "Starting GATK MarkDuplicates"
log "Sample            : ${SAMPLE}"
log "Input BAM         : ${INPUT_BAM}"
log "Output BAM        : ${FINAL_BAM}"
log "Remove duplicates : ${REMOVE_DUPLICATES}"
log "Java heap         : ${JAVA_HEAP}"

START_TIME="${SECONDS}"

if "${GATK_BIN}" \
    --java-options "-Xmx${JAVA_HEAP} -Djava.io.tmpdir=${TASK_TMP}" \
    MarkDuplicates \
    --INPUT "${INPUT_BAM}" \
    --OUTPUT "${PARTIAL_BAM}" \
    --METRICS_FILE "${PARTIAL_METRICS}" \
    --REMOVE_DUPLICATES "${REMOVE_DUPLICATES}" \
    --CREATE_INDEX false \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE \
        "${OPTICAL_DUPLICATE_PIXEL_DISTANCE}" \
    --TMP_DIR "${TASK_TMP}" \
    > "${GATK_LOG}" 2>&1
then
    :
else
    status=$?
    die "GATK MarkDuplicates failed with exit status ${status}; see ${GATK_LOG}"
fi

# ---------- Validate and publish ----------

[[ -s "${PARTIAL_BAM}" ]] ||
    die "MarkDuplicates did not create a BAM: ${PARTIAL_BAM}"

[[ -s "${PARTIAL_METRICS}" ]] ||
    die "MarkDuplicates did not create metrics: ${PARTIAL_METRICS}"

"${SAMTOOLS_BIN}" quickcheck -v "${PARTIAL_BAM}" ||
    die "MarkDuplicates BAM failed quickcheck: ${PARTIAL_BAM}"

mv -- "${PARTIAL_METRICS}" "${FINAL_METRICS}"
PARTIAL_METRICS=""

mv -- "${PARTIAL_BAM}" "${FINAL_BAM}"
PARTIAL_BAM=""

"${SAMTOOLS_BIN}" index \
    -@ "${THREADS}" \
    "${FINAL_BAM}"

[[ -s "${FINAL_BAI}" ]] ||
    die "BAM index was not created: ${FINAL_BAI}"

log "GATK MarkDuplicates completed: ${SAMPLE}"
log "Elapsed time: $((SECONDS - START_TIME)) seconds"
