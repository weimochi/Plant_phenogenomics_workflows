#!/usr/bin/env bash
#SBATCH --job-name=bwa-mem2-map
#SBATCH --partition=CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=2-00:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -Eeuo pipefail

# One paired-end sample per array task.
# MANIFEST format (no header): sample<TAB>R1<TAB>R2

MANIFEST="${MANIFEST:-}"
REFERENCE="${REFERENCE:-}"
OUTPUT_DIR="${OUTPUT_DIR:-}"

BWA_MEM2_BIN="${BWA_MEM2_BIN:-bwa-mem2}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}"

TOTAL_THREADS="${SLURM_CPUS_PER_TASK:-16}"
BWA_THREADS="${BWA_THREADS:-12}"
SORT_THREADS="${SORT_THREADS:-3}"
# samtools sort -m is an approximate limit per sorting thread, not a total.
SORT_MEMORY="${SORT_MEMORY:-2G}"

RG_PLATFORM="${RG_PLATFORM:-ILLUMINA}"
RG_LIBRARY_PREFIX="${RG_LIBRARY_PREFIX:-lib}"
RG_CENTER="${RG_CENTER:-}"

log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2
}

die() {
    log "ERROR: $*"
    exit 1
}

require_positive_integer() {
    local name="$1"
    local value="$2"
    [[ "${value}" =~ ^[1-9][0-9]*$ ]] ||
        die "${name} must be a positive integer: ${value}"
}

PARTIAL_BAM=""
cleanup_partial_output() {
    if [[ -n "${PARTIAL_BAM}" && -f "${PARTIAL_BAM}" ]]; then
        rm -f -- "${PARTIAL_BAM}"
    fi
}
trap cleanup_partial_output EXIT

[[ -n "${MANIFEST}" ]] || die "MANIFEST is required"
[[ -s "${MANIFEST}" ]] || die "MANIFEST does not exist or is empty: ${MANIFEST}"
[[ -n "${REFERENCE}" ]] || die "REFERENCE is required"
[[ -s "${REFERENCE}" ]] || die "Reference FASTA does not exist or is empty: ${REFERENCE}"
[[ -n "${OUTPUT_DIR}" ]] || die "OUTPUT_DIR is required"
[[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || die "Submit this script as a Slurm array job"
[[ "${SLURM_ARRAY_TASK_ID}" =~ ^[0-9]+$ ]] ||
    die "Invalid SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"

require_positive_integer "TOTAL_THREADS" "${TOTAL_THREADS}"
require_positive_integer "BWA_THREADS" "${BWA_THREADS}"
require_positive_integer "SORT_THREADS" "${SORT_THREADS}"

if (( BWA_THREADS + SORT_THREADS + 1 > TOTAL_THREADS )); then
    die "Thread allocation exceeds Slurm CPUs: BWA=${BWA_THREADS}, sort=${SORT_THREADS}, total=${TOTAL_THREADS}"
fi

command -v "${BWA_MEM2_BIN}" >/dev/null 2>&1 ||
    die "BWA-MEM2 executable not found: ${BWA_MEM2_BIN}"
command -v "${SAMTOOLS_BIN}" >/dev/null 2>&1 ||
    die "samtools executable not found: ${SAMTOOLS_BIN}"

# Modern BWA-MEM2 indexes (v2.1+). Build these before submitting the array.
BWA_INDEX_SUFFIXES=(
    ".0123"
    ".amb"
    ".ann"
    ".bwt.2bit.64"
    ".pac"
)
for suffix in "${BWA_INDEX_SUFFIXES[@]}"; do
    [[ -s "${REFERENCE}${suffix}" ]] ||
        die "Missing BWA-MEM2 reference index: ${REFERENCE}${suffix}"
done

LINE_NUMBER=$((SLURM_ARRAY_TASK_ID + 1))
LINE="$(sed -n "${LINE_NUMBER}p" "${MANIFEST}")"
LINE="${LINE%$'\r'}"
[[ -n "${LINE}" ]] ||
    die "No manifest entry at line ${LINE_NUMBER}: ${MANIFEST}"

IFS=$'\t' read -r SAMPLE R1 R2 EXTRA <<< "${LINE}"
[[ -n "${SAMPLE}" ]] || die "Missing sample name at manifest line ${LINE_NUMBER}"
[[ -n "${R1}" ]] || die "Missing R1 at manifest line ${LINE_NUMBER}"
[[ -n "${R2}" ]] || die "Missing R2 at manifest line ${LINE_NUMBER}"
[[ -z "${EXTRA:-}" ]] || die "Manifest line ${LINE_NUMBER} has more than three columns"
[[ "${SAMPLE}" =~ ^[A-Za-z0-9._-]+$ ]] ||
    die "Sample name contains unsupported characters: ${SAMPLE}"
[[ -f "${R1}" && -r "${R1}" ]] || die "R1 is not readable: ${R1}"
[[ -f "${R2}" && -r "${R2}" ]] || die "R2 is not readable: ${R2}"
[[ "${R1}" != "${R2}" ]] || die "R1 and R2 point to the same file"

BAM_DIR="${OUTPUT_DIR}/bam"
LOG_DIR="${OUTPUT_DIR}/logs"
mkdir -p "${BAM_DIR}" "${LOG_DIR}"

FINAL_BAM="${BAM_DIR}/${SAMPLE}.sorted.bam"
FINAL_BAI="${FINAL_BAM}.bai"

# Recover cleanly when mapping succeeded previously but indexing did not.
if [[ -e "${FINAL_BAM}" ]]; then
    [[ -s "${FINAL_BAM}" ]] || die "Existing BAM is empty: ${FINAL_BAM}"
    "${SAMTOOLS_BIN}" quickcheck -v "${FINAL_BAM}" ||
        die "Existing BAM failed quickcheck: ${FINAL_BAM}"

    if [[ -s "${FINAL_BAI}" ]]; then
        log "Valid BAM and index already exist; skipping: ${SAMPLE}"
        exit 0
    fi

    log "Valid BAM exists but its index is missing; rebuilding index"
    "${SAMTOOLS_BIN}" index -@ "${SORT_THREADS}" "${FINAL_BAM}"
    [[ -s "${FINAL_BAI}" ]] || die "Failed to rebuild BAM index: ${FINAL_BAI}"
    log "BAM index rebuilt successfully: ${SAMPLE}"
    exit 0
fi

[[ ! -e "${FINAL_BAI}" ]] ||
    die "BAM index exists without its BAM: ${FINAL_BAI}"

JOB_TAG="${SLURM_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID}"
PARTIAL_BAM="${BAM_DIR}/.${SAMPLE}.${JOB_TAG}.partial.bam"
BWA_LOG="${LOG_DIR}/${SAMPLE}.bwa-mem2.log"
SORT_LOG="${LOG_DIR}/${SAMPLE}.samtools-sort.log"
VERSION_LOG="${LOG_DIR}/${SAMPLE}.tool-versions.log"

TMP_ROOT="${SLURM_TMPDIR:-${OUTPUT_DIR}/tmp}"
TASK_TMP="${TMP_ROOT}/mapping_${JOB_TAG}"
SORT_PREFIX="${TASK_TMP}/${SAMPLE}.sort"
mkdir -p "${TASK_TMP}"

RG_LIBRARY="${RG_LIBRARY_PREFIX}.${SAMPLE}"
READ_GROUP="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${RG_LIBRARY}\tPL:${RG_PLATFORM}"
if [[ -n "${RG_CENTER}" ]]; then
    READ_GROUP+="\tCN:${RG_CENTER}"
fi

{
    printf 'sample\t%s\n' "${SAMPLE}"
    printf 'reference\t%s\n' "${REFERENCE}"
    printf 'r1\t%s\n' "${R1}"
    printf 'r2\t%s\n' "${R2}"
    printf 'bwa_threads\t%s\n' "${BWA_THREADS}"
    printf 'sort_threads\t%s\n' "${SORT_THREADS}"
    printf 'sort_memory_per_thread\t%s\n' "${SORT_MEMORY}"
    printf '\n[bwa-mem2]\n'
    "${BWA_MEM2_BIN}" version 2>&1 || true
    printf '\n[samtools]\n'
    "${SAMTOOLS_BIN}" --version 2>&1 || true
} > "${VERSION_LOG}"

log "Starting BWA-MEM2 mapping: ${SAMPLE}"
log "Reference: ${REFERENCE}"
log "Output BAM: ${FINAL_BAM}"
START_TIME="${SECONDS}"

if "${BWA_MEM2_BIN}" mem \
    -t "${BWA_THREADS}" \
    -R "${READ_GROUP}" \
    "${REFERENCE}" "${R1}" "${R2}" \
    2> "${BWA_LOG}" |
    "${SAMTOOLS_BIN}" sort \
        -@ "${SORT_THREADS}" \
        -m "${SORT_MEMORY}" \
        -T "${SORT_PREFIX}" \
        -O BAM \
        -o "${PARTIAL_BAM}" \
        - \
        2> "${SORT_LOG}"
then
    :
else
    status=$?
    die "BWA-MEM2/samtools pipeline failed with exit status ${status}"
fi

[[ -s "${PARTIAL_BAM}" ]] || die "Mapping did not create a BAM: ${PARTIAL_BAM}"
"${SAMTOOLS_BIN}" quickcheck -v "${PARTIAL_BAM}" ||
    die "BAM failed quickcheck: ${PARTIAL_BAM}"

mv -- "${PARTIAL_BAM}" "${FINAL_BAM}"
PARTIAL_BAM=""

"${SAMTOOLS_BIN}" index -@ "${SORT_THREADS}" "${FINAL_BAM}"
[[ -s "${FINAL_BAI}" ]] || die "BAM index was not created: ${FINAL_BAI}"

log "Mapping completed: ${SAMPLE}"
log "Elapsed time: $((SECONDS - START_TIME)) seconds"
