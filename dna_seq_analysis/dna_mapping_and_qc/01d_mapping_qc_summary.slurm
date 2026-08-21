#!/usr/bin/env bash
#SBATCH --job-name=mapping-qc-summary
#SBATCH --partition=CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -Eeuo pipefail

# ==============================================================================
# 01d_mapping_qc_summary.slurm
#
# Combine per-sample mapping QC metrics produced by:
#   01c_mapping_qc.slurm
#
# Expected input structure:
# QC_ROOT/
# ├── sample01/sample01.mapping_qc_metrics.tsv
# ├── sample02/sample02.mapping_qc_metrics.tsv
# └── sample03/sample03.mapping_qc_metrics.tsv
# ==============================================================================

QC_ROOT="${QC_ROOT:-}"
OUTPUT_TSV="${OUTPUT_TSV:-}"

log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2
}

die() {
    log "ERROR: $*"
    exit 1
}

PARTIAL_OUTPUT=""
TEMP_ROWS=""

cleanup_partial_outputs() {
    if [[ -n "${PARTIAL_OUTPUT}" && -f "${PARTIAL_OUTPUT}" ]]; then
        rm -f -- "${PARTIAL_OUTPUT}"
    fi

    if [[ -n "${TEMP_ROWS}" && -f "${TEMP_ROWS}" ]]; then
        rm -f -- "${TEMP_ROWS}"
    fi
}

trap cleanup_partial_outputs EXIT

# ---------- Validate configuration ----------

[[ -n "${QC_ROOT}" ]] ||
    die "QC_ROOT is required"

[[ -d "${QC_ROOT}" ]] ||
    die "QC_ROOT does not exist: ${QC_ROOT}"

if [[ -z "${OUTPUT_TSV}" ]]; then
    OUTPUT_TSV="${QC_ROOT}/mapping_qc_summary.tsv"
fi

SUMMARY_DIR="$(dirname "${OUTPUT_TSV}")"
mkdir -p "${SUMMARY_DIR}"

[[ ! -d "${OUTPUT_TSV}" ]] ||
    die "OUTPUT_TSV points to a directory: ${OUTPUT_TSV}"

# ---------- Collect per-sample metrics ----------

mapfile -d '' -t METRIC_FILES < <(
    find "${QC_ROOT}" \
        -mindepth 2 \
        -maxdepth 2 \
        -type f \
        -name '*.mapping_qc_metrics.tsv' \
        -print0 |
        sort -z
)

(( ${#METRIC_FILES[@]} > 0 )) ||
    die "No per-sample mapping QC metrics found under: ${QC_ROOT}"

log "Found ${#METRIC_FILES[@]} per-sample metrics files"

# ---------- Temporary output paths ----------

JOB_TAG="${SLURM_JOB_ID:-manual}"
PARTIAL_OUTPUT="${SUMMARY_DIR}/.mapping_qc_summary.${JOB_TAG}.partial.tsv"
TEMP_ROWS="${SUMMARY_DIR}/.mapping_qc_rows.${JOB_TAG}.tmp.tsv"

: > "${TEMP_ROWS}"

# ---------- Validate and collect rows ----------

EXPECTED_HEADER=""
declare -A SEEN_SAMPLES=()

for metrics_file in "${METRIC_FILES[@]}"; do
    log "Reading: ${metrics_file}"

    # Each per-sample metrics file must contain exactly:
    #   one header row
    #   one sample row
    NONEMPTY_LINES="$(
        awk 'NF > 0 { count++ } END { print count + 0 }' "${metrics_file}"
    )"

    [[ "${NONEMPTY_LINES}" -eq 2 ]] ||
        die "Expected exactly two non-empty lines: ${metrics_file}"

    HEADER="$(awk 'NF > 0 { print; exit }' "${metrics_file}")"
    ROW="$(awk 'NF > 0 { count++; if (count == 2) print }' "${metrics_file}")"

    [[ -n "${HEADER}" ]] ||
        die "Missing header: ${metrics_file}"

    [[ -n "${ROW}" ]] ||
        die "Missing sample row: ${metrics_file}"

    if [[ -z "${EXPECTED_HEADER}" ]]; then
        EXPECTED_HEADER="${HEADER}"

        FIRST_COLUMN="${EXPECTED_HEADER%%$'\t'*}"
        [[ "${FIRST_COLUMN}" == "sample" ]] ||
            die "First summary column must be named 'sample': ${metrics_file}"
    elif [[ "${HEADER}" != "${EXPECTED_HEADER}" ]]; then
        die "Header differs from other metrics files: ${metrics_file}"
    fi

    SAMPLE="${ROW%%$'\t'*}"

    [[ -n "${SAMPLE}" ]] ||
        die "Sample field is empty: ${metrics_file}"

    [[ "${SAMPLE}" =~ ^[A-Za-z0-9._-]+$ ]] ||
        die "Invalid sample name in metrics file: ${SAMPLE}"

    if [[ -n "${SEEN_SAMPLES[${SAMPLE}]:-}" ]]; then
        die "Duplicate sample in QC metrics: ${SAMPLE}"
    fi

    SEEN_SAMPLES["${SAMPLE}"]=1

    printf '%s\n' "${ROW}" >> "${TEMP_ROWS}"
done

# ---------- Build sorted summary ----------

{
    printf '%s\n' "${EXPECTED_HEADER}"
    sort -t $'\t' -k1,1 "${TEMP_ROWS}"
} > "${PARTIAL_OUTPUT}"

EXPECTED_ROWS=$(( ${#METRIC_FILES[@]} + 1 ))

ACTUAL_ROWS="$(
    awk 'NF > 0 { count++ } END { print count + 0 }' "${PARTIAL_OUTPUT}"
)"

[[ "${ACTUAL_ROWS}" -eq "${EXPECTED_ROWS}" ]] ||
    die "Summary row count is incorrect: expected ${EXPECTED_ROWS}, got ${ACTUAL_ROWS}"

[[ -s "${PARTIAL_OUTPUT}" ]] ||
    die "Summary output is empty: ${PARTIAL_OUTPUT}"

# ---------- Atomically publish summary ----------

mv -- "${PARTIAL_OUTPUT}" "${OUTPUT_TSV}"
PARTIAL_OUTPUT=""

rm -f -- "${TEMP_ROWS}"
TEMP_ROWS=""

log "Mapping QC summary completed"
log "Samples: ${#METRIC_FILES[@]}"
log "Output : ${OUTPUT_TSV}"