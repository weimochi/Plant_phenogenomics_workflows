#!/usr/bin/env bash
#SBATCH --job-name=glnexus-joint
#SBATCH --partition=CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=2-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -Eeuo pipefail
export LC_ALL=C

# ==============================================================================
# 02_glnexus_joint_calling.slurm
#
# Joint-call all completed DeepVariant gVCFs in one non-array Slurm job.
# Expected input layout:
#   DEEPVARIANT_OUTPUT_DIR/samples/<sample>/<sample>.deepvariant.g.vcf.gz
# ==============================================================================

DEEPVARIANT_OUTPUT_DIR="${DEEPVARIANT_OUTPUT_DIR:-}"
OUTPUT_DIR="${OUTPUT_DIR:-}"
COHORT_NAME="${COHORT_NAME:-cohort}"

GVCF_PATTERN="${GVCF_PATTERN:-*.deepvariant.g.vcf.gz}"
MIN_SAMPLES="${MIN_SAMPLES:-2}"

BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

GLNEXUS_CONFIG="${GLNEXUS_CONFIG:-DeepVariantWGS}"
GLNEXUS_VERSION="${GLNEXUS_VERSION:-1.4.1}"
GLNEXUS_IMAGE="${GLNEXUS_IMAGE:-docker://ghcr.io/dnanexus-rnd/glnexus:v${GLNEXUS_VERSION}}"
GLNEXUS_BIN_IN_CONTAINER="${GLNEXUS_BIN_IN_CONTAINER:-/usr/local/bin/glnexus_cli}"
CONTAINER_RUNTIME="${CONTAINER_RUNTIME:-}"
EXTRA_BINDS="${EXTRA_BINDS:-}"

# RocksDB performs heavy random I/O. Use node-local scratch whenever possible.
# If the cluster does not define SLURM_TMPDIR, GLNEXUS_DB_ROOT is required.
GLNEXUS_DB_ROOT="${GLNEXUS_DB_ROOT:-${SLURM_TMPDIR:-}}"

log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2
}

die() {
    log "ERROR: $*"
    exit 1
}

STAGING_DIR=""
DB_DIR=""

cleanup_partial_outputs() {
    if [[ -n "${STAGING_DIR}" && -d "${STAGING_DIR}" ]]; then
        rm -rf -- "${STAGING_DIR}"
    fi

    if [[ -n "${DB_DIR}" &&
          -d "${DB_DIR}" &&
          "${DB_DIR}" == */glnexus_db_* ]]
    then
        rm -rf -- "${DB_DIR}"
    fi
}

trap cleanup_partial_outputs EXIT

has_variant_index() {
    local variant_file="$1"
    [[ -s "${variant_file}.tbi" || -s "${variant_file}.csi" ]]
}

# ---------- Validate configuration ----------

[[ -n "${DEEPVARIANT_OUTPUT_DIR}" ]] || die "DEEPVARIANT_OUTPUT_DIR is required"
[[ -d "${DEEPVARIANT_OUTPUT_DIR}" ]] ||
    die "DEEPVARIANT_OUTPUT_DIR does not exist: ${DEEPVARIANT_OUTPUT_DIR}"
[[ -n "${OUTPUT_DIR}" ]] || die "OUTPUT_DIR is required"
[[ -n "${GLNEXUS_DB_ROOT}" ]] ||
    die "Set GLNEXUS_DB_ROOT to a sufficiently large scratch directory"
[[ -d "${GLNEXUS_DB_ROOT}" ]] ||
    die "GLNEXUS_DB_ROOT does not exist: ${GLNEXUS_DB_ROOT}"
[[ -w "${GLNEXUS_DB_ROOT}" ]] ||
    die "GLNEXUS_DB_ROOT is not writable: ${GLNEXUS_DB_ROOT}"
[[ "${THREADS}" =~ ^[1-9][0-9]*$ ]] ||
    die "THREADS must be a positive integer: ${THREADS}"
[[ "${MIN_SAMPLES}" =~ ^[1-9][0-9]*$ ]] ||
    die "MIN_SAMPLES must be a positive integer: ${MIN_SAMPLES}"
[[ "${COHORT_NAME}" =~ ^[A-Za-z0-9._-]+$ ]] ||
    die "COHORT_NAME contains unsupported characters: ${COHORT_NAME}"

case "${GLNEXUS_CONFIG}" in
    DeepVariant|DeepVariantWGS|DeepVariantWES|DeepVariant_unfiltered)
        ;;
    *)
        die "Unsupported GLNEXUS_CONFIG: ${GLNEXUS_CONFIG}"
        ;;
esac

command -v "${BCFTOOLS_BIN}" >/dev/null 2>&1 ||
    die "bcftools executable not found: ${BCFTOOLS_BIN}"

DEEPVARIANT_OUTPUT_DIR="$(cd "${DEEPVARIANT_OUTPUT_DIR}" && pwd -P)"
mkdir -p "${OUTPUT_DIR}"
OUTPUT_DIR="$(cd "${OUTPUT_DIR}" && pwd -P)"
GLNEXUS_DB_ROOT="$(cd "${GLNEXUS_DB_ROOT}" && pwd -P)"

# ---------- Locate Apptainer or Singularity ----------

if [[ -z "${CONTAINER_RUNTIME}" ]]; then
    if command -v apptainer >/dev/null 2>&1; then
        CONTAINER_RUNTIME="apptainer"
    elif command -v singularity >/dev/null 2>&1; then
        CONTAINER_RUNTIME="singularity"
    else
        if command -v module >/dev/null 2>&1; then
            module load apptainer 2>/dev/null ||
                module load singularity/3.10.4 2>/dev/null ||
                true
        fi

        if command -v apptainer >/dev/null 2>&1; then
            CONTAINER_RUNTIME="apptainer"
        elif command -v singularity >/dev/null 2>&1; then
            CONTAINER_RUNTIME="singularity"
        else
            die "Neither Apptainer nor Singularity was found"
        fi
    fi
fi

command -v "${CONTAINER_RUNTIME}" >/dev/null 2>&1 ||
    die "Container runtime not found: ${CONTAINER_RUNTIME}"

if [[ "${GLNEXUS_IMAGE}" != *"://"* ]]; then
    [[ -s "${GLNEXUS_IMAGE}" ]] || die "GLnexus image does not exist: ${GLNEXUS_IMAGE}"
fi

# ---------- Discover gVCFs ----------

GVCF_ROOT="${DEEPVARIANT_OUTPUT_DIR}/samples"
[[ -d "${GVCF_ROOT}" ]] || die "DeepVariant samples directory is missing: ${GVCF_ROOT}"

mapfile -d '' -t GVCF_FILES < <(
    find "${GVCF_ROOT}" \
        -mindepth 2 \
        -maxdepth 2 \
        -type f \
        -name "${GVCF_PATTERN}" \
        -print0 |
        sort -z
)

N_GVCFS="${#GVCF_FILES[@]}"
(( N_GVCFS >= MIN_SAMPLES )) ||
    die "Found ${N_GVCFS} gVCFs; at least ${MIN_SAMPLES} are required"

# ---------- Output and staging paths ----------

FINAL_DIR="${OUTPUT_DIR}/${COHORT_NAME}"

FINAL_BCF="${FINAL_DIR}/${COHORT_NAME}.deepvariant.glnexus.bcf"
FINAL_VCF="${FINAL_DIR}/${COHORT_NAME}.deepvariant.glnexus.vcf.gz"
FINAL_STATS="${FINAL_DIR}/${COHORT_NAME}.deepvariant.glnexus.bcftools.stats.txt"

if [[ -d "${FINAL_DIR}" ]]; then
    if [[ -s "${FINAL_BCF}" &&
          -s "${FINAL_BCF}.csi" &&
          -s "${FINAL_VCF}" &&
          -s "${FINAL_STATS}" ]] &&
       has_variant_index "${FINAL_VCF}"
    then
        "${BCFTOOLS_BIN}" view --header-only "${FINAL_BCF}" >/dev/null ||
            die "Existing cohort BCF failed validation: ${FINAL_BCF}"
        "${BCFTOOLS_BIN}" view --header-only "${FINAL_VCF}" >/dev/null ||
            die "Existing cohort VCF failed validation: ${FINAL_VCF}"
        log "Valid GLnexus cohort outputs already exist; skipping: ${COHORT_NAME}"
        exit 0
    fi

    die "Incomplete existing cohort directory: ${FINAL_DIR}"
fi

JOB_TAG="${SLURM_JOB_ID:-manual}"
STAGING_DIR="${OUTPUT_DIR}/.${COHORT_NAME}.${JOB_TAG}.partial"
[[ ! -e "${STAGING_DIR}" ]] || die "Staging directory already exists: ${STAGING_DIR}"
mkdir -p "${STAGING_DIR}"

GVCF_LIST="${STAGING_DIR}/${COHORT_NAME}.gvcf.list"
SAMPLE_LIST="${STAGING_DIR}/${COHORT_NAME}.samples.txt"
BCF="${STAGING_DIR}/${COHORT_NAME}.deepvariant.glnexus.bcf"
VCF="${STAGING_DIR}/${COHORT_NAME}.deepvariant.glnexus.vcf.gz"
STATS="${STAGING_DIR}/${COHORT_NAME}.deepvariant.glnexus.bcftools.stats.txt"
CONFIG_LOG="${STAGING_DIR}/${COHORT_NAME}.glnexus.config.tsv"
GLNEXUS_LOG="${STAGING_DIR}/${COHORT_NAME}.glnexus.log"

DB_DIR="${GLNEXUS_DB_ROOT}/glnexus_db_${JOB_TAG}"
[[ ! -e "${DB_DIR}" ]] || die "GLnexus database directory already exists: ${DB_DIR}"
mkdir -p "${DB_DIR}"

# ---------- Validate gVCFs and sample identities ----------

declare -A SEEN_SAMPLES=()
: > "${GVCF_LIST}"
: > "${SAMPLE_LIST}"

for gvcf in "${GVCF_FILES[@]}"; do
    [[ -s "${gvcf}" ]] || die "gVCF is empty: ${gvcf}"
    has_variant_index "${gvcf}" || die "gVCF index is missing: ${gvcf}"

    "${BCFTOOLS_BIN}" view --header-only "${gvcf}" >/dev/null ||
        die "gVCF failed validation: ${gvcf}"
    "${BCFTOOLS_BIN}" index --stats "${gvcf}" >/dev/null ||
        die "gVCF index failed validation: ${gvcf}"

    mapfile -t GVCF_SAMPLES < <("${BCFTOOLS_BIN}" query -l "${gvcf}")
    [[ "${#GVCF_SAMPLES[@]}" -eq 1 ]] ||
        die "Expected exactly one sample in gVCF: ${gvcf}"

    sample="${GVCF_SAMPLES[0]}"
    [[ -n "${sample}" ]] || die "gVCF has an empty sample name: ${gvcf}"

    if [[ -n "${SEEN_SAMPLES[${sample}]:-}" ]]; then
        die "Duplicate sample '${sample}' found in multiple gVCFs"
    fi

    SEEN_SAMPLES["${sample}"]=1
    printf '%s\n' "${gvcf}" >> "${GVCF_LIST}"
    printf '%s\n' "${sample}" >> "${SAMPLE_LIST}"
done

[[ "$(wc -l < "${GVCF_LIST}")" -eq "${N_GVCFS}" ]] ||
    die "gVCF list row count does not match discovered inputs"

# ---------- Construct container command ----------

STAGING_ABS="$(cd "${STAGING_DIR}" && pwd -P)"
DB_ABS="$(cd "${DB_DIR}" && pwd -P)"

CONTAINER_ARGS=(
    exec
    --cleanenv
    --bind "${GVCF_ROOT}:${GVCF_ROOT}"
    --bind "${STAGING_ABS}:${STAGING_ABS}"
    --bind "${DB_ABS}:${DB_ABS}"
)

if [[ -n "${EXTRA_BINDS}" ]]; then
    CONTAINER_ARGS+=(--bind "${EXTRA_BINDS}")
fi

{
    printf 'cohort_name\t%s\n' "${COHORT_NAME}"
    printf 'sample_count\t%s\n' "${N_GVCFS}"
    printf 'glnexus_config\t%s\n' "${GLNEXUS_CONFIG}"
    printf 'glnexus_version\t%s\n' "${GLNEXUS_VERSION}"
    printf 'glnexus_image\t%s\n' "${GLNEXUS_IMAGE}"
    printf 'container_runtime\t%s\n' "${CONTAINER_RUNTIME}"
    printf 'threads\t%s\n' "${THREADS}"
    printf 'glnexus_db_root\t%s\n' "${GLNEXUS_DB_ROOT}"
    printf '\n[db_filesystem]\n'
    df -h "${GLNEXUS_DB_ROOT}" 2>&1 || true
} > "${CONFIG_LOG}"

# ---------- Run GLnexus ----------

log "Starting GLnexus joint calling"
log "Cohort : ${COHORT_NAME}"
log "Samples: ${N_GVCFS}"
log "Preset : ${GLNEXUS_CONFIG}"
log "Threads: ${THREADS}"
log "DB root: ${GLNEXUS_DB_ROOT}"
df -h "${GLNEXUS_DB_ROOT}" >&2 || true

START_TIME="${SECONDS}"

if "${CONTAINER_RUNTIME}" \
    "${CONTAINER_ARGS[@]}" \
    "${GLNEXUS_IMAGE}" \
    "${GLNEXUS_BIN_IN_CONTAINER}" \
    --config "${GLNEXUS_CONFIG}" \
    --threads "${THREADS}" \
    --dir "${DB_DIR}" \
    --list "${GVCF_LIST}" \
    > "${BCF}" \
    2> "${GLNEXUS_LOG}"
then
    :
else
    status=$?
    die "GLnexus failed with exit status ${status}; see ${GLNEXUS_LOG}"
fi

# ---------- Validate BCF and create VCF ----------

[[ -s "${BCF}" ]] || die "GLnexus did not create a cohort BCF: ${BCF}"
"${BCFTOOLS_BIN}" view --header-only "${BCF}" >/dev/null ||
    die "GLnexus cohort BCF failed validation: ${BCF}"

"${BCFTOOLS_BIN}" index --threads "${THREADS}" --csi "${BCF}"
[[ -s "${BCF}.csi" ]] || die "Cohort BCF index was not created: ${BCF}.csi"

"${BCFTOOLS_BIN}" view \
    --threads "${THREADS}" \
    --output-type z \
    --output-file "${VCF}" \
    "${BCF}"

[[ -s "${VCF}" ]] || die "Cohort VCF was not created: ${VCF}"
"${BCFTOOLS_BIN}" view --header-only "${VCF}" >/dev/null ||
    die "Cohort VCF failed validation: ${VCF}"

"${BCFTOOLS_BIN}" index --threads "${THREADS}" --tbi "${VCF}" ||
    "${BCFTOOLS_BIN}" index --threads "${THREADS}" --csi "${VCF}"
has_variant_index "${VCF}" || die "Cohort VCF index was not created: ${VCF}"

"${BCFTOOLS_BIN}" stats --threads "${THREADS}" "${VCF}" > "${STATS}"
[[ -s "${STATS}" ]] || die "Cohort bcftools stats output is empty: ${STATS}"

# Confirm the output contains the expected number of samples.
OUTPUT_SAMPLE_COUNT="$("${BCFTOOLS_BIN}" query -l "${VCF}" | wc -l)"
[[ "${OUTPUT_SAMPLE_COUNT}" -eq "${N_GVCFS}" ]] ||
    die "Cohort output has ${OUTPUT_SAMPLE_COUNT} samples; expected ${N_GVCFS}"

# ---------- Publish completed cohort ----------

mv -- "${STAGING_DIR}" "${FINAL_DIR}"
STAGING_DIR=""

rm -rf -- "${DB_DIR}"
DB_DIR=""

log "GLnexus joint calling completed"
log "Cohort BCF: ${FINAL_BCF}"
log "Cohort VCF: ${FINAL_VCF}"
log "Samples   : ${N_GVCFS}"
log "Elapsed time: $((SECONDS - START_TIME)) seconds"
