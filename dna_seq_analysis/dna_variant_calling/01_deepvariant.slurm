#!/usr/bin/env bash
#SBATCH --job-name=deepvariant
#SBATCH --partition=GPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=48G
#SBATCH --time=2-00:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -Eeuo pipefail
export LC_ALL=C

# ==============================================================================
# 01_deepvariant.slurm
#
# Run DeepVariant on one analysis-ready BAM per Slurm array task.
# BAMs are discovered from INPUT_BAM_DIR and sorted deterministically by path.
# Default input pattern: *.markdup.bam
# ==============================================================================

INPUT_BAM_DIR="${INPUT_BAM_DIR:-}"
REFERENCE="${REFERENCE:-}"
OUTPUT_DIR="${OUTPUT_DIR:-}"

BAM_SUFFIX="${BAM_SUFFIX:-.markdup.bam}"
BAM_PATTERN="${BAM_PATTERN:-*${BAM_SUFFIX}}"

SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}"
BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"

MODEL_TYPE="${MODEL_TYPE:-WGS}"
NUM_SHARDS="${NUM_SHARDS:-${SLURM_CPUS_PER_TASK:-8}}"

# USE_GPU=1 requires a GPU allocation and a DeepVariant GPU image.
USE_GPU="${USE_GPU:-1}"
DEEPVARIANT_VERSION="${DEEPVARIANT_VERSION:-1.10.0}"

if [[ "${USE_GPU}" == "1" ]]; then
    DEFAULT_IMAGE="docker://google/deepvariant:${DEEPVARIANT_VERSION}-gpu"
else
    DEFAULT_IMAGE="docker://google/deepvariant:${DEEPVARIANT_VERSION}"
fi

# For production runs, prefer a pre-pulled local SIF over docker://.
DEEPVARIANT_IMAGE="${DEEPVARIANT_IMAGE:-${DEFAULT_IMAGE}}"
CONTAINER_RUNTIME="${CONTAINER_RUNTIME:-}"
EXTRA_BINDS="${EXTRA_BINDS:-}"

REGIONS="${REGIONS:-}"
CUSTOM_MODEL="${CUSTOM_MODEL:-}"
VCF_STATS_REPORT="${VCF_STATS_REPORT:-true}"
MAKE_EXAMPLES_EXTRA_ARGS="${MAKE_EXAMPLES_EXTRA_ARGS:-}"
CALL_VARIANTS_EXTRA_ARGS="${CALL_VARIANTS_EXTRA_ARGS:-}"
POSTPROCESS_VARIANTS_EXTRA_ARGS="${POSTPROCESS_VARIANTS_EXTRA_ARGS:-}"

log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2
}

die() {
    log "ERROR: $*"
    exit 1
}

STAGING_DIR=""
TASK_TMP=""

cleanup_partial_outputs() {
    if [[ -n "${STAGING_DIR}" && -d "${STAGING_DIR}" ]]; then
        rm -rf -- "${STAGING_DIR}"
    fi

    if [[ -n "${TASK_TMP}" &&
          -d "${TASK_TMP}" &&
          "${TASK_TMP}" == */deepvariant_* ]]
    then
        rm -rf -- "${TASK_TMP}"
    fi
}

trap cleanup_partial_outputs EXIT

has_variant_index() {
    local variant_file="$1"
    [[ -s "${variant_file}.tbi" || -s "${variant_file}.csi" ]]
}

# ---------- Validate configuration ----------

[[ -n "${INPUT_BAM_DIR}" ]] || die "INPUT_BAM_DIR is required"
[[ -d "${INPUT_BAM_DIR}" ]] || die "INPUT_BAM_DIR does not exist: ${INPUT_BAM_DIR}"
[[ -n "${REFERENCE}" ]] || die "REFERENCE is required"
[[ -s "${REFERENCE}" ]] || die "Reference FASTA does not exist or is empty: ${REFERENCE}"
[[ -s "${REFERENCE}.fai" ]] || die "Reference FASTA index is missing: ${REFERENCE}.fai"
[[ -n "${OUTPUT_DIR}" ]] || die "OUTPUT_DIR is required"
[[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || die "Submit this script as a Slurm array job"
[[ "${SLURM_ARRAY_TASK_ID}" =~ ^[0-9]+$ ]] ||
    die "Invalid SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
[[ "${NUM_SHARDS}" =~ ^[1-9][0-9]*$ ]] ||
    die "NUM_SHARDS must be a positive integer: ${NUM_SHARDS}"
[[ "${USE_GPU}" == "0" || "${USE_GPU}" == "1" ]] ||
    die "USE_GPU must be 0 or 1"
[[ "${VCF_STATS_REPORT}" == "true" || "${VCF_STATS_REPORT}" == "false" ]] ||
    die "VCF_STATS_REPORT must be true or false"

TOTAL_CPUS="${SLURM_CPUS_PER_TASK:-1}"
if (( NUM_SHARDS > TOTAL_CPUS )); then
    die "NUM_SHARDS=${NUM_SHARDS} exceeds allocated CPUs=${TOTAL_CPUS}"
fi

case "${MODEL_TYPE}" in
    WGS|WES|PACBIO|ONT_R104|HYBRID_PACBIO_ILLUMINA)
        ;;
    *)
        die "Unsupported MODEL_TYPE: ${MODEL_TYPE}"
        ;;
esac

command -v "${SAMTOOLS_BIN}" >/dev/null 2>&1 ||
    die "samtools executable not found: ${SAMTOOLS_BIN}"
command -v "${BCFTOOLS_BIN}" >/dev/null 2>&1 ||
    die "bcftools executable not found: ${BCFTOOLS_BIN}"

# Resolve important directories before passing paths into the container.
INPUT_BAM_DIR="$(cd "${INPUT_BAM_DIR}" && pwd -P)"
REFERENCE_DIR="$(cd "$(dirname "${REFERENCE}")" && pwd -P)"
REFERENCE="${REFERENCE_DIR}/$(basename "${REFERENCE}")"

mkdir -p "${OUTPUT_DIR}"
OUTPUT_DIR="$(cd "${OUTPUT_DIR}" && pwd -P)"

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

if [[ "${DEEPVARIANT_IMAGE}" != *"://"* ]]; then
    [[ -s "${DEEPVARIANT_IMAGE}" ]] ||
        die "DeepVariant image does not exist: ${DEEPVARIANT_IMAGE}"
fi

if [[ -n "${CUSTOM_MODEL}" ]]; then
    [[ -e "${CUSTOM_MODEL}" || -e "${CUSTOM_MODEL}.index" ]] ||
        die "CUSTOM_MODEL does not exist: ${CUSTOM_MODEL}"
fi

# ---------- Discover and select one BAM ----------

mapfile -d '' -t BAM_FILES < <(
    find "${INPUT_BAM_DIR}" \
        -maxdepth 1 \
        -type f \
        -name "${BAM_PATTERN}" \
        -print0 |
        sort -z
)

N_BAMS="${#BAM_FILES[@]}"
(( N_BAMS > 0 )) ||
    die "No BAM files matching '${BAM_PATTERN}' found in ${INPUT_BAM_DIR}"

IDX="${SLURM_ARRAY_TASK_ID}"
if (( IDX >= N_BAMS )); then
    die "Array index ${IDX} is out of range; found ${N_BAMS} BAM files"
fi

INPUT_BAM="${BAM_FILES[${IDX}]}"
BAM_BASENAME="$(basename "${INPUT_BAM}")"

[[ "${BAM_BASENAME}" == *"${BAM_SUFFIX}" ]] ||
    die "BAM filename does not end with ${BAM_SUFFIX}: ${BAM_BASENAME}"

SAMPLE="${BAM_BASENAME%"${BAM_SUFFIX}"}"
[[ -n "${SAMPLE}" ]] || die "Could not derive sample name from ${BAM_BASENAME}"
[[ "${SAMPLE}" =~ ^[A-Za-z0-9._-]+$ ]] ||
    die "Derived sample name contains unsupported characters: ${SAMPLE}"
[[ -s "${INPUT_BAM}" && -r "${INPUT_BAM}" ]] ||
    die "Input BAM is not readable: ${INPUT_BAM}"

"${SAMTOOLS_BIN}" quickcheck -v "${INPUT_BAM}" ||
    die "Input BAM failed quickcheck: ${INPUT_BAM}"

BAM_STEM="${INPUT_BAM%.bam}"
if [[ ! -s "${INPUT_BAM}.bai" &&
      ! -s "${BAM_STEM}.bai" &&
      ! -s "${INPUT_BAM}.csi" &&
      ! -s "${BAM_STEM}.csi" ]]
then
    die "Input BAM index is missing: ${INPUT_BAM}"
fi

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

# ---------- Output and temporary paths ----------

SAMPLE_ROOT="${OUTPUT_DIR}/samples"
LOG_DIR="${OUTPUT_DIR}/logs"
FINAL_DIR="${SAMPLE_ROOT}/${SAMPLE}"

mkdir -p "${SAMPLE_ROOT}" "${LOG_DIR}"

FINAL_VCF="${FINAL_DIR}/${SAMPLE}.deepvariant.vcf.gz"
FINAL_GVCF="${FINAL_DIR}/${SAMPLE}.deepvariant.g.vcf.gz"

# ---------- Resume protection ----------

if [[ -d "${FINAL_DIR}" ]]; then
    if [[ -s "${FINAL_VCF}" && -s "${FINAL_GVCF}" ]] &&
       has_variant_index "${FINAL_VCF}" &&
       has_variant_index "${FINAL_GVCF}"
    then
        "${BCFTOOLS_BIN}" view --header-only "${FINAL_VCF}" >/dev/null ||
            die "Existing VCF failed validation: ${FINAL_VCF}"
        "${BCFTOOLS_BIN}" view --header-only "${FINAL_GVCF}" >/dev/null ||
            die "Existing gVCF failed validation: ${FINAL_GVCF}"
        "${BCFTOOLS_BIN}" index --stats "${FINAL_VCF}" >/dev/null ||
            die "Existing VCF index failed validation: ${FINAL_VCF}"
        "${BCFTOOLS_BIN}" index --stats "${FINAL_GVCF}" >/dev/null ||
            die "Existing gVCF index failed validation: ${FINAL_GVCF}"

        log "Valid DeepVariant outputs already exist; skipping: ${SAMPLE}"
        exit 0
    fi

    die "Incomplete existing DeepVariant directory: ${FINAL_DIR}"
fi

JOB_TAG="${SLURM_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID}"
STAGING_DIR="${SAMPLE_ROOT}/.${SAMPLE}.${JOB_TAG}.partial"

TMP_ROOT="${SLURM_TMPDIR:-${OUTPUT_DIR}/tmp}"
TASK_TMP="${TMP_ROOT}/deepvariant_${JOB_TAG}"
INTERMEDIATE_DIR="${TASK_TMP}/intermediate"

[[ ! -e "${STAGING_DIR}" ]] || die "Staging directory already exists: ${STAGING_DIR}"
mkdir -p "${STAGING_DIR}" "${INTERMEDIATE_DIR}"

VCF="${STAGING_DIR}/${SAMPLE}.deepvariant.vcf.gz"
GVCF="${STAGING_DIR}/${SAMPLE}.deepvariant.g.vcf.gz"
DEEPVARIANT_LOG="${LOG_DIR}/${SAMPLE}.deepvariant.log"
CONFIG_LOG="${LOG_DIR}/${SAMPLE}.deepvariant.config.tsv"

# ---------- Construct container arguments ----------

TASK_TMP_ABS="$(cd "${TASK_TMP}" && pwd -P)"
STAGING_ABS="$(cd "${STAGING_DIR}" && pwd -P)"

CONTAINER_ARGS=(
    exec
    --cleanenv
    --bind "${REFERENCE_DIR}:${REFERENCE_DIR}"
    --bind "${INPUT_BAM_DIR}:${INPUT_BAM_DIR}"
    --bind "${STAGING_ABS}:${STAGING_ABS}"
    --bind "${TASK_TMP_ABS}:${TASK_TMP_ABS}"
)

if (( USE_GPU == 1 )); then
    CONTAINER_ARGS+=(--nv)
fi

if [[ -n "${EXTRA_BINDS}" ]]; then
    CONTAINER_ARGS+=(--bind "${EXTRA_BINDS}")
fi

DEEPVARIANT_ARGS=(
    /opt/deepvariant/bin/run_deepvariant
    --model_type="${MODEL_TYPE}"
    --ref="${REFERENCE}"
    --reads="${INPUT_BAM}"
    --output_vcf="${VCF}"
    --output_gvcf="${GVCF}"
    --sample_name="${SAMPLE}"
    --num_shards="${NUM_SHARDS}"
    --vcf_stats_report="${VCF_STATS_REPORT}"
    --intermediate_results_dir="${INTERMEDIATE_DIR}"
    --logging_dir="${STAGING_DIR}/stage_logs"
)

if [[ -n "${REGIONS}" ]]; then
    DEEPVARIANT_ARGS+=(--regions="${REGIONS}")

    if [[ -e "${REGIONS}" ]]; then
        REGIONS_DIR="$(cd "$(dirname "${REGIONS}")" && pwd -P)"
        REGIONS="${REGIONS_DIR}/$(basename "${REGIONS}")"
        CONTAINER_ARGS+=(--bind "${REGIONS_DIR}:${REGIONS_DIR}")
        DEEPVARIANT_ARGS[${#DEEPVARIANT_ARGS[@]}-1]="--regions=${REGIONS}"
    fi
fi

if [[ -n "${CUSTOM_MODEL}" ]]; then
    CUSTOM_MODEL_DIR="$(cd "$(dirname "${CUSTOM_MODEL}")" && pwd -P)"
    CUSTOM_MODEL="${CUSTOM_MODEL_DIR}/$(basename "${CUSTOM_MODEL}")"
    CONTAINER_ARGS+=(--bind "${CUSTOM_MODEL_DIR}:${CUSTOM_MODEL_DIR}")
    DEEPVARIANT_ARGS+=(--customized_model="${CUSTOM_MODEL}")
fi

if [[ -n "${MAKE_EXAMPLES_EXTRA_ARGS}" ]]; then
    DEEPVARIANT_ARGS+=(--make_examples_extra_args="${MAKE_EXAMPLES_EXTRA_ARGS}")
fi

if [[ -n "${CALL_VARIANTS_EXTRA_ARGS}" ]]; then
    DEEPVARIANT_ARGS+=(--call_variants_extra_args="${CALL_VARIANTS_EXTRA_ARGS}")
fi

if [[ -n "${POSTPROCESS_VARIANTS_EXTRA_ARGS}" ]]; then
    DEEPVARIANT_ARGS+=(--postprocess_variants_extra_args="${POSTPROCESS_VARIANTS_EXTRA_ARGS}")
fi

# ---------- Record configuration ----------

{
    printf 'sample\t%s\n' "${SAMPLE}"
    printf 'array_index\t%s\n' "${IDX}"
    printf 'bam_count\t%s\n' "${N_BAMS}"
    printf 'bam\t%s\n' "${INPUT_BAM}"
    printf 'reference\t%s\n' "${REFERENCE}"
    printf 'model_type\t%s\n' "${MODEL_TYPE}"
    printf 'deepvariant_version\t%s\n' "${DEEPVARIANT_VERSION}"
    printf 'image\t%s\n' "${DEEPVARIANT_IMAGE}"
    printf 'container_runtime\t%s\n' "${CONTAINER_RUNTIME}"
    printf 'use_gpu\t%s\n' "${USE_GPU}"
    printf 'num_shards\t%s\n' "${NUM_SHARDS}"
    printf 'regions\t%s\n' "${REGIONS:-whole_genome}"
    printf 'custom_model\t%s\n' "${CUSTOM_MODEL:-default}"
} > "${CONFIG_LOG}"

# ---------- Run DeepVariant ----------

log "Starting DeepVariant"
log "Array task: ${IDX}/${N_BAMS}"
log "Sample    : ${SAMPLE}"
log "BAM       : ${INPUT_BAM}"
log "Reference : ${REFERENCE}"
log "Model     : ${MODEL_TYPE}"
log "Image     : ${DEEPVARIANT_IMAGE}"
log "GPU mode  : ${USE_GPU}"
log "Shards    : ${NUM_SHARDS}"

START_TIME="${SECONDS}"

if "${CONTAINER_RUNTIME}" \
    "${CONTAINER_ARGS[@]}" \
    "${DEEPVARIANT_IMAGE}" \
    "${DEEPVARIANT_ARGS[@]}" \
    > "${DEEPVARIANT_LOG}" 2>&1
then
    :
else
    status=$?
    die "DeepVariant failed with exit status ${status}; see ${DEEPVARIANT_LOG}"
fi

# ---------- Validate outputs ----------

[[ -s "${VCF}" ]] || die "DeepVariant VCF is missing: ${VCF}"
[[ -s "${GVCF}" ]] || die "DeepVariant gVCF is missing: ${GVCF}"
has_variant_index "${VCF}" || die "DeepVariant VCF index is missing: ${VCF}"
has_variant_index "${GVCF}" || die "DeepVariant gVCF index is missing: ${GVCF}"

"${BCFTOOLS_BIN}" view --header-only "${VCF}" >/dev/null ||
    die "DeepVariant VCF failed validation: ${VCF}"
"${BCFTOOLS_BIN}" view --header-only "${GVCF}" >/dev/null ||
    die "DeepVariant gVCF failed validation: ${GVCF}"
"${BCFTOOLS_BIN}" index --stats "${VCF}" >/dev/null ||
    die "DeepVariant VCF index failed validation: ${VCF}"
"${BCFTOOLS_BIN}" index --stats "${GVCF}" >/dev/null ||
    die "DeepVariant gVCF index failed validation: ${GVCF}"

# Publish the complete sample directory in one rename operation.
mv -- "${STAGING_DIR}" "${FINAL_DIR}"
STAGING_DIR=""

rm -rf -- "${TASK_TMP}"
TASK_TMP=""

log "DeepVariant completed: ${SAMPLE}"
log "VCF : ${FINAL_VCF}"
log "gVCF: ${FINAL_GVCF}"
log "Elapsed time: $((SECONDS - START_TIME)) seconds"
