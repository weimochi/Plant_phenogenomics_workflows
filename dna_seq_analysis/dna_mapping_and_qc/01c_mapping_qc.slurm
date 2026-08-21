#!/usr/bin/env bash
#SBATCH --job-name=mapping-qc
#SBATCH --partition=CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -Eeuo pipefail

# ==============================================================================
# 01c_mapping_qc.slurm
#
# Run mapping QC on one markdup BAM per Slurm array task.
#
# MANIFEST format, without header:
# sample<TAB>markdup_bam
#
# Outputs per sample:
#   flagstat.txt
#   idxstats.txt
#   samtools.stats.txt
#   mosdepth outputs
#   mapping_qc_metrics.tsv
# ==============================================================================

MANIFEST="${MANIFEST:-}"
OUTPUT_DIR="${OUTPUT_DIR:-}"

SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}"
MOSDEPTH_BIN="${MOSDEPTH_BIN:-mosdepth}"

THREADS="${SLURM_CPUS_PER_TASK:-4}"

# mosdepth filters reads with any of these SAM flags.
# 1796 = unmapped + secondary + QC-fail + duplicate
MOSDEPTH_EXCLUDE_FLAG="${MOSDEPTH_EXCLUDE_FLAG:-1796}"
MOSDEPTH_MAPQ="${MOSDEPTH_MAPQ:-0}"
MOSDEPTH_PRECISION="${MOSDEPTH_PRECISION:-5}"

# 1 enables --fast-mode; 0 uses full CIGAR/mate-overlap processing.
MOSDEPTH_FAST_MODE="${MOSDEPTH_FAST_MODE:-1}"

# Coverage levels extracted from mosdepth's cumulative global distribution.
COVERAGE_LEVELS="${COVERAGE_LEVELS:-1,5,10,20,30}"

log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2
}

die() {
    log "ERROR: $*"
    exit 1
}

require_nonnegative_integer() {
    local name="$1"
    local value="$2"

    [[ "${value}" =~ ^[0-9]+$ ]] ||
        die "${name} must be a non-negative integer: ${value}"
}

STAGING_DIR=""

cleanup_staging() {
    if [[ -n "${STAGING_DIR}" && -d "${STAGING_DIR}" ]]; then
        rm -rf -- "${STAGING_DIR}"
    fi
}

trap cleanup_staging EXIT

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

require_nonnegative_integer "MOSDEPTH_EXCLUDE_FLAG" "${MOSDEPTH_EXCLUDE_FLAG}"
require_nonnegative_integer "MOSDEPTH_MAPQ" "${MOSDEPTH_MAPQ}"
require_nonnegative_integer "MOSDEPTH_PRECISION" "${MOSDEPTH_PRECISION}"

[[ "${MOSDEPTH_FAST_MODE}" == "0" || "${MOSDEPTH_FAST_MODE}" == "1" ]] ||
    die "MOSDEPTH_FAST_MODE must be 0 or 1"

[[ "${COVERAGE_LEVELS}" =~ ^[0-9]+(,[0-9]+)*$ ]] ||
    die "COVERAGE_LEVELS must be comma-separated integers: ${COVERAGE_LEVELS}"

command -v "${SAMTOOLS_BIN}" >/dev/null 2>&1 ||
    die "samtools executable not found: ${SAMTOOLS_BIN}"

command -v "${MOSDEPTH_BIN}" >/dev/null 2>&1 ||
    die "mosdepth executable not found: ${MOSDEPTH_BIN}"

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

# ---------- Check BAM sort order and index ----------

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

BAM_STEM="${INPUT_BAM%.bam}"

if [[ ! -s "${INPUT_BAM}.bai" &&
      ! -s "${BAM_STEM}.bai" &&
      ! -s "${INPUT_BAM}.csi" &&
      ! -s "${BAM_STEM}.csi" ]]
then
    die "Input BAM index is missing: ${INPUT_BAM}"
fi

# ---------- Output and staging directories ----------

QC_ROOT="${OUTPUT_DIR}/samples"
FINAL_DIR="${QC_ROOT}/${SAMPLE}"

JOB_TAG="${SLURM_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID}"
STAGING_DIR="${QC_ROOT}/.${SAMPLE}.${JOB_TAG}.partial"

mkdir -p "${QC_ROOT}"

# ---------- Resume protection ----------

REQUIRED_OUTPUTS=(
    "${FINAL_DIR}/${SAMPLE}.flagstat.txt"
    "${FINAL_DIR}/${SAMPLE}.idxstats.txt"
    "${FINAL_DIR}/${SAMPLE}.samtools.stats.txt"
    "${FINAL_DIR}/${SAMPLE}.mosdepth.mosdepth.summary.txt"
    "${FINAL_DIR}/${SAMPLE}.mosdepth.mosdepth.global.dist.txt"
    "${FINAL_DIR}/${SAMPLE}.mapping_qc_metrics.tsv"
)

if [[ -d "${FINAL_DIR}" ]]; then
    outputs_complete=1

    for output in "${REQUIRED_OUTPUTS[@]}"; do
        if [[ ! -s "${output}" ]]; then
            outputs_complete=0
            break
        fi
    done

    if (( outputs_complete == 1 )); then
        log "Complete mapping QC already exists; skipping: ${SAMPLE}"
        exit 0
    fi

    die "Incomplete existing QC directory found: ${FINAL_DIR}"
fi

[[ ! -e "${STAGING_DIR}" ]] ||
    die "Staging directory already exists: ${STAGING_DIR}"

mkdir -p "${STAGING_DIR}"

FLAGSTAT="${STAGING_DIR}/${SAMPLE}.flagstat.txt"
IDXSTATS="${STAGING_DIR}/${SAMPLE}.idxstats.txt"
SAMTOOLS_STATS="${STAGING_DIR}/${SAMPLE}.samtools.stats.txt"

MOSDEPTH_PREFIX="${STAGING_DIR}/${SAMPLE}.mosdepth"
MOSDEPTH_SUMMARY="${MOSDEPTH_PREFIX}.mosdepth.summary.txt"
MOSDEPTH_GLOBAL_DIST="${MOSDEPTH_PREFIX}.mosdepth.global.dist.txt"

METRICS="${STAGING_DIR}/${SAMPLE}.mapping_qc_metrics.tsv"
SAMTOOLS_LOG="${STAGING_DIR}/${SAMPLE}.samtools.log"
MOSDEPTH_LOG="${STAGING_DIR}/${SAMPLE}.mosdepth.log"
VERSION_LOG="${STAGING_DIR}/${SAMPLE}.tool_versions.txt"

# ---------- Record configuration ----------

{
    printf 'sample\t%s\n' "${SAMPLE}"
    printf 'input_bam\t%s\n' "${INPUT_BAM}"
    printf 'mosdepth_exclude_flag\t%s\n' "${MOSDEPTH_EXCLUDE_FLAG}"
    printf 'mosdepth_mapq\t%s\n' "${MOSDEPTH_MAPQ}"
    printf 'mosdepth_fast_mode\t%s\n' "${MOSDEPTH_FAST_MODE}"
    printf 'coverage_levels\t%s\n' "${COVERAGE_LEVELS}"

    printf '\n[samtools]\n'
    "${SAMTOOLS_BIN}" --version 2>&1 || true

    printf '\n[mosdepth]\n'
    "${MOSDEPTH_BIN}" --version 2>&1 || true
} > "${VERSION_LOG}"

# ---------- Run samtools QC ----------

log "Starting mapping QC: ${SAMPLE}"
log "Input BAM: ${INPUT_BAM}"

START_TIME="${SECONDS}"

{
    log "Running samtools flagstat"
    "${SAMTOOLS_BIN}" flagstat -@ "${THREADS}" "${INPUT_BAM}" > "${FLAGSTAT}"

    log "Running samtools idxstats"
    "${SAMTOOLS_BIN}" idxstats "${INPUT_BAM}" > "${IDXSTATS}"

    log "Running samtools stats"
    "${SAMTOOLS_BIN}" stats -@ "${THREADS}" "${INPUT_BAM}" > "${SAMTOOLS_STATS}"
} 2> "${SAMTOOLS_LOG}"

# ---------- Run mosdepth ----------

MOSDEPTH_ARGS=(
    --threads "${THREADS}"
    --no-per-base
    --flag "${MOSDEPTH_EXCLUDE_FLAG}"
    --mapq "${MOSDEPTH_MAPQ}"
)

if (( MOSDEPTH_FAST_MODE == 1 )); then
    MOSDEPTH_ARGS+=(--fast-mode)
fi

log "Running mosdepth"

if MOSDEPTH_PRECISION="${MOSDEPTH_PRECISION}" \
    "${MOSDEPTH_BIN}" \
    "${MOSDEPTH_ARGS[@]}" \
    "${MOSDEPTH_PREFIX}" \
    "${INPUT_BAM}" \
    > "${MOSDEPTH_LOG}" 2>&1
then
    :
else
    status=$?
    die "mosdepth failed with exit status ${status}; see ${MOSDEPTH_LOG}"
fi

# ---------- Validate QC outputs ----------

[[ -s "${FLAGSTAT}" ]] ||
    die "flagstat output is missing or empty"

[[ -s "${IDXSTATS}" ]] ||
    die "idxstats output is missing or empty"

[[ -s "${SAMTOOLS_STATS}" ]] ||
    die "samtools stats output is missing or empty"

[[ -s "${MOSDEPTH_SUMMARY}" ]] ||
    die "mosdepth summary is missing: ${MOSDEPTH_SUMMARY}"

[[ -s "${MOSDEPTH_GLOBAL_DIST}" ]] ||
    die "mosdepth global distribution is missing: ${MOSDEPTH_GLOBAL_DIST}"

# ---------- Build per-sample metrics ----------

MAPPED_PERCENT="$(
    awk -F '[()%]' '
        / mapped \(/ {
            value = $2
            sub(/%.*/, "", value)
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", value)
            print value
            exit
        }
    ' "${FLAGSTAT}"
)"

MEAN_DEPTH="$(
    awk '
        $1 == "total" {
            print $4
            exit
        }
    ' "${MOSDEPTH_SUMMARY}"
)"

[[ -n "${MAPPED_PERCENT}" ]] ||
    die "Could not extract mapped percentage from ${FLAGSTAT}"

[[ -n "${MEAN_DEPTH}" ]] ||
    die "Could not extract total mean depth from ${MOSDEPTH_SUMMARY}"

{
    printf 'sample\tmapped_pct\tmean_depth'

    IFS=',' read -r -a LEVEL_ARRAY <<< "${COVERAGE_LEVELS}"

    for level in "${LEVEL_ARRAY[@]}"; do
        printf '\tcov_%sx_pct' "${level}"
    done

    printf '\n'
    printf '%s\t%s\t%s' "${SAMPLE}" "${MAPPED_PERCENT}" "${MEAN_DEPTH}"

    for level in "${LEVEL_ARRAY[@]}"; do
        proportion="$(
            awk -v target="${level}" '
                $1 == "total" && $2 == target {
                    print $3
                    exit
                }
            ' "${MOSDEPTH_GLOBAL_DIST}"
        "

        [[ -n "${proportion}" ]] ||
            die "Coverage level ${level} was not found in ${MOSDEPTH_GLOBAL_DIST}"

        percentage="$(
            awk -v value="${proportion}" 'BEGIN { printf "%.5f", value * 100 }'
        )"

        printf '\t%s' "${percentage}"
    done

    printf '\n'
} > "${METRICS}"

[[ -s "${METRICS}" ]] ||
    die "Per-sample metrics file was not created"

# ---------- Atomically publish sample QC directory ----------

mv -- "${STAGING_DIR}" "${FINAL_DIR}"
STAGING_DIR=""

log "Mapping QC completed: ${SAMPLE}"
log "Results: ${FINAL_DIR}"
log "Elapsed time: $((SECONDS - START_TIME)) seconds"