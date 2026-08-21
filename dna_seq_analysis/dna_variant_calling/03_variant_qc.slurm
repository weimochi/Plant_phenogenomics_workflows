#!/usr/bin/env bash
#SBATCH --job-name=variant-qc
#SBATCH --partition=CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -Eeuo pipefail
export LC_ALL=C

# ==============================================================================
# 03_variant_qc.slurm
#
# Generate baseline QC for one joint-called cohort VCF/BCF before filtering.
# This script does not modify or filter variants.
# ==============================================================================

INPUT_VCF="${INPUT_VCF:-}"
OUTPUT_DIR="${OUTPUT_DIR:-}"
COHORT_NAME="${COHORT_NAME:-cohort}"

BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2
}

die() {
    log "ERROR: $*"
    exit 1
}

STAGING_DIR=""

cleanup_staging() {
    if [[ -n "${STAGING_DIR}" && -d "${STAGING_DIR}" ]]; then
        rm -rf -- "${STAGING_DIR}"
    fi
}

trap cleanup_staging EXIT

has_variant_index() {
    local variant_file="$1"
    [[ -s "${variant_file}.tbi" || -s "${variant_file}.csi" ]]
}

# ---------- Validate configuration ----------

[[ -n "${INPUT_VCF}" ]] || die "INPUT_VCF is required"
[[ -s "${INPUT_VCF}" ]] || die "INPUT_VCF does not exist or is empty: ${INPUT_VCF}"
[[ -n "${OUTPUT_DIR}" ]] || die "OUTPUT_DIR is required"
[[ "${COHORT_NAME}" =~ ^[A-Za-z0-9._-]+$ ]] ||
    die "COHORT_NAME contains unsupported characters: ${COHORT_NAME}"
[[ "${THREADS}" =~ ^[1-9][0-9]*$ ]] ||
    die "THREADS must be a positive integer: ${THREADS}"

command -v "${BCFTOOLS_BIN}" >/dev/null 2>&1 ||
    die "bcftools executable not found: ${BCFTOOLS_BIN}"

has_variant_index "${INPUT_VCF}" ||
    die "Input VCF/BCF index is missing: ${INPUT_VCF}"

"${BCFTOOLS_BIN}" view --header-only "${INPUT_VCF}" >/dev/null ||
    die "Input VCF/BCF failed header validation: ${INPUT_VCF}"

"${BCFTOOLS_BIN}" index --stats "${INPUT_VCF}" >/dev/null ||
    die "Input VCF/BCF index failed validation: ${INPUT_VCF}"

INPUT_DIR="$(cd "$(dirname "${INPUT_VCF}")" && pwd -P)"
INPUT_VCF="${INPUT_DIR}/$(basename "${INPUT_VCF}")"

mkdir -p "${OUTPUT_DIR}"
OUTPUT_DIR="$(cd "${OUTPUT_DIR}" && pwd -P)"

# ---------- Output paths ----------

FINAL_DIR="${OUTPUT_DIR}/${COHORT_NAME}"

FINAL_STATS="${FINAL_DIR}/${COHORT_NAME}.bcftools.stats.txt"
FINAL_SITE_SUMMARY="${FINAL_DIR}/${COHORT_NAME}.site_summary.tsv"
FINAL_SAMPLE_QC="${FINAL_DIR}/${COHORT_NAME}.sample_qc.tsv"
FINAL_FILTER_SUMMARY="${FINAL_DIR}/${COHORT_NAME}.filter_summary.tsv"
FINAL_CONTIG_SUMMARY="${FINAL_DIR}/${COHORT_NAME}.contig_summary.tsv"
FINAL_SAMPLE_LIST="${FINAL_DIR}/${COHORT_NAME}.samples.txt"

REQUIRED_OUTPUTS=(
    "${FINAL_STATS}"
    "${FINAL_SITE_SUMMARY}"
    "${FINAL_SAMPLE_QC}"
    "${FINAL_FILTER_SUMMARY}"
    "${FINAL_CONTIG_SUMMARY}"
    "${FINAL_SAMPLE_LIST}"
)

if [[ -d "${FINAL_DIR}" ]]; then
    complete=1

    for output in "${REQUIRED_OUTPUTS[@]}"; do
        if [[ ! -s "${output}" ]]; then
            complete=0
            break
        fi
    done

    if (( complete == 1 )); then
        log "Complete variant QC already exists; skipping: ${COHORT_NAME}"
        exit 0
    fi

    die "Incomplete existing variant QC directory: ${FINAL_DIR}"
fi

JOB_TAG="${SLURM_JOB_ID:-manual}"
STAGING_DIR="${OUTPUT_DIR}/.${COHORT_NAME}.${JOB_TAG}.variant_qc.partial"
[[ ! -e "${STAGING_DIR}" ]] || die "Staging directory already exists: ${STAGING_DIR}"
mkdir -p "${STAGING_DIR}"

STATS="${STAGING_DIR}/${COHORT_NAME}.bcftools.stats.txt"
SITE_SUMMARY="${STAGING_DIR}/${COHORT_NAME}.site_summary.tsv"
SAMPLE_QC="${STAGING_DIR}/${COHORT_NAME}.sample_qc.tsv"
FILTER_SUMMARY="${STAGING_DIR}/${COHORT_NAME}.filter_summary.tsv"
CONTIG_SUMMARY="${STAGING_DIR}/${COHORT_NAME}.contig_summary.tsv"
SAMPLE_LIST="${STAGING_DIR}/${COHORT_NAME}.samples.txt"
CONFIG_LOG="${STAGING_DIR}/${COHORT_NAME}.variant_qc.config.tsv"
BCFTOOLS_LOG="${STAGING_DIR}/${COHORT_NAME}.bcftools.log"

# ---------- Record configuration ----------

{
    printf 'cohort_name\t%s\n' "${COHORT_NAME}"
    printf 'input_vcf\t%s\n' "${INPUT_VCF}"
    printf 'threads\t%s\n' "${THREADS}"
    printf '\n[bcftools]\n'
    "${BCFTOOLS_BIN}" --version 2>&1 || true
} > "${CONFIG_LOG}"

# ---------- Run baseline bcftools QC ----------

log "Starting cohort variant QC"
log "Cohort: ${COHORT_NAME}"
log "Input : ${INPUT_VCF}"

START_TIME="${SECONDS}"

if "${BCFTOOLS_BIN}" stats \
    --threads "${THREADS}" \
    --samples - \
    "${INPUT_VCF}" \
    > "${STATS}" \
    2> "${BCFTOOLS_LOG}"
then
    :
else
    status=$?
    die "bcftools stats failed with exit status ${status}; see ${BCFTOOLS_LOG}"
fi

[[ -s "${STATS}" ]] || die "bcftools stats output is empty: ${STATS}"

# ---------- Sample list ----------

"${BCFTOOLS_BIN}" query -l "${INPUT_VCF}" > "${SAMPLE_LIST}"
[[ -s "${SAMPLE_LIST}" ]] || die "No samples found in cohort VCF/BCF"

SAMPLE_COUNT="$(wc -l < "${SAMPLE_LIST}")"
[[ "${SAMPLE_COUNT}" -gt 0 ]] || die "Cohort sample count is zero"

UNIQUE_SAMPLE_COUNT="$(sort -u "${SAMPLE_LIST}" | wc -l)"
[[ "${UNIQUE_SAMPLE_COUNT}" -eq "${SAMPLE_COUNT}" ]] ||
    die "Duplicate sample names found in cohort VCF/BCF"

# ---------- Overall site summary ----------

{
    printf 'metric\tvalue\n'

    awk -F '\t' '
        $1 == "SN" {
            metric = $3
            sub(/:$/, "", metric)
            gsub(/[[:space:]]+/, "_", metric)
            print metric "\t" $4
        }
    ' "${STATS}"
} > "${SITE_SUMMARY}"

[[ "$(awk 'END { print NR + 0 }' "${SITE_SUMMARY}")" -gt 1 ]] ||
    die "Could not extract site summary from bcftools stats"

NUMBER_OF_RECORDS="$(
    awk -F '\t' '$1 == "SN" && $3 == "number of records:" { print $4; exit }' "${STATS}"
)"

[[ "${NUMBER_OF_RECORDS}" =~ ^[0-9]+$ ]] ||
    die "Could not extract number of records from bcftools stats"

# ---------- Per-sample QC from PSC records ----------

{
    printf 'sample\tn_ref_hom_snp\tn_alt_hom_snp\tn_het_snp'
    printf '\tsnp_het_fraction\tn_transitions\tn_transversions\tts_tv'
    printf '\tn_indels\tmean_depth\tn_singletons\tn_missing\tmissing_pct\n'

    awk -F '\t' -v records="${NUMBER_OF_RECORDS}" '
        BEGIN { OFS = "\t" }

        $1 == "PSC" {
            sample = $3
            ref_hom = $4 + 0
            alt_hom = $5 + 0
            het = $6 + 0
            ts = $7 + 0
            tv = $8 + 0
            indels = $9 + 0
            mean_depth = $10
            singletons = $11 + 0
            missing = ($14 == "" ? 0 : $14 + 0)

            snp_called = ref_hom + alt_hom + het
            het_fraction = snp_called > 0 ? het / snp_called : 0
            ts_tv = tv > 0 ? ts / tv : 0
            missing_pct = records > 0 ? 100 * missing / records : 0

            printf "%s\t%d\t%d\t%d\t%.6f\t%d\t%d\t%.6f\t%d\t%s\t%d\t%d\t%.6f\n", \
                sample, ref_hom, alt_hom, het, het_fraction, ts, tv, ts_tv, \
                indels, mean_depth, singletons, missing, missing_pct
        }
    ' "${STATS}"
} > "${SAMPLE_QC}"

SAMPLE_QC_ROWS=$(( $(wc -l < "${SAMPLE_QC}") - 1 ))
[[ "${SAMPLE_QC_ROWS}" -eq "${SAMPLE_COUNT}" ]] ||
    die "Per-sample QC has ${SAMPLE_QC_ROWS} rows; expected ${SAMPLE_COUNT}"

# ---------- FILTER distribution ----------

{
    printf 'filter\tn_records\n'

    "${BCFTOOLS_BIN}" query -f '%FILTER\n' "${INPUT_VCF}" |
        sort |
        uniq -c |
        awk 'BEGIN { OFS = "\t" } { count = $1; $1 = ""; sub(/^ +/, ""); print $0, count }'
} > "${FILTER_SUMMARY}"

# ---------- Per-contig record counts ----------

{
    printf 'contig\tlength\tn_records\n'
    "${BCFTOOLS_BIN}" index --stats "${INPUT_VCF}"
} > "${CONTIG_SUMMARY}"

# ---------- Final validation ----------

for output in \
    "${SITE_SUMMARY}" \
    "${SAMPLE_QC}" \
    "${FILTER_SUMMARY}" \
    "${CONTIG_SUMMARY}" \
    "${SAMPLE_LIST}" \
    "${CONFIG_LOG}"
do
    [[ -s "${output}" ]] || die "Expected QC output is missing or empty: ${output}"
done

# Publish the complete QC directory in one rename operation.
mv -- "${STAGING_DIR}" "${FINAL_DIR}"
STAGING_DIR=""

log "Variant QC completed"
log "Samples: ${SAMPLE_COUNT}"
log "Records: ${NUMBER_OF_RECORDS}"
log "Output : ${FINAL_DIR}"
log "Elapsed time: $((SECONDS - START_TIME)) seconds"
