#!/usr/bin/env bash
#SBATCH --job-name=filtered-variant-qc
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
# 05_filtered_variant_qc.slurm
#
# Compare a joint-called cohort VCF before and after variant filtering.
# This script does not modify either input VCF/BCF.
# ==============================================================================

RAW_VCF="${RAW_VCF:-}"
FILTERED_VCF="${FILTERED_VCF:-}"
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

validate_variant_file() {
    local label="$1"
    local variant_file="$2"

    [[ -s "${variant_file}" ]] || die "${label} does not exist or is empty: ${variant_file}"
    has_variant_index "${variant_file}" || die "${label} index is missing: ${variant_file}"
    "${BCFTOOLS_BIN}" view --header-only "${variant_file}" >/dev/null ||
        die "${label} failed header validation: ${variant_file}"
    "${BCFTOOLS_BIN}" index --stats "${variant_file}" >/dev/null ||
        die "${label} index failed validation: ${variant_file}"
}

extract_site_summary() {
    local stats_file="$1"
    local output_file="$2"

    {
        printf 'metric\tvalue\n'

        awk -F '\t' '
            $1 == "SN" {
                metric = $3
                sub(/:$/, "", metric)
                gsub(/[[:space:]]+/, "_", metric)
                print metric "\t" $4
            }
        ' "${stats_file}"
    } > "${output_file}"

    [[ "$(awk 'END { print NR + 0 }' "${output_file}")" -gt 1 ]] ||
        die "Could not extract site summary from ${stats_file}"
}

extract_sample_qc() {
    local stats_file="$1"
    local output_file="$2"
    local number_of_records

    number_of_records="$(
        awk -F '\t' '$1 == "SN" && $3 == "number of records:" { print $4; exit }' "${stats_file}"
    )"

    [[ "${number_of_records}" =~ ^[0-9]+$ ]] ||
        die "Could not extract number of records from ${stats_file}"

    {
        printf 'sample\tn_ref_hom_snp\tn_alt_hom_snp\tn_het_snp'
        printf '\tsnp_het_fraction\tn_transitions\tn_transversions\tts_tv'
        printf '\tn_indels\tmean_depth\tn_singletons\tn_missing\tmissing_pct\n'

        awk -F '\t' -v records="${number_of_records}" '
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
        ' "${stats_file}"
    } > "${output_file}"
}

# ---------- Validate configuration ----------

[[ -n "${RAW_VCF}" ]] || die "RAW_VCF is required"
[[ -n "${FILTERED_VCF}" ]] || die "FILTERED_VCF is required"
[[ -n "${OUTPUT_DIR}" ]] || die "OUTPUT_DIR is required"
[[ "${COHORT_NAME}" =~ ^[A-Za-z0-9._-]+$ ]] ||
    die "COHORT_NAME contains unsupported characters: ${COHORT_NAME}"
[[ "${THREADS}" =~ ^[1-9][0-9]*$ ]] ||
    die "THREADS must be a positive integer: ${THREADS}"

command -v "${BCFTOOLS_BIN}" >/dev/null 2>&1 ||
    die "bcftools executable not found: ${BCFTOOLS_BIN}"

validate_variant_file "RAW_VCF" "${RAW_VCF}"
validate_variant_file "FILTERED_VCF" "${FILTERED_VCF}"

RAW_DIR="$(cd "$(dirname "${RAW_VCF}")" && pwd -P)"
FILTERED_DIR="$(cd "$(dirname "${FILTERED_VCF}")" && pwd -P)"
RAW_VCF="${RAW_DIR}/$(basename "${RAW_VCF}")"
FILTERED_VCF="${FILTERED_DIR}/$(basename "${FILTERED_VCF}")"

[[ "${RAW_VCF}" != "${FILTERED_VCF}" ]] ||
    die "RAW_VCF and FILTERED_VCF resolve to the same file"

mkdir -p "${OUTPUT_DIR}"
OUTPUT_DIR="$(cd "${OUTPUT_DIR}" && pwd -P)"

# ---------- Output paths and resume protection ----------

FINAL_DIR="${OUTPUT_DIR}/${COHORT_NAME}"

REQUIRED_NAMES=(
    "${COHORT_NAME}.raw.bcftools.stats.txt"
    "${COHORT_NAME}.filtered.bcftools.stats.txt"
    "${COHORT_NAME}.raw.site_summary.tsv"
    "${COHORT_NAME}.filtered.site_summary.tsv"
    "${COHORT_NAME}.raw.sample_qc.tsv"
    "${COHORT_NAME}.filtered.sample_qc.tsv"
    "${COHORT_NAME}.pre_post_site_comparison.tsv"
    "${COHORT_NAME}.pre_post_sample_comparison.tsv"
    "${COHORT_NAME}.filtered.filter_summary.tsv"
    "${COHORT_NAME}.filtered.contig_summary.tsv"
    "${COHORT_NAME}.samples.txt"
)

if [[ -d "${FINAL_DIR}" ]]; then
    complete=1

    for name in "${REQUIRED_NAMES[@]}"; do
        if [[ ! -s "${FINAL_DIR}/${name}" ]]; then
            complete=0
            break
        fi
    done

    if (( complete == 1 )); then
        log "Complete filtered variant QC already exists; skipping: ${COHORT_NAME}"
        exit 0
    fi

    die "Incomplete existing filtered variant QC directory: ${FINAL_DIR}"
fi

JOB_TAG="${SLURM_JOB_ID:-manual}"
STAGING_DIR="${OUTPUT_DIR}/.${COHORT_NAME}.${JOB_TAG}.filtered_variant_qc.partial"
[[ ! -e "${STAGING_DIR}" ]] || die "Staging directory already exists: ${STAGING_DIR}"
mkdir -p "${STAGING_DIR}"

RAW_STATS="${STAGING_DIR}/${COHORT_NAME}.raw.bcftools.stats.txt"
FILTERED_STATS="${STAGING_DIR}/${COHORT_NAME}.filtered.bcftools.stats.txt"
RAW_SITE_SUMMARY="${STAGING_DIR}/${COHORT_NAME}.raw.site_summary.tsv"
FILTERED_SITE_SUMMARY="${STAGING_DIR}/${COHORT_NAME}.filtered.site_summary.tsv"
RAW_SAMPLE_QC="${STAGING_DIR}/${COHORT_NAME}.raw.sample_qc.tsv"
FILTERED_SAMPLE_QC="${STAGING_DIR}/${COHORT_NAME}.filtered.sample_qc.tsv"
SITE_COMPARISON="${STAGING_DIR}/${COHORT_NAME}.pre_post_site_comparison.tsv"
SAMPLE_COMPARISON="${STAGING_DIR}/${COHORT_NAME}.pre_post_sample_comparison.tsv"
FILTER_SUMMARY="${STAGING_DIR}/${COHORT_NAME}.filtered.filter_summary.tsv"
CONTIG_SUMMARY="${STAGING_DIR}/${COHORT_NAME}.filtered.contig_summary.tsv"
RAW_SAMPLE_LIST="${STAGING_DIR}/.raw.samples.tmp"
FILTERED_SAMPLE_LIST="${STAGING_DIR}/${COHORT_NAME}.samples.txt"
CONFIG_LOG="${STAGING_DIR}/${COHORT_NAME}.filtered_variant_qc.config.tsv"
RAW_LOG="${STAGING_DIR}/${COHORT_NAME}.raw.bcftools.log"
FILTERED_LOG="${STAGING_DIR}/${COHORT_NAME}.filtered.bcftools.log"

# ---------- Record configuration ----------

{
    printf 'cohort_name\t%s\n' "${COHORT_NAME}"
    printf 'raw_vcf\t%s\n' "${RAW_VCF}"
    printf 'filtered_vcf\t%s\n' "${FILTERED_VCF}"
    printf 'threads\t%s\n' "${THREADS}"
    printf '\n[bcftools]\n'
    "${BCFTOOLS_BIN}" --version 2>&1 || true
} > "${CONFIG_LOG}"

# ---------- Run identical QC before and after filtering ----------

log "Starting pre/post-filter variant QC"
log "Raw VCF     : ${RAW_VCF}"
log "Filtered VCF: ${FILTERED_VCF}"

START_TIME="${SECONDS}"

if "${BCFTOOLS_BIN}" stats \
    --threads "${THREADS}" \
    --samples - \
    "${RAW_VCF}" \
    > "${RAW_STATS}" \
    2> "${RAW_LOG}"
then
    :
else
    status=$?
    die "bcftools stats failed for RAW_VCF with exit status ${status}"
fi

if "${BCFTOOLS_BIN}" stats \
    --threads "${THREADS}" \
    --samples - \
    "${FILTERED_VCF}" \
    > "${FILTERED_STATS}" \
    2> "${FILTERED_LOG}"
then
    :
else
    status=$?
    die "bcftools stats failed for FILTERED_VCF with exit status ${status}"
fi

[[ -s "${RAW_STATS}" ]] || die "Raw bcftools stats output is empty"
[[ -s "${FILTERED_STATS}" ]] || die "Filtered bcftools stats output is empty"

# ---------- Confirm sample identity and order ----------

"${BCFTOOLS_BIN}" query -l "${RAW_VCF}" > "${RAW_SAMPLE_LIST}"
"${BCFTOOLS_BIN}" query -l "${FILTERED_VCF}" > "${FILTERED_SAMPLE_LIST}"

[[ -s "${RAW_SAMPLE_LIST}" ]] || die "RAW_VCF contains no samples"
[[ -s "${FILTERED_SAMPLE_LIST}" ]] || die "FILTERED_VCF contains no samples"

cmp -s "${RAW_SAMPLE_LIST}" "${FILTERED_SAMPLE_LIST}" ||
    die "Sample names or sample order changed during filtering"

SAMPLE_COUNT="$(wc -l < "${FILTERED_SAMPLE_LIST}")"
UNIQUE_SAMPLE_COUNT="$(sort -u "${FILTERED_SAMPLE_LIST}" | wc -l)"
[[ "${UNIQUE_SAMPLE_COUNT}" -eq "${SAMPLE_COUNT}" ]] ||
    die "Duplicate sample names found in cohort"

# ---------- Extract normalized summaries ----------

extract_site_summary "${RAW_STATS}" "${RAW_SITE_SUMMARY}"
extract_site_summary "${FILTERED_STATS}" "${FILTERED_SITE_SUMMARY}"
extract_sample_qc "${RAW_STATS}" "${RAW_SAMPLE_QC}"
extract_sample_qc "${FILTERED_STATS}" "${FILTERED_SAMPLE_QC}"

RAW_SAMPLE_ROWS=$(( $(wc -l < "${RAW_SAMPLE_QC}") - 1 ))
FILTERED_SAMPLE_ROWS=$(( $(wc -l < "${FILTERED_SAMPLE_QC}") - 1 ))
[[ "${RAW_SAMPLE_ROWS}" -eq "${SAMPLE_COUNT}" ]] ||
    die "Raw per-sample QC row count does not match sample count"
[[ "${FILTERED_SAMPLE_ROWS}" -eq "${SAMPLE_COUNT}" ]] ||
    die "Filtered per-sample QC row count does not match sample count"

# ---------- Site-level pre/post comparison ----------

awk -F '\t' '
    BEGIN { OFS = "\t" }

    NR == FNR {
        if (FNR > 1) {
            before[$1] = $2
            order[++n] = $1
        }
        next
    }

    FNR > 1 { after[$1] = $2 }

    END {
        print "metric", "before", "after", "retained_pct"

        for (i = 1; i <= n; i++) {
            metric = order[i]
            b = before[metric]
            a = (metric in after ? after[metric] : "NA")

            if (b ~ /^[0-9]+([.][0-9]+)?$/ &&
                a ~ /^[0-9]+([.][0-9]+)?$/ && b > 0) {
                retained = 100 * a / b
                printf "%s\t%s\t%s\t%.6f\n", metric, b, a, retained
            } else {
                print metric, b, a, "NA"
            }
        }
    }
' "${RAW_SITE_SUMMARY}" "${FILTERED_SITE_SUMMARY}" > "${SITE_COMPARISON}"

# ---------- Per-sample pre/post comparison ----------

awk -F '\t' '
    BEGIN { OFS = "\t" }

    NR == FNR {
        if (FNR > 1) {
            sample = $1
            raw_called[sample] = $2 + $3 + $4
            raw_het[sample] = $5
            raw_tstv[sample] = $8
            raw_depth[sample] = $10
            raw_missing[sample] = $13
        }
        next
    }

    FNR == 1 {
        print "sample", "raw_called_snp_gt", "filtered_called_snp_gt", \
              "called_snp_gt_retained_pct", "raw_snp_het_fraction", \
              "filtered_snp_het_fraction", "delta_snp_het_fraction", \
              "raw_ts_tv", "filtered_ts_tv", "raw_mean_depth", \
              "filtered_mean_depth", "raw_missing_pct", \
              "filtered_missing_pct", "delta_missing_pct"
        next
    }

    {
        sample = $1

        if (!(sample in raw_called)) {
            print "ERROR: sample missing from raw QC: " sample > "/dev/stderr"
            exit 2
        }

        filtered_called = $2 + $3 + $4
        retained = raw_called[sample] > 0 ? 100 * filtered_called / raw_called[sample] : 0

        printf "%s\t%d\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%s\t%s\t%.6f\t%.6f\t%.6f\n", \
            sample, raw_called[sample], filtered_called, retained, \
            raw_het[sample], $5, $5 - raw_het[sample], \
            raw_tstv[sample], $8, raw_depth[sample], $10, \
            raw_missing[sample], $13, $13 - raw_missing[sample]
    }
' "${RAW_SAMPLE_QC}" "${FILTERED_SAMPLE_QC}" > "${SAMPLE_COMPARISON}"

# ---------- Filter and contig summaries for the filtered VCF ----------

{
    printf 'filter\tn_records\n'

    "${BCFTOOLS_BIN}" query -f '%FILTER\n' "${FILTERED_VCF}" |
        sort |
        uniq -c |
        awk 'BEGIN { OFS = "\t" } { count = $1; $1 = ""; sub(/^ +/, ""); print $0, count }'
} > "${FILTER_SUMMARY}"

{
    printf 'contig\tlength\tn_records\n'
    "${BCFTOOLS_BIN}" index --stats "${FILTERED_VCF}"
} > "${CONTIG_SUMMARY}"

# ---------- Final validation and publish ----------

for name in "${REQUIRED_NAMES[@]}"; do
    [[ -s "${STAGING_DIR}/${name}" ]] ||
        die "Expected QC output is missing or empty: ${STAGING_DIR}/${name}"
done

[[ -s "${CONFIG_LOG}" ]] || die "Configuration log is missing"

rm -f -- "${RAW_SAMPLE_LIST}"

mv -- "${STAGING_DIR}" "${FINAL_DIR}"
STAGING_DIR=""

RAW_RECORDS="$(
    awk -F '\t' '$1 == "number_of_records" { print $2; exit }' "${FINAL_DIR}/${COHORT_NAME}.raw.site_summary.tsv"
)"
FILTERED_RECORDS="$(
    awk -F '\t' '$1 == "number_of_records" { print $2; exit }' "${FINAL_DIR}/${COHORT_NAME}.filtered.site_summary.tsv"
)"

log "Filtered variant QC completed"
log "Samples         : ${SAMPLE_COUNT}"
log "Raw records     : ${RAW_RECORDS:-NA}"
log "Filtered records: ${FILTERED_RECORDS:-NA}"
log "Output          : ${FINAL_DIR}"
log "Elapsed time: $((SECONDS - START_TIME)) seconds"
