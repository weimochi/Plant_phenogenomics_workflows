#!/usr/bin/env bash
#SBATCH --job-name=fastp-trim
#SBATCH --partition=CPU
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

# ======================
# User-configurable parameters
# ======================
INPUT_DIR="${INPUT_DIR:-}"          # Auto-pair mode: folder with fastq.gz
MANIFEST="${MANIFEST:-}"            # Manifest mode: TSV with sample \t R1 \t R2
OUTPUT_DIR="${OUTPUT_DIR:-}"        # required
TRIM_FRONT="${TRIM_FRONT:-15}"      # cut front bp for both reads
PATTERN_R1="${PATTERN_R1:-*_R1_*.fastq.gz}"   
R2_TOKEN_FROM_R1="${R2_TOKEN_FROM_R1:-_R1_}"  
R2_TOKEN="${R2_TOKEN:-_R2_}"                 

# Adapters (Illumina standard by default, configurable via --export)
AD1="${AD1:-AGATCGGAAGAGCACACGTCTGAACTCCAGTCA}"
AD2="${AD2:-AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT}"

# fastp options and executable
FASTP_BIN="${FASTP_BIN:-fastp}"
FASTP_OPTS="${FASTP_OPTS:--5 -3 -W 16 -M 20 -q 20 -u 40 -l 80 --trim_poly_g --trim_poly_x}"
OUT_SUFFIX="${OUT_SUFFIX:-trim${TRIM_FRONT}}"

# ======================
# Validation & Setup
# ======================
if [[ -z "${OUTPUT_DIR}" ]]; then
  echo "ERROR: OUTPUT_DIR is required. Please set via --export=OUTPUT_DIR=..." >&2
  exit 2
fi

mkdir -p "logs" "${OUTPUT_DIR}/logs"

: "${SLURM_ARRAY_TASK_ID:?Submit with --array=0-(N-1)}"

# Build filelist
LIST="${OUTPUT_DIR}/filelist_r1.txt"

if [[ -n "${MANIFEST}" ]]; then
  [[ -s "${MANIFEST}" ]] || { echo "ERROR: MANIFEST not found: ${MANIFEST}" >&2; exit 2; }
  MODE="manifest"
  N_FILES=$(wc -l < "${MANIFEST}")
else
  [[ -n "${INPUT_DIR}" ]] || { echo "ERROR: Set MANIFEST or INPUT_DIR" >&2; exit 2; }
  [[ -d "${INPUT_DIR}" ]] || { echo "ERROR: INPUT_DIR not found: ${INPUT_DIR}" >&2; exit 2; }
  MODE="auto"
  find -L "${INPUT_DIR}" -maxdepth 1 -type f -name "${PATTERN_R1}" | sort > "${LIST}"
  N_FILES=$(wc -l < "${LIST}")
fi

if [[ "${SLURM_ARRAY_TASK_ID}" -ge "${N_FILES}" ]]; then
  echo "ERROR: Array index ${SLURM_ARRAY_TASK_ID} out of range (N=${N_FILES})" >&2
  exit 2
fi

# Get sample files
if [[ "${MODE}" == "manifest" ]]; then
  LINE="$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "${MANIFEST}")"
  SAMPLE="$(echo "$LINE" | cut -f1)"
  R1="$(echo "$LINE" | cut -f2)"
  R2="$(echo "$LINE" | cut -f3)"
else
  R1="$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "${LIST}")"
  R2="${R1/${R2_TOKEN_FROM_R1}/${R2_TOKEN}}"
  SAMPLE="$(basename "${R1}")"
  SAMPLE="${SAMPLE%%${R2_TOKEN_FROM_R1}*}"
fi

[[ -r "${R1}" ]] || { echo "ERROR: R1 not readable: ${R1}" >&2; exit 2; }
[[ -r "${R2}" ]] || { echo "ERROR: R2 not readable: ${R2}" >&2; exit 2; }

# Output paths
OUT_R1="${OUTPUT_DIR}/${SAMPLE}_R1.${OUT_SUFFIX}.fq.gz"
OUT_R2="${OUTPUT_DIR}/${SAMPLE}_R2.${OUT_SUFFIX}.fq.gz"
JSON="${OUTPUT_DIR}/${SAMPLE}.${OUT_SUFFIX}.fastp.json"
HTML="${OUTPUT_DIR}/${SAMPLE}.${OUT_SUFFIX}.fastp.html"
LOG="${OUTPUT_DIR}/logs/${SAMPLE}.${OUT_SUFFIX}.log"

# Resume protection
if [[ -s "${OUT_R1}" && -s "${OUT_R2}" && -s "${HTML}" ]]; then
  echo "Already done: ${SAMPLE}"
  exit 0
fi

echo "[$(date)] fastp trim_front=${TRIM_FRONT} | sample=${SAMPLE}"

SECONDS=0
{
  ${FASTP_BIN} \
    -i "${R1}" -I "${R2}" \
    -o "${OUT_R1}" -O "${OUT_R2}" \
    --adapter_sequence "${AD1}" \
    --adapter_sequence_r2 "${AD2}" \
    --trim_front1 "${TRIM_FRONT}" --trim_front2 "${TRIM_FRONT}" \
    ${FASTP_OPTS} \
    --thread "${SLURM_CPUS_PER_TASK}" \
    -j "${JSON}" -h "${HTML}"
} &> "${LOG}"

echo "[$(date)] Done: ${SAMPLE} (${SECONDS}s)"