#!/usr/bin/env bash
#SBATCH --job-name=gsnap-map
#SBATCH --partition=CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=40G
#SBATCH --time=3-12:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

: "${SLURM_ARRAY_TASK_ID:?Submit with --array=0-(N-1)}"

# ======================
# User-configurable variables
# ======================
READS="${READS:-}"
BASE="${BASE:-}"
DB_DIR="${DB_DIR:-}"
DB_NAME="${DB_NAME:-}"
REF_TAG="${REF_TAG:-$DB_NAME}"
RUN_TAG="${RUN_TAG:-${REF_TAG}_gsnap}"

if [[ -z "$READS" || -z "$BASE" || -z "$DB_DIR" || -z "$DB_NAME" ]]; then
  echo "[ERR] Required env vars missing: READS, BASE, DB_DIR, DB_NAME" >&2
  exit 2
fi

# Performance & Environment Settings
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-20}}"
SORT_THREADS="${SORT_THREADS:-4}"
SORT_MEM="${SORT_MEM:-3G}"
TMPDIR="${TMPDIR:-/tmp/${USER}_gsnap_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}}"
mkdir -p "$TMPDIR" "logs" "$BASE"

# GSNAP Specific Alignment Knobs
LOCAL_SPLICE="${LOCAL_SPLICE:-5000}"
MAX_INTRON_MID="${MAX_INTRON_MID:-20000}"
MAX_INTRON_END="${MAX_INTRON_END:-5000}"
PAIRMAX_RNA="${PAIRMAX_RNA:-1000}"
NPATHS="${NPATHS:-10}"
NOVEL="${NOVEL:-1}"
ONLY_CONCORDANT="${ONLY_CONCORDANT:-1}"

# Optional: Load HPC modules if they exist
module purge 2>/dev/null || true
module load gmap/2024-05-20 2>/dev/null || true
module load samtools 2>/dev/null || true

# Executables
GSNAP_BIN="${GSNAP_BIN:-gsnap}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}"

# Output Layout
OUT="${BASE}/${RUN_TAG}"
OUT_BAM="${OUT}/bam"
OUT_LOG="${OUT}/logs"
OUT_QC="${OUT}/qc"
mkdir -p "$OUT_BAM" "$OUT_LOG" "$OUT_QC"

LIST="${OUT}/r1_list.txt"

# Helper functions
pair_r2 () {
  local r1="$1"
  local d="$(dirname "$r1")"
  local b="$(basename "$r1")"
  local cand="${b/R1/R2}"
  if [[ "$cand" != "$b" ]] && [[ -r "${d}/${cand}" ]]; then
      echo "${d}/${cand}"
      return 0
  fi
  return 1
}

sample_from_r1 () {
  local b="$1"
  b="${b%.fastq.gz}"; b="${b%.fq.gz}"
  b="$(echo "$b" | sed -E 's/(_R?1(_001)?|[._]R1|[._]1)([._-].*)?$//')"
  echo "$b" | sed -E 's/[-._]$//'
}

# Build R1 list
if [[ ! -s "$LIST" ]]; then
  find -L "$READS" -maxdepth 1 -type f \( -name "*.fq.gz" -o -name "*.fastq.gz" \) | grep -E '(_R?1[._-]|_R?1$|[._]1[._-])' | sort > "$LIST"
fi

R1="$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$LIST")"
[[ -n "$R1" ]] || { echo "[ERR] Empty R1 at index $SLURM_ARRAY_TASK_ID" >&2; exit 4; }

# Safe R2 pairing and validation
R2="$(pair_r2 "$R1" || true)"
[[ -n "$R2" ]] || { echo "[ERR] Cannot find matching R2 for R1: $R1" >&2; exit 5; }
[[ -r "$R2" ]] || { echo "[ERR] R2 not readable: $R2" >&2; exit 5; }

SAMPLE="$(sample_from_r1 "$(basename "$R1")")"
BAM="${OUT_BAM}/${SAMPLE}.bam"
LOG="${OUT_LOG}/${SAMPLE}.gsnap.log"

if [[ -s "$BAM" && -s "${BAM}.bai" ]]; then
  echo "Already done: $SAMPLE"
  rm -rf "$TMPDIR"
  exit 0
fi

# Execution
help_txt="$("$GSNAP_BIN" --help 2>&1 || true)"
MID_FLAG="$(grep -q -- "--max-middle-intronlength" <<< "$help_txt" && echo "--max-middle-intronlength" || echo "--max-intronlength-middle")"
END_FLAG="$(grep -q -- "--max-end-intronlength" <<< "$help_txt" && echo "--max-end-intronlength" || echo "--max-intronlength-ends")"

cmd=("$GSNAP_BIN" -D "$DB_DIR" -d "$DB_NAME" --gunzip --nthreads "${THREADS}" -A sam \
      --localsplicedist "${LOCAL_SPLICE}" --pairmax-rna "${PAIRMAX_RNA}" \
      --novelsplicing "${NOVEL}" --npaths "${NPATHS}" --quiet-if-excessive \
      --read-group-id "${SAMPLE}" --read-group-name "${SAMPLE}" --read-group-platform "ILLUMINA")

[[ "$ONLY_CONCORDANT" == "1" ]] && cmd+=(--only-concordant)
[[ -n "$MID_FLAG" ]] && cmd+=("${MID_FLAG}=${MAX_INTRON_MID}")
[[ -n "$END_FLAG" ]] && cmd+=("${END_FLAG}=${MAX_INTRON_END}")
cmd+=("$R1" "$R2")

set +e
"${cmd[@]}" 2> "$LOG" | "$SAMTOOLS_BIN" sort -@ "${SORT_THREADS}" -m "${SORT_MEM}" -T "${TMPDIR}/${SAMPLE}_sorttmp" -o "$BAM" -
p0=${PIPESTATUS[0]}
p1=${PIPESTATUS[1]}
set -e

if [[ $p0 -ne 0 || $p1 -ne 0 ]]; then
  echo "Mapping failed for $SAMPLE"
  rm -rf "$TMPDIR"
  exit 1
fi

"$SAMTOOLS_BIN" index "$BAM"

# Basic QC & TSV Generation
tot=$(zcat "$R1" | awk 'END{printf "%.0f", NR/4}')
map_primary_pp=$("$SAMTOOLS_BIN" view -c -f 2 -F 3852 "$BAM" || echo 0)
map_all=$("$SAMTOOLS_BIN" view -c -F 4 "$BAM" || echo 0)
rate_primary_pp=$(awk -v m="$map_primary_pp" -v t="$tot" 'BEGIN{printf (t? "%.2f":"0"), 100*m/t}')
rate_all=$(awk -v m="$map_all" -v t="$tot" 'BEGIN{printf (t? "%.2f":"0"), 100*m/t}')

"$SAMTOOLS_BIN" flagstat "$BAM" > "${OUT_QC}/${SAMPLE}.flagstat.txt"
{
  echo -e "sample_id\treads_total_R1\tmapped_primary_properly\tmapped_rate_primary_properly(%)\tmapped_all\tmapped_rate_all(%)"
  echo -e "${SAMPLE}\t${tot}\t${map_primary_pp}\t${rate_primary_pp}\t${map_all}\t${rate_all}"
} > "${OUT_QC}/${SAMPLE}.${REF_TAG}.mapping.tsv"

# Cleanup
rm -rf "$TMPDIR"
echo "Finished mapping: $SAMPLE"