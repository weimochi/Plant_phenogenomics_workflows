#!/usr/bin/env bash
#SBATCH --job-name=star-map
#SBATCH --partition=CPU
#SBATCH --cpus-per-task=24
#SBATCH --mem=40G
#SBATCH --time=3-12:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail
: "${SLURM_ARRAY_TASK_ID:?Submit with --array=0-(N-1)}"

# ============================================================
# User Configurable Variables
# ============================================================
READS="${READS:-}"
BASE="${BASE:-}"
GENOME_DIR="${GENOME_DIR:-}"
REF_TAG="${REF_TAG:-$(basename "$GENOME_DIR")}"
RUN_TAG="${RUN_TAG:-${REF_TAG}_star}"

if [[ -z "$READS" || -z "$BASE" || -z "$GENOME_DIR" ]]; then
  echo "[ERR] Missing required vars: READS / BASE / GENOME_DIR. Set via --export" >&2
  exit 2
fi

# Performance & Params
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-20}}"
SORT_THREADS="${SORT_THREADS:-4}"
TWOPASS_MODE="${TWOPASS_MODE:-Basic}"
BAMSORT_RAM="${BAMSORT_RAM:-20000000000}"

# STAR Specific Alignment Knobs
OUT_FILTER_MULTIMAP_NMAX="${OUT_FILTER_MULTIMAP_NMAX:-10}"
ALIGN_INTRON_MAX="${ALIGN_INTRON_MAX:-20000}"
ALIGN_MATES_GAP_MAX="${ALIGN_MATES_GAP_MAX:-20000}"
STAR_EXTRA="${STAR_EXTRA:-}"

# Executables & Modules
module purge 2>/dev/null || true
module load STAR 2>/dev/null || true
module load samtools 2>/dev/null || true
STAR_BIN="${STAR_BIN:-STAR}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}"

# Directories
OUT="${BASE}/${RUN_TAG}"
OUT_BAM="${OUT}/bam"
OUT_LOG="${OUT}/logs"
OUT_QC="${OUT}/qc"
OUT_TMP="${OUT}/tmp"
mkdir -p "$OUT_BAM" "$OUT_LOG" "$OUT_QC" "$OUT_TMP" "logs"

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

# Build file list dynamically
if [[ ! -s "$LIST" ]]; then
  find -L "$READS" -maxdepth 1 -type f \( -name "*.fq.gz" -o -name "*.fastq.gz" \) | grep -E '(_R?1[._-]|_R?1$|[._]1[._-])' | sort > "$LIST"
fi

R1="$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$LIST")"
[[ -n "$R1" ]] || { echo "[ERR] Empty R1 at index $SLURM_ARRAY_TASK_ID" >&2; exit 4; }

R2="$(pair_r2 "$R1" || true)"
[[ -n "$R2" ]] || { echo "[ERR] Cannot find matching R2 for R1: $R1" >&2; exit 5; }
[[ -r "$R2" ]] || { echo "[ERR] R2 not readable: $R2" >&2; exit 5; }

SAMPLE="$(basename "$R1" | sed -E 's/(_R?1(_001)?|[._]R1)(\.fq\.gz|\.fastq\.gz)?//')"
BAM="${OUT_BAM}/${SAMPLE}.bam"
STAR_PREFIX="${OUT_BAM}/${SAMPLE}."
STAR_TMP="${OUT_TMP}/${SAMPLE}_star_tmp"
LOG="${OUT_LOG}/${SAMPLE}.star.log"

if [[ -s "$BAM" && -s "${BAM}.bai" ]]; then
  echo "⏭️ already done: $SAMPLE"
  exit 0
fi

# Execution
rm -f "$BAM" "${BAM}.bai" "${STAR_PREFIX}"* || true
rm -rf "$STAR_TMP" || true

set +e
"$STAR_BIN" \
  --runThreadN "${THREADS}" \
  --genomeDir "$GENOME_DIR" \
  --readFilesIn "$R1" "$R2" \
  --readFilesCommand zcat \
  --outFileNamePrefix "$STAR_PREFIX" \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattrRGline ID:"$SAMPLE" SM:"$SAMPLE" LB:lib1 PL:ILLUMINA \
  --twopassMode "${TWOPASS_MODE}" \
  --limitBAMsortRAM "${BAMSORT_RAM}" \
  --outFilterMultimapNmax "${OUT_FILTER_MULTIMAP_NMAX}" \
  --alignIntronMax "${ALIGN_INTRON_MAX}" \
  --alignMatesGapMax "${ALIGN_MATES_GAP_MAX}" \
  --outTmpDir "$STAR_TMP" \
  ${STAR_EXTRA} \
  > "$LOG" 2>&1
p0=$?
set -e

if [[ $p0 -ne 0 ]]; then
  echo "STAR FAILED for $SAMPLE"
  rm -rf "$STAR_TMP" || true
  exit 1
fi

mv "${STAR_PREFIX}Aligned.sortedByCoord.out.bam" "$BAM"
"$SAMTOOLS_BIN" index -@ "$SORT_THREADS" "$BAM"

# Clean up STAR specific logs to QC folder
for f in Log.final.out Log.out Log.progress.out SJ.out.tab; do
  [[ -f "${STAR_PREFIX}${f}" ]] && mv "${STAR_PREFIX}${f}" "${OUT_QC}/${SAMPLE}.${f}"
done

# QC Stats
tot=$(zcat "$R1" | awk 'END{printf "%.0f", NR/4}')
map_primary=$("$SAMTOOLS_BIN" view -c -f 2 -F 3852 "$BAM" || echo 0)
rate=$(awk -v m="$map_primary" -v t="$tot" 'BEGIN{printf (t? "%.2f":"0"), 100*m/t}')

echo -e "sample_id\ttotal_R1\tmapped_properly\trate(%)\n${SAMPLE}\t${tot}\t${map_primary}\t${rate}" > "${OUT_QC}/${SAMPLE}.${REF_TAG}.tsv"
"$SAMTOOLS_BIN" flagstat "$BAM" > "${OUT_QC}/${SAMPLE}.flagstat.txt"

rm -rf "$STAR_TMP" || true
echo "Done: $SAMPLE"