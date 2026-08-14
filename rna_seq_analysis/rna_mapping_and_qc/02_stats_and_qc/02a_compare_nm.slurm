#!/usr/bin/env bash
#SBATCH --job-name=compare_nm
#SBATCH --partition=CPU
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail
LC_ALL=C

# ============================================================
# NM Mismatch & Coverage Multi-Genome Competitive Analysis
# ============================================================
# User configurable variables (passed via --export)
GENOMES="${GENOMES:-}"         # Format: "Ref1:/path/1 Ref2:/path/2 Ref3:/path/3 ..."
OUTDIR="${OUTDIR:-}"           # Output directory for the TSV report
MAPQ="${MAPQ:-0}"              # Minimum MAPQ filter
THREADS="${THREADS:-2}"        # samtools threads per job
JOBS="${JOBS:-4}"              # Parallel sample jobs (JOBS * THREADS <= SLURM_CPUS_PER_TASK)
BAM_SUFFIX="${BAM_SUFFIX:-.bam}"

if [[ -z "$GENOMES" || -z "$OUTDIR" ]]; then
  echo "[ERR] Missing required vars: GENOMES, OUTDIR. Set via --export" >&2
  exit 2
fi

mkdir -p "$OUTDIR" "logs"
module load samtools 2>/dev/null || true

# Parse GENOMES string into arrays
declare -A GENOME_DIRS
GENOME_ORDER=()
for g in $GENOMES; do
  gname="${g%%:*}"
  gdir="${g#*:}"
  [[ -d "$gdir" ]] || { echo "[ERR] Genome dir not found: $gdir" >&2; exit 2; }
  GENOME_DIRS["$gname"]="$gdir"
  GENOME_ORDER+=("$gname")
done

N_GENOMES=${#GENOME_ORDER[@]}
STAMP=$(date +%Y%m%d_%H%M)
LABEL=$(IFS=_; echo "${GENOME_ORDER[*]}")
OUT="${OUTDIR}/${LABEL}_compare_${STAMP}.tsv"
TMPDIR="${OUTDIR}/.nm_tmp_${STAMP}"
mkdir -p "$TMPDIR"

echo "[INFO] Genomes (${N_GENOMES}) : ${GENOME_ORDER[*]}"
echo "[INFO] Output          : $OUT"

# ============================================================
# Worker Function (Exported for parallel execution)
# ============================================================
process_sample() {
  local sample="$1" mapq="$2" thr="$3" suffix="$4" tmpdir="$6" n_genomes="$7"
  shift 7
  local genome_order=()
  declare -A genome_dirs
  while [[ $# -ge 2 ]]; do
    genome_order+=("$1")
    genome_dirs["$1"]="$2"
    shift 2
  done

  _nm_rate() {
    local bam="$1"
    [[ -s "$bam" ]] || { echo "NA"; return; }
    samtools view -@ "$thr" -F 0x904 -q "$mapq" "$bam" | awk '{
        L=length($10)
        for(i=12;i<=NF;i++){
          if($i ~ /^(NM|nM):i:/){ split($i,a,":"); nm += a[3]; tb += L; break }
        }
      } END{ if(tb>0) printf "%.4f\n", 100*nm/tb; else print "NA" }'
  }

  _cov_stats() {
    local bam="$1"
    [[ -s "$bam" ]] || { printf "NA\tNA\tNA\tNA"; return; }
    samtools coverage --min-MQ "$mapq" "$bam" 2>/dev/null | awk 'NR==1{next}
           { numreads+=$4; totlen+=($3-$2); wcov+=$6*($3-$2); wdepth+=$7*($3-$2); wmapq+=$9*($3-$2) }
           END{ if(totlen==0){ print "NA\tNA\tNA\tNA" }
                else{ printf "%d\t%.4f\t%.4f\t%.2f\n", numreads, wcov/totlen, wdepth/totlen, wmapq/totlen } }'
  }

  local ROW="$sample"
  declare -A nm_val
  local best_ref="NA"
  local min_nm=999999

  # 1. 依序對每個 Reference 計算 NM 與 Coverage
  for g in "${genome_order[@]}"; do
    local bam="${genome_dirs[$g]}/${sample}${suffix}"
    local nm; nm=$(_nm_rate "$bam")
    nm_val["$g"]="$nm"
    local stats; stats=$(_cov_stats "$bam")
    ROW+="\t${nm}\t${stats}"

    # 動態尋找 NM 最低（最佳匹配）的 Reference
    if [[ "$nm" != "NA" ]]; then
      # 透過 awk 比較浮點數大小
      local is_smaller
      is_smaller=$(awk -v a="$nm" -v b="$min_nm" 'BEGIN{print (a < b) ? 1 : 0}')
      if [[ "$is_smaller" -eq 1 ]]; then
        min_nm="$nm"
        best_ref="$g"
      fi
    fi
  done

  # 2. 附加最佳匹配結果（Best Hit）
  ROW+="\t${best_ref}"

  # 3. 附加各基因體的實際 BAM 路徑
  for g in "${genome_order[@]}"; do 
    ROW+="\t${genome_dirs[$g]}/${sample}${suffix}"
  done

  printf '%s\n' "$ROW" > "${tmpdir}/${sample}.row"
}
export -f process_sample

# ============================================================
# Execution
# ============================================================
GENOME_ARGS=()
for g in "${GENOME_ORDER[@]}"; do GENOME_ARGS+=("$g" "${GENOME_DIRS[$g]}"); done

REF_DIR="${GENOME_DIRS[${GENOME_ORDER[0]}]}"
mapfile -t SAMPLES < <(find "$REF_DIR" -maxdepth 1 -name "*${BAM_SUFFIX}" | LC_ALL=C sort | xargs -I{} basename {} "${BAM_SUFFIX}")

echo "[INFO] Processing ${#SAMPLES[@]} samples..."

# Use GNU parallel if available, otherwise fallback to background job loop
if command -v parallel >/dev/null 2>&1 && parallel --version 2>&1 | grep -q "GNU parallel"; then
  printf '%s\n' "${SAMPLES[@]}" | parallel --jobs "$JOBS" --bar \
    process_sample {} "$MAPQ" "$THREADS" "$BAM_SUFFIX" "" "$TMPDIR" "$N_GENOMES" "${GENOME_ARGS[@]}"
else
  run_count=0
  for sample in "${SAMPLES[@]}"; do
    process_sample "$sample" "$MAPQ" "$THREADS" "$BAM_SUFFIX" "" "$TMPDIR" "$N_GENOMES" "${GENOME_ARGS[@]}" &
    ((++run_count % JOBS == 0)) && wait
  done
  wait
fi

# Assemble Header & Results
HEADER="sample_id"
for g in "${GENOME_ORDER[@]}"; do 
  HEADER+="\t${g}_NM(%)\t${g}_numreads\t${g}_coverage\t${g}_meandepth\t${g}_meanmapq"
done

# 支援任意 N 個基因體的總結欄位：直接輸出 Best Match
HEADER+="\tbest_match"

for g in "${GENOME_ORDER[@]}"; do 
  HEADER+="\t${g}_bam"
done

{
  echo -e "$HEADER"
  for sample in "${SAMPLES[@]}"; do cat "${TMPDIR}/${sample}.row"; done
} > "$OUT"

rm -rf "$TMPDIR"
echo "Analysis complete: $OUT"