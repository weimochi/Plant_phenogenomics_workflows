#!/usr/bin/env bash
#SBATCH --job-name=mosdepth-cov
#SBATCH --partition=CPU
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail
shopt -s nullglob

# ============================================================
# Mosdepth Coverage & Gene-Body Analysis
# ============================================================
BAMDIR="${BAMDIR:-}"
OUT="${OUT:-}"
GENEBED="${GENEBED:-}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
Q="${Q:-10}"

if [[ -z "$BAMDIR" || -z "$OUT" || -z "$GENEBED" ]]; then
  echo "[ERR] Required env vars missing: BAMDIR, OUT, GENEBED. Set via --export" >&2
  exit 1
fi

[[ -d "$BAMDIR" ]] || { echo "[ERR] BAMDIR not found: $BAMDIR" >&2; exit 2; }
[[ -s "$GENEBED" ]] || { echo "[ERR] GENEBED not found or empty: $GENEBED" >&2; exit 2; }
command -v mosdepth >/dev/null 2>&1 || { echo "[ERR] mosdepth not found in PATH" >&2; exit 3; }

mkdir -p "$(dirname "$OUT")" "logs"

mapfile -t BAMLIST < <(find "$BAMDIR" -maxdepth 1 -type f -name "*.bam" ! -name "*.bai" | sort)
[[ ${#BAMLIST[@]} -gt 0 ]] || { echo "[ERR] No .bam found under: $BAMDIR" >&2; exit 4; }

echo -e "sample\tmean_depth\tcoverage_1x\tgene_body_mean_depth\tbam" > "$OUT"

for bam in "${BAMLIST[@]}"; do
  sample="$(basename "$bam" .bam)"
  prefix="${BAMDIR}/${sample}"

  echo "[INFO] Running mosdepth for: $sample"

  # 1) Whole-genome / global coverage
  mosdepth -t "$THREADS" -n -x -m -Q "$Q" "$prefix" "$bam"

  summary="${prefix}.mosdepth.summary.txt"
  global="${prefix}.mosdepth.global.dist.txt"
  [[ -s "$summary" ]] || { echo "[ERR] Missing summary: $summary" >&2; exit 5; }
  [[ -s "$global"  ]] || { echo "[ERR] Missing global : $global"  >&2; exit 6; }

  mean_depth=$(awk '$1=="total"{print $4}' "$summary")
  cov1=$(awk 'BEGIN{cov=0;total=0}{total+=$2;if($1>=1)cov+=$2}END{printf(total? "%.6f":"0", cov/total)}' "$global")

  # 2) Gene-body mean depth (using BED regions)
  gb_prefix="${prefix}.genebody"
  mosdepth -t "$THREADS" -n -x -m -Q "$Q" --by "$GENEBED" "$gb_prefix" "$bam"

  gb_regions="${gb_prefix}.regions.bed.gz"
  [[ -s "$gb_regions" ]] || { echo "[ERR] Missing regions: $gb_regions" >&2; exit 7; }

  gene_body_mean_depth=$(
    zcat "$gb_regions" | awk '
      BEGIN{sum=0; len=0}
      {L=$3-$2; if(L>0){sum += L*$5; len += L}}
      END{ if(len>0) printf("%.6f", sum/len); else print "0" }'
  )

  echo -e "${sample}\t${mean_depth}\t${cov1}\t${gene_body_mean_depth}\t${bam}" >> "$OUT"
done

echo "Coverage & Gene-body analysis complete: $OUT"