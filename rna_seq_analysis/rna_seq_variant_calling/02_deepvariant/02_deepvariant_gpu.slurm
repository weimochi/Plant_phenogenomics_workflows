#!/usr/bin/env bash
#SBATCH --job-name=dv_gpu
#SBATCH --partition=GPU
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --time=4-00:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -Eeuo pipefail
: "${SLURM_ARRAY_TASK_ID:?Submit with --array=0-(N-1)}"

# =====================================================
# 1. User Configurable Variables (Passed via --export)
# =====================================================
BAMDIR="${BAMDIR:-}"
REF="${REF:-}"
OUTDIR="${OUTDIR:-}"
TAG="${TAG:-}"                               # 預設為空字串，即不使用 TAG
MODEL="${MODEL:-WES}"
BIN_VERSION="${BIN_VERSION:-1.9.0}"
FAST_STORAGE="${FAST_STORAGE:-/tmp}"         # 高速暫存路徑，預設使用系統 /tmp
EXTRA_BINDS="${EXTRA_BINDS:-}"               # 自訂額外的容器掛載點 (例如 "-B /nas -B /data")

if [[ -z "$BAMDIR" || -z "$REF" || -z "$OUTDIR" ]]; then
  echo "[ERR] Missing required vars: BAMDIR, REF, OUTDIR" >&2
  exit 1
fi

mkdir -p "$OUTDIR" "logs"

# =====================================================
# 2. Modules & Environment
# =====================================================
module purge || true
module load singularity/3.10.4 2>/dev/null || module load apptainer 2>/dev/null || true
CNT="$(command -v singularity || command -v apptainer)"
module load samtools/1.16.1 2>/dev/null || true

# =====================================================
# 3. SSD Temp Directories Setup (Crucial for GPU IO)
# =====================================================
# 使用 JOB_ID 與 TASK_ID 建立獨立暫存，避免陣列任務互相覆蓋
MY_SSD_TMP="${FAST_STORAGE}/${USER}_dv_tmp_${SLURM_JOB_ID:-$$}_${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$MY_SSD_TMP"
export TMPDIR="$MY_SSD_TMP"
export TEMP="$MY_SSD_TMP"
export TMP="$MY_SSD_TMP"
export SINGULARITY_TMPDIR="$MY_SSD_TMP"
export APPTAINER_TMPDIR="$MY_SSD_TMP"

# =====================================================
# 4. Target Selection via Array
# =====================================================
[[ -s "${REF}.fai" ]] || samtools faidx "$REF"

# 動態處理 TAG 綴詞：若 TAG 有值則為 ".TAG"，否則為空字串
TAG_SUFFIX="${TAG:+.${TAG}}"

LIST="${OUTDIR}/bam_list${TAG_SUFFIX}.txt"
if [[ ! -s "$LIST" ]]; then
  # 【修正 1】防競爭危害 (Race Condition) 的安全寫入法
  # 找尋 *.bam 確保無論 Step 1 有無 md TAG 都能抓到
  TMP_LIST="${LIST}.${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.tmp"
  find -L "$BAMDIR" -maxdepth 1 -type f -name "*.bam" | sort > "$TMP_LIST"
  mv -n "$TMP_LIST" "$LIST" 2>/dev/null || rm -f "$TMP_LIST"
  sleep 2 # 給網路硬碟 (NFS) 一點時間同步
fi

BAM="$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$LIST")"
[[ -n "$BAM" ]] || { echo "[ERR] Empty BAM at index $SLURM_ARRAY_TASK_ID" >&2; exit 2; }

# 【修正 2】無敵檔名剝除法，確保無論副檔名是什麼都能切得乾乾淨淨
SAMPLE="$(basename "$BAM")"
SAMPLE="${SAMPLE%%.md.bam}"
SAMPLE="${SAMPLE%%.bam}"

VCF="${OUTDIR}/${SAMPLE}${TAG_SUFFIX}.vcf.gz"
GVCF="${OUTDIR}/${SAMPLE}${TAG_SUFFIX}.g.vcf.gz"

if [[ -s "$VCF" && -s "$GVCF" ]]; then
  echo "Already done: $SAMPLE"
  rm -rf "$MY_SSD_TMP"
  exit 0
fi

[[ -s "${BAM}.bai" ]] || samtools index -@ 4 "$BAM"

# =====================================================
# 5. DeepVariant Execution
# =====================================================
if [[ -n "$TAG" ]]; then
  echo "[RUN][GPU] $SAMPLE | Tag=$TAG | Model=$MODEL"
else
  echo "[RUN][GPU] $SAMPLE | Model=$MODEL"
fi

IMG="docker://google/deepvariant:${BIN_VERSION}-gpu"

# 動態組裝掛載參數：必定掛載 FAST_STORAGE，如果有給 EXTRA_BINDS 就接在後面
BIND_ARGS="-B ${FAST_STORAGE}"
[[ -n "$EXTRA_BINDS" ]] && BIND_ARGS="${BIND_ARGS} ${EXTRA_BINDS}"

"$CNT" run --nv $BIND_ARGS \
  "$IMG" /opt/deepvariant/bin/run_deepvariant \
    --model_type="${MODEL}" \
    --ref="$REF" \
    --reads="$BAM" \
    --output_vcf="$VCF" \
    --output_gvcf="$GVCF" \
    --sample_name="$SAMPLE" \
    --vcf_stats_report=true \
    --num_shards="${SLURM_CPUS_PER_TASK:-20}" \
    --intermediate_results_dir "$MY_SSD_TMP"

# 清理 SSD 暫存
rm -rf "$MY_SSD_TMP"
echo "Finished: $SAMPLE -> $OUTDIR"