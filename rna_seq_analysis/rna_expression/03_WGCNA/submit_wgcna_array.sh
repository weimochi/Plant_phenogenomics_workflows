#!/bin/bash
#SBATCH --job-name=WGCNA_CC_Array
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --export=ALL

# ==============================================================================
# [Super Matrix] 6 Powers × 2 Merge Thresholds = 12 Parallel Array Tasks!
# ==============================================================================
#SBATCH --array=1-12
#SBATCH --output=wgcna_array_%A_%a.log

. /etc/profile
. ~/.bashrc
module load R 2>/dev/null || module load r 2>/dev/null || echo "[INFO] No module system found, trying path..."

# ==============================================================================
# [IMPORTANT] Please update these paths to match your current project directory!
# ==============================================================================
INPUT_RDS="/path/to/your/output_dir_from_step01/datExpr_ready_for_WGCNA.rds"
OUT_DIR="/path/to/your/desired_wgcna_output_dir"
SCRIPT_PATH="/path/to/your/scripts_dir/02_run_wgcna_modules.R"

# 1. Define all Powers and Merge thresholds you want to test
POWERS=(14 14 16 16 18 18 20 20 22 22 24 24)
MERGES=(0.25 0.20 0.25 0.20 0.25 0.20 0.25 0.20 0.25 0.20 0.25 0.20)

# 2. Fetch the corresponding values based on the current Slurm array task ID (1 to 12)
# (Since Bash arrays are 0-indexed, we subtract 1)
IDX=$((SLURM_ARRAY_TASK_ID - 1))

CURRENT_POWER=${POWERS[$IDX]}
CURRENT_MERGE=${MERGES[$IDX]}

echo "==== [Slurm Task $SLURM_ARRAY_TASK_ID Started] ===="
echo "Execution Time: $(date)"
echo "Array Index:    Position $IDX"
echo "Parameters:     Power = $CURRENT_POWER , Merge Cut = $CURRENT_MERGE"
echo "Input RDS:      $INPUT_RDS"
echo "Output Dir:     $OUT_DIR"
echo "================================================="

# 3. Run the R Pipeline
Rscript $SCRIPT_PATH $INPUT_RDS $CURRENT_POWER $CURRENT_MERGE $OUT_DIR

echo "==== [Slurm Task $SLURM_ARRAY_TASK_ID Finished]: $(date) ===="