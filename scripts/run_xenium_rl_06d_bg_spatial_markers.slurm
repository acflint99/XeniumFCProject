#!/bin/bash
# Submit with a manifest-derived --array; see submit_xenium_rl_bg_spatial_markers.sh.
#SBATCH --job-name=Xen_RL_BG_diag
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/%x_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/%x_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --partition=express
#SBATCH --constraint=intel
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND
set -euo pipefail
module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject || exit 1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
: "${SLURM_ARRAY_TASK_ID:?Submit using the manifest-derived array helper}"
driver=scripts/xenium_rl_06d_bg_spatial_markers.R
Rscript --max-ppsize=500000 "$driver" --dry-run "$SLURM_ARRAY_TASK_ID"
driver_options=("$SLURM_ARRAY_TASK_ID")
if [[ "${RL_BG_DIAGNOSTICS_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi
Rscript --max-ppsize=500000 "$driver" "${driver_options[@]}"
