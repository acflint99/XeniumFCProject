#!/bin/bash
#SBATCH --mail-type=FAIL,INVALID_DEPEND
#SBATCH --mail-user=acflint@uab.edu
# Defaults match the existing pooled VZ/RL plotting resource envelope.
# submit_refresh_current_plots.sh sets group-specific array/resource overrides.
#SBATCH --job-name=Xen_PaletteRefresh
#SBATCH --array=1-3%1
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_PaletteRefresh_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_PaletteRefresh_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --partition=short
#SBATCH --constraint=intel

set -euo pipefail
module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
plot_group="${PLOT_REFRESH_GROUP:-pooled}"
plot_tag="${PLOT_REFRESH_TAG:-palette_20261019}"
Rscript --max-ppsize=500000 scripts/refresh_current_plots.R \
  "--group=${plot_group}" "--tag=${plot_tag}" --dry-run "${SLURM_ARRAY_TASK_ID}"
Rscript --max-ppsize=500000 scripts/refresh_current_plots.R \
  "--group=${plot_group}" "--tag=${plot_tag}" "${SLURM_ARRAY_TASK_ID}"
