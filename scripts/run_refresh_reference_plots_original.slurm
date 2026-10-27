#!/bin/bash
#SBATCH --mail-type=FAIL,INVALID_DEPEND
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --job-name=Ref_OriginalPlots
#SBATCH --array=1-3%1
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Ref_OriginalPlots_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Ref_OriginalPlots_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --partition=express
#SBATCH --constraint=intel
set -euo pipefail
module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
Rscript --max-ppsize=500000 scripts/refresh_reference_plots_original.R --dry-run "${SLURM_ARRAY_TASK_ID}"
if [[ "${REFERENCE_PLOTS_OVERWRITE:-false}" != true ]]; then
  echo "Set REFERENCE_PLOTS_OVERWRITE=true for the requested replacement of existing plots." >&2
  exit 2
fi
Rscript --max-ppsize=500000 scripts/refresh_reference_plots_original.R --overwrite "${SLURM_ARRAY_TASK_ID}"
