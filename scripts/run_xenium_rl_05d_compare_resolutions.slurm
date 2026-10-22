#!/bin/bash
#SBATCH --job-name=RL_ResolutionCompare
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/RL_ResolutionCompare_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/RL_ResolutionCompare_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --partition=short
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
# Project .Rprofile activates renv. No package installation in production jobs.
driver_options=("--run-id=${RL_RESOLUTION_COMPARE_RUN_ID:-review1}")
if [[ "${RL_RESOLUTION_COMPARE_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
else
  Rscript --max-ppsize=500000 scripts/xenium_rl_05d_compare_resolutions.R \
    --run-id="${RL_RESOLUTION_COMPARE_RUN_ID:-review1}" --dry-run
fi
Rscript --max-ppsize=500000 scripts/xenium_rl_05d_compare_resolutions.R "${driver_options[@]}"
