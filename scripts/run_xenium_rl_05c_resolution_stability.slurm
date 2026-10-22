#!/bin/bash
#SBATCH --job-name=RL_ResolutionStability
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/RL_ResolutionStability_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/RL_ResolutionStability_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=600G
#SBATCH --partition=medium
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

# .Rprofile activates the existing project renv; never install inside this job.
driver_options=(--region=RL "--repeats=${STABILITY_REPEATS:-5}"
  "--fraction=${STABILITY_FRACTION:-0.8}" "--run-id=${STABILITY_RUN_ID:-pilot5}")
if [[ "${STABILITY_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
else
  Rscript --max-ppsize=500000 scripts/xenium_regional_05c_resolution_stability.R \
    --region=RL --repeats="${STABILITY_REPEATS:-5}" \
    --fraction="${STABILITY_FRACTION:-0.8}" --run-id="${STABILITY_RUN_ID:-pilot5}" --dry-run
fi
Rscript --max-ppsize=500000 scripts/xenium_regional_05c_resolution_stability.R "${driver_options[@]}"
