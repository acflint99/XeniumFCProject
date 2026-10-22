#!/bin/bash
#SBATCH --job-name=VZ_ResolutionCompare
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/VZ_ResolutionCompare_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/VZ_ResolutionCompare_%j.err
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
# Keep one worker: VZ parallel object duplication previously caused OOMs.
# The project .Rprofile activates the existing renv; no installation in jobs.
driver_options=("--run-id=${VZ_RESOLUTION_COMPARE_RUN_ID:-review1}")
if [[ "${VZ_RESOLUTION_COMPARE_RESUME:-false}" == "true" ]]; then
  driver_options+=(--resume-checkpoint)
fi
if [[ "${VZ_RESOLUTION_COMPARE_TIFF:-false}" == "true" ]]; then
  driver_options+=(--tiff)
fi
if [[ "${VZ_RESOLUTION_COMPARE_OVERWRITE:-false}" == "true" ]]; then
  Rscript --max-ppsize=500000 scripts/xenium_vz_05d_compare_resolutions.R --list
  driver_options+=(--overwrite)
else
  Rscript --max-ppsize=500000 scripts/xenium_vz_05d_compare_resolutions.R "${driver_options[@]}" --dry-run
fi
Rscript --max-ppsize=500000 scripts/xenium_vz_05d_compare_resolutions.R "${driver_options[@]}"
