#!/bin/bash
#SBATCH --job-name=Science_res3_top20
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Science_res3_top20_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Science_res3_top20_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --partition=short
#SBATCH --constraint=intel
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND

# Full Science reprocessing and marker testing: estimated 1-8 h, excluding queue.
# Very low confidence: 70,731 cells in the prior report; no measured full-run history.
# Keep 256 GB from the existing Science diagnostic envelope; collect MaxRSS to refine.
# Time limit includes headroom; short allows 12 h per current Cheaha hardware docs
# https://docs.rc.uab.edu/cheaha/hardware/ (checked 2026-10-21).
# One Science reference object, not an array. No sample/task mapping.
set -eo pipefail
module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject || exit 1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

driver_options=()
case "${SCIENCE_RES3_OVERWRITE:-false}" in
  false) ;;
  true) driver_options+=(--overwrite) ;;
  *) echo "SCIENCE_RES3_OVERWRITE must be true or false." >&2; exit 2 ;;
esac

echo "Science resolution 3 / top 20: job=${SLURM_JOB_ID:-manual}; start=$(date -Is)"
# Project-root startup activates existing renv; expose loading errors, never install.
Rscript --max-ppsize=500000 -e 'cat(R.version.string, "\n"); print(.libPaths()); library(Seurat); library(harmony); library(dplyr); library(ggplot2)'
Rscript --max-ppsize=500000 'scripts/fetal cerebellum_science_ACF.R' --resolution=3 --dry-run
Rscript --max-ppsize=500000 'scripts/fetal cerebellum_science_ACF.R' --resolution=3 "${driver_options[@]}"
echo "Science resolution 3 / top 20 completed: end=$(date -Is)"
