#!/bin/bash
#SBATCH --job-name=Science_eCN_Markers
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Science_eCN_Markers_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Science_eCN_Markers_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --partition=express
#SBATCH --constraint=intel
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND

# Expected runtime: 2-20 min including RDS load, excluding queue wait.
# Low confidence: static inspection; no measured Science eCN job history.
# Retain the existing run_xenium_diagnostic_cn_markers.slurm resource envelope.
# Related VZ stage-6 job 40966278: 31.52 min / 152.09 GB MaxRSS;
# this is supporting context, not a measurement of this Science report.
# One reference object per job; no array and no Xenium sample/task mapping.
set -eo pipefail
module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject || exit 1
# Start R at the project root for the established .Rprofile/renv activation.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

driver_options=(--ecn-only)
case "${SCIENCE_ECN_OVERWRITE:-false}" in
  false) ;;
  true) driver_options+=(--overwrite) ;;
  *) echo "SCIENCE_ECN_OVERWRITE must be true or false." >&2; exit 2 ;;
esac

echo "Science raw-cluster eCN report: job=${SLURM_JOB_ID:-manual}; start=$(date -Is)"
# Show the actual loading error if Seurat or an imported dependency is broken.
# Never install packages or restore the environment within this job.
Rscript --max-ppsize=500000 -e 'cat(R.version.string, "\n"); print(.libPaths()); library(Seurat)'
Rscript --max-ppsize=500000 'scripts/fetal cerebellum_science_ACF.R' --ecn-only --dry-run
Rscript --max-ppsize=500000 'scripts/fetal cerebellum_science_ACF.R' "${driver_options[@]}"
echo "Science raw-cluster eCN report completed: end=$(date -Is)"
