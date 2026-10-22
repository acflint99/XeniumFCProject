#!/bin/bash
#SBATCH --job-name=Xen_CN_Markers
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_CN_Markers_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_CN_Markers_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --partition=express
#SBATCH --constraint=intel
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND

# Expected merged report: 15-90 minutes, excluding queue wait; low confidence.
# Resource basis: same-merge expression export 41005351 took 3.45 min / 22.86 GB;
# processed merged plotting 41006727 reached 24.88 min / 182.34 GB before failure.
# These are different workloads; retain 256 GB with headroom, no marker tests.
# One combined cohort per job; no sample array or task-to-sample mapping.
set -eo pipefail
module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject || exit 1
# Starting R here lets the project .Rprofile activate its existing renv.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

branch="${CN_MARKER_BRANCH:-nonspatial}"
case "${branch}" in
  nonspatial|vz|rl|residual) ;;
  *) echo "CN_MARKER_BRANCH must be nonspatial, vz, rl, or residual." >&2; exit 2 ;;
esac
driver_options=("--branch=${branch}")
case "${CN_MARKER_OVERWRITE:-false}" in
  false) ;;
  true) driver_options+=(--overwrite) ;;
  *) echo "CN_MARKER_OVERWRITE must be true or false." >&2; exit 2 ;;
esac

echo "CN marker report: branch=${branch}; job=${SLURM_JOB_ID:-manual}; start=$(date -Is)"
# Expose the underlying package-loading error if Seurat or an import is broken.
# No installation, restore, or environment changes are performed.
Rscript --max-ppsize=500000 -e 'cat(R.version.string, "\n"); print(.libPaths()); library(Seurat)'
Rscript --max-ppsize=500000 scripts/xenium_diagnostic_cn_markers.R "--branch=${branch}" --dry-run
Rscript --max-ppsize=500000 scripts/xenium_diagnostic_cn_markers.R "${driver_options[@]}"
echo "CN marker report completed: branch=${branch}; end=$(date -Is)"
