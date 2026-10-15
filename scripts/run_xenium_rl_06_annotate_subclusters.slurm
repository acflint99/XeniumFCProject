#!/bin/bash
#SBATCH --job-name=Xen_RL_Subclusters
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_RL_Subclusters_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_RL_Subclusters_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --partition=short
#SBATCH --constraint=intel
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND

module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject || exit 1

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

driver_options=()
subcluster_resolution="${RL_SUBCLUSTER_RESOLUTION:-0.5}"
case "${subcluster_resolution}" in
  0.3|0.5|0.8)
    driver_options+=("--resolution=${subcluster_resolution}")
    ;;
  *)
    echo "RL_SUBCLUSTER_RESOLUTION must be 0.3, 0.5, or 0.8." >&2
    exit 2
    ;;
esac
if [[ "${RL_SUBCLUSTER_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

case "${RL_SUBCLUSTER_MODE:-review}" in
  review)
    driver_options+=(--review-only)
    ;;
  apply)
    driver_options+=(--apply-reviewed-labels)
    ;;
  refresh-plots)
    driver_options+=(--plots-only)
    ;;
  *)
    echo "RL_SUBCLUSTER_MODE must be 'review', 'apply', or 'refresh-plots'." >&2
    exit 2
    ;;
esac

Rscript --max-ppsize=500000 scripts/xenium_rl_06_annotate_subclusters.R "${driver_options[@]}"
