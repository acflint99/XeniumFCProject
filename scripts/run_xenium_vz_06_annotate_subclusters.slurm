#!/bin/bash
#SBATCH --job-name=Xen_VZ_Subclusters
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_VZ_Subclusters_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_VZ_Subclusters_%j.err
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
if [[ "${VZ_SUBCLUSTER_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

case "${VZ_SUBCLUSTER_MODE:-review}" in
  review)
    driver_options+=(--review-only)
    ;;
  apply)
    driver_options+=(--apply-reviewed-labels)
    ;;
  *)
    echo "VZ_SUBCLUSTER_MODE must be 'review' or 'apply'." >&2
    exit 2
    ;;
esac

Rscript scripts/xenium_vz_06_annotate_subclusters.R "${driver_options[@]}"
