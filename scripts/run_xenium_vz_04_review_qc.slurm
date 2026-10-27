#!/bin/bash
#SBATCH --job-name=Xen_VZ_QC
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_VZ_QC_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_VZ_QC_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
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
if [[ "${VZ_QC_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

case "${VZ_QC_MODE:-review}" in
  review)
    driver_options+=(--qc-only)
    ;;
  apply|crosswalk)
    if [[ -z "${VZ_QC_REMOVE_CLUSTERS:-}" ]]; then
      echo "VZ_QC_REMOVE_CLUSTERS must explicitly name the newly reviewed clusters or 'none' (current reviewed decision: 23,24,25,27,28)." >&2
      exit 2
    fi
    if [[ "${VZ_QC_MODE}" == "crosswalk" ]]; then
      driver_options+=(--crosswalk-only)
      if [[ -z "${VZ_QC_ADDITIONAL_POSTQC_CLUSTERS:-}" ]]; then
        echo "Crosswalk mode requires VZ_QC_ADDITIONAL_POSTQC_CLUSTERS." >&2
        exit 2
      fi
    fi
    remove_clusters="${VZ_QC_REMOVE_CLUSTERS}"
    echo "Applying reviewed VZ QC decision: remove cluster(s) ${remove_clusters}."
    driver_options+=("--remove-clusters=${remove_clusters}")
    if [[ -n "${VZ_QC_ADDITIONAL_POSTQC_CLUSTERS:-}" ]]; then
      driver_options+=("--additional-postqc-clusters=${VZ_QC_ADDITIONAL_POSTQC_CLUSTERS}")
    fi
    ;;
  *)
    echo "VZ_QC_MODE must be 'review', 'apply', or 'crosswalk'." >&2
    exit 2
    ;;
esac

Rscript --max-ppsize=500000 \
  scripts/xenium_vz_04_review_qc.R "${driver_options[@]}"
