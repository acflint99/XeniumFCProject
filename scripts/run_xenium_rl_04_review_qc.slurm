#!/bin/bash
#SBATCH --job-name=Xen_RL_QC
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_RL_QC_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_RL_QC_%j.err
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
if [[ "${RL_QC_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

# Reviewed 2026-10-02 from the resolution-0.8 cluster QC, sample-composition,
# and top-50 marker evidence. Review mode remains the launcher default; this
# decision is used only when RL_QC_MODE=apply is explicitly requested.
reviewed_removal_clusters="6,15,24,26,27"

case "${RL_QC_MODE:-review}" in
  review)
    driver_options+=(--qc-only)
    ;;
  apply)
    selected_removal_clusters="${RL_QC_REMOVE_CLUSTERS:-${reviewed_removal_clusters}}"
    echo "Applying reviewed RL removal clusters: ${selected_removal_clusters}"
    driver_options+=("--remove-clusters=${selected_removal_clusters}")
    ;;
  *)
    echo "RL_QC_MODE must be 'review' or 'apply'." >&2
    exit 2
    ;;
esac

Rscript --max-ppsize=500000 \
  scripts/xenium_rl_04_review_qc.R "${driver_options[@]}"
