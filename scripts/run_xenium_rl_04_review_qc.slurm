#!/bin/bash
#SBATCH --job-name=Xen_RL_QC
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_RL_QC_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_RL_QC_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
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

case "${RL_QC_MODE:-review}" in
  review)
    driver_options+=(--qc-only)
    ;;
  apply)
    if [[ -z "${RL_QC_REMOVE_CLUSTERS:-}" ]]; then
      echo "RL_QC_REMOVE_CLUSTERS is required when RL_QC_MODE=apply." >&2
      echo "Use a comma-separated cluster list or the literal value 'none'." >&2
      exit 2
    fi
    driver_options+=("--remove-clusters=${RL_QC_REMOVE_CLUSTERS}")
    ;;
  *)
    echo "RL_QC_MODE must be 'review' or 'apply'." >&2
    exit 2
    ;;
esac

Rscript scripts/xenium_rl_04_review_qc.R "${driver_options[@]}"
