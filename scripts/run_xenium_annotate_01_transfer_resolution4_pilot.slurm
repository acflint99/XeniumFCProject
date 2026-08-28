#!/bin/bash
# Label-transfer pilot for the three resolution-4.0 whole-tissue samples.
# Submit once for each REFERENCE: Aldinger, Sepp, and Science.

#SBATCH --job-name=Xen_res4_ABT
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/%x_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/%x_%A_%a.err
#SBATCH --array=1-3%3
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --partition=short
#SBATCH --constraint=intel
#SBATCH --export=ALL

module purge
module load rc-base
module load GCC/13.2.0
module load R/4.4.1-foss-2023b

cd /home/acflint/R/Projects/XeniumFCProject || exit 1

export RENV_CONFIG_SANDBOX_ENABLED=FALSE
export RENV_CONFIG_NAMESPACES_CHECK=FALSE
export RENV_CONFIG_SYNCHRONIZED_CHECK=FALSE
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

case "${REFERENCE:-}" in
  Aldinger|Sepp|Science) ;;
  *)
    echo "REFERENCE must be Aldinger, Sepp, or Science." >&2
    exit 2
    ;;
esac

driver_options=(--pilot-res4)
if [[ "${PILOT_RES4_ABT_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

Rscript --max-ppsize=500000 scripts/xenium_annotate_01_label_transfer_rpca.R \
  "${driver_options[@]}" "${REFERENCE}" "${SLURM_ARRAY_TASK_ID}"
