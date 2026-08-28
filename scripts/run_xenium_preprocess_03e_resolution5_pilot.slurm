#!/bin/bash
# Resolution-5.0 whole-tissue pilot for the same three manifest samples.
# Reclusters the saved resolution-1.5 SNN graph and writes only beneath
# outputs/xenium/preprocess/03e_resolution5_pilot/.

#SBATCH --job-name=Xen_res5_pilot
#SBATCH --array=1-3%3
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_res5_pilot_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_res5_pilot_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=short
#SBATCH --constraint=intel
#SBATCH --export=ALL

module purge
module load rc-base
module load R/4.4.1-foss-2023b

cd /home/acflint/R/Projects/XeniumFCProject || exit 1

export RENV_CONFIG_SANDBOX_ENABLED=FALSE
export RENV_CONFIG_NAMESPACES_CHECK=FALSE
export RENV_CONFIG_SYNCHRONIZED_CHECK=FALSE
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

driver_options=(--resolution5)
if [[ "${RES5_PILOT_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

Rscript scripts/xenium_preprocess_03b_resolution2_pilot.R \
  "${driver_options[@]}" "${SLURM_ARRAY_TASK_ID}"
