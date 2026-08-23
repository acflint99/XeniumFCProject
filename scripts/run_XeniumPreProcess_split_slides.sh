#!/bin/bash
# Submit all 15 manually separated biological samples. At most three heavy
# Xenium jobs run concurrently to limit memory and shared-filesystem pressure.

#SBATCH --job-name=Xen_split_preprocess
#SBATCH --array=1-15%3
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/split_preprocess_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/split_preprocess_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --partition=medium
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

driver_options=()
if [[ "${SPLIT_PREPROCESS_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

Rscript scripts/XeniumPreProcess_split_slides.R \
  "${driver_options[@]}" "${SLURM_ARRAY_TASK_ID}"
