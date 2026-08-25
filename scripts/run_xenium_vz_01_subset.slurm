#!/bin/bash
#SBATCH --job-name=Xenium_VZSubsets
#SBATCH --array=1-34%3
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/VZsubset_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/VZsubset_%a.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G             
#SBATCH --partition=short
#SBATCH --constraint=intel

# 1. Load Modules
module purge
module load rc-base
module load R/4.4.1-foss-2023b

# 2. Move into project directory (CRITICAL for here() and renv)
cd ~/R/Projects/XeniumFCProject

# 3. Disable renv sandbox (Smart move—prevents slow symlinking on HPC)
export RENV_CONFIG_SANDBOX_ENABLED=FALSE

driver_options=()
if [[ "${VZ_SUBSET_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

Rscript scripts/xenium_vz_01_subset_samples.R \
  "${driver_options[@]}" "${SLURM_ARRAY_TASK_ID}"
