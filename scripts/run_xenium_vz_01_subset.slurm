#!/bin/bash
#SBATCH --job-name=Xenium_VZSubsets
#SBATCH --array=1-41%3
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/VZsubset_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/VZsubset_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=express
#SBATCH --constraint=intel
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND

# 1. Load Modules
module purge
module load rc-base
module load R/4.4.1-foss-2023b

# 2. Move into project directory (CRITICAL for here() and renv)
cd /home/acflint/R/Projects/XeniumFCProject || exit 1

# 3. Disable renv sandbox (Smart move—prevents slow symlinking on HPC)
export RENV_CONFIG_SANDBOX_ENABLED=FALSE

driver_options=()
if [[ "${VZ_SUBSET_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

Rscript --max-ppsize=500000 scripts/xenium_vz_01_subset_samples.R \
  "${driver_options[@]}" "${SLURM_ARRAY_TASK_ID}"
