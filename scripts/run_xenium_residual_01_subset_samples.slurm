#!/bin/bash
#SBATCH --job-name=Xen_residual_subset
#SBATCH --array=1-41%3
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/%x_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/%x_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=express
#SBATCH --account=acflint
#SBATCH --qos=normal
#SBATCH --constraint=intel
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=FAIL
module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject || exit 1
export RENV_CONFIG_SANDBOX_ENABLED=FALSE
driver_options=()
if [[ "${RESIDUAL_SUBSET_OVERWRITE:-false}" == "true" ]]; then driver_options+=(--overwrite); fi
Rscript --max-ppsize=500000 \
  scripts/xenium_residual_01_subset_samples.R \
  --dry-run "${SLURM_ARRAY_TASK_ID}" || exit 2
Rscript --max-ppsize=500000 \
  scripts/xenium_residual_01_subset_samples.R \
  "${driver_options[@]}" "${SLURM_ARRAY_TASK_ID}"
