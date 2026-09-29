#!/bin/bash
# Validate one complete developmental age/stage field in each reference before
# any stratified label-transfer task is released.

#SBATCH --job-name=Xen_res5_strata_preflight
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/%x_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --partition=express
#SBATCH --account=acflint
#SBATCH --qos=normal
#SBATCH --constraint=intel
#SBATCH --export=ALL
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=FAIL

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

Rscript scripts/xenium_annotate_01d_stratified_sampling_sensitivity.R \
  --inspect-strata --dry-run || exit 2

driver_options=(--inspect-strata)
if [[ "${STRATIFIED_SAMPLING_PREFLIGHT_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

Rscript scripts/xenium_annotate_01d_stratified_sampling_sensitivity.R \
  "${driver_options[@]}"
