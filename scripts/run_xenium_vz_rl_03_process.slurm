#!/bin/bash
#SBATCH --job-name=Xen_CombProcess
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_CombProcess_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_CombProcess_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --partition=medium
#SBATCH --constraint=intel

module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject || exit 1

# One R process is used. Permit its compiled PCA/Harmony/UMAP routines to use
# the CPUs allocated by Slurm without spawning full copies of the Seurat object.
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

driver_options=(--process-only)
if [[ "${COMBINED_PROCESS_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi
Rscript 'scripts/xenium_vz_rl_03_process_and_plot.R' "${driver_options[@]}"
