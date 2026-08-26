#!/bin/bash
#SBATCH --job-name=Xen_RL_Counts
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_RL_Counts_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_RL_Counts_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --partition=short
#SBATCH --constraint=intel

module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject || exit 1

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

driver_options=()
if [[ "${RL_COUNT_PLOTS_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

Rscript scripts/xenium_rl_09_plot_counts.R "${driver_options[@]}"
