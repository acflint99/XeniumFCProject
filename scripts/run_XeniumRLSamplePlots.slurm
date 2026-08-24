#!/bin/bash
#SBATCH --job-name=Xen_RL_Plots
#SBATCH --array=1-34%3
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_RL_Plots_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_RL_Plots_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --partition=short
#SBATCH --constraint=intel

module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject || exit 1

driver_options=()
if [[ "${RL_PLOTS_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi
Rscript scripts/Xenium_RL_SamplePlots_Loop_res1.5.R \
  "${driver_options[@]}" "${SLURM_ARRAY_TASK_ID}"
