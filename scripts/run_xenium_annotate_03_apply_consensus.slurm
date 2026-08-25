#!/bin/bash
#SBATCH --job-name=Xen_Consensus
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_Consensus_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_Consensus_%A_%a.err
#SBATCH --array=1-34%3
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --partition=short
#SBATCH --constraint=intel

module purge
module load rc-base
module load R/4.4.1-foss-2023b

cd /home/acflint/R/Projects/XeniumFCProject || exit 1

driver_options=()
if [[ "${CONSENSUS_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi
if [[ "${CONSENSUS_DOTPLOT_ONLY:-false}" == "true" ]]; then
  driver_options+=(--dotplot-only)
fi

Rscript scripts/xenium_annotate_03_apply_consensus.R \
  "${driver_options[@]}" "${SLURM_ARRAY_TASK_ID}"
