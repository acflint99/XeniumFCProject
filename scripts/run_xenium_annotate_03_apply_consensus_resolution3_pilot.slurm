#!/bin/bash
# Apply resolution-3.0 consensus labels to the three pilot samples.

#SBATCH --job-name=Xen_res3_consensus
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_res3_consensus_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_res3_consensus_%A_%a.err
#SBATCH --array=1-3%3
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --partition=short
#SBATCH --constraint=intel

module purge
module load rc-base
module load R/4.4.1-foss-2023b

cd /home/acflint/R/Projects/XeniumFCProject || exit 1

driver_options=(--pilot-res3)
if [[ "${PILOT_RES3_CONSENSUS_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi
if [[ "${PILOT_RES3_CONSENSUS_DOTPLOT_ONLY:-false}" == "true" ]]; then
  driver_options+=(--dotplot-only)
fi

Rscript scripts/xenium_annotate_03_apply_consensus.R \
  "${driver_options[@]}" "${SLURM_ARRAY_TASK_ID}"
