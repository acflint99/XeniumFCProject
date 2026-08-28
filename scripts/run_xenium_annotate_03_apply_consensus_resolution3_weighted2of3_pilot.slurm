#!/bin/bash
# Apply the weighted-reference, 2-of-3 consensus method to the three
# completed resolution-3.0 pilot samples. Legacy consensus outputs are untouched.

#SBATCH --job-name=Xen_res3_w2of3
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_res3_w2of3_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_res3_w2of3_%A_%a.err
#SBATCH --array=1-3%3
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=short
#SBATCH --constraint=intel

module purge
module load rc-base
module load R/4.4.1-foss-2023b

cd /home/acflint/R/Projects/XeniumFCProject || exit 1

driver_options=(--pilot-res3 --weighted-2of3)
if [[ "${PILOT_RES3_W2OF3_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

Rscript scripts/xenium_annotate_03_apply_consensus.R \
  "${driver_options[@]}" "${SLURM_ARRAY_TASK_ID}"
