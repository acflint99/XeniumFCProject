#!/bin/bash
# Apply weighted 2-of-3 consensus and create final plots for one selected
# resolution-5 sample.

#SBATCH --job-name=Xen_res5_selected_w2of3
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_res5_selected_w2of3_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_res5_selected_w2of3_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=express
#SBATCH --constraint=intel
#SBATCH --export=ALL

module purge
module load rc-base
module load R/4.4.1-foss-2023b

cd /home/acflint/R/Projects/XeniumFCProject || exit 1

if [[ -z "${SELECTED_SAMPLE_ID:-}" ]]; then
  echo "SELECTED_SAMPLE_ID is required." >&2
  exit 2
fi

driver_options=(--pilot-res5 --selected-sample --weighted-2of3)
if [[ "${RES5_SELECTED_W2OF3_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

Rscript scripts/xenium_annotate_03_apply_consensus.R \
  "${driver_options[@]}" 1
