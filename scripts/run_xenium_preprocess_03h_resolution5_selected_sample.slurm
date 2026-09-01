#!/bin/bash
# Resolution-5 reclustering for one manifest-validated selected sample.

#SBATCH --job-name=Xen_res5_selected
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_res5_selected_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_res5_selected_%j.err
#SBATCH --time=01:00:00
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

driver_options=(--resolution5 --selected-sample)
if [[ "${RES5_SELECTED_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

Rscript scripts/xenium_preprocess_03b_resolution2_pilot.R \
  "${driver_options[@]}" 1
