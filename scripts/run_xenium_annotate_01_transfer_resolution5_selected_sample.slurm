#!/bin/bash
# One reference transfer for a selected resolution-5 sample. Submit once per
# REFERENCE: Aldinger, Sepp, and Science.

#SBATCH --job-name=Xen_res5_selected_ABT
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/%x_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/%x_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --partition=short
#SBATCH --constraint=intel
#SBATCH --export=ALL

module purge
module load rc-base
module load GCC/13.2.0
module load R/4.4.1-foss-2023b

cd /home/acflint/R/Projects/XeniumFCProject || exit 1

if [[ -z "${SELECTED_SAMPLE_ID:-}" ]]; then
  echo "SELECTED_SAMPLE_ID is required." >&2
  exit 2
fi
case "${REFERENCE:-}" in
  Aldinger|Sepp|Science) ;;
  *)
    echo "REFERENCE must be Aldinger, Sepp, or Science." >&2
    exit 2
    ;;
esac

driver_options=(--pilot-res5 --selected-sample)
if [[ "${RES5_SELECTED_ABT_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

Rscript --max-ppsize=500000 scripts/xenium_annotate_01_label_transfer_rpca.R \
  "${driver_options[@]}" "${REFERENCE}" 1
