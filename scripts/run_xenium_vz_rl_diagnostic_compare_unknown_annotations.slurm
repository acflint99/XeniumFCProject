#!/bin/bash
#SBATCH --job-name=Xen_Unknown_VZRL
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_Unknown_VZRL_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_Unknown_VZRL_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --partition=short
#SBATCH --constraint=intel
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND

module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject || exit 1

driver_options=()
if [[ "${UNKNOWN_COMPARISON_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi
Rscript --max-ppsize=500000 scripts/xenium_vz_rl_diagnostic_compare_unknown_annotations.R "${driver_options[@]}"
