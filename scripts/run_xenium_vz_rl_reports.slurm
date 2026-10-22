#!/bin/bash
#SBATCH --job-name=Xen_VZRL_Reports
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_VZRL_Reports_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_VZRL_Reports_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --partition=express
#SBATCH --constraint=intel
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND

set -euo pipefail
module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject

driver_options=("--labels=${VZRL_REPORT_LABELS:-both}" "--point-size=${VZRL_REPORT_POINT_SIZE:-0.6}")
Rscript --max-ppsize=500000 scripts/xenium_vz_rl_reports.R "${driver_options[@]}" --dry-run
if [[ "${VZRL_REPORT_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi
Rscript --max-ppsize=500000 scripts/xenium_vz_rl_reports.R "${driver_options[@]}"
