#!/bin/bash
# Assemble the existing Aldinger resolution-5 label-transfer TIFFs into one
# 34-page PDF without recomputing the label-transfer analysis.

#SBATCH --job-name=Xen_res5_Ald_pdf
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_res5_Ald_pdf_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_res5_Ald_pdf_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=express
#SBATCH --constraint=intel
#SBATCH --export=ALL

module purge
module load rc-base
module load R/4.4.1-foss-2023b

cd /home/acflint/R/Projects/XeniumFCProject || exit 1

Rscript --max-ppsize=500000 \
  scripts/xenium_annotate_01b_build_aldinger_resolution5_plot_report.R \
  --dry-run || exit 2

report_options=()
if [[ "${ALDINGER_RES5_REPORT_OVERWRITE:-false}" == "true" ]]; then
  report_options+=(--overwrite)
fi

Rscript --max-ppsize=500000 \
  scripts/xenium_annotate_01b_build_aldinger_resolution5_plot_report.R \
  "${report_options[@]}"
