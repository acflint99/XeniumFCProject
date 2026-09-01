#!/bin/bash
# Merge 34 resolution-5 weighted 2-of-3 report pages into one PDF.

#SBATCH --job-name=Xen_res5_w2of3_merge
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_res5_w2of3_merge_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_res5_w2of3_merge_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=express
#SBATCH --constraint=intel
#SBATCH --export=ALL

module purge
module load rc-base
module load R/4.4.1-foss-2023b

cd /home/acflint/R/Projects/XeniumFCProject || exit 1

report_options=(--all-samples-res5 --weighted-2of3 --merge-pages)
if [[ "${RES5_W2OF3_REPORT_OVERWRITE:-false}" == "true" ]]; then
  report_options+=(--overwrite)
fi

Rscript scripts/xenium_annotate_03d_plot_report.R "${report_options[@]}"
