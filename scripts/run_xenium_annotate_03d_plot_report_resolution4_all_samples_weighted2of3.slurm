#!/bin/bash
# Render PCW-ordered resolution-4 weighted 2-of-3 report pages in parallel.

#SBATCH --job-name=Xen_res4_w2of3_pages
#SBATCH --array=1-34%6
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_res4_w2of3_pages_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_res4_w2of3_pages_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --partition=express
#SBATCH --constraint=intel
#SBATCH --export=ALL

module purge
module load rc-base
module load R/4.4.1-foss-2023b

cd /home/acflint/R/Projects/XeniumFCProject || exit 1

report_options=(--all-samples-res4 --weighted-2of3 --render-page)
if [[ "${RES4_W2OF3_REPORT_OVERWRITE:-false}" == "true" ]]; then
  report_options+=(--overwrite)
fi

Rscript scripts/xenium_annotate_03d_plot_report.R "${report_options[@]}"
