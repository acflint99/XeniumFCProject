#!/bin/bash
#SBATCH --job-name=Xen_Consensus_Report
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_Consensus_Report_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_Consensus_Report_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --partition=medium
#SBATCH --constraint=intel

module purge
module load rc-base
module load R/4.4.1-foss-2023b

cd /home/acflint/R/Projects/XeniumFCProject || exit 1

report_options=()
if [[ "${CONSENSUS_REPORT_OVERWRITE:-false}" == "true" ]]; then
  report_options+=(--overwrite)
fi

Rscript scripts/xenium_annotate_03d_plot_report.R "${report_options[@]}"
