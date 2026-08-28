#!/bin/bash
# Audit relative QC of Unknown consensus clusters in the three resolution-3
# weighted-2-of-3 pilot objects. No Seurat objects are modified.

#SBATCH --job-name=Xen_res3_unknownQC
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_res3_unknownQC_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_res3_unknownQC_%A_%a.err
#SBATCH --array=1-3%3
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=short
#SBATCH --constraint=intel

module purge
module load rc-base
module load R/4.4.1-foss-2023b

cd /home/acflint/R/Projects/XeniumFCProject || exit 1

driver_options=(--pilot-res3 --weighted-2of3)
if [[ "${UNKNOWN_QC_AUDIT_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

Rscript scripts/xenium_annotate_04_audit_unknown_qc.R \
  "${driver_options[@]}" "${SLURM_ARRAY_TASK_ID}"
