#!/bin/bash
# Combine the 27 completed intermediate-cap tasks with the existing cap1000 and
# full-reference results into five-mode convergence and abstention summaries.

#SBATCH --job-name=Xen_res5_cap_curve_combine
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/%x_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/%x_%j.err
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=express
#SBATCH --account=acflint
#SBATCH --qos=normal
#SBATCH --constraint=intel
#SBATCH --export=ALL
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND

module purge
module load rc-base
module load R/4.4.1-foss-2023b

cd /home/acflint/R/Projects/XeniumFCProject || exit 1

export RENV_CONFIG_SANDBOX_ENABLED=FALSE
export RENV_CONFIG_NAMESPACES_CHECK=FALSE
export RENV_CONFIG_SYNCHRONIZED_CHECK=FALSE

Rscript scripts/xenium_annotate_01e_reference_cap_convergence.R \
  --combine --dry-run || exit 2

driver_options=(--combine)
if [[ "${REFERENCE_CAP_CONVERGENCE_COMBINE_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

Rscript scripts/xenium_annotate_01e_reference_cap_convergence.R \
  "${driver_options[@]}"
