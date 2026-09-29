#!/bin/bash
# Combine the 18 completed k.anchor sensitivity cluster-score tables.

#SBATCH --job-name=Xen_res5_kanchor_combine
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/%x_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/%x_%j.err
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=express
#SBATCH --constraint=intel
#SBATCH --export=ALL

module purge
module load rc-base
module load R/4.4.1-foss-2023b

cd /home/acflint/R/Projects/XeniumFCProject || exit 1

export RENV_CONFIG_SANDBOX_ENABLED=FALSE
export RENV_CONFIG_NAMESPACES_CHECK=FALSE
export RENV_CONFIG_SYNCHRONIZED_CHECK=FALSE

Rscript scripts/xenium_annotate_01a_k_anchor_sensitivity.R \
  --combine --dry-run || exit 2

driver_options=(--combine)
if [[ "${K_ANCHOR_SENSITIVITY_COMBINE_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

Rscript scripts/xenium_annotate_01a_k_anchor_sensitivity.R \
  "${driver_options[@]}"
