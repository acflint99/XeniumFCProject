#!/bin/bash
# Combine the 18 completed reference-sampling/abstention sensitivity tasks into
# reference-vote, equal 2-of-3 consensus, reviewed-cluster, and Unknown summaries.

#SBATCH --job-name=Xen_res5_sampling_abstain_combine
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

Rscript scripts/xenium_annotate_01c_reference_sampling_abstention_sensitivity.R \
  --combine --dry-run || exit 2

driver_options=(--combine)
if [[ "${REFERENCE_SAMPLING_ABSTENTION_COMBINE_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

Rscript scripts/xenium_annotate_01c_reference_sampling_abstention_sensitivity.R \
  "${driver_options[@]}"
