#!/bin/bash
#SBATCH --job-name=Xen_unknown_evidence
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_unknown_evidence_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_unknown_evidence_%j.err
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=express
#SBATCH --constraint=intel
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND

module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject || exit 1
export RENV_CONFIG_SANDBOX_ENABLED=FALSE

driver_options=()
if [[ "${UNKNOWN_EVIDENCE_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi
Rscript scripts/xenium_unknown_04a_compile_reference_evidence.R \
  "${driver_options[@]}"
