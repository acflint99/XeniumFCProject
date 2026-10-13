#!/bin/bash
#SBATCH --job-name=Xen_VZ_Harmony_audit
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/%x_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/%x_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --partition=express
#SBATCH --account=acflint
#SBATCH --qos=normal
#SBATCH --constraint=intel
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND
module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject || exit 1
export RENV_CONFIG_SANDBOX_ENABLED=FALSE
driver_options=()
if [[ "${VZ_HARMONY_AUDIT_OVERWRITE:-false}" == "true" ]]; then driver_options+=(--overwrite); fi
Rscript scripts/xenium_vz_03a_audit_harmony.R --dry-run || exit 2
Rscript scripts/xenium_vz_03a_audit_harmony.R "${driver_options[@]}"
