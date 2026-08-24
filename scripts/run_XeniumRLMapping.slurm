#!/bin/bash
#SBATCH --job-name=Xen_RL_Map
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_RL_Map_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_RL_Map_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --partition=short
#SBATCH --constraint=intel

module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject || exit 1

driver_options=()
if [[ "${RL_MAPPING_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi
Rscript scripts/Xenium_RL_Mapping_res1.5.R "${driver_options[@]}"
