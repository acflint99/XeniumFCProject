#!/bin/bash
# Build 34 weighted 2-of-3 consensus tables after all resolution-5 transfers.

#SBATCH --job-name=Xen_res5_build_w2of3
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_res5_build_w2of3_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_res5_build_w2of3_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=express
#SBATCH --constraint=intel
#SBATCH --export=ALL

module purge
module load rc-base
module load R/4.4.1-foss-2023b

cd /home/acflint/R/Projects/XeniumFCProject || exit 1

Rscript scripts/xenium_annotate_02_build_consensus.R \
  --all-samples-res5 --weighted-2of3 --dry-run || exit 2

driver_options=(--all-samples-res5 --weighted-2of3)
if [[ "${RES5_ALL_CONSENSUS_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

Rscript scripts/xenium_annotate_02_build_consensus.R "${driver_options[@]}"
