#!/bin/bash
# Regenerate all final consensus plots from existing resolution-4 weighted
# 2-of-3 consensus RDS objects. No RDS or consensus table is rewritten.

#SBATCH --job-name=Xen_res4_w2of3_replot
#SBATCH --array=1-34%4
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_res4_w2of3_replot_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_res4_w2of3_replot_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=express
#SBATCH --constraint=intel
#SBATCH --export=ALL

module purge
module load rc-base
module load R/4.4.1-foss-2023b

cd /home/acflint/R/Projects/XeniumFCProject || exit 1

Rscript scripts/xenium_annotate_03_apply_consensus.R \
  --all-samples-res4 --weighted-2of3 --plots-only "${SLURM_ARRAY_TASK_ID}"
