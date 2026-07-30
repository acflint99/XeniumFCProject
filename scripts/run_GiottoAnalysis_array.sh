#!/bin/bash
#SBATCH --job-name=Giotto_Array
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Giotto_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Giotto_%A_%a.err
#SBATCH --time=12:00:00        # Dropped to 2 hours since each job only does 1 sample!
#SBATCH --cpus-per-task=8
#SBATCH --mem=400G
#SBATCH --partition=medium
#SBATCH --array=1-2           # Launch 16 identical jobs

# Load R module
module purge
module load rc-base
module load R/4.4.1-foss-2023b

cd ~/R/Projects/XeniumFCProject
export R_LIBS_USER=~/R/Projects/XeniumFCProject/renv/library/linux-rhel-7.9/R-4.4/x86_64-pc-linux-gnu

# Define the remaining 15 samples in a bash array
SAMPLES=(
  #"GZFB_9_X_G_3" "FB328_1_X_G" 
  #"FB330_1_X_G" 
  #"FB78_X_G" "GZFB5_X_G" "GZFB_12_X_G_1" 
  #"GZFB_12_X_G_2" "GZFB_12_X_G_3" "GZFB_12_X_G_4" "GZFB_12_X_G_5" 
  #"GZFB_9_X_G_1" "GZFB_9_X_G_2" 
  #"FB124_X_G" 
  "FB198_X_G" "GZFB4_X_G"
   #"GZFB_1_X_G"
)

# Bash arrays are 0-indexed, but Slurm arrays start at 1. We subtract 1 to get the right sample.
CURRENT_SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID - 1]}

echo "Launching Giotto Pipeline for Sample: $CURRENT_SAMPLE"

# Pass the sample name directly to the R script
Rscript scripts/Xenium_Giotto_Analysis_array.R "$CURRENT_SAMPLE"