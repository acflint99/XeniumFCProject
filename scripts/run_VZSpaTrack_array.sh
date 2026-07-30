#!/bin/bash
#SBATCH --job-name=SpaTrack_Array
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/spatrack_%A_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/spatrack_%A_%a.err
#SBATCH --array=1
#SBATCH --time=12:00:00 #10 for most samples
#SBATCH --cpus-per-task=4
#SBATCH --mem=500G #200G for most samples
#SBATCH --partition=medium

# Load R module
module purge
module load rc-base
module load R/4.4.1-foss-2023b

cd ~/R/Projects/XeniumFCProject
export R_LIBS_USER=~/R/Projects/XeniumFCProject/renv/library/linux-rhel-7.9/R-4.4/x86_64-pc-linux-gnu

# NO COMMAS in the list below
SAMPLES=(
  "GZFB4_X_G" #"FB124_X_G" 
  #"FB198_X_G" "FB328_1_X_G"
  #"FB330_1_X_G" "FB78_X_G" "GZFB5_X_G" "GZFB_12_X_G_1"
  #"GZFB_12_X_G_2" "GZFB_12_X_G_3" "GZFB_12_X_G_4" "GZFB_12_X_G_5"
  #"GZFB_1_X_G" "GZFB_9_X_G_1" "GZFB_9_X_G_2" "GZFB_9_X_G_3"
)

INDEX=$(($SLURM_ARRAY_TASK_ID - 1))
CURRENT_SAMPLE=${SAMPLES[$INDEX]}

echo "Launching VZ SpaTrack Pipeline for Sample: $CURRENT_SAMPLE"

Rscript scripts/Xenium_VZRL_Subclusters_VZ_SpaTrack.R $CURRENT_SAMPLE