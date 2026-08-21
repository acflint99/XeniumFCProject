#!/bin/bash
#SBATCH --job-name=Xen_res1.5_ABT
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/xen_ABT_Sci3_%a.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/xen_ABT_Sci3_%a.err
#SBATCH --array=1-3                 
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8              
#SBATCH --mem=128G                  
#SBATCH --partition=short
#SBATCH --export=ALL

# 1. Load Modules
module purge
module load rc-base
module load GCC/13.2.0    
module load R/4.4.1-foss-2023b         

# 2. Move into project directory
cd ~/R/Projects/XeniumFCProject

# 3. Disable renv sandbox & silence warnings for parallel workers
export RENV_CONFIG_SANDBOX_ENABLED=FALSE
export RENV_CONFIG_NAMESPACES_CHECK=FALSE
export RENV_CONFIG_SYNCHRONIZED_CHECK=FALSE

# 4. OPTIMIZATION: Prevent C++ OpenMP from thread-thrashing your Slurm allocation
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# --------------------------
# Run your R script using the project library
# --------------------------
Rscript --max-ppsize=500000 scripts/AnchorBasedTransfer_RPCA_res1.5.R $SLURM_ARRAY_TASK_ID
