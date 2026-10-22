#!/bin/bash
#SBATCH --job-name=Xen_GeneMatrix
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_GeneMatrix_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_GeneMatrix_%j.err
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --partition=express
#SBATCH --constraint=intel
#SBATCH --mail-user=acflint@uab.edu
#SBATCH --mail-type=END,FAIL,INVALID_DEPEND

# Expected runtime: 5-30 minutes, low confidence; excludes queue wait.
# Provisional 128G memory: 24G merged RDS plus deserialization and count-layer
# temporaries. No completed export MaxRSS measurement is available yet.
# Collect sacct after completion before refining this request.
set -euo pipefail
module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject || exit 1

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# No overwrite flag: existing or partial outputs require review.
Rscript --max-ppsize=500000 "scripts/Gene expression matrix.R" --dry-run
Rscript --max-ppsize=500000 "scripts/Gene expression matrix.R"
