#!/bin/bash
#SBATCH --job-name=Science_R_job
#SBATCH --output=/home/acflint/R/Projects/XeniumFCProject/logs/Science_R_job_%j.out
#SBATCH --error=/home/acflint/R/Projects/XeniumFCProject/logs/Science_R_job_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=express

# Load R module first
module purge
module load rc-base
module load R/4.4.1-foss-2023b

# Move into the project directory (CRITICAL)
cd /home/acflint/R/Projects/XeniumFCProject || exit 1

# The project .Rprofile activates renv. Run renv::restore() manually before
# submission when the lockfile changes; never mutate the library inside a job.

# Inspect paths without loading either RDS before the protected rerun.
Rscript scripts/science_02_subset_gene_panel.R --dry-run || exit 2

driver_options=()
if [[ "${SCIENCE_REFERENCE_OVERWRITE:-false}" == "true" ]]; then
  driver_options+=(--overwrite)
fi

# Run the analysis using the existing project library.
Rscript scripts/science_02_subset_gene_panel.R "${driver_options[@]}"
