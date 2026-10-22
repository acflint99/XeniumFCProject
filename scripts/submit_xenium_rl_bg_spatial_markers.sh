#!/bin/bash
# Lists and preflights by default. --submit explicitly submits the reviewed array.
set -euo pipefail
if [[ $# -gt 1 || ( $# -eq 1 && "$1" != "--submit" ) ]]; then
  echo "Usage: bash scripts/submit_xenium_rl_bg_spatial_markers.sh [--submit]" >&2
  exit 2
fi
module purge
module load rc-base
module load R/4.4.1-foss-2023b
cd /home/acflint/R/Projects/XeniumFCProject || exit 1
driver=scripts/xenium_rl_06d_bg_spatial_markers.R
task_map=$(Rscript --max-ppsize=500000 "$driver" --list)
echo "$task_map"
task_count=$(echo "$task_map" | awk 'NR>1 {n++} END {print n+0}')
if [[ "$task_count" -lt 1 ]]; then echo "Empty manifest task mapping." >&2; exit 2; fi
Rscript --max-ppsize=500000 "$driver" --dry-run
echo "Array: 1-${task_count}%3; expected 5-30 min/task, about 1-7 h for 41 samples at concurrency 3, excluding queue wait."
echo "Requested: 2 hours/task, 64 GB/task, express; low-confidence estimate pending first-run sacct."
if [[ "${1:-}" == "--submit" ]]; then
  sbatch --array="1-${task_count}%3" scripts/run_xenium_rl_06d_bg_spatial_markers.slurm
fi
