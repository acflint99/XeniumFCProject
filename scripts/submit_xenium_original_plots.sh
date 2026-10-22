#!/bin/bash
# In-place figure redraw only. No annotation, integration, or RDS-writing jobs.
set -euo pipefail
cd /home/acflint/R/Projects/XeniumFCProject
mode=dry-run
selected=consensus,vz,rl,residual,combined,vz-pooled,rl-pooled,residual-pooled,counts-vz,counts-rl
reports=true
for option in "$@"; do
  case "$option" in
    --dry-run) mode=dry-run ;;
    --submit) mode=submit ;;
    --groups=*) selected="${option#--groups=}" ;;
    --no-reports) reports=false ;;
    *) echo "Use [--dry-run|--submit] [--groups=STAGE,...] [--no-reports]." >&2; exit 2 ;;
  esac
done
if [[ ! "$selected" =~ ^[a-z]+(-[a-z]+)?(,[a-z]+(-[a-z]+)?)*$ ]]; then echo "Invalid group list." >&2; exit 2; fi
IFS=',' read -r -a groups <<< "$selected"
seen=,
for stage in "${groups[@]}"; do
  case "$stage" in consensus|vz|rl|residual|combined|vz-pooled|rl-pooled|residual-pooled|counts-vz|counts-rl|nonspatial|spatial) ;;
    *) echo "Unknown native plot stage: ${stage}" >&2; exit 2 ;;
  esac
  if [[ "$seen" == *",${stage},"* ]]; then echo "Repeated stage: ${stage}" >&2; exit 2; fi
  seen+="${stage},"
done

# Inspect every selected input set before submitting anything. A path preflight
# cannot certify that saved biological annotations are current.
failed=false
for stage in "${groups[@]}"; do
  if ! Rscript --max-ppsize=500000 scripts/xenium_original_plot_plan.R "--stage=${stage}" --dry-run; then failed=true; fi
done
if [[ "$failed" == true ]]; then echo "Preflight failed; no jobs submitted." >&2; exit 2; fi

# Resolve every bound before the first sbatch. Ignore renv startup messages.
counts=()
for stage in "${groups[@]}"; do
  mapping="$(Rscript --max-ppsize=500000 scripts/xenium_original_plot_plan.R "--stage=${stage}" --list)"
  if ! count="$(awk -F '\t' '$1~/^[1-9][0-9]*$/ && NF>=2 {if(seen[$1]++)bad=1;n++;if($1>max)max=$1} END {if(bad||n<1||n!=max)exit 2;print n}' <<< "$mapping")"; then
    echo "Invalid/noncontiguous task mapping: ${stage}" >&2; exit 2
  fi
  counts+=("$count")
done
submit_job() {
  local job_stage="$1" array="$2" memory="$3" limit="$4" partition="$5" dependency="${6:-}"
  local result
  if [[ -n "$dependency" ]]; then
    result="$(sbatch --parsable --dependency="afterok:${dependency}" --array="$array" --mem="$memory" --time="$limit" --partition="$partition" \
      --export="ALL,XENIUM_ORIGINAL_STAGE=${job_stage},XENIUM_ORIGINAL_OVERWRITE=true" scripts/run_xenium_original_plots.slurm)"
  else
    result="$(sbatch --parsable --array="$array" --mem="$memory" --time="$limit" --partition="$partition" \
      --export="ALL,XENIUM_ORIGINAL_STAGE=${job_stage},XENIUM_ORIGINAL_OVERWRITE=true" scripts/run_xenium_original_plots.slurm)"
  fi
  result="${result%%;*}"
  if [[ ! "$result" =~ ^[0-9]+$ ]]; then echo "Unrecognized sbatch job ID: ${result}" >&2; return 2; fi
  printf '%s\n' "$result"
}
for index in "${!groups[@]}"; do
  stage="${groups[$index]}";count="${counts[$index]}"
  concurrency=3;memory=64G;limit=02:00:00;partition=express
  case "$stage" in
    residual-pooled) concurrency=1;memory=32G;limit=01:00:00;partition=express ;;
    *-pooled) concurrency=1;memory=256G;limit=12:00:00;partition=short ;;
    counts-*) concurrency=1;memory=128G ;;
    nonspatial|spatial) concurrency=1;memory=600G;limit=24:00:00;partition=medium ;;
  esac
  echo "${stage}: ${count} tasks, concurrency ${concurrency}, ${memory}, limit ${limit}; original plot paths"
  if [[ "$mode" == submit ]]; then
    job="$(submit_job "$stage" "1-${count}%${concurrency}" "$memory" "$limit" "$partition")"
    echo "Submitted ${stage}: ${job}"
    if [[ "$reports" == true ]]; then
      case "$stage" in
        vz|rl|residual)
          report="$(submit_job "${stage}-report" '1-1%1' 64G 02:00:00 express "$job")"
          echo "Submitted ${stage} combined PDF after ${job}: ${report}" ;;
        combined)
          report="$(submit_job combined-report '1-1%1' 64G 02:00:00 express "$job")"
          report085="$(submit_job combined-report085 '1-1%1' 64G 02:00:00 express "$job")"
          echo "Submitted combined/corrected-broad PDFs at both point sizes after ${job}: ${report}, ${report085}" ;;
        consensus)
          pages="$(submit_job consensus-pages "1-${count}%3" 64G 02:00:00 express "$job")"
          report="$(submit_job consensus-report '1-1%1' 64G 02:00:00 express "$pages")"
          echo "Submitted consensus pages after ${job}: ${pages}; merged reports after ${pages}: ${report}" ;;
      esac
    fi
  fi
done
if [[ "$mode" == dry-run ]]; then echo "No jobs submitted. --submit explicitly replaces original plot files; no snapshot directory is used."; fi
