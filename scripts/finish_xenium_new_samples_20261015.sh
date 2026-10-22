#!/bin/bash
# Apply the user-approved keep decisions and queue only the new regional splits.
# --dry-run: inspect paths/CSV files only. --submit: back up and refresh shared
# tables/reports, then submit dependency-linked jobs. No consensus RDS replacement.
set -euo pipefail
if [[ $# -ne 1 ]]; then
  echo "Usage: bash scripts/finish_xenium_new_samples_20261015.sh [--dry-run|--submit]" >&2
  exit 2
fi
mode=$1
case "$mode" in
  --dry-run|--submit|--combine-qc|--verify) ;;
  *) echo "Unknown mode: $mode" >&2; exit 2 ;;
esac
cd /home/acflint/R/Projects/XeniumFCProject
cpu_vendor=$(lscpu | awk -F: '/Vendor ID/ {gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2); print $2; exit}')
[[ $cpu_vendor == GenuineIntel ]] || {
  echo "Use an Intel session for the project R library." >&2; exit 2;
}
module purge
module load rc-base
module load R/4.4.1-foss-2023b
export RENV_CONFIG_SANDBOX_ENABLED=FALSE
export RENV_CONFIG_NAMESPACES_CHECK=FALSE
export RENV_CONFIG_SYNCHRONIZED_CHECK=FALSE

# Several legacy dry-runs print FAIL but return zero: check their text as well.
run_dry() {
  local output
  if ! output=$(Rscript --max-ppsize=500000 "$@" 2>&1); then
    printf '%s\n' "$output" >&2; exit 2
  fi
  printf '%s\n' "$output"
  if [[ "$output" == *"DRY-RUN FAIL"* || "$output" != *"DRY-RUN PASS"* ]]; then
    echo "Dry-run did not pass; stopping." >&2; exit 2
  fi
}

if [[ $mode == --combine-qc ]]; then
  run_dry scripts/xenium_annotate_04b_apply_unknown_qc.R --combine --dry-run
  Rscript --max-ppsize=500000 scripts/xenium_annotate_04b_apply_unknown_qc.R --combine --overwrite
  exit 0
fi

if [[ $mode == --verify ]]; then
  Rscript --max-ppsize=500000 -e '
source("scripts/R/config.R")
cfg <- load_pipeline_config()
ids <- readLines(resolve_config_path("config/new_samples_20261014.txt",cfg))
stopifnot(length(ids)==7L,!anyDuplicated(ids))
paths <- unlist(lapply(c("vz","rl"),function(branch) {
  file.path(here::here(cfg$project$outputs_dir),"xenium",branch,"01_subsets","rds",
            paste0(ids,"_",toupper(branch),"subset.rds"))
}))
stopifnot(length(paths)==14L,all(file.exists(paths)),all(file.info(paths)$size>0))
cat("PASS: all 14 new VZ/RL subset files exist and are nonempty.\n")
cat(paste(paths,collapse="\n"),"\n")
'
  exit 0
fi

[[ -f reviewed_unknown_qc_APPROVED.csv ]] || {
  echo "Upload reviewed_unknown_qc_APPROVED.csv to the project root first." >&2; exit 2;
}

# Use the HPC review as authoritative. Import only the approved cohort additions,
# preserving existing rows/notes/extra columns rather than replacing from a local copy.
# This same code is read-only during preflight and writes only after backup.
qc_review_code=$(cat <<'RSCRIPT'
source("scripts/R/config.R")
source("scripts/R/consensus_labels.R")
source("scripts/R/xenium_unknown_qc_review.R")
source("scripts/R/xenium_broad_reference.R")
cfg <- load_pipeline_config()
samples <- load_sample_manifest(cfg)
ids <- readLines(resolve_config_path("config/new_samples_20261014.txt",cfg))
stopifnot(length(ids)==7L,!anyDuplicated(ids),all(ids %in% samples$sample_id))
old <- load_reviewed_unknown_qc(cfg)
uploaded <- read.csv("reviewed_unknown_qc_APPROVED.csv",check.names=FALSE,stringsAsFactors=FALSE,na.strings=character())
validate_unknown_qc_review(old,samples$sample_id,require_complete=TRUE)
validate_unknown_qc_review(uploaded,samples$sample_id,require_complete=TRUE)
required <- c("sample_id","consensus_label","decision","review_notes")
additions <- uploaded[uploaded$sample_id %in% ids,required,drop=FALSE]
if(nrow(additions)!=19L || any(additions$decision!="keep")) {
  stop("Expected exactly 19 approved keep decisions for the seven new samples.")
}
old_keys <- unknown_qc_review_key(old)
addition_keys <- unknown_qc_review_key(additions)
unexpected <- setdiff(old_keys[old$sample_id %in% ids],addition_keys)
if(length(unexpected)) stop("Existing new-sample QC keys are absent from the approved additions: ",paste(unexpected,collapse=", "))
hit <- match(addition_keys,old_keys)
existing <- !is.na(hit)
conflict <- existing
conflict[existing] <- tolower(trimws(old$decision[hit[existing]]))!="keep"
if(any(conflict)) stop("Existing HPC decisions conflict with approved keep additions: ",paste(addition_keys[conflict],collapse=", "))
missing <- additions[!existing,,drop=FALSE]
new <- old
if(nrow(missing)) {
  extra_columns <- setdiff(names(old),names(missing))
  for(column in extra_columns) missing[[column]] <- old[[column]][rep(NA_integer_,nrow(missing))]
  missing <- missing[,names(old),drop=FALSE]
  new <- rbind(old,missing)
}
if(!identical(old,new[seq_len(nrow(old)),,drop=FALSE])) stop("Internal error: existing HPC rows changed.")
validate_unknown_qc_review(new,samples$sample_id,require_complete=TRUE)
if(nrow(new)!=91L) stop("Expected 91 total decisions after merge; found ",nrow(new),". Inspect the active HPC review.")
if(any(tolower(trimws(new$decision))!="keep")) {
  stop("Existing remove decisions require separate handling to preserve filtered RDS files. No review was changed.")
}
invisible(load_reviewed_broad_thresholds(cfg,require_approved=TRUE))
root <- here::here(cfg$project$outputs_dir)
filtered <- file.path(root,"xenium","annotation","04_unknown_qc_reviewed","rds",
                      paste0(samples$sample_id,"_Consensus_UnknownQC_reviewed.rds"))
if(any(file.exists(filtered))) stop("A filtered QC RDS exists; inspect before an all-keep refresh: ",paste(filtered[file.exists(filtered)],collapse=", "))
paths <- unlist(lapply(c("vz","rl"),function(branch) {
  file.path(root,"xenium",branch,"01_subsets","rds",paste0(ids,"_",toupper(branch),"subset.rds"))
}))
if(any(file.exists(paths))) stop("New subset outputs already exist; inspect before retrying: ",paste(paths[file.exists(paths)],collapse=", "))
inputs <- file.path(get_production_output_paths(cfg)$consensus_rds_dir,
                    paste0(samples$sample_id,"_Consensus_annotated.rds"))
if(any(!file.exists(inputs))) stop("Missing consensus RDS inputs: ",paste(inputs[!file.exists(inputs)],collapse=", "))
cat("PASS: preserved ",nrow(old)," existing HPC rows; ",nrow(missing)," approved additions to append; 91 keep decisions total.\n",sep="")
cat("PASS: consensus inputs present; all 14 new subset outputs absent.\n")
if(identical(Sys.getenv("APPROVED_QC_WRITE"),"yes")) {
  review_path <- resolve_config_path(cfg$manifests$reviewed_unknown_qc,cfg)
  if(nrow(missing)) {
    write.csv(new,review_path,row.names=FALSE,na="")
    saved <- read.csv(review_path,check.names=FALSE,stringsAsFactors=FALSE,na.strings=character())
    for(column in names(old)) {
      expected <- as.character(old[[column]])
      observed <- as.character(saved[[column]][seq_len(nrow(old))])
      if(!identical(expected,observed)) stop("Written review changed existing column: ",column,"; restore the backup.")
    }
    message("Appended only the approved new-sample decisions; existing HPC review values preserved.")
  } else {
    message("All approved decisions already exist; leaving the review file unchanged.")
  }
}
RSCRIPT
)
export APPROVED_QC_WRITE=no
Rscript --max-ppsize=500000 -e "$qc_review_code"
Rscript --max-ppsize=500000 scripts/validate_config.R --check-files
run_dry scripts/xenium_unknown_02_compile_qc_review.R --dry-run
run_dry scripts/xenium_broad_03_compile_calibration.R --dry-run

cohort=config/new_samples_20261014.txt
[[ -f $cohort ]] || cohort=scripts/config/new_samples_20261014.txt
task_map=$(Rscript --max-ppsize=500000 scripts/xenium_unknown_03_apply_qc_review.R --list)
all_ids=$(printf '%s\n' "$task_map" | awk '$1 ~ /^[0-9]+$/ {print $1}' | paste -sd, -)
new_ids=$(printf '%s\n' "$task_map" | awk 'NR==FNR {if(NF) ids[$1]=1; next} $1 ~ /^[0-9]+$/ && $2 in ids {print $1}' "$cohort" - | paste -sd, -)
[[ -n $all_ids && -n $new_ids ]] || { echo "Task mapping is empty." >&2; exit 2; }
echo "QC provenance refresh: $all_ids at concurrency 3"
echo "New VZ/RL splits: $new_ids at concurrency 3 per branch"
[[ $mode == --submit ]] || exit 0

mkdir -p logs outputs/review_gate_backups
backup_root=$(mktemp -d outputs/review_gate_backups/approved_new_samples_20261015_XXXXXX)
cp -p config/reviewed_unknown_qc.csv "$backup_root/reviewed_unknown_qc_before.csv"
cp -p config/reviewed_broad_reference_thresholds.csv "$backup_root/reviewed_broad_reference_thresholds_before.csv"
for path in \
  outputs/xenium/annotation/04_unknown_qc_review/tables \
  outputs/xenium/annotation/04_unknown_qc_review/plots/Xenium_Unknown_QC_All_Samples.pdf \
  outputs/xenium/annotation/04_unknown_qc_reviewed/tables \
  outputs/validation/broad_reference_workflow/03_calibration/tables \
  outputs/xenium/annotation/05_broad_reference_consensus/tables \
  outputs/xenium/annotation/04_unknown_reference_routing/tables; do
  if [[ -e $path ]]; then
    mkdir -p "$backup_root/$(dirname "$path")"
    cp -a "$path" "$backup_root/$path"
  fi
done
export APPROVED_QC_WRITE=yes
Rscript --max-ppsize=500000 -e "$qc_review_code"
cp -p config/reviewed_unknown_qc.csv "$backup_root/reviewed_unknown_qc_applied.csv"
echo "Backups: $backup_root"

submit_job() {
  local job_id
  job_id=$(sbatch --parsable "$@")
  job_id=${job_id%%;*}
  [[ $job_id =~ ^[0-9]+$ ]] || { echo "Invalid sbatch response: $job_id" >&2; exit 2; }
  printf '%s\n' "$job_id"
}
qc_compile=$(submit_job --export=ALL,UNKNOWN_QC_COMPILE_OVERWRITE=true scripts/run_xenium_unknown_02_compile_qc_review.slurm)
echo "QC report refresh: $qc_compile"
broad_compile=$(submit_job --export=ALL,BROAD_CALIBRATION_OVERWRITE=true scripts/run_xenium_broad_03_compile_calibration.slurm)
echo "Broad calibration refresh: $broad_compile"
qc_apply=$(submit_job --array="${all_ids}%3" --dependency="afterok:$qc_compile" --export=ALL,UNKNOWN_QC_APPLY_OVERWRITE=true scripts/run_xenium_unknown_03_apply_qc_review.slurm)
echo "QC provenance refresh: $qc_apply"
broad_apply=$(submit_job --dependency="afterok:$broad_compile" --export=ALL,BROAD_CONSENSUS_OVERWRITE=true scripts/run_xenium_broad_04_apply_reviewed_thresholds.slurm)
echo "Existing broad thresholds: $broad_apply"
combine=$(submit_job --job-name=Xen_unknown_QC_combine --dependency="afterok:$qc_apply" \
  --partition=express --constraint=intel --cpus-per-task=1 --mem=4G --time=00:15:00 \
  --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_unknown_QC_combine_%j.out \
  --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_unknown_QC_combine_%j.err \
  --export=ALL --wrap='bash scripts/finish_xenium_new_samples_20261015.sh --combine-qc')
echo "Combined QC routing: $combine"
routes=$(submit_job --dependency="afterok:$combine:$broad_apply" --export=ALL,UNKNOWN_ROUTE_OVERWRITE=true scripts/run_xenium_unknown_04_build_reference_routes.slurm)
echo "Reference routes: $routes"
vz=$(submit_job --array="${new_ids}%3" --dependency="afterok:$routes" --export=ALL,VZ_SUBSET_OVERWRITE=false scripts/run_xenium_vz_01_subset.slurm)
echo "New VZ subsets: $vz"
rl=$(submit_job --array="${new_ids}%3" --dependency="afterok:$routes" --export=ALL,RL_SUBSET_OVERWRITE=false scripts/run_xenium_rl_01_subset.slurm)
echo "New RL subsets: $rl"
verify=$(submit_job --job-name=Xen_new_split_verify --dependency="afterok:$vz:$rl" \
  --partition=express --constraint=intel --cpus-per-task=1 --mem=4G --time=00:15:00 \
  --output=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_new_split_verify_%j.out \
  --error=/home/acflint/R/Projects/XeniumFCProject/logs/Xen_new_split_verify_%j.err \
  --export=ALL --wrap='bash scripts/finish_xenium_new_samples_20261015.sh --verify')
echo "Final 14-file check: $verify"
printf '%s\n' "$qc_compile" "$broad_compile" "$qc_apply" "$broad_apply" "$combine" "$routes" "$vz" "$rl" "$verify" > "$backup_root/job_ids.txt"
echo "Submitted. This helper has finished; the queued jobs run independently of this terminal."
echo "If submission stops partway through, inspect printed job IDs before retrying."
