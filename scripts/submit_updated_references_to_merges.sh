#!/usr/bin/env bash
# Four dependency chains, separated only by scientific review gates.
# Dry-runs inspect paths/CSV tables/task maps, never large RDS objects.
# Runtime estimates and required manual decisions: UPDATED_REFERENCES_TO_MERGES.md.
set -euo pipefail
if [[ $# -ne 2 || ! $2 =~ ^--(dry-run|submit)$ ]]; then
  echo "Usage: bash scripts/submit_updated_references_to_merges.sh {annotation|regions|postqc|merges} --dry-run|--submit" >&2
  exit 2
fi
phase=$1
mode=$2
case "$phase" in annotation|regions|postqc|merges) ;; *) echo "Unknown phase: $phase" >&2; exit 2 ;; esac
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd -- "$script_dir/.."
run_id=${REFERENCE_REFRESH_RUN_ID:-updated_reference_labels_20261022}
[[ $run_id =~ ^[A-Za-z0-9_-]+$ ]] || { echo "Invalid REFERENCE_REFRESH_RUN_ID" >&2; exit 2; }
run_dir="outputs/reference_label_refresh_runs/$run_id"
ledger="$run_dir/jobs.tsv"

# Use the established HPC R environment loaded by the setup commands in the guide.
# Clear stale modes that could otherwise reuse an old embedding or QC crosswalk.
unset VZ_QC_ADDITIONAL_POSTQC_CLUSTERS
export VZ_INTEGRATE_RESUME_CHECKPOINT=false RL_INTEGRATE_RESUME_CHECKPOINT=false
export VZ_INTEGRATE_PLOTS_ONLY=false RL_INTEGRATE_PLOTS_ONLY=false
export VZ_INTEGRATE_INCLUDE_UNCORRECTED=false RL_INTEGRATE_INCLUDE_UNCORRECTED=false
export VZ_POST_QC_RESUME_CHECKPOINT=false VZ_INTEGRATE_CLUSTER_WORKERS=1 RL_INTEGRATE_CLUSTER_WORKERS=1
export VZ_POST_QC_CLUSTER_WORKERS=1 RL_POST_QC_CLUSTER_WORKERS=1
export VZ_SUBCLUSTER_WORKERS=1 RL_SUBCLUSTER_WORKERS=1

check_r() {
  local output
  if ! output=$(Rscript --max-ppsize=500000 "$@" 2>&1); then
    printf '%s\n' "$output" >&2; return 1
  fi
  printf '%s\n' "$output"
  if [[ $output == *"DRY-RUN FAIL"* ]]; then return 1; fi
}
task_count() {
  local map
  map=$(Rscript --max-ppsize=500000 "$@" --list)
  printf '%s\n' "$map" | awk -F '\t' '
    $1 ~ /^[0-9]+$/ {n++; if ($1 != n || NF < 2 || $2 == "") bad=1}
    END {if (!n || bad) exit 2; print n}'
}
require_previous_phase() {
  local previous=$1 ids accounting
  [[ -s $ledger ]] || { echo "Missing earlier submission ledger: $ledger" >&2; exit 2; }
  ids=$(awk -F '\t' -v phase="$previous" '$1==phase {print $3}' "$ledger" | paste -sd, -)
  [[ -n $ids ]] || { echo "No recorded $previous jobs." >&2; exit 2; }
  accounting=$(sacct -n -X -P -j "$ids" --format=JobID,State,ExitCode)
  # Check each allocation/array element, not .batch/.extern steps.
  while IFS=$'\t' read -r recorded_phase label job expected; do
    [[ $recorded_phase == "$previous" ]] || continue
    if ! printf '%s\n' "$accounting" | awk -F '|' -v id="$job" -v expected="$expected" '
      $1==id || index($1,id "_")==1 {
        n++; if ($2!="COMPLETED" || $3!="0:0") bad=1
      }
      END {if (n < expected || bad) exit 1}'; then
      echo "Wait for/inspect $label ($job); all prior jobs/tasks must be COMPLETED with 0:0." >&2
      exit 2
    fi
  done < "$ledger"
}
queue_stage() {
  local label=$1 script=$2 dependency=$3 array=$4 exports=$5 expected=1 result job
  local submission=(sbatch --parsable)
  [[ -z $dependency ]] || submission+=("--dependency=afterok:$dependency")
  if [[ -n $array ]]; then
    submission+=("--array=$array")
    expected=${array#1-}; expected=${expected%%%*}
  fi
  if [[ -n $exports ]]; then submission+=("--export=ALL,$exports"); else submission+=(--export=ALL); fi
  submission+=("scripts/$script")
  if [[ $mode == --dry-run ]]; then
    printf 'PLAN %s: ' "$label" >&2; printf '%q ' "${submission[@]}" >&2; printf '\n' >&2
    printf 'plan_%s' "$label"
  else
    result=$("${submission[@]}"); job=${result%%;*}
    [[ $job =~ ^[0-9]+$ ]] || { echo "Invalid sbatch result: $result" >&2; exit 2; }
    printf '%s\t%s\t%s\t%s\n' "$phase" "$label" "$job" "$expected" >> "$ledger"
    printf '%s: %s\n' "$label" "$job" >&2
    printf '%s' "$job"
  fi
}

check_r scripts/validate_config.R
check_r -e 'source("scripts/R/config.R"); source("scripts/R/xenium_broad_reference.R"); c <- load_pipeline_config(); validate_taxonomy_against_config(load_reference_broad_taxonomy(c),c); stopifnot("eCN" %in% unlist(c$regional_subsets$rl_broad_labels),"Vascular" %in% unlist(c$regional_subsets$nonregional_broad_labels)); cat("Active configuration:",attr(c,"config_path"),"\n")'
active_paths=$(Rscript --max-ppsize=500000 -e 'source("scripts/R/config.R"); c <- load_pipeline_config(); cat("CONFIG",attr(c,"config_path"),sep="\t"); cat("\nMANIFEST",resolve_config_path(c$manifests$samples,c),sep="\t"); cat("\nRESIDUAL",resolve_config_path("config/reviewed_residual_qc_removals.csv",c),sep="\t"); cat("\nSCIENCE",here::here(c$project$outputs_dir,"references","science","reprocessed","res2","reviewed_retained","Science_res2_reviewed_retained_reprocessed.rds"),sep="\t"); cat("\n")')
config_path=$(printf '%s\n' "$active_paths" | awk -F '\t' '$1=="CONFIG" {print $2}')
manifest_path=$(printf '%s\n' "$active_paths" | awk -F '\t' '$1=="MANIFEST" {print $2}')
residual_decision_path=$(printf '%s\n' "$active_paths" | awk -F '\t' '$1=="RESIDUAL" {print $2}')
science_reviewed_path=$(printf '%s\n' "$active_paths" | awk -F '\t' '$1=="SCIENCE" {print $2}')
[[ -n $science_reviewed_path ]] || { echo "Cannot resolve reviewed Science source path." >&2; exit 2; }
science_needs_reviewed=false
[[ -s $config_path && -s $manifest_path ]] || { echo "Cannot resolve active configuration/manifest." >&2; exit 2; }
if [[ $phase != annotation ]]; then
  cmp -s "$manifest_path" "$run_dir/samples_annotation.csv" || { echo "Sample manifest differs from this run's annotation snapshot." >&2; exit 2; }
  cmp -s "$config_path" "$run_dir/config_annotation.yml" || { echo "Configuration differs from this run's annotation snapshot; review before mixing phases." >&2; exit 2; }
fi
n=$(task_count scripts/xenium_annotate_01_label_transfer_rpca.R --all-samples-res5)
b=$(task_count scripts/xenium_broad_02_score_clusters.R)
[[ $b -eq $((3*n)) ]] || { echo "Broad task map must contain three references per sample." >&2; exit 2; }
printf 'Manifest: %s biological samples; %s reference-score tasks.\n' "$n" "$b"
case "$phase" in
  annotation)
    # Check source inputs, not future standardized/panel outputs.
    if [[ ! -f $science_reviewed_path ]]; then
      science_needs_reviewed=true
      echo "Reviewed Science RDS is missing; prepend approved resolution-2 mapping, removal and reprocessing (estimated 1-8 h; very low confidence)."
      for required in run_science_map_reviewed_res2.slurm science_03_map_reviewed_res2.R R/science_reviewed_mapping.R R/science_ecn_marker_evidence.R R/dotplot_helpers.R; do
        [[ -s scripts/$required ]] || { echo "Missing required script: scripts/$required" >&2; exit 2; }
      done
      check_r 'scripts/fetal cerebellum_science_ACF.R' --resolution=2 --map-reviewed --dry-run
    else
      check_r scripts/science_01_standardize_cell_types.R --dry-run
    fi
    check_r -e 'source("scripts/R/config.R"); c <- load_pipeline_config(); ids <- load_sample_manifest(c)$sample_id; p <- get_production_output_paths(c); inputs <- c(here::here(c$project$outputs_dir,"references","aldinger","rds","Aldinger_seurat_updated.rds"),file.path(c$inputs$published_data_root,"SeppFC","Sepp_hum_sce_final.rds"),resolve_config_path(c$inputs$xenium_panel_genes,c),file.path(p$preprocessing_rds_dir,paste0(ids,"_whole_tissue_Res",p$resolution_tag,".rds"))); ok <- compact_dry_run("Reference sources and existing whole-tissue clusters (no RDS loading)",inputs); if(!ok) quit(status=1)'
    ;;
  regions)
    require_previous_phase annotation
    check_r scripts/xenium_broad_04_apply_reviewed_thresholds.R --dry-run
    check_r scripts/xenium_unknown_03_apply_qc_review.R --dry-run 1
    ;;
  postqc)
    require_previous_phase regions
    : "${VZ_QC_REMOVE_CLUSTERS:?Set freshly reviewed VZ cluster IDs or none; never reuse old IDs automatically}"
    : "${RL_QC_REMOVE_CLUSTERS:?Set freshly reviewed RL cluster IDs or none; never reuse launcher defaults}"
    : "${RESIDUAL_QC_ACTION:?Set remove after updating the reviewed residual decision, or none after explicit no-removal review}"
    [[ $VZ_QC_REMOVE_CLUSTERS =~ ^(none|[0-9]+(,[0-9]+)*)$ && $RL_QC_REMOVE_CLUSTERS =~ ^(none|[0-9]+(,[0-9]+)*)$ ]] || { echo "Invalid QC cluster selection." >&2; exit 2; }
    export VZ_QC_REMOVE_CLUSTERS RL_QC_REMOVE_CLUSTERS
    check_r scripts/xenium_vz_04_review_qc.R "--remove-clusters=$VZ_QC_REMOVE_CLUSTERS" --dry-run
    check_r scripts/xenium_rl_04_review_qc.R "--remove-clusters=$RL_QC_REMOVE_CLUSTERS" --dry-run
    case "$RESIDUAL_QC_ACTION" in
      remove) check_r scripts/xenium_residual_04b_apply_qc.R --dry-run ;;
      none) echo "Reviewed residual no-removal decision: prior decision/filtered outputs will be archived, not reused." ;;
      *) echo "RESIDUAL_QC_ACTION must be remove or none." >&2; exit 2 ;;
    esac
    ;;
  merges)
    require_previous_phase postqc
    check_r scripts/xenium_vz_06_annotate_subclusters.R --resolution=0.8 --apply-reviewed-labels --dry-run
    check_r scripts/xenium_rl_06_annotate_subclusters.R --resolution=0.8 --apply-reviewed-labels --dry-run
    check_r scripts/xenium_residual_05_apply_reviewed_labels.R --dry-run
    ;;
esac
if [[ -f $ledger ]] && awk -F '\t' -v phase="$phase" '$1==phase {found=1} END {exit !found}' "$ledger"; then
  echo "This phase already has recorded jobs; inspect $ledger and recover failed jobs individually." >&2
  exit 2
fi
if [[ $mode == --submit ]]; then
  active_names=$(squeue -h -u acflint -o '%j')
  conflicts=$(printf '%s\n' "$active_names" | awk '/^(Aldinger_|Sepp_|Science_|Xen_(res5_|broad_|unknown_|VZ_|RL_|residual_|Combine|CorrectBroad|CleanMerge|SpatialMerge))/')
  if [[ -n $conflicts ]]; then
    printf 'Wait for these potentially conflicting jobs before submitting a new phase:\n%s\n' "$conflicts" >&2
    exit 2
  fi
  echo "Replacing this phase's existing derived outputs through protected overwrite flags. Source Xenium clusters are reused."
  mkdir -p "$run_dir" logs
  cp -p "$config_path" "$run_dir/config_${phase}.yml"
  cp -p "$manifest_path" "$run_dir/samples_${phase}.csv"
fi

case "$phase" in
  annotation)
    a=$(queue_stage Aldinger_standardize run_aldinger_02_standardize_cell_types.slurm '' '' '')
    ap=$(queue_stage Aldinger_panel run_aldinger_03_subset_gene_panel.slurm "$a" '' '')
    s=$(queue_stage Sepp_standardize run_sepp_01_standardize_cell_types.slurm '' '' '')
    sp=$(queue_stage Sepp_panel run_sepp_02_subset_gene_panel.slurm "$s" '' '')
    science_dependency=''
    if [[ $science_needs_reviewed == true ]]; then
      science_dependency=$(queue_stage Science_reviewed_reprocess run_science_map_reviewed_res2.slurm '' '' 'SCIENCE_REVIEWED_OVERWRITE=true')
    fi
    c=$(queue_stage Science_standardize run_science_01_standardize_cell_types.slurm "$science_dependency" '' 'SCIENCE_REFERENCE_OVERWRITE=true')
    cp=$(queue_stage Science_panel run_science_02_subset_gene_panel.slurm "$c" '' 'SCIENCE_REFERENCE_OVERWRITE=true')
    v=$(queue_stage reference_validation run_xenium_broad_01_validate_references.slurm "$ap:$sp:$cp" '' 'BROAD_REFERENCE_VALIDATION_OVERWRITE=true')
    transfers=()
    for reference in Aldinger Sepp Science; do
      t=$(queue_stage "${reference}_transfer" run_xenium_annotate_01_transfer_resolution5_all_samples.slurm "$v" "1-$n%3" "REFERENCE=$reference,RES5_ALL_ABT_OVERWRITE=true")
      transfers+=("$t")
    done
    dependencies=$(IFS=:; printf '%s' "${transfers[*]}")
    cbuild=$(queue_stage consensus_build run_xenium_annotate_02_build_consensus_resolution5_all_samples_weighted2of3.slurm "$dependencies" '' 'RES5_ALL_CONSENSUS_OVERWRITE=true')
    apply=$(queue_stage consensus_apply run_xenium_annotate_03_apply_consensus_resolution5_all_samples_weighted2of3.slurm "$cbuild" "1-$n%6" 'RES5_ALL_W2OF3_OVERWRITE=true')
    scores=$(queue_stage broad_scores run_xenium_broad_02_score_clusters.slurm "$apply" "1-$b%6" 'BROAD_SCORE_OVERWRITE=true')
    queue_stage broad_calibration run_xenium_broad_03_compile_calibration.slurm "$scores" '' 'BROAD_CALIBRATION_OVERWRITE=true' >/dev/null
    qc=$(queue_stage unknown_qc run_xenium_unknown_01_audit_qc.slurm "$apply" "1-$n%3" 'UNKNOWN_QC_AUDIT_OVERWRITE=true')
    queue_stage unknown_qc_compile run_xenium_unknown_02_compile_qc_review.slurm "$qc" '' 'UNKNOWN_QC_COMPILE_OVERWRITE=true' >/dev/null
    echo "STOP: review broad thresholds and Unknown QC; see the guide before regions."
    ;;
  regions)
    thresholds=$(queue_stage broad_apply run_xenium_broad_04_apply_reviewed_thresholds.slurm '' '' 'BROAD_CONSENSUS_OVERWRITE=true')
    qc=$(queue_stage unknown_qc_apply run_xenium_unknown_03_apply_qc_review.slurm '' "1-$n%3" 'UNKNOWN_QC_APPLY_OVERWRITE=true')
    routes=$(queue_stage routes run_xenium_unknown_04_build_reference_routes.slurm "$thresholds:$qc" '' 'UNKNOWN_ROUTE_OVERWRITE=true')
    for region in vz rl; do
      upper=$(printf '%s' "$region" | tr '[:lower:]' '[:upper:]')
      subset=$(queue_stage "${region}_subset" "run_xenium_${region}_01_subset.slurm" "$routes" "1-$n%3" "${upper}_SUBSET_OVERWRITE=true")
      merge=$(queue_stage "${region}_merge" "run_xenium_${region}_02_merge_samples.slurm" "$subset" '' "${upper}_MERGE_OVERWRITE=true")
      integrated=$(queue_stage "${region}_integrate" "run_xenium_${region}_03_integrate.slurm" "$merge" '' "${upper}_INTEGRATE_OVERWRITE=true")
      queue_stage "${region}_qc_review" "run_xenium_${region}_04_review_qc.slurm" "$integrated" '' "${upper}_QC_MODE=review,${upper}_QC_OVERWRITE=true" >/dev/null
      queue_stage "${region}_harmony_audit" "run_xenium_${region}_03a_audit_harmony.slurm" "$integrated" '' "${upper}_HARMONY_AUDIT_OVERWRITE=true" >/dev/null
    done
    subset=$(queue_stage residual_subset run_xenium_residual_01_subset_samples.slurm "$routes" "1-$n%3" 'RESIDUAL_SUBSET_OVERWRITE=true')
    merge=$(queue_stage residual_merge run_xenium_residual_02_merge_samples.slurm "$subset" '' 'RESIDUAL_MERGE_OVERWRITE=true')
    integrated=$(queue_stage residual_integrate run_xenium_residual_03_integrate.slurm "$merge" '' 'RESIDUAL_INTEGRATE_OVERWRITE=true')
    queue_stage residual_qc_review run_xenium_residual_04a_review_qc.slurm "$integrated" '' 'RESIDUAL_QC_OVERWRITE=true' >/dev/null
    echo "STOP: review fresh VZ/RL/residual QC; do not reuse previous removal IDs."
    ;;
  postqc)
    for region in vz rl; do
      upper=$(printf '%s' "$region" | tr '[:lower:]' '[:upper:]')
      qc=$(queue_stage "${region}_qc_apply" "run_xenium_${region}_04_review_qc.slurm" '' '' "${upper}_QC_MODE=apply,${upper}_QC_OVERWRITE=true")
      integrated=$(queue_stage "${region}_postqc" "run_xenium_${region}_05_reintegrate_post_qc.slurm" "$qc" '' "${upper}_POST_QC_OVERWRITE=true")
      # Current review mode builds/reuses provenance-matched replicate evidence itself.
      queue_stage "${region}_marker_review" "run_xenium_${region}_06_annotate_subclusters.slurm" "$integrated" '' "${upper}_SUBCLUSTER_MODE=review,${upper}_SUBCLUSTER_RESOLUTION=0.8,${upper}_SUBCLUSTER_OVERWRITE=true" >/dev/null
    done
    residual_dependency=''
    if [[ $RESIDUAL_QC_ACTION == remove ]]; then
      residual_dependency=$(queue_stage residual_qc_apply run_xenium_residual_04b_apply_qc.slurm '' '' 'RESIDUAL_QC_APPLY_OVERWRITE=true')
    elif [[ $mode == --submit ]]; then
      archive="$run_dir/residual_prior_qc"
      mkdir "$archive"
      for old in "$residual_decision_path" config/reviewed_residual_qc_removals.csv scripts/config/reviewed_residual_qc_removals.csv outputs/xenium/residual/04_qc_filtered; do
        if [[ -e $old ]]; then
          destination="$archive/${old//\//__}"
          printf '%s\t%s\n' "$old" "$destination" >> "$archive/moves.tsv"
          mv -- "$old" "$destination"
        fi
      done
      printf 'Reviewed no-removal decision for refreshed residual object.\n' > "$archive/no_removal.txt"
    fi
    queue_stage residual_marker_review run_xenium_residual_04_review_subclusters.slurm "$residual_dependency" '' 'RESIDUAL_REVIEW_OVERWRITE=true' >/dev/null
    echo "STOP: approve NEW VZ/RL maps, input-count guards/validator, residual labels and any new fine-to-broad identities."
    ;;
  merges)
    maps=()
    for region in vz rl; do
      upper=$(printf '%s' "$region" | tr '[:lower:]' '[:upper:]')
      labels=$(queue_stage "${region}_labels" "run_xenium_${region}_06_annotate_subclusters.slurm" '' '' "${upper}_SUBCLUSTER_MODE=apply,${upper}_SUBCLUSTER_RESOLUTION=0.8,${upper}_SUBCLUSTER_OVERWRITE=true")
      map=$(queue_stage "${region}_map" "run_xenium_${region}_07_map_to_samples.slurm" "$labels" '' "${upper}_MAPPING_OVERWRITE=true")
      maps+=("$map")
    done
    labels=$(queue_stage residual_labels run_xenium_residual_05_apply_reviewed_labels.slurm '' '' 'RESIDUAL_LABEL_OVERWRITE=true')
    map=$(queue_stage residual_map run_xenium_residual_06_map_to_samples.slurm "$labels" '' 'RESIDUAL_MAPPING_OVERWRITE=true')
    maps+=("$map")
    dependencies=$(IFS=:; printf '%s' "${maps[*]}")
    combined=$(queue_stage combined_labels run_xenium_vz_rl_01_combine_labels.slurm "$dependencies" '' 'COMBINE_OVERWRITE=true')
    broad=$(queue_stage corrected_broad run_xenium_vz_rl_02_correct_broad_labels.slurm "$combined" '' 'CORRECTED_BROAD_OVERWRITE=true')
    queue_stage nonspatial_merge run_xenium_vz_rl_03_nonspatial_merge.slurm "$broad" '' 'CLEAN_MERGE_OVERWRITE=true' >/dev/null
    queue_stage spatial_merge run_xenium_vz_rl_03_spatial_merge.slurm "$broad" '' 'SPATIAL_MERGE_OVERWRITE=true' >/dev/null
    echo "Endpoint: refreshed nonspatial and spatial merged objects; merge processing/integration is separate."
    ;;
esac
[[ $mode != --dry-run ]] || echo "DRY-RUN: no Slurm jobs submitted and no outputs changed."
