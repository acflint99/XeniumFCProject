# Cohort-wide outputs after adding seven samples

Scope: production resolution-5 whole-tissue outputs through Unknown review,
broad routing, and the requested new VZ/RL subsets. The manifest now has 41
samples. Scientific settings, output names, and the original 34 IDs are retained.
This inventory comes from current scripts, not a live HPC filesystem scan.
Run the included read-only audit on Cheaha to inspect actual file coverage.

## Upload and refresh the consensus PDFs

Upload `xenium_merged_reports_update_20261015.tar.gz` to
`/home/acflint/R/Projects/XeniumFCProject/`. In the Intel terminal:

```bash
cd /home/acflint/R/Projects/XeniumFCProject
source_backup="scripts_backup_merged_reports_$(date +%Y%m%d_%H%M%S)"
mkdir "$source_backup"
tar -tzf xenium_merged_reports_update_20261015.tar.gz | while IFS= read -r path; do
  if [[ -f "$path" ]]; then
    mkdir -p "$source_backup/$(dirname "$path")"
    cp -p "$path" "$source_backup/$path"
  fi
done
tar -xzf xenium_merged_reports_update_20261015.tar.gz
bash scripts/refresh_xenium_consensus_reports_20261015.sh --dry-run
```

Preflight should take 1–3 minutes: metadata, paths, TIFF presence, and small
CSV inventories only; no Seurat RDS objects are loaded. It prints the PCW-sorted
41-page map, validates both report builders, and inventories other merged files.
PDF page counts require the existing `pdfinfo` command. CSV coverage in
conflict/unrouted reports can legitimately be below 41; inspect their intended
population before interpreting a coverage count as a missing-sample error.

After all report preflight checks pass:

```bash
bash scripts/refresh_xenium_consensus_reports_20261015.sh --submit
```

This archives the previous three PDFs and both consensus PNG page-cache
directories under `outputs/report_backups/`, then submits fresh page rendering,
a dependent consensus merge, and an independent Aldinger merge. It preserves
all analysis RDS files, input TIFFs, and overwrite protection on fresh reports.
If submission fails partway through, inspect the printed job IDs before retrying.

Consensus rendering loads one consensus RDS per array task and builds both
variants: estimate 3–15 minutes per task, 21–105 minutes for 41 tasks at
concurrency six, plus 1–5 minutes for the merge. These are low-confidence
estimates from static inspection and the new samples' observed 1:26–2:49
consensus-object/plot jobs, which are a related but different stage. The renderer
requests 2 hours/64 GB per task; merge requests 30 minutes/16 GB. The Aldinger
TIFF-only report is estimated at 3–15 minutes, low confidence from 369 input
TIFF reads and raster downsizing; requested limit remains 1 hour/16 GB.
All jobs retain the Intel constraint and existing module/renv setup. Queue wait
is excluded. Collect Elapsed/MaxRSS after these jobs to refine estimates.

Both consensus PDFs must contain **41 PCW-ordered pages**. The Aldinger report
must contain **41 manifest-ordered pages**. Rebuilding only global tasks 35–41
is incorrect for the PCW-sorted report: its task IDs are page numbers, and new
PCW12/PCW17 samples change page positions and cache names for older samples.

## Outputs that must include the new samples

Paths below are relative to `outputs/`.

| Output family | Main outputs | How to refresh / readiness |
|---|---|---|
| Consensus PDFs | `xenium/annotation/03_consensus_labels_weighted_2of3/plots/Xenium_Consensus_Res5_Weighted2of3_All_Samples_Report.pdf` and `..._All_Samples_KnownOnly_Report.pdf` | Included in the new report helper; render 41 pages for each variant, then merge. |
| Aldinger combined plots | `xenium/annotation/01_label_transfer/aldinger/plots/Xenium_Aldinger_Res5_All_Samples_Plot_Report.pdf` | Included in the new report helper. The prior 34-sample guard is replaced by validated manifest membership/count. No label-transfer rerun. No equivalent Sepp/Science merged-report builder was found in this local script set. |
| Unknown QC PDF and combined tables | `xenium/annotation/04_unknown_qc_review/plots/Xenium_Unknown_QC_All_Samples.pdf`; `tables/Xenium_Unknown_QC_Evidence.csv`, `Xenium_Unknown_QC_Review_Template.csv`, `Xenium_Unknown_QC_Audit_manifest.csv` | `run_xenium_unknown_02_compile_qc_review.slurm`, protected `UNKNOWN_QC_COMPILE_OVERWRITE=true`. Expected 41 audit-manifest rows and all 91 reviewed population keys in the evidence. Do not replace the approved review with the generated blank template. |
| Broad calibration | Six CSVs in `validation/broad_reference_workflow/03_calibration/tables/`, including `Xenium_Broad_All_Cluster_Scores.csv` and `Xenium_Broad_Reference_Threshold_Sensitivity.csv` | `run_xenium_broad_03_compile_calibration.slurm`, protected `BROAD_CALIBRATION_OVERWRITE=true`; compile 123 sample/reference inputs. Existing reviewed thresholds stay unchanged. |
| Broad-reference consensus | Four combined/provenance CSVs plus 41 per-sample CSVs under `xenium/annotation/05_broad_reference_consensus/tables/` | After calibration, `run_xenium_broad_04_apply_reviewed_thresholds.slurm`, `BROAD_CONSENSUS_OVERWRITE=true`. Included in the approved-split helper. |
| Reviewed QC input manifest | `xenium/annotation/04_unknown_qc_reviewed/tables/Xenium_Unknown_QC_routing_manifest.csv` and per-sample routing/removal CSVs | Refresh all 41 QC-application tasks after appending the approved review, then combine. All decisions are keep, so consensus RDS files are reused. Included in the approved-split helper. |
| Reference routing | `xenium/annotation/04_unknown_reference_routing/tables/Xenium_Unknown_Reference_Routes.csv`, `Xenium_Unknown_Reviewed_Input_manifest.csv`, validation and summary CSVs | After QC application and broad consensus, `run_xenium_unknown_04_build_reference_routes.slurm`, `UNKNOWN_ROUTE_OVERWRITE=true`. Included in the approved-split helper. |
| Unknown reference/marker evidence | `xenium/annotation/04_unknown_reference_routing/tables/Xenium_Unknown_Reference_Score_Marker_Evidence.csv` | Optional report, after fresh routes. Commands below. |
| Unrouted Unknown summaries | Three CSVs under `xenium/annotation/05_unrouted_unknown_report/tables/` | After fresh routes; reports only retained Other/unresolved populations, so not every sample necessarily appears. Commands below. |
| Consensus validation and counts | `validation/resolution5_through_consensus/reference_validation.csv`, `sample_validation.csv`, `consensus_label_counts.csv`, `validation_summary.txt` | `validate_outputs_through_consensus.R` derives sample membership from the manifest. It reads many large RDS files; estimated 1–5 hours for 41 samples, low confidence; schedule on HPC, not as local preflight. No built-in overwrite protection, so back up the report directory before intentionally rerunning. Not included in the lightweight refresh helper. |
| Whole-tissue cell proportions | `xenium/annotation/03_consensus_labels_weighted_2of3/plots/proportions/XeniumConsensus_Res5_Weighted2of3_ClusterPropPlot.pdf` and `.tif` | `xenium_annotate_03c_plot_proportions.R` discovers all consensus RDS files, rather than validating the manifest. It has no dry-run/overwrite guard. Flagged for a focused guarded update before automatic refresh. It pools cells by PCW within RL/UBC/Granule/Purkinje/GABA; the new PCW12/PCW17 cells change those pooled proportions without changing that definition. Estimated 10–60 minutes, low confidence from 41 RDS reads. Not run. |

The initial preprocessing drivers and transfer/consensus builders write their
main outputs per sample. The seven new samples already completed those stages;
there is no separate production all-sample preprocessing RDS merge to rebuild
before the regional split. Reference preparation, gene-panel inputs, calibration
pilots, k.anchor/k.weight/reference-cap sensitivity cohorts, and legacy
resolution comparisons are not automatically expanded by adding samples.

Regional all-sample spatial reports, regional count/proportion reports, combined
VZ/RL/residual label outputs, and merged/integrated regional RDS files belong
to later reviewed branch analyses. The new samples have not gone through
regional integration/subcluster annotation/mapping yet. Including them there
requires those stages and fresh cluster review; it is not a PDF-only refresh.
Giotto/SpaTrack and historical directories remain outside this audit's analysis
scope. The current local script set contains additional dated 34-sample wording
or guards in downstream regional helpers; those are not changed by this update.

## Optional routing-dependent merged reports

Wait until the approved-split helper's route job has succeeded. These dry-runs
inspect CSV/workbook paths and small tables (no large Seurat RDS objects),
expected under 1 minute each:

```bash
Rscript --max-ppsize=500000 scripts/xenium_unknown_04a_compile_reference_evidence.R --dry-run
Rscript --max-ppsize=500000 scripts/xenium_unknown_05_report_unrouted.R --dry-run
```

After both pass, preserve prior versions and submit the protected updates:

```bash
routing_report_backup="outputs/report_backups/routing_reports_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$routing_report_backup"
if [[ -f outputs/xenium/annotation/04_unknown_reference_routing/tables/Xenium_Unknown_Reference_Score_Marker_Evidence.csv ]]; then
  cp -p outputs/xenium/annotation/04_unknown_reference_routing/tables/Xenium_Unknown_Reference_Score_Marker_Evidence.csv "$routing_report_backup/"
fi
if [[ -d outputs/xenium/annotation/05_unrouted_unknown_report/tables ]]; then
  cp -a outputs/xenium/annotation/05_unrouted_unknown_report/tables "$routing_report_backup/unrouted_tables"
fi
sbatch --export=ALL,UNKNOWN_EVIDENCE_OVERWRITE=true scripts/run_xenium_unknown_04a_compile_reference_evidence.slurm
sbatch --export=ALL,UNKNOWN_REPORT_OVERWRITE=true scripts/run_xenium_unknown_05_report_unrouted.slurm
```

Expected active runtime: 1–5 minutes for evidence, under 1 minute for unrouted
summaries, low confidence from small-table/workbook processing. Both request
15 minutes, 4 GB, one CPU on express/Intel. These reports do not change routes,
labels, thresholds, or QC decisions. They are not part of the consensus-only helper.

## Check actual merged-output coverage on Cheaha

```bash
Rscript --max-ppsize=500000 scripts/audit_xenium_merged_outputs.R --dry-run
```

This read-only inventory should take under 1 minute. It prints file existence,
CSV sample coverage where applicable, PDF page counts, and modification times.
It does not load RDS objects or write reports. Missing files and older coverage
identify work still outstanding; PRESENT alone is not proof of current contents.
Collect job states and exit codes after merges finish, and verify that both
consensus PDFs and the Aldinger PDF have 41 pages.

The Aldinger R script was syntax-parsed, all changed shell launchers were
checked with bash -n, and task mappings were checked locally. HPC preflight
and report rendering have not been run by Codex. Review the diff and make a
separate Git checkpoint when ready; no package changes or renv snapshot are needed.
