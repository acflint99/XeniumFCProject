# Approved labels through combined regional labels, updated 2026-10-19

Run from `/home/acflint/R/Projects/XeniumFCProject`. Active manifest: 41 samples. These are cohort jobs, not sample arrays. Queue wait is excluded from all estimates.

## Files to sync

Required current changes, under `scripts/`: `xenium_vz_06_annotate_subclusters.R`, `color_palette.R`, `xenium_vz_plot_sample.R`, `validate_xenium_regional_mapping.R`, and `run_xenium_vz_rl_01_combine_labels.slurm`. Also ensure the previously updated `R/dotplot_helpers.R` is synced; `workflow_dotplot_grob` must accept `strip_colors`. Existing VZ/RL/residual mapping launchers and drivers are reused without changes.

VZ Stage 5 and the new Stage-6 REVIEW are complete. Do not repeat QC, reintegration or raw marker review for this label-only update. The exact approved map is:

```r
c(
  "0"="Patterning PCs", "1"="BG", "2"="Maturing PCs", "3"="MLIs",
  "4"="RG Progenitors", "5"="Early-born PCs", "6"="Late-born PCs",
  "7"="eCN", "8"="Pre-OPCs", "9"="GABA Progenitors",
  "10"="Cycling Glia", "11"="iCN", "12"="Astrocytes/Ependyma",
  "13"="Cycling Glia", "14"="RL VZ", "15"="Maturing PCs",
  "16"="Golgi Cells", "17"="Differentiating GCs", "18"="OPCs", "19"="VZPs",
  "20"="GABA Interneurons", "21"="iCN", "22"="Unresolved",
  "23"="Patterning PCs", "24"="Maturing PCs"
)
```

This is a small in-memory vector, effectively instantaneous to construct. The label-application script checks the reviewed 1,447,477-cell cohort as well as existing exact-input provenance and exact cluster-set guards. No cells, sample IDs, expression values, integration settings or clustering parameters change. The tiny cluster 22 remains present as Unresolved; 23 is explicitly user-assigned Patterning PCs rather than independently supported as a subtype.

## 1. Apply VZ labels and regenerate annotated plots

Path/small-CSV preflight, under one minute, no large RDS load:

```bash
Rscript --max-ppsize=500000 scripts/xenium_vz_06_annotate_subclusters.R \
  --apply-reviewed-labels --resolution=0.8 --dry-run
```

If checks pass and application outputs do not already exist:

```bash
VZ_SUBCLUSTER_MODE=apply \
VZ_SUBCLUSTER_RESOLUTION=0.8 \
VZ_SUBCLUSTER_OVERWRITE=false \
sbatch --export=ALL scripts/run_xenium_vz_06_annotate_subclusters.slurm
```

Expected 20–60 minutes, low confidence from the previous label application reaching plotting at about 12 minutes and the large current cohort. Existing request: 256G, one CPU, 12-hour limit. This writes the master annotated RDS, approved label manifest, QC summary, annotated UMAPs, compact DotPlots and heatmap. Raw diagnostic review files remain separate. If existing application outputs are reported, preserve/review those before explicitly using the supported overwrite flag; do not repeatedly submit after a partial failure.

Production VZ plots use exactly the 20 observed approved label names. Cluster 17 is now provisionally Differentiating GCs based on review of the VZ marker evidence; its former Mature UBCs identity is removed from the VZ taxonomy and Other DotPlot. The VZ Granule panel uses NEUROD1/CBLN1/RELN, with the broad Granule green header. The independent RL and residual Mature UBCs labels remain unchanged. Neuroblasts/Immature UBCs are absent from the current VZ taxonomy; legacy shared colors remain available to other branches. SP9 remains in the GABA panel. No cycling-marker columns are added. GABA Interneurons appears in GABA/lineage spatial highlights and facets; Unresolved is grey in global/Other plots and has no forced marker panel or lineage. eCN retains its subcluster orange-red. Per-sample production panels continue to filter absent marker groups and identity levels. Raw diagnostic cycling markers remain intentionally available for QC review.

## 2. Map the new VZ labels to individual samples

Wait for the application job to finish COMPLETED/0:0. Preserve stale Stage-7 outputs if present (a move, not a new large copy):

```bash
if [[ -d outputs/xenium/vz/07_mapped ]]; then
  mv outputs/xenium/vz/07_mapped \
    "outputs/xenium/vz/07_mapped_before_approved_map_$(date +%Y%m%d_%H%M%S)"
fi
Rscript --max-ppsize=500000 scripts/xenium_vz_07_map_to_samples.R --dry-run
```

After PASS:

```bash
VZ_MAPPING_OVERWRITE=false VZ_MAPPING_WORKERS=1 \
sbatch --export=ALL scripts/run_xenium_vz_07_map_to_samples.slurm
```

Expected 10–25 minutes, based on the previous 13:29 VZ mapping run (70.77 GiB peak) with a slightly larger cohort. Existing request: 128G, one CPU, four-hour limit. The driver checks per-sample mapped counts and total master-cell coverage; every mapped VZ cell retains its approved label, including Unresolved.

After success, lightweight manifest validation, under one minute; no large RDS load:

```bash
Rscript --max-ppsize=500000 scripts/validate_xenium_regional_mapping.R VZ
```

## 3. RL: reuse current approved mapping when valid

The RL labels and QC were not changed by the VZ decision. Previous RL sample reports completed across all 41 tasks, so first check existing mapping rather than automatically repeating it:

```bash
Rscript --max-ppsize=500000 scripts/validate_xenium_regional_mapping.R RL
```

Expected under one minute, reading CSVs and file metadata only. A PASS means the sample IDs/order, output paths, nonempty file sizes, mapped counts and source/output timestamps are consistent. It does not inspect RDS contents or prove label validity independently.

If labels are already applied but the mapping is missing/stale, preserve existing mapping outputs and rerun mapping only:

```bash
if [[ -d outputs/xenium/rl/07_mapped ]]; then
  mv outputs/xenium/rl/07_mapped \
    "outputs/xenium/rl/07_mapped_before_remap_$(date +%Y%m%d_%H%M%S)"
fi
Rscript --max-ppsize=500000 scripts/xenium_rl_07_map_to_samples.R --dry-run
```

After PASS:

```bash
RL_MAPPING_OVERWRITE=false RL_MAPPING_WORKERS=1 \
sbatch --export=ALL scripts/run_xenium_rl_07_map_to_samples.slurm
```

Expected 10–30 minutes, low confidence from the analogous sequential VZ mapping job and 41 whole-tissue object reads/writes. Existing request: 128G, one CPU, four-hour limit.

ONLY if the approved RL master/label manifest is missing and its current Stage-6 review is ready, use the existing approved RL map:

```bash
Rscript --max-ppsize=500000 scripts/xenium_rl_06_annotate_subclusters.R \
  --apply-reviewed-labels --resolution=0.8 --dry-run
```

After PASS and confirming no existing application outputs need replacement:

```bash
RL_SUBCLUSTER_MODE=apply RL_SUBCLUSTER_RESOLUTION=0.8 RL_SUBCLUSTER_OVERWRITE=false \
sbatch --export=ALL scripts/run_xenium_rl_06_annotate_subclusters.slurm
```

Expected 15–45 minutes, based on the previous 18:56 approved RL application; request 256G/one CPU/12 hours. Wait for completion before the RL mapping preflight/submission. Current approved RL identities are unchanged by this update. If its input was reclustered, obtain a new review rather than applying an old cluster map.

## 4. Residual: already complete, validate and reuse

Recorded successful jobs: label application 40969781 (00:00:40), mapping 40969788 (00:22:06). The resource history records prior verification of all 41 output paths, 9,356 mapped cells and 977 QC exclusions. No residual label changes were made here.

```bash
Rscript --max-ppsize=500000 scripts/validate_xenium_regional_mapping.R RESIDUAL
```

Expected under one minute, no RDS load. If PASS, skip residual application and mapping.

For a genuinely missing/stale residual mapping, after preserving existing outputs and confirming current approved Stage-5 labels:

```bash
if [[ -d outputs/xenium/residual/06_mapped ]]; then
  mv outputs/xenium/residual/06_mapped \
    "outputs/xenium/residual/06_mapped_before_remap_$(date +%Y%m%d_%H%M%S)"
fi
Rscript --max-ppsize=500000 scripts/xenium_residual_06_map_to_samples.R --dry-run
```

After PASS:

```bash
RESIDUAL_MAPPING_OVERWRITE=false RESIDUAL_MAPPING_WORKERS=1 \
sbatch --export=ALL scripts/run_xenium_residual_06_map_to_samples.slurm
```

Expected 15–35 minutes from measured 22:06; request 128G/one CPU/four hours. If the residual Stage-5 master is genuinely absent (not merely mapping), the existing approved residual map can be applied first:

```bash
Rscript --max-ppsize=500000 scripts/xenium_residual_05_apply_reviewed_labels.R --dry-run
```

After PASS and ensuring existing label outputs are preserved:

```bash
RESIDUAL_LABEL_OVERWRITE=false \
sbatch --export=ALL scripts/run_xenium_residual_05_apply_reviewed_labels.slurm
```

Expected 1–5 minutes from the measured 40-second run; request 96G/one CPU/one hour. Wait for completion before the mapping preflight/submission.

## 5. Confirm all three mappings before combining branches

```bash
Rscript --max-ppsize=500000 scripts/validate_xenium_regional_mapping.R VZ RL RESIDUAL
```

Expected under one minute; output should show 41/41 for each branch. It additionally checks VZ's exact approved 0–24 label manifest and 1,447,477 mapped-cell total. It deliberately does not read any large RDS. VZ sample reports/counts and joint combined VZ/RL/residual analyses remain stale until regenerated. The shared palette/master order now recognizes GABA Interneurons and Unresolved for the combined stage.

## 6. Combine the three regional mappings

Preserve the old combined output tree before submitting the replacement. Submit the combined-label launcher only after the VZ map job completes; a Slurm `afterok` dependency can enforce this. The updated launcher validates all three mapping manifests and requires a `DRY-RUN PASS` from the combined driver before starting its large RDS work. Those checks are expected to take under one minute and inspect paths, CSVs, and file metadata without loading RDS objects.

```bash
if [[ -d outputs/xenium/vz_rl/01_combined_labels ]]; then
  mv outputs/xenium/vz_rl/01_combined_labels \
    "outputs/xenium/vz_rl/01_combined_labels_before_approved_map_$(date +%Y%m%d_%H%M%S)"
fi
sbatch --dependency=afterok:VZ_MAP_JOB_ID --export=ALL,COMBINE_OVERWRITE=false \
  scripts/run_xenium_vz_rl_01_combine_labels.slurm
```

Combined-label job 40971846 completed in 44:57 for the 41-sample cohort. A similar rerun is estimated at 40–90 minutes, low confidence from one observed run. The launcher requests 24 hours, 128G, and one CPU on `medium`. Queue wait is excluded. Existing combined results and downstream outputs are stale when the VZ labels change; retain the archived tree for rollback. Collect `MaxRSS` before future resource changes.

## 7. Correct broad labels in new outputs before merging (approved 2026-10-19)

Combined-label job 40971846 completed in 44:57. Its existing RDS files and
plots remain unchanged. The approved fine-to-broad map is in
`R/corrected_broad_labels.R`; it applies to every retained cell, preserves
`consensus_label`, and adds `corrected_broad_label`. The new Stage 01c writes
41 corrected sample RDS files, 82 new broad spatial plots, per-sample count
and transition tables, and a 41-row provenance manifest under
`outputs/xenium/vz_rl/02_corrected_broad/`. PC states map to Purkinje;
Pre-OPCs and OPCs map to OPC; iCN maps to GABA; RL VZ/SVZ map to RL; RG/BG,
cycling glial, and astrocyte/ependymal states map to Glia. eCN, Pericytes,
Cycling Cells, and Unresolved remain separate broad identities.

```bash
Rscript --max-ppsize=500000 scripts/xenium_vz_rl_02_correct_broad_labels.R --dry-run
```

This path-only preflight takes under one minute and does not load RDS objects.
After PASS with `inputs=41/41` and `outputs-existing=0/206`, submit:

```bash
sbatch scripts/run_xenium_vz_rl_02_correct_broad_labels.slurm
```

Job 40973926 completed in 12:56 for 41 samples. A repeat is estimated at
10–30 minutes, low confidence from one observed run; batch MaxRSS was 16.44
GiB. The post-job validator passed for 41/41 samples, 2,355,487 cells,
82 plots, and the count/transition tables without loading RDS contents.
The request is 12
hours/128G/one CPU on `short`; queue wait is excluded. Check the new manifest,
transition tables, plots, and `sacct` before the regular or spatial merge.
The merge launchers now require the corrected-broad manifest and objects and
will fail early if they are missing or stale. No new QC, subclustering, or
original combined-label run is required.

```bash
Rscript --max-ppsize=500000 scripts/validate_xenium_vz_rl_02_correct_broad.R
```

This audits all 41 corrected sample paths, 82 nonempty new plots, small
count/transition tables, and mapping provenance without loading RDS objects.
Estimate under two minutes, depending on filesystem latency. Inspect at least
the FB124_X_G corrected-broad spatial plot alongside its original consensus
plot before proceeding to `xenium_vz_rl_03_nonspatial_merge.R --dry-run`.

After each important job, collect resource evidence:

```bash
sacct -j JOB_NUMBER --format=JobID,Elapsed,State,MaxRSS,ReqMem,AllocCPUS,ExitCode
```

Local checks: changed R scripts parse; approved 25-entry map exactly matches the 20-name display order and complete palettes; current assay evidence contains every production marker; active-panel filtering excludes unused groups and cycling-marker columns. A small artificial Seurat DotPlot was inspected for layout. The local Seurat/ggplot2 pairing required applying facets through the current ggplot API in the test harness; no package/environment or production faceting code was changed. HPC production rendering must still be confirmed. Supported local preflights reported unavailable HPC inputs, not successful cohort readiness. No full analysis or output replacement ran locally. No packages or renv.lock changed.

Rollback is the prior scripts plus matching prior annotated/mapped outputs, kept as coherent sets. This local directory is not a Git checkout; in the actual project checkout, review git status/diff and validation, then commit this logical change when ready. No commit or push was performed.
