# Fetal Cerebellum Xenium Analysis

This directory contains the R and Slurm scripts used to analyze 10x Genomics Xenium spatial-transcriptomics data from the developing human cerebellum. The workflow crops each Xenium field to cerebellar tissue, performs cell-level quality control and Seurat clustering, transfers cell-type labels from published single-cell references, and then performs focused analyses of the ventricular zone (VZ) and rhombic lip (RL) lineages.

The downstream scripts generate spatial maps, cell-type proportions, cross-study comparisons, spatially variable genes, neighborhood enrichment, ligand–receptor results, gene-expression gradients, and trajectory-ready `h5ad` files.

> [!IMPORTANT]
> This is an active research workflow, not a packaged command-line pipeline. Intentional focused cohorts, Slurm array bounds, cluster selections, and some absolute HPC paths must be checked before each run. The execution order below is inferred from script inputs and outputs.

## Biological scope

The project integrates Xenium measurements with three fetal cerebellum single-cell or single-nucleus references:

- **Aldinger**
- **Sepp**
- **Science**

The transferred broad annotations are subsequently used to isolate and refine two developmental compartments:

- **VZ-related cells:** explicit VZ cells, glia, GABAergic cells, Purkinje cells, and OPCs
- **RL-related cells:** rhombic-lip cells, granule-lineage cells, and unipolar brush cells (UBCs)

Exact labels are controlled by the active identities and metadata in each Seurat object. The shared cell-type order and plotting colors are defined in `color_palette.R`.

For recoloring plots in their original output directories with the original
UMAP formatting, see [RESTORE_ORIGINAL_PLOT_FORMAT_20261020.md](RESTORE_ORIGINAL_PLOT_FORMAT_20261020.md).
The reference redraw reads saved objects without rewriting RDS files; other
stages use their native plot-only drivers. The earlier generic snapshot
exporter is optional historical tooling, not the recommended redraw route.
For replacing the Xenium originals, including dependent PDF reports, use
[XENIUM_ORIGINAL_PLOTS_20261020.md](XENIUM_ORIGINAL_PLOTS_20261020.md).
Its default path-only preflight takes about 1–3 minutes; plotting and Slurm
array estimates are documented in that guide.

## Workflow overview

Current broad-label convention: **Vascular = Endothelial + Pericytes**.
Current reference preparation and downstream broad annotation use Vascular;
Pericytes remains available as a fine subcluster label. See
[VASCULAR_LABEL_UPDATE_20261022.md](VASCULAR_LABEL_UPDATE_20261022.md) for the
source-label exceptions, required regeneration, validation and runtime estimates.

For the 2026-10-22 Science/Sepp/Aldinger label refresh, use
[`UPDATED_REFERENCES_TO_MERGES.md`](UPDATED_REFERENCES_TO_MERGES.md) and
`submit_updated_references_to_merges.sh`. It rebuilds reference preparations,
all manifest-defined transfers/consensus and regional branches through both
merges plus post-merge reference validation, stopping for fresh threshold/QC/subcluster-label decisions. eCN enters
RL and Vascular remains nonregional, as explicitly approved on 2026-10-22.
It derives 41-sample/123-task bounds; older 34-sample submitters are not used.

```text
Xenium import -> cerebellum crop -> cell QC -> PCA/UMAP/neighbors
                 |
                 v
Production resolution-5 clusters + panel-matched Aldinger / Sepp / Science references
                 |
                 v
RPCA transfer -> weighted 2-of-3 consensus -> reviewed Unknown QC / broad routing
                 |
        +--------+---------+
        |        |         |
       VZ       RL      Residual
        |        |         |
  merge / integration -> QC review and approved application
        |        |         |
  clustree + subsampling stability -> marker review -> approved fine labels
        |        |         |
        +--------+---------+
                 |
      map labels to samples -> combine / corrected broad labels
                 |
        +--------+---------+
        |                  |
  nonspatial merge    spatial merge
        |
  Spearman + MetaNeighborUS validation against original references
```

See [WORKFLOW_ALIGNMENT_20261022.md](WORKFLOW_ALIGNMENT_20261022.md) for the
current script/diagram alignment, residual clustree entry point, protected
outputs, runtime estimates and validation limits. All three branches use the
same review sequence. Established numerical settings remain branch-specific;
VZ/RL post-QC processing and residual's retained embedding are preserved.


## Expected project layout

The scripts assume they are stored in a project whose root can be found by the `here` package:

```text
XeniumFCProject/
├── config/
│   ├── config.yml
│   ├── samples.csv
│   └── slides.csv
├── data/
│   └── FCXeniumProject/
│       └── <sample>/
│           ├── <standard Xenium output files>
│           └── cerebellum_cells_stats.csv
├── metadata/
│   └── samples_meta.xlsx
├── scripts/
│   └── <this directory>
├── outputs/
├── logs/
├── renv/
└── renv.lock
```

For each Xenium sample, `cerebellum_cells_stats.csv` must contain a `Cell ID` column. An area column named `cell_area`, `Area`, or another name containing `area` is preferred.

The files in `config/` record project paths, plotting defaults, biological
samples, the six multi-sample slide directories, and the approved resolution-5
whole-tissue production source used by downstream VZ/RL work. Not every
analysis script reads configuration yet. From the project root, validate its
structure and referenced input files with:

```bash
Rscript scripts/validate_config.R --check-files
```

## Main analysis stages

### Atlas coverage/QC reporting and manuscript panels

See [ATLAS_COVERAGE_QC.md](ATLAS_COVERAGE_QC.md) for the additive, read-only
atlas reporting layer, companion donor/section/anatomy templates, and figure
exports. `xenium_atlas_qc.R` supports manifest-derived `--list`, path-only
`--dry-run`, per-sample `--task`, and small-table `--compile`.
`xenium_atlas_plot_qc.R` exports coverage, retention, QC and donor-support
panels when evidence exists, explicitly recording blocked panels otherwise.
The 41 manifest entries correspond to 12 working cerebellum groups under
user-supplied grouping, including two provisional FB/GZFB merges. These counts
must not be substituted for independently verified donor replication.
Compilation/plotting each take approximately 1–5 minutes. HPC extraction is
estimated at 2–10 minutes/task, 28–140 minutes for 41 tasks at concurrency three,
excluding queue wait (low confidence, saved-object/raw-table I/O). No production
filtering, identities, objects or stable output names are changed.

### Downstream resource updates 2026-10-17

The current manifest has 41 samples. The implemented resource and coverage
updates are documented in [DOWNSTREAM_RESOURCE_UPDATES_20261017.md](DOWNSTREAM_RESOURCE_UPDATES_20261017.md),
including expected runtimes, requested limits, validation results, and HPC
submission examples. This implementation supersedes resource recommendations
in the dated [review](DOWNSTREAM_RESOURCE_REVIEW_20261016.md).

Clustering, stage-6 markers, and per-sample mapping default to one R worker.
Whole-cohort combined/spatial launchers provisionally request 600 GB; these
requests still need measured peak-memory evidence. Scientific parameters,
identities, seeds, output names, overwrite protections, and the intentional
nine-sample stage-9 cohort are preserved. Array submission should use
`submit_xenium_downstream_array.sh` to derive task bounds from the R driver's
manifest map (currently 41 samples; concurrency three).

### New samples added 2026-10-14

The manifest now contains 41 analyzed samples (25 single-slide and 16 split-slide
components). The original 34 row/task IDs are preserved. See
[`NEW_SAMPLES_20261014.md`](NEW_SAMPLES_20261014.md) for the targeted seven-sample
commands; historical 34-sample run examples below are not submission instructions
for this cohort. Only section 2 is included from `GZFB124_1_X_G_3__2__1`.

For cohort-wide report coverage, see
[`MERGED_OUTPUTS_20261015.md`](MERGED_OUTPUTS_20261015.md). Both consensus PDFs
and the Aldinger combined PDF require 41 pages. The dedicated report-refresh
helper archives old PDFs/caches before rebuilding them; its path-only preflight
is estimated at 1–3 minutes. Consensus rendering is estimated at 3–15 minutes
per task, 21–105 minutes for 41 tasks at concurrency six, plus a 1–5-minute
merge, excluding queue wait (low confidence). It does not repeat label transfer.

### 1. Xenium gene panel

For an optional, independent sample sex-marker check after whole-tissue QC,
see [SEX_MARKER_EXPRESSION.md](SEX_MARKER_EXPRESSION.md). It summarizes SRY,
DDX3Y, RPS4Y1, KDM5D, UTY, EIF1AY, USP9Y, and XIST in a manifest-derived
Slurm array, reports marker coverage, and produces provisional
male-like/female-like/unresolved evidence without updating sample metadata.
Expected runtime is 1–5 minutes per sample, approximately 15–70 minutes for
41 samples at concurrency three plus 1–3 minutes to combine, excluding queue
wait (low confidence; based on saved-object I/O and static inspection).

`xenium_extract_gene_panel.R` extracts target genes from a Xenium panel JSON file. The resulting gene list is used to restrict the reference datasets to genes measurable by Xenium.

### 2. Reference preparation

The reference workflows harmonize cell-type labels, prepare Seurat objects, restrict them to the Xenium panel, and recompute dimensional reductions.

| Reference | Initial preparation | Xenium-panel subset | Additional analysis |
|---|---|---|---|
| Aldinger | `aldinger_01_update_seurat_object.R`, `aldinger_02_standardize_cell_types.R` | `aldinger_03_subset_gene_panel.R` | `aldinger_03b_plot_summary.R`, `aldinger_04_subset_vz.R` |
| Sepp | `sepp_01_standardize_cell_types.R` | `sepp_02_subset_gene_panel.R` | `sepp_03_subset_vz.R` |
| Science | `science_01_standardize_cell_types.R` | `science_02_subset_gene_panel.R` | — |

The corresponding `run_*_subset_gene_panel.slurm` files submit the panel-subsetting jobs.

Reference objects and figures use study-specific directories beneath
`outputs/references/`: `aldinger/`, `sepp/`, and `science/` each contain
`rds/` and `plots/`, while reference-specific VZ objects use `vz/rds/` and
`vz/plots/`. Cross-study figures are written to
`outputs/references/cross_study/plots/`. The updated Aldinger object created by
`aldinger_01_update_seurat_object.R` is stored with the other Aldinger project
outputs; the published source object under `/data/` remains unchanged.

### 3. Xenium import, cerebellum cropping, and QC

The two manifest-driven preprocessing drivers source the same ordered steps:

1. `xenium_preprocess_01_crop_cerebellum.R`
2. `xenium_preprocess_02_qc_cells.R`
3. `xenium_preprocess_03_normalize_cluster.R`

`xenium_preprocess_01_crop_cerebellum.R`:

- loads a sample with `Seurat::LoadXenium()`;
- retains cells listed in `cerebellum_cells_stats.csv`;
- adds cell-area metadata;
- saves a cropped spatial feature plot; and
- writes `<sample>_CB.rds` to
  `outputs/xenium/preprocess/01_cropped/rds/`.

For combined-slide data, `xenium_preprocess_01_crop_cerebellum.R` accepts the shared input
directory and sample-specific cell-stat CSV separately. The manifest-driven
`xenium_preprocess_split_slides.R` driver uses this interface automatically.
The earlier hard-coded combined-slide crop and driver scripts are archived in
`scripts/OLD/`.

`xenium_preprocess_02_qc_cells.R` calculates:

- `nCount_Xenium`;
- `nFeature_Xenium`;
- cell area; and
- the percentage of BlankCodeword, ControlCodeword, ControlProbe, and GenomicControl counts.

The default filters are:

| Metric | Lower threshold | Upper threshold |
|---|---:|---:|
| Xenium transcripts | max(20, median − 3 MAD) | 99th percentile |
| detected Xenium genes | max(10, median − 3 MAD) | 99th percentile |
| cell area | 15 | 99th percentile |
| control percentage | — | 5% |

QC thresholds are stored in `object@misc$QC_thresholds`. Reports are written
to `outputs/xenium/preprocess/02_qc/reports/`, and QC-filtered objects are
written to `outputs/xenium/preprocess/02_qc/rds/`.

There are two active manifest-driven drivers based on input layout:

- `xenium_preprocess_single_slides.R` – 19 biological samples from slides
  containing one sample; and
- `xenium_preprocess_split_slides.R` – 15 manually separated biological
  samples from six multi-sample slides.

Both drivers support `--list`, path-only `--dry-run`, protected `--overwrite`,
and one sample per Slurm array task. Their task maps come directly from
`config/samples.csv`.

### 4. Initial normalization and clustering

`xenium_preprocess_03_normalize_cluster.R` defines `process_xenium_clusters()`. It performs:

- log normalization using the median transcript count as the scale factor;
- 2,000 variable features;
- scaling and 50-component PCA;
- UMAP and Annoy nearest neighbors using PCs 1–50;
- Louvain clustering (`algorithm = 1`) at resolution 1.5; and
- UMAP, global spatial, and faceted spatial cluster plots.

Processed objects are saved as
`outputs/xenium/preprocess/03_clustered/rds/<sample>_CB_QC_cluster.rds`, with
plots in `outputs/xenium/preprocess/03_clustered/plots/`.

#### Resolution-2.0 whole-tissue pilot

Before changing the production clustering resolution, use
`xenium_preprocess_03b_resolution2_pilot.R` to compare resolution 2.0 with the
existing resolution-1.5 clusters on the same saved SNN graph. The three pilot
samples are `FB328_1_X_G` (PCW12), `GZFB_9_X_G_2` (PCW11), and
`GZFB_22_X_G_4` (PCW19). The script preserves cell IDs and all existing
metadata, validates the baseline cluster column, and does not repeat QC,
normalization, PCA, UMAP, or neighbor finding.

Inspect the task map and a path-only task before submission:

```bash
Rscript scripts/xenium_preprocess_03b_resolution2_pilot.R --list
Rscript scripts/xenium_preprocess_03b_resolution2_pilot.R --dry-run 1
```

Submit all three tasks with:

```bash
sbatch scripts/run_xenium_preprocess_03b_resolution2_pilot.slurm
```

Pilot RDS objects, side-by-side UMAP and spatial plots, faceted resolution-2.0
spatial plots, cluster counts, cluster-transition tables, and provenance are
written beneath `outputs/xenium/preprocess/03b_resolution2_pilot/`. Existing
production resolution-1.5 outputs are never replaced. Pilot outputs are also
protected; use the launcher's `RES2_PILOT_OVERWRITE=true` switch only after
reviewing an intentional replacement.

Expected runtime is approximately 15–60 minutes per task and 15–60 minutes
total at the configured concurrency of three, excluding queue wait. This is a
low-confidence estimate based on static inspection because no matching `sacct`
history was found; timing one task will improve it. The requested Slurm time
limit is 4 hours per task, safely above the expected runtime.

#### Label transfer for the resolution-2.0 pilot

After all three pilot RDS files have been reviewed for completeness,
`xenium_annotate_01_label_transfer_rpca.R --pilot-res2` runs the unchanged
RPCA label-transfer settings against Aldinger, Sepp, or Science. Pilot mode
uses the same three-task mapping from
`config/resolution2_pilot_samples.csv` and refuses to read the production
resolution-1.5 inputs.

Inspect the pilot mapping and all nine reference/sample input combinations
before submission. These dry-runs inspect paths only and do not load RDS files:

```bash
Rscript scripts/xenium_annotate_01_label_transfer_rpca.R --pilot-res2 --list
for reference in Aldinger Sepp Science; do
  for task_id in 1 2 3; do
    Rscript scripts/xenium_annotate_01_label_transfer_rpca.R \
      --pilot-res2 --dry-run "${reference}" "${task_id}"
  done
done
```

Submit one three-sample array for each reference:

```bash
sbatch --job-name=Xen_res2_Aldinger \
  --export=ALL,REFERENCE=Aldinger \
  scripts/run_xenium_annotate_01_transfer_resolution2_pilot.slurm
sbatch --job-name=Xen_res2_Sepp \
  --export=ALL,REFERENCE=Sepp \
  scripts/run_xenium_annotate_01_transfer_resolution2_pilot.slurm
sbatch --job-name=Xen_res2_Science \
  --export=ALL,REFERENCE=Science \
  scripts/run_xenium_annotate_01_transfer_resolution2_pilot.slurm
```

Each reference writes its annotated RDS objects, majority/weighted comparison
tables, prediction-count tables, score histograms, UMAPs, global and faceted
spatial plots, and broad-marker DotPlots beneath
`outputs/xenium/preprocess/03b_resolution2_pilot/annotation/01_label_transfer/<reference>/`.
No production annotation output is replaced. Pilot annotation outputs are
protected; set `PILOT_ABT_OVERWRITE=true` only after reviewing an intentional
replacement.

Observed runtime for the nine pilot label-transfer tasks was approximately
8–22 minutes per task. Each three-sample reference array completed in about
21–22 minutes, and all three arrays completed in about 22 minutes when they ran
concurrently, excluding queue wait. Peak memory was approximately 9–25 GB per
task. These measurements are from Slurm jobs 39831999, 39832000, and 39832014;
the requested time limit remains 4 hours per task.

#### Consensus labels for the resolution-2.0 pilot

Build the three cluster-level consensus tables only after all nine label-transfer
tables are present. Pilot mode requires Aldinger, Sepp, and Science results for
all three manifest samples, validates one-to-one cluster joins and identical
cluster sets across references, and refuses to overwrite an existing consensus
table.

```bash
Rscript scripts/xenium_annotate_02_build_consensus.R --pilot-res2 --list
Rscript scripts/xenium_annotate_02_build_consensus.R --pilot-res2 --dry-run
Rscript scripts/xenium_annotate_02_build_consensus.R --pilot-res2
```

The dry-run inspects paths only and should finish in under 10 seconds. The CSV
consensus merge should finish in under 1 minute. It writes the three tables to
`outputs/xenium/preprocess/03b_resolution2_pilot/annotation/02_consensus/tables/`.
If reviewed tables intentionally need replacement, rerun the final command with
`--overwrite`.

Inspect the three-task apply mapping and path-only inputs, then submit the three
samples as one array:

```bash
Rscript scripts/xenium_annotate_03_apply_consensus.R --pilot-res2 --list
for task_id in 1 2 3; do
  Rscript scripts/xenium_annotate_03_apply_consensus.R \
    --pilot-res2 --dry-run "${task_id}"
done

sbatch scripts/run_xenium_annotate_03_apply_consensus_resolution2_pilot.slurm
```

The apply stage reads the pilot Aldinger object as the expression/spatial
container, verifies that `seurat_clusters` still exactly matches the
resolution-2.0 clusters, adds all reference labels and `consensus_label`, and
writes annotated RDS objects, UMAPs, global and faceted spatial maps, and marker
DotPlots beneath
`outputs/xenium/preprocess/03b_resolution2_pilot/annotation/03_consensus_labels/`.
It never writes into the production annotation tree. Existing pilot consensus
outputs are protected; set `PILOT_CONSENSUS_OVERWRITE=true` only after reviewing
an intentional replacement.

Expected runtime is approximately 3–15 minutes per task and 3–15 minutes total
at the configured concurrency of three, excluding queue wait. This is a
moderate-confidence estimate based on the 8–22 minute, 9–25 GB pilot
label-transfer history and the apply stage's lighter static workload. The
launcher requests 1 hour, 1 CPU, and 64 GB per task.

#### Resolution-3.0 pilot for the same three samples

Resolution 3.0 uses the same three-sample pilot manifest but writes to a new
`outputs/xenium/preprocess/03c_resolution3_pilot/` tree. It reclusters the
saved production SNN graph and preserves the resolution-1.5 and resolution-2.0
results. Faceted spatial plots made in resolution-3 mode use point size 0.03;
global spatial plots and earlier pilot figures are unchanged.

Inspect and submit clustering:

```bash
Rscript scripts/xenium_preprocess_03b_resolution2_pilot.R --resolution3 --list
for task_id in 1 2 3; do
  Rscript scripts/xenium_preprocess_03b_resolution2_pilot.R \
    --resolution3 --dry-run "${task_id}"
done
sbatch scripts/run_xenium_preprocess_03c_resolution3_pilot.slurm
```

After all three clustering tasks complete, inspect and submit all three
references:

```bash
Rscript scripts/xenium_annotate_01_label_transfer_rpca.R --pilot-res3 --list
for reference in Aldinger Sepp Science; do
  for task_id in 1 2 3; do
    Rscript scripts/xenium_annotate_01_label_transfer_rpca.R \
      --pilot-res3 --dry-run "${reference}" "${task_id}"
  done
done

for reference in Aldinger Sepp Science; do
  sbatch --job-name="Xen_res3_${reference}" \
    --export=ALL,REFERENCE="${reference}" \
    scripts/run_xenium_annotate_01_transfer_resolution3_pilot.slurm
done
```

After all nine transfers complete, build and apply consensus labels:

```bash
Rscript scripts/xenium_annotate_02_build_consensus.R --pilot-res3 --dry-run
Rscript scripts/xenium_annotate_02_build_consensus.R --pilot-res3

for task_id in 1 2 3; do
  Rscript scripts/xenium_annotate_03_apply_consensus.R \
    --pilot-res3 --dry-run "${task_id}"
done
sbatch scripts/run_xenium_annotate_03_apply_consensus_resolution3_pilot.slurm
```

To apply the weighted-reference 2-of-3 method to completed resolution-3 label
transfers without replacing the legacy consensus, run:

```bash
Rscript scripts/xenium_annotate_02_build_consensus.R \
  --pilot-res3 --weighted-2of3 --dry-run
Rscript scripts/xenium_annotate_02_build_consensus.R \
  --pilot-res3 --weighted-2of3

for task_id in 1 2 3; do
  Rscript scripts/xenium_annotate_03_apply_consensus.R \
    --pilot-res3 --weighted-2of3 --dry-run "${task_id}"
done
sbatch \
  scripts/run_xenium_annotate_03_apply_consensus_resolution3_weighted2of3_pilot.slurm
```

New tables and objects are written beneath
`02_consensus_weighted_2of3/` and `03_consensus_labels_weighted_2of3/` inside
the resolution-3 annotation tree.

To test whether resolution-3 weighted-2-of-3 `Unknown` clusters have relatively
low QC among the cells that already passed the original QC filters, run:

```bash
Rscript scripts/xenium_annotate_04_audit_unknown_qc.R \
  --pilot-res3 --weighted-2of3 --list

for task_id in 1 2 3; do
  Rscript scripts/xenium_annotate_04_audit_unknown_qc.R \
    --pilot-res3 --weighted-2of3 --dry-run "${task_id}"
done

sbatch \
  scripts/run_xenium_annotate_04_audit_unknown_qc_resolution3_weighted2of3_pilot.slurm
```

The audit is descriptive and does not filter or relabel cells. It writes
within-sample Unknown-versus-Known QC summaries, rank-biserial effect sizes,
Known-cell tail fractions for each Unknown cluster, stored QC thresholds, and a
two-page PDF beneath `03_consensus_labels_weighted_2of3/unknown_qc_audit/`.
Cell-level P values are not treated as independent biological replication.

For completed resolution-4 weighted-2-of-3 consensus objects, use the same
audit driver with `--pilot-res4`, or submit:

```bash
sbatch \
  scripts/run_xenium_annotate_04_audit_unknown_qc_resolution4_weighted2of3_pilot.slurm
```

Expected clustering runtime is approximately 3–15 minutes per task and total
at concurrency three (moderate confidence from the same saved SNN graphs and
sample sizes). Label transfer is expected to remain about 8–25 minutes per task
and approximately 20–30 minutes per three-sample reference array at concurrency
two, based on the measured resolution-2 jobs; all three reference arrays may
finish in about 20–30 minutes if scheduled together. The CSV merge should take
under 1 minute, and consensus application approximately 3–15 minutes per task
and total at concurrency three. Estimates exclude queue wait. Each resolution-3
launcher requests a 1-hour time limit.

#### Resolution-4.0 pilot for the same three samples

Resolution 4.0 reuses the three-sample pilot manifest and writes only beneath
`outputs/xenium/preprocess/03d_resolution4_pilot/`. It preserves production,
resolution-2, and resolution-3 results. As in resolution-3 mode, all faceted
spatial plots use point size 0.03 while global spatial plots are unchanged.

Inspect and submit clustering:

```bash
Rscript scripts/xenium_preprocess_03b_resolution2_pilot.R --resolution4 --list
for task_id in 1 2 3; do
  Rscript scripts/xenium_preprocess_03b_resolution2_pilot.R \
    --resolution4 --dry-run "${task_id}"
done
sbatch scripts/run_xenium_preprocess_03d_resolution4_pilot.slurm
```

After clustering completes, inspect and submit the three label transfers:

```bash
for reference in Aldinger Sepp Science; do
  for task_id in 1 2 3; do
    Rscript scripts/xenium_annotate_01_label_transfer_rpca.R \
      --pilot-res4 --dry-run "${reference}" "${task_id}"
  done
done

for reference in Aldinger Sepp Science; do
  sbatch --job-name="Xen_res4_${reference}" \
    --export=ALL,REFERENCE="${reference}" \
    scripts/run_xenium_annotate_01_transfer_resolution4_pilot.slurm
done
```

After all nine transfers complete, build and apply consensus labels:

```bash
Rscript scripts/xenium_annotate_02_build_consensus.R \
  --pilot-res4 --weighted-2of3 --dry-run
Rscript scripts/xenium_annotate_02_build_consensus.R \
  --pilot-res4 --weighted-2of3

for task_id in 1 2 3; do
  Rscript scripts/xenium_annotate_03_apply_consensus.R \
    --pilot-res4 --weighted-2of3 --dry-run "${task_id}"
done
sbatch \
  scripts/run_xenium_annotate_03_apply_consensus_resolution4_weighted2of3_pilot.slurm
```

The weighted 2-of-3 method gives each reference one vote using its weighted
cluster label. Majority/weighted agreement is retained for auditing but is not
counted as a second vote. A final non-Unknown label requires agreement from at
least two references. Its tables and annotated outputs are isolated beneath
`02_consensus_weighted_2of3/` and `03_consensus_labels_weighted_2of3/`, so the
legacy six-vote consensus remains available for comparison.

Based on the measured resolution-3 pilot, clustering should take approximately
4–8 minutes per task and total at concurrency three. Label transfer should take
approximately 9–25 minutes per task and per reference array at concurrency
three; all three arrays may finish in about 9–25 minutes if all nine tasks are
scheduled together. The CSV merge should take under 1 minute, and consensus application
approximately 1–3 minutes per task and total at concurrency three. Estimates
exclude queue wait. Clustering requests 1 hour, 1 CPU, and 32 GB; label transfer
requests 1 hour, 1 CPU, and 64 GB; consensus requests 30 minutes, 1 CPU, and
32 GB. Three simultaneously submitted reference arrays can reserve up to 576 GB
cluster-wide at their maximum combined concurrency of nine.

#### Resolution-5.0 pilot for the same three samples

Resolution 5.0 writes only beneath
`outputs/xenium/preprocess/03e_resolution5_pilot/` and preserves all lower-resolution
results. It reuses the saved SNN graph, the three approved references, the
weighted 2-of-3 consensus method, and faceted spatial point size 0.01.

Run clustering:

```bash
Rscript scripts/xenium_preprocess_03b_resolution2_pilot.R --resolution5 --list
for task_id in 1 2 3; do
  Rscript scripts/xenium_preprocess_03b_resolution2_pilot.R \
    --resolution5 --dry-run "${task_id}"
done
sbatch scripts/run_xenium_preprocess_03e_resolution5_pilot.slurm
```

After all three clustering tasks complete, run the nine label transfers:

```bash
for reference in Aldinger Sepp Science; do
  for task_id in 1 2 3; do
    Rscript scripts/xenium_annotate_01_label_transfer_rpca.R \
      --pilot-res5 --dry-run "${reference}" "${task_id}"
  done
done

for reference in Aldinger Sepp Science; do
  sbatch --job-name="Xen_res5_${reference}" \
    --export=ALL,REFERENCE="${reference}" \
    scripts/run_xenium_annotate_01_transfer_resolution5_pilot.slurm
done
```

After all nine transfers complete, build and apply weighted consensus labels:

```bash
Rscript scripts/xenium_annotate_02_build_consensus.R \
  --pilot-res5 --weighted-2of3 --dry-run
Rscript scripts/xenium_annotate_02_build_consensus.R \
  --pilot-res5 --weighted-2of3

for task_id in 1 2 3; do
  Rscript scripts/xenium_annotate_03_apply_consensus.R \
    --pilot-res5 --weighted-2of3 --dry-run "${task_id}"
done
sbatch \
  scripts/run_xenium_annotate_03_apply_consensus_resolution5_weighted2of3_pilot.slurm
```

The optional Unknown-cluster QC audit is:

```bash
sbatch \
  scripts/run_xenium_annotate_04_audit_unknown_qc_resolution5_weighted2of3_pilot.slurm
```

Expected runtimes and requested resources are the same as the resolution-4
pilot estimates above. Resolution 5 may create smaller clusters; compare cluster
sizes, marker coherence, spatial coherence, and Unknown frequency with the
resolution-3 and resolution-4 pilots before selecting it for all 34 samples.

#### Resolution-5.0 weighted 2-of-3 pilot for one selected sample

Use this isolated mode to inspect one manifest-validated sample without rerunning
the original three-sample pilot. Set the sample explicitly; the example below is
`GZFB_22_X_G_5`. All results are protected beneath
`outputs/xenium/preprocess/03h_resolution5_selected_sample/`.

First verify the mapping and submit clustering:

```bash
export SELECTED_SAMPLE_ID=GZFB_22_X_G_5

Rscript scripts/xenium_preprocess_03b_resolution2_pilot.R \
  --resolution5 --selected-sample --list
Rscript scripts/xenium_preprocess_03b_resolution2_pilot.R \
  --resolution5 --selected-sample --dry-run 1

cluster_job=$(sbatch --parsable \
  --export=ALL,SELECTED_SAMPLE_ID="${SELECTED_SAMPLE_ID}" \
  scripts/run_xenium_preprocess_03h_resolution5_selected_sample.slurm)
echo "${cluster_job}"
```

After clustering completes successfully, verify and submit all three reference
transfers. The three jobs may run concurrently:

```bash
for reference in Aldinger Sepp Science; do
  Rscript scripts/xenium_annotate_01_label_transfer_rpca.R \
    --pilot-res5 --selected-sample --dry-run "${reference}" 1
done

for reference in Aldinger Sepp Science; do
  sbatch --job-name="Xen_res5_${reference}_${SELECTED_SAMPLE_ID}" \
    --export=ALL,SELECTED_SAMPLE_ID="${SELECTED_SAMPLE_ID}",REFERENCE="${reference}" \
    scripts/run_xenium_annotate_01_transfer_resolution5_selected_sample.slurm
done
```

After all three transfers complete, build the weighted 2-of-3 consensus table
and submit the final annotation and plot job:

```bash
Rscript scripts/xenium_annotate_02_build_consensus.R \
  --pilot-res5 --selected-sample --weighted-2of3 --dry-run
Rscript scripts/xenium_annotate_02_build_consensus.R \
  --pilot-res5 --selected-sample --weighted-2of3

Rscript scripts/xenium_annotate_03_apply_consensus.R \
  --pilot-res5 --selected-sample --weighted-2of3 --dry-run 1
final_job=$(sbatch --parsable \
  --export=ALL,SELECTED_SAMPLE_ID="${SELECTED_SAMPLE_ID}" \
  scripts/run_xenium_annotate_03_apply_consensus_resolution5_selected_sample_weighted2of3.slurm)
echo "${final_job}"
```

The requested faceted plot is written to
`outputs/xenium/preprocess/03h_resolution5_selected_sample/annotation/03_consensus_labels_weighted_2of3/plots/GZFB_22_X_G_5_Consensus_FacetSpatial.tif`.
Its point size is 0.01.

Expected active runtimes, excluding queue wait, are approximately 5–15 minutes
for clustering; 15 minutes to 3 hours per label-transfer job and about the same
total wall time when all three run concurrently; under 2 minutes for consensus
table construction; and 10–60 minutes for final annotation and plotting. These
moderate-confidence estimates use the measured resolution-3/4 pilot jobs and the
observed 91,577-cell label-transfer outlier. The requested Slurm limits are 1,
4, and 2 hours for clustering, each transfer, and final plotting, respectively.

#### Consensus-label consistency across resolutions 3, 4, and 5

After weighted-2-of-3 consensus objects exist for all three resolutions, the
resolution-consistency audit validates identical cell sets and compares labels
and cluster membership without changing any Seurat object. It writes per-cell
agreement status, pairwise label agreement and Cohen's kappa, cluster ARI/NMI,
per-label Jaccard/retention, label and cluster transition tables, UMAP/spatial
agreement plots, and PDF reports beneath
`outputs/validation/resolution_selection/weighted_2of3/`.

Inspect and submit the three-sample array:

```bash
Rscript scripts/xenium_annotate_05_compare_consensus_across_resolutions.R --list
for task_id in 1 2 3; do
  Rscript scripts/xenium_annotate_05_compare_consensus_across_resolutions.R \
    --dry-run "${task_id}"
done
sbatch scripts/run_xenium_annotate_05_resolution_consistency_pilot.slurm
```

After all three tasks complete, combine their summaries:

```bash
Rscript scripts/xenium_annotate_05_compare_consensus_across_resolutions.R \
  --combine --dry-run
Rscript scripts/xenium_annotate_05_compare_consensus_across_resolutions.R \
  --combine
```

The `cross_resolution_modal_label` column is descriptive and is never written
back to a Seurat object. Treat agreement across resolutions as a sensitivity
analysis, not independent biological replication.

#### Resolution-4 analysis for all 34 samples

The approved all-sample resolution-4 workflow reuses each saved production SNN
graph. Resolution-4 clustering objects, plots, and tables are written beneath
`outputs/xenium/preprocess/03g_resolution4_all_samples/`; label transfer and
consensus outputs are written beneath
`outputs/xenium/annotation/resolution4_all_samples/`. It does not overwrite
resolution-1.5 production objects or any pilot. All broad-label DotPlots use
the same scaled-expression range (-2.5 to 2.5), percent-expressed range
(0% to 100%), cell-type order, genes-on-x orientation, and vertical right-side
legends. Faceted spatial point size is 0.01.

Inspect and submit clustering:

```bash
Rscript scripts/xenium_preprocess_03b_resolution2_pilot.R \
  --all-samples-res4 --list

for task_id in $(seq 1 34); do
  Rscript scripts/xenium_preprocess_03b_resolution2_pilot.R \
    --all-samples-res4 --dry-run "${task_id}"
done

sbatch scripts/run_xenium_preprocess_03g_resolution4_all_samples.slurm
```

After all 34 clustering tasks complete, inspect and submit all three reference
transfers:

```bash
for reference in Aldinger Sepp Science; do
  for task_id in $(seq 1 34); do
    Rscript scripts/xenium_annotate_01_label_transfer_rpca.R \
      --all-samples-res4 --dry-run "${reference}" "${task_id}"
  done
done

for reference in Aldinger Sepp Science; do
  sbatch --job-name="Xen_res4_all_${reference}" \
    --export=ALL,REFERENCE="${reference}" \
    scripts/run_xenium_annotate_01_transfer_resolution4_all_samples.slurm
done
```

After all 102 transfer tasks complete, build weighted 2-of-3 consensus tables:

```bash
Rscript scripts/xenium_annotate_02_build_consensus.R \
  --all-samples-res4 --weighted-2of3 --dry-run
Rscript scripts/xenium_annotate_02_build_consensus.R \
  --all-samples-res4 --weighted-2of3
```

Inspect and submit consensus application and plots:

```bash
for task_id in $(seq 1 34); do
  Rscript scripts/xenium_annotate_03_apply_consensus.R \
    --all-samples-res4 --weighted-2of3 --dry-run "${task_id}"
done

sbatch \
  scripts/run_xenium_annotate_03_apply_consensus_resolution4_all_samples_weighted2of3.slurm
```

To regenerate all final consensus plots after a visualization-only change,
without rewriting any consensus RDS or table, inspect and submit the dedicated
plots-only array:

```bash
for task_id in $(seq 1 34); do
  Rscript scripts/xenium_annotate_03_apply_consensus.R \
    --all-samples-res4 --weighted-2of3 --plots-only --dry-run "${task_id}"
done
sbatch \
  scripts/run_xenium_annotate_03_regenerate_plots_resolution4_all_samples_weighted2of3.slurm
```

The plots-only mode intentionally replaces the five plot files for each sample
but reads the existing consensus object as its sole analysis input and never
rewrites RDS or consensus tables. Expected runtime is approximately 10-60
minutes per task and 1.5-8 hours total at concurrency four (low confidence,
based on the measured plotting time for the largest sample); each task requests
2 hours.

After all 34 consensus objects and per-sample plots complete, inspect and submit
the PCW-ordered merged report:

```bash
Rscript scripts/xenium_annotate_03d_plot_report.R \
  --all-samples-res4 --weighted-2of3 --list
Rscript scripts/xenium_annotate_03d_plot_report.R \
  --all-samples-res4 --weighted-2of3 --dry-run
page_job=$(sbatch --parsable \
  scripts/run_xenium_annotate_03d_plot_report_resolution4_all_samples_weighted2of3.slurm)
merge_job=$(sbatch --parsable --dependency="afterok:${page_job}" \
  scripts/run_xenium_annotate_03d_plot_report_resolution4_all_samples_weighted2of3_merge.slurm)
echo "Page array: ${page_job}; dependent merge: ${merge_job}"
```

The protected 34-page output is
`outputs/xenium/annotation/resolution4_all_samples/03_consensus_labels_weighted_2of3/plots/Xenium_Consensus_Res4_Weighted2of3_All_Samples_Report.pdf`.
Each page validates resolution-4 identities and `weighted_2of3` provenance,
uses the shared fixed DotPlot scales and orientation, and uses point size 0.01
in the faceted spatial panel. Pages are rendered independently at concurrency
six, then a dependent single job verifies and merges the 34 ordered PNG pages.
Expected page-rendering runtime is approximately 10-45 minutes per task and
1-4 hours total at concurrency six; the merge should take under 10 minutes.
These are low-confidence estimates based on static inspection and the measured
cost of plotting the largest sample. Page tasks request 2 hours each and the
merge requests 30 minutes.

Expected clustering runtime is approximately 4-10 minutes per task and
35-90 minutes total at concurrency four. Each reference transfer is expected
to take approximately 8-30 minutes per task and 1.5-6 hours per 34-task array
at concurrency three; the three arrays can run simultaneously if cluster
capacity permits. Consensus-table construction should take under 2 minutes,
and consensus application approximately 1-4 minutes per task and 10-25 minutes
total at concurrency six. Estimates exclude queue wait and are based on the
measured three-sample pilots; estimates for unmeasured larger samples are
moderate confidence.

#### Resolution-5 analysis for all 34 samples, start to finish

The current production workflow uses the complete Aldinger, Sepp, and Science
references. Each reference contributes one cluster-level weighted vote only
when its winning mean score is at least 0.4 and exceeds the runner-up mean score
by at least 0.05. Consensus remains an equal, unique two-of-three vote; a
cluster without that support becomes `Unknown-##`. Seed 42, RPCA dimensions
1:30, `k.anchor = 20`, `k.score = 30`, `k.weight = 50`, marker sets, and fixed
DotPlot scales are unchanged. These approved production settings are recorded
under `label_transfer.resolution5_all_samples` in `config/config.yml`.
Clustering outputs are written beneath
`outputs/xenium/preprocess/04_resolution5_clustered/`; annotation and report
outputs are written directly beneath `outputs/xenium/annotation/` because
resolution 5 is the only active whole-tissue annotation workflow.
Existing resolution-1.5, resolution-4, and pilot results are not replaced.
The `production` section of `config/config.yml` designates these resolution-5
weighted two-of-three consensus objects as the sole whole-tissue source for
downstream spatial plotting, VZ extraction, RL extraction, and output
validation.

The apparently resolution-neutral directories `01_cropped`, `02_qc`, and
`03_clustered`, plus the current Aldinger, Sepp, and Science reference objects,
remain required to reproduce the resolution-5 outputs. They must not be
removed as part of legacy-resolution cleanup.

For installations that still use the historical
`03i_resolution5_all_samples` directory name, validate and perform the one-time
rename before running any current production or cleanup script:

```bash
cd /home/acflint/R/Projects/XeniumFCProject
bash scripts/rename_resolution5_preprocess_output.sh --dry-run
bash scripts/rename_resolution5_preprocess_output.sh --move
```

The utility checks every sample in `config/samples.csv`, requires exactly 34
nonempty resolution-5 clustering RDS files, refuses to merge with an existing
target, and then renames the directory to `04_resolution5_clustered`. The move
does not rewrite any RDS, table, or plot. On the same filesystem it should take
seconds; this estimate excludes time spent checking directory entries and is
high confidence for a normal metadata-only rename.

To replace the existing resolution-5 annotation while reusing the completed
34 clustering objects, first run the path-only preflight:

```bash
cd /home/acflint/R/Projects/XeniumFCProject
bash \
  scripts/submit_xenium_resolution5_all_samples_full_reference_margin005.sh \
  --dry-run
```

After reviewing the 34-row mapping and existing output counts, explicitly
authorize replacement of the resolution-5 transfer, consensus, per-sample
plot, report-page, and merged-report outputs:

```bash
bash \
  scripts/submit_xenium_resolution5_all_samples_full_reference_margin005.sh \
  --overwrite
```

This annotation-only submission does not rewrite the resolution-5 clustering
RDS files. It submits three 34-task full-reference arrays followed by consensus
construction, consensus application, report pages, and report merge with
`afterok` dependencies. No production launcher overwrites by default; the
helper propagates overwrite permission only when invoked with `--overwrite`.
Each reference table now includes the pre-abstention winner, runner-up, mean
scores, top-two margin, failed-rule flags, sampling mode, transfer parameters,
and seed. Matching per-task reference-class balance and transfer-provenance
CSVs are written beside the established comparison and cell-count tables.

After reviewing the 34-row mapping, this single shell block submits the complete
workflow with `afterok` dependencies:

```bash
set -euo pipefail
cd /home/acflint/R/Projects/XeniumFCProject

Rscript scripts/xenium_preprocess_03b_resolution2_pilot.R \
  --all-samples-res5 --list

for task_id in $(seq 1 34); do
  Rscript scripts/xenium_preprocess_03b_resolution2_pilot.R \
    --all-samples-res5 --dry-run "${task_id}"
done

cluster_submit=$(sbatch --parsable \
  scripts/run_xenium_preprocess_04_resolution5_clustered.slurm)
cluster_job=${cluster_submit%%;*}

transfer_jobs=()
for reference in Aldinger Sepp Science; do
  transfer_submit=$(sbatch --parsable \
    --dependency="afterok:${cluster_job}" \
    --array="1-34%3" \
    --job-name="Xen_res5_all_${reference}" \
    --export=ALL,REFERENCE="${reference}" \
    scripts/run_xenium_annotate_01_transfer_resolution5_all_samples.slurm)
  transfer_jobs+=("${transfer_submit%%;*}")
done
transfer_dependency=$(IFS=:; echo "${transfer_jobs[*]}")

build_submit=$(sbatch --parsable \
  --dependency="afterok:${transfer_dependency}" \
  scripts/run_xenium_annotate_02_build_consensus_resolution5_all_samples_weighted2of3.slurm)
build_job=${build_submit%%;*}

apply_submit=$(sbatch --parsable \
  --dependency="afterok:${build_job}" \
  scripts/run_xenium_annotate_03_apply_consensus_resolution5_all_samples_weighted2of3.slurm)
apply_job=${apply_submit%%;*}

pages_submit=$(sbatch --parsable \
  --dependency="afterok:${apply_job}" \
  scripts/run_xenium_annotate_03d_plot_report_resolution5_all_samples_weighted2of3.slurm)
pages_job=${pages_submit%%;*}

merge_submit=$(sbatch --parsable \
  --dependency="afterok:${pages_job}" \
  scripts/run_xenium_annotate_03d_plot_report_resolution5_all_samples_weighted2of3_merge.slurm)
merge_job=${merge_submit%%;*}

printf 'clustering=%s\ntransfers=%s\nbuild=%s\napply=%s\npages=%s\nmerge=%s\n' \
  "${cluster_job}" "${transfer_dependency}" "${build_job}" \
  "${apply_job}" "${pages_job}" "${merge_job}"
```

Every production launcher repeats its relevant path-only `--dry-run` after its
upstream dependency succeeds. No launcher enables overwrite by default. The
all-cell final report is
`outputs/xenium/annotation/03_consensus_labels_weighted_2of3/plots/Xenium_Consensus_Res5_Weighted2of3_All_Samples_Report.pdf`.
A matched report with every `Unknown-##` cell removed is written as
`Xenium_Consensus_Res5_Weighted2of3_All_Samples_KnownOnly_Report.pdf` in the
same directory.
Report facet panels use point size 0.01, a one-pixel glyph, and alpha 0.55 to
make reduced overplotting visible after rasterization.

Validate the completed resolution-5 production chain with:

```bash
Rscript scripts/validate_outputs_through_consensus.R
```

The validator reads the 34 resolution-5 clustering objects, three annotated
objects per sample, final consensus objects, and their tables and plots. It
writes reports under `outputs/validation/resolution5_through_consensus/`.
Expect approximately 1-4 hours on HPC, excluding queue wait (low confidence,
based on static inspection and the volume of RDS data read sequentially).

After the resolution-5 replacement set is complete, inspect the exact legacy
cleanup set without moving anything:

```bash
bash scripts/cleanup_non_resolution5_outputs.sh --dry-run
```

The script requires all 34 resolution-5 clustering objects, all three complete
34-sample transfer sets and their provenance tables, all 34 consensus tables,
all 34 consensus objects, and the current three reference objects. Only after
that preflight passes can `--quarantine` move the enumerated resolution-1.5,
resolution-2/3/4, superseded `_old`, smoke-test, historical, and pilot paths to
`quarantine_non_resolution5_2026-09-29/`. It deliberately retains
`01_cropped`, `02_qc`, `03_clustered`, current references, the full 34-sample
resolution-5 tree, and the small cross-resolution comparison summary. A
separate confirmation token is required for permanent quarantine deletion.

Completed resolution-selection and method-development evidence is kept outside
the production trees. Preview or apply the one-time structural move with:

```bash
bash scripts/reorganize_resolution5_validation_outputs.sh --dry-run
bash scripts/reorganize_resolution5_validation_outputs.sh --move
```

The move places the cross-resolution comparison beneath
`outputs/validation/resolution_selection/` and the transfer-threshold,
reference-sampling, k-anchor, and transfer-diagnostic studies beneath
`outputs/validation/resolution5_method_development/`. Production preprocessing
then ends at `04_resolution5_clustered`, while the annotation directory
contains only numbered production stages 01 through 03.

After that method-development move, remove the redundant resolution-5 wrapper
without rewriting any analysis file:

```bash
bash scripts/flatten_resolution5_annotation.sh --dry-run
bash scripts/flatten_resolution5_annotation.sh --move
```

The utility requires all 34 transfer, consensus-table, and consensus-RDS
outputs, refuses collisions or unexpected wrapper contents, moves the three
production stage directories to `outputs/xenium/annotation/`, and records the
old-to-new mappings under `outputs/validation/`.

Final insufficient-support clusters are retained as distinct downstream
identities. They are numbered in numeric `seurat_clusters` order within each
sample as `Unknown-01`, `Unknown-02`, and so on; a sample with one such cluster
uses `Unknown-01`. Existing consensus UMAP, global spatial, faceted spatial,
and marker DotPlot filenames are the all-cell versions and show detailed
Unknown identities with deterministic grey shades. Each also has a matched
`_Consensus_KnownOnly_...` output made after removing all Unknown cells.

Every sample with at least one Unknown cluster also receives
`tables/unknown_markers/<sample>_Unknown_cluster_top20_markers.xlsx`, with one
sheet per `Unknown-##`. Markers are positive one-versus-rest Wilcoxon results
from the normalized `Xenium` assay (`only.pos=TRUE`, `logfc.threshold=0.1`,
`min.pct=0.01`), ranked by average log2 fold change, adjusted P value, and gene,
then limited to 20 rows. A sample with only `Unknown-01` therefore receives a
one-sheet workbook; only samples with no Unknown clusters receive no workbook.
During an explicit overwrite rerun, a pre-existing workbook for a sample that
now has no Unknown clusters is renamed with a
`.superseded_no_unknown_<timestamp>_<process-id>` suffix so it cannot be
mistaken for a current result and remains recoverable.
This conditional output requires `writexl` in the active project `renv`; the
production script never installs it.

The resolution-5 all-sample application and parallel report-page launchers
start R with `--max-ppsize=500000`, the R-supported maximum
pointer-protection-stack size, because marker calculations and report rendering
for larger samples can overflow R's default stack. This runtime safeguard does
not alter the marker comparison, report content, or outputs.

The saved Aldinger resolution-5 transfer figures can be assembled without
recomputing label transfer by running
`xenium_annotate_01b_build_aldinger_resolution5_plot_report.R`. It validates all
34 manifest samples and all nine unique TIFFs per sample, then produces one
3-by-3 landscape page per sample at
`01_label_transfer/aldinger/plots/Xenium_Aldinger_Res5_All_Samples_Plot_Report.pdf`.
The TIFF versions are used so the three PDF companions are not duplicated.
The report requires `tiff` in the project `renv`; the script and Slurm job never
install packages. Inspect paths with `--list` and `--dry-run`, then submit
`run_xenium_annotate_01b_build_aldinger_resolution5_plot_report.slurm`.

Expected runtime per clustering task is 4-10 minutes and approximately 35-90
minutes total at concurrency four. In full-reference sensitivity job 40474781,
the nine full-reference tasks completed in 9.1-12.8 minutes and used at most
24.65 GB, but those selected-sample tasks did not render the complete
production output set. Across the prior 170 capped production transfers,
observed times were 7.3-47.3 minutes and maximum RSS was 49.28 GB for manifest
task 6 (`GZFB4_X_G`). Combining those measurements, full-reference production
is estimated at 15-75 minutes per task and approximately 3-15 hours per
34-task reference array at concurrency three, excluding queue wait. This is a
moderate-confidence estimate. Each transfer requests 4 hours, 1 CPU, and 96 GB
on `short`; all three arrays can reserve at most 864 GB when their nine tasks
run simultaneously.
Consensus-table construction should take under 2 minutes. Consensus
application, paired plot sets, and conditional Unknown marker tests are
expected to take 15-90 minutes per task and approximately 1.5-9 hours total at
concurrency six. Paired report pages are expected to take 20-90 minutes per
task and approximately 2-9 hours total at concurrency six; merging both reports
should take under 15 minutes. These revised estimates are low confidence and
come from prior single-version plotting history plus static inspection of the
new doubled plotting and marker workload; the next `sacct` results should be
used to refine them. Requested limits remain 1 hour for clustering, 4 hours per
transfer, 30 minutes for table construction, 2 hours per
consensus/report-page task, and 30 minutes for the merge.

Because broad consensus labels determine VZ and RL cell inclusion, existing
VZ/RL subset, integration, QC, subclustering, and relabeling outputs become
stale after this production annotation replacement. Review the new broad-label
summary and Unknown clusters before regenerating either regional branch.

##### Compact audit of reviewed resolution-5 cluster mislabels

`config/reviewed_annotation_clusters.csv` records the five currently reviewed
cluster/expected-label pairs. The non-mutating audit
`xenium_annotate_00c_audit_reviewed_clusters.R` reads the existing Aldinger,
Sepp, and Science transfer RDS files for those samples and writes one row per
reviewed cluster/reference combination. The table reports reference-level
majority and weighted labels, mean RL/Granule/Glia/VZ scores, the top-two
margin, expected-versus-current score difference, cell-level prediction
fractions, and the fraction above the legacy 0.4 maximum-score threshold. It
never writes an RDS or modifies a label, threshold, or consensus table.

Inspect the reviewed rows and nine required RDS paths without loading an
analysis object:

```bash
Rscript scripts/xenium_annotate_00c_audit_reviewed_clusters.R --list
Rscript scripts/xenium_annotate_00c_audit_reviewed_clusters.R --dry-run
```

Submit the sequential audit on Cheaha with:

```bash
sbatch scripts/run_xenium_annotate_00c_audit_reviewed_clusters.slurm
```

The compact result is written to
`outputs/validation/resolution5_method_development/04_transfer_threshold_calibration_audit/reviewed_clusters/tables/Reviewed_cluster_reference_score_audit.csv`.
The path-only dry-run should take under 10 seconds. The full audit is expected
to take approximately 20-90 minutes because it sequentially loads nine
full-size Seurat objects; this is a low-confidence estimate based on RDS I/O
and the observed 7.3-47.3-minute complete transfer-task history. The launcher
requests 2 hours and 64 GB on `express`; collect `Elapsed` and `MaxRSS` with
`sacct` after the first run to refine both requests.

##### Focused reference-sampling and abstention sensitivity at resolution 5

`xenium_annotate_01c_reference_sampling_abstention_sensitivity.R` isolates two
questions without changing production annotations:

- Does using each complete reference instead of the production cap of 1,000
  cells per `clusters_refined` identity change the reference or consensus label?
- How many additional clusters become unresolved when a reference must meet a
  score-margin or cell-support requirement before casting its equal vote?

The three query samples are derived from the unique sample IDs in
`config/reviewed_annotation_clusters.csv`. Each is mapped to Aldinger, Sepp,
and Science under `cap1000` and `full`, producing 18 isolated tasks. Both
sampling modes retain the production RPCA settings: seed 42, 3,000 reference
variable features, dimensions 1:30, `k.anchor = 20`, `k.score = 30`,
`k.weight = 50`, and the 0.4 winning-score threshold. The only transfer-stage
change is the set of reference cells retained. Seed 42 is reset before
reference sampling, each PCA, and anchor finding so stochastic state does not
become a second difference between the two sampling modes.

Each transfer is evaluated with four prespecified sensitivity rules:

1. `production_score`: winning mean class score at least 0.4;
2. `score_margin_0.05`: production score plus a top-two mean-score margin of
   at least 0.05;
3. `score_margin_0.10`: production score plus a margin of at least 0.10; and
4. `score_margin_0.05_cell_support_0.60`: production score, margin at least
   0.05, and at least 60% of cluster cells individually predicted as the
   weighted winner.

These candidate cutoffs are sensitivity conditions, not assumed published
standards. The reviewed expected labels are never used during transfer or to
construct a vote; they are joined only into the final review table. Each
reference still contributes at most one vote, and every combined non-Unknown
label still requires agreement from at least two references.

On Cheaha, inspect all mappings and input paths without loading an RDS object:

```bash
cd /home/acflint/R/Projects/XeniumFCProject

Rscript \
  scripts/xenium_annotate_01c_reference_sampling_abstention_sensitivity.R \
  --list

bash \
  scripts/submit_xenium_annotate_01c_reference_sampling_abstention_sensitivity.sh \
  --dry-run
```

Submit the array and its `afterok` combine job with:

```bash
bash \
  scripts/submit_xenium_annotate_01c_reference_sampling_abstention_sensitivity.sh
```

Compact task tables and combined summaries are written beneath
`outputs/validation/resolution5_method_development/reference_sampling_abstention_sensitivity/`.
The most direct review files in `combined_tables/` are:

- `sampling_abstention_consensus_results.csv`, containing every equal
  two-of-three result;
- `sampling_abstention_unknown_summary.csv`, reporting cluster and cell
  coverage under each abstention rule;
- `sampling_abstention_reviewed_cluster_results.csv`, showing the outcomes for
  the five reviewed clusters; and
- `sampling_abstention_sampling_comparison.csv`, identifying labels changed by
  full-reference sampling.

No RDS object, production transfer table, consensus label, or downstream VZ/RL
output is written or replaced. Based on the completed capped-reference jobs,
`cap1000` tasks are expected to take approximately 8-50 minutes each. The full
reference tasks are estimated at approximately 20-120 minutes each, with low
confidence because their memory and anchor costs have not yet been measured.
At array concurrency three, all 18 tasks should take approximately 2-8 hours
total, excluding queue wait. The task launcher requests 4 hours, 1 CPU, and
128 GB on `short`; the combine job should take under 2 minutes but requests 15
minutes and 4 GB on `express`. After completion, collect `JobID`, `Elapsed`,
`State`, `MaxRSS`, `ReqMem`, `AllocCPUS`, and `ExitCode` with `sacct` before
changing these requests or applying either sensitivity condition to production.

##### Age/stage-stratified reference sampling sensitivity at resolution 5

`xenium_annotate_01d_stratified_sampling_sensitivity.R` adds a third,
isolated sampling mode to the completed `cap1000` versus `full` experiment.
`stratified_cap1000` retains at most 1,000 cells per `clusters_refined`
identity, but distributes that identity-specific allowance as evenly as
possible across developmental age/stage levels. If a small stratum is
exhausted, its unused allowance is redistributed among the remaining strata.
Ties and cell selection are reproducible with seed 42.

The test uses the same three reviewed query samples, Aldinger/Sepp/Science
references, RPCA settings, four abstention rules, and equal two-of-three
consensus as the preceding sensitivity analysis. It therefore changes only
which cells are retained within each capped reference identity. It does not
write an RDS or modify a production label, threshold, consensus result, or
the existing `combined_tables/` results.

A required preflight loads all three references and selects one complete,
informative developmental field for each reference. It prioritizes the known
Aldinger `age` and Sepp `Stage` fields and common Science PCW/age/stage field
names. A field must have no missing values and 2-50 distinct levels. The run
stops before label transfer if any reference cannot meet those safeguards.
The selected fields and level counts are recorded in
`stratified_preflight/reference_age_stage_strata.csv`.

On Cheaha, inspect the task mapping and reference paths without loading an RDS
object or submitting a job:

```bash
cd /home/acflint/R/Projects/XeniumFCProject

bash \
  scripts/submit_xenium_annotate_01d_stratified_sampling_sensitivity.sh \
  --dry-run
```

Submit the preflight, nine-task array, and dependent three-mode combine job:

```bash
bash \
  scripts/submit_xenium_annotate_01d_stratified_sampling_sensitivity.sh
```

New task tables are written under `stratified_cap1000/`. The most useful
three-mode summaries are written to `stratified_combined_tables/`:

- `three_mode_sampling_comparison.csv` compares the consensus label from
  `cap1000`, `stratified_cap1000`, and `full` for every cluster and rule;
- `three_mode_reviewed_cluster_results.csv` isolates the five reviewed
  clusters and reports agreement with the reviewed expected label;
- `three_mode_unknown_summary.csv` compares Unknown cluster and cell fractions;
  and
- `stratified_reference_strata_balance.csv` records the full and retained cell
  count for every reference class-by-stage combination.

The preflight is expected to take approximately 5-30 minutes and requests 1
hour, 1 CPU, and 64 GB on `express`; this estimate is low confidence because
reference RDS read times have not been measured separately. Based on the
observed 7.3-47.3-minute capped production transfers and 44.5-49.3 GB peak
memory, each stratified task is expected to take approximately 8-50 minutes
and requests 2 hours, 1 CPU, and 64 GB on `express`. At concurrency three, the
nine tasks should take approximately 25-150 minutes total, excluding queue
wait. The combine step should take under 2 minutes and requests 15 minutes and
4 GB on `express`. Collect `sacct` evidence after completion before changing
these requests or considering this sampling mode for production.

##### Reference-cap convergence screen at resolution 5

`xenium_annotate_01e_reference_cap_convergence.R` tests whether an
intermediate per-identity reference cap can reproduce the completed
full-reference results without requiring full-reference computation for every
sample. It runs caps of 2,500, 5,000, and 10,000 cells for the same three
reviewed query samples and Aldinger/Sepp/Science references. The transfer
settings, seed 42, four abstention rules, and equal two-of-three consensus are
unchanged from the preceding sampling analysis. The combine step reuses its
completed `cap1000` and `full` tables, yielding a five-mode comparison.

This is an isolated sensitivity screen. It writes compact CSV tables only and
does not change a production cap, threshold, reference, consensus label, RDS
object, or downstream VZ/RL result. The intended decision is the smallest cap
whose labels and confidence measures have plateaued toward the full-reference
result. Before adopting that cap, repeat the selected cap with additional
random seeds to confirm that the result is not specific to seed 42.

On Cheaha, first inspect all paths and the exact 27-task mapping without
loading the large RDS objects:

```bash
cd /home/acflint/R/Projects/XeniumFCProject

bash scripts/submit_xenium_annotate_01e_reference_cap_convergence.sh \
  --dry-run
```

Submit the 27-task array and its dependent combine job:

```bash
bash scripts/submit_xenium_annotate_01e_reference_cap_convergence.sh
```

Intermediate task tables are written under `cap2500/`, `cap5000/`, and
`cap10000/` within
`outputs/validation/resolution5_method_development/reference_sampling_abstention_sensitivity/`.
The combined decision tables are written to `cap_convergence_screen_tables/`:

- `cap_convergence_to_full_summary.csv` reports cluster- and cell-weighted
  agreement with full-reference consensus, Unknown burden, consensus support,
  and reviewed-cluster accuracy for every cap and abstention rule;
- `cap_convergence_cluster_comparison.csv` identifies the exact clusters and
  labels that differ from full-reference results;
- `cap_convergence_reviewed_cluster_results.csv` shows all five reviewed
  clusters across the five sampling modes and four rules;
- `cap_convergence_reference_vote_summary.csv` determines whether remaining
  disagreement is concentrated in Aldinger, Sepp, or Science; and
- `cap_convergence_unknown_summary.csv` and
  `cap_convergence_reference_class_balance.csv` document coverage and the
  exact number of reference cells retained.

Based on the observed 7.3-47.3-minute capped transfers and increasing retained
reference sizes, estimated runtimes are 10-60 minutes per `cap2500` task,
15-90 minutes per `cap5000` task, and 20-120 minutes per `cap10000` task. These
are low-confidence estimates until this screen supplies direct `sacct`
measurements. At concurrency three, the 27 tasks are expected to take roughly
2.25-13.5 hours total, excluding queue wait. Each task requests 4 hours, 1 CPU,
and 128 GB on `short`; the combine job should take under 2 minutes but requests
15 minutes and 4 GB on `express`. After completion, collect `JobID`, `Elapsed`,
`State`, `MaxRSS`, `ReqMem`, `AllocCPUS`, and `ExitCode` with `sacct`.

##### Focused `k.anchor` sensitivity at resolution 5

`xenium_annotate_01a_k_anchor_sensitivity.R` runs only
`GZFB_9_X_G_1`, `GZFB_20_X_G_1`, and `GZFB_22_X_G_3` against Aldinger,
Sepp, and Science with `k.anchor = 5` and `k.anchor = 10`. This is an
18-task isolated sensitivity analysis. It preserves the production RPCA
settings (`k.score = 30`, `k.weight = 50`, dimensions 1:30, 3,000 variable
features, at most 1,000 reference cells per identity, seed 42, and the legacy
0.4 reporting threshold) and changes only `k.anchor`. It writes compact score,
prediction-count, and provenance tables; it does not write an RDS or modify
production transfer or consensus outputs.

Inspect the exact task mapping before submission:

```bash
Rscript scripts/xenium_annotate_01a_k_anchor_sensitivity.R --list
for task_id in $(seq 1 18); do
  Rscript scripts/xenium_annotate_01a_k_anchor_sensitivity.R \
    --dry-run "${task_id}"
done
```

The dry-run checks paths and task/manifest mapping without loading the large
RDS objects. On Cheaha, submit the array and then its dependent combine job:

```bash
array_submit=$(sbatch --parsable \
  scripts/run_xenium_annotate_01a_k_anchor_sensitivity.slurm)
array_job=${array_submit%%;*}
sbatch --dependency="afterok:${array_job}" \
  scripts/run_xenium_annotate_01a_k_anchor_sensitivity_combine.slurm
```

Per-task outputs are written below
`outputs/validation/resolution5_method_development/k_anchor_sensitivity/k_anchor_<5-or-10>/<reference>/tables/`.
The compact comparison table is
`outputs/validation/resolution5_method_development/k_anchor_sensitivity/combined_tables/k_anchor_5_10_combined_cluster_score_summary.csv`.
Based on the existing 7.3-47.3-minute full transfer-task history and the
omission of UMAP, plots, and RDS output here, expect approximately 5-40 minutes
per task and approximately 15-120 minutes total at concurrency six. This is a
moderate-to-low-confidence estimate; the launcher requests 2 hours per task,
and the combine job should take under 1 minute while requesting 15 minutes.
Use `sacct` to collect `Elapsed` and `MaxRSS` after completion.

##### Minimal transfer-source diagnostic at resolution 5

`xenium_annotate_01b_transfer_diagnostics.R` tests the highest-priority
non-reference explanations for the reviewed RL/Granule/Glia/VZ discrepancies
without changing production labels. It uses only `GZFB_9_X_G_1`,
`GZFB_20_X_G_1`, and `GZFB_22_X_G_3`, all three references, resolution 5,
`k.anchor = 10`, `k.score = 30`, dimensions 1:30, at most 1,000 reference
cells per broad identity, seed 42, and the existing 0.4 reporting threshold.
It compares:

- `baseline_rpca`: the normalized data already stored in each object;
- `norm10k_rpca`: reference and query both re-normalized with LogNormalize and
  `scale.factor = 10000`;
- `norm10k_cca`: the same matched 10k normalization followed by CCA instead of
  RPCA.

Each AnchorSet is reused for `k.weight = 25` and 50. The resulting tables test
normalization compatibility, cross-platform reduction choice, transfer-neighbor
smoothing, anchor-label composition and coverage, RPCA mapping quality,
reference class/age/donor/subtype balance, marker-feature availability, and
cell-level score entropy and margins. CCA mapping scores are intentionally `NA`
because Seurat's `MappingScore` diagnostic is RPCA-specific. The focused
clusters and reasons for inclusion are declared in
`config/transfer_diagnostic_clusters.csv`; this manifest records review targets
but does not force an expected identity. The driver prefers this project-root
location on Cheaha and retains `scripts/config/` as a local-layout fallback.

Inspect the exact 27-task mapping and run path-only checks with:

```bash
Rscript scripts/xenium_annotate_01b_transfer_diagnostics.R --list
for task_id in $(seq 1 27); do
  Rscript scripts/xenium_annotate_01b_transfer_diagnostics.R \
    --dry-run "${task_id}" || exit 1
done
```

The dry-run checks input/output paths, sample-manifest membership, target
clusters, and fixed settings without loading the large RDS objects; expect
under 10 seconds per invocation. Submit the two arrays and one combine job on
Cheaha with:

```bash
rpca_submit=$(sbatch --parsable \
  scripts/run_xenium_annotate_01b_transfer_diagnostics_rpca.slurm)
rpca_job=${rpca_submit%%;*}

cca_submit=$(sbatch --parsable \
  scripts/run_xenium_annotate_01b_transfer_diagnostics_cca.slurm)
cca_job=${cca_submit%%;*}

sbatch --dependency="afterok:${rpca_job}:${cca_job}" \
  scripts/run_xenium_annotate_01b_transfer_diagnostics_combine.slurm
```

All three launchers use the literal email `acflint@uab.edu`. Arrays email on
failure; the dependent combine job emails on completion, failure, or an
invalid dependency. They also explicitly request Slurm account `acflint` and
QoS `normal` because those authorized values may not currently be inferred as
defaults. Outputs are isolated below
`outputs/validation/resolution5_method_development/transfer_diagnostics_k10/`.
The fastest review table is
`combined_tables/diagnostic_target_cluster_cross_reference_summary.csv`, with
per-cell cross-reference scores in
`combined_tables/diagnostic_target_cells_cross_reference.csv`.

Based on the observed 7.3-47.3-minute RPCA transfer history, expect roughly
10-60 minutes per RPCA task and 30 minutes to 2 hours per CCA task. At the
configured concurrency, the 18-task RPCA array should take about 30 minutes to
3 hours and the 9-task CCA array about 1.5-6 hours; because the arrays can run
concurrently, the expected end-to-end wall time is approximately 1.5-6 hours,
excluding queue wait. The RPCA estimate is moderate confidence; the CCA
estimate is low confidence because there is no project-specific CCA history.
The launchers request 2 hours/64 GB on `express` for RPCA, 4 hours/128 GB on
`short` for CCA, and 15 minutes/4 GB on `express` for the combine step. After
completion, collect `JobID`, `Elapsed`, `State`, `MaxRSS`, `ReqMem`,
`AllocCPUS`, and `ExitCode` with `sacct` to refine these requests.

To propagate a reviewed Sepp panel-reference update without rerunning Xenium
clustering or the unchanged Aldinger and Science transfers, use the dedicated
submitter. It uses the shared six-task resolution-5 transfer concurrency:

```bash
bash scripts/submit_sepp_reference_resolution5_consensus.sh --dry-run
bash scripts/submit_sepp_reference_resolution5_consensus.sh
```

If `sepp_02_subset_gene_panel.R` has already completed successfully, skip it
and start at the six-concurrent-task Sepp transfer:

```bash
bash scripts/submit_sepp_reference_resolution5_consensus.sh \
  --downstream-only --dry-run
bash scripts/submit_sepp_reference_resolution5_consensus.sh \
  --downstream-only
```

For a focused review after `sepp_02_subset_gene_panel.R` has already been
rerun, propagate the updated Sepp reference through only `GZFB_9_X_G_1`,
`GZFB_20_X_G_1`, and `GZFB_22_X_G_3`:

```bash
bash scripts/submit_sepp_reference_resolution5_three_sample_consensus.sh --dry-run
bash scripts/submit_sepp_reference_resolution5_three_sample_consensus.sh
```

The submitter resolves and verifies manifest task IDs 9, 17, and 31, reruns
only those three Sepp transfers, rebuilds only their three weighted 2-of-3
consensus tables, and replaces only their consensus objects and per-sample
plots. Existing Aldinger and Science transfers are reused. It does not rebuild
the 34-sample report PDFs. Expect approximately 8-41 minutes per Sepp transfer
and about the same transfer wall time at concurrency three, under 2 minutes for
the consensus merge, and approximately 1-4 minutes per consensus-application
task. The expected end-to-end active wall time is about 10-50 minutes,
excluding queue wait; these estimates are based on completed resolution-5
transfer history and the measured all-sample consensus workflow.

The dry-run validates the 34-row manifest and all Sepp transfer input paths but
does not load large RDS objects. The full workflow explicitly replaces the 34
Sepp transfer outputs, weighted 2-of-3 tables, consensus objects, paired plots,
Unknown marker workbooks, report pages, and both final PDFs. The Sepp panel script itself writes its stable RDS
and plot names directly and has no overwrite flag, so run the full mode only
after reviewing those existing reference outputs. The completed Sepp array
ran in 8.2-40.4 minutes per task with a mean of 19.0 minutes. At concurrency
six, expect approximately 2-3 hours for the array and roughly 6-23 hours for
the complete downstream chain, excluding queue wait. The transfer estimate is
high confidence from all 34 samples; downstream plotting remains lower
confidence.

Each preflight line begins with `DRY-RUN PASS` or `DRY-RUN FAIL`. The submitter
also ends with a boxed `DRY-RUN SUCCESS` or `DRY-RUN FAILED` summary that names
the failed stage and confirms whether any jobs were submitted. The
`outputs-existing` field is informational: it counts protected outputs already
on disk and does not itself indicate failure. The submitter suppresses repeated
`renv` synchronization notices during its 34 path checks without changing the
environment; inspect it once with `renv::status()` before production.
The preflight also requires `writexl` in the active project library before any
Slurm jobs are submitted because consensus application uses it for all
Unknown-cluster marker workbooks.

To propagate a reviewed Science-reference update without rerunning unchanged
Xenium clustering or the unchanged Aldinger and Sepp transfers, run this
protected dependency submitter from the project root on an allocated compute
node:

```bash
bash scripts/submit_science_reference_resolution5_consensus.sh --dry-run
bash scripts/submit_science_reference_resolution5_consensus.sh
```

If the two Science reference stages have already completed and only downstream
stages 3-7 need to be rerun, use:

```bash
bash scripts/submit_science_reference_resolution5_consensus.sh \
  --downstream-only --dry-run
bash scripts/submit_science_reference_resolution5_consensus.sh \
  --downstream-only
```

The dry-run checks the manifest, Science reference, and all 34 label-transfer
inputs, then submits nothing. Checks for the consensus and report stages are
deferred until their upstream dependency jobs have produced the required
inputs; each production launcher repeats its own path-only dry-run before
running. Active R drivers report dry-runs in a compact format such as
`DRY-RUN PASS | stage | inputs=2/2 | outputs-existing=0/15`; missing inputs and
failed checks are still named without printing every expected output path.
The full second command repeats the initial checks, rebuilds the two Science
reference stages, reruns all 34 resolution-5 Science transfers, rebuilds and
applies the weighted 2-of-3 consensus, renders both versions of all 34 report
pages, and merges both final PDFs. With `--downstream-only`, the two reference stages are skipped
and the existing Science reference output is used. Both submission modes
explicitly enable each submitted stage's protected overwrite switch; neither
submits preprocessing, clustering, VZ/RL, Giotto, or SpaTrack jobs.
`RLSVZ` remains mapped to `RL`.

If the two Science reference stages have already completed successfully, skip
them and start the protected dependency chain at the 34-sample Science label
transfer:

```bash
bash scripts/submit_science_reference_resolution5_consensus.sh \
  --downstream-only --dry-run
bash scripts/submit_science_reference_resolution5_consensus.sh \
  --downstream-only
```

The downstream-only dry-run verifies the existing panel-restricted Science
reference through the transfer input checks before submitting anything.

Observed Science-reference jobs took about 2-3 minutes each in the available
Slurm history, so their 1-hour requested limits are conservative. Three
completed 34-task Science transfer arrays ran in 7.7-47.3 minutes per task,
with means of 19.1-19.8 minutes. Expect about 2-3 hours per array at concurrency
six. Consensus construction should take under 2 minutes; consensus application,
paired plots, and marker tests should take 15-90 minutes per task and about
1.5-9 hours total at concurrency six; paired report pages should take 20-90
minutes per task and about 2-9 hours total at concurrency six; merging both
reports should take under 15 minutes. Approximate active wall time is 6-23
hours, excluding queue wait. Transfer estimates are high confidence; the
revised downstream estimates are low confidence and should be refined with
`sacct` after completion.

### 5. Annotation implementation and retired resolution-1.5 default mode

The implementation described here also powers resolution-5 production, but
the unflagged/default commands and `outputs/xenium/annotation/legacy_resolution1_5/` paths
below are retired resolution-1.5 workflow details. Do not run those default
launchers for new production work. Use the explicit resolution-5 commands in
the preceding section. The generic resolution-1.5 outputs may be removed once
the resolution-5 validation report passes.

`xenium_annotate_01_label_transfer_rpca.R` is the primary label-transfer implementation. It:

- uses the complete reference for `--all-samples-res5`; legacy/default modes
  retain their prior cap of at most 1,000 cells per reference identity;
- finds genes shared by the reference and Xenium object;
- builds reciprocal-PCA transfer anchors over 30 dimensions;
- transfers `clusters_refined` labels;
- applies a winner mean-score threshold of 0.4 and, for resolution-5
  production, a top-two mean-score margin of 0.05; and
- derives cluster-level majority and weighted-vote labels.

The reference and sample are selected explicitly from the command line. The
sample task ID maps to all 34 rows in `config/samples.csv`:

```bash
Rscript scripts/xenium_annotate_01_label_transfer_rpca.R --list
Rscript scripts/xenium_annotate_01_label_transfer_rpca.R --dry-run Aldinger 1
```

`run_xenium_annotate_01_transfer.slurm` submits a 34-task array, capped at two concurrent jobs. To
submit all three reference analyses:

```bash
sbatch --job-name=Xen_ABT_Aldinger --export=ALL,REFERENCE=Aldinger scripts/run_xenium_annotate_01_transfer.slurm
sbatch --job-name=Xen_ABT_Sepp --export=ALL,REFERENCE=Sepp scripts/run_xenium_annotate_01_transfer.slurm
sbatch --job-name=Xen_ABT_Science --export=ALL,REFERENCE=Science scripts/run_xenium_annotate_01_transfer.slurm
```

The driver refuses to replace any existing annotation output unless called
with `--overwrite`, or submitted with `ABT_OVERWRITE=true` after review.
For each reference, RDS objects, plots, and tables are written beneath
`outputs/xenium/annotation/legacy_resolution1_5/01_label_transfer/<reference>/`.
Combined voting tables are written to
`outputs/xenium/annotation/legacy_resolution1_5/02_consensus/tables/`.

Related scripts:

- `xenium_annotate_02_build_consensus.R` – merges the three reference comparisons, calculates one cluster-level `consensus_label`, and numbers insufficient-support clusters as `Unknown-##`;
- `xenium_annotate_03_apply_consensus.R` – adds consensus labels and PCW metadata to individual objects, makes consensus labels the active identities, writes all-cell and known-only consensus figures, and conditionally exports Unknown-cluster marker workbooks;
- `run_xenium_annotate_03_apply_consensus.slurm` – submits all 34 consensus-label tasks, capped at three concurrent jobs;
- `xenium_annotate_03b_plot_spatial.R` – additional consensus spatial maps;
- `xenium_annotate_03c_plot_proportions.R` – attaches sample metadata and plots consensus cell-type proportions;
- `xenium_annotate_03d_plot_report.R` – writes matched all-cell and known-only PCW-ordered PDF reports, with each page combining the global spatial map, faceted spatial map, UMAP, and marker DotPlot;
- `run_xenium_annotate_03d_plot_report.slurm` – runs the combined report as one sequential job so concurrent array tasks cannot write to the same PDF;
- `xenium_vz_rl_optional_consensus_spatial.R` – consensus maps after regional refinement/QC.

The consensus vote ignores missing and reference-level `Unknown` labels and selects the
most frequent label across Aldinger, Sepp, and Science majority/weighted
results. Final insufficient-support clusters are numbered `Unknown-01`,
`Unknown-02`, and so on in numeric cluster order. Ties retain the existing
comparison-table order. Original reference
labels remain in the object; `consensus_label` is added without overwriting
them. Inspect the task mapping and one task before running:

```bash
Rscript scripts/xenium_annotate_03_apply_consensus.R --list
Rscript scripts/xenium_annotate_03_apply_consensus.R --dry-run 1
```

Consensus objects are written to
`outputs/xenium/annotation/legacy_resolution1_5/03_consensus_labels/rds/<sample>_Consensus_annotated.rds`.
Consensus UMAP, spatial, and marker plots are written together under
`outputs/xenium/annotation/legacy_resolution1_5/03_consensus_labels/plots/`; only the separate
proportion-summary analysis retains its `plots/proportions/` subdirectory.
The stable filenames contain all cells, while matched `_Consensus_KnownOnly_`
files remove every `Unknown-##` cell. Unknown marker workbooks are written under
`03_consensus_labels/tables/unknown_markers/` whenever a sample has at least one
Unknown cluster.
The driver refuses to overwrite existing objects or plots unless explicitly
given `--overwrite` or submitted with `CONSENSUS_OVERWRITE=true`. PCW is read
from `metadata/samples_meta.xlsx` during this stage, so downstream objects
inherit it without a second full-size PCW RDS copy. Consensus UMAP and spatial
plots are TIFF-only; the marker DotPlot is saved as both TIFF and Cairo PDF.
The DotPlot displays the consensus identities in reverse y-axis order. To
regenerate only the DotPlot files for existing validated consensus objects,
without rewriting the RDS or spatial plots, submit with
`CONSENSUS_DOTPLOT_ONLY=true`:

```bash
sbatch --array=1-34%3 --export=ALL,CONSENSUS_DOTPLOT_ONLY=true \
  scripts/run_xenium_annotate_03_apply_consensus.slurm
```

After all consensus-label array tasks have completed, inspect the 34 report
inputs and PCW-defined page order without loading any RDS objects:

```bash
Rscript scripts/xenium_annotate_03d_plot_report.R --list
Rscript scripts/xenium_annotate_03d_plot_report.R --dry-run
```

Submit the report with:

```bash
sbatch scripts/run_xenium_annotate_03d_plot_report.slurm
```

The job writes
`outputs/xenium/annotation/legacy_resolution1_5/03_consensus_labels/plots/Xenium_Consensus_All_Samples_Report.pdf`
and the matched
`Xenium_Consensus_All_Samples_KnownOnly_Report.pdf`.
Each of its 34 landscape pages contains one sample's global spatial map,
faceted spatial map, UMAP, and marker DotPlot. Pages are ordered by numeric PCW,
with `config/samples.csv` order breaking ties, and each page title includes the
sample ID and PCW. The script validates each consensus object's PCW against the
configured sample metadata before plotting. Pages are rasterized at 150 DPI
before PDF embedding to keep the many spatial points from producing an
impractically large vector file. Existing reports are protected; after review,
an intentional replacement can be submitted with
`--export=ALL,CONSENSUS_REPORT_OVERWRITE=true`. The expected runtime is about
2–6 hours for all 34 samples (low confidence, based on static inspection and
sequentially loading/plotting one full consensus object at a time); the Slurm
launcher requests 8 hours. A measured per-page time from the first run would
substantially improve this estimate.

### 6. VZ analysis

The VZ branch now reads all 34 whole-tissue inputs from
`outputs/xenium/annotation/03_consensus_labels_weighted_2of3/rds/`,
as selected by the `production` configuration. This changes the upstream broad
identities from the retired resolution-1.5 annotation to the approved
resolution-5 weighted consensus. Consequently, any VZ/RL objects previously
derived from resolution-1.5 whole-tissue inputs are stale and must be
regenerated.

Regional artifacts now use semantic names such as `Xenium_VZ_integrated.rds`,
`Xenium_VZ_postQC.rds`, and `Xenium_VZ_subclusters.rds`. The old `Ald`,
`XenAld`, and `Res1.5` tokens are not part of the active VZ/RL filename
contract. Numerical `Res_0.5`, `Res_0.8`, and `Res_1` plot names still describe
the regional clustering resolutions actually evaluated by the scripts.

Before either regional branch starts, sample-qualified Unknown populations pass
through the independent post-consensus workflow below:

#### Broad-reference calibration and consensus

Broad routing is now evaluated at the VZ/RL/Other level for every cluster,
not only for populations already called Unknown. Original reference identities
and the established detailed `consensus_label` are preserved. `Other` is a
routing category defined by `config/reference_broad_taxonomy.csv`; reference
objects are not relabelled to a literal cell type named Other.

The protected calibration workflow is:

```bash
Rscript scripts/xenium_broad_01_validate_references.R --dry-run
sbatch scripts/run_xenium_broad_01_validate_references.slurm

Rscript scripts/xenium_broad_02_score_clusters.R --list
Rscript scripts/xenium_broad_02_score_clusters.R --dry-run 1
sbatch scripts/run_xenium_broad_02_score_clusters.slurm

Rscript scripts/xenium_broad_03_compile_calibration.R --dry-run
sbatch scripts/run_xenium_broad_03_compile_calibration.slurm
```

Stage 1 confirms that the production `_5k.rds` references contain no genes
outside `inputs/xenium_5k_genes.rds` and that every `clusters_refined` identity
maps exactly once to VZ, RL, or Other. Stage 2 reads the existing 102 annotated
reference/sample RDS objects, sums fine-label prediction probabilities within
each broad category at the cell level, and then reports cluster-level winner,
margin, winner-cell fraction, and entropy. It does not rerun label transfer.
Stage 3 evaluates prespecified score, margin, and cell-support grids separately
for each reference. Agreement of the other two raw broad winners is an internal
concordance target, not biological ground truth. It also audits all established
known clusters so threshold selection is not optimized only on Unknowns.

Review `Xenium_Broad_Reference_Threshold_Sensitivity.csv`, known-label
conflicts, and the current-numeric-threshold diagnostic. Then copy explicitly
chosen values and review notes to
`config/reviewed_broad_reference_thresholds.csv`, set all three `approved`
values to `TRUE`, and run:

```bash
Rscript scripts/xenium_broad_04_apply_reviewed_thresholds.R --dry-run
sbatch scripts/run_xenium_broad_04_apply_reviewed_thresholds.slurm
```

The apply stage refuses unapproved or incomplete thresholds. It writes a
versioned broad consensus beneath `05_broad_reference_consensus/`, preserves
the detailed consensus, and reports known-label conflicts without silently
relabeling them. Do not choose thresholds merely because they minimize the
Unknown count.

Expected runtime is about 5–20 minutes for reference validation; 2–15 minutes
per broad-score task and about 35 minutes–4.5 hours total at concurrency six;
2–10 minutes for calibration compilation; and under 2–10 minutes for reviewed
threshold application. These are low-confidence static estimates until matching
`sacct` observations are collected. Requested limits are 45 minutes, 2 hours,
30 minutes, and 30 minutes, respectively.

1. `xenium_unknown_01_audit_qc.R` audits one production sample without
   filtering or relabeling it. It joins `cell_area`, `nucleus_area`, and
   `nucleus_count` from the sample's raw Xenium `cells.csv.gz` or `cells.csv`
   by exact cell ID, verifies that raw and retained-object cell areas agree
   within the 0.00501-square-micrometer tolerance required by the retained
   table's two-decimal area rounding, and derives cell:nucleus area ratio and
   transcripts per cell area. The validation tolerance and maximum observed
   difference are written to each sample's provenance table.
2. `xenium_unknown_02_compile_qc_review.R` combines all 34 audits and writes the
   review template, then merges all sample QC PDFs in `config/samples.csv`
   order into `Xenium_Unknown_QC_All_Samples.pdf`.
3. `xenium_unknown_03_apply_qc_review.R` applies explicit whole-population
   `keep`/`remove` decisions one sample at a time.
4. `xenium_unknown_04_build_reference_routes.R` consumes the reviewed
   broad-reference votes and requires agreement from at least two references
   for VZ, RL, or Other. A reference abstains when its reviewed broad score,
   margin, or winner-cell-fraction requirement fails.
   - Optional: `xenium_unknown_04a_compile_reference_evidence.R` is a
     read-only review diagnostic. It joins each Unknown population to all
     three references' continuous winner/runner-up scores, score margins,
     threshold distances, and abstention reasons plus its whole-tissue top-20
     marker summary. It does not alter QC decisions, routes, labels, or RDS
     objects and is not an input to either regional branch.
5. `xenium_unknown_05_report_unrouted.R` reports retained `route_other` and
   `route_unresolved` populations without clustering or relabeling them.

These are separate, single-purpose entry points; production users do not select
review/apply/combine modes within one driver. Only standard `--list`,
`--dry-run`, and protected `--overwrite` options are exposed where relevant.
The older multi-mode annotation audit remains available for the established
resolution pilot launchers but is not the production entry point.

The production order is:

```bash
Rscript scripts/xenium_unknown_01_audit_qc.R --list
Rscript scripts/xenium_unknown_01_audit_qc.R --dry-run 1
sbatch scripts/run_xenium_unknown_01_audit_qc.slurm

Rscript scripts/xenium_unknown_02_compile_qc_review.R --dry-run
sbatch scripts/run_xenium_unknown_02_compile_qc_review.slurm

# Review the generated evidence, then populate config/reviewed_unknown_qc.csv.
Rscript scripts/xenium_unknown_03_apply_qc_review.R --dry-run 1
sbatch scripts/run_xenium_unknown_03_apply_qc_review.slurm

Rscript scripts/xenium_unknown_04_build_reference_routes.R --dry-run
sbatch scripts/run_xenium_unknown_04_build_reference_routes.slurm

# Optional score/margin and marker evidence table.
Rscript scripts/xenium_unknown_04a_compile_reference_evidence.R --dry-run
sbatch scripts/run_xenium_unknown_04a_compile_reference_evidence.slurm

Rscript scripts/xenium_unknown_05_report_unrouted.R --dry-run
sbatch scripts/run_xenium_unknown_05_report_unrouted.slurm
```

The first production QC audit (job 40905015), before raw Xenium morphology was
added, completed in 42–75 seconds per sample with a peak task MaxRSS of about
11.6 GiB. Reading the selected columns from `cells.csv.gz` should add seconds
to a few minutes per task, so a repeat is expected to take approximately 1–3
minutes per task and 12–36 minutes total at the launcher's three-task
concurrency, excluding queue time. This updated estimate is low-confidence
until the first morphology-enabled run; the requested limit remains 30 minutes
and 32 GiB per task. The first table compilation (job 40905016) completed in
28 seconds. Compilation plus
the new PDF merge should normally finish in under 1–2 minutes, while its
launcher requests 15 minutes. Applying reviewed QC decisions should take about
1–5 minutes per sample and 12–60 minutes total at concurrency three, with a
30-minute requested limit; this estimate remains low-confidence until its first
production run. Reference routing, optional reference-score/marker evidence
compilation, and the unrouted report should each finish in under 1–5 minutes
and request 15 minutes. The evidence compilation reads small CSV tables and
marker workbooks only; it does not load RDS objects.

QC tail fractions and cell-level tests are descriptive evidence only. They do
not automatically remove a population. Every removal requires an exact
`sample_id + consensus_label` decision with review notes. The default review
decision should be `keep` unless multiple QC measures support an artifact.
Upper-tail comparisons use values strictly above the within-sample q90, and a
metric with no variation among any audited cells is recorded as non-informative
with blank tail fractions rather than classifying tied zero values as high QC
tails.

Morphology evidence is also descriptive. The combined evidence reports median
nucleus area, cell:nucleus area ratio, transcript density, and nucleus count;
the fractions below the known-cell q10 for nucleus area and transcript density;
both q10 and q90 tails for the cell:nucleus ratio; missing or zero nucleus area;
zero nuclei; and nucleus counts other than one. Nucleus-area summaries and
tails use finite measurements only, while missing area is retained as a
separate fraction and is not imputed. These additions do not change the
original QC thresholds or automatically remove cells. Raw morphology joins and
source-file provenance are validated before outputs are written.

Reference routing uses the reviewed per-reference broad thresholds. Fine-label
scores are summed within VZ, RL, and Other before broad score and margin gates
are evaluated. `VZ`, `Glia`, `GABA`, `Purkinje`, and `OPC` map to `route_vz`;
`RL`, `Granule`, and `UBC` map to `route_rl`; `Meninges`, `Endothelial`, and
`Immune` map to `route_other`. A failed broad vote abstains. Populations without
two-reference broad agreement receive `route_unresolved`.

The optional Stage-4a evidence table is written to
`04_unknown_reference_routing/tables/Xenium_Unknown_Reference_Score_Marker_Evidence.csv`.
It contains one row per sample-qualified Unknown population. Marker adjusted
P-value counts are descriptive: no adjusted-significant positive top-20
markers means that the population lacks a distinctive program in that
one-versus-rest contrast, but does not by itself establish failed cell QC or
authorize removal. Score and margin distances are also diagnostic only and do
not relax the configured 0.4 winner-score or 0.05 top-two-margin thresholds.

Regional Stage 1 includes configured branch labels plus only retained Unknown
populations assigned to that branch. `route_other` and `route_unresolved`
cells remain in reviewed whole-tissue and combined outputs as `Unknown`, but do
not enter VZ or RL integration. The subsets carry `VZ_subset_reason` or
`RL_subset_reason`, the broad route, and reference-support count. Merge
manifests report configured-label and reference-routed Unknown cells
separately.

`Unknown-##` numbering is local to each sample and therefore must not imply
cross-sample identity. During both regional merges, the original
`consensus_label` is preserved for exact joins back to the whole-tissue object,
while `regional_consensus_label` appends the sample ID to every Unknown label
(for example, `Unknown-01_FB124_X_G`) and becomes the active identity. The
companion metadata columns `regional_unknown_label` and
`regional_consensus_label_pooled` support, respectively, an Unknown-only UMAP
with sample-qualified labels and an overview in which all Unknown cells are
pooled into one grey category. Known consensus labels are unchanged in all
three regional metadata columns. Each merge manifest records the exact
sample-qualified Unknown identities in `sample_qualified_unknown_labels`.

Because exclusive reference routing changes both regional cell sets, all
existing VZ and RL outputs from stage 1 onward are stale, even when an
early-stage filename did not change. Do not mix old and new stages. First run
the Stage-1 `--list` and `--dry-run` checks, then deliberately regenerate the
two Stage-1 arrays.
Existing stage-1 subsets and stage-2 merges are protected, so replacement
requires the branch's explicit `*_SUBSET_OVERWRITE=true` and
`*_MERGE_OVERWRITE=true` launcher variables after the old files have been
reviewed or archived. VZ and RL stages 1 through 7 can then proceed as
independent branches. Their Stage-6 review evidence and label application must
be regenerated because routed Unknown cells can change regional neighbors,
clusters, marker rankings, and therefore the scientifically appropriate
cluster-to-label mapping. Combined and spatial VZ/RL outputs are stale until
both refreshed stage-7 mappings exist.

The VZ branch generally follows this order:

For the 2026-10-15 expansion to **41 samples**, use
[`REBUILD_REGIONS_20261015.md`](REBUILD_REGIONS_20261015.md) and
`rebuild_xenium_regions_20261015.sh` after stage-1 subsets exist. The helper
archives old stage-2 through stage-9 results, queues manifest-defined merges,
fresh integration, QC review, and Harmony audits for both branches, and stops
before applying old removal IDs or annotation mappings. References to 34 samples
below describe the earlier cohort; active merge and integration drivers use the
current manifest. Preflight is estimated at 1-3 minutes without loading RDS
objects. Low-confidence active-time estimates through QC review are 5.5-21 hours
for VZ and 3-12 hours for RL, excluding queue waits; the linked guide gives the
per-stage basis and unchanged resource requests.

1. `xenium_vz_01_subset_samples.R` – extracts configured consensus identities per sample.
2. `xenium_vz_02_merge_samples.R` – verifies and incrementally merges all 34 sample subsets after removing spatial overhead.
3. `xenium_vz_03_integrate.R` – validates the complete merge, then runs normalization, PCA, Harmony integration, UMAP, and clustering.
   `xenium_vz_03a_audit_harmony.R` separately quantifies sample mixing,
   biological-label preservation, developmental-stage neighborhood change,
   and pre/post neighbor retention without changing the integrated object.
4. `xenium_vz_04_review_qc.R` – 4a reviews merged-cluster QC; 4b records the
   decision and conditionally materializes filtered sample objects.
5. `xenium_vz_05_reintegrate_post_qc.R` – reuses Stage-3 results after a
   zero-removal decision, or removes failed cells and recomputes after removals.
   Before selecting the Stage-6 partition, run the protected Stage-5c
   [clustree and subsampling stability review](REGIONAL_RESOLUTION_REVIEW.md)
   with `run_xenium_vz_05c_resolution_stability.slurm`. Existing completed
   Stage-5 objects can be reviewed directly; integration need not be rerun.
   For an efficient full-data comparison of 0.5/0.6/0.7/0.8, use
   [VZ resolution comparison](VZ_RESOLUTION_COMPARISON.md) and
   `run_xenium_vz_05d_compare_resolutions.slurm`. This reuses the saved graph
   and 0.5/0.8 assignments, calculates only 0.6/0.7, and supports a validated
   partition checkpoint for retries. Production labels remain unchanged.
6. `xenium_vz_06a_build_replicate_marker_evidence.R` – builds sample-aware
   positive and anti-marker candidates; `xenium_vz_06_annotate_subclusters.R`
   retains cell-level Wilcoxon rankings as descriptive review evidence and
   applies explicitly reviewed labels.
7. `xenium_vz_07_map_to_samples.R` – maps refined VZ labels back to all 34 configured samples.
8. `xenium_vz_plot_sample.R` and `xenium_vz_08_plot_samples.R` – manifest-defined per-sample spatial reports.
9. `xenium_vz_09_plot_counts.R` – VZ subcluster percent abundance within each plotted lineage.

The VZ driver reads all 34 samples from `config/samples.csv` and supports
`--list`, `--dry-run`, and protected `--overwrite` operation. Submit
per-sample extraction with `run_xenium_vz_01_subset.slurm`. Submit
merging with `run_xenium_vz_02_merge_samples.slurm`, integration with
`run_xenium_vz_03_integrate.slurm`, post-QC integration with
`run_xenium_vz_05_reintegrate_post_qc.slurm`, and subcluster annotation with
`run_xenium_vz_06_annotate_subclusters.slurm`.
`run_xenium_vz_07_map_to_samples.slurm` maps the refined identities, and
`run_xenium_vz_08_plot_samples.slurm` launches 34 protected plotting tasks.
Every VZ Stage-8 per-sample filename begins `<sample>_VZ_`; the corresponding
RL outputs begin `<sample>_RL_`, so copied artifacts remain branch-identifiable
outside their original directories. Combined reports use the same branch-aware
one-page report names as inputs.

After Stage 3, run the quantitative Harmony audit independently:

```bash
Rscript scripts/xenium_vz_03a_audit_harmony.R --dry-run
sbatch scripts/run_xenium_vz_03a_audit_harmony.slurm
Rscript scripts/xenium_rl_03a_audit_harmony.R --dry-run
sbatch scripts/run_xenium_rl_03a_audit_harmony.slurm
```

Each audit uses a fixed-seed stratified sample of at most 100,000 cells and 30
neighbors. Increased cross-sample mixing is interpreted only alongside
retention of biological-label purity, PCW neighborhood structure, and pre/post
neighbor overlap. Expected runtime is approximately 10–60 minutes per branch,
low confidence; each launcher requests 2 hours and 256 GB on `express`.
The 2026-10-16 audit fix releases the integrated Seurat object before neighbor
searches, retains only the required metadata and sampled embeddings, and logs
each major phase. The memory increase supplies headroom for the initial RDS
load after VZ job 40951436 segfaulted at 61.20 GB MaxRSS under a 64-GB request;
this observation suggests memory pressure but is not a confirmed Slurm OOM.
The age parser accepts both numeric values and production `PCW12` labels.
Dimensions, neighbor count, sample cap, seed, and metric definitions are
unchanged. See `HARMONY_AUDIT_RETRY_20261016.md` to retry only the audits.

After the post-QC objects exist, run the combined resolution-0.8 review.
Both VZ and RL review build missing replicate-aware evidence first and then
create cluster review outputs in the same job. Matching evidence is reused;
a separate Stage-06a submission is no longer required. The path-only
preflights take under one minute per branch:

```bash
Rscript --max-ppsize=500000 scripts/xenium_vz_06_annotate_subclusters.R \
  --review-only --resolution=0.8 --dry-run
Rscript --max-ppsize=500000 scripts/xenium_rl_06_annotate_subclusters.R \
  --review-only --resolution=0.8 --dry-run
```

After the selected branch passes preflight, submit its existing launcher:

```bash
sbatch scripts/run_xenium_vz_06_annotate_subclusters.slurm
sbatch scripts/run_xenium_rl_06_annotate_subclusters.slurm
```

Each defaults to review at resolution 0.8, one CPU, 256 GB and a 12-hour
limit. Combined runtime is estimated at 1–10 hours for VZ and 0.75–7 hours
for RL (low confidence from static inspection and preceding separate-stage
estimates; queue wait excluded). Missing evidence is planned output in the
preflight; partial/stale evidence is protected until explicit overwrite is
authorized. Standalone Stage-06a scripts/launchers remain available for
compatibility and evidence-only work. See [STAGE6_MERGE_20261017.md](STAGE6_MERGE_20261017.md).

These stages compare each raw subcluster with the remaining branch cells within
each biological sample, then summarize direction, magnitude, detection
contrast, and a sample-level sign test. They output separately ranked positive
and anti-marker candidates. Cell-level Wilcoxon tables remain useful for
descriptive ranking, but their P-values are not treated as biological-replicate
inference. The production resolution-0.8 VZ/RL review and label-application
steps verify these five evidence files and their input-RDS, cluster-column,
size, and modification-time provenance. Review builds them when absent;
label application still requires matching completed evidence and review.
The optional evidence-only launchers retain four-hour limits.

Each VZ Stage-8 global or lineage-highlight plot renders opaque grey context
cells at size 0.5 first and colored subcluster points at size 0.75 second, so
context cannot cover the highlighted cells. Faceted plots omit context cells
and display only the colored cells assigned to each named cluster. It
writes both the all-inclusive global map and a `Global_Spatial_NoOther` map,
a lineage-faceted TIFF, a separate Other-faceted TIFF, marker DotPlots, and a
one-page `Spatial_Report.pdf` containing both global maps, the lineage facets,
and marker DotPlot. Every individual-sample DotPlot and its report page use
the shared `vz_all_compact_markers` panel: the same current marker pool as the
Stage-6 production DotPlot, restricted to marker groups represented by that
sample. Unused groups and cycling-marker columns are omitted. The plotting code constructs all figures before writing any
output, requires every requested marker in the Xenium assay, and uses that assay
explicitly for DotPlot features, so a plotting error cannot silently leave a
report that appears complete.
Standalone TIFFs draw points directly on the final 600-dpi raster device so
their requested point size is preserved without a blurry intermediate raster.
Only spatial panels embedded in the merged PDF are pre-rasterized with
Seurat's imported `scattermore` dependency; `ggrastr` is retained as a
secondary raster backend. The report deliberately does not fall back to
millions of vector PDF points because that can overflow R's protection stack
on full Xenium samples.
DotPlots select mapped identities directly from the source object rather than
subsetting its FOV and centroid structures, avoiding an unnecessary spatial
object copy for each sample.
Stage 8 opens `CairoTIFF` directly in inch units at 600 dpi; the array driver
rejects TIFF outputs smaller than 10 KB so an incorrect device dimension
cannot be mistaken for a completed report.
Global and lineage-highlight TIFF names record the native 0.75 colored-point
size; faceted TIFF names record the native 0.6 colored-point size. In global
and highlight plots, opaque grey context cells use size 0.5 and are rendered
before the larger colored cells.

VZ outputs follow the script order under `outputs/xenium/vz/`:
`01_subsets`, `02_merged`, `03_integrated`, `04_qc`, `05_post_qc`,
`06_subclusters`, `07_mapped`, `08_sample_reports`, and
`09_cluster_counts`. Artifact types such as `rds`, `plots`, and `tables` are
separated within each stage.

Before merging, inspect completeness without loading any Seurat objects:

```bash
Rscript scripts/xenium_vz_02_merge_samples.R --dry-run
```

The merge requires exactly the 34 manifest-defined subset files and writes
`outputs/xenium/vz/02_merged/rds/Xenium_Merged_VZSubsets.rds` plus a CSV
manifest under `02_merged/tables/` with each sample's input path, total cell
count, configured-label count, reference-routed Unknown count, and PCW.

After that merge exists, inspect processing readiness without loading it:

```bash
Rscript scripts/xenium_vz_03_integrate.R --dry-run
```

Processing writes the stable object
`outputs/xenium/vz/03_integrated/rds/Xenium_VZ_integrated.rds`. The QC and
post-QC scripts consume that name. Existing processing RDS or UMAP outputs
require an explicit `--overwrite` rerun after review.

The critical stage-3 path computes only the corrected Harmony UMAP. The
pre-Harmony diagnostic UMAP is deferred because it does not affect Harmony,
neighbors, clustering, QC, or downstream labels. Corrected UMAP and clustering
use separate worker settings: the corrected UMAP can use all four allocated
CPUs, while clustering defaults to one process after three forked clustering
workers reached 422 GB MaxRSS and failed in VZ job 40676165. The launcher maps
`VZ_INTEGRATE_CLUSTER_WORKERS` to `XENIUM_STAGE3_CLUSTER_WORKERS`; leave it
unset at the safe default of 1 unless new memory evidence supports increasing
it. The launcher requests 600 GB on `medium`, retaining four CPUs for the
pre-clustering UMAP while clustering remains single-worker. This resource-only
increase follows two separate October 13 node losses during the same VZ
clustering phase: full Stage-3 job 40925744 lost `c0194` after 2 hours 43
minutes, and precluster-resume job 40930904 lost `c0191` after 58 minutes.
Both jobs requested 256 GB, both nodes became non-responsive before MaxRSS
could be returned, and neither failure supplied evidence that the scientific
parameters were invalid. The larger request reserves headroom for the large
SNN graph and reduces co-tenancy on the 754-GB Intel nodes without changing
dimensions, UMAP neighbors, graph construction, resolutions, seed, or labels.

Before clustering, stage 3 writes the provenance-checked restart object
`Xenium_VZ_precluster_checkpoint.rds`; immediately after clustering it writes
`Xenium_VZ_postcluster_checkpoint.rds` before creating any plots. Resume mode
automatically prefers the post-clustering checkpoint when it exists, otherwise
it uses the pre-clustering checkpoint. Checkpoints are written to a temporary
file in the destination directory and atomically renamed only after a nonempty
write completes, preventing a node loss during `saveRDS()` from replacing the
last valid checkpoint with a partial file. The launcher explicitly permits
requeue, appends logs across attempts, and automatically selects checkpoint
resume when `SLURM_RESTART_COUNT` is positive and a checkpoint exists. Slurm
therefore does not repeat normalization, PCA, Harmony, UMAP, and neighbor
construction after a node-level restart. Each stage-3 script also parses its
complete source at startup so a damaged later line fails before loading data or
running analysis. If a job fails after either checkpoint, inspect and resume:

```bash
Rscript scripts/xenium_vz_03_integrate.R --resume-checkpoint --dry-run
sbatch --export=ALL,VZ_INTEGRATE_RESUME_CHECKPOINT=true,VZ_INTEGRATE_OVERWRITE=true \
  scripts/run_xenium_vz_03_integrate.slurm
```

A precluster resume is expected to take approximately 1-4 hours, excluding
queue wait; this is low confidence because both recent attempts lost their
nodes during clustering before reporting final resource use. The launcher
retains a 24-hour limit to provide margin. After the first successful 600-GB
run, collect `JobID`, `Elapsed`, `State`, `MaxRSS`, `ReqMem`, `AllocCPUS`, and
`ExitCode` with `sacct` before reducing the request.

Stage 3 writes
`Xenium_VZ_UnknownClusters_BySample_UMAP.tif`, where each independently
numbered Unknown cluster includes its sample ID, and
`Xenium_VZ_ConsensusLabel_UMAP.tif`, where all Unknown cells are pooled as
`Unknown` and colored grey. These are integration summaries, not QC evidence.
The sample-qualified Unknown legend is restricted to two columns on a taller
canvas. To regenerate the corrected stage-3 UMAPs from the existing integrated
object without rerunning normalization, Harmony, clustering, or rewriting the
RDS, inspect and submit plot-only mode with:

```bash
Rscript scripts/xenium_vz_03_integrate.R --plots-only --dry-run
sbatch --export=ALL,VZ_INTEGRATE_PLOTS_ONLY=true,VZ_INTEGRATE_OVERWRITE=true \
  scripts/run_xenium_vz_03_integrate.slurm
```

Generate the deferred pre-/post-Harmony batch comparison separately with:

```bash
Rscript scripts/xenium_vz_03_integrate.R --plots-only --include-uncorrected --dry-run
sbatch --export=ALL,VZ_INTEGRATE_PLOTS_ONLY=true,VZ_INTEGRATE_INCLUDE_UNCORRECTED=true,VZ_INTEGRATE_OVERWRITE=true \
  scripts/run_xenium_vz_03_integrate.slurm
```

Until that deferred job completes, any older
`Xenium_VZ_Batch_Comp_UMAP.tif` should be treated as stale; full stage 3 warns
when it detects that file but does not delete or overwrite it.

The QC violin PDF, resolution-0.8 cluster summary, and top-10 marker table are
created only by the explicit stage-4 review job described below, under
`outputs/xenium/vz/04_qc/plots/` and `outputs/xenium/vz/04_qc/tables/`.

VZ QC is an explicit two-part review gate. Stage 4a first inspects readiness,
then generates the
resolution-0.8 QC summary, combined violin PDF, marker table, sample-qualified
Unknown evidence, sample-by-cluster QC/composition evidence, and review
manifest:

```bash
Rscript scripts/xenium_vz_04_review_qc.R --qc-only --dry-run
sbatch scripts/run_xenium_vz_04_review_qc.slurm
```

VZ stage 4 also writes
`Xenium_VZ_RawSubcluster_BySample_QC_Summary.csv`, containing count and feature
quartiles, within-sample median ratios, sample and cluster fractions, and
VZ/Glia/GABA/Purkinje/OPC/Unknown composition for every observed sample/cluster
pair. `Xenium_VZ_RawSubcluster_BySample_QC_Heatmaps.pdf` shows the log2 median
count and feature ratios relative to all VZ-branch cells from the same sample,
while `Xenium_VZ_RawSubcluster_SampleComposition_Heatmap.pdf` shows the fraction
of each expression cluster contributed by each sample. These descriptive
outputs mirror the RL stage-4 checks and do not impose a new threshold or
remove cells automatically.

The sample-qualified Unknown evidence includes:

- `Xenium_VZ_SampleQualifiedUnknown_ClusterComposition.csv`, an observed-cell
  cross-tabulation of resolution-0.8 expression cluster, sample, and qualified
  Unknown population, with both cluster and population fractions;
- `Xenium_VZ_SampleQualifiedUnknown_QC_Summary.csv`, including count/feature
  quartiles, sample-matched configured-cell medians, median ratios, and the
  fraction of each population below those sample-matched medians;
- `Xenium_VZ_SampleQualifiedUnknown_QC_Violins.pdf`, with count and feature
  distributions for each sample's Unknown populations beside configured
  regional cells from that same sample;
- `VZ_SampleQualifiedUnknown_top10_Markers.csv`, containing positive markers
  from within-sample one-versus-all-other-regional-cell tests; and
- `VZ_SampleQualifiedUnknown_Marker_Status.csv`, recording target/reference
  cell counts, completed comparisons, empty results, and populations skipped
  because either comparison group has fewer than three cells.

These outputs are descriptive review evidence only. They do not impose a new
QC threshold or remove cells automatically.

Review outputs created before the sample-by-cluster checks were added are
incomplete for the current decision gate. After confirming that the existing
review evidence may be replaced, regenerate the complete Stage-4a evidence
with:

```bash
VZ_QC_MODE=review \
  VZ_QC_OVERWRITE=true \
  sbatch --export=ALL \
  scripts/run_xenium_vz_04_review_qc.slurm
```

Stage 4b records the reviewed decision explicitly. Because the VZ branch now
includes cells whose consensus label is explicitly `VZ`, cluster membership
and numbering from the earlier branch are stale. The launcher therefore has
no apply-mode default: `VZ_QC_REMOVE_CLUSTERS` must explicitly contain the
newly reviewed comma-separated cluster IDs or the literal value `none`.
The script validates all manifest-defined whole-tissue inputs (currently 41), protects existing outputs,
and records per-sample removal counts. Post-QC processing reads that manifest
and cannot independently choose a different cluster. After inspecting the
regenerated Stage-4a evidence, the current approved VZ decision is to remove
pre-QC resolution-0.8 clusters `10`, `23`, `24`, `25`, `27`, and `28`
(updated 2026-10-18). The original five-cluster decision from 2026-10-16 removed
10 cells; it is superseded by adding pre-QC cluster 10. The new post-QC review
contains 1,447,477 cells. These pre-QC cluster IDs must not be confused with
new post-QC IDs. All other pre-QC clusters are retained by this decision.

Before submission, run the path-only preflight (expected under 1 minute; it
does not load RDS objects):

```bash
Rscript --max-ppsize=500000 scripts/xenium_vz_04_review_qc.R \
  --remove-clusters=10,23,24,25,27,28 --dry-run
```

Submit the approved decision with:

```bash
VZ_QC_MODE=apply \
  VZ_QC_REMOVE_CLUSTERS='10,23,24,25,27,28' \
  sbatch --export=ALL \
  scripts/run_xenium_vz_04_review_qc.slurm
```

Use the no-removal form only after a future review explicitly replaces the
current decision:

```bash
VZ_QC_MODE=apply \
  VZ_QC_REMOVE_CLUSTERS=none \
  VZ_QC_OVERWRITE=true \
  sbatch --export=ALL \
  scripts/run_xenium_vz_04_review_qc.slurm
```

`VZ_QC_OVERWRITE=true` replaces existing Stage-4b decision outputs after the
refreshed review; it does not overwrite Stage-4a evidence in apply mode.

Do not place a comma-separated cluster list directly inside Slurm's
`--export=...` argument: Slurm uses commas to delimit exported variables and
would pass only the first cluster ID. Set the shell environment variable as
shown above and submit with `--export=ALL` so the entire list reaches the
launcher.

The launcher defaults to Stage 4a review mode, so it cannot apply a decision
unless `VZ_QC_MODE=apply` and an explicit `VZ_QC_REMOVE_CLUSTERS` value are
supplied. It requests 12 hours, 1 CPU, and 256 GB on `short`. Stage 4a is expected to take
approximately 2-8 hours because it runs one sample-matched marker comparison
per qualified Unknown population. A Stage-4b `none` decision writes only a
routing manifest and should finish in under 1 minute; the reviewed nonzero
removal materializes all 41 manifest-defined filtered objects sequentially and
is expected to take approximately 20-90 minutes (low confidence). The resource
history records successful sequential VZ job 40758964 at 15.42 minutes for 34
samples, with 61.51 GB MaxRSS. The current cohort has 41 samples and a roughly
79-GB integrated object; object sizes and filesystem throughput make linear
scaling uncertain. Even a ten-cell removal
requires reading and writing every whole-tissue object. After completion,
check that the removal manifest contains all 41 samples, the approved cluster
set, and a total of 10 removed cells. Regenerate stage 5 and downstream VZ
annotations, mapped labels, and analyses after applying this decision. Preserve
the stage-3 object and review evidence to allow a reviewed rollback. These
estimates exclude queue wait. The 12-hour limit and 256-GB memory request are
retained conservatively; collect current-run elapsed time and MaxRSS before
revising them.

The Stage-4 launcher starts the main R process with `--max-ppsize=500000`, the
R-supported maximum pointer-protection-stack size. Materialization is now
deliberately sequential in that same process, matching RL, and the integrated
object is released before individual samples are loaded. These safeguards were
added after the earlier cluster-15-only VZ apply job 40697743 failed with
`protect(): protection stack overflow` after 11 minutes 45 seconds at 61.53 GB
MaxRSS, well below its 256 GB request, before writing any filtered RDS or
removal manifest. They change only runtime behavior; they do not change the
reviewed clusters, filtering logic, cell identifiers, or stable output names.

Inspect post-QC integration without loading the merged Seurat object:

```bash
Rscript scripts/xenium_vz_05_reintegrate_post_qc.R --dry-run
```

When Stage 4b records zero removals, VZ stage 5 reuses the Stage-3 PCA,
Harmony integration, SNN graph, UMAP, and resolution-0.5/0.8/1 assignments.
It does not calculate resolution 0.3. Any nonzero removal retains the full
PCA/Harmony/UMAP/neighbors/reclustering path at resolutions 0.3/0.5/0.8/1.
The full path uses 30 PCA/Harmony dimensions, 30 UMAP neighbors, and
`k.param=30`; resolution 0.8 is the user-selected Stage-6 review partition,
while other Stage-5 partitions remain available as sensitivity checks. The original-label
UMAP uses `regional_consensus_label_pooled` so independently numbered Unknown
populations do not compress the plotting panel with a sample-qualified legend.
The zero-removal path is expected to take approximately 10-45 minutes,
dominated by loading, joining, plotting, and saving.

After the Stage-5 dry run passes, start the first run with the new parameters
from the Stage-3 input rather than an old checkpoint:

```bash
VZ_POST_QC_OVERWRITE=true \
  sbatch --export=ALL scripts/run_xenium_vz_05_reintegrate_post_qc.slurm
```

Full VZ reintegration saves
`Xenium_VZ_postQC_precluster_checkpoint.rds` immediately after PCA, Harmony,
UMAP, and SNN construction, then saves
`Xenium_VZ_postQC_postcluster_checkpoint.rds` before plotting. Resume mode
validates the merged input, review and removal manifests, sample set, cell
counts, seed, resolutions, and the configured PCA/UMAP/graph parameters before
skipping completed work. Checkpoints from the earlier 50-dimension,
100-UMAP-neighbor run are deliberately rejected and cannot be used for this
first rerun:

```bash
Rscript scripts/xenium_vz_05_reintegrate_post_qc.R \
  --resume-checkpoint --dry-run
VZ_POST_QC_RESUME_CHECKPOINT=true \
  VZ_POST_QC_OVERWRITE=true \
  sbatch --export=ALL scripts/run_xenium_vz_05_reintegrate_post_qc.slurm
```

The full-reintegration launcher requests 36 hours, 1 CPU, and 600 GB on
`medium`. Clustering defaults to one sequential worker after job 40958650
failed with OUT_OF_MEMORY at 600 GB using four fork workers. Forked memory
sharing did not bound transient clustering allocations. The scientific graph,
Louvain settings, resolutions, and random seed are unchanged. A pre-clustering
checkpoint resume is expected to take 8–24 hours; a fresh full run may take
10–28 hours (low confidence, based on the recent OOM and older sequential
12-hour timeout under different parameters). Exclude queue wait. Use
`VZ_POST_QC_RESUME_CHECKPOINT=true` only after a matching checkpoint preflight.
Leave overwrite disabled unless replacing reviewed outputs is authorized.
The 600-GB request is provisional, not a validated peak-memory bound.

This stage writes the stable object
`outputs/xenium/vz/05_post_qc/rds/Xenium_VZ_postQC.rds` and a provenance
manifest. Stage 6 generates post-QC review evidence without assigning
subcluster identities. Resolution 0.8 is the default; select one of the
Stage-5 partitions with `--resolution=0.3`, `--resolution=0.5`, or
`--resolution=0.8`. For a protected resolution-0.8 review, inspect paths first
and then explicitly replace any stale resolution-0.8 evidence:

```bash
Rscript scripts/xenium_vz_06_annotate_subclusters.R \
  --review-only --resolution=0.8 --dry-run
sbatch --export=ALL,VZ_SUBCLUSTER_MODE=review,VZ_SUBCLUSTER_RESOLUTION=0.8,VZ_SUBCLUSTER_OVERWRITE=true \
  scripts/run_xenium_vz_06_annotate_subclusters.slurm
```

Review the resolution-specific raw-cluster UMAP, stable top-10 table, expanded
top-50 table, review manifest, and focused lineage-marker DotPlot beneath
`outputs/xenium/vz/06_subclusters/`. The DotPlot contrasts VZ progenitor,
cell-cycle, Purkinje, GABAergic, iCN, radial/Bergmann glial,
astrocyte/ependymal, pre-OPC/OPC, RL-VZ, immature/mature UBC, eCN, microglial,
granule, and endothelial-contrast programs. The deliberately broad panel is
intended for review and can be reduced after inspecting the first production
DotPlot. Its matching marker-availability CSV records every requested gene and
whether it exists in the active Xenium panel; absent genes are reported and
omitted rather than causing a plotting failure. Both the focused panel and the
shared `vz_markers` list from `color_palette.R` are included. Palette genes
already present in the focused panel retain their biological facet, while new
palette-only genes are appended without duplicate dots. Both marker tables are
derived from the same `FindAllMarkers()` result, so changing
the review resolution does not change `min.pct`, log-fold-change, sampling, or
statistical-test settings. The UMAP coordinates are reused from `umap_clean`;
only the displayed cluster assignments change.

The user-approved resolution-0.8 VZ map, updated 2026-10-19 for cluster 17,
covers the new
1,447,477-cell object's clusters 0–24 after pre-QC removal of 10,23,24,25,27,28.
See [VZ_RL_RESIDUAL_LABEL_MAPPING_20261018.md](VZ_RL_RESIDUAL_LABEL_MAPPING_20261018.md)
for the exact map and commands through all three regional per-sample mappings.
Cluster 20 is GABA Interneurons; 22 remains Unresolved; 23 matches cluster 0
as Patterning PCs; 18 is OPCs; 17 is provisionally Differentiating GCs rather
than Mature UBCs. The former 2026-10-17 map is historical.
Resolutions 0.3 and 0.5 remain review-only until separately approved.

After syncing the updated driver, palette, and sample-plot helper, preflight
label application (under one minute, paths/small CSVs only):

```bash
Rscript --max-ppsize=500000 scripts/xenium_vz_06_annotate_subclusters.R \
  --apply-reviewed-labels --resolution=0.8 --dry-run
```

After PASS and checking existing outputs, submit without enabling overwrite:

```bash
VZ_SUBCLUSTER_MODE=apply \
  VZ_SUBCLUSTER_RESOLUTION=0.8 \
  sbatch --export=ALL scripts/run_xenium_vz_06_annotate_subclusters.slurm
```

Expected application runtime is 15–120 minutes (low confidence from an older
11.53-minute RL apply run and the larger 1.50-million-cell VZ object). The
launcher requests 256 GB, one CPU and 12 hours; queue wait excluded. All
mapped VZ labels and downstream combined labels, reports/counts and spatial
outputs become stale after applying this map and must be regenerated in
workflow order. No stage-5 reintegration is needed for annotation-only changes.
Existing outputs remain protected until explicitly approved for replacement.

The current production DotPlots use the approved VZ taxonomy with 20 observed
labels and 34 unique marker columns. Purkinje markers are ordered DAB1,
FOXP1, FOXP4, EBF2, EBF1, EPHA4; GABA markers are TFAP2B, NXPH2, SLC6A5,
GRIK3, SP9, SOX14, DMBX1; glial markers are SOX2, TOP2A, TNC, AQP4, FOXJ1.
GABA identities are displayed as GABA Progenitors, Golgi Cells, MLIs,
GABA Interneurons, iCN in legends and DotPlots. HES6, SOX10, and C1QL3
are excluded from these production display panels. The removed Neuroblast
and UBC panels remain absent. The VZ Granule panel uses NEUROD1,
CBLN1, and RELN for the provisional Differentiating GCs identity.
Unresolved is
retained in global/Other plots without a forced marker panel or lineage.
The Granule header uses the broad Granule color; eCN uses its subcluster color.
Per-sample panels omit absent groups and unused identity levels. Raw diagnostic
review panels remain independent and unchanged. The current implementation
supersedes the earlier panel in VZ_DOTPLOT_OPTIMIZATION_20261018.md.
Marker panels and identity membership are defined in `color_palette.R`. Stage 6
refuses to write plots if a requested gene is absent from the Xenium assay, and
its final validation fails unless every required output exists and is nonempty.

The 2026-10-19 cluster-17 correction changes only the VZ resolution-0.8 fine
label from Mature UBCs to provisional Differentiating GCs. The original
`consensus_label` remains unchanged, while a rerun of corrected-broad Stage 02
assigns these cells to Granule. RL and residual Mature UBCs remain unchanged.
After all existing VZ/RL jobs finish and the numbered folder migration is
complete, rerun VZ Stage 06 application and Stage 07 mapping, then combined
Stage 01, corrected-broad Stage 02, and both Stage 03 merge branches. Existing
outputs require each launcher's explicit overwrite switch for this user-approved
replacement. Downstream non-spatial processing/plots, spatial integration,
combined reports, and VZ per-sample reports are stale until regenerated.
Raw marker review, VZ QC, and VZ post-QC integration do not need repeating.

```bash
Rscript scripts/xenium_vz_06_annotate_subclusters.R \
  --apply-reviewed-labels --resolution=0.8 --dry-run
sbatch --export=ALL,VZ_SUBCLUSTER_MODE=apply,VZ_SUBCLUSTER_RESOLUTION=0.8 \
  scripts/run_xenium_vz_06_annotate_subclusters.slurm
```

Application is refused unless the review outputs are nonempty and their
provenance matches the exact post-QC RDS, cell count, selected cluster column
and cluster set, random seed, marker-table settings, and focused marker panel.
Stage 6 then validates complete cluster-to-label coverage
before writing
`outputs/xenium/vz/06_subclusters/rds/Xenium_VZ_subclusters.rds`.
After an annotated VZ object exists, `--umaps-only` (launcher mode
`VZ_SUBCLUSTER_MODE=refresh-umaps`) safely redraws only the all-inclusive and
VZ-lineage UMAP TIFFs. The lineage-only output is named
`Xenium_VZ_Lineage_Subcluster_UMAP.tif`; it replaces the obsolete
`rmGCP,eCN,Cyc` filename. Every VZ label uses a white box with line-free repulsion;
the annotated RDS, label manifest, tables, and other plots are not rewritten.
To regenerate all production tables and plots, including the merged and split
DotPlots, from the existing annotated RDS without rewriting its labels, use:

```bash
Rscript scripts/xenium_vz_06_annotate_subclusters.R \
  --plots-only --resolution=0.8 --dry-run
sbatch --export=ALL,VZ_SUBCLUSTER_MODE=refresh-plots,VZ_SUBCLUSTER_RESOLUTION=0.8,VZ_SUBCLUSTER_OVERWRITE=true \
  scripts/run_xenium_vz_06_annotate_subclusters.slurm
```

Review job 40874352 completed in 19 minutes 50 seconds with 112.00 GB MaxRSS.
Plot refresh is expected to take approximately 30-90 minutes because it also
regenerates the production marker table and heatmap; this is a low-confidence
estimate until a refresh completes. The launcher requests 12 hours, 4 CPUs,
and 256 GB on `short`; the higher memory request follows the earlier 128-GB
segmentation failure at 125.73 GB MaxRSS.

Before mapping refined labels back to the whole-tissue objects, inspect all
expected inputs and outputs without loading Seurat objects:

```bash
Rscript scripts/xenium_vz_07_map_to_samples.R --dry-run
```

The mapping stage obtains exactly one input per sample from the Stage-4b
decision manifest. A `reuse_original` row points to the original consensus RDS;
a `filtered_copy` row points to a Stage-4 filtered object. File size and mtime
must still match the recorded decision. Legacy manifests remain supported and
retain the prior strict Stage-4-directory checks. The stage refuses to replace
existing outputs unless given `--overwrite`, verifies that every cell in the
master VZ object maps exactly once, and writes
`Xenium_VZ_Mapping_manifest.csv` with per-sample cell counts, mapped
reference-routed Unknown counts, and PCW. Per-sample mapped objects use the neutral name
`<sample>_VZ_QC_Subclusters.rds`.

Inspect the VZ plotting task map and one sample without creating plots:

```bash
Rscript scripts/xenium_vz_08_plot_samples.R --list
Rscript scripts/xenium_vz_08_plot_samples.R --dry-run 1
```

Stage 9 intentionally compares a focused nine-sample, one-sample-per-PCW
developmental cohort spanning every observed PCW from 10 through 19; it does
not use all 34 samples. PCW19 is represented by `GZFB_22_X_G_4`. Inspect its protected outputs with
`Rscript scripts/xenium_vz_09_plot_counts.R --dry-run` and submit it with
`run_xenium_vz_09_plot_counts.slurm`.

VZ and RL Stage-9 figures show percent abundance: each subcluster count is
divided by the total cells in that plotted lineage within the same sample,
then multiplied by 100. Each age bar totals 100%; these figures describe
lineage composition rather than abundance relative to whole tissue. A sample
with no cells in a plotted lineage is an error because its percentage is
undefined. The optional combined-lineage script uses the same denominator,
pooling cell counts across samples at each age before normalization.
Existing `09_cluster_counts` directories and `ClusterCountPlot` filenames are
retained for compatibility. Existing figures need regeneration; RDS objects
and downstream annotations remain valid. Overwrite protection is unchanged.
Expected Stage-9 runtime is 5–30 minutes per VZ or RL job, excluding queue
wait (low confidence, based on static inspection of one RDS load and 2–3
TIFF/PDF figure pairs); the launchers request 2 hours each. The path-only
`--dry-run` takes under one minute and does not load RDS objects.

The shared `refresh_current_plots.R` helper also displays cluster abundance
as a percentage of all labeled cells in each object/view, excluding unmapped
background cells. Raw counts remain in `Cluster_counts.csv`. Both legacy
`Counts_cells` and `Counts_percent_of_labeled` figure stems are retained and
now display the same percentages with explicit axis labels. Choose a new
snapshot tag to regenerate these views without replacing previous snapshots.

### 7. RL analysis

The RL branch mirrors the VZ branch:

1. `xenium_rl_01_subset_samples.R` – extracts configured consensus identities per sample.
2. `xenium_rl_02_merge_samples.R` – verifies and incrementally merges all manifest-defined sample subsets.
3. `xenium_rl_03_integrate.R` – validates the complete merge before processing.
4. `xenium_rl_04_review_qc.R` – 4a reviews QC; 4b records the decision and
   conditionally materializes filtered sample objects.
5. `xenium_rl_05_reintegrate_post_qc.R` – reuses Stage-3 results after a
   zero-removal decision, or reintegrates after removals.
   Before selecting the Stage-6 partition, run the protected Stage-5c
   [clustree and subsampling stability review](REGIONAL_RESOLUTION_REVIEW.md)
   with `run_xenium_rl_05c_resolution_stability.slurm`.
   To compare intermediate resolutions without rerunning integration, use the
   [RL 0.5/0.6/0.7/0.8 comparison](RL_RESOLUTION_COMPARISON.md), launched with
   `run_xenium_rl_05d_compare_resolutions.slurm`. Saved 0.5/0.8 assignments stay
   fixed; only 0.6/0.7 are computed on the saved SNN, and cluster 6 is displayed
   as identity unresolved in this separate diagnostic.
6. `xenium_rl_06_annotate_subclusters.R`
7. `xenium_rl_07_map_to_samples.R` – maps refined RL labels back to all manifest-defined configured samples.
8. `xenium_rl_plot_sample.R` and `xenium_rl_08_plot_samples.R` – manifest-defined per-sample reports.
9. `xenium_rl_09_plot_counts.R`

A separate read-only diagnostic plots the same RL-VZ marker panel across the
RL and VZ resolution-0.8 post-QC partitions. The expanded panel includes
RL-VZ/roof-plate core genes (`LMX1A`, `LMX1B`, `WLS`, `WNT1`, `PAX6`, `SOX2`,
`CALCB`, `CRYAB`, `COL2A1`, `MSX1`, `MSX2`, `GDF7`, `BMP4`, and `BMP7`), RL-VZ
proliferation, RL-SVZ/early-output, BG/radial-glia, canonical cerebellar
VZ/GABA, and ependymal/choroid programs. These matched positive and contrast
programs allow glial-like RL populations to be evaluated without changing a
graph, cluster assignment, or RDS. Genes absent from the active Xenium panel
are recorded in the availability table and omitted from the plot. Task 1 loads
RL and task 2 loads VZ, so the two large objects are handled independently:

```bash
Rscript scripts/xenium_vz_rl_diagnostic_rl_vz_marker_dotplot.R --list
Rscript scripts/xenium_vz_rl_diagnostic_rl_vz_marker_dotplot.R \
  --branch=RL --dry-run
Rscript scripts/xenium_vz_rl_diagnostic_rl_vz_marker_dotplot.R \
  --branch=VZ --dry-run
sbatch scripts/run_xenium_vz_rl_diagnostic_rl_vz_marker_dotplots.slurm
```

Matched TIFF/PDF DotPlots, marker-availability tables, cluster summaries, and
provenance manifests are written beneath
`outputs/xenium/vz_rl/diagnostics/rl_vz_marker_dotplots/`. Expected runtime is
approximately 15-45 minutes per task and 15-45 minutes total at the configured
array concurrency of two, excluding queue wait. This moderate-confidence
estimate is based on the prior 16-minute RL Stage-6 review, while allowing for
the larger VZ object and high-resolution rendering. The launcher requests two
hours, one CPU, and 256 GB per task on `express`.

A separate Purkinje diagnostic uses the existing Stage-6 annotated VZ object
(`06_subclusters/rds/Xenium_VZ_subclusters.rds`) and plots Maturing PCs,
Early-born PCs, Late-born PCs, and Patterning PCs, in that top-to-bottom order.
It uses the 15 fully legible markers from the user's cropped reference image:
VSTM2L, COL5A1, ITPR1, CNTNAP4, NDRG1, KCNC3, NEFL, NXPH4, PCDH10,
BCL11A, NPTX2, RFTN1, EBF2, EBF1, and RORB. Additional cropped markers
are not inferred. These are the current project identities; no PC1–PC4,
PC-TNC, or PC-MAB assignments are inferred from the reference.

From the project root on Cheaha:

```bash
Rscript --max-ppsize=500000 scripts/xenium_vz_diagnostic_purkinje_dotplot.R --dry-run
sbatch scripts/run_xenium_vz_diagnostic_purkinje_dotplot.slurm
```

The preflight checks paths and plotting dependencies without reading RDS
objects (expected under one minute). The single pooled plotting job is
expected to take 10–45 minutes, excluding queue wait; confidence is low,
based on loading the 1,447,477-cell VZ object and prior Stage-6 reviews of
31–46 minutes, with this job restricted to a small PC marker panel and no
marker tests. The launcher requests two hours, one CPU, and 256 GB on
`express`, matching the existing diagnostic resource envelope; prior VZ
Stage-6 review reached 152 GB. Measure this job's Elapsed and MaxRSS to
refine the estimate.

Outputs beneath `outputs/xenium/vz/diagnostics/purkinje_reference_markers/`
include PDF and 600-dpi TIFF plots, requested-marker availability, plotted
expression summaries, sample/source-cluster cell counts, and input provenance.
Genes absent from the Xenium panel are explicitly reported and omitted;
no available genes, missing expected PC identities, or non-finite expression
summaries cause an error. Dot size shows percent expressed; color uses the
project's grey-to-red gene-wise scaled mean expression across the four PC
groups, with fixed -2.5 to 2.5 limits. This is a pooled descriptive plot,
not a replicate-adjusted comparison. The script uses existing normalized
Xenium data and never saves a modified RDS or changes labels. Existing
diagnostic outputs are protected; after explicit review, replacement requires
`--overwrite` (launcher: `PURKINJE_DOTPLOT_OVERWRITE=true`).

#### Optional canonical iCN/eCN marker evidence (after the nonspatial merge)

`xenium_diagnostic_cn_markers.R` now defaults to the complete merged nonspatial
object at `outputs/xenium/vz_rl/nonspatial/03_merged/rds/Xenium_VZRL_merged.rds`.
It selects **eCN, iCN, MLIs, Golgi Cells, GABA Interneurons, and GABA
Progenitors**, in that display order, grouping only by existing `comb_subcluster` labels. VZ, RL, and Residual
cells are pooled within each label in the plot, expression means, coexpression,
and sample-level summaries. Other identities (including the separate Purkinje
lineage) are excluded from this report only. It requires iCN, eCN, and Residual-derived eCN,
all active manifest samples, and agreement between combined and source branch
labels. It does not read or modify the source regional RDS objects.

The combined-label stage does not propagate original numeric regional cluster
IDs. Therefore merged iCN pools former VZ clusters 11 and 21; the script
never interprets a whole-tissue clustering column as a regional cluster ID.
For a regional per-numeric-cluster view, explicitly select `--branch=vz`
(launcher: `CN_MARKER_BRANCH=vz`); `rl` and `residual` also remain available.

Submit the Slurm job from `scripts/` on Cheaha:

```bash
sbatch run_xenium_diagnostic_cn_markers.slurm
```

The dry-run checks paths and dependencies without loading RDS objects (expected
under one minute). Full merged execution is expected to take **15–90 minutes**,
excluding queue wait, with low confidence based on static inspection and
related workloads. Same-merge expression export job 41005351 took 3.45 minutes
and 22.86 GB MaxRSS; processed-merge plotting job 41006727 reached 24.88 minutes
and 182.34 GB before an unrelated marker-heatmap failure. Neither measures this
diagnostic. Optional regional reports retain the earlier 10–60-minute estimate
per branch. This is an expected runtime, not a requested
Slurm time limit. The launcher requests **two hours, one CPU, and 256 GB on
express**, matching the existing VZ diagnostic envelope. VZ Stage-6 job 40966278
took 31.52 minutes with 152.09 GB MaxRSS; this is supporting evidence, not a
measurement of the new diagnostic. Cheaha's published partition table confirms
express allows two hours (https://docs.rc.uab.edu/cheaha/hardware/, checked
2026-10-21). Loading the full object dominates memory. The launcher loads the
established R/4.4.1 module, starts R at the project root for existing renv
activation, exposes Seurat load errors, then runs the path-only dry-run before
the analysis. It never installs or restores packages. Default branch is
`nonspatial`; `CN_MARKER_BRANCH=vz`, `rl`, or `residual` explicitly selects an
regional report using the expanded marker panel. This is one pooled job, not an array.
After explicit output review, `CN_MARKER_OVERWRITE=true` enables the driver's
protected overwrite mode. Logs use `logs/Xen_CN_Markers_<job_id>.out` and `.err`.
After completion, collect `JobID`, `Elapsed`, `State`, `MaxRSS`, `ReqMem`,
`AllocCPUS`, and `ExitCode` with `sacct` to refine future resource requests.

Outputs are protected under
`outputs/xenium/vz_rl/nonspatial/diagnostics/canonical_cn_markers_gaba_ecn_mli/`: an unscaled PDF
DotPlot, gene availability, pooled and sample-by-subcluster expression,
equal-sample-weight expression means, within-cell marker-pair detection,
sample/cluster counts, manifest sample coverage, provenance, session information,
and interpretation notes with literature sources. The script reads the current
sample manifest, rather than assuming a fixed sample count. Missing samples in
the merged cohort are errors; regional-only reports record absent samples.
Missing metadata, unknown sample IDs, duplicate
cell IDs, incomplete layers, or conflicting cluster-to-label mappings are errors.
Genes absent from the panel are reported as unavailable, never zero expression.
An unavailable coexpression pair has `available=FALSE` and missing detection
values; other pairs and genes remain available for review. Merged result rows
do not split or group by `regional_assignment`. Sample coverage records source
and selected cell counts; provenance records the input/selected cell totals and
requested/evaluated labels. Absent non-CN requested labels are logged and omitted,
not displayed as zero-expression populations. Previous origin-separated and
regional reports remain in their existing directories and are never overwritten
by this focused report. Optional regional reports also use the new
`canonical_cn_markers_gaba_ecn_mli` directory name but retain their original numeric
cluster grouping and all regional labels.

The expanded panel requests 49 genes. It distinguishes SOX14/DMBX1 iCN
candidates; shared iCN-supporting ZFHX3, ZFHX4, TOX, SLC24A4, CHRM2, MEIS1,
and TAC1; GABA transmission/lineage genes; SLC17A6/MEIS2 eCN evidence;
subtype/developmental TBR1, LHX9, LMX1A, BARHL1, and NEUROD6; shared
eCN-supporting NTNG1, ZNF804A, NR4A2, NRP1, C1QL3, NTNG2, CREG2, and
NEUROD2; and MLI, GABA subtype context, and pan-neuronal genes.
The Purkinje contrast panel is replaced by PVALB, SORCS3, NXPH1, KCNA2, GRIK1,
KIT, ERBB4, GRM8, and CACNA1E. Genes moved out of GABA subtype context are not
duplicated across panels. These MLI markers have maturation/subtype dependence
and overlap with other GABA populations; SORCS3/NXPH1 transcriptional types
should not be equated directly with basket/stellate morphology. Sources:
https://www.nature.com/articles/s41586-021-03220-z and
https://www.nature.com/articles/s41467-022-30977-2. Added shared CN candidates are informed by the human developmental and
organoid comparison at https://doi.org/10.1242/dmm.052290 and developmental
subtypes by https://www.nature.com/articles/s41586-023-06884-x. These are not
all CN-specific or universally expressed: species, stage, and subtype matter.
Mouse evidence is distinguished from human developmental context in `interpretation.txt`.
SOX14 is not a universal requirement for all inhibitory CN subtypes. Generic
neuronal markers or GABA transmission alone do not establish CN identity.
Detection means raw counts greater than zero; colors show log1p of the mean
linearized existing normalized Xenium expression. No normalization is rerun.
The main `CN_markers_DotPlot.pdf` retains this unscaled, cell-pooled view.
A second `CN_identity_Focused_DotPlot.pdf` highlights 24 requested CN/GABA/MLI
genes, with eCN and iCN first, using equal-sample mean detection for dot size
and a gene-wise z-score of log1p(sample-balanced mean linear expression) for
color. Its six short panels are eCN core/support, iCN core/support, GABA, and
MLI. The MLI panel uses the same nine genes as the main plot: PVALB, SORCS3,
NXPH1, KCNA2, GRIK1, KIT, ERBB4, GRM8, and CACNA1E. The focused figure widens
with the number of available genes to preserve label readability.
Constant genes have z-score zero. This presentation highlights relative marker
enrichment and can amplify weak signals: inspect dot sizes and unscaled
values/sample counts in `CN_identity_Focused_Expression.csv` before drawing an
identity conclusion. Missing panel genes are omitted and recorded. The source
RDS, clustering, label definitions, and other analyses are unchanged.
Sample-balanced means include represented samples only, with no absent
sample-cluster zero imputation, minimum-cell threshold, or donor-level inference.
Inspect counts alongside these means, especially for small sample-cluster groups.
Same-cell detection can reflect segmentation contamination as well as biology.
These are descriptive identity checks, without significance tests or automatic
pass/fail assignments. Spatial nuclear localization and independent reference
evidence are needed for stronger identity claims. The only cell selection is
the requested report labels; QC thresholds, labels, cell IDs, and RDS objects
remain unchanged. No downstream analysis outputs become stale. Existing
report files require explicit review before rerunning with `--overwrite`.

The RL driver also reads all manifest-defined samples from `config/samples.csv` and supports
the same inspection and overwrite protections. `run_xenium_rl_01_subset.slurm`
submits RL extraction. Use `run_xenium_rl_02_merge_samples.slurm`,
`run_xenium_rl_03_integrate.slurm`,
`run_xenium_rl_05_reintegrate_post_qc.slurm`, and
`run_xenium_rl_06_annotate_subclusters.slurm` for the corresponding merged
stages.

RL outputs mirror the same numbered layout under `outputs/xenium/rl/`, from
`01_subsets` through `09_cluster_counts`. This makes corresponding VZ and RL
stages directly comparable with parallel, region-specific output filenames.

Each RL Stage-8 task uses the same back-to-front spatial rendering and writes
both all-inclusive and `Global_Spatial_NoOther` maps. Its principal faceted
TIFF places RL SVZ, Immature UBCs, and Mature UBCs in the first row, followed
by the five granule-lineage clusters in the second row. BG, Maturing PCs, eCN,
GABA Progenitors, Microglia, and Endothelial are written to a separate
Other-faceted TIFF. The per-sample `Spatial_Report.pdf` combines both global
maps, principal lineage facets, and Xenium marker DotPlot on one landscape
page. Standalone global and highlight TIFFs draw native 0.75-size ggplot2
points and use `_pt0.75` filenames. Both faceted files omit grey context cells,
draw only the requested cluster's cells at native size 0.6, and use `_pt0.6`
filenames. Global-plot grey context cells use size 0.5 and are drawn before the
larger size-0.75 colored cells. RL SVZ uses vivid magenta (`#FF4FA3`) and
Immature UBCs use vivid orange (`#FFA600`) to keep these adjacent lineage
colors visually distinct; the Immature UBC color is shared with the VZ
palette. Only the copies embedded in the merged PDF are pre-rasterized, which
avoids enlarging a lower-resolution point layer in the high-resolution TIFFs.

After all manifest-defined one-page Stage-8 reports for a branch are complete, the existing
RL or VZ Stage-8 driver can run in `--combine` mode. That mode validates the
reports and combines them into one manifest-sized PDF (currently 41 pages). Pages are ordered by numeric
PCW from `metadata/samples_meta.xlsx`, with `config/samples.csv` order breaking
ties, matching the established merged-report behavior.
It requires the Poppler `pdfinfo` and `pdfunite` commands, verifies that every
input has one page, verifies that the final output has one page per manifest sample, and writes a
page-to-sample CSV manifest containing PCW and original manifest position
beside the PDF. It never installs dependencies
inside the job. Inspect RL readiness and submit the same RL launcher in combine
mode after the RL array succeeds:

```bash
Rscript scripts/xenium_rl_08_plot_samples.R --combine --dry-run
rl_combined=$(sbatch --parsable \
  --array=1 \
  --dependency=afterok:${rl_stage8} \
  --export=ALL,RL_PLOTS_MODE=combine \
  scripts/run_xenium_rl_08_plot_samples.slurm)
```

The RL outputs are
`outputs/xenium/rl/08_sample_reports/Xenium_RL_AllSamples_Spatial_Report.pdf`
and its `_manifest.csv`. The `--array=1` override is required so only one array
task attempts the merge. For VZ, use
`xenium_vz_08_plot_samples.R --combine --dry-run` and submit
`run_xenium_vz_08_plot_samples.slurm` with `--array=1` and
`--export=ALL,VZ_PLOTS_MODE=combine`. Combining existing PDFs should take
approximately 1-5 minutes per branch. The shared Stage-8 launchers retain their
two-hour, 64-GB request because they also run the per-sample plotting tasks;
combine mode is expected to use far less time and memory.

Inspect RL merge readiness with:

```bash
Rscript scripts/xenium_rl_02_merge_samples.R --dry-run
```

The stable merged output is
`outputs/xenium/rl/02_merged/rds/Xenium_Merged_RLSubsets.rds`, accompanied by
the corresponding input, configured-label, reference-routed Unknown,
total-cell, and PCW
manifest under `02_merged/tables/`.
Both merge scripts refuse partial input sets, unexpected top-level RDS files,
and accidental output replacement.

Inspect processing readiness before stage 3:

```bash
Rscript scripts/xenium_rl_03_integrate.R --dry-run
```

Stage 3 writes the stable object
`outputs/xenium/rl/03_integrated/rds/Xenium_RL_integrated.rds`; RL QC and post-QC
processing consume that path.

As in VZ, the critical RL path defers the pre-Harmony diagnostic UMAP, allows
the corrected UMAP to use all allocated CPUs, and defaults clustering to one
process. The launcher maps `RL_INTEGRATE_CLUSTER_WORKERS` to
`XENIUM_STAGE3_CLUSTER_WORKERS`; leave it unset at 1 unless new memory evidence
supports increasing it. Stage 3 saves `Xenium_RL_precluster_checkpoint.rds`
before clustering and `Xenium_RL_postcluster_checkpoint.rds` immediately after
clustering, before plotting. Resume mode automatically uses the most advanced
available validated checkpoint. The RL script also performs the same complete
up-front source parse. Resume a failed job with:

```bash
Rscript scripts/xenium_rl_03_integrate.R --resume-checkpoint --dry-run
sbatch --export=ALL,RL_INTEGRATE_RESUME_CHECKPOINT=true,RL_INTEGRATE_OVERWRITE=true \
  scripts/run_xenium_rl_03_integrate.slurm
```

It also writes `Xenium_RL_UnknownClusters_BySample_UMAP.tif`, with sample IDs
appended to the independently numbered Unknown labels, and
`Xenium_RL_ConsensusLabel_UMAP.tif`, with every Unknown pooled into one grey
category. The sample-qualified Unknown legend is restricted to two columns on
a taller canvas. To regenerate the corrected stage-3 UMAPs from the existing
integrated object without rerunning normalization, Harmony, clustering, or
rewriting the RDS, inspect and submit plot-only mode with:

```bash
Rscript scripts/xenium_rl_03_integrate.R --plots-only --dry-run
sbatch --export=ALL,RL_INTEGRATE_PLOTS_ONLY=true,RL_INTEGRATE_OVERWRITE=true \
  scripts/run_xenium_rl_03_integrate.slurm
```

Generate the deferred pre-/post-Harmony batch comparison separately with:

```bash
Rscript scripts/xenium_rl_03_integrate.R --plots-only --include-uncorrected --dry-run
sbatch --export=ALL,RL_INTEGRATE_PLOTS_ONLY=true,RL_INTEGRATE_INCLUDE_UNCORRECTED=true,RL_INTEGRATE_OVERWRITE=true \
  scripts/run_xenium_rl_03_integrate.slurm
```

Until that deferred job completes, any older
`Xenium_RL_Batch_Comp_UMAP.tif` should be treated as stale; full stage 3 warns
when it detects that file but does not delete or overwrite it.

Marker testing and QC violins intentionally begin at stage 4, not stage 3;
those outputs are written under `outputs/xenium/rl/04_qc/plots/` and
`outputs/xenium/rl/04_qc/tables/`.

RL QC uses the same explicit Stage-4a review gate:

```bash
Rscript scripts/xenium_rl_04_review_qc.R --qc-only --dry-run
sbatch scripts/run_xenium_rl_04_review_qc.slurm
```

Its expanded review outputs mirror VZ under `outputs/xenium/rl/04_qc/`:
`Xenium_RL_SampleQualifiedUnknown_ClusterComposition.csv`,
`Xenium_RL_SampleQualifiedUnknown_QC_Summary.csv`,
`Xenium_RL_SampleQualifiedUnknown_QC_Violins.pdf`,
`RL_SampleQualifiedUnknown_top10_Markers.csv`, and
`RL_SampleQualifiedUnknown_Marker_Status.csv`. Marker comparisons use each
Unknown population against all other regional cells from the same sample, and
no removal is inferred from these descriptive outputs.

RL stage 4 also writes
`Xenium_RL_RawSubcluster_BySample_QC_Summary.csv`, containing count and feature
quartiles, within-sample median ratios, sample and cluster fractions, and
RL/Granule/UBC/Unknown composition for every observed sample/cluster pair.
Its raw-cluster marker table is named
`Xenium_RL_RawSubcluster_top10_Markers_Res0.8.csv`, matching the `Xenium_RL`
prefix used by the neighboring Stage-4 outputs. The equivalent VZ table is
`Xenium_VZ_RawSubcluster_top10_Markers_Res0.8.csv`.
`Xenium_RL_RawSubcluster_BySample_QC_Heatmaps.pdf` shows the log2 median count
and feature ratios relative to all RL-branch cells from the same sample, while
`Xenium_RL_RawSubcluster_SampleComposition_Heatmap.pdf` shows the fraction of
each expression cluster contributed by each sample. These descriptive outputs
are intended to distinguish reproducible low-QC clusters from sample- or
age-specific composition before any removal decision.

In Stage 4b, the user-approved 2026-10-16 resolution-0.8 decision removes only
clusters `19,20,21,22,23` from the current 41-sample RL object. Each contains
two cells: 10 of 827,193 cells total (approximately 0.0012%). All clusters
0-18 are retained. This supersedes the older `6,15,24,26,27` decision; cluster
numbers must not be carried between different integrated objects without review.
The launcher stores the current set as its apply-mode default; review remains
the default mode, so removal still requires `RL_QC_MODE=apply`. An explicit
`RL_QC_REMOVE_CLUSTERS` value can override the reviewed set if a later review
changes the decision. Do not place a comma-separated override directly inside
Slurm's `--export` list: Slurm treats those commas as separators between
environment variables. Export the variable in the shell first, or omit it to
use the reviewed launcher default. RL QC starts independently from the approved
whole-tissue consensus objects, not from VZ-mapped objects. Therefore VZ and RL
review/application can proceed independently after their respective stage 3
jobs complete. Inspect the decision and submit it with:

```bash
Rscript scripts/xenium_rl_04_review_qc.R \
  --remove-clusters=19,20,21,22,23 --dry-run
sbatch --export=ALL,RL_QC_MODE=apply,RL_QC_OVERWRITE=true \
  scripts/run_xenium_rl_04_review_qc.slurm
```

`RL_QC_OVERWRITE=true` is required when replacing the earlier no-removal
Stage-4b RDS files and decision manifest. It does not overwrite the Stage-4a
QC evidence. The approved removal contains 10 cells in the reviewed object;
Stage 5 therefore takes the full reintegration path, and all
existing Stage-5 and downstream RL results become stale. As in the VZ Stage-4
launcher, RL Stage 4 starts R with `--max-ppsize=500000`. RL materialization is
also deliberately sequential so each spatial Seurat object is subset in that
main process. Jobs 40698160 and 40698252 showed that the command-line stack
setting was insufficient for the `multisession` path: both failed with
`protect(): protection stack overflow` below 39 GB MaxRSS. Job 40698252 also
demonstrated the Slurm comma-parsing hazard by receiving only cluster `6` from
a comma-separated `--export` override.

The launcher defaults to Stage 4a review mode and requests 12 hours, 1 CPU,
and 256 GB on `short`. Stage 4b records one consistent decision for all 41
samples. A `none` decision writes only a routing manifest pointing to the
original consensus objects and should finish in under 1 minute. A nonzero
removal writes 41 filtered copies sequentially and is expected to take
approximately 25-60 minutes (low confidence), extrapolated from job 40698364,
which completed the earlier 34-sample cohort in 24 minutes 46 seconds with
36.83 GB MaxRSS. Timing the first 41-sample run will improve this estimate.
Expanded Stage 4a review remains approximately 2-8 hours. These estimates
exclude queue wait. The
12-hour limit and 256-GB memory request remain conservative pending additional
successful runs.

Inspect the protected post-QC and subcluster stages with:

```bash
Rscript scripts/xenium_rl_05_reintegrate_post_qc.R --dry-run
Rscript scripts/xenium_rl_06_annotate_subclusters.R --review-only --dry-run
```

When a nonzero Stage-4 removal makes the Stage-3 and Stage-5 UMAPs look
substantially different, run the protected diagnostic rather than replacing
the reviewed QC decision or the Stage-5 object. It performs four controlled
checks: (1) plots only retained cells on the unchanged Stage-3 UMAP, (2)
highlights reviewed removal clusters 6 and 15 separately on that same UMAP,
(3) exports their pooled-label, sample, QC, and existing Stage-4 top-10 marker
evidence, and (4) uses the exact Stage-5 retained cells and existing
`harmony_clean` coordinates to compare dimensions 1-30 with `k.param=30`
versus `k.param=50`. The controlled comparison keeps UMAP neighbors at 30,
the random seed at the configured value, and resolutions at 0.3, 0.5, and
0.8, so neighbor `k` is the only graph parameter that changes. It writes no
RDS and does not change the removal manifest, Stage-3 object, Stage-5 object,
or Stage-6 review outputs.

Inspect and submit it with:

```bash
Rscript scripts/xenium_rl_05a_diagnose_post_qc.R --dry-run
sbatch scripts/run_xenium_rl_05a_diagnose_post_qc.slurm
```

Outputs are protected under
`outputs/xenium/rl/05_post_qc/diagnostics/{plots,tables}/`. If a failed or
previously completed diagnostic left final outputs that you explicitly intend
to replace, submit with
`RL_POST_QC_DIAGNOSTIC_OVERWRITE=true`. Partial results from a failed job stay
inside a job-specific hidden staging directory and are never treated as final
diagnostic outputs. Based on completed RL Stage-5 job 40699023 (3 hours 32
minutes, 136.39 GB MaxRSS for the larger full reintegration) and Stage-6 job
40792498 (16 minutes, 69.09 GB MaxRSS), this diagnostic is expected to take
approximately 1.5-4 hours with four CPUs and remain below 256 GB. This is a
moderate-confidence estimate; it excludes queue wait. The launcher requests
12 hours on `short` to provide margin for two neighbor graphs and six matched
clusterings.

To test whether smaller expression neighborhoods separate the glial-like RL
population without changing any production object, run the protected Stage-5b
sensitivity analysis. It holds the existing 30-dimensional Harmony embedding,
`umap_clean` coordinates, SNN pruning, Louvain algorithm, marker thresholds,
and random seed fixed. It rebuilds graphs with `k.param=10` and `20` and
clusters each graph at resolutions 0.8 and 1.0. Outputs include matched UMAPs,
production-cluster-6 highlight plots, cell-level assignments, cluster/sample/
input-label composition tables, cross-parameter contingencies, top-50 marker
tables, and focused RL DotPlots. Results are staged before publication and no
RDS, QC decision, production graph, or Stage-6 review output is changed.

Inspect and submit with:

```bash
Rscript scripts/xenium_rl_05b_k_resolution_sensitivity.R --dry-run
sbatch scripts/run_xenium_rl_05b_k_resolution_sensitivity.slurm
```

Outputs are written beneath
`outputs/xenium/rl/05_post_qc/diagnostics/k_resolution_sensitivity/`. To replace
an earlier completed sensitivity run after review, submit with
`RL_K_RES_SENSITIVITY_OVERWRITE=true`. Expected runtime is approximately
45-120 minutes with four CPUs and 256 GB, excluding queue wait. This estimate
is moderate-to-low confidence: it is based on the 16-minute single-partition
Stage-6 marker review plus two graph builds and four marker/DotPlot reviews,
but this exact four-partition diagnostic has not yet been measured. The Slurm
launcher requests four hours on `short`.

When the complete stage-4 removal manifest records `clusters_removed=none`,
zero removed cells, and unchanged input/remaining counts for all manifest-defined samples,
RL stage 5 takes a protected zero-removal fast path. It preserves the stage-3
PCA, Harmony integration, SNN graph, and resolution-0.5/0.8 assignments; copies
the reviewed `umap_harmony` coordinates to the stage-5 `umap_clean` reduction;
and writes only the resolution-0.5/0.8 fast-path plots. It does not calculate a
new resolution-0.3 partition. Any nonzero removal uses the existing full
PCA/Harmony/UMAP/neighbors/reclustering path at resolutions 0.3/0.5/0.8. After
the controlled post-QC diagnostic, the full path uses 30 PCA/Harmony
dimensions, 30 UMAP neighbors, and `k.param=30`. Resolution 0.8 is the
user-selected default Stage-6 review partition; resolutions 0.3 and 0.5 remain
available through explicit overrides for sensitivity checks. The Stage-5 provenance manifest records these
parameters, the selected processing mode, and which reductions and resolutions
were reused or computed. The original-label UMAP uses
`regional_consensus_label_pooled`, preventing sample-qualified Unknown labels
from creating an oversized legend that compresses the plotting panel.

The zero-removal fast path is expected to take approximately 10-45 minutes,
dominated by loading, joining, plotting, and saving the large RL object. Full
post-removal reintegration for the current 41-sample cohort is estimated at
3-8 hours (low confidence). Job 40958628 reached clustering after approximately
1 hour 40 minutes but was killed by the 256-GB cgroup limit when four
multisession workers started; it did not measure clustering completion time.
The earlier 34-sample job 40699023 finished in 3 hours 32 minutes. Estimates
exclude queue wait. The
launcher requests 24 hours, 1 CPU, and 600 GB on `medium`; clustering defaults to one worker.
The 2026-10-16 fix replaces multisession object serialization with fork-based
clustering and, at the user's request, defaults to three concurrent clustering
workers so the three resolutions can run together. Override with
`RL_POST_QC_CLUSTER_WORKERS=1` for sequential clustering. Worker count is capped
at allocated CPUs and the three resolutions. Unsupported fork environments
fail before loading the large object and instruct the user to select one worker.
This runtime change preserves the configured dimensions, neighbors, resolutions,
seed, cell selection, output names, and overwrite protection. The provenance
manifest now records the clustering worker count and backend. The 600-GB request
provides a conservative retry budget after the 256-GB OOM; it is not a measured
requirement or a guarantee of completion. The failed job reported 173.50 GiB
MaxRSS for an individual process, which does not measure concurrent job memory.
Successful peak memory remains unmeasured; collect sacct after the retry to
refine the request. Current Cheaha hardware documentation, checked 2026-10-16,
lists Intel short-partition nodes with up to 754 GB and a 12-hour limit:
https://docs.rc.uab.edu/cheaha/hardware/ . The request stays on `short`, retaining
the Intel constraint and existing 12-hour limit. The 3-8-hour estimate remains
low confidence because neither the two- nor three-fork-worker retry has completed.
There is no RL Stage-5 restart checkpoint, so this retry repeats PCA, Harmony,
UMAP, and neighbor finding before clustering.
Historically, the memory request was raised from 128 to 256 GB after full-reintegration job
40698736 was killed while loading or validating the merged RL object at
124.25 GB MaxRSS, before any PCA, Harmony, UMAP, neighbor, or clustering
progress was logged. This resource-only change does not alter filtering,
integration, clustering, or output parameters.

After confirming the Stage-5 dry run, explicitly replace the stale Stage-5
outputs from the earlier parameterization with:

```bash
RL_POST_QC_OVERWRITE=true \
  sbatch --export=ALL scripts/run_xenium_rl_05_reintegrate_post_qc.slurm
```

For uncertain resolution-0.8 identities after Stage-6 review, run the read-only
identity diagnostics before applying labels:

```bash
Rscript --max-ppsize=500000 scripts/xenium_rl_06b_diagnose_identities.R --dry-run
sbatch scripts/run_xenium_rl_06b_diagnose_identities.slurm
```

The diagnostic checks the current Stage-6 review provenance against the Stage-5
RDS and preserves every cell, sample, cluster, and original prediction label.
It summarizes cluster sizes/counts/features, within-sample QC ratios, and
within-cell raw-count detection for targeted marker pairs. Detection means
counts greater than zero, a descriptive rule rather than a QC or identity
threshold. Missing panel genes are recorded as unavailable. All clusters receive
aggregate diagnostics; clusters 5, 6, 10, 12, 18, and 19 also receive compressed
cell-level marker-count and reference-prediction CSVs when present.
Original Aldinger, Sepp, and Science annotation objects are read sequentially
for every manifest sample. Each retained RL cell must match exactly once using
the Stage-2 sample prefix; missing/duplicate IDs, blank labels, invalid scores,
and partial sample/reference coverage are errors. Saved `predicted.id` and
`prediction.score.max` values are reported without confidence filtering,
recalibration, label transfer, consensus voting, or changes to analysis RDS files.
Scores should be interpreted within each reference. Label-composition tables
are descriptive cell fractions, not independent-donor statistical inference.

Outputs are protected under
`outputs/xenium/rl/06_subclusters/diagnostics/identity_review/{tables,plots}/`.
They include three PDF summaries, complete per-cluster/per-sample CSVs, the
focus-cell CSVs, an exact reference-join/input audit, and provenance. Existing
diagnostic outputs require explicit `--overwrite` after review, exposed through
`RL_IDENTITY_DIAGNOSTICS_OVERWRITE=true` in the launcher. The path-only preflight
checks the post-QC object, review manifest, and all reference paths without
loading RDS files; it should take under 1 minute. Expected full runtime is
1-6 hours (low confidence from static inspection of one approximately 43-GB
RL-object load and 123 sequential reference loads for 41 samples; measure the
first completed run to refine this). The job requests 12 hours, one CPU, 256 GB,
and Intel nodes on `short`; queue wait is excluded. Memory is budgeted above
the completed RL review's 80.32-GiB peak, with only one reference object resident
at a time after the large RL object is released.

#### Optional BG spatial and within-cell marker diagnostics (after RL Stage 7)

`xenium_rl_06d_bg_spatial_markers.R` reads the same manifest-defined, whole-tissue
Stage-7 mapped objects used by the sample spatial reports. It analyzes the
existing **BG Progenitors** label without renaming cells, filtering, reclustering,
or writing analysis RDS files. Run this on the HPC. From the project root:

```bash
cd /home/acflint/R/Projects/XeniumFCProject
# Lists manifest tasks and checks every input path; no RDS loads or submission.
bash scripts/submit_xenium_rl_bg_spatial_markers.sh
# Repeat the preflight, then submit the manifest-derived array at concurrency 3.
bash scripts/submit_xenium_rl_bg_spatial_markers.sh --submit
```

Upload the driver, `R/rl_bg_spatial_diagnostics.R`, launcher
`run_xenium_rl_06d_bg_spatial_markers.slurm`, and submission helper together.
The helper preserves the existing R/4.4.1 module setup. There are currently 41
manifest samples; task IDs follow manifest order, not PDF page order. Inspect
`Rscript --max-ppsize=500000 scripts/xenium_rl_06d_bg_spatial_markers.R --list`
before selecting a single task. A single task accepts `--dry-run TASK_ID` or
`TASK_ID`. The all-sample `--dry-run` checks paths only and returns nonzero when
inputs are missing. Preflight/listing should take under 1 minute.

Outputs are separate, protected diagnostic files under
`outputs/xenium/rl/06_subclusters/diagnostics/bg_spatial_markers/<sample_id>/`:

- A nine-page PDF: BG-only map, RL context map, six EOMES-pair maps
  (SLC1A3/BCAN/SOX2/NES/PTPRZ1/FABP7), and raw-count marker detection across
  all mapped RL identities. All cells are displayed; there is no random sampling.
- Compressed cell-level counts, centroids, IDs, labels, age and QC metadata for
  every mapped RL cell, plus marker availability/detection and pair tables.
  Pair tables report the four detection states, fractions of measured cells,
  conditional fractions, and denominators at both >=1 and >=2 raw transcripts.
  These are descriptive sensitivity checks, not identity or QC thresholds.
- Per-BG-cell nearest-centroid distances to RL SVZ, UBCs (immature + mature),
  and granule identities, with per-sample medians. Coordinates use the report's
  x/y swap; Euclidean distance is unchanged. Units remain native Xenium centroid
  units and are not independently calibrated to microns here.
- A sample summary and a completion/provenance CSV written last. The input's
  size/mtime, configuration/manifest hashes, helper hash and package versions
  are recorded. Partial outputs remain protected until reviewed.

Sample identity is checked against the Stage-7 RL mapping manifest's exact
sample/output-path match, whole-tissue and mapped-RL cell counts, and PCW.
The original per-sample `orig.ident` may be a legacy Seurat project identifier;
it is preserved in the cell table and provenance rather than rewritten or
required to equal the manifest sample ID. Multiple/blank origin identifiers,
an origin naming a different active sample, or conflicting explicit `sample_id`
metadata remain errors. Preflight reads paths and the small mapping CSV only,
with no RDS loads. The submission helper counts numeric task rows only, so
startup messages such as a `renv` status notice cannot create extra tasks.

NA RL labels are expected whole-tissue background, not missing RL assignments.
All object cells must have exactly one finite centroid; mapped RL cells must
have exactly one raw-count layer. Missing genes are unavailable (NA), never
zero expression. A sample without BG cells is explicitly `BG_absent`, and an
absent distance-comparison identity produces NA distances with an explicit status.
The cohort summary rejects an entirely absent BG label. These diagnostics do
not inspect transcript-to-boundary assignments, establish lineage, or determine
whether proximity reflects anatomy, cell density, or section geometry.
Samples may include sections from the same donor; no cell-level or sample-level
inferential p-values are calculated.

Expected runtime is **5-30 minutes per task**, or approximately **1-7 hours for
41 samples at concurrency 3**, excluding queue wait (low confidence from static
inspection: one mapped object load, targeted counts, nearest-neighbor searches
and nine plots). The launcher requests **2 hours, 64 GB and one CPU per task
on express**, following the existing Stage-8 per-sample launcher; this new job
has no measured memory peak yet. The historical per-sample transfer jobs peaked
around 44-49 GB, but that is only an indirect comparator. Collect first-run
`sacct` JobID, Elapsed, State, MaxRSS, ReqMem, AllocCPUS and ExitCode to refine
these estimates. No major partition/resource change is introduced.

After every array task completes, combine the small tables (normally under
1 minute; no RDS objects loaded):

```bash
Rscript --max-ppsize=500000 scripts/xenium_rl_06d_bg_spatial_markers.R --summarize --dry-run
Rscript --max-ppsize=500000 scripts/xenium_rl_06d_bg_spatial_markers.R --summarize
```

Combination requires complete outputs for every current manifest sample and
matching input/config/manifest/helper provenance; it preserves sample rows
rather than pooling cells for donor inference. Existing diagnostic outputs
require explicit `--overwrite` after review. For arrays only, the equivalent
launcher opt-in is `RL_BG_DIAGNOSTICS_OVERWRITE=true`. No existing production
outputs become stale from running these read-only diagnostics.

After stage 5, submit `run_xenium_rl_06_annotate_subclusters.slurm` in its
default review mode. Resolution-0.8 review generates the five sample-aware
marker-evidence tables within the same job before generating the review plots
and descriptive marker tables. Complete evidence with matching input provenance
is reused. Partial or stale evidence is protected and requires explicit
`--overwrite` after review; label application still requires completed review
outputs and current evidence. A missing evidence set is an output to build,
not a failed prerequisite in review-mode `--dry-run`. The preflight reads paths
and small provenance CSVs only and takes under 1 minute.
The evidence builder releases its object before the review loads it, avoiding
two simultaneously resident copies. Existing stable output paths are preserved.
The combined first RL review is estimated at 0.75–7 hours (low confidence, combining
the 20-120-minute evidence estimate with the earlier marker-review timing and
allowance for the 41-sample object); the launcher retains its 12-hour limit,
one CPU, and 256 GB. A separate RL 06a job is no longer needed.
Both the R driver and launcher default to resolution 0.8, while
`--resolution=0.3|0.5|0.8` and the launcher's
`RL_SUBCLUSTER_RESOLUTION` select another Stage-5 partition. Each review writes
a resolution-specific raw UMAP, stable top-10 marker table, expanded top-50
marker table, provenance manifest, focused RL lineage-marker DotPlot in TIFF
and PDF, and marker-availability CSV. The focused panel covers RL VZ/RL SVZ,
UBC, GCP/GC maturation, eCN, glial, Purkinje/GABAergic, OPC/ependymal, immune,
and endothelial programs. Missing Xenium-panel genes are recorded and omitted
from the plot. The shared `rl_markers` list from `color_palette.R` is also
included; overlapping genes retain their focused biological facet and only
palette-only additions receive a separate facet, preventing duplicate dots.
Marker testing retains the same scientific thresholds, and the
UMAP plot reuses `umap_clean` coordinates. Earlier resolution-specific evidence
from a previous Stage-5 object is stale and must be explicitly overwritten.
Label approval is separate from evidence generation. Review now writes a blank
`Xenium_RL_Subcluster_Label_Approval_Template_Res0.8.csv` under the Stage-6 tables
directory. After reviewing the current fit, fill its labels, notes and approval
flags and save a copy as `scripts/reviews/RL_subcluster_labels_Res0.8.csv`.
Do not edit its source/review/partition fingerprints. Apply requires this
explicit decision and checks it against both the fitted cell-to-cluster
assignments and the configured `reviewed_label_maps`; matching integer cluster
IDs alone cannot authorize reuse. The small-table preflight stays under one
minute and loads no RDS. Full apply fingerprints the already loaded assignments
(estimated under 10 seconds extra at ~827k cells, low confidence); plot-only
refresh is unaffected. Synthetic approval checks take under 10 seconds.
The existing user-approved 2026-10-17 resolution-0.8 map remains in the script
for the historical 827,183-cell object's clusters 0-19:
`0/2/10/19` Differentiating GCs; `1/18` Mature UBCs; `3/5` Immature UBCs;
`4/13` Migrating GCs; `6` BG Progenitors; `7/9/14/15` Prolif GCPs;
`8/17` RL SVZ; `11` Mature GCs; `12` Maturing PCs; `16` Maturing GCPs.
This is a marker-based annotation approved by the user, with sparse 18/19
assigned tentatively using the UMAP and focused DotPlot. Assignments for
5/6/10/12 remain biologically provisional despite approval of the label names.
BG Progenitors is registered in RL and master orders with the existing BG
violet color and unchanged BG marker panel; VZ BG labels are unchanged.
All retained cells, including 18/19, remain in the annotated RDS. The lineage-only
UMAP excludes BG Progenitors and Maturing PCs for display only. Shared taxonomy and palette checks
prevent any reviewed label from becoming `NA` during factor conversion. The
combined count plots report `eCN` as its own lineage and place `Microglia` with
the other non-neural/off-lineage populations rather than counting either as
UBC. Production RL-branch plots use only the 10 identities present in that
approved map; the RL-branch's obsolete `RL VZ`, `Intermediate Progenitors`,
and `Ependymal Cells` levels are not displayed. This does not remove the
separately reviewed `RL VZ` identity from the VZ resolution-0.8 taxonomy. The
production output retains the merged DotPlot and adds separate compact RL
lineage and Other DotPlots. The lineage panel uses three RL and three UBC
markers plus the six user-selected granule genes (see Shared DotPlot
presentation). The Other panel displays only the marker groups for observed
identities: currently BG Progenitors (SLC1A3/SOX2/BCAN) and Maturing PCs
(CALB1/FOXP4/DAB1). Inactive eCN, GABA Progenitors, Microglia, and Endothelial
facets are omitted. Nonspecific PAX6 and assay-absent GRIA4 remain omitted
from these production panels.
Both RL UMAPs use centered white boxes for every cluster label, without leader
lines or changes to UMAP coordinates or identities.

Applying this map requires no Stage-5 reintegration or Stage-6 marker-review
rerun. It changes the annotated RL RDS, label manifest, and label-dependent
Stage-6 plots/tables, making earlier mapped per-sample RL and combined-label
outputs stale. Preserve prior outputs/label CSVs for rollback and regenerate
downstream stages only after reviewing the new annotation. Existing outputs
remain protected unless `--overwrite` is explicitly supplied. Label application
is estimated at 10-30 minutes (low confidence based on the earlier 11m32s job
40865299 and the larger current object), excluding queue wait; the launcher
retains its 12-hour limit and 256-GB request.

```bash
Rscript scripts/xenium_rl_06_annotate_subclusters.R \
  --review-only --resolution=0.8 --dry-run
sbatch --export=ALL,RL_SUBCLUSTER_MODE=review,RL_SUBCLUSTER_RESOLUTION=0.8,RL_SUBCLUSTER_OVERWRITE=true \
  scripts/run_xenium_rl_06_annotate_subclusters.slurm
```

Inspect and submit the reviewed resolution-0.8 application with:

```bash
Rscript scripts/xenium_rl_06_annotate_subclusters.R \
  --apply-reviewed-labels --resolution=0.8 --dry-run
sbatch --export=ALL,RL_SUBCLUSTER_MODE=apply,RL_SUBCLUSTER_RESOLUTION=0.8 \
  scripts/run_xenium_rl_06_annotate_subclusters.slurm
```

Application is refused if the post-QC input changed after review. Review mode
with already current replicate evidence is estimated at 30-120 minutes for
the 41-sample object (low confidence). The earlier resolution-0.8 marker-review
job 40792498 completed in 16 minutes 9 seconds with 69.09 GB MaxRSS; it did not
include the current combined evidence build. Allow 1-4 hours for the first
combined review, as described above. Full resolution-0.8 annotation is
expected to take approximately 10-20 minutes based on completed job 40865299,
which ran in 11 minutes 32 seconds with 102.77 GB MaxRSS. The launcher requests
12 hours, 1 CPU, and 256 GB on `short`.

To refresh only the Stage-6 production tables and plots from the existing
annotated RDS, without rewriting the label manifest or RDS, inspect and then
submit plot-only mode:

```bash
Rscript scripts/xenium_rl_06_annotate_subclusters.R \
  --plots-only --resolution=0.8 --dry-run
sbatch --export=ALL,RL_SUBCLUSTER_MODE=refresh-plots,RL_SUBCLUSTER_RESOLUTION=0.8,RL_SUBCLUSTER_OVERWRITE=true \
  scripts/run_xenium_rl_06_annotate_subclusters.slurm
```

Plot-only refresh is expected to take approximately 10-20 minutes, with
moderate confidence based on completed full Stage-6 job 40865299. It reuses the
launcher's 12-hour, 4-CPU, 128-GB request and does not rewrite
`Xenium_RL_subclusters.rds`.

After reviewed RL labels exist, inspect the complete mapping inputs and
protected outputs:

```bash
Rscript scripts/xenium_rl_07_map_to_samples.R --dry-run
```

The RL mapping stage obtains its inputs from the Stage-4b routing manifest,
validates their recorded file provenance, and applies the same completeness,
overwrite, barcode-matching, and manifest checks as VZ. Legacy manifests still
route to existing Stage-4 `*_RL_QC.rds` files. It writes independent
`<sample>_RL_QC_Subclusters.rds` objects and
`Xenium_RL_Mapping_manifest.csv`; VZ labels are joined later by exact cell ID.

Use `run_xenium_rl_07_map_to_samples.slurm` for mapping and
`run_xenium_rl_08_plot_samples.slurm` for the 41-task plotting array at concurrency
three, matching the current manifest. The sample plotting helper includes
BG Progenitors in Other facets and uses the six selected granule markers and
broad-colored DotPlot headers. Merged and per-sample production DotPlots share
`rl_production_marker_panel()`: only marker groups for observed identities are
shown. Samples with the same identities as the merged object therefore use the
same ordered markers; absent identity groups are omitted rather than retaining
the old generic Other markers. Raw review markers/signatures are unchanged.
Inspect one
plotting task first with:

```bash
Rscript scripts/xenium_rl_08_plot_samples.R --dry-run 1
```

After the sample array completes, inspect and submit the single-task combine:

```bash
Rscript --max-ppsize=500000 scripts/xenium_rl_08_plot_samples.R --combine --dry-run
sbatch --array=1 --export=ALL,RL_PLOTS_MODE=combine scripts/run_xenium_rl_08_plot_samples.slurm
```

Stage-8 sample reports are estimated at 5-20 minutes per task and approximately
70-280 minutes for 41 tasks at concurrency three, excluding queue wait (low
confidence from static inspection of per-sample RDS reads, spatial TIFF/PDF
rendering, and DotPlots; no matching sacct measurement yet). Each task requests
2 hours, one CPU, and 64 GB on express. Path-only sample preflights take under
1 minute each. Combining the 41 one-page reports is estimated at 1-5 minutes,
with the same 2-hour requested limit; its preflight checks PDF files/page counts
and the PCW page map without loading RDS objects and should take under 1 minute.
Expected combined coverage is 41 pages in current manifest/PCW order. Existing
sample and combined outputs remain protected by the launcher's RL_PLOTS_OVERWRITE
switch. No analysis RDS files or mapped labels are modified by Stage 8.

Optional RL UMAP comparison: `xenium_rl_06c_compare_umap.R` reads the approved
Stage-6 `Xenium_RL_subclusters.rds` and fits four visualization alternatives
using its existing `harmony_clean` coordinates: neighbors 30/15 crossed with
minimum distance 0.1/0.05. Seed and dimensions come from configuration; spread
remains 1 and the metric is cosine. It never changes clustering, approved labels,
or production objects. All cells are fitted; a fixed display subset of up to
2,000 cells per identity reduces crowding in the five-page PDF, whose first
page shows the existing embedding. Full coordinates, display cell IDs, parameters,
input provenance, and session information are saved under
`outputs/xenium/rl/06c_umap_comparison/`. Existing comparison outputs are protected
unless explicitly replaced with `--overwrite`.

```bash
Rscript --max-ppsize=500000 scripts/xenium_rl_06c_compare_umap.R --list
Rscript --max-ppsize=500000 scripts/xenium_rl_06c_compare_umap.R --dry-run
sbatch scripts/run_xenium_rl_06c_compare_umap.slurm
```

List/path-only preflight should take under one minute and do not load RDS data.
The four sequential fits are estimated at 3-7 hours total, excluding queue wait,
based on the earlier approximately 42-minute UMAP fit for 827,183 cells plus
loading/rendering overhead (low confidence). The launcher requests 8 hours,
256 GB and one CPU on short; no array or parallel R workers are used. Review the
PDF before adopting any alternative; production plots remain current because
this optional comparison does not replace their coordinates. After completion,
collect Elapsed, State, MaxRSS, ReqMem, AllocCPUS and ExitCode with sacct.

User-selected RL UMAP (2026-10-20): 30 neighbors, minimum distance 0.05, spread 1.
The RL-only `regional_subsets.rl_umap_min_dist_post_qc` setting controls future
Stage-5 full recomputations; it does not alter clustering or the zero-removal
reuse path. Existing production RDS reductions are retained. To display every
cell from the already computed `nn30_md0.05` comparison, use:

```bash
Rscript --max-ppsize=500000 scripts/xenium_rl_06c_compare_umap.R --render-selected --dry-run
sbatch --export=ALL,RL_UMAP_OVERWRITE=true scripts/run_xenium_rl_06c_render_selected_umap.slurm
```

The selected-render mode checks source provenance and renders all saved cell
coordinates using the original Stage-6 DimPlot formatting. It loads the approved
Seurat object, substitutes coordinates in memory, and never fits UMAP again. It
replaces the two stable Stage-6 full/lineage-only UMAP TIFFs, adds a two-page
`Xenium_RL_Subcluster_UMAP_AllCells.pdf`, and records a selected-UMAP manifest.
The full plot includes every cell; the companion filtered plot retains the
established lineage-only exclusion. Estimated runtime is 10-25 minutes including
the large object read (low confidence, based on earlier Stage-6 plot refresh timings). The launcher
requests 2 hours, 256 GB and one CPU on express. Preflight is path-only, under one minute.
Subsequent Stage-6 plot-only refreshes honor the selected manifest and redraw
these UMAPs from the same comparison. Changed source-object provenance requires
a new reviewed comparison. Stage-7 mapped labels, Stage-8 spatial reports, and
marker results remain valid and do not need regeneration. No existing Seurat
RDS is overwritten; their stored UMAP coordinates are unchanged.

To refresh only the five upstream RL UMAPs with the same selected coordinates:

```bash
Rscript --max-ppsize=500000 scripts/xenium_rl_06c_compare_umap.R --render-selected --upstream-only --dry-run
sbatch --export=ALL,RL_UMAP_SCOPE=upstream,RL_UMAP_OVERWRITE=true scripts/run_xenium_rl_06c_render_selected_umap.slurm
```

This replaces the four Stage-5 UMAP TIFFs (pooled input labels and resolutions
0.3/0.5/0.8) and the Stage-6 raw resolution-0.8 UMAP TIFF. It reads the approved
annotated object, which retains the upstream grouping columns, substitutes the
selected full-cell coordinates in memory, and uses the original Stage-5/review
plot settings verbatim. No marker tests, cluster assignments, labels, graphs,
RDS objects, production annotated UMAPs, or other plots are changed. A separate
`Xenium_RL_Upstream_UMAP_refresh_manifest.csv` records all five output paths.
Expected runtime is 15-35 minutes, excluding queue wait (low confidence, based
on earlier Stage-6 object loading/rendering times); requested resources remain
2 hours, 256 GB and one CPU. Path-only preflight takes under one minute.

Stages 5 and 6 write the stable objects `Xenium_RL_postQC.rds` and
`Xenium_RL_subclusters.rds` in their numbered stage directories. Stage
9 intentionally uses the same focused nine-sample, one-sample-per-PCW cohort
as VZ, including `GZFB_22_X_G_4` for PCW19. Inspect it
with `Rscript scripts/xenium_rl_09_plot_counts.R --dry-run` and submit it with
`run_xenium_rl_09_plot_counts.slurm`.

### 7a. Residual Other/unresolved branch

Residual post-QC processing now matches VZ/RL: refit after approved removals,
or reuse the original analysis after an explicit no-removal decision. See
[RESIDUAL_POST_QC_PROCESSING.md](RESIDUAL_POST_QC_PROCESSING.md) for the new stage,
fresh label-review requirement, protected outputs and validation. Expected
runtime is 15–120 minutes after removals or 1–10 minutes for reuse (low confidence,
~10,000 cells; no array), versus the inherited 12-hour requested limit.


Retained high-quality Unknown populations assigned `route_other` or
`route_unresolved` enter an independent residual branch instead of being
forced into VZ or RL. This branch is deliberately non-spatial through label
review; spatial plots remain subsequent validation.

```bash
Rscript scripts/xenium_residual_01_subset_samples.R --list
Rscript scripts/xenium_residual_01_subset_samples.R --dry-run 1
sbatch scripts/run_xenium_residual_01_subset_samples.slurm

Rscript scripts/xenium_residual_02_merge_samples.R --dry-run
sbatch scripts/run_xenium_residual_02_merge_samples.slurm

Rscript scripts/xenium_residual_03_integrate.R --dry-run
sbatch scripts/run_xenium_residual_03_integrate.slurm

# Independent QC evidence; can run after Stage 3 without rerunning markers.
Rscript --max-ppsize=500000 scripts/xenium_residual_04a_review_qc.R --dry-run
sbatch scripts/run_xenium_residual_04a_review_qc.slurm

# Run marker review after QC decisions are applied, if any are approved.
Rscript --max-ppsize=500000 scripts/xenium_residual_04_review_subclusters.R --dry-run
sbatch scripts/run_xenium_residual_04_review_subclusters.slurm
```

Stage 4a adds descriptive QC before label approval. It reads the integrated
resolution-0.8 object and the merge manifest, checks coverage against the current
sample manifest, and joins raw Xenium morphology by exact sample-qualified cell
IDs. Only samples declared nonempty in the merge manifest need raw cell tables;
every manifest sample, including empty samples, appears in the coverage table.
The raw/retained cell-area validation uses the production audit's existing
0.00501-square-micrometer rounding tolerance. Missing raw tables, duplicate or
unmatched IDs, inconsistent sample counts, or invalid measurements stop the report.

Outputs go to `outputs/xenium/residual/04_qc_review/`: eight metrics cover
transcripts, detected genes, control percentage, cell area, nucleus area/count,
cell:nucleus area ratio, and transcripts per cell area. The report includes
per-cell measurements, pooled cluster quantiles, per-sample/cluster quantiles,
within-sample target-versus-other-residual-cluster median differences/ratios,
sample-balanced summaries, segmentation diagnostics, sample coverage, provenance,
and a PDF with pooled and per-sample boxplots for every metric and cluster.
Like the RL/VZ QC review, the PDF also includes transcript/gene median-ratio
heatmaps relative to all residual cells in each sample, plus a sample-composition
heatmap. The detailed table records these sample-relative ratios separately from
comparisons against other clusters, so their denominators are explicit.
The cluster-5 question can be assessed from these same outputs without hard-coding
that cluster as low quality. Its identity may change after reclustering.

Sample-balanced comparisons reuse marker review's minimum of 10 measured target
and 10 measured background cells per biological sample, recording excluded small
groups in the detailed table. Each eligible sample has equal weight. This is
comparison eligibility, not a cell-removal rule. Comparisons use other residual
cells, not the full tissue or known-label populations; cell-type differences can
therefore affect QC distributions. Zero/missing nucleus area makes the area ratio
undefined, with missing denominators and nucleus-count anomalies reported
explicitly. An absent finite measurement yields an explicit missing value, not
zero evidence of a QC problem. No hypothesis tests, QC pass/fail classification,
filtering, relabeling, or RDS writes are performed. No QC thresholds change.

Stage 4a's dry run reads only the small merge CSV and checks paths; it does not
load an RDS or scan raw Xenium tables. Expect 5–30 seconds. A full QC report is
estimated at 2–10 minutes for the current 10,333 cells across 11 nonempty samples,
excluding queue wait. This is a low-confidence static estimate dominated by one
integrated RDS read and 11 sequential raw CSV reads. The launcher requests
30 minutes, 32 GiB, one CPU and `express`, following the existing morphology-audit
resource pattern. Collect `sacct` JobID, Elapsed, State, MaxRSS, ReqMem, AllocCPUS,
and ExitCode after the first run to refine the estimate. No integration or marker
rerun is required solely to add this report, and existing outputs are protected
unless `--overwrite` is explicitly supplied after review.

Stage 4b applies an explicitly approved whole-cluster exclusion recorded in
`config/reviewed_residual_qc_removals.csv`. The 2026-10-18 decision removes the
977 cells of resolution-0.8 cluster `5` from the reviewed 10,333-cell source,
leaving 9,356 cells. The source RDS size and modification time, cluster column,
and expected counts bind this decision to that snapshot. A changed source must
receive a new review rather than reusing the old cluster number. The reason is
conservative exclusion after ambiguous markers and consistently lower within-sample
RNA complexity, not an assertion that all excluded cells are technically invalid.

After reviewing Stage 4a and recording approval, run:

```bash
Rscript --max-ppsize=500000 scripts/xenium_residual_04b_apply_qc.R --dry-run
sbatch scripts/run_xenium_residual_04b_apply_qc.slurm
# Once removal completes, run post-QC processing before marker review:
Rscript --max-ppsize=500000 scripts/xenium_residual_05_reintegrate_post_qc.R --dry-run
sbatch --export=ALL,RESIDUAL_QC_ACTION=remove scripts/run_xenium_residual_05_reintegrate_post_qc.slurm
# Wait for successful completion, then regenerate marker review:
Rscript --max-ppsize=500000 scripts/xenium_residual_04_review_subclusters.R --dry-run
sbatch scripts/run_xenium_residual_04_review_subclusters.slurm
```

Stage 4b preserves the original integrated RDS and pre-QC reports. It writes
`04_qc_filtered/rds/Xenium_Residual_QC_Filtered.rds`, an exact excluded-cell
list, a removal-count manifest covering every current sample, and provenance.
The surviving cell IDs, metadata and cluster numbers remain unchanged. Existing
embedding coordinates are retained; Harmony, neighbors, UMAP and clustering are
not recomputed. These coordinates describe the pre-QC embedding, not a new
post-removal integration. Existing QC-removal files require protected `--overwrite`.

Stage 5 now reads that filtered snapshot and refits on retained cells; an
explicit no-removal decision reuses the integrated source. Marker review and
label application require `05_post_qc/rds/Xenium_Residual_postQC.rds` and its
current provenance. Fresh evidence and a source-bound blank approval template
go to `05_post_qc/marker_review/tables/`, preserving the original evidence. Marker backgrounds change after exclusion,
so review these tables and the new cluster partitions before approving fresh labels. Stages 5–6
and combined outputs from before the removal become stale and must be archived
or replaced explicitly before regeneration. Label application checks marker-input
provenance; mapping refuses a master object without the current QC decision.

Stage 6 keeps whole-tissue mapping copies intact and records each excluded cell
as `residual_qc_status=removed` with no residual label. Expected mapped counts are
the Stage-1 counts minus the exact approved exclusions. The combined-label stage
then drops those explicitly listed cells, verifies their sample/cell identities,
and records `Removed_By_Residual` counts and the decision fingerprint. It still
requires a label for every other retained residual target. VZ/RL decisions,
reference routing, consensus definitions and label-transfer settings are unchanged.

Estimated Stage-4b runtime is 1–5 minutes, followed by approximately 1–10 minutes
for post-QC marker evidence, excluding queue wait. These low-confidence estimates
use the current 10,333-cell object and earlier residual marker-review observations.
The requested limits are 30 minutes/32 GiB/one CPU for removal and the existing
4 hours/96 GiB/one CPU for marker review. Dry runs read small CSVs/check paths and
file metadata only, typically 5–30 seconds; they never load an RDS. Collect the
first `sacct` records to refine the estimates. Later mapping still reads/writes
whole-tissue objects and retains its separately documented 5–30-minute estimate.

The lightweight `test_xenium_residual_qc_removal.R` is expected to take under
one minute. It uses tiny isolated synthetic objects to verify the decision
binding, exact exclusions, preserved IDs/metadata/embedding, empty samples,
protected outputs, new marker/label provenance, and downstream exclusion joins.
Rollback consists of preserving the post-QC outputs, restoring the earlier
approved-decision state (including absence of the decision CSV if appropriate),
and restoring the corresponding prior label/mapping/combined outputs together.
Never mix outputs across decision versions.

Stage 1 records a status table for every manifest sample and writes an RDS only
when retained residual cells are present. Stage 2 validates all manifest-defined statuses and
merges only declared nonempty subsets. Stage 3 reuses the approved regional
30-dimensional Harmony and resolution 0.5/0.8/1.0 settings. Stage 4 generates
replicate-aware positive and anti-marker evidence plus a blank label-review
template. Copy the completed reviewed rows to
`config/reviewed_residual_subcluster_labels.csv`; every cluster requires an
explicit decision, confidence, evidence summaries, notes, and `approved=TRUE`.
`final_label=Unknown` is allowed when no defensible winner exists. Then run:

The user-approved 2026-10-18 post-QC resolution-0.8 residual map assigns `0`
to `eCN`, `1` to `BG`, `2` to `Immature UBCs`, `3` to `Microglia`, `4` to
`Pericytes`, `6` to `Mature UBCs`, `7` to `GABA Progenitors`, `8` to
`Cycling Cells`, and `9` to `OPCs`. Cluster `5` was excluded by the separate
approved QC decision. All nine surviving clusters have explicit approved rows
in `config/reviewed_residual_subcluster_labels.csv`; the previous map is
preserved in `config/backups/reviewed_residual_subcluster_labels_before_post_qc_approval_20261018.csv`.
OPCs is explicitly approved despite only one eligible biological sample and
mixed oligodendroglial/myeloid signals; those limitations remain in the evidence
notes and the annotation confidence is low. Approval is a labeling decision and
does not establish adjusted sample-level significance. The current map applies
to the 9,356-cell QC-filtered object, not the historical 34-sample result or any
future reclustering. Existing Stage-5 annotations, Stage-6 mapping, and combined
outputs based on an older label CSV become stale and require explicit regeneration.

```bash
Rscript scripts/xenium_residual_05_apply_reviewed_labels.R --dry-run
sbatch scripts/run_xenium_residual_05_apply_reviewed_labels.slurm

Rscript scripts/xenium_residual_06_map_to_samples.R --dry-run
sbatch scripts/run_xenium_residual_06_map_to_samples.slurm
```

For a pooled residual DotPlot matching the RL/VZ subcluster presentation, run
the separate plotting-only Stage-5 driver after label application:

```bash
Rscript --max-ppsize=500000 scripts/xenium_residual_05_plot_subclusters.R --dry-run
sbatch scripts/run_xenium_residual_05_plot_subclusters.slurm
```

It reads `05_subclusters/rds/Xenium_Residual_Subclusters.rds` and writes
`05_subclusters/plots/Xenium_Residual_SubclusterMarker_DotPlot.pdf` and `.tif`,
plus marker availability/source and input/review provenance CSVs under `tables/`.
This pools all retained residual cells across samples; it does not concatenate
sample DotPlots or represent replicate-adjusted differential expression.
It reuses the same shared VZ/RL panels and residual-only markers as Stage 7,
omits absent identities, checks labels against the approved cluster map, and
uses existing normalized Xenium data. No RDS is written and mapping need not
be rerun. Existing plot outputs require explicit `--overwrite` (launcher:
`RESIDUAL_DOTPLOT_OVERWRITE=true`).

Expected runtime is 1–5 minutes excluding queue wait, low confidence, based on
9,356 cells and the completed 1–4-minute sample plotting tasks. The launcher
requests one CPU, 32 GiB and one hour; the small labeled input previously
completed label application in 40 seconds, but pooled plotting memory has not
yet been measured. Preflight checks paths and dependencies without loading RDS,
typically under one minute. Syntax parsing takes under ten seconds; shell
validation under five seconds. Collect `sacct` elapsed/state/MaxRSS/ReqMem/
AllocCPUS/ExitCode after the first run to refine resources.

Stage 6 maps the reviewed residual identities back to all manifest-defined reviewed
whole-tissue inputs by exact cell ID. It writes one
`<sample>_Residual_Subclusters.rds` object for every manifest sample, including
samples with zero residual cells, plus
`06_mapped/tables/Xenium_Residual_Mapping_Manifest.csv`. The mapping fails if
the Stage-1 input provenance or cell counts changed, if a residual identity is
assigned to a non-Unknown cell, or if any master cell cannot be mapped exactly
once. Existing Stage-6 outputs are protected unless `--overwrite` is supplied.

Observed runtimes were about 1–2 minutes per Stage-1 task (roughly 12–20
minutes total at concurrency three), 2.5 minutes for merge, 2 minutes for
integration, 1 minute for marker review, and 1 minute for label application.
Stage 6 is expected to take 5–30 minutes; this estimate is low confidence and
is based on four concurrent whole-tissue RDS reads/writes plus exact metadata
joins, so collect its first `sacct` record to refine it. Requested limits are
30 minutes, 2 hours, 12 hours, 4 hours, 1 hour, and 4 hours, respectively.

Stage 7 (`xenium_residual_07_plot_samples.R`) reads the Stage-6 whole-tissue
sample objects and writes reports under `residual/07_sample_reports/<sample>/`.
Each task creates global residual spatial TIFFs at sizes 0.6 and 0.85, residual
identity facets at size 0.6, marker DotPlots as TIFF and Cairo PDF, a marker
availability CSV, residual counts/percentages, a provenance/status CSV, and a
single-page spatial/marker PDF containing both global sizes, facets, and DotPlot.
Residual identities use the approved labels and shared master colors. Grey
context represents other cells; the already-approved residual QC exclusions
are omitted from the plotted context. No filtering decision, cell identifier,
metadata, or source RDS is changed. Only identities present in the sample are
faceted or used for marker groups. Shared identities reuse the existing production
panels in `color_palette.R`: RL Other for eCN/GABA Progenitors/Microglia,
RL production UBC for Immature UBCs, the retained CBLN1/RELN display pair
for residual Mature UBCs,
VZ Glia for BG, and VZ OPCs for OPCs. These reuse complete panel definitions,
including the wider Glia and OPC groups, rather than independently selecting
residual genes for those labels. Residual-only identities (currently Pericytes
and Cycling Cells) use the first three unique semicolon-separated genes from
their approved `positive_marker_support`. Marker-availability CSVs record each
gene's panel source. A marker-panel signature in the status file prevents
combining outdated DotPlots/reports after the shared panel definitions change.
Assay-absent genes are recorded and omitted. A
nonempty residual identity with no available panel genes is an error. DotPlots
use existing normalized Xenium data and residual cells only; no normalization
is rerun. Samples with zero mapped residual cells still receive grey-context
plots, explicit empty-report/DotPlot panels, and empty count/marker tables.

After Stage 6 completes, preview and submit both the sample array and its
dependent combined report with:

```bash
bash scripts/submit_xenium_residual_reports.sh --dry-run
bash scripts/submit_xenium_residual_reports.sh
```

The preflight takes under one minute and checks all mapped paths, the mapping
manifest, approved review table, plotting packages, and `pdfinfo`; it does not
load large RDS files. Array bounds come from the active manifest via `--list`.
At 41 samples and concurrency three, estimated runtime is 5–20 minutes per
task and 70–280 minutes for the full array, excluding queue wait (low confidence,
based on the existing VZ/RL reporting workflow and the residual mapping run's
11.30-GiB peak memory). Both sample and combine jobs request a two-hour limit,
64G, and one CPU on `express`. Combining the 41 single-page PDFs is expected
to take 1–5 minutes. Collect `sacct` JobID, Elapsed, State, MaxRSS, ReqMem,
AllocCPUS, and ExitCode after the first run to refine these estimates.

The combined output is
`outputs/xenium/residual/07_sample_reports/Xenium_Residual_AllSamples_Spatial_Report.pdf`
and its adjacent `_manifest.csv`, ordered by numeric PCW and then sample
manifest order. The report audit checks all sample statuses, input sizes/mtimes,
approved review signature, counts, saved output sizes, and PDF page counts
before combining. To audit afterward without loading RDS contents (estimated
under two minutes), run:

```bash
Rscript --max-ppsize=500000 scripts/xenium_residual_07_plot_samples.R --validate
```

Separate submissions are also supported by
`submit_xenium_downstream_array.sh residual-plots`; combine mode uses
`sbatch --array=1 --export=ALL,RESIDUAL_PLOTS_MODE=combine scripts/run_xenium_residual_07_plot_samples.slurm`.
Existing report outputs are protected. Use `--overwrite` or the launcher export
`RESIDUAL_PLOTS_OVERWRITE=true` only after reviewing outputs intended for replacement.
After validation, review the script changes with `git status`/`git diff` and
checkpoint the logical unit when ready; no package environment changes or
`renv.lock` updates are needed for this reporting stage.

### 8. Combined VZ/RL/residual objects

The refined regional labels are combined and analyzed in this order. The two
stage-3 branches can run independently after stage 2.

- `xenium_vz_rl_01_combine_labels.R` – combines refined labels in all manifest-defined configured individual samples;
- `xenium_vz_rl_02_correct_broad_labels.R` – writes separate corrected-broad sample objects, plots, and audit tables;
- `xenium_vz_rl_03_nonspatial_merge.R` – cleans and merges per-sample objects;
- `xenium_vz_rl_04_05_nonspatial_process_and_plot.R` – processes and visualizes the combined object;
- `xenium_vz_rl_03_spatial_merge.R` – memory-conscious spatial merge;
- `xenium_vz_rl_04_spatial_integrate.R` – Seurat v5 sketch/integration workflow on the merged data;
- `xenium_vz_rl_optional_cluster_counts.R` – optional lineage percent-abundance plots, still reading the original combined objects rather than corrected-broad objects.

`xenium_vz_rl_reports.R` builds separate all-sample spatial PDF reports from
the completed combined-subcluster and corrected-broad TIFFs, without loading
or rewriting RDS objects. It follows the same PCW page ordering as the VZ/RL
Stage-8 reports, with one full spatial TIFF (including its legend) per sample
page. By default it uses both label sets and the existing point-size-0.6 plots;
`--labels=combined|corrected-broad|both` selects the label sets, and
`--point-size=0.85` selects the alternative existing TIFFs. No TIFFs are regenerated.
The driver recognizes both the renamed layout and the previous
`01_combined_labels/plots/combined/` and `01_corrected_broad/plots/` locations.
It uses one location per label set and refuses duplicate or split sets across
layouts. Missing files still fail the preflight; no paths or plots are moved.

```bash
Rscript --max-ppsize=500000 scripts/xenium_vz_rl_reports.R --dry-run
sbatch scripts/run_xenium_vz_rl_reports.slurm
```

The path/metadata preflight should take under one minute and does not decode
TIFF pixels or load RDS files. It checks the manifest-defined TIFF set,
PCW metadata, `tiff`, and Poppler `pdfinfo`. The report job
is estimated at 5–20 minutes for the current 41-sample cohort (low confidence,
based on sequentially decoding 82 existing 600-dpi TIFFs); its requested limit
is two hours with 64G and one CPU on `express`, following the existing Stage-8
report launcher resource class. Queue wait is excluded; collect its first
`sacct` record to refine the estimate.
Outputs are `reports/Xenium_VZRL_AllSamples_CombinedSubclusters_Spatial_Report_pt0.6.pdf`
under `01_combined_labels/` and
`reports/Xenium_VZRL_AllSamples_CorrectedBroad_Spatial_Report_pt0.6.pdf`
under the detected corrected-broad stage (`02_corrected_broad/` or the previous
`01_corrected_broad/`), each with an adjacent `_manifest.csv` recording
page/sample/PCW, label column, point size, TIFF path/size/mtime, and PDF path.
The driver verifies the final page count and refuses existing PDFs or manifests
unless `--overwrite` is supplied. For the launcher use
`VZRL_REPORT_LABELS`, `VZRL_REPORT_POINT_SIZE`, and `VZRL_REPORT_OVERWRITE=true`
as the corresponding exported controls. Run this after the sample TIFF jobs
finish; it can run alongside later merge jobs that only read those inputs.

Combined outputs are organized under `outputs/xenium/vz_rl/`:

- `01_combined_labels/` contains all manifest-defined labelled sample objects, original consensus and fine-subcluster spatial plots directly in `plots/`, and count tables. Runs made before this layout change placed those four TIFFs per sample in `plots/combined/`; move them into `plots/` only after checking for filename collisions. The RDS and tables do not need to be regenerated;
- `02_corrected_broad/` contains separate fine-derived broad-label sample objects, spatial plots, count tables, and a manifest; original combined files remain intact;
- `nonspatial/03_merged/` contains the cleaned sample objects, merged RDS, and merge manifest;
- `nonspatial/04_processed/` contains the processed object and provenance manifest; `nonspatial/05_plots/` contains its summary plots;
- `spatial/03_merged/` and `spatial/04_integrated/` contain the spatial-preserving branch.

For an existing completed run, first sync the renamed scripts to Cheaha, then
run `python3 scripts/migrate_vz_rl_layout.py` to preview the moves and
`python3 scripts/migrate_vz_rl_layout.py --apply` after all VZ/RL jobs have
left the queue. This migration moves directories on the same filesystem,
flattens `01_combined_labels/plots/combined/`, updates absolute path columns
in the small provenance manifests, and saves their originals with a
`.pre_vzrl_layout_20261019.csv` suffix. It archives older script names in
`scripts/legacy_vz_rl_names_20261019/` after checking that their renamed
replacements are present. It does not read or regenerate RDS
objects. Check the corrected-broad validator and both merge `--dry-run` results
afterward; the migration refuses a partial stage, a destination collision, or
an active VZ/RL job.

The earlier `xenium_vz_rl_diagnostic_compare_unknown_annotations.R` dual-branch
comparison and `config/reviewed_unknown_relabels.csv` are retained only as
historical/sensitivity tools. They are not production prerequisites after
exclusive reference routing. Do not run the comparison unless a reviewer asks
for the prior include-in-both-branches sensitivity analysis.

The spatial branch has a separate validated merge because it preserves the
FOV/image data removed from the smaller non-spatial merge. Inspect both stages
before submission:

```bash
Rscript --max-ppsize=500000 scripts/xenium_vz_rl_03_spatial_merge.R --dry-run
Rscript --max-ppsize=500000 scripts/xenium_vz_rl_04_spatial_integrate.R --dry-run
```

The merge requires all manifest-defined manifest inputs and writes
`spatial/03_merged/rds/Xenium_VZRL_spatial_merged.rds` plus a cell/image and
regional-assignment count manifest under `spatial/03_merged/tables/`.
It merges ordered groups of six samples, then merges those groups. Set
`SPATIAL_MERGE_BATCH_SIZE=2` through `10` when submitting the launcher to
compare future runs; the default is six. The manifest records the group for
each sample. The merge checks the exact cell IDs and order, total cell and
image counts, and per-sample assignment counts before saving. Batch merging
has not yet been benchmarked on the full 41-sample spatial dataset; the
24-hour Slurm limit and 600-GiB request remain conservative pending a run.
The merge uses one CPU because its Seurat merge calls run sequentially.
Sketch-based Harmony processing
writes `spatial/04_integrated/rds/Xenium_VZRL_spatial_integrated.rds`. UMAPs
are TIFF-only. Giotto, h5ad export, and SpaTrack scripts remain outside the
current scope and still require path review before they are reactivated.
Submit these stages with `run_xenium_vz_rl_03_spatial_merge.slurm` and then
`run_xenium_vz_rl_04_spatial_integrate.slurm`.

Before creating the combined per-sample RDS files, spatial TIFFs, and count
tables, inspect input completeness and existing outputs with:

```bash
Rscript --max-ppsize=500000 scripts/xenium_vz_rl_01_combine_labels.R --dry-run
```

This stage reads sample IDs from `config/samples.csv`, requires the reviewed
whole-tissue input manifest, the reference-route manifest, and all manifest-defined
independently mapped VZ, RL, and residual objects. It refuses to replace any
existing output unless given `--overwrite`. Exact cell IDs and unchanged
consensus labels validate every join. Known regional cells receive their VZ or
RL subcluster. A routed Unknown population can receive a fine label from only
its supported branch. Retained `route_other` and `route_unresolved`
populations receive their explicitly reviewed residual label; the approved
residual `Unknown` cluster remains `Unknown` rather than being forced. The
original `consensus_label`, `unknown_reference_route`, residual reference
route, confidence, and annotation status preserve provenance.
`regional_assignment` records VZ, RL, `Residual`, or
`unassigned_consensus`, and a per-sample routing table records the resulting
cell counts. Any overlap among VZ, RL, and residual assignments is an error.
Submit it with `run_xenium_vz_rl_01_combine_labels.slurm`.
The first reported 41-sample combined-label run completed in 44:57 (job
40971846, 2026-10-19; output set and MaxRSS still require independent checks).
For a repeat of the same cohort, expect approximately 40–90 minutes, based on
that single observed runtime; the estimate remains low confidence until more
runs and output checks are available.
The launcher checks the VZ/RL/residual mapping manifests and combined-driver
dry run before processing large RDS objects. It requests 24 hours and 128G on
`medium`; queue wait is excluded.
Collect `MaxRSS` with `sacct` before changing the requested memory.

Before either merge, apply the user-approved 2026-10-19 fine-to-broad map in
`R/corrected_broad_labels.R` to all retained cells. This adds
`corrected_broad_label` while preserving the original `consensus_label` and
every cell ID. PC subclusters become `Purkinje`; Pre-OPCs/OPCs become `OPC`;
`iCN` becomes `GABA`; RL VZ/SVZ become `RL`; radial/Bergmann/cycling glial
subclusters become `Glia`; `eCN`, `Pericytes`, `Cycling Cells`, and `Unresolved`
are distinct broad identities. Granule, UBC, and other reviewed subclusters
follow the corresponding broad lineage. Run the new read-only preflight:

```bash
Rscript --max-ppsize=500000 scripts/xenium_vz_rl_02_correct_broad_labels.R --dry-run
```

It checks 41 combined inputs and 206 expected new outputs without loading
large RDS objects (usually under one minute). After `DRY-RUN PASS` and
`outputs-existing=0/206`, submit
`sbatch scripts/run_xenium_vz_rl_02_correct_broad_labels.slurm`.
The first reported 41-sample run completed in 12:56 (job 40973926,
2026-10-19; 16.44 GiB batch MaxRSS). The post-job validator passed for 41/41
samples, 2,355,487 cells, and 82 plots; its CSV/path checks do not load RDS
contents. For the
same cohort, expect about 10–30 minutes, low confidence from one observed run.
The launcher requests 12 hours, 128G, and one CPU on `short`.
Check the 41-row corrected-broad manifest and per-sample
`*_consensus_to_corrected_broad.csv` tables before merging. Existing
`01_combined_labels/` objects and plots are not replaced.
After job completion, run
`Rscript --max-ppsize=500000 scripts/validate_xenium_vz_rl_02_correct_broad.R`
to check all 41 paths, 82 nonempty new plots, counts and transition tables,
and mapping provenance without loading RDS objects (estimated under two
minutes; filesystem latency may vary). The later merge preflight checks that
the corrected-broad manifest and objects still match the original combined
inputs and the approved mapping script.

Applying or changing the reviewed Unknown QC decisions, reference-route
manifest, or reviewed residual labels makes every existing combined Stage-1
output and all downstream VZ/RL/residual merge, processed, spatial, count, and
trajectory-derived outputs stale. Regenerate them only after reviewing the
Unknown-QC, route, and residual-label tables. A rollback is the prior reviewed
QC and residual-label CSVs plus the corresponding protected Stage-1 outputs;
the whole-tissue consensus and mapped branch objects are never modified.

The subsequent clean/merge stage requires all manifest-defined corrected-broad
objects and their current 41-row provenance manifest. It removes
spatial images, scale data, and control assays before merging:

```bash
Rscript --max-ppsize=500000 scripts/xenium_vz_rl_03_nonspatial_merge.R --dry-run
```

It writes the stable merged object
`nonspatial/03_merged/rds/Xenium_VZRL_merged.rds` and a per-sample cell-count, PCW, and
regional-assignment manifest under `nonspatial/03_merged/tables/`. Combined integration and plotting are
separate operations:

```bash
Rscript --max-ppsize=500000 scripts/xenium_vz_rl_04_05_nonspatial_process_and_plot.R --dry-run
Rscript --max-ppsize=500000 scripts/xenium_vz_rl_04_05_nonspatial_process_and_plot.R --process-only
Rscript --max-ppsize=500000 scripts/xenium_vz_rl_04_05_nonspatial_process_and_plot.R --plots-only
```

The stable processed object is
`nonspatial/04_processed/rds/Xenium_VZRL_processed.rds`; plots are written to
`nonspatial/05_plots/`.
Plot-only runs verify that the input merge has not changed, and they never
repeat normalization, PCA, Harmony, or UMAP. UMAP outputs remain TIFF-only;
DotPlots, the marker heatmap, and violin plots are saved as TIFF and Cairo PDF.
The corresponding launchers are `run_xenium_vz_rl_03_nonspatial_merge.slurm`,
`run_xenium_vz_rl_04_nonspatial_process.slurm`, and `run_xenium_vz_rl_05_nonspatial_plots.slurm`.
For an isolated rerun from an archived merge, export
`COMBINED_MERGE_INPUT_DIR` to the directory containing its `rds/` and `tables/`,
and `COMBINED_NONSPATIAL_OUTPUT_DIR` to a fresh destination. Paths may be
project-relative. Use both variables unchanged for processing and plotting;
the destination receives `04_processed/` and `05_plots/`. The merge backup is
read-only, existing outputs remain protected, and unset variables retain the
production paths. Run `--process-only --dry-run` before processing and
`--plots-only --dry-run` after processing; these inspect paths without loading
RDS objects (normally under one minute). Allow roughly 4–24 hours for processing
and 0.5–4 hours for plotting, excluding queue wait; these are low-confidence
static estimates for the full cohort, not measured job runtimes. The existing
launchers request 48 and 24 hours respectively.
The plot-only stage writes 24 figure files plus
`Xenium_Merged_plot_manifest.csv` and `Xenium_Merged_Heatmap_marker_summary.csv`. All eight plot families have an all-cell
version and a `NoUnknown` version. UMAPs use original broad consensus labels,
corrected broad labels, or combined subcluster labels; the batch comparison
uses the sample manifest order and one shared legend below its two panels.
Original consensus assignments remain available in the original-label UMAP;
the redundant original-consensus marker DotPlot is retired and replaced by
the corrected-broad marker DotPlot. Existing legacy files are never deleted.

The marker heatmap and signaling violins now group cells by
`corrected_broad_label`. Heatmap marker-test thresholds, per-identity test cap
(500), and display cap (100 cells) remain unchanged. This changes which cells
contribute to each plotted broad group and its top marker genes; regenerate
these figures, but no merge, integration, labels, or saved RDS needs rebuilding.
Both heatmaps share the full-cohort marker tests/scaling and the same sampled
known cells. Unresolved may have no qualifying positive markers: retain its
heatmap columns, report zero markers in the marker-summary CSV and plot notes,
and display a subtitle explaining the absent marker block. Missing markers for
other identities and entirely empty marker results remain errors. Thresholds
are unchanged. This explicitly supported exception fixes job 41006727, which
failed at this check after 24 minutes 53 seconds with peak batch RSS 182.34 GiB. Cluster columns, dotplot rows/section boxes, and figure legends
follow the broad order, with eCN immediately after Granule and Pericytes after
Endothelial. Dotplot section boxes use their registered cluster colors.
The redundant heatmap Identity key is removed; the expression colorbar remains.
UMAP labels use repulsion and expanded plotting space, with ordered keys below.
Their final placement still needs visual inspection on the full cohort.

NoUnknown UMAPs filter using their displayed label column; batch plots,
corrected-broad DotPlots and heatmaps use corrected broad labels. The filter
recognizes Unknown and numbered Unknown labels, retains Unresolved, and never
changes saved cells or embeddings. Purkinje-only and selected signaling panels
already exclude Unknown identities, so their paired files are intentionally
identical. The plot manifest records filenames, grouping, displayed/total
cell counts, seed, input RDS path, and these exceptions. All-cell filenames
retain their established names; new variants insert `_NoUnknown` before
`_UMAP` or append it to the other stems.

Use `COMBINED_PLOT_OUTPUT_DIR` to render a fresh review folder while reusing the
existing processed object. This overrides only the figure destination; the
existing merge/processed directory overrides still work independently. See
`scripts/NONSPATIAL_PLOT_REFRESH.md` for the safe HPC dry-run and submission.
Allow approximately 30–90 minutes for the paired plotting job, excluding queue
wait (low confidence: failed job 41006727 reached the end of marker testing in
about 25 minutes; remaining heatmap/violin rendering time is unmeasured). The unchanged
launcher requests 24 hours, one CPU and 600 GB. Refine this estimate using the
completed job's accounting. The path-only dry-run should take under one minute
and does not load RDS objects. The spatial merge also consumes the
corrected-broad objects, and spatial integration adds a separate
corrected-broad UMAP.

### 9. Spatial and trajectory analyses

#### Giotto

`Xenium_Giotto_Analysis.R` and `Xenium_Giotto_Analysis_array.R` convert annotated Seurat objects to Giotto objects and calculate:

- Delaunay spatial networks;
- spatially variable genes with `binSpect()`;
- cell-proximity enrichment; and
- spatial ligand–receptor interactions using Xenium-detectable pairs from the human CellChat database.

`Xenium_Giotto_BroadCluster_Analysis_array.R` runs the corresponding broad-cluster analysis. Plotting helpers include `Xenium_Giotto_SVGPlots.R`, `Xenium_Giotto_SVGCompPlots.R`, and `Xenium_Giotto_TimePlots.R`.

Launchers are `run_GiottoAnalysis.slurm`, `run_GiottoAnalysis_array.slurm`, and
`run_GiottoBroadAnalysis_array.slurm`.

#### SpaTrack and CellRank export

- `Xenium_VZRL_Subclusters_VZ_h5adExport.R` exports VZ-lineage expression, metadata, embeddings, and spatial coordinates to `h5ad`.
- `Xenium_VZRL_Subclusters_VZ_SpaTrack.R` prepares/runs the per-sample SpaTrack workflow through `reticulate`.
- `Xenium_VZRL_Subclusters_VZ_SpaTrack_Plotting.R` imports pseudotime and creates spatial plots.
- `run_VZSpaTrack_array.slurm` submits per-sample trajectory jobs.

These scripts require configured Python/Conda environments and currently contain environment-specific paths.

#### Other downstream analyses

- `Xenium_GeneExpGradient_Analysis.R` – spatial expression gradients across selected lineages or regions.
- `xenium_annotate_03c_plot_proportions.R` – proportions by sample metadata such as post-conception week.
- `xenium_vz_compare_aldinger.R` and `xenium_vz_compare_sepp.R` – average-expression correlations between Xenium VZ clusters and reference VZ clusters.

### 10. Cross-study reference comparisons

- `cross_study_cluster_correlations.R` computes pseudobulk cluster-expression correlations across Aldinger, Sepp, and Science.
- `cross_study_cluster_marker_comparison.R` compares cluster marker programs and produces heatmaps.
- `cross_study_module_scores.R` transfers dataset-specific marker modules among references.
- `run_cross_study_cluster_marker_comparison.slurm` submits the marker comparison.

For all-pairs expression matching of corrected Xenium broad populations and
subclusters to original reference data and published annotations, see
[`XENIUM_REFERENCE_SIMILARITY.md`](XENIUM_REFERENCE_SIMILARITY.md) and
`xenium_compare_reference_similarity.R`. It reports Spearman correlations,
population-specificity correlations, competing matches, and conditional
sample-bootstrap stability without changing annotations. Original panel JSON Ensembl IDs
align Sepp and Science genes with an auditable one-to-one match. The lightweight dry run
should take under one minute; a full HPC comparison is estimated at 1
to 6 hours or longer (low confidence, based on static inspection; queue wait excluded).
Original Aldinger/Sepp inputs are supported. Science reads original HDF5 counts
in bounded chunks, using the user-approved cerebellum/9–15 PCW cohort and
original published broad CellClass names (numeric source IDs remain audited).
Sepp cells lacking original named annotations are excluded and audited.
Heatmaps include numeric legends and separate raw/specificity color scales.
Its metadata-only `--science-plan` preflight is
estimated at 1–3 minutes from the observed 67-second metadata audit.
Copy `R/science_original_reference.R` and `R/reference_similarity_plots.R` with
the driver. `xenium_redraw_reference_similarity.R` redraws revised correlation
CSVs without expression loading (typically <1–2 minutes). Both plot paths use
master_subcluster_order and celltype_order from color_palette.R, top to bottom.
Cached reference ranks
and one-pass query aggregation speed up full comparisons. The approved bootstrap
default is 100 iterations (seed 42). The guide explains why a
nearest match alone does not establish that a subtype is present.

For a separate cell-level cluster-replicability comparison using the original
references, see [XENIUM_METANEIGHBOR.md](XENIUM_METANEIGHBOR.md),
`xenium_compare_metaneighbor.R` and `run_xenium_compare_metaneighbor.slurm`.
The three-task reference array runs both Xenium broad and subcluster labels
against each original reference's published labels, permitting coarser names.
It exports directed MetaNeighborUS AUROCs, mean scores, competing matches and
reciprocal best hits without changing annotations. The official low-memory
MetaNeighbor package is required in project renv. Path/package preflights take
under one minute each; a synthetic engine smoke test takes roughly 10–60 seconds.
The first full run is provisionally estimated at 30 minutes–4 hours per reference,
or 1.5–12 hours total at concurrency one, excluding queue wait (low confidence,
static estimate; no completed MetaNeighbor run measured). The launcher requests
12 hours and 256 GiB per task; see the guide for provenance and limitations.

## Software requirements

Every current shell launcher and submission helper starts R with
`Rscript --max-ppsize=500000`, including list, dry-run, and inline-expression
calls. This raises the pointer-protection stack limit used during large object
imports and processing. It does not change analysis parameters or Slurm memory
requests, and it cannot guarantee that every protection-stack overflow will be
avoided. Keep the option in new launchers as required by `AGENTS.md`.

Check this rule from `scripts/` with
`python3 update_r_launchers_stack_limit.py --check`. The checker also runs
`bash -n` on current `.sh` and `.slurm` files. It should take under 1 minute
and runs no R analysis. To update an existing HPC checkout, copy the updater
into its `scripts/` directory, use `--dry-run`, then `--apply`; it preserves
previous file versions in a timestamped backup directory. Historical directories
are excluded. Existing stage runtime estimates and requested limits still apply.

The configured Slurm environment uses R 4.4.1 from the
`R/4.4.1-foss-2023b` module. Frequently used R packages include:

- `Seurat`
- `harmony`
- `future` and `future.apply`
- `dplyr`, `tidyr`, `purrr`, and `tibble`
- `ggplot2`, `patchwork`, `randomcoloR`, and `Cairo`
- `ComplexHeatmap`, `pheatmap`, `circlize`, and `viridis`
- `Giotto`
- `reticulate`
- `readxl` and `jsonlite`
- `AnnotationDbi` and `org.Hs.eg.db`

Some scripts call `source("renv/activate.R")`, while project-root jobs also
activate `renv` through `.Rprofile`. Active launchers do not install, restore,
refresh, or snapshot dependencies. Run `renv::restore()` manually once before
submission when the lockfile changes; compute jobs must treat the library as
read-only.

Giotto and trajectory scripts additionally require the corresponding Python/Conda environments. Giotto downloads the human CellChat database at runtime, so compute nodes need network access or the database must be cached locally.

## Running the workflow

Run R scripts from the **project root**, not from `scripts/`, so `here()` resolves the expected data and output paths.

```bash
cd /path/to/XeniumFCProject
Rscript scripts/xenium_preprocess_single_slides.R --list
```

For a script that consumes a one-based task ID:

```bash
Rscript scripts/xenium_vz_01_subset_samples.R 1
```

On the configured cluster, submit the matching launcher:

```bash
sbatch scripts/run_xenium_vz_01_subset.slurm
```

To inspect both preprocessing task maps and dry-run one task from each:

```bash
Rscript scripts/xenium_preprocess_single_slides.R --list
Rscript scripts/xenium_preprocess_single_slides.R --dry-run 1
Rscript scripts/xenium_preprocess_split_slides.R --list
Rscript scripts/xenium_preprocess_split_slides.R --dry-run 1
```

To submit all tasks after reviewing the dry runs:

```bash
sbatch scripts/run_xenium_preprocess_single_slides.slurm
sbatch scripts/run_xenium_preprocess_split_slides.slurm
```

The single-slide launcher submits 19 tasks and the split-slide launcher submits
15; each runs at most three tasks concurrently. Every task resolves its input
directory and cell-stat CSV from `config/samples.csv`. The drivers stop if any
expected output already exists. An intentional rerun requires `--overwrite`,
`SINGLE_PREPROCESS_OVERWRITE=true`, or `SPLIT_PREPROCESS_OVERWRITE=true`, as
appropriate. Both launchers request Intel nodes because the project's native R
packages are compiled with Broadwell flags; loading `yaml` on the AMD nodes
tested here causes an illegal-instruction failure.

Before submission, verify:

1. the script's `--list` task mapping;
2. the Slurm `--array` range;
3. input and output RDS filenames;
4. requested memory, time, CPUs, and partition;
5. the project and Python paths; and
6. that the launcher invokes the intended R script.

### Cheaha partition selection

`config/cheaha_slurm_partitions.csv` records the Cheaha partition table
provided on 2026-09-01. It is a dated planning reference, not a guarantee of
future availability. Recheck current Cheaha policy before changing or
submitting important jobs.

For ordinary CPU batch jobs, choose the highest-priority partition that safely
contains the evidence-based requested time limit:

- `express`: at most 2 hours; 52 nodes; unlimited nodes per researcher;
  priority tier 20;
- `short`: over 2 and at most 12 hours; 52 nodes; 44 nodes per researcher;
  priority tier 16;
- `medium`: over 12 and at most 50 hours; 52 nodes; 44 nodes per researcher;
  priority tier 12;
- `long`: over 50 and at most 150 hours; 52 nodes; 5 nodes per researcher;
  priority tier 8.

Use `interactive` only for interactive work. Use `intel-dcb`, `amd-hdr100`, GPU,
or large-memory partitions only when the workload explicitly needs that
hardware or cannot fit ordinary CPU nodes. The project's Intel constraint is
still required for launchers using native R packages compiled with Broadwell
flags. A higher priority tier does not override fair-share, matching-node
memory, hardware constraints, or per-researcher node limits.

Measured job summaries belong in `slurm_resource_history.psv`; raw Slurm logs
remain outside Git. Before creating or editing a launcher, review that file and
recent `sacct` evidence, then keep `--time` safely above expected runtime while
placing the job on the narrowest valid partition. Do not reduce a justified
time limit solely to qualify for a higher-priority partition.

## Output conventions

Most outputs are written beneath `outputs/` in stage-specific directories. Common artifact types are:

- uncompressed Seurat `.rds` files for large intermediate objects;
- 600-DPI TIFF figures written with `Cairo::CairoTIFF` (using Cairo's default
  TIFF compression behavior);
- Cairo PDF companions for histograms, dot plots, bar graphs, and other
  non-spatial graphs;
- TIFF-only output for UMAPs and spatial plots, except for spatial QC pages
  retained in the combined QC PDF report;
- CSV marker, QC, spatial-gene, and comparison tables;
- Giotto objects and interaction statistics in `.rds` format; and
- AnnData `.h5ad` files for Python trajectory tools.

Some reference-preparation and optional comparison artifacts still use manually
assigned dates. The primary Xenium preprocessing, annotation, and VZ/RL branch
contracts use stable stage names.

### Shared DotPlot presentation

The 2026-10-17 production RL/VZ Stage-6 DotPlots render facet backgrounds with
their broad identity colors from `cluster_colors`: RL pink, UBC orange, Granule
green, VZ magenta, Purkinje blue, GABA teal, Glia light blue, and OPC yellow.
Groups without a corresponding broad identity remain neutral. The shared
`workflow_dotplot_grob()` colors the final plot strips for both PDF and TIFF;
it does not change expression values, dot sizes, color scales, or identities.
No additional plotting package is required.

The RL production Granule panel uses the six user-selected genes, in order:
SFRP1, HEY1, NTF3, CNTN2, GAP43, and GABRA6. RL and VZ Other DotPlots use a
6-inch minimum width, 0.22 inches per gene plus label space, and vertically
stacked top legends; merged and lineage plots keep their existing dimensions. The BG Progenitors panel uses SLC1A3/SOX2/BCAN;
Maturing PCs uses CALB1/FOXP4/DAB1. Other-panel facets are restricted to observed
identities, and unused factor levels and colors are dropped from the RL heatmap.
Raw Stage-6 review marker panels/signatures and the approved identity map are
unchanged. The heatmap's display-only cell sample explicitly uses the configured
seed. Existing stable output names are preserved.

Use the two Stage-6 `--plots-only --resolution=0.8 --dry-run` modes before
submitting `RL_SUBCLUSTER_MODE=refresh-plots` or
`VZ_SUBCLUSTER_MODE=refresh-plots` with their protected overwrite switches.
Preflights inspect paths only and take under 1 minute. Plot refreshes are
estimated at 15-40 minutes for RL and 20-60 minutes for VZ, excluding queue wait
(low confidence based on the 18m56s RL annotation job and larger VZ object).
Both launchers request 12 hours. Refreshes regenerate production figures and
descriptive marker/QC tables without saving an RDS or rerunning integration,
clustering, or annotation. The local R installation lacks Cairo shared libraries;
synthetic layout validation used base PDF, while production keeps Cairo exports.

Current reference, label-transfer, consensus, VZ/RL, and focused-marker DotPlots
use `R/dotplot_helpers.R` (normally sourced through `color_palette.R`). Dot size
always maps 0–100% to the same radius range (0–6); expression colors always use
light grey to red with limits −2.5 to 2.5. Both legends appear above the plot,
with Percent Expressed first. Gene and identity labels use 12-point text;
cluster/subcluster labels use plain text. VZ/RL standalone figure dimensions
follow the displayed gene and identity counts to avoid oversized spacing.

The shared color range is a display convention, not a pooled normalization:
Seurat still scales each gene within the plotted dataset. Existing single-identity
plots retain `scale = FALSE` and explicitly label the legend as log-normalized;
their expression values must not be interpreted as z-scores. Values outside the
shared color limits saturate at the endpoints. No normalization, QC, clustering,
annotation, identity selection, or saved analysis object changes are implied.

Existing TIFF/PDF figures and reports must be regenerated to show this styling;
RDS outputs remain valid. Review protected plot-only modes and output paths
before authorizing replacement. Local syntax checks should take under 10 seconds
per script and path-only dry-runs under 1 minute (static-inspection estimates);
full figure regeneration depends on object size and should be timed on the HPC.

## Reproducibility notes and known caveats

- Several scripts contain absolute paths under `/home/acflint`, `/data/user/acflint`, or `~/R/Projects/XeniumFCProject`. Replace these when running elsewhere.
- General preprocessing, consensus labelling, and VZ/RL stages use
  `config/samples.csv`; a few targeted plotting analyses still define selected
  scientific cohorts explicitly.
- VZ and RL post-QC merged processing use the matching
  `run_xenium_vz_05_reintegrate_post_qc.slurm` and `run_xenium_rl_05_reintegrate_post_qc.slurm`
  launchers.
- Large objects can require 64–600 GB RAM depending on the stage. Prefer
  fork-based shared-memory workers for the largest read-mostly graphs; do not
  serialize those graphs into multiple PSOCK workers.
- Random seeds are set in major Seurat stages. VZ and RL subcluster identities remain manually reviewed mappings whose exact cluster coverage is validated before output.
- The ignored `OLD/` directory contains superseded or experimental versions
  and is not part of the active workflow. This includes
  `AnchorBasedTransfer_RCTD_res1.5.R`, `XeniumABT_seq.R`, `SeppPreProcess.R`,
  `SaveRawClusterPlots.R`, `XeniumPreProcess_seq.R`, the two former
  `xenium_preprocess_selected_samples*.R` drivers, and the broken
  `run_NormCluster1.5.sh` launcher.

## Reviewed Science resolution-2 mapping, removal, and retained-cell reprocessing

The mapping approved on 2026-10-22 is recorded in
`reviews/Science_res2_reviewed_labels_20261022.csv`. It covers all 45 clusters in
the reviewed 70,731-cell, 16-sample Science run. Labels include eCN for 26/44,
GABA for 23/32, Glia for 24, VZ for 28/40, and Vascular for 39/42. There is no
Cycling label. Brainstem (34), Other (36), and RBC (37) are explicitly removed:
**1,366 cells removed; 69,365 retained**. The review includes expected cluster
counts and ordered top-20 gene evidence; mismatches are errors to prevent using
the mapping on another run's cluster IDs.

Submit the reviewed mapping/removal/reprocessing job from `scripts/`:

```bash
sbatch run_science_map_reviewed_res2.slurm
```

The launcher calls `fetal cerebellum_science_ACF.R --resolution=2 --map-reviewed`,
which dispatches to `science_03_map_reviewed_res2.R`. It reads the existing saved
resolution-2 RDS and marker evidence from `outputs/references/science/reprocessed/res2/`;
it does not repeat the original raw clustering before mapping. The path-only
`--dry-run` checks paths, dependencies and review-table structure, without loading
RDS objects (expected under one minute). Cell alignment, cluster counts, marker
evidence and sample retention are checked when the real job loads the source.

The retained cells repeat the established originalexp normalization, cell-cycle
scoring and regression of S.Score/G2M.Score, PCA with 50 components, Harmony on
sample_id and dimensions 1:40, neighbors and resolution-2 clustering, and UMAP on
Harmony dimensions 1:30. Seed 42 is set. Removing cells changes these fitted
representations and can change cluster boundaries and numeric IDs. Reviewed
`science_broad_label` and original `science_cluster_res2` stay attached to each
retained cell; new numeric IDs are saved in `science_cluster_res2_postreview`.
No old numeric-cluster label map is applied to the new clusters. All retained
cell/sample IDs and unrelated metadata are checked for preservation.

Expected runtime is **1–8 hours excluding queue wait**, very low confidence from
static inspection and 69,365 retained cells, including marker testing. The job
retains **12 hours, one CPU, 256 GB, short** from the established full reprocessing
launcher. No full Science reprocessing history is available in
`slurm_resource_history.psv`; collect JobID, Elapsed, State, MaxRSS, ReqMem,
AllocCPUS and ExitCode with `sacct` after this run. Logs use
`logs/Science_res2_reviewed_<job_id>.out` and `.err`.

New protected outputs go under
`outputs/references/science/reprocessed/res2/reviewed_retained/`, prefix
`Science_res2_reviewed`:

- `_UMAP_before_removal.pdf`: all reviewed broad labels, original coordinates.
- `_UMAP_after_removal_original_coordinates.pdf`: retained cells, same axis limits
  and original coordinates, to show exactly what was removed.
- `_UMAP_after_reprocessing.pdf`: reviewed labels on the newly fitted UMAP.
- `_raw_clusters_after_reprocessing_UMAP.pdf`: newly generated numeric clusters.
- `_broad_marker_DotPlot.pdf`: broad-marker evidence after reprocessing.
- `_retained_reprocessed.rds`: retained, reprocessed object with active broad identities.
- `_cell_mapping.csv`: all original cells, source clusters, labels, removal actions
  and post-review cluster IDs (explicitly missing for removed cells).
- `_cluster_retention.csv` and `_sample_label_counts.csv`: reconciled retention counts.
- `_postreview_top20_markers.csv` and `_postreview_marker_counts.csv`: markers for new
  numeric clusters, keeping only.pos=TRUE, min.pct=0.25 and logfc.threshold=0.25.
- `_provenance.csv` and `_sessionInfo.txt`: inputs, review fingerprint, settings and environment.

Existing source RDS, marker evidence and earlier plots remain unchanged. Replacing
reviewed outputs requires explicit `--overwrite` or
`SCIENCE_REVIEWED_OVERWRITE=true` after review. A full rerun does not automatically
switch production reference paths or regenerate downstream label transfers;
those require separate adoption and validation of this new reference. Review
`git status`/`git diff` in the tracked checkout and checkpoint this logical unit
after validation when ready; package versions and renv.lock are unchanged.

## Reprocessed Science resolution-2 clusters, top-20 markers, and eCN evidence

`fetal cerebellum_science_ACF.R` now runs the existing 9/12/13/15-week Science
reprocessing workflow at **resolution 2**, preserving raw identities in
`science_cluster_res2`. Other normalization, cell-cycle regression, Harmony,
and dimensionality settings are unchanged. Seed 42 is set for reproducibility.
Changing resolution changes cluster membership and numeric identities;
the legacy resolution-1.2 broad-label map is retained in the script but is
not applied at resolution 2. Review the new markers before assigning labels.

Resolution 2 is the default. Both resolutions 2 and 3 remain supported with
`--resolution=2` or `--resolution=3`; each writes to its own resolution directory.
The existing `run_science_reprocess_res3.slurm` explicitly pins resolution 3.
Resolution-2 results require fresh markers, plots, and eCN summaries; existing
resolution-3 results remain available for comparison.

Submit the full reprocessing and marker-export job from `scripts/`:

```bash
sbatch run_science_reprocess_res2.slurm
```

Expected runtime is **1–8 hours**, excluding queue wait, with very low confidence
from static inspection and 70,731 cells in the prior report. Marker testing can
dominate runtime. The job requests **12 hours, one CPU, 256 GB on short**; this
keeps the existing Science memory envelope with time headroom. Cheaha's current
[hardware table](https://docs.rc.uab.edu/cheaha/hardware/) confirms short allows
12 hours (checked 2026-10-21). There is no array. Logs use
`logs/Science_res2_top20_<job_id>.out` and `.err`. Collect elapsed time and MaxRSS
after completion to refine these estimates.

Outputs are protected under `outputs/references/science/reprocessed/res2/`:

- `Science_subset_ccscore_harmony_clusters_res2.rds`
- `Science_clusters_res2_top20_markers.csv`
- `Science_clusters_res2_marker_counts.csv`
- `Science_clusters_res2_numeric_cluster_marker_DotPlot.pdf`
- `Science_clusters_res2_UMAP.pdf`
- `Science_clusters_res2_provenance.csv`

`FindAllMarkers` runs on raw clusters in `originalexp`, using the previously
documented `only.pos=TRUE`, `min.pct=0.25`, and `logfc.threshold=0.25` settings.
Exports select up to 20 qualifying markers per cluster by descending average
log2 fold change, breaking ties by gene name. Zero qualifying markers for any
cluster is an error; fewer than 20 is explicitly warned and recorded in the
marker-count table. The old RDS and reports are preserved. Replacing new outputs
requires explicit `--overwrite`, or the full launcher's
`SCIENCE_RES2_OVERWRITE=true`, after review. Resolution-2 markers, raw plots,
and eCN summaries must be regenerated together; production reference paths and
annotations are not switched to these new results automatically.

`fetal cerebellum_science_ACF.R` exports an additional eCN report for **every**
preserved `science_cluster_res2` cluster in `originalexp`, before broad-label
assignment. It uses the existing canonical eCN core, developmental, and shared
supporting panel, with GABA and other excitatory context. It does not assign eCN
labels or alter filtering or the reviewed broad-label map.

To inspect the newly saved `Science_subset_ccscore_harmony_clusters_res2.rds`
without repeating preprocessing, submit from `scripts/` on Cheaha:

```bash
sbatch run_science_ecn_markers.slurm
```

The launcher loads `R/4.4.1-foss-2023b`, starts at the project root for existing
renv activation, exposes Seurat loading errors, and runs the path-only dry-run
before the diagnostic. It requests **two hours, one CPU, 256 GB, express**, retaining
the existing CN diagnostic resource envelope. These are conservative resource
requests, not measured Science requirements. There is no array. Logs use
`logs/Science_eCN_Markers_<job_id>.out` and `.err`. After reviewing existing
diagnostic outputs, `sbatch --export=ALL,SCIENCE_ECN_OVERWRITE=true run_science_ecn_markers.slurm`
enables protected overwrite. No packages are installed/restored by the job.
After completion, collect `JobID`, `Elapsed`, `State`, `MaxRSS`, `ReqMem`,
`AllocCPUS`, and `ExitCode` with `sacct` to refine these resources.

For an existing interactive HPC allocation with the same module and project R
environment, the equivalent direct commands from `scripts/` are:

```bash
Rscript --max-ppsize=500000 'fetal cerebellum_science_ACF.R' --ecn-only --dry-run
Rscript --max-ppsize=500000 'fetal cerebellum_science_ACF.R' --ecn-only
```

The dry-run only inspects paths/dependencies and should take under one minute;
it never loads an RDS. Diagnostic-only execution is estimated at **2–20 minutes**,
including 1–10 minutes to load and 1–10 minutes to summarize/plot. This is a
low-confidence estimate from static inspection, not a Slurm time limit; record
cell count and elapsed time to refine it. Full reprocessing still runs the
existing scaling/PCA/Harmony/neighbors/UMAP pipeline and is much more expensive
(roughly 1–8 hours, very low confidence without measured full-run history).
Use diagnostic-only execution for this marker check.

New outputs go under the project-relative directory
`outputs/references/science/diagnostics/raw_cluster_ecn/`, with prefix
`Science_res2_eCN_`: an unscaled DotPlot,
gene availability, cluster and sample-cluster expression/detection, sample-cluster
SLC17A6/MEIS2 same-cell coexpression, provenance, and interpretation notes.
Existing diagnostic files cause an error unless `--overwrite` is explicitly
provided after review. Diagnostic-only execution does not write RDS objects or
the existing broad-marker/UMAP plots. Missing genes are recorded as unavailable;
missing sample-cluster groups are not treated as zero expression. Expression and
coexpression provide descriptive evidence, with no identity threshold or
statistical test. Inspect detection rates, sample counts, and supporting/context
markers together before reviewing a cluster's identity.

## Recommended maintenance

For future runs, the highest-impact improvements would be:

1. migrate remaining hard-coded sample definitions and paths to the existing
   configuration files;
2. validate Slurm array lengths automatically;
3. eliminate absolute user-specific paths;
4. pin all R and Python dependencies; and
5. encode the stages in a workflow manager such as `targets`, Snakemake, or Nextflow.
