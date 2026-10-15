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

## Workflow overview

```text
Raw Xenium output + cerebellum cell list
                 |
                 v
Crop cerebellum -> cell QC -> normalization/PCA/UMAP -> clustering (resolution 1.5)
                 |
                 v
Prepare Aldinger / Sepp / Science reference objects and subset to Xenium panel genes
                 |
                 v
RPCA label transfer -> compare references -> consensus-labelled samples
                 |
                 v
Unknown QC audit/review -> apply approved removals -> broad 2-of-3 reference routing
                 |
          +------+------+
          |             |
          v             v
      VZ subset      RL subset
          |             |
          v             v
 merge/integrate/QC/recluster/annotate/map labels back to individual samples
          |             |
          +------+------+
                 |
                 v
combined spatial objects, plots, proportions, Giotto, gradients, SpaTrack/CellRank
```

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

### New samples added 2026-10-14

The manifest now contains 41 analyzed samples (25 single-slide and 16 split-slide
components). The original 34 row/task IDs are preserved. See
[`NEW_SAMPLES_20261014.md`](NEW_SAMPLES_20261014.md) for the targeted seven-sample
commands; historical 34-sample run examples below are not submission instructions
for this cohort. Only section 2 is included from `GZFB124_1_X_G_3__2__1`.

### 1. Xenium gene panel

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
- `xenium_vz_rl_01b_plot_spatial.R` – consensus maps after regional refinement/QC.

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
6. `xenium_vz_06a_build_replicate_marker_evidence.R` – builds sample-aware
   positive and anti-marker candidates; `xenium_vz_06_annotate_subclusters.R`
   retains cell-level Wilcoxon rankings as descriptive review evidence and
   applies explicitly reviewed labels.
7. `xenium_vz_07_map_to_samples.R` – maps refined VZ labels back to all 34 configured samples.
8. `xenium_vz_plot_sample.R` and `xenium_vz_08_plot_samples.R` – manifest-defined per-sample spatial reports.
9. `xenium_vz_09_plot_counts.R` – VZ cell/subcluster counts.

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
neighbor overlap. Expected runtime is approximately 10–45 minutes per branch,
low confidence; each launcher requests 2 hours and 64 GB on `express`.

After the post-QC objects exist and before final label application, generate
sample-aware marker evidence:

```bash
Rscript scripts/xenium_vz_06a_build_replicate_marker_evidence.R --dry-run
sbatch scripts/run_xenium_vz_06a_build_replicate_marker_evidence.slurm
Rscript scripts/xenium_rl_06a_build_replicate_marker_evidence.R --dry-run
sbatch scripts/run_xenium_rl_06a_build_replicate_marker_evidence.slurm
```

These stages compare each raw subcluster with the remaining branch cells within
each biological sample, then summarize direction, magnitude, detection
contrast, and a sample-level sign test. They output separately ranked positive
and anti-marker candidates. Cell-level Wilcoxon tables remain useful for
descriptive ranking, but their P-values are not treated as biological-replicate
inference. The production resolution-0.8 VZ/RL review and label-application
steps now require these five evidence files and verify their input-RDS,
cluster-column, size, and modification-time provenance before proceeding.
Expected runtime is approximately 20–120 minutes per branch,
low confidence; the launchers request 4 hours.

Each VZ Stage-8 global or lineage-highlight plot renders opaque grey context
cells at size 0.5 first and colored subcluster points at size 0.75 second, so
context cannot cover the highlighted cells. Faceted plots omit context cells
and display only the colored cells assigned to each named cluster. It
writes both the all-inclusive global map and a `Global_Spatial_NoOther` map,
a lineage-faceted TIFF, a separate Other-faceted TIFF, marker DotPlots, and a
one-page `Spatial_Report.pdf` containing both global maps, the lineage facets,
and marker DotPlot. Every individual-sample DotPlot and its report page use
the shared `vz_all_compact_markers` panel: exactly the same 42 genes and program
headers as the Stage-6 production DotPlot, with no additional broad-panel
genes. The plotting code constructs all figures before writing any
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
it. Before clustering, stage 3 writes the provenance-checked restart object
`Xenium_VZ_precluster_checkpoint.rds`; immediately after clustering it writes
`Xenium_VZ_postcluster_checkpoint.rds` before creating any plots. Resume mode
automatically prefers the post-clustering checkpoint when it exists, otherwise
it uses the pre-clustering checkpoint. Each stage-3 script also parses its
complete source at startup so a damaged later line fails before loading data or
running analysis. If a job fails after either checkpoint, inspect and resume:

```bash
Rscript scripts/xenium_vz_03_integrate.R --resume-checkpoint --dry-run
sbatch --export=ALL,VZ_INTEGRATE_RESUME_CHECKPOINT=true,VZ_INTEGRATE_OVERWRITE=true \
  scripts/run_xenium_vz_03_integrate.slurm
```

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
The script validates all 34 whole-tissue inputs, protects existing outputs,
and records per-sample removal counts. Post-QC processing reads that manifest
and cannot independently choose a different cluster. After inspecting the
regenerated Stage-4a evidence, the current reviewed VZ decision is to remove
resolution-0.8 clusters `9`, `25`, and `26`. Submit that approved decision with:

```bash
VZ_QC_MODE=apply \
  VZ_QC_REMOVE_CLUSTERS='9,25,26' \
  VZ_QC_OVERWRITE=true \
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
removal materializes 34 filtered objects sequentially and is expected to take
approximately 30-120 minutes. This is a low-confidence estimate based on the
11-minute-45-second failed VZ run and the RL materialization history; the first
successful sequential VZ run will provide the useful measurement. These
estimates exclude queue wait. The 12-hour limit and 256-GB memory request are
retained conservatively until that run completes.

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
`k.param=30`; resolution 0.5 is the preferred Stage-6 review partition, while
0.3, 0.8, and 1 remain available as sensitivity checks. The original-label
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

The full-reintegration launcher requests 16 hours, 4 CPUs, and 600 GB on
`medium`. Its four resolution jobs use `future::multicore`, so the large SNN
graph is shared by forked workers rather than serialized into independent
PSOCK processes. The graph, Louvain algorithm, ten starts, singleton grouping,
resolutions, and random seed are unchanged. Fresh full reintegration is
expected to take approximately 4-8 hours, excluding queue wait; a checkpoint
resume is shorter according to the saved stage. This estimate is low
confidence because clustering still evaluates four resolutions and no complete
run with the new parameters has been measured.

This design follows two observed failure modes. Four PSOCK workers in job
40761237 failed after 2 hours 22 minutes while exporting a 62.89-GiB closure,
at 225.95 GB MaxRSS. A one-worker retry, job 40774031, remained stable at
225.08 GB but timed out after 12 hours during sequential multi-resolution
clustering. Earlier VZ Stage-3 job 40676165 reached approximately 422 GB with
three forked clustering workers, so 600 GB provides headroom for four workers.
Current Cheaha documentation lists up to 754 GB per Intel `medium` node.

This stage writes the stable object
`outputs/xenium/vz/05_post_qc/rds/Xenium_VZ_postQC.rds` and a provenance
manifest. Stage 6 generates post-QC review evidence without assigning
subcluster identities. Resolution 0.5 remains the default; select one of the
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

The approved resolution-0.8 map assigns clusters `11` to VZPs, `12/23` to
Cycling VZPs, `7` to Cycling Glia, `16` to Pre-OPCs, `18` to RL VZ, `19/21`
to Immature/Mature UBCs, and `20` to Microglia. The remaining reviewed labels
cover the Purkinje, GABAergic, glial, OPC, eCN, iCN, and generic cycling
populations. `GCPs` is no longer part of the VZ taxonomy. The filtered VZ UMAP
removes RL VZ, both UBC states, eCN, Microglia, and generic Cycling Cells while
retaining lineage-resolved Cycling VZPs and Cycling Glia. Resolutions 0.3 and
0.5 remain review-only until separate maps are approved.

The main production DotPlot uses the unique union of the two approved compact
marker lists and no additional genes. It displays the 42 genes beneath VZ,
Purkinje, GABA, Glia, OPC, Cycling, RL-VZ, UBC, eCN, and Microglial
program headers. `MKI67` and `TOP2A` occur once under Cycling and are also read
across the Cycling Cells identity row. The VZ-lineage plot groups its selected
genes into VZ, Purkinje, GABA, Glia, OPC, and Cycling programs while retaining
the approved lineage identities as rows. The separate Other plot groups RL-VZ,
combined immature/mature UBC, eCN, generic Cycling Cell, and Microglial markers.
Marker panels and identity membership are defined in `color_palette.R`. Stage 6
refuses to write plots if a requested gene is absent from the Xenium assay, and
its final validation fails unless every required output exists and is nonempty.

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

### 7. RL analysis

The RL branch mirrors the VZ branch:

1. `xenium_rl_01_subset_samples.R` – extracts configured consensus identities per sample.
2. `xenium_rl_02_merge_samples.R` – verifies and incrementally merges all 34 sample subsets.
3. `xenium_rl_03_integrate.R` – validates the complete merge before processing.
4. `xenium_rl_04_review_qc.R` – 4a reviews QC; 4b records the decision and
   conditionally materializes filtered sample objects.
5. `xenium_rl_05_reintegrate_post_qc.R` – reuses Stage-3 results after a
   zero-removal decision, or reintegrates after removals.
6. `xenium_rl_06_annotate_subclusters.R`
7. `xenium_rl_07_map_to_samples.R` – maps refined RL labels back to all 34 configured samples.
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
Rscript scripts/xenium_vz_rl_plot_rl_vz_marker_dotplot.R --list
Rscript scripts/xenium_vz_rl_plot_rl_vz_marker_dotplot.R \
  --branch=RL --dry-run
Rscript scripts/xenium_vz_rl_plot_rl_vz_marker_dotplot.R \
  --branch=VZ --dry-run
sbatch scripts/run_xenium_vz_rl_rl_vz_marker_dotplots.slurm
```

Matched TIFF/PDF DotPlots, marker-availability tables, cluster summaries, and
provenance manifests are written beneath
`outputs/xenium/vz_rl/diagnostics/rl_vz_marker_dotplots/`. Expected runtime is
approximately 15-45 minutes per task and 15-45 minutes total at the configured
array concurrency of two, excluding queue wait. This moderate-confidence
estimate is based on the prior 16-minute RL Stage-6 review, while allowing for
the larger VZ object and high-resolution rendering. The launcher requests two
hours, four CPUs, and 128 GB per task on `express`.

The RL driver also reads all 34 samples from `config/samples.csv` and supports
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

After all 34 one-page Stage-8 reports for a branch are complete, the existing
RL or VZ Stage-8 driver can run in `--combine` mode. That mode validates the
reports and combines them into one 34-page PDF. Pages are ordered by numeric
PCW from `metadata/samples_meta.xlsx`, with `config/samples.csv` order breaking
ties, matching the established merged-report behavior.
It requires the Poppler `pdfinfo` and `pdfunite` commands, verifies that every
input has one page, verifies that the final output has 34 pages, and writes a
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

In Stage 4b, the reviewed resolution-0.8 decision removes clusters
`6,15,24,26,27`. Clusters 6 and 15 have consistently depressed count and
feature distributions, while clusters 24, 26, and 27 contain only 3, 2, and 2
cells. The launcher stores this set as its apply-mode default; review remains
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
  --remove-clusters=6,15,24,26,27 --dry-run
sbatch --export=ALL,RL_QC_MODE=apply,RL_QC_OVERWRITE=true \
  scripts/run_xenium_rl_04_review_qc.slurm
```

`RL_QC_OVERWRITE=true` is required when replacing the earlier no-removal
Stage-4b RDS files and decision manifest. It does not overwrite the Stage-4a
QC evidence. The approved removal contains 69,962 cells in the reviewed object
(about 8.8%); Stage 5 therefore takes the full reintegration path, and all
existing Stage-5 and downstream RL results become stale. As in the VZ Stage-4
launcher, RL Stage 4 starts R with `--max-ppsize=500000`. RL materialization is
also deliberately sequential so each spatial Seurat object is subset in that
main process. Jobs 40698160 and 40698252 showed that the command-line stack
setting was insufficient for the `multisession` path: both failed with
`protect(): protection stack overflow` below 39 GB MaxRSS. Job 40698252 also
demonstrated the Slurm comma-parsing hazard by receiving only cluster `6` from
a comma-separated `--export` override.

The launcher defaults to Stage 4a review mode and requests 12 hours, 1 CPU,
and 256 GB on `short`. Stage 4b records one consistent decision for all 34
samples. A `none` decision writes only a routing manifest pointing to the
original consensus objects and should finish in under 1 minute. A nonzero
removal writes 34 filtered copies sequentially and is expected to take
approximately 20-40 minutes based on job 40698364, which completed all 34
samples in 24 minutes 46 seconds with 36.83 GB MaxRSS. Expanded Stage 4a review
remains approximately 2-8 hours. These estimates exclude queue wait. The
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
zero removed cells, and unchanged input/remaining counts for all 34 samples,
RL stage 5 takes a protected zero-removal fast path. It preserves the stage-3
PCA, Harmony integration, SNN graph, and resolution-0.5/0.8 assignments; copies
the reviewed `umap_harmony` coordinates to the stage-5 `umap_clean` reduction;
and writes only the resolution-0.5/0.8 fast-path plots. It does not calculate a
new resolution-0.3 partition. Any nonzero removal uses the existing full
PCA/Harmony/UMAP/neighbors/reclustering path at resolutions 0.3/0.5/0.8. After
the controlled post-QC diagnostic, the full path uses 30 PCA/Harmony
dimensions, 30 UMAP neighbors, and `k.param=30`. Resolution 0.5 is the preferred
Stage-6 review partition; resolutions 0.3 and 0.8 remain available for coarse
and fine sensitivity checks. The Stage-5 provenance manifest records these
parameters, the selected processing mode, and which reductions and resolutions
were reused or computed. The original-label UMAP uses
`regional_consensus_label_pooled`, preventing sample-qualified Unknown labels
from creating an oversized legend that compresses the plotting panel.

The zero-removal fast path is expected to take approximately 10-45 minutes,
dominated by loading, joining, plotting, and saving the large RL object. Full
post-removal reintegration with the reviewed 30-dimension settings is expected
to take approximately 2.5-4 hours. This moderate-confidence estimate is based
on job 40699023, which completed the 50-dimension path in 3 hours 32 minutes,
and the completed controlled diagnostic. Estimates exclude queue wait. The
launcher requests 12 hours, 4 CPUs, and 256 GB on `short`.
The memory request was raised from 128 GB after full-reintegration job
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

After stage 5, submit `run_xenium_rl_06_annotate_subclusters.slurm` in its
default review mode. Resolution 0.5 remains the default, while
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
The old identity map is deliberately not reused. The reviewed resolution-0.8
map is now configured in Stage 6 from the current post-QC UMAP, top-50 marker
table, focused lineage DotPlot, and RL-VZ diagnostic DotPlot. The map groups
clusters `8/12` as cycling RL-SVZ/UBC progenitors, `9/15` as proliferating
GCPs, `3/5/10` as migrating GCs, `4/7/14` as differentiating GCs, and
`13/21/22` as maturing PCs. It retains the reviewed off-lineage identities
`eCN` (17), `GABA Progenitors` (18), `Microglia` (19), and `Endothelial` (20)
instead of forcing them into an RL lineage. Shared taxonomy and palette checks
prevent any reviewed label from becoming `NA` during factor conversion. The
combined count plots report `eCN` as its own lineage and place `Microglia` with
the other non-neural/off-lineage populations rather than counting either as
UBC. Production RL-branch plots use only the 14 identities present in that
approved map; the RL-branch's obsolete `RL VZ`, `Intermediate Progenitors`,
and `Ependymal Cells` levels are not displayed. This does not remove the
separately reviewed `RL VZ` identity from the VZ resolution-0.8 taxonomy. The
production output retains the merged DotPlot and adds separate compact RL
lineage and Other DotPlots. The lineage panel uses three existing production
markers for each RL, UBC, and granule program across the corresponding lineage
identities. The Other panel uses three markers each for BG, Maturing PCs, eCN,
GABA Progenitors, Microglia, and Endothelial. `HES1`, `BCAN`, nonspecific `PAX6`, weak per-sample
Maturing-PC marker `FOXP2`, and assay-absent `GRIA4` remain omitted.
Both RL UMAPs use centered white boxes for every cluster label, without leader
lines or changes to UMAP coordinates or identities.

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

Application is refused if the post-QC input changed after review. Review mode,
including the focused DotPlot, is expected to take approximately 15-45 minutes
based on job 40792498, which completed the earlier resolution-0.8 marker review
in 16 minutes 9 seconds with 69.09 GB MaxRSS. This estimate is moderate
confidence; the added DotPlot should be small relative to marker testing but
has not yet been timed on the full object. Full resolution-0.8 annotation is
expected to take approximately 10-20 minutes based on completed job 40865299,
which ran in 11 minutes 32 seconds with 102.77 GB MaxRSS. The launcher requests
12 hours, 4 CPUs, and 128 GB on `short`.

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
`run_xenium_rl_08_plot_samples.slurm` for the 34-task plotting array. Inspect one
plotting task first with:

```bash
Rscript scripts/xenium_rl_08_plot_samples.R --dry-run 1
```

Stages 5 and 6 write the stable objects `Xenium_RL_postQC.rds` and
`Xenium_RL_subclusters.rds` in their numbered stage directories. Stage
9 intentionally uses the same focused nine-sample, one-sample-per-PCW cohort
as VZ, including `GZFB_22_X_G_4` for PCW19. Inspect it
with `Rscript scripts/xenium_rl_09_plot_counts.R --dry-run` and submit it with
`run_xenium_rl_09_plot_counts.slurm`.

### 7a. Residual Other/unresolved branch

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

Rscript scripts/xenium_residual_04_review_subclusters.R --dry-run
sbatch scripts/run_xenium_residual_04_review_subclusters.slurm
```

Stage 1 records a status table for every manifest sample and writes an RDS only
when retained residual cells are present. Stage 2 validates all 34 statuses and
merges only declared nonempty subsets. Stage 3 reuses the approved regional
30-dimensional Harmony and resolution 0.5/0.8/1.0 settings. Stage 4 generates
replicate-aware positive and anti-marker evidence plus a blank label-review
template. Copy the completed reviewed rows to
`config/reviewed_residual_subcluster_labels.csv`; every cluster requires an
explicit decision, confidence, evidence summaries, notes, and `approved=TRUE`.
`final_label=Unknown` is allowed when no defensible winner exists. Then run:

The approved resolution-0.8 residual map assigns clusters `0` to `eCN`, `1`
to `BG`, `2` to `Immature UBCs`, `3` to `Microglia`, `4` to `Pericytes`, `5`
to `Mature UBCs`, `6` to `Unknown`, `7` to `GABA Progenitors`, and `8` to
`Cycling Cells`. The retained `Unknown` records the absence of a defensible
canonical lineage winner rather than forcing a label from generic immature
neuronal markers.

```bash
Rscript scripts/xenium_residual_05_apply_reviewed_labels.R --dry-run
sbatch scripts/run_xenium_residual_05_apply_reviewed_labels.slurm

Rscript scripts/xenium_residual_06_map_to_samples.R --dry-run
sbatch scripts/run_xenium_residual_06_map_to_samples.slurm
```

Stage 6 maps the reviewed residual identities back to all 34 reviewed
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

### 8. Combined VZ/RL/residual objects

The refined regional labels are combined and analyzed with:

- `xenium_vz_rl_01_combine_labels.R` – combines refined labels in all 34 configured individual samples;
- `xenium_vz_rl_02_merge_samples.R` – cleans and merges per-sample objects;
- `xenium_vz_rl_03_process_and_plot.R` – processes and visualizes the combined object;
- `xenium_vz_rl_spatial_01_merge.R` – memory-conscious spatial merge;
- `xenium_vz_rl_spatial_02_integrate.R` – Seurat v5 sketch/integration workflow on the merged data;
- `xenium_vz_rl_03b_plot_cluster_counts.R` – combined lineage count plots.

Combined outputs are organized under `outputs/xenium/vz_rl/`:

- `01_combined_labels/` contains all 34 labelled sample objects, per-sample spatial plots, and count tables;
- `02_merged/` contains the cleaned sample objects, merged RDS, and merge manifest;
- `03_processed/` contains the integrated object, provenance manifest, summary plots, and cluster-count plots;
- `spatial/01_merged/` and `spatial/02_integrated/` contain the spatial-preserving branch.

The earlier `xenium_vz_rl_00_compare_unknown_annotations.R` dual-branch
comparison and `config/reviewed_unknown_relabels.csv` are retained only as
historical/sensitivity tools. They are not production prerequisites after
exclusive reference routing. Do not run the comparison unless a reviewer asks
for the prior include-in-both-branches sensitivity analysis.

The spatial branch has a separate validated merge because it preserves the
FOV/image data removed from the smaller non-spatial merge. Inspect both stages
before submission:

```bash
Rscript scripts/xenium_vz_rl_spatial_01_merge.R --dry-run
Rscript scripts/xenium_vz_rl_spatial_02_integrate.R --dry-run
```

The merge requires all 34 manifest inputs and writes
`spatial/01_merged/rds/Xenium_VZRL_spatial_merged.rds` plus a cell/image and
regional-assignment count manifest under `spatial/01_merged/tables/`.
Sketch-based Harmony processing
writes `spatial/02_integrated/rds/Xenium_VZRL_spatial_integrated.rds`. UMAPs
are TIFF-only. Giotto, h5ad export, and SpaTrack scripts remain outside the
current scope and still require path review before they are reactivated.
Submit these stages with `run_xenium_vz_rl_spatial_01_merge.slurm` and then
`run_xenium_vz_rl_spatial_02_integrate.slurm`.

Before creating the combined per-sample RDS files, spatial TIFFs, and count
tables, inspect input completeness and existing outputs with:

```bash
Rscript scripts/xenium_vz_rl_01_combine_labels.R --dry-run
```

This stage reads sample IDs from `config/samples.csv`, requires the reviewed
whole-tissue input manifest, the reference-route manifest, and all 34
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
The combined-label stage is expected to take approximately 4–12 hours for all
manifest samples (low confidence, based on static inspection of 34 sequential
three-object joins, 136 high-resolution spatial plots, and 34 uncompressed RDS
writes; no matching `sacct` record is present in `slurm_resource_history.psv`).
The existing launcher requests 24 hours on `medium`; queue wait is excluded.
Collecting `Elapsed` and `MaxRSS` from the first reviewed run will substantially
improve this estimate.

Applying or changing the reviewed Unknown QC decisions, reference-route
manifest, or reviewed residual labels makes every existing combined Stage-1
output and all downstream VZ/RL/residual merge, processed, spatial, count, and
trajectory-derived outputs stale. Regenerate them only after reviewing the
Unknown-QC, route, and residual-label tables. A rollback is the prior reviewed
QC and residual-label CSVs plus the corresponding protected Stage-1 outputs;
the whole-tissue consensus and mapped branch objects are never modified.

The subsequent clean/merge stage also requires all 34 inputs and removes
spatial images, scale data, and control assays before merging:

```bash
Rscript scripts/xenium_vz_rl_02_merge_samples.R --dry-run
```

It writes the stable merged object
`02_merged/rds/Xenium_VZRL_merged.rds` and a per-sample cell-count, PCW, and
regional-assignment manifest under `02_merged/tables/`. Combined integration and plotting are
separate operations:

```bash
Rscript scripts/xenium_vz_rl_03_process_and_plot.R --dry-run
Rscript scripts/xenium_vz_rl_03_process_and_plot.R --process-only
Rscript scripts/xenium_vz_rl_03_process_and_plot.R --plots-only
```

The stable processed object is
`03_processed/rds/Xenium_VZRL_processed.rds`.
Plot-only runs verify that the input merge has not changed, and they never
repeat normalization, PCA, Harmony, or UMAP. UMAP outputs remain TIFF-only;
DotPlots, the marker heatmap, and violin plots are saved as TIFF and Cairo PDF.
The corresponding launchers are `run_xenium_vz_rl_02_merge_samples.slurm`,
`run_xenium_vz_rl_03_process.slurm`, and `run_xenium_vz_rl_03_plots_only.slurm`.

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

## Software requirements

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

## Recommended maintenance

For future runs, the highest-impact improvements would be:

1. migrate remaining hard-coded sample definitions and paths to the existing
   configuration files;
2. validate Slurm array lengths automatically;
3. eliminate absolute user-specific paths;
4. pin all R and Python dependencies; and
5. encode the stages in a workflow manager such as `targets`, Snakemake, or Nextflow.
