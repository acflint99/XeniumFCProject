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

- **VZ-related cells:** glia, GABAergic cells, Purkinje cells, and OPCs
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
samples, and the six multi-sample slide directories. They are currently
validated metadata scaffolding; not every analysis script reads them yet. From
the project root, validate their structure and referenced input files with:

```bash
Rscript scripts/validate_config.R --check-files
```

## Main analysis stages

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
`outputs/xenium/preprocess/03f_resolution_consistency/weighted_2of3/`.

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

This isolated workflow changes only the whole-tissue clustering resolution from
4 to 5. It retains the saved SNN graph, three reference-transfer settings,
weighted 2-of-3 consensus rules, seeds, marker sets, and fixed DotPlot scales.
Clustering outputs are written beneath
`outputs/xenium/preprocess/03i_resolution5_all_samples/`; annotation and report
outputs are written beneath `outputs/xenium/annotation/resolution5_all_samples/`.
Existing resolution-1.5, resolution-4, and pilot results are not replaced.

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
  scripts/run_xenium_preprocess_03i_resolution5_all_samples.slurm)
cluster_job=${cluster_submit%%;*}

transfer_jobs=()
for reference in Aldinger Sepp Science; do
  transfer_submit=$(sbatch --parsable \
    --dependency="afterok:${cluster_job}" \
    --array="1-34%6" \
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
`outputs/xenium/annotation/resolution5_all_samples/03_consensus_labels_weighted_2of3/plots/Xenium_Consensus_Res5_Weighted2of3_All_Samples_Report.pdf`.
A matched report with every `Unknown-##` cell removed is written as
`Xenium_Consensus_Res5_Weighted2of3_All_Samples_KnownOnly_Report.pdf` in the
same directory.
Report facet panels use point size 0.01, a one-pixel glyph, and alpha 0.55 to
make reduced overplotting visible after rasterization.

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
This conditional output requires `writexl` in the active project `renv`; the
production script never installs it.

Expected runtime per clustering task is 4-10 minutes and approximately 35-90
minutes total at concurrency four. Across 170 completed resolution-5 transfer
tasks from jobs 39875880-39875882, 39912562, and 39943655, observed task times
were 7.3-47.3 minutes with per-array means of 18.7-19.8 minutes. Aldinger,
Sepp, and Science are each capped at six concurrent tasks; each array should
finish in approximately 2-3 hours and all three may do so concurrently when
all 18 tasks can be scheduled. At maximum simultaneous use, the transfer
arrays can reserve up to 1,152 GB cluster-wide. Observed maximum RSS was
49.28 GB, consistently peaking for manifest task 6 (`GZFB4_X_G`), so the
64-GB per-task request is retained.
Consensus-table construction should take under 2 minutes. Consensus
application, paired plot sets, and conditional Unknown marker tests are
expected to take 15-90 minutes per task and approximately 1.5-9 hours total at
concurrency six. Paired report pages are expected to take 20-90 minutes per
task and approximately 2-9 hours total at concurrency six; merging both reports
should take under 15 minutes. These revised estimates are low confidence and
come from prior single-version plotting history plus static inspection of the
new doubled plotting and marker workload; the next `sacct` results should be
used to refine them. Requested limits remain 1 hour for
clustering, 4 hours per transfer, 30
minutes for table construction, 2 hours per consensus/report-page task, and 30
minutes for the merge.

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

### 5. Reference-based annotation

`xenium_annotate_01_label_transfer_rpca.R` is the primary label-transfer implementation. It:

- balances the reference by downsampling to at most 1,000 cells per reference identity;
- finds genes shared by the reference and Xenium object;
- builds reciprocal-PCA transfer anchors over 30 dimensions;
- transfers `clusters_refined` labels;
- applies a default prediction-score threshold of 0.4; and
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
`outputs/xenium/annotation/01_label_transfer/<reference>/`. Combined voting
tables are written to `outputs/xenium/annotation/02_consensus/tables/`.

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
`outputs/xenium/annotation/03_consensus_labels/rds/<sample>_Consensus_annotated.rds`.
Consensus UMAP, spatial, and marker plots are written together under
`outputs/xenium/annotation/03_consensus_labels/plots/`; only the separate
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
`outputs/xenium/annotation/03_consensus_labels/plots/Xenium_Consensus_All_Samples_Report.pdf`
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

The VZ branch generally follows this order:

1. `xenium_vz_01_subset_samples.R` – extracts configured consensus identities per sample.
2. `xenium_vz_02_merge_samples.R` – verifies and incrementally merges all 34 sample subsets after removing spatial overhead.
3. `xenium_vz_03_integrate.R` – validates the complete merge, then runs normalization, PCA, Harmony integration, UMAP, and clustering.
4. `xenium_vz_04_review_qc.R` – merged-cluster QC summaries and cell flags.
5. `xenium_vz_05_reintegrate_post_qc.R` – removes failed cells and recomputes integration and reductions.
6. `xenium_vz_06_annotate_subclusters.R` – VZ subcluster annotation and marker analysis.
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
manifest under `02_merged/tables/` with each sample's input path, cell count,
and PCW.

After that merge exists, inspect processing readiness without loading it:

```bash
Rscript scripts/xenium_vz_03_integrate.R --dry-run
```

Processing writes the stable object
`outputs/xenium/vz/03_integrated/rds/Xenium_VZ_Res1.5.rds`. The QC and
post-QC scripts consume that name. Existing processing RDS or UMAP outputs
require an explicit `--overwrite` rerun after review.

VZ QC is an explicit review gate. First inspect readiness, then generate the
QC summary, combined violin PDF, marker table, and review manifest:

```bash
Rscript scripts/xenium_vz_04_review_qc.R --dry-run
Rscript scripts/xenium_vz_04_review_qc.R --qc-only
```

After reviewing those outputs, record the decision explicitly. For example,
remove clusters 7 and 8 with `--remove-clusters=7,8`, or preserve every cell
with `--remove-clusters=none`. The script validates all 34 whole-tissue inputs,
protects existing outputs, and records per-sample removal counts. Post-QC
processing reads that manifest and cannot independently choose a different
cluster.

Inspect post-QC integration without loading the merged Seurat object:

```bash
Rscript scripts/xenium_vz_05_reintegrate_post_qc.R --dry-run
```

This stage writes the stable object
`outputs/xenium/vz/05_post_qc/rds/Xenium_VZ_postQC_Res1.5.rds` and a provenance
manifest. Stage 6 validates the complete cluster-to-label mapping before it
writes `outputs/xenium/vz/06_subclusters/rds/Xenium_VZ_subclusters_Res1.5.rds`:

```bash
Rscript scripts/xenium_vz_06_annotate_subclusters.R --dry-run
```

Before mapping refined labels back to the whole-tissue objects, inspect all
expected inputs and outputs without loading Seurat objects:

```bash
Rscript scripts/xenium_vz_07_map_to_samples.R --dry-run
```

The mapping stage requires exactly one input for every sample in
`config/samples.csv`, rejects unexpected RDS files in its input directory, and
refuses to replace existing outputs unless given `--overwrite`. It verifies
that every cell in the master VZ object maps exactly once and writes
`Xenium_VZ_Mapping_manifest.csv` with per-sample cell counts and PCW.

Inspect the VZ plotting task map and one sample without creating plots:

```bash
Rscript scripts/xenium_vz_08_plot_samples.R --list
Rscript scripts/xenium_vz_08_plot_samples.R --dry-run 1
```

Stage 9 intentionally compares a focused eight-sample developmental cohort;
it does not use all 34 samples. Inspect its protected outputs with
`Rscript scripts/xenium_vz_09_plot_counts.R --dry-run` and submit it with
`run_xenium_vz_09_plot_counts.slurm`.

### 7. RL analysis

The RL branch mirrors the VZ branch:

1. `xenium_rl_01_subset_samples.R` – extracts configured consensus identities per sample.
2. `xenium_rl_02_merge_samples.R` – verifies and incrementally merges all 34 sample subsets.
3. `xenium_rl_03_integrate.R` – validates the complete merge before processing.
4. `xenium_rl_04_review_qc.R`
5. `xenium_rl_05_reintegrate_post_qc.R`
6. `xenium_rl_06_annotate_subclusters.R`
7. `xenium_rl_07_map_to_samples.R` – maps refined RL labels back to all 34 configured samples.
8. `xenium_rl_plot_sample.R` and `xenium_rl_08_plot_samples.R` – manifest-defined per-sample reports.
9. `xenium_rl_09_plot_counts.R`

The RL driver also reads all 34 samples from `config/samples.csv` and supports
the same inspection and overwrite protections. `run_xenium_rl_01_subset.slurm`
submits RL extraction. Use `run_xenium_rl_02_merge_samples.slurm`,
`run_xenium_rl_03_integrate.slurm`,
`run_xenium_rl_05_reintegrate_post_qc.slurm`, and
`run_xenium_rl_06_annotate_subclusters.slurm` for the corresponding merged
stages.

RL outputs mirror the same numbered layout under `outputs/xenium/rl/`, from
`01_subsets` through `09_cluster_counts`. This makes corresponding VZ and RL
stages directly comparable without changing any output filenames.

Inspect RL merge readiness with:

```bash
Rscript scripts/xenium_rl_02_merge_samples.R --dry-run
```

The stable merged output is
`outputs/xenium/rl/02_merged/rds/Xenium_Merged_RLSubsets.rds`, accompanied by
the corresponding input/cell-count/PCW manifest under `02_merged/tables/`.
Both merge scripts refuse partial input sets, unexpected top-level RDS files,
and accidental output replacement.

RL QC uses the same explicit review gate:

```bash
Rscript scripts/xenium_rl_04_review_qc.R --dry-run
Rscript scripts/xenium_rl_04_review_qc.R --qc-only
```

After reviewing the RL evidence, use `--remove-clusters=<IDs>` or
`--remove-clusters=none`. The recorded decision is applied consistently to all
34 individual objects and the merged post-QC object.

Inspect processing readiness with:

```bash
Rscript scripts/xenium_rl_03_integrate.R --dry-run
```

Inspect the complete RL mapping inputs and protected outputs with:

```bash
Rscript scripts/xenium_rl_07_map_to_samples.R --dry-run
```

The RL mapping stage applies the same completeness, unexpected-file,
overwrite, barcode-matching, and manifest checks as the VZ mapping stage. It
also requires the VZ subcluster metadata inherited from the preceding branch
and writes `Xenium_RL_Mapping_manifest.csv`.

Use `run_xenium_rl_07_map_to_samples.slurm` for mapping and
`run_xenium_rl_08_plot_samples.slurm` for the 34-task plotting array. Inspect one
plotting task first with:

```bash
Rscript scripts/xenium_rl_08_plot_samples.R --dry-run 1
```

The stable processed object is
`outputs/xenium/rl/03_integrated/rds/Xenium_RL_Res1.5.rds`; RL QC and post-QC
processing consume that path.

Inspect the protected post-QC and subcluster stages with:

```bash
Rscript scripts/xenium_rl_05_reintegrate_post_qc.R --dry-run
Rscript scripts/xenium_rl_06_annotate_subclusters.R --dry-run
```

They write the stable objects `Xenium_RL_postQC_Res1.5.rds` and
`Xenium_RL_subclusters_Res1.5.rds` in their numbered stage directories. Stage
9 intentionally uses the same focused eight-sample cohort as VZ. Inspect it
with `Rscript scripts/xenium_rl_09_plot_counts.R --dry-run` and submit it with
`run_xenium_rl_09_plot_counts.slurm`.

### 8. Combined VZ/RL objects

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

The spatial branch has a separate validated merge because it preserves the
FOV/image data removed from the smaller non-spatial merge. Inspect both stages
before submission:

```bash
Rscript scripts/xenium_vz_rl_spatial_01_merge.R --dry-run
Rscript scripts/xenium_vz_rl_spatial_02_integrate.R --dry-run
```

The merge requires all 34 manifest inputs and writes
`spatial/01_merged/rds/XenAld_VZRL_spatial_merged.rds` plus a cell/image-count
manifest under `spatial/01_merged/tables/`. Sketch-based Harmony processing
writes `spatial/02_integrated/rds/XenAld_VZRL_spatial_integrated.rds`. UMAPs
are TIFF-only. Giotto, h5ad export, and SpaTrack scripts remain outside the
current scope and still require path review before they are reactivated.
Submit these stages with `run_xenium_vz_rl_spatial_01_merge.slurm` and then
`run_xenium_vz_rl_spatial_02_integrate.slurm`.

Before creating the combined per-sample RDS files, spatial TIFFs, and count
tables, inspect input completeness and existing outputs with:

```bash
Rscript scripts/xenium_vz_rl_01_combine_labels.R --dry-run
```

This stage reads sample IDs from `config/samples.csv`, requires all 34 mapped
inputs, and refuses to replace any existing output unless given `--overwrite`.
Combined identities use VZ subclusters first, RL subclusters second, and the
broad consensus label for all remaining cells.
Submit it with `run_xenium_vz_rl_01_combine_labels.slurm`.

The subsequent clean/merge stage also requires all 34 inputs and removes
spatial images, scale data, and control assays before merging:

```bash
Rscript scripts/xenium_vz_rl_02_merge_samples.R --dry-run
```

It writes the stable merged object
`02_merged/rds/XenAld_VZRL_clean_merge.rds` and a per-sample cell-count/PCW
manifest under `02_merged/tables/`. Combined integration and plotting are
separate operations:

```bash
Rscript scripts/xenium_vz_rl_03_process_and_plot.R --dry-run
Rscript scripts/xenium_vz_rl_03_process_and_plot.R --process-only
Rscript scripts/xenium_vz_rl_03_process_and_plot.R --plots-only
```

The stable processed object is
`03_processed/rds/XenAld_VZRL_clean_merge_processed.rds`.
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
- Large objects can require 64–500 GB RAM depending on the stage. Do not copy the largest objects across too many `future` workers.
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
