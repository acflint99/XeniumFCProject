# Fetal Cerebellum Xenium Analysis

This directory contains the R and Slurm scripts used to analyze 10x Genomics Xenium spatial-transcriptomics data from the developing human cerebellum. The workflow crops each Xenium field to cerebellar tissue, performs cell-level quality control and Seurat clustering, transfers cell-type labels from published single-cell references, and then performs focused analyses of the ventricular zone (VZ) and rhombic lip (RL) lineages.

The downstream scripts generate spatial maps, cell-type proportions, cross-study comparisons, spatially variable genes, neighborhood enrichment, ligand–receptor results, gene-expression gradients, and trajectory-ready `h5ad` files.

> [!IMPORTANT]
> This is an active research workflow, not a packaged command-line pipeline. Sample lists, dated RDS filenames, Slurm array bounds, cluster selections, and some absolute HPC paths must be checked before each run. The execution order below is inferred from script inputs and outputs.

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

The corresponding `run_*GenePanelSubset.slurm` files submit the panel-subsetting jobs.

### 3. Xenium import, cerebellum cropping, and QC

`xenium_preprocess_selected_samples.R` is the main orchestration script. It sources:

1. `xenium_preprocess_01_crop_cerebellum.R`
2. `xenium_preprocess_02_qc_cells.R`
3. `xenium_preprocess_03_normalize_cluster.R`

`xenium_preprocess_01_crop_cerebellum.R`:

- loads a sample with `Seurat::LoadXenium()`;
- retains cells listed in `cerebellum_cells_stats.csv`;
- adds cell-area metadata;
- saves a cropped spatial feature plot; and
- writes `<sample>_CB.rds` to `outputs/XeniumRDS/`.

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

QC thresholds are stored in `object@misc$QC_thresholds`. Reports are written to `outputs/XeniumQCPlots/`.

There are separate driver variants for different resource strategies:

- `xenium_preprocess_selected_samples.R` – configurable primary driver;
- `xenium_preprocess_selected_samples_parallel.R` – multisample parallel processing;
- `XeniumPreProcess_seq.R` – one-sample sequential processing;
- `xenium_preprocess_split_slides.R` – manifest-driven Slurm-array driver for
  all manually separated biological samples.

### 4. Initial normalization and clustering

`xenium_preprocess_03_normalize_cluster.R` defines `process_xenium_clusters()`. It performs:

- log normalization using the median transcript count as the scale factor;
- 2,000 variable features;
- scaling and 50-component PCA;
- UMAP and Annoy nearest neighbors using PCs 1–50;
- Louvain clustering (`algorithm = 1`) at resolution 1.5; and
- UMAP, global spatial, and faceted spatial cluster plots.

Processed objects are saved as `outputs/Xenium_Res1.5_RDS/<sample>_CB_QC_cluster.rds`.

`SaveRawClusterPlots.R` regenerates raw cluster plots for a selected object.

### 5. Reference-based annotation

`AnchorBasedTransfer_RPCA_res1.5.R` is the primary label-transfer implementation. It:

- balances the reference by downsampling to at most 1,000 cells per reference identity;
- finds genes shared by the reference and Xenium object;
- builds reciprocal-PCA transfer anchors over 30 dimensions;
- transfers `clusters_refined` labels;
- applies a default prediction-score threshold of 0.4; and
- derives cluster-level majority and weighted-vote labels.

The reference and sample are selected explicitly from the command line. The
sample task ID maps to all 34 rows in `config/samples.csv`:

```bash
Rscript scripts/AnchorBasedTransfer_RPCA_res1.5.R --list
Rscript scripts/AnchorBasedTransfer_RPCA_res1.5.R --dry-run Aldinger 1
```

`run_ABT_res1.5.slurm` submits a 34-task array, capped at two concurrent jobs. To
submit all three reference analyses:

```bash
sbatch --job-name=Xen_ABT_Aldinger --export=ALL,REFERENCE=Aldinger scripts/run_ABT_res1.5.slurm
sbatch --job-name=Xen_ABT_Sepp --export=ALL,REFERENCE=Sepp scripts/run_ABT_res1.5.slurm
sbatch --job-name=Xen_ABT_Science --export=ALL,REFERENCE=Science scripts/run_ABT_res1.5.slurm
```

The driver refuses to replace any existing annotation output unless called
with `--overwrite`, or submitted with `ABT_OVERWRITE=true` after review.

Related scripts:

- `XeniumABT_seq.R` – sequential annotation driver;
- `Xenium_ABT_res1.5_comp.R` – merges the three reference comparisons and calculates one cluster-level `consensus_label`;
- `Xenium_ABT_ConsensusLabels.R` – adds consensus labels and PCW metadata to individual objects, makes consensus labels the active identities, and writes consensus UMAP, spatial, and marker DotPlot figures;
- `run_XeniumConsensusLabels.slurm` – submits all 34 consensus-label tasks, capped at three concurrent jobs;
- `Xenium_ABT_GlobalPlot.R` – additional consensus spatial maps;
- `Xenium_ABT_PropGraph.R` – attaches sample metadata and plots consensus cell-type proportions;
- `Xenium_ABT_GlobalPlot_postQC.R` – consensus maps after regional refinement/QC.

The consensus calculation ignores missing and `Unknown` labels and selects the
most frequent label across Aldinger, Sepp, and Science majority/weighted
results. Ties retain the existing comparison-table order. Original reference
labels remain in the object; `consensus_label` is added without overwriting
them. Inspect the task mapping and one task before running:

```bash
Rscript scripts/Xenium_ABT_ConsensusLabels.R --list
Rscript scripts/Xenium_ABT_ConsensusLabels.R --dry-run 1
```

Consensus objects are written to
`outputs/Xenium_ConsensusABT_Res1.5_RDS/<sample>_Consensus_annotated.rds`.
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
  scripts/run_XeniumConsensusLabels.slurm
```

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
post-QC merged processing with `run_xenium_vz_05_reintegrate_post_qc.slurm`.
`run_xenium_vz_07_map_to_samples.slurm` maps the refined identities, and
`run_xenium_vz_08_plot_samples.slurm` launches 34 protected plotting tasks.

Before merging, inspect completeness without loading any Seurat objects:

```bash
Rscript scripts/xenium_vz_02_merge_samples.R --dry-run
```

The merge requires exactly the 34 manifest-defined subset files and writes the
stable output `Merged/Xenium_Merged_VZSubsets.rds` plus a CSV manifest with
each sample's input path, cell count, and PCW.

After that merge exists, inspect processing readiness without loading it:

```bash
Rscript scripts/xenium_vz_03_integrate.R --dry-run
```

Processing writes the stable object
`outputs/XenAld_VZ_Res1.5_RDS/Xenium_VZ_Res1.5.rds`. The QC and post-QC
scripts consume that name. Existing processing RDS or UMAP outputs require an
explicit `--overwrite` rerun after review.

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
submits RL extraction, and
`run_xenium_rl_05_reintegrate_post_qc.slurm` submits post-QC merged processing.

Inspect RL merge readiness with:

```bash
Rscript scripts/xenium_rl_02_merge_samples.R --dry-run
```

The stable merged output is `Merged/Xenium_Merged_RLSubsets.rds`, accompanied
by the corresponding input/cell-count/PCW manifest. Both merge scripts refuse
partial input sets, unexpected top-level RDS files, and accidental output
replacement.

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
`outputs/XenAld_RL_Res1.5_RDS/Xenium_RL_Res1.5.rds`; RL QC and post-QC
processing now consume that path.

### 8. Combined VZ/RL objects

The refined regional labels are combined and analyzed with:

- `Xenium_Combine_Subclusters.R` – combines refined labels in all 34 configured individual samples;
- `Xenium_CombSubclusters_Clean&Merge.R` – cleans and merges per-sample objects;
- `Xenium_CombSubclusters_Process&Plot.R` – processes and visualizes the combined object;
- `Xenium_VZRL_Subclusters_Spatial_Merge.R` – memory-conscious spatial merge;
- `Xenium_VZRL_Subclusters_Spatial_Merge_Processing.R` – Seurat v5 sketch/integration workflow on the merged data;
- `Xenium_VZ&RL_CountPlot_res1.5.R` – combined lineage count plots.

The spatial branch has a separate validated merge because it preserves the
FOV/image data removed from the smaller non-spatial merge. Inspect both stages
before submission:

```bash
Rscript scripts/Xenium_VZRL_Subclusters_Spatial_Merge.R --dry-run
Rscript scripts/Xenium_VZRL_Subclusters_Spatial_Merge_Processing.R --dry-run
```

The merge requires all 34 manifest inputs and writes
`XenAld_VZRL_spatial_merged.rds` plus a cell/image-count manifest. Sketch-based
Harmony processing writes `XenAld_VZRL_spatial_integrated.rds`; this stable name
is also used by the h5ad export and SpaTrack scripts. UMAPs are TIFF-only.
Submit these stages with `run_XeniumVZRLSpatialMerge.slurm` and then
`run_XeniumVZRLSpatialProcess.slurm`.

Before creating the combined per-sample RDS files, spatial TIFFs, and count
tables, inspect input completeness and existing outputs with:

```bash
Rscript scripts/Xenium_Combine_Subclusters.R --dry-run
```

This stage reads sample IDs from `config/samples.csv`, requires all 34 mapped
inputs, and refuses to replace any existing output unless given `--overwrite`.
Combined identities use VZ subclusters first, RL subclusters second, and the
broad consensus label for all remaining cells.
Submit it with `run_XeniumCombineSubclusters.slurm`.

The subsequent clean/merge stage also requires all 34 inputs and removes
spatial images, scale data, and control assays before merging. Because its
filename contains `&`, quote it in the shell:

```bash
Rscript 'scripts/Xenium_CombSubclusters_Clean&Merge.R' --dry-run
```

It writes the stable merged object `XenAld_VZRL_clean_merge.rds` and a
per-sample cell-count/PCW manifest. Combined integration and plotting are
separate operations:

```bash
Rscript 'scripts/Xenium_CombSubclusters_Process&Plot.R' --dry-run
Rscript 'scripts/Xenium_CombSubclusters_Process&Plot.R' --process-only
Rscript 'scripts/Xenium_CombSubclusters_Process&Plot.R' --plots-only
```

The stable processed object is `XenAld_VZRL_clean_merge_processed.rds`.
Plot-only runs verify that the input merge has not changed, and they never
repeat normalization, PCA, Harmony, or UMAP. UMAP outputs remain TIFF-only;
DotPlots, the marker heatmap, and violin plots are saved as TIFF and Cairo PDF.
The corresponding launchers are `run_XeniumCombinedCleanMerge.slurm`,
`run_XeniumCombinedProcess.slurm`, and `run_XeniumCombinedPlots.slurm`.

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
- `Xenium_ABT_PropGraph.R` – proportions by sample metadata such as post-conception week.
- `XenAld_VZ_Compare.R` and `XenSepp_VZ_Compare.R` – average-expression correlations between Xenium VZ clusters and reference VZ clusters.

### 10. Cross-study reference comparisons

- `CrossStudyClusterCorr.R` computes pseudobulk cluster-expression correlations across Aldinger, Sepp, and Science.
- `CrossStudyClusterMarkerComp.R` compares cluster marker programs and produces heatmaps.
- `CrossStudyModuleScoreTransfer.R` transfers dataset-specific marker modules among references.
- `run_MarkerComp.slurm` submits the marker comparison.

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
Rscript scripts/XeniumPreProcess_seq.R
```

For a script that consumes a one-based task ID:

```bash
Rscript scripts/xenium_vz_01_subset_samples.R 1
```

On the configured cluster, submit the matching launcher:

```bash
sbatch scripts/run_xenium_vz_01_subset.slurm
```

To inspect the split-slide task mapping and dry-run one task:

```bash
Rscript scripts/xenium_preprocess_split_slides.R --list
Rscript scripts/xenium_preprocess_split_slides.R --dry-run 1
```

To submit all tasks after reviewing the dry runs:

```bash
sbatch scripts/run_xenium_preprocess_split_slides.slurm
```

The launcher submits 15 tasks and runs at most three concurrently. Each task
processes one biological sample using its shared slide directory and unique
cell-stat CSV from `config/samples.csv`. The driver stops if any expected
output already exists. An intentional rerun requires the explicit
`--overwrite` option, or `SPLIT_PREPROCESS_OVERWRITE=true` when submitting the
launcher. The launcher requests Intel nodes because the project's native R
packages are compiled with Broadwell flags; loading `yaml` on the AMD nodes
tested here causes an illegal-instruction failure.

Before submission, verify:

1. the script's `--list` task mapping;
2. the Slurm `--array` range;
3. input and output RDS filenames;
4. requested memory, time, CPUs, and partition;
5. the project and Python paths; and
6. that the launcher invokes the intended R script.

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

Many intermediate filenames contain manually assigned dates. Downstream scripts often refer to those exact filenames, so renaming or regenerating an object requires updating its consumers.

## Reproducibility notes and known caveats

- Several scripts contain absolute paths under `/home/acflint`, `/data/user/acflint`, or `~/R/Projects/XeniumFCProject`. Replace these when running elsewhere.
- Some legacy scripts still select samples by commenting entries in or out.
  Consensus labelling and VZ/RL subsetting instead use `config/samples.csv`.
- VZ and RL post-QC merged processing use the matching
  `run_xenium_vz_05_reintegrate_post_qc.slurm` and `run_xenium_rl_05_reintegrate_post_qc.slurm`
  launchers.
- Large objects can require 64–500 GB RAM depending on the stage. Do not copy the largest objects across too many `future` workers.
- Random seeds are set in major Seurat stages, but manual cluster renaming and dated file selection remain part of the workflow.
- The ignored `OLD/` directory contains superseded or experimental versions
  and is not part of the active workflow. This includes
  `AnchorBasedTransfer_RCTD_res1.5.R`, `SeppPreProcess.R`, and the broken
  `run_NormCluster1.5.sh` launcher.

## Recommended maintenance

For future runs, the highest-impact improvements would be:

1. migrate remaining hard-coded sample definitions and paths to the existing
   configuration files;
2. replace dated input filenames with stable stage names or a manifest;
3. validate Slurm array lengths automatically;
4. eliminate absolute user-specific paths;
5. pin all R and Python dependencies; and
6. encode the stages in a workflow manager such as `targets`, Snakemake, or Nextflow.
