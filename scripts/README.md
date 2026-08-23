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
RPCA label transfer -> broad annotated samples
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

`XeniumGenePanel.R` extracts target genes from a Xenium panel JSON file. The resulting gene list is used to restrict the reference datasets to genes measurable by Xenium.

### 2. Reference preparation

The reference workflows harmonize cell-type labels, prepare Seurat objects, restrict them to the Xenium panel, and recompute dimensional reductions.

| Reference | Initial preparation | Xenium-panel subset | Additional analysis |
|---|---|---|---|
| Aldinger | `AldingerPreProcess.R`, `AldingerPreProcess2.R` | `AldingerGenePanelSubset.R` | `AldingerGenePanelSubset_Plots.R`, `Aldinger_VZ_Subset.R` |
| Sepp | `SeppPreProcess_v2.R` | `SeppGenePanelSubset.R` | `Sepp_VZ_Subset.R` |
| Science | `SciencePreProcess.R` | `ScienceGenePanelSubset.R` | — |

The corresponding `run_*GenePanelSubset.slurm` files submit the panel-subsetting jobs.

### 3. Xenium import, cerebellum cropping, and QC

`XeniumPreProcess.R` is the main orchestration script. It sources:

1. `XeniumCropCerebellum.R`
2. `XeniumQC.R`
3. `XeniumNormCluster_res1.5.R`

`XeniumCropCerebellum.R`:

- loads a sample with `Seurat::LoadXenium()`;
- retains cells listed in `cerebellum_cells_stats.csv`;
- adds cell-area metadata;
- saves a cropped spatial feature plot; and
- writes `<sample>_CB.rds` to `outputs/XeniumRDS/`.

For combined-slide data, `XeniumCropCerebellum.R` accepts the shared input
directory and sample-specific cell-stat CSV separately. The manifest-driven
`XeniumPreProcess_split_slides.R` driver uses this interface automatically.
`XeniumCropCerebellum_duo.R` remains as the earlier manual version.

`XeniumQC.R` calculates:

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

- `XeniumPreProcess.R` – configurable primary driver;
- `XeniumPreProcess_parallel.R` – multisample parallel processing;
- `XeniumPreProcess_seq.R` – one-sample sequential processing;
- `XeniumPreProcess_seq_duo.R` – earlier manual combined-slide driver; and
- `XeniumPreProcess_split_slides.R` – manifest-driven Slurm-array driver for
  all manually separated biological samples.

### 4. Initial normalization and clustering

`XeniumNormCluster_res1.5.R` defines `process_xenium_clusters()`. It performs:

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

The sample to process is selected by a one-based command-line task ID. `run_ABT_res1.5.sh` submits the current array.

Related scripts:

- `XeniumABT_seq.R` – sequential annotation driver;
- `Xenium_ABT_res1.5_comp.R` – merges comparison tables from the three reference runs;
- `Xenium_ABT_GlobalPlot.R` – broad annotation maps;
- `Xenium_ABT_PropGraph.R` – attaches sample metadata and plots cell-type proportions;
- `Xenium_ABT_GlobalPlot_postQC.R` – annotation maps after regional refinement/QC.

### 6. VZ analysis

The VZ branch generally follows this order:

1. `Xenium_VZ_Subset_res1.5.R` – extracts target broad identities per sample.
2. `Xenium_VZ_Merge_res1.5.R` – merges sample subsets after removing unnecessary spatial overhead.
3. `Xenium_VZ_Merge_Processing_res1.5.R` – normalization, PCA, Harmony integration, UMAP, and clustering.
4. `Xenium_VZ_Merge_QC_res1.5.R` – merged-cluster QC summaries and cell flags.
5. `Xenium_VZ_Merge_postQC_Processing_res1.5.R` – removes failed cells and recomputes integration and reductions.
6. `Xenium_VZ_Analysis_res1.5.R` – VZ subcluster annotation and marker analysis.
7. `Xenium_VZ_Mapping_res1.5.R` – maps refined VZ labels back to individual samples.
8. `Xenium_VZ_SamplePlots_res1.5.R` and `XeniumVZSamplePlots_Loop_res1.5.R` – per-sample spatial reports.
9. `Xenium_VZ_CountPlot_res1.5.R` – VZ cell/subcluster counts.

Submit per-sample extraction with `run_XeniumVZSubset_res1.5.sh`. `run_VZSamplePlots.slurm` launches the plotting stage.

### 7. RL analysis

The RL branch mirrors the VZ branch:

1. `Xenium_RL_Subset_res1.5.R`
2. `Xenium_RL_Merge_res1.5.R`
3. `Xenium_RL_Merge_Processing_res1.5.R`
4. `Xenium_RL_Merge_QC_res1.5.R`
5. `Xenium_RL_Merge_postQC_Processing_res1.5.R`
6. `Xenium_RL_Analysis_res1.5.R`
7. `Xenium_RL_Mapping_res1.5.R`
8. `Xenium_RL_SamplePlots_res1.5.R` and `Xenium_RL_SamplePlots_Loop_res1.5.R`
9. `Xenium_RL_CountPlot_res1.5.R`

`run_XeniumRLSubset_res1.5.sh` submits RL extraction. `run_XeniumRLAnalysis.slurm` currently invokes a post-QC processing script and should be checked against the intended stage before submission.

### 8. Combined VZ/RL objects

The refined regional labels are combined and analyzed with:

- `Xenium_Combine_Subclusters.R` – combines refined labels in individual samples;
- `Xenium_CombSubclusters_Clean&Merge.R` – cleans and merges per-sample objects;
- `Xenium_CombSubclusters_Process&Plot.R` – processes and visualizes the combined object;
- `Xenium_VZRL_Subclusters_Spatial_Merge.R` – memory-conscious spatial merge;
- `Xenium_VZRL_Subclusters_Spatial_Merge_Processing.R` – Seurat v5 sketch/integration workflow on the merged data;
- `Xenium_VZ&RL_CountPlot_res1.5.R` – combined lineage count plots.

`run_VZRL_Subcluster_Merge.sh` submits the merge stage.

### 9. Spatial and trajectory analyses

#### Giotto

`Xenium_Giotto_Analysis.R` and `Xenium_Giotto_Analysis_array.R` convert annotated Seurat objects to Giotto objects and calculate:

- Delaunay spatial networks;
- spatially variable genes with `binSpect()`;
- cell-proximity enrichment; and
- spatial ligand–receptor interactions using Xenium-detectable pairs from the human CellChat database.

`Xenium_Giotto_BroadCluster_Analysis_array.R` runs the corresponding broad-cluster analysis. Plotting helpers include `Xenium_Giotto_SVGPlots.R`, `Xenium_Giotto_SVGCompPlots.R`, and `Xenium_Giotto_TimePlots.R`.

Launchers are `run_GiottoAnalysis.sh`, `run_GiottoAnalysis_array.sh`, and `run_GiottoBroadAnalysis_array.sh`.

#### SpaTrack and CellRank export

- `Xenium_VZRL_Subclusters_VZ_h5adExport.R` exports VZ-lineage expression, metadata, embeddings, and spatial coordinates to `h5ad`.
- `Xenium_VZRL_Subclusters_VZ_SpaTrack.R` prepares/runs the per-sample SpaTrack workflow through `reticulate`.
- `Xenium_VZRL_Subclusters_VZ_SpaTrack_Plotting.R` imports pseudotime and creates spatial plots.
- `run_VZSpaTrack_array.sh` submits per-sample trajectory jobs.

These scripts require configured Python/Conda environments and currently contain environment-specific paths.

#### Other downstream analyses

- `Xenium_GeneExpGradient_Analysis.R` – spatial expression gradients across selected lineages or regions.
- `Xenium_ABT_PropGraph.R` – proportions by sample metadata such as post-conception week.
- `XenAld_VZ_Compare.R` and `XenSepp_VZ_Compare.R` – average-expression correlations between Xenium VZ clusters and reference VZ clusters.

### 10. Cross-study reference comparisons

- `CrossStudyClusterCorr.R` computes pseudobulk cluster-expression correlations across Aldinger, Sepp, and Science.
- `CrossStudyClusterMarkerComp.R` compares cluster marker programs and produces heatmaps.
- `CrossStudyModuleScoreTransfer.R` transfers dataset-specific marker modules among references.
- `run_MarkerComp.sh` submits the marker comparison.

## Software requirements

The configured Slurm environment uses R 4.4.1 from the
`R/4.4.1-foss-2023b` module. Frequently used R packages include:

- `Seurat`
- `harmony`
- `future` and `future.apply`
- `dplyr`, `tidyr`, `purrr`, and `tibble`
- `ggplot2`, `patchwork`, `ggh4x`, `randomcoloR`, and `Cairo`
- `ComplexHeatmap`, `pheatmap`, `circlize`, and `viridis`
- `Giotto`
- `reticulate`
- `readxl` and `jsonlite`
- `AnnotationDbi` and `org.Hs.eg.db`

Some scripts call `source("renv/activate.R")`, and several launchers run `renv::restore()`. Restore the environment once before launching large arrays; simultaneous array jobs should not compete for an `renv` library lock.

Giotto and trajectory scripts additionally require the corresponding Python/Conda environments. Giotto downloads the human CellChat database at runtime, so compute nodes need network access or the database must be cached locally.

## Running the workflow

Run R scripts from the **project root**, not from `scripts/`, so `here()` resolves the expected data and output paths.

```bash
cd /path/to/XeniumFCProject
Rscript scripts/XeniumPreProcess_seq.R
```

For a script that consumes a one-based task ID:

```bash
Rscript scripts/Xenium_VZ_Subset_res1.5.R 1
```

On the configured cluster, submit the matching launcher:

```bash
sbatch scripts/run_XeniumVZSubset_res1.5.sh
```

To inspect the split-slide task mapping and dry-run one task:

```bash
Rscript scripts/XeniumPreProcess_split_slides.R --list
Rscript scripts/XeniumPreProcess_split_slides.R --dry-run 1
```

To submit all tasks after reviewing the dry runs:

```bash
sbatch scripts/run_XeniumPreProcess_split_slides.sh
```

The launcher submits 15 tasks and runs at most three concurrently. Each task
processes one biological sample using its shared slide directory and unique
cell-stat CSV from `config/samples.csv`. The driver stops if any expected
output already exists. An intentional rerun requires the explicit
`--overwrite` option, or `SPLIT_PREPROCESS_OVERWRITE=true` when submitting the
launcher.

Before submission, verify:

1. the uncommented R sample list;
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
- Active sample lists are commonly selected by commenting entries in or out. Treat the Slurm array bounds and R sample vectors as a coupled configuration.
- `run_NormCluster1.5.sh` calls `XeniumNormCluster_res1.5.R`, but that file currently defines a function and does not read the array task ID. Use a preprocessing driver or update the launcher before running it.
- `run_XeniumRLAnalysis.slurm` invokes `Xenium_RL_Merge_postQC_Processing_res1.5.R`, despite its broader job name.
- Some scripts restore or refresh `renv` during a job, while others assume the library is already available.
- Large objects can require 64–500 GB RAM depending on the stage. Do not copy the largest objects across too many `future` workers.
- Random seeds are set in major Seurat stages, but manual cluster renaming and dated file selection remain part of the workflow.
- The ignored `OLD/` directory contains superseded or experimental versions
  and is not part of the active workflow. This includes
  `AnchorBasedTransfer_RCTD_res1.5.R` and `SeppPreProcess.R`.

## Recommended maintenance

For future runs, the highest-impact improvements would be:

1. migrate remaining hard-coded sample definitions and paths to the existing
   configuration files;
2. replace dated input filenames with stable stage names or a manifest;
3. validate Slurm array lengths automatically;
4. eliminate absolute user-specific paths;
5. pin all R and Python dependencies; and
6. encode the stages in a workflow manager such as `targets`, Snakemake, or Nextflow.
