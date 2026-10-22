# Reference refresh: legacy Seurat assay compatibility

The reference failures in jobs 40984267_1, _2 and _6 used a legacy `Assay`.
For that class, `Layers(assay, search="^data($|[.])")` looks for an exact slot
name rather than applying a regular expression. The refresh therefore falsely
rejected an existing `data` slot. This behavior was reproduced locally.

The fix lists available layers with `search=NA`, then filters their names using
`grep`. It supports both legacy `Assay` and `Assay5`. Split Assay5 normalized
layers are joined only in memory, with their existing values preserved. No
normalization, package upgrade, PCA/UMAP recomputation, or RDS writing is added.
Truly missing normalized data still fails rather than falling back to counts.

The refresh also restores the established broad/reference behavior of omitting
marker groups with zero available genes, while recording all requested genes
in Marker_availability.csv. Missing required VZ marker genes remain errors.

The archive includes the corrected numeric task-map parser from the preceding
array fix. It counts 6 reference tasks even when renv/here startup messages
precede the TSV header. Colors, marker definitions, metadata and identities
are not changed by this patch.

After uploading the archive into the HPC project's scripts directory:

```bash
cd /home/acflint/R/Projects/XeniumFCProject/scripts
tar -xzf plot_palette_refresh_20261020_reference_fix.tar.gz
cd /home/acflint/R/Projects/XeniumFCProject

bash scripts/submit_refresh_current_plots.sh --dry-run \
  --groups=references --tag=palette_20261020_refs

# Run after the path-only preflight passes.
bash scripts/submit_refresh_current_plots.sh --submit \
  --groups=references --tag=palette_20261020_refs
```

Use the existing project R/4.4.1 module and renv library. Do not reinstall
packages or rerun reference-building stages for this compatibility fix.

The fresh output tag protects the partial previous outputs. New files appear
under `outputs/plot_refresh/palette_20261020_refs/references/`, with full/panel
subdirectories for each reference. Each directory contains Marker_DotPlot.pdf,
Plot_Report.pdf, TIFFs, marker availability and provenance files.

Expected runtime remains low confidence: 5-20 minutes per reference task,
6 tasks at concurrency 1 gives approximately 30-120 minutes excluding queue
wait. The requested Slurm limit is 1 hour and 64 GB per task; it is not an
elapsed-time estimate. The previous failed tasks stopped too early to measure
the complete plotting workload. Path-only preflight and completion validation
are expected under 1 minute and do not load RDS objects.

After completion:

```bash
Rscript --max-ppsize=500000 scripts/refresh_current_plots.R \
  --group=references --tag=palette_20261020_refs --validate
```

Local validation: changed R parses; legacy and modern layer tests; split-layer
value/cell-ID preservation; genuinely missing-layer detection; and a small
legacy-reference rendering test. The rendering test used base-PDF surrogate
devices because local Cairo/X11 is unavailable. Production Cairo/TIFF output
must be inspected on Cheaha. No project RDS files were loaded locally.

The pooled VZ/sample errors are a separate stale-taxonomy issue: those saved
objects still contain Mature UBCs, which the current VZ taxonomy no longer
includes. This patch intentionally does not rename those cells. Reconcile
their reviewed-label and mapping workflow before retrying those snapshot
groups. RL/residual Mature UBCs remain valid and unchanged.
