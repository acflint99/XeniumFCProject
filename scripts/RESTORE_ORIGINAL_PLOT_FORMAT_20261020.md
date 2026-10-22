# Restore original UMAP style and output locations

Use the native plotting stages rather than the generic snapshot exporter.
The generic `submit_refresh_current_plots.sh` workflow is no longer the
recommended redraw route: it writes a separate plot_refresh tree and uses a
generic UMAP layout. Existing snapshot files are left alone.

## References

`refresh_reference_plots_original.R` reads the saved reference objects and
redraws their original files under:

```
outputs/references/aldinger/plots/
outputs/references/sepp/plots/
outputs/references/science/plots/
```

UMAPs match the reference source scripts: 7 by 6 inches at 600 dpi for TIFF,
the original saved reduction, Seurat's default point size and rasterization,
the original label behavior, and forward broad legend order. No explicit
1200-pixel raster override or 11-by-8-inch canvas is used. The Aldinger panel
title and existing .tif summary aliases are preserved where applicable.

The full standardized and final Xenium-panel UMAPs are always redrawn.
Earlier standardized/pre-recomputed UMAPs are also redrawn if their original
plot files exist. Their saved source objects must be present; missing input
objects cause a preflight failure rather than recalculating an embedding.
The original raw published-label diagnostics are separate and are not
reclassified by this plot-only command.

Panel-reference DotPlots retain the updated dynamic dimensions, colored
headers, italic marker labels, and shared legends/scales. They replace the
original marker TIFF/PDF filenames. Broad markers still come from the shared
palette; only available genes are drawn. Duplicate shared genes retain their
first display group, consistent with the snapshot convention. The supplied
DotPlot helper applies this display convention to the native plot stages too.
Residual display marker lists allow repeated genes so the shared helper can
handle TOP2A under both glial/cycling panels. No marker-test or annotation
parameters are changed.

No RDS is written, no normalization/PCA/UMAP is run, and no cell is removed or
relabeled. The current palette supplies the approved colors. The new reference
manifest records input paths/timestamps, reduction names, output filenames
and palette hash within each original plots directory.

Upload the code archive into the HPC scripts folder, then:

```bash
cd /home/acflint/R/Projects/XeniumFCProject/scripts
tar -xzf original_plot_format_20261020.tar.gz
cd /home/acflint/R/Projects/XeniumFCProject
module purge
module load rc-base
module load R/4.4.1-foss-2023b

Rscript --max-ppsize=500000 scripts/refresh_reference_plots_original.R --list
Rscript --max-ppsize=500000 scripts/refresh_reference_plots_original.R --dry-run

# Explicit replacement of existing reference plot files, as requested.
# The verified task map is 1 Aldinger, 2 Sepp, 3 Science.
sbatch --export=ALL,REFERENCE_PLOTS_OVERWRITE=true \
  scripts/run_refresh_reference_plots_original.slurm
```

The preflight is path-only and should take under 1 minute. Expected plotting
runtime is 10-40 minutes per reference task, 3 tasks at concurrency 1, giving
30-120 minutes overall excluding queue wait. This is a low-confidence static
estimate for reading/redrawing multiple saved reference views, not a complete
HPC measurement. The requested limit is 1 hour and 64 GB per task. No package
installation or renv lockfile changes are required by this formatting update.

## Other production figures

Their native plot-only drivers retain their original UMAP dimensions,
point settings, boxed labels and layout while reading the updated shared
palette. They also write into their established stage directories. Do not
launch duplicate work while existing annotation/mapping/report jobs are still
running, and do not rerun label application or integration just to recolor plots.
The earlier stale VZ Mature-UBCs labels must be resolved in the annotation and
mapping workflow before trying to redraw those objects.

For pooled VZ/RL, the supported native plot-refresh modes are:

```bash
Rscript --max-ppsize=500000 scripts/xenium_vz_06_annotate_subclusters.R \
  --plots-only --resolution=0.8 --dry-run
Rscript --max-ppsize=500000 scripts/xenium_rl_06_annotate_subclusters.R \
  --plots-only --resolution=0.8 --dry-run

# Only after their inputs are current and the corresponding preflight passes:
sbatch --export=ALL,VZ_SUBCLUSTER_MODE=refresh-plots,VZ_SUBCLUSTER_RESOLUTION=0.8,VZ_SUBCLUSTER_OVERWRITE=true \
  scripts/run_xenium_vz_06_annotate_subclusters.slurm
sbatch --export=ALL,RL_SUBCLUSTER_MODE=refresh-plots,RL_SUBCLUSTER_RESOLUTION=0.8,RL_SUBCLUSTER_OVERWRITE=true \
  scripts/run_xenium_rl_06_annotate_subclusters.slurm
```

These dry-runs inspect paths, not RDS contents, and should take under 1 minute.
Expected merged refresh runtime is approximately VZ 30-90 minutes and RL
5-20 minutes, each one task requesting 12 hours and 256 GB. Estimates are low
confidence; the prior RL refresh took 9.33 minutes, while the VZ workload is
larger and includes marker-table/heatmap regeneration. These modes retain the
annotated RDS and its labels, but replace the stage's existing production
plot/table outputs under the explicit overwrite flag.

Native per-sample reports remain `xenium_vz_08_plot_samples.R`,
`xenium_rl_08_plot_samples.R`, and `xenium_residual_07_plot_samples.R` with their
existing Slurm launchers and protected overwrite options. Their bounds must
continue to match --list/config/samples.csv (currently 41). Current VZ/RL
per-sample report estimates are 5-30 minutes/task, concurrency 3, approximately
1.2-7 hours per 41-task array; requested limits are 2 hours/task. Residual
completed reports previously took 0.95-4 minutes, but the unresolved task-12
outlier means a complete-array elapsed estimate is still low confidence.

## Verification

All changed R files parse; shell syntax and the required R startup-stack check
pass. A small synthetic legacy-assay test confirmed original output names,
7-by-6-inch UMAP page sizes, and identical source-RDS hashes before/after.
Local graphics testing used base-PDF surrogate devices because local Cairo/X11
is unavailable; production TIFF/Cairo rendering still needs visual inspection
on Cheaha. The local production dry-run could not verify HPC RDS files because
they are not available locally. No project analysis objects were loaded.

After reviewing the code and completed figures, checkpoint this logical change
in the tracked Git checkout if ready. No commit, push, renv update or production
job submission was performed here.
