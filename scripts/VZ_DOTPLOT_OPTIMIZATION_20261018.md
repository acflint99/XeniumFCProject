# VZ production DotPlot optimization — 2026-10-18

Scope: main, lineage, and Other production VZ DotPlots, plus per-sample standalone DotPlots and embedded sample-report DotPlots. Source review: VZ Subclusters 10-18-26, approved cluster-label and QC tables, prior raw top-50 marker tables and replicate-aware marker evidence. No annotation, raw cluster assignment, sample/cell ID, QC threshold, statistical marker test, normalization, scaling definition or routing decision was changed. No existing attached plot, RDS, or output manifest was overwritten locally.

## Marker panels

| Group | Markers |
|---|---|
| VZ | PRDM13, ASCL1, HES6 |
| Neuroblasts | DCX, STMN2, SOX11, MAP2 |
| Purkinje | FOXP2, ITPR1, CALB1, FOXP1, PTGER3, EBF1, FOXP4, EBF2, NDNF, TRPC3 |
| GABA | PAX2, TFAP2B, GAD1, LRRTM1, GRIK3, NXPH2, SLC6A5, SOX14, DMBX1 |
| Glia | SOX2, NES, PTPRZ1, SLC1A3, TNC, FAM107A, AQP4, GJA1, FOXJ1 |
| OPCs | EGFR, OLIG1, OLIG2, PDGFRA, SOX10 |
| RL VZ | CALCB, OTX2, RSPO3 |
| UBCs | EOMES, ZIC2, CBLN1, RELN, TRH |
| eCN | MEIS2, SLC17A6, C1QL3 |

All 51 unique genes occur in the supplied 5,001-gene replicate evidence table. No MKI67, TOP2A, PCNA, TYMS, RRM2, or CCND2 columns remain in production panels. Existing broad raw-review marker lists retain their diagnostic cell-cycle panels and signatures. Marker selection above is for display and does not recompute or change biological labels. Low expression of the new Neuroblasts panel remains informative; adding these markers does not strengthen the assignment by itself.

Purkinje changes improve the display of the assigned four states: core FOXP2/ITPR1/CALB1, FOXP1/PTGER3, EBF1/FOXP4, EBF2/NDNF, and TRPC3. GABA markers cover the assigned progenitor, MLI, Golgi and iCN groups. Glial markers contrast progenitor/radial signals with TNC/FAM107A/SLC1A3 and AQP4/GJA1/FOXJ1. Neither low neuroblast complexity nor tentative fine-subtype interpretations were relabeled by these presentation changes.

## Ordering and colors

The current VZ order contains exactly the 20 observed approved labels. Neuroblasts follows VZPs and precedes the four Purkinje labels. Unused Cycling VZPs, generic Cycling Cells and Microglia are omitted from the current VZ order and production panels. Assigned Cycling Glia stays present. Legacy palette colors and other branches' labels remain available in the shared lookup; no cells were deleted. The master combined order places Neuroblasts by VZPs without dropping RL/residual identities.

The production lineage DotPlot includes Neuroblasts immediately after VZPs. This is a presentation-only partition: biological lineage/Other membership, cell filtering and VZ/RL routing are unchanged. The Other DotPlot has four rows: RL VZ, Immature UBCs, Mature UBCs and eCN.

One combined UBCs header uses broad cluster_colors['UBC'] (#FF7F00), per the user's correction; no separate immature/mature headers are introduced. eCN uses vz_palette['eCN'] (#FF4500). RL VZ uses its existing subcluster color (#E00065), and Neuroblasts uses its existing color (#FF66B2). Remaining headers retain their broad group colors. VZ-specific colors are passed explicitly to the shared renderer so RL/other callers' default color behavior is unchanged.

Marker groups are selected according to identities actually present in each merged/sample plot. Empty group headers and their marker columns are omitted. Sample-report embedded DotPlots use the same colored rendered object as standalone exports.

## Validation

All four changed R files parsed successfully (<1 second combined). Palette invariants confirmed unique genes, all selected genes in the supplied panel, exact observed-label coverage, Neuroblasts in the requested position, cycling-marker exclusion and exact UBC/eCN color equality. Missing-group selection was checked for a Neuroblasts-only and a mature-UBC-only example.

Synthetic layout checks rendered main, lineage and Other panels to PNG and were visually inspected. An actual Seurat DotPlot test on 100 synthetic cells passed through the same rendering helper and an embedded patchwork report. These were explicitly illustrative layout fixtures, not scientific result plots. No new PDF artifact was authored and no real RDS was loaded. Project startup-stack audit passed. The two supported real-project path-only dry-runs executed but reported missing local HPC input RDS files, so full-data regeneration remains an HPC step.

Raw review panels/signatures and the approved 0–24 label map were checked to remain unchanged. These display-only edits do not require rerunning integration, marker testing, label application or sample mapping. Regenerate production plots and any desired existing sample report pages to show the new panel.

## HPC refresh

Sync these files together:

- color_palette.R
- R/dotplot_helpers.R (new optional strip_colors argument)
- xenium_vz_06_annotate_subclusters.R
- xenium_vz_plot_sample.R

The palette fails early if an old renderer without strip_colors was synced. No new dependencies or renv changes are required.

From /home/acflint/R/Projects/XeniumFCProject, run the path-only preflight (<1 minute; no large RDS load):

```bash
Rscript --max-ppsize=500000 scripts/xenium_vz_06_annotate_subclusters.R \
  --plots-only --resolution=0.8 --dry-run
```

After PASS, refresh existing production plots/QC summary using the protected plotting mode:

```bash
VZ_SUBCLUSTER_MODE=refresh-plots \
VZ_SUBCLUSTER_RESOLUTION=0.8 \
VZ_SUBCLUSTER_OVERWRITE=true \
  sbatch --export=ALL scripts/run_xenium_vz_06_annotate_subclusters.slurm
```

This explicitly replaces the existing stage-6 plotting outputs and QC summary, preserving the annotated RDS and label manifest. Expected runtime 10–90 minutes, low confidence from static inspection and the earlier partial apply/refresh workflow; queue wait excluded. Requested limit remains 12 hours, one CPU and 256 GB. Plot refresh also regenerates existing UMAP/QC/heatmap products through the established mode; it does not rerun integration or apply labels.

Per-sample report refresh, when desired, uses the existing 41-task stage-8 launcher/submission helper after sample mapping is complete. Expected 5–30 minutes/task and approximately 70–420 minutes for 41 tasks at concurrency three (low confidence; plot sizes and I/O dependent), versus the two-hour limit and 64 GB/task. An existing sample report refresh explicitly replaces that report's established outputs. First inspect --list and the supported path-only dry-run (<1 minute each); no such report overwrite or job submission occurred locally.

For rollback, retain the existing production plots until ready to replace them, and restore the prior display-panel/renderer files if needed. Analysis objects remain unchanged. In the actual repository review git status/diff and checkpoint this logical plotting update when ready. This local folder is not a Git checkout; no commit, push, package install, or renv snapshot was performed.
