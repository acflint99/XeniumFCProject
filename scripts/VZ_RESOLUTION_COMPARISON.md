# Efficient VZ comparison: resolutions 0.5, 0.6, 0.7 and 0.8

Run `xenium_vz_05d_compare_resolutions.R` on the existing post-QC VZ object.
This is a protected full-data comparison after Stage 5, not a new subsampling
run. No production labels, QC decisions, cell/sample IDs, integration, graph
parameters, marker thresholds or existing outputs are changed.

## Speed and memory improvements

- Load the approximately 77.6-GB input file once per job.
- Preserve the saved 0.5 and 0.8 assignments exactly. Compute only 0.6/0.7 on
  the existing `Xenium_snn`, using the saved Louvain parameters and original seed.
  No PCA, Harmony, neighbor construction or UMAP fitting is performed.
- Save a partition-only checkpoint immediately after clustering. Explicit resume
  can skip both new clusterings after a downstream failure. It still loads the
  production object for expression summaries and marker tests.
- Release graphs after clustering. After joining Xenium layers, slim only the
  in-memory working object to Xenium counts/data, all metadata/cells and the saved
  UMAP. Drop unused reductions, scaled data and unused assays. Synthetic checks
  require identical normalized values, counts, metadata, UMAP and marker results.
- Produce PDF marker plots once by default, using R's standard PDF device.
  Optional 600-dpi TIFFs require `--tiff`; no biological calculation depends on
  this display-only option. The PDF device avoids an additional external Cairo
  dependency for the default review. Existing production callers of the shared
  plotting helper retain their original TIFF and Cairo-PDF behavior.
- Keep marker testing unchanged: cell-level Wilcoxon, positive markers,
  `min.pct=0.25`, `logfc.threshold=0.25`, up to 1000 cells per identity, configured
  marker seed. Do not pre-subset the whole object or reduce this cap for speed:
  that could change fold changes, background comparisons and marker selection.
- Record whether optional `presto` is already installed. Seurat may use its
  existing accelerated backend; this script never installs or changes packages.
- Keep one worker. Earlier VZ parallel object duplication caused OOM failures.
  There is no claimed measured speedup until this new comparison is timed.

## Files to synchronize to HPC

Copy these into matching `scripts/` locations:

1. `xenium_vz_05d_compare_resolutions.R`
2. `run_xenium_vz_05d_compare_resolutions.slurm`
3. `R/vz_resolution_efficiency.R`
4. Updated `R/subcluster_review_markers.R` (adds optional TIFF omission and PDF
   device selection; existing callers retain their defaults).

Also keep the shared helpers from the earlier RL comparison/stability review
current: `R/regional_resolution_comparison.R`, `R/regional_cluster_stability.R`
(including the explicit ggraph edge guide), `R/dotplot_helpers.R`, `R/config.R`
and `color_palette.R`.

## Run on Cheaha

Preflight should take under one minute per command. `--list` reads the authoritative
manifest (currently 41 samples). `--dry-run` checks configuration, paths, output
protection and checkpoint presence only; it never loads the large input.
`--check-packages` also renders a tiny synthetic clustree PDF and removes it.
Actual execution verifies cell/sample coverage, ordered IDs, input provenance,
graph and UMAP values, and saved clustering settings before running the comparison.

```bash
cd /home/acflint/R/Projects/XeniumFCProject
module purge
module load rc-base
module load R/4.4.1-foss-2023b

if Rscript --max-ppsize=500000 scripts/xenium_vz_05d_compare_resolutions.R --check-packages &&
   Rscript --max-ppsize=500000 scripts/xenium_vz_05d_compare_resolutions.R --list &&
   Rscript --max-ppsize=500000 scripts/xenium_vz_05d_compare_resolutions.R --dry-run
then
  sbatch scripts/run_xenium_vz_05d_compare_resolutions.slurm
fi
```

Estimated active runtime: **3–10 hours**, excluding queue wait, for one pooled
job. Confidence is low. Basis: the 1,447,477-cell VZ pilot needed 158–173 minutes
per repeat for an 80% graph rebuild and four partitions. This job computes two
full-data partitions without rebuilding the graph, then four marker/plot reviews.
The current VZ Stage-6 review took 46m19s; the earlier review measured 152.09 GiB
MaxRSS, but stages and evidence reuse differ. The new comparison has not been
benchmarked, and the improvements do not justify a guaranteed speedup factor.

Requested resources: **12 hours, one CPU, 256 GB, short partition**, matching
the established VZ Stage-6 review launcher. This allocation is for the sequential
diagnostic, not the 600-GB full reintegration. Initial loading/clustering remains
the peak-memory uncertainty; `timings.csv` records individual load, clustering,
layer-join, slimming, marker and focused-plot times to refine future requests.
The pipeline's previous parallel worker failures are a reason to retain sequential
execution rather than add CPUs without a measured benefit.

Collect actual resource history after completion (normally under one minute):

```bash
sacct -j <JOB_ID> --format=JobID,Elapsed,State,MaxRSS,ReqMem,AllocCPUS,ExitCode
```

## Retry without recalculating completed partitions

After both new partitions are computed, a label-only checkpoint is published at:

```text
outputs/xenium/vz/05_post_qc/diagnostics/resolution_comparison_05_06_07_08/.review1_partitions.rds
```

It contains assignments, IDs and provenance, not expression matrices or a Seurat
analysis object. It is uncompressed for fast I/O, so it can still occupy hundreds
of MB. It is copied into a successful report directory as `partitions_checkpoint.rds`.
Do not commit diagnostic outputs to Git.

If marker testing or plotting fails after checkpoint creation, retain the files
and resume explicitly:

```bash
sbatch --export=ALL,VZ_RESOLUTION_COMPARE_RESUME=true \
  scripts/run_xenium_vz_05d_compare_resolutions.slurm
```

Resume checks input path/size/time, Stage-5 manifest, configuration and sample
manifest checksums, driver/helper checksums, Seurat/SeuratObject versions, saved
clustering parameters, exact ordered cell/sample IDs and original 0.5/0.8 anchors.
A mismatch is an error. A script/input change may therefore require a fresh run ID.
Resume saves the two clustering steps but repeats input loading, marker tests and
plots; expected runtime is provisionally **1–5 hours**, low confidence from the
four review tasks. If a failure happens before the checkpoint is complete, use
a new run ID rather than request resume from a missing checkpoint.

Existing outputs and checkpoints are protected. A fresh run uses
`--export=ALL,VZ_RESOLUTION_COMPARE_RUN_ID=review2`. Explicit reviewed replacement
uses `VZ_RESOLUTION_COMPARE_OVERWRITE=true`; previous diagnostics/checkpoints are
archived before replacement. A killed job may leave a `.lock_<run-id>` directory;
verify that no job still owns it before manually removing an abandoned lock.
R-handled errors release the lock automatically and retain partial staged outputs.

Optional TIFF export adds `VZ_RESOLUTION_COMPARE_TIFF=true` and extra rendering
time/I/O. Standard PDF review uses all cells; there is no display subsampling.

## Outputs and interpretation

Default complete directory:

```text
/home/acflint/R/Projects/XeniumFCProject/outputs/xenium/vz/05_post_qc/diagnostics/resolution_comparison_05_06_07_08/review1/
```

- `clustree.pdf`: four-resolution full-data membership relationships.
- `res<resolution>_UMAP.pdf`: all cells on the identical saved coordinates.
- `res<resolution>_focused_markers.pdf`: established VZ/PC/GABA, glial, RL VZ,
  UBC and other lineage panels; no identities are automatically assigned.
- `res<resolution>_top50_markers.csv`, `marker_status.csv` and
  `marker_availability.csv`: descriptive markers and explicit rare-cluster/panel
  limitations. Clusters with fewer than three cells are retained but not tested.
- `production0.8_transitions.csv`: exact transitions from each current 0.8 cluster
  into 0.5/0.6/0.7, with source-retention and destination-composition fractions.
- `focus_population_transitions.csv`: the same evidence for current clusters
  8/12/17/19/20/21/24, chosen from the weak/variable pilot populations.
- `all_pair_transitions.csv`, `assignments.csv.gz`, `cluster_summary.csv`,
  `sample_composition.csv`, `input_label_composition.csv`: complete cell-ID
  assignments and cross-resolution/sample/input-label summaries.
- `partitions_checkpoint.rds`, `timings.csv`, `provenance.csv`, `sessionInfo.txt`,
  `METHODS.txt`, `outputs.csv`: checkpoint, measured steps, provenance and checksums.

Treat current working labels as hypotheses, not truth for choosing resolution.
Cell-level marker P-values do not establish independent biological-replicate
significance. DotPlot colours are gene-wise scaled within each dataset, not absolute
cross-resolution expression changes. This review does not replace replicate-aware
Stage-6 marker evidence. No downstream production output becomes stale from the
diagnostic alone. Approving a different production partition would require a new
Stage-6 label review and regeneration of affected Stage-7 mappings and dependent
combined/spatial/count/export outputs.

## Lightweight validation and checkpoint

Syntax checks take under 10 seconds per changed R file. The synthetic efficiency
test should take 5–30 seconds on 120 artificial cells and verifies unchanged
counts, normalized data, metadata, UMAP, and exactly equal marker results after
slimming, plus checkpoint rejection and PDF-only rendering. Shell syntax takes
under 5 seconds; the launcher audit under one minute. These do not validate the
scientific interpretation or establish full-size performance.

```bash
Rscript --max-ppsize=500000 -e 'for (f in c("scripts/xenium_vz_05d_compare_resolutions.R", "scripts/R/vz_resolution_efficiency.R", "scripts/R/subcluster_review_markers.R", "scripts/test_vz_resolution_efficiency.R")) parse(f)'
Rscript --max-ppsize=500000 scripts/test_vz_resolution_efficiency.R
bash -n scripts/run_xenium_vz_05d_compare_resolutions.slurm
cd scripts
python3 update_r_launchers_stack_limit.py --check
```

After reviewing changes in the actual Git repository with `git status` and
`git diff` and running validation, consider a commit for this logical unit when
ready. No commit/push is performed here. No package environment or `renv.lock`
change is made; a snapshot is not needed solely for these code edits.
