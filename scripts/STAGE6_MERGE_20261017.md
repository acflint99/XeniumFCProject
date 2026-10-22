# Combined VZ/RL stage-6 review — 2026-10-17

The existing stage-6 review drivers now include the replicate-aware marker prerequisite at production resolution 0.8. One existing launcher submission builds missing evidence and then generates cluster review outputs. Both driver/launcher defaults are resolution 0.8. The standalone 06a scripts/launchers remain optional evidence-only entry points using the same scientific builder. Alternate-resolution review, explicit label application, and plot-refresh behavior are preserved.

Evidence handling: absent files are built; all five nonempty files with matching input/cluster provenance are reused. Partial, empty, or stale evidence is blocked unless --overwrite was explicitly supplied after review. Existing review-output overwrite protection runs before evidence generation. Builder failure prevents review. Evidence is validated again after generation. The evidence builder releases its large object and garbage collection runs before review loads its object. Two sequential object reads remain; this avoids changing the builder's scientific object interface or holding two phase objects concurrently.

Expected fresh combined runtime: VZ 1–10 hours; RL 0.75–7 hours, low confidence from static inspection and earlier separate-stage estimates. Both launchers retain 256 GB, one CPU and a 12-hour limit. Queue wait excluded. Actual elapsed/peak-memory measurements remain necessary; these requests are not guarantees.

After syncing the changed scripts/shared helper to the HPC, run from the project root:

```bash
Rscript --max-ppsize=500000 scripts/xenium_vz_06_annotate_subclusters.R \
  --review-only --resolution=0.8 --dry-run
```

The preflight checks paths and small provenance CSVs only; expected under one minute. When it passes and existing outputs have been checked:

```bash
sbatch scripts/run_xenium_vz_06_annotate_subclusters.slurm
```

RL uses the corresponding rl driver and launcher. Missing evidence is shown as planned output rather than a missing input. Explicit overwrite remains necessary for replacing stale/partial evidence or existing review outputs. No overwrite was enabled by default. Inspect review plots/tables before applying labels.

Scientific settings, sample/cell IDs, ranking definitions, sample-level inference, marker caps/thresholds, seeds, eight review output names and five evidence output names remain unchanged. Rebuilding evidence intentionally changes its provenance signature, so downstream reviewed label application must use the matching regenerated review. The post-QC object remains unchanged for rollback.

Validation: all three edited R files parsed; both launchers passed bash -n; pointer-stack audit passed. Lightweight evidence-planning checks covered absent, matching, partial, stale and explicit overwrite cases. Both actual drivers passed path-only integration preflights using a present-input/absent-evidence fixture; partial evidence correctly blocked preflight without writes. No fixture was loaded as an RDS. Actual-project preflights still lack local HPC RDS inputs, so real readiness and scientific output remain unverified. No analysis, package installation, configuration/sample-manifest edit or job submission was performed.

Changed files: R/xenium_replicate_marker_evidence.R (planning helper only; scientific builder unchanged), both xenium_{vz,rl}_06_annotate_subclusters.R drivers, both matching stage-6 launchers, README.md, and this note. Sync the shared R helper with the drivers. Earlier resource notes describing a mandatory separate 06a job are superseded by this merged order. Review git status/diff and checkpoint this logical unit in the actual repository when ready; no commit or push was made, and this local scripts folder is not a Git checkout. No package changes require renv::snapshot().
