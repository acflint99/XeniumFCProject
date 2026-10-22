# Preserve FOVs across the two-level spatial merge

Job 40985456 failed after 06:18:51 with exit 1:0 at the final image-count
validation. Its batch MaxRSS was 233350536K (about 222.5 GiB), with 600G
requested. All seven batches and their final merge completed; validation
stopped the job before the final saveRDS and manifest write. Existing output
files therefore cannot be treated as products of this failed run.

Seurat's duplicate-image suffixes can collide with existing image names when
merging already-merged batches. A synthetic 41-sample/seven-batch fixture with
one FOV per sample reproduced 36 surviving FOVs. The fixed path retains 41.
This fixture establishes a reproducible mechanism; the failed HPC log does
not report the actual surviving names or count.

The driver now assigns sample-qualified FOV names before any batching, such
as FB124_X_G__fov. Only names of image-list entries change: cell IDs, sample
IDs, FOV contents, coordinates, counts, metadata, and analysis parameters are
preserved. The merge helper refuses overlapping image names before merging
and checks exact image-name preservation after each merge. Existing cell,
image-count, and regional-assignment checks remain. Two additive manifest
columns record input and merged image names. Output filenames are unchanged.

Files: xenium_vz_rl_03_spatial_merge.R, R/vz_rl_batch_merge.R,
test_vz_rl_spatial_merge.R. No launcher, package, or configuration change.

## HPC checks before a replacement job

Copy those files to their matching scripts/ paths, or extract the code-only
spatial_fov_merge_fix_20261020.tar.gz archive in the HPC scripts directory.
From /home/acflint/R/Projects/XeniumFCProject, using the established R/4.4.1
module and renv environment:

```bash
Rscript --max-ppsize=500000 scripts/test_vz_rl_spatial_merge.R
Rscript --max-ppsize=500000 scripts/xenium_vz_rl_03_spatial_merge.R --dry-run
ls -lh outputs/xenium/vz_rl/spatial/03_merged/rds/Xenium_VZRL_spatial_merged.rds \
  outputs/xenium/vz_rl/spatial/03_merged/tables/Xenium_VZRL_spatial_merged_manifest.csv
squeue -u acflint -o '%.24i %.36j %.10T %.12M %.180E %R'
```

The synthetic test should take 15–90 seconds (low-confidence cross-host
estimate); it uses tiny objects and writes no analysis outputs. The driver
dry-run is path/CSV/provenance-only, with no RDS loading, estimated under one
minute. Missing files in ls are information, not a reason to create them.

Review/preserve any existing merge outputs before explicit overwrite. Do not
remove the validation guard or rerun the whole annotation/plot workflow.
There is no persisted batch checkpoint in this driver: a new merge must reread
the inputs. Expected replacement merge runtime is roughly 6–9 hours excluding
queue wait, low confidence based on the 6:18:51 failed run reaching the final
validation plus time to serialize the large successful result. The existing
request remains 24 hours, 600G, one CPU, medium.

Job 40985459 depends on the failed merge and needs a replacement dependency
chain. The VZ palette installer 40985724 also included the failed merge and
spatial processing in its dependency snapshot. Inspect its full current
dependency list and states before reconnecting it: retain all other unfinished
prerequisites. Do not simply release either job unconditionally. Downstream
spatial integration must use the validated replacement merge. The independent
nonspatial branch and successful sample redraws do not require rerunning.

## Validation and rollback

All three changed/added R files parse. The synthetic fixture checks batched
versus pairwise merging for all 41 samples, including FOV names/memberships,
coordinates, cell IDs/order, expression counts, and sample/label metadata.
It also rejects overlapping image names. Local checks used the available
R 4.5.1 / SeuratObject 5.2.0, not the HPC R 4.4.1/renv environment; run the
same tiny fixture on HPC before submission. The local production preflight
cannot certify HPC input readiness because production paths are unavailable.

No production outputs were overwritten and no HPC jobs were submitted here.
Retain the previous scripts and corresponding prior valid merge/integration
outputs together for rollback; a failed run does not provide a new valid
merge. Inspect the new manifest and collect sacct JobID, Elapsed, State,
MaxRSS, ReqMem, AllocCPUS, and ExitCode after completion. In the tracked HPC
checkout, review git status/diff and validation, then commit this logical fix
when ready. No package environment change calls for renv::snapshot().
