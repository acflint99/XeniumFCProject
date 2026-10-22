# Rebuild VZ and RL with all 41 samples

The seven new samples have reached the stage-1 split. Retain all 82 existing
stage-1 subset inputs. Rebuild each complete branch from stage 2; do not merge
seven new integrated objects into the old embedding or reuse old cluster IDs.
The current drivers obtain sample membership from `config/samples.csv`.

## Upload and run

Upload `xenium_regions_rebuild_20261015.tar.gz` to
`/home/acflint/R/Projects/XeniumFCProject/`, then run these commands in an Intel
interactive terminal. The archive contains source code only, including the
approved palette; it does not replace configuration, sample metadata, or RDS
inputs. Preserve the current source before extraction:

```bash
cd /home/acflint/R/Projects/XeniumFCProject
source_backup="scripts_backup_regions_$(date +%Y%m%d_%H%M%S)"
mkdir "$source_backup"
tar -tzf xenium_regions_rebuild_20261015.tar.gz | while IFS= read -r path; do
  if [[ -f "$path" ]]; then
    mkdir -p "$source_backup/$(dirname "$path")"
    cp -p "$path" "$source_backup/$path"
  fi
done
tar -xzf xenium_regions_rebuild_20261015.tar.gz
bash scripts/rebuild_xenium_regions_20261015.sh --dry-run
```

Expected preflight: approximately 1-3 minutes, low confidence from namespace
startup and filesystem checks. It checks dependencies, all 41 manifest IDs,
82 nonempty subset files, 41 consensus inputs for QC, full configuration file
checks, merge dry-runs, launcher syntax, and active regional jobs. It does not
load any large RDS objects or write analysis outputs. Existing stage-2 outputs
are expected before archiving; their presence does not invalidate a dry-run.
Do not bypass missing, unexpected, or empty inputs. The local code checks do
not establish that the current HPC data are scientifically valid.

After preflight passes:

```bash
bash scripts/rebuild_xenium_regions_20261015.sh --submit
```

This rechecks readiness, then moves previous `02_*` through `09_*` directories
for both branches into a unique directory below
`outputs/regional_rebuild_backups/`. This preserves old RDS/checkpoint bytes
without making another complete copy on the same filesystem. Those old results
are no longer available at their original paths while rebuilding. All `01_subsets`
and whole-tissue inputs stay in place. A TSV records every archived directory.

The helper queues eight jobs using existing launchers' resource, log, mail,
module, and R settings:

```text
VZ merge -> VZ integration -> VZ QC review
                         -> VZ Harmony audit
RL merge -> RL integration -> RL QC review
                         -> RL Harmony audit
```

Each child starts only after its parent exits successfully. Before its actual
work it runs the relevant path/manifest dry-run on the newly created parent
outputs; a printed `DRY-RUN FAIL` is treated as a failure even if the R driver
returns zero. Fresh output overwrite protection stays enabled. The archived
34-sample checkpoints are not used. The VZ launcher's existing automatic
checkpoint resume can operate on new-run checkpoints after a Slurm requeue;
it remains subject to provenance and overwrite checks.

The helper explicitly clears inherited resume/plots-only/overwrite choices,
pins clustering to the existing single-worker default, and forces QC **review**
mode. It does not apply either branch's previous removal list. No dimensions,
resolutions, seed, integration method, scientific thresholds, or labels change.
The only R-driver edits in this update replace misleading fixed `34` messages
with the manifest sample count.

Submission returns after queuing. Leave the files unchanged while jobs run.
The persistent `outputs/regional_rebuild_20261015.run` file points to the backup
directory; `job_ids.tsv` there records every successful submission immediately.
The helper refuses a duplicate submission. If submission stops partway through,
inspect the record and job states before recovery; do not remove the record or
rerun the helper blindly. There is no automatic cancellation or deletion.

## Expected runtime and unchanged resource requests

These are single jobs per stage, not arrays. Both branches can overlap if
resources permit. Queue waits and Slurm restarts are excluded.

| Stage | Expected active runtime | Existing requested limit / memory |
|---|---|---|
| Merge, each branch | 30-120 minutes | 12 hours / 256 GB / 1 CPU / short |
| VZ integration | 4-16 hours | 24 hours / 600 GB / 4 CPUs / medium |
| RL integration | 2-8 hours | 24 hours / 256 GB / 4 CPUs / medium |
| VZ QC review | 1-3 hours | 12 hours / 256 GB / 1 CPU / short |
| RL QC review | 30-120 minutes | 12 hours / 256 GB / 1 CPU / short |
| Harmony audit, each branch | 10-45 minutes | 2 hours / 64 GB / 1 CPU / express |

All estimates are low confidence for the enlarged cohort. Merge estimates come
from static inspection of 41 incremental RDS reads/merges, not a timed full run.
For VZ, the prior full stage-3 job reached a precluster checkpoint and lost its
node after 2h43m; another resume lost a node after 59m. Older sequential post-QC
clustering reached a 12-hour timeout. For RL, prior 34-sample post-QC
reintegration completed in 3h32m with 136 GB MaxRSS; stage 3 and the new cell
count differ. Prior VZ expanded QC took 76 minutes at 153 GB MaxRSS. RL review
and Harmony estimates are based on similar workloads and code inspection.
New sample/cell counts, completed elapsed times, and MaxRSS will refine these
ranges. Existing resource requests are retained, not increased by this helper.

Allow roughly **5.5-21 hours active time for VZ** and **3-12 hours for RL**
through QC review, with substantial additional queue time possible. The audit
can overlap QC after integration. A requested time limit is not a runtime
prediction. Do not start a second concurrent regional rebuild to accelerate it.

## Completion and review gate

```bash
squeue -u acflint
cat outputs/regional_rebuild_20261015.run
```

Read `job_ids.tsv` in that backup directory. For the eight recorded IDs:

```bash
sacct -j JOBID1,JOBID2,JOBID3,JOBID4,JOBID5,JOBID6,JOBID7,JOBID8 \
  --format=JobID,JobName%28,Elapsed,State,MaxRSS,ReqMem,AllocCPUS,ExitCode \
  --units=G
```

Expect completed jobs with `0:0`. Missing success dependencies require inspection
of the failed parent's logs; do not resubmit the entire workflow. The new merge
manifests under each branch's `02_merged/tables/` must contain all 41 sample IDs
and positive per-sample cell counts. The integration drivers validate manifest
membership and cell totals against the loaded object.

Review these outputs in **both** branches:

- `03_integrated/plots/`: corrected UMAPs at resolutions 0.5, 0.8 and 1,
  consensus-label and sample-qualified Unknown summaries.
- `03_integrated/tables/harmony_audit/`: mixing, biological-label retention,
  developmental-stage neighborhood change, and pre/post neighbor overlap.
- `04_qc/plots/`: raw-subcluster QC violins, per-sample QC and sample-composition
  heatmaps, and Unknown QC evidence.
- `04_qc/tables/`: raw-subcluster and per-sample QC summaries, top-10 markers,
  Unknown evidence, and the new QC review manifest.

**Stop here to review new cluster IDs.** Then explicitly choose no removal or
identify newly reviewed removal clusters. Only afterward should stages 4b and
5 proceed, followed by regenerated sample-aware marker evidence, descriptive
markers, fresh subcluster-label review, mapping, and spatial/count reports.
Do not reuse old QC-removal cluster numbers or old stage-6 mappings. The prior
stage-5 through stage-9 artifacts are archived because they belong to the old
cohort. Existing combined VZ/RL/residual and spatial outputs elsewhere are also
stale; this helper does not regenerate them before both new branches are mapped.

## Rollback and source checkpoint

Retain the backup until the new analysis is reviewed. To roll back, first stop
or wait for all recorded regional jobs, preserve any new partial results under
another unique archive, and restore complete old stage directories together
using `archive_manifest.tsv`. Never mix an old checkpoint, QC decision, or
annotation object with the new merge. Restore source from the source backup
only if required. No automatic rollback is performed.

Before a source Git checkpoint on Cheaha, inspect `git status --short` and
`git diff`, parse the changed R scripts (under 10 seconds each), run `bash -n`
on the helper and included launchers (under 5 seconds each), and confirm the
HPC preflight passed. Commit only reviewed source changes; no generated RDS,
reports, backup directories, logs, or run records belong in the commit. Package
versions were not changed, so no renv snapshot is needed for this update.
