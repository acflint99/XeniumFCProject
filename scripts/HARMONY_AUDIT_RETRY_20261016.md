# Retry only the VZ and RL Harmony audits

VZ audit 40951436 failed after 3m34s with a segmentation fault, `sacct` exit
11:0, 61.20 GB MaxRSS, and a 64-GB request. The stderr did not include an R
traceback, and Slurm did not classify this as OUT_OF_MEMORY. Memory pressure
is a working hypothesis, not a confirmed root cause. The RL failure was
reported without its job ID or accounting, so its exact failure is unverified.

## Changes

- Both audit launchers request 256 GB rather than 64 GB, retaining one CPU,
  the Intel constraint, express partition, two-hour time limit, module setup,
  log names, project environment, and protection-stack option. The larger
  request provides initial-load headroom: prior regional QC used 153 GB and
  prior RL reintegration used 136 GB, although they are different workloads.
  [Current Cheaha hardware documentation](https://docs.rc.uab.edu/cheaha/hardware/)
  lists express with a two-hour limit and Intel nodes with 754 GB. This does
  not guarantee that a 256-GB job starts immediately or that memory caused
  the prior failure.
- The shared audit logs loading, extraction, each neighbor search, and writing.
  After extracting the sampled embeddings and required metadata, it releases
  the full Seurat object before the searches. It also releases embeddings
  after each corresponding search. Initial `readRDS()` still requires enough
  RAM for the full input; this patch does not make the input a streaming file.
- Numeric/finite matrix, dimension-count, and neighbor-count checks fail
  informatively before native nearest-neighbor work on invalid inputs.
- The audit age parser accepts numeric ages and stored `PCW12` / `PCW17`
  labels. This fixes an additional input-format problem; it does not explain
  the observed segmentation fault or change ages stored in the objects.

The 100,000-cell cap, seed, sampling algorithm, 30 dimensions, 30 neighbors,
and metric definitions are unchanged. No input RDS, clustering, integration,
QC decision, subcluster annotation, or regional mapping is rewritten.
Only the failed audit tables need regeneration.

## Upload and retry

Upload `xenium_harmony_audit_fix_20261016.tar.gz` to the project root. The
archive contains the shared audit helper, two audit launchers, retry helper,
and this guide; it contains no data or configuration. Preserve current source:

```bash
cd /home/acflint/R/Projects/XeniumFCProject
source_backup="scripts_backup_harmony_audit_$(date +%Y%m%d_%H%M%S)"
mkdir "$source_backup"
tar -tzf xenium_harmony_audit_fix_20261016.tar.gz | while IFS= read -r path; do
  if [[ -f "$path" ]]; then
    mkdir -p "$source_backup/$(dirname "$path")"
    cp -p "$path" "$source_backup/$path"
  fi
done
tar -xzf xenium_harmony_audit_fix_20261016.tar.gz
bash scripts/retry_xenium_harmony_audits_20261016.sh --dry-run
```

Preflight is expected to take under one minute, based on two path-only driver
checks and scheduler accounting queries. It does not load RDS objects. It
prints the original rebuild job outcomes and requires both recorded
integration jobs to be COMPLETED with 0:0. It checks that no audit is already
active, verifies launcher syntax and the 256-GB request, and checks both input
paths. A missing original run record, accounting failure, unsuccessful
integration, or missing input blocks submission. Diagnose that upstream job
instead of rerunning the rebuild helper. QC status is printed for separate
review; an audit retry does not repair a failed QC job.

After preflight passes:

```bash
bash scripts/retry_xenium_harmony_audits_20261016.sh --submit
```

The helper rechecks readiness, archives existing/partial audit directories
under `outputs/report_backups/harmony_audit_retry_20261016_*/`, and submits
only the two audit jobs with overwrite protection enabled. It records each
submitted ID immediately in `job_ids.tsv`. If submission stops partway
through, inspect these IDs and active jobs before retrying. Do not rerun
`rebuild_xenium_regions_20261015.sh` to repair this audit failure.

Expected audit runtime: **10-60 minutes per branch**, low confidence based on
loading a large regional RDS and two searches capped at 100,000 cells in 30
dimensions. The failed 3m34s is not a successful-runtime measurement. Both
jobs can run concurrently if memory is available: approximately 10-60 minutes
combined active elapsed time if simultaneous, or 20-120 minutes if serial.
Queue wait is excluded. Each job still requests two hours and 256 GB.
Peak RSS and new stage messages are needed to refine the request or locate
a persistent native-code failure. Do not run the full audits interactively.

## Verify the retry

Use the new IDs printed by the helper:

```bash
sacct -j VZ_JOB_ID,RL_JOB_ID \
  --format=JobID,JobName%28,Elapsed,State,MaxRSS,ReqMem,AllocCPUS,ExitCode \
  --units=G
```

Expected outputs in each branch's `03_integrated/tables/harmony_audit/`:
`Xenium_<VZ|RL>_Harmony_Audit_Summary.csv`,
`Xenium_<VZ|RL>_Harmony_Audit_Paired_Metrics.csv`, and
`Xenium_<VZ|RL>_Harmony_Audit_Provenance.csv`.
Require both jobs to finish COMPLETED/0:0, all three nonempty tables per
branch, and provenance naming the new integrated input. If another crash
occurs, retain both stdout and stderr plus sacct; the new messages identify
the last phase entered. Do not assume another memory increase is sufficient.

Local checks: changed R syntax parsed; launcher/helper shell syntax and
protection-stack rule checked. A 240-cell synthetic Seurat fixture produced
identical numeric-age audit summaries and paired metrics before and after
the memory change, with the same seed and sampled cells. PCW-prefixed ages
produced those same metrics after parsing. Overwrite and non-finite input
guards were exercised. Path-only dry-runs and simulated Slurm retry gates
were checked; no full HPC audit was run by Codex.

For a source checkpoint, review `git status --short` and `git diff`, then
commit the audit fix after validation if ready. Do not add generated audit
tables, logs, or backups. No packages or package versions changed, and this
fix does not require a renv snapshot. Preserve the source and audit-output
backups for rollback; never restore failed partial tables as successful ones.
