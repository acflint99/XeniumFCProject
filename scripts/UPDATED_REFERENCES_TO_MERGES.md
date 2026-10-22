# Updated Science, Sepp and Aldinger labels through both merges

Run the four phases below from the project root on Cheaha. Each phase submits
its complete dependency chain and stops at the next scientific review gate.
No Xenium import/QC or resolution-5 clustering is repeated. The existing
resolution-5 RPCA settings, full-reference sampling, 0.4 score / 0.05 margin
abstention, weighted 2-of-3 consensus and regional resolution 0.8 are preserved.
User-approved routing: **eCN -> RL; Vascular -> Other/nonregional**. This changes
regional cell membership; all previous regional cluster maps are stale.

The updated helper derives **41 samples / 123 reference-score tasks** from the
active manifest/driver task maps. Do not use the older helpers that require 34.
All R calls use `--max-ppsize=500000`.

## Upload and prepare once

Sync the updated scripts, the new submitter/launcher, and the modified configuration:

- `submit_updated_references_to_merges.sh`
- `run_aldinger_02_standardize_cell_types.slurm`
- `science_01_standardize_cell_types.R`
- `fetal cerebellum_science_ACF.R`, `science_03_map_reviewed_res2.R`
- `R/science_reviewed_mapping.R`, `R/science_ecn_marker_evidence.R`, `R/dotplot_helpers.R`
- `reviews/Science_res2_reviewed_labels_20261022.csv`
- `run_science_map_reviewed_res2.slurm`
- `R/corrected_broad_labels.R`
- `config/config.yml` and `config/reference_broad_taxonomy.csv`

`color_palette.R` must include the already reviewed eCN and Vascular colors/order.
Keep the user's updated Aldinger/Sepp standardization scripts. The Science source
is the newly reviewed retained object, not the old resolution-1.2 file:
`outputs/references/science/reprocessed/res2/reviewed_retained/Science_res2_reviewed_retained_reprocessed.rds`.
If this RDS is missing, the annotation chain automatically queues the approved
resolution-2 mapping/removal/reprocessing job before Science standardization.
Its preflight requires the existing resolution-2 RDS and top-20 marker CSV in
`outputs/references/science/reprocessed/res2/`. The mapping validates all 45
reviewed source clusters and their marker evidence, removes Brainstem/RBC/Other,
and reprocesses the retained cells. Existing partial reviewed outputs are replaced
using the protected overwrite flag. The raw resolution-2 source is preserved.

The configuration loader prefers project-root `config/` over `scripts/config/`.
Update the **active** files reported by preflight; avoid leaving a stale root
configuration shadowing the uploaded one. Review the diff when copying config;
do not replace unrelated sample manifests, QC decisions or scientific settings.
Required reference sources are the existing Aldinger updated Seurat RDS, original
Sepp SCE RDS, reviewed Science RDS (or its raw resolution-2 RDS and top-20 marker CSV),
and `inputs/xenium_5k_genes.rds`.

Use an existing Intel compute allocation with the established project environment:

```bash
cd /home/acflint/R/Projects/XeniumFCProject
module purge
module load rc-base
module load GCC/13.2.0
module load R/4.4.1-foss-2023b
export REFERENCE_REFRESH_RUN_ID=updated_reference_labels_20261022
```

Do not run a competing reference/annotation/regional rebuild against these same
output paths. Review or snapshot existing derived outputs before `--submit`:
submission explicitly enables their protected overwrite switches. Sepp and
Aldinger reference scripts have legacy direct writes; those reference RDS/plots
will also be replaced. No analysis outputs are modified by `--dry-run`.
Dry-runs inspect paths, dependencies, small CSVs and task maps, not large RDS files.
Allow roughly **1–3 minutes per phase preflight**, low confidence from repeated R
startup and filesystem latency. Missing files, printed `DRY-RUN FAIL`, or failed
dependencies stop submission. Future-output checks run in their downstream jobs;
the job's actual driver still validates its input even where a legacy launcher
does not provide a separate dry-run.

## Phase 1: prepare references, annotate, and generate review evidence

```bash
bash scripts/submit_updated_references_to_merges.sh annotation --dry-run
bash scripts/submit_updated_references_to_merges.sh annotation --submit
```

This queues three standardization -> panel-subset chains. Aldinger starts at
`aldinger_02_standardize_cell_types.R`; the unchanged
`Aldinger_seurat_updated.rds` must exist. It does not rerun
`aldinger_01_update_seurat_object.R`. Sepp starts at `sepp_01_standardize_cell_types.R`;
Science starts at `science_01_standardize_cell_types.R` when its reviewed retained
RDS exists; otherwise the chain first builds that object using
`run_science_map_reviewed_res2.slurm`. Aldinger/Sepp can run while Science is
reprocessed. After cluster removal,
Aldinger/Sepp recompute normalization, variable features, scaling, PCA, neighbors
and UMAP while preserving labels; their panel scripts repeat that processing
after gene restriction. They do not rerun Harmony or recluster.

Then: validate final references/panel/taxonomy -> all three transfer arrays
(41 tasks each, concurrency 3 each) -> consensus table build -> consensus apply
(41 tasks, concurrency 6). Two evidence paths follow in parallel:

- broad scores (123 tasks, concurrency 6) -> threshold calibration;
- Unknown QC audit (41 tasks, concurrency 3) -> combined QC tables/PDF.

**Manual gate 1:** Review the newly generated threshold sensitivity, known-label
conflicts and Unknown QC evidence. Update/reapprove the active copies of:

- `config/reviewed_broad_reference_thresholds.csv` (all three references approved;
  do not assume old thresholds remain appropriate merely because approved=TRUE);
- `config/reviewed_unknown_qc.csv` (every current sample-qualified Unknown needs an
  explicit keep/remove decision and notes; old Unknown numbers can change).

Calibration tables are under
`outputs/validation/broad_reference_workflow/03_calibration/tables/`.
The Unknown workflow prints its combined review-template and PDF paths in the
compilation log. Do not populate decisions from old Unknown IDs.

## Phase 2: apply decisions, route cells, rebuild all regional branches

After every Phase-1 job has completed and the decisions above are reviewed:

```bash
bash scripts/submit_updated_references_to_merges.sh regions --dry-run
bash scripts/submit_updated_references_to_merges.sh regions --submit
```

This queues reviewed broad-threshold application and Unknown QC application,
then builds fresh routes. VZ and RL each run subset arrays -> merge -> fresh
integration -> QC review and Harmony audit. Residual runs subset array -> merge
-> integration -> morphology/QC review. All subset arrays use 41 tasks at
concurrency 3. No old QC removal list is applied.

**Manual gate 2:** Inspect each branch's newly generated QC evidence, markers,
sample composition and embeddings. Decide fresh VZ/RL removal IDs, or explicitly
`none`. For residual, choose a new approved removal decision bound to the new
source object, or explicitly decide no removal.

For residual removal, update `config/reviewed_residual_qc_removals.csv` using the
new QC per-cell/provenance tables. Cluster IDs, expected input/removal cell
counts, source size/mtime, review notes and approved status must match the new
`03_integrated/rds/Xenium_Residual_Integrated.rds`. The old 10,333-cell decision
and cluster 5 must not be reused automatically.

## Phase 3: post-QC processing and new subcluster-label review

Set these variables to the decisions from **this** review, then run the next
chain. Do not paste old cluster IDs. Use literal `none` for an approved VZ/RL
no-removal decision. Comma-separated values are exported through the environment
so Slurm's `--export` parser does not truncate the list.

```bash
# Replace each value with the newly reviewed decision before running.
export VZ_QC_REMOVE_CLUSTERS='REVIEWED_VZ_IDS_OR_none'
export RL_QC_REMOVE_CLUSTERS='REVIEWED_RL_IDS_OR_none'
export RESIDUAL_QC_ACTION='remove_OR_none'

bash scripts/submit_updated_references_to_merges.sh postqc --dry-run
bash scripts/submit_updated_references_to_merges.sh postqc --submit
```

The placeholder values above intentionally fail validation until replaced.
For `RESIDUAL_QC_ACTION=remove`, the newly approved residual CSV must exist.
For `none`, submission preserves any previous residual decision and
`outputs/xenium/residual/04_qc_filtered/` under this run's `residual_prior_qc/`
and records the explicit no-removal decision. Both possible decision locations
are archived so the loader cannot fall back to a stale scripts/config copy.
Paths and destinations are recorded in `moves.tsv`; no files are deleted.
Marker review then reads the new unfiltered integrated residual object.

VZ/RL apply QC -> post-QC processing -> resolution-0.8 marker review. The review
drivers build fresh replicate-aware evidence when needed. A zero-removal decision
uses the driver's supported reuse path instead of forcing reintegration.
Residual applies approved removals when requested, then regenerates marker
evidence and a blank label-review template. Its filtered embedding retains
pre-removal coordinates, as defined by the existing residual workflow.

**Manual gate 3:** Review new marker evidence and approve fresh cluster labels.

- Update `reviewed_label_maps` in `xenium_vz_06_annotate_subclusters.R` and
  `xenium_rl_06_annotate_subclusters.R` for the new resolution-0.8 cluster sets.
  The old VZ 1,447,477-cell guard is a source-specific review safeguard: replace
  it with a guard for the newly reviewed input, rather than deleting validation.
- Update the corresponding VZ expected map/cell total in
  `validate_xenium_regional_mapping.R`. Do not bypass the validator.
- Populate and approve `config/reviewed_residual_subcluster_labels.csv` from the
  newly generated residual template, including marker/sample evidence.
- If review introduces new fine labels, explicitly approve their broad map in
  `R/corrected_broad_labels.R` and add display colors/orders. Existing Vascular
  already maps to Vascular; other fine-label definitions remain unchanged.
- Parse every edited R script (normally under 10 seconds each), run structural
  `Rscript --max-ppsize=500000 scripts/validate_config.R` (under one minute), then
  run the next path/CSV preflight. Full file checking on HPC uses `--check-files`;
  typically 1–5 minutes, filesystem-dependent.

## Phase 4: apply approved labels and produce nonspatial and spatial merges

```bash
bash scripts/submit_updated_references_to_merges.sh merges --dry-run
bash scripts/submit_updated_references_to_merges.sh merges --submit
```

This queues VZ/RL label application and mapping to all samples, plus residual
label application and mapping. After all three maps complete: combine labels
-> corrected broad labels -> **both merges in parallel**. The combined launcher
validates all three mapping manifests; corrected-broad validation also runs
before merge jobs begin.

Final merged objects:

- `outputs/xenium/vz_rl/nonspatial/03_merged/rds/Xenium_VZRL_merged.rds`
- `outputs/xenium/vz_rl/spatial/03_merged/rds/Xenium_VZRL_spatial_merged.rds`

This is the requested merge endpoint. Subsequent nonspatial processing and
spatial integration are separate stages and are not silently included.

## Job records and completion checks

`outputs/reference_label_refresh_runs/$REFERENCE_REFRESH_RUN_ID/jobs.tsv` records
each successful submission immediately, including expected array-task counts.
The next phase requires every prior phase allocation/task to be COMPLETED with
ExitCode 0:0. A phase with recorded submissions cannot be blindly submitted again.
On failure, inspect its logs and recover the failed stage; do not start a second
overlapping rebuild. The configuration and manifest are snapshotted for each phase.

```bash
cat "outputs/reference_label_refresh_runs/$REFERENCE_REFRESH_RUN_ID/jobs.tsv"
job_ids=$(awk -F '\t' '{print $3}' \
  "outputs/reference_label_refresh_runs/$REFERENCE_REFRESH_RUN_ID/jobs.tsv" | paste -sd, -)
squeue -u acflint
sacct -j "$job_ids" \
  --format=JobID,JobName%30,Elapsed,State,MaxRSS,ReqMem,AllocCPUS,ExitCode --units=G
```

After Phase 4, validate regional/corrected-broad manifests (CSV/path checks,
under 1–2 minutes each, no large RDS loads):

```bash
Rscript --max-ppsize=500000 scripts/validate_xenium_regional_mapping.R VZ RL RESIDUAL
Rscript --max-ppsize=500000 scripts/validate_xenium_vz_rl_02_correct_broad.R
```

## Runtime estimates, excluding queue waits and human review

Existing launcher limits/memory/CPUs are reused. These are expected runtimes,
not Slurm limits. Arrays below use the explicit helper concurrency.

| Stage | Expected per task/job | Approximate array wall time | Basis/confidence |
|---|---|---|---|
| Conditional Science reviewed mapping/removal/reprocessing | 1–8 h | Adds 1–8 h before Science standardization if reviewed RDS is missing | Very low confidence; static inspection of 69,365 retained cells; requested limit 12 h |
| Reference standardization / panel preparation | 5–60 min per stage | Three two-stage chains overlap: 10–120 min | Static inspection; Science previously 2–3 min, others low confidence |
| Reference validation | 5–20 min | Single job | Static inspection, low |
| Label transfer | 8–50 min/task | 41 tasks/reference at 3: about 2–12 h/reference; three references may overlap | Prior 34-task arrays averaged 19–20 min/task; moderate for new labels |
| Consensus build | Under 2 min | Single job | Earlier workflow/history |
| Consensus apply, paired figures and markers | 15–90 min/task | 41 tasks at 6: 1.75–10.5 h | Revised figures/static inspection, low |
| Broad-score tables | 2–15 min/task | 123 tasks at 6: 0.7–5.25 h | Static inspection, low |
| Broad calibration / threshold application | 2–10 min each | Single jobs | Static inspection, low |
| Unknown QC audit | 1–3 min/task | 41 tasks at 3: 14–42 min | Earlier 42–75 s/task plus morphology I/O, low–moderate |
| Unknown compilation | 1–2 min | Single job | Earlier 28 s plus PDF assembly |
| Unknown QC application | 1–5 min/task | 41 tasks at 3: 14–70 min | Static inspection, low |
| Routing | 1–5 min | Single job | Small CSV tables, low |
| Regional subset | 1–5 min/task | 41 tasks/branch at 3: 14–70 min | RDS I/O, low |
| VZ/RL merge | 30–120 min each | Single jobs may overlap | 41 sequential RDS reads, low |
| VZ / RL integration | 4–16 h / 2–8 h | Single jobs may overlap | Related history; new routed cell sets, low |
| VZ / RL QC review | 1–3 h / 0.5–2 h | Single jobs may overlap | Prior VZ 76 min, related workloads |
| Harmony audit | 10–60 min/branch | Single jobs may overlap | Static inspection, low |
| Residual merge/integration/QC | 5–30 min / 15–120 min / 2–20 min | Single jobs | Prior ~10k-cell branch; new size unknown, low |
| VZ/RL QC application | 5–30 min each | Single jobs | Prior 13–25 min, moderate |
| VZ / RL post-QC processing after removals | 6–16 h / 3–8 h | Single jobs may overlap | VZ 373.5 min; RL 211.9 min; new cells, low |
| No-removal post-QC reuse | 5–30 min/branch | Single jobs | Supported reuse path, low |
| VZ/RL marker review including replicate evidence | 0.5–3 h/branch | Single jobs | Prior combined reviews 12–46 min; new sets, low |
| Residual QC apply / marker review | 1–10 min / 5–60 min | Single jobs | Smaller branch; new size unknown, low |
| VZ/RL label application / mapping | 10–60 min per job | Single jobs | Prior 11–19 min application, 13.5 min VZ mapping |
| Residual label application / mapping | 1–10 min / 10–40 min | Single jobs | Prior 0.67 / 22.1 min |
| Combined labels / corrected broad labels | 40–90 min / 10–30 min | Single jobs | Prior 44.95 / 12.93 min, low–moderate |
| Nonspatial / spatial merge | 0.5–4 h / 1–12 h | Single jobs may overlap | Large-object merge, spatial batching unbenchmarked, low |

Allow several days of active computation across all phases; human review, queue
waits, node caps and memory availability can dominate elapsed calendar time.
Do not shorten justified time limits to target a higher-priority partition.
Use the new sacct observations to refine future requests.

Before a source checkpoint, review `git status --short` and `git diff`, the
validation results, and only commit the completed logical unit when ready.
Keep generated outputs/logs/run ledgers out of Git. No package changes were made;
no renv snapshot is required solely for these script changes.
