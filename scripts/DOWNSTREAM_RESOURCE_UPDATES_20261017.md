# Downstream resource updates — 2026-10-17

## Implemented changes

The current config/samples.csv has 41 samples. This implementation supersedes the recommendations in DOWNSTREAM_RESOURCE_REVIEW_20261016.md; that report remains a historical record of the pre-change code.

- VZ stage 5 defaults to one clustering worker in both the R driver and launcher, requests one CPU / 600 GB / 36 hours on medium, and retains explicit checkpoint-resume mode. No automatic resume or overwrite was enabled.
- RL stage 5 defaults to one clustering worker in both the R driver and launcher, requesting one CPU / 600 GB / 24 hours on medium. Its existing driver does not implement VZ-style checkpoint resume; none was added in this focused resource update.
- Stage-6 VZ/RL marker review and apply jobs default to one R worker and request one CPU / 256 GB / 12 hours on short. Raising allocated CPUs alone no longer silently enables parallel marker exports. Explicit VZ_SUBCLUSTER_WORKERS / RL_SUBCLUSTER_WORKERS overrides are bounded by allocated CPUs.
- VZ/RL/residual mapping default to one R worker at 128 GB / four-hour limits. Explicit VZ_MAPPING_WORKERS, RL_MAPPING_WORKERS, and RESIDUAL_MAPPING_WORKERS overrides are bounded by allocated CPUs. Sequential future evaluation is used for one worker.
- VZ/RL stage-6a evidence requests 256/128 GB respectively, retaining one CPU and four-hour limits.
- Combined merge/process/plot and spatial merge/integration provisionally request 600 GB. These requests provide headroom for full-cohort objects and transient copies but are not measured memory guarantees. Existing CPU counts and time limits remain unchanged.
- Combined and spatial integration use sequential future evaluation, including UMAP, so CPU allocation does not start independent full-object R sessions. Compiled routines may still use their existing thread configuration. UMAP may take longer than with parallel workers; no UMAP scientific parameters were changed.
- Combined plotting releases the completed Purkinje subset, completed plot objects, and the scaled heatmap temporary object before the next memory-intensive step. Marker selection and the order of scaling/downsampling are unchanged.
- Combined/spatial merge and processing log object size and R-managed memory after loads/merges and major preprocessing steps. The new R/downstream_resource_logging.R is required alongside those scripts. Logged R-managed memory is not Slurm RSS, peak allocation, or enforced cgroup usage; use sacct and logs together.
- The two-branch marker DotPlot array retains array=1-2%2, requesting one CPU / 256 GB per task. This array maps branches, not biological samples.
- RL diagnostic 5a now uses one CPU / 256 GB / 24 hours on medium; sensitivity 5b uses one CPU / 256 GB / 12 hours on short. Explicit 34-sample requirements were replaced with nonempty, unique, nonblank manifest checks. Cluster-6 sensitivity focus remains unchanged and needs biological review after new clustering.
- Residual mapping preflight no longer rejects a valid 41-sample manifest.
- VZ/RL stage-8 plotting and residual stage-1 subset launchers now default to array=1-41%3, matching the current manifest without reordering task IDs. A new submission helper derives bounds from each R driver's --list output to avoid future fixed-count drift. Concurrency remains three. Report-combine mode remains a separate single-task submission.

No changes were made under OLD/OLD-2, backups, Giotto or SpaTrack. No packages were installed or updated, no renv lockfile was changed, no configuration/sample manifest was edited, and no scientific analysis was run. Resource changes only affect future submissions; running/queued jobs do not inherit edited local scripts automatically.

## Expected runtimes versus requested limits

Estimates exclude queue wait and remain low confidence unless measured. Basis: static code inspection, 41 manifest samples, observed 1,501,411 VZ / 827,183 RL retained cells, failed parallel jobs, older sequential timeouts, and stage-4 input handling history. Sample count is not a substitute for cell/object size. Report elapsed time and peak memory after the first representative successful run before reducing requests or adding workers.

| Updated workload | Expected wall time | Requested limit | Memory / R workers by default |
|---|---|---|---|
| VZ stage-5 precluster resume | 8–24 hours | 36 hours | 600 GB / 1 |
| VZ stage-5 fresh run | 10–28 hours | 36 hours | 600 GB / 1 |
| RL stage-5 fresh run | 4–16 hours | 24 hours | 600 GB / 1 |
| RL diagnostic 5a | 3–16 hours | 24 hours | 256 GB / 1 |
| RL sensitivity 5b | 2–8 hours | 12 hours | 256 GB / 1 |
| VZ stage-6 review | 1–8 hours | 12 hours | 256 GB / 1 |
| RL stage-6 review | 0.5–6 hours | 12 hours | 256 GB / 1 |
| VZ stage-6 apply/plot modes | 15–120 minutes | 12 hours | 256 GB / 1 |
| RL stage-6 apply/plot modes | 10–90 minutes | 12 hours | 256 GB / 1 |
| VZ / RL stage-6a evidence | 15–90 / 10–60 minutes | 4 hours each | 256 / 128 GB, 1 each |
| VZ / RL / residual mapping | 20–120 minutes per complete branch | 4 hours each | 128 GB / 1 each |
| Combined non-spatial merge | 1–6 hours | 24 hours | 600 GB / 1 |
| Combined process | 8–36 hours | 48 hours | 600 GB / sequential futures, 8 CPUs |
| Combined plots | 1–6 hours | 24 hours | 600 GB / 1 |
| Spatial merge | 1–8 hours | 24 hours | 600 GB / 1 |
| Spatial integration | 12–40 hours | 48 hours | 600 GB / sequential futures, 8 CPUs |
| VZ/RL marker DotPlots | 10–60 min/task; 10–60 min total at concurrency 2 | 2 hours/task | 256 GB/task, 1 each |
| VZ/RL spatial report arrays | 5–30 min/task; 70–420 min total for 41 tasks at concurrency 3 | 2 hours/task | unchanged 64 GB/task, 1 each |
| Residual subset array | 2–15 min/task; 28–210 min total for 41 tasks at concurrency 3 | unchanged 30 min/task | unchanged 32 GB/task, 1 each |

Unchanged residual integration/review/label jobs still need residual cell/object-size measurements before their requests can be validated. No arbitrary RAM increase was applied without cohort-size evidence. Stage-9's intentional nine-sample developmental cohort remains unchanged; it still loads the entire branch object and needs memory monitoring (expected 5–30 minutes versus two-hour limits).

Cheaha partition policy was verified against current official documentation, https://docs.rc.uab.edu/cheaha/hardware/ (checked 2026-10-17): medium supports up to 50 hours; Intel nodes list 754 GB. Requested 600 GB / 24–48-hour jobs fit these documented limits; availability/fair-share and measured job peak still matter.

## Submission and safety

Sync the changed local files, including both new helpers, to the HPC project before using these settings. Keep the established project directory and R/4.4.1-foss-2023b environment. Do not rerun completed jobs solely because resource defaults changed.

From /home/acflint/R/Projects/XeniumFCProject, inspect sample task mappings without submitting anything:

```bash
bash scripts/submit_xenium_downstream_array.sh vz-plots --dry-run
bash scripts/submit_xenium_downstream_array.sh rl-plots --dry-run
bash scripts/submit_xenium_downstream_array.sh residual-subset --dry-run
```

Each helper preflight is expected to take under one minute and only reads manifest/configuration data via --list. It validates contiguous task IDs and unique nonblank sample IDs, prints the exact task map, and prints the proposed submission. Once inputs are ready, omit --dry-run to submit the selected stage. The helper does not validate spatial RDS contents or authorize overwrites. Missing input preflight/output manifests still need to be checked by the selected stage. Existing outputs are protected unless an explicit overwrite setting is supplied.

VZ recovery preflight (under one minute; paths and manifests only, no RDS load):

```bash
Rscript --max-ppsize=500000 scripts/xenium_vz_05_reintegrate_post_qc.R \
  --resume-checkpoint --dry-run
```

After HPC preflight passes and existing final outputs have been reviewed, a checkpoint recovery without overwriting final outputs is:

```bash
VZ_POST_QC_RESUME_CHECKPOINT=true \
  sbatch --export=ALL scripts/run_xenium_vz_05_reintegrate_post_qc.slurm
```

The real resume validates checkpoint input/parameter provenance when it loads the object. No scientific parameter, cluster-removal set, marker test/cap/threshold, seed, identity definition, label map, join cardinality, output filename, or overwrite rule was changed. Exact scientific call expressions were compared with the pre-edit code. Changing parallel execution can still expose stochastic/numerical differences; inspect resulting cluster/QC manifests before downstream labeling.

The five approved VZ removal clusters remain 23,24,25,27,28. The count-cohort definitions remain intact. New sample coverage requires completing the added samples' reports/subsets and rebuilding their combined summaries after validation; it does not require deleting or overwriting existing reports without review. Existing stage-3 objects and decision manifests remain the rollback basis. Keep the matched VZ precluster checkpoint for recovery.

## Validation performed locally

- All 14 new/changed R files parsed successfully (combined check under one second). Parse checks do not establish scientific correctness.
- Parsed call-expression comparisons confirmed that normalization, variable-feature selection, scaling, PCA, UMAP, Harmony, neighbor construction, clustering, marker tests, sketching and projection calls were unchanged across the 13 edited drivers.
- All 20 changed launchers plus the new submission helper passed bash -n (under one second combined).
- Stack-limit audit passed: 131 shell files and 200 Rscript invocations, each with --max-ppsize=500000 exactly once where required.
- All three array-helper dry runs produced exactly 41 task IDs with sample IDs in the authoritative manifest order and concurrency three. No job was submitted. The two-branch marker array was verified to remain 1-2%2.
- The memory logger ran successfully on a tiny in-memory object; no large object was loaded.
- All 13 changed R drivers executed their supported path-only dry runs (combined processing and plotting were both inspected). They reported missing HPC-only inputs, as expected. Local missing analysis RDS/review/route manifests prevent a readiness PASS. The RL sensitivity manifest checks passed 3/3; diagnostic/residual manifest validity no longer fails merely because there are 41 samples. Several existing preflights print FAIL but return exit zero; inspect the printed status as well as exit status. No change to that established behavior was made in this update.
- No configuration or sample-manifest change required validate_config.R. No package/environment change required renv::snapshot().

Git status could not be collected because this local scripts directory is not a Git checkout. In the actual repository, review git status and git diff, run the HPC path preflights, and checkpoint this logical update when ready. No commit or push was made.

## Changed files

- `run_xenium_residual_01_subset_samples.slurm`
- `run_xenium_residual_06_map_to_samples.slurm`
- `run_xenium_rl_05_reintegrate_post_qc.slurm`
- `run_xenium_rl_05a_diagnose_post_qc.slurm`
- `run_xenium_rl_05b_k_resolution_sensitivity.slurm`
- `run_xenium_rl_06_annotate_subclusters.slurm`
- `run_xenium_rl_06a_build_replicate_marker_evidence.slurm`
- `run_xenium_rl_07_map_to_samples.slurm`
- `run_xenium_rl_08_plot_samples.slurm`
- `run_xenium_vz_05_reintegrate_post_qc.slurm`
- `run_xenium_vz_06_annotate_subclusters.slurm`
- `run_xenium_vz_06a_build_replicate_marker_evidence.slurm`
- `run_xenium_vz_07_map_to_samples.slurm`
- `run_xenium_vz_08_plot_samples.slurm`
- `run_xenium_vz_rl_02_merge_samples.slurm`
- `run_xenium_vz_rl_03_plots_only.slurm`
- `run_xenium_vz_rl_03_process.slurm`
- `run_xenium_vz_rl_rl_vz_marker_dotplots.slurm`
- `run_xenium_vz_rl_spatial_01_merge.slurm`
- `run_xenium_vz_rl_spatial_02_integrate.slurm`
- `xenium_residual_06_map_to_samples.R`
- `xenium_rl_05_reintegrate_post_qc.R`
- `xenium_rl_05a_diagnose_post_qc.R`
- `xenium_rl_05b_k_resolution_sensitivity.R`
- `xenium_rl_06_annotate_subclusters.R`
- `xenium_rl_07_map_to_samples.R`
- `xenium_vz_05_reintegrate_post_qc.R`
- `xenium_vz_06_annotate_subclusters.R`
- `xenium_vz_07_map_to_samples.R`
- `xenium_vz_rl_02_merge_samples.R`
- `xenium_vz_rl_03_process_and_plot.R`
- `xenium_vz_rl_spatial_01_merge.R`
- `xenium_vz_rl_spatial_02_integrate.R`
- `R/downstream_resource_logging.R` (new)
- `submit_xenium_downstream_array.sh` (new)
- `README.md`
- `DOWNSTREAM_RESOURCE_UPDATES_20261017.md` (this file)

### Local preflight outcomes

| Driver | Exit | First status line |
|---|---:|---|
| xenium_vz_05_reintegrate_post_qc.R | 0 | DRY-RUN FAIL  /  VZ post-QC checkpoint resume  /  inputs=0/4  /  outputs-existing=0/7  /  checks=1/3 |
| xenium_rl_05_reintegrate_post_qc.R | 0 | DRY-RUN FAIL  /  RL post-QC reintegration  /  inputs=0/3  /  outputs-existing=0/6  /  checks=1/3 |
| xenium_rl_05a_diagnose_post_qc.R | 0 | DRY-RUN FAIL  /  RL stage-3 versus stage-5 diagnostic  /  inputs=0/5  /  outputs-existing=0/17  /  checks=1/3 |
| xenium_rl_05b_k_resolution_sensitivity.R | 0 | DRY-RUN FAIL  /  RL k.param/resolution sensitivity  /  inputs=0/2  /  outputs-existing=0/31  /  checks=3/3 |
| xenium_residual_06_map_to_samples.R | 1 | DRY-RUN FAIL  /  Residual label mapping  /  inputs=0/84  /  outputs-existing=0/42  /  checks=1/2 |
| xenium_vz_06_annotate_subclusters.R | 1 | DRY-RUN FAIL  /  VZ post-QC subcluster review  /  inputs=0/1  /  outputs-existing=0/8 |
| xenium_rl_06_annotate_subclusters.R | 1 | DRY-RUN FAIL  /  RL post-QC subcluster review  /  inputs=0/1  /  outputs-existing=0/8 |
| xenium_vz_07_map_to_samples.R | 0 | DRY-RUN FAIL  /  VZ label mapping  /  inputs=0/43  /  outputs-existing=0/42  /  checks=1/3 |
| xenium_rl_07_map_to_samples.R | 0 | DRY-RUN FAIL  /  RL label mapping  /  inputs=0/43  /  outputs-existing=0/42  /  checks=1/3 |
| xenium_vz_rl_spatial_02_integrate.R | 0 | DRY-RUN FAIL  /  Spatial regional integration  /  inputs=0/2  /  outputs-existing=0/5 |
| xenium_vz_rl_03_process_and_plot.R | 0 | DRY-RUN FAIL  /  Combined regional analysis [processing]  /  inputs=0/2  /  outputs-existing=0/12 |
| xenium_vz_rl_02_merge_samples.R | 0 | DRY-RUN FAIL  /  Combined regional merge  /  inputs=0/41  /  outputs-existing=0/43  /  checks=1/1 |
| xenium_vz_rl_spatial_01_merge.R | 0 | DRY-RUN FAIL  /  Spatial regional merge  /  inputs=0/41  /  outputs-existing=0/2  /  checks=1/1 |
