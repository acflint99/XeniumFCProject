# Pipeline configuration

The 2026-10-22 reference-label refresh explicitly approves `eCN -> RL` and
`Vascular -> Other/nonregional`. Both `regional_subsets` in `config.yml` and
`reference_broad_taxonomy.csv` must agree. See
`../UPDATED_REFERENCES_TO_MERGES.md` for staged submission and manual review gates.

These files provide validated configuration metadata for the fetal cerebellum
Xenium pipeline. Split-slide preprocessing, annotation, consensus labelling,
and VZ/RL subset drivers read the sample manifest directly; most other analysis
scripts have not yet been migrated.

The consensus-label stage also reads `metadata/samples_meta.xlsx` and stores
the matched `PCW` value directly in every consensus object.

When this directory is moved to the HPC project root, the intended layout is:

```text
/home/acflint/R/Projects/XeniumFCProject/
├── config/
│   ├── cheaha_slurm_partitions.csv
│   ├── config.yml
│   ├── reviewed_unknown_qc.csv
│   ├── samples.csv
│   └── slides.csv
├── data/
├── metadata/
├── outputs/
├── scripts/
└── logs/
```

After moving it, scripts should resolve the files with:

```r
config_path <- here::here("config", "config.yml")
samples_path <- here::here("config", "samples.csv")
```

## Scope

All 41 biological samples are intended to undergo preprocessing, annotation,
VZ analysis, and RL analysis. The manifest therefore does not contain
stage-specific inclusion columns. Completion is determined from expected
outputs rather than stored as changing state in the configuration.

`samples.csv` has one row per biological sample. `slides.csv` has one row per
physical Xenium input directory and records which biological samples were
included from that directory. For the GZFB124 three-section slide, only section 2
is included; its declared count of one is the analyzed count. Each manually separated sample has its own
`cell_stats_file` within the shared input directory.

`config.yml` records split-slide preprocessing paths, runtime defaults, the
confirmed v2 reference RDS paths, and the consensus broad labels used for VZ
and RL subsetting. `regional_subsets.unknown_routing_method` records the
approved rule that a sample-qualified `Unknown-##` population enters only the
VZ or RL branch supported by at least two of the three existing reference
votes. Nonregional and unresolved populations remain outside both branches.
The configured label lists define how existing reference votes collapse to VZ,
RL, or nonregional routing categories.
VZ/RL integration settings still differ among scripts and should not be
unified without checking the intended scientific workflow.

`reviewed_unknown_qc.csv` is the explicit population-level QC decision
manifest. Its key is `sample_id` plus `consensus_label`; every row requires a
`keep` or `remove` decision and nonblank review notes. Removal applies to the
entire sample-qualified Unknown population. The conservative default during
review is `keep`.

`reviewed_unknown_relabels.csv` is retained for provenance of the earlier
dual-branch comparison workflow. It is not used by the reference-routed
production workflow because an Unknown population can enter at most one
regional branch.

`cheaha_slurm_partitions.csv` is a dated scheduling-policy snapshot rather
than pipeline configuration. It records the partition classes, node counts,
per-researcher node limits, maximum wall times, and priority tiers shown in the
Cheaha "Available Slurm Partitions" table supplied on 2026-09-01. Check the
current Cheaha information before relying on it for a later submission.

The configured v2 reference objects are stored below
`outputs/references/<study>/rds/`. Reference plots and VZ-specific reference
outputs use the parallel `plots/` and `vz/` directories.

## Validation

From the current local `scripts/` directory, run structural validation with:

```bash
Rscript validate_config.R
```

After moving `config/` to the HPC project root, run full input-file validation
from that project root with:

```bash
Rscript scripts/validate_config.R --check-files
```
