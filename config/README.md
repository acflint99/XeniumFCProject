# Pipeline configuration

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
│   ├── config.yml
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

All 34 biological samples are intended to undergo preprocessing, annotation,
VZ analysis, and RL analysis. The manifest therefore does not contain
stage-specific inclusion columns. Completion is determined from expected
outputs rather than stored as changing state in the configuration.

`samples.csv` has one row per biological sample. `slides.csv` has one row per
physical Xenium input directory and records which biological samples were
captured together. Each manually separated sample has its own
`cell_stats_file` within the shared input directory.

`config.yml` records split-slide preprocessing paths, runtime defaults, the
confirmed v2 reference RDS paths, and the consensus broad labels used for VZ
and RL subsetting. It is not yet the source of truth for the full pipeline.
VZ/RL integration settings still differ among scripts and should not be
unified without checking the intended scientific workflow.

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
