# Pipeline configuration

These files provide validated configuration metadata for the fetal cerebellum
Xenium pipeline. The split-slide preprocessing driver reads them directly;
most other analysis scripts have not yet been migrated.

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

## Review required

All stage inclusion columns in `samples.csv` are intentionally blank. They
must be reviewed before they control stage-specific job submission. Blank
means “not decided,” not `FALSE`. The split-slide preprocessing driver selects
rows by `input_layout` and does not use these inclusion columns.

`samples.csv` has one row per biological sample. `slides.csv` has one row per
physical Xenium input directory and records which biological samples were
captured together. Each manually separated sample has its own
`cell_stats_file` within the shared input directory.

`config.yml` is the source of truth for the split-slide preprocessing paths
and runtime defaults. It is not yet the source of truth for the full pipeline.
In particular, VZ/RL integration settings differ among scripts and should not
be unified without checking the intended scientific workflow.

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
