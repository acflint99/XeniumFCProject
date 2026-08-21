# Pipeline configuration

These files are the first, non-operative configuration draft for the fetal
cerebellum Xenium pipeline. No analysis script reads them yet.

When this directory is moved to the HPC project root, the intended layout is:

```text
/home/acflint/R/Projects/XeniumFCProject/
├── config/
│   ├── config.yml
│   └── samples.csv
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
must be reviewed before the manifest controls any job submission or analysis.
Blank means “not decided,” not `FALSE`.

`samples.csv` has one row per biological sample. `slides.csv` has one row per
physical Xenium input directory and records which biological samples were
captured together. The `cell_stats_file` field is intentionally blank where a
manually separated sample-specific cell list still needs to be identified.

`config.yml` records current defaults observed in the scripts. It is not yet a
source of truth. In particular, VZ/RL integration settings differ among
scripts and should not be unified without checking the intended scientific
workflow.

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
