#!/usr/bin/env Rscript
rm(list = ls())
suppressPackageStartupMessages(library(here))
source(here("scripts", "R", "config.R"))
source(here("scripts", "R", "xenium_harmony_audit.R"))
config <- load_pipeline_config()
args <- commandArgs(trailingOnly = TRUE)
valid <- c("--dry-run", "--overwrite")
if (any(!args %in% valid)) stop("Usage: Rscript scripts/xenium_vz_03a_audit_harmony.R [--dry-run|--overwrite]")
dry_run <- "--dry-run" %in% args; overwrite <- "--overwrite" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")
input <- here(config$project$outputs_dir, "xenium", "vz", "03_integrated", "rds", "Xenium_VZ_integrated.rds")
out <- here(config$project$outputs_dir, "xenium", "vz", "03_integrated", "tables", "harmony_audit")
outputs <- file.path(out, c("Xenium_VZ_Harmony_Audit_Summary.csv", "Xenium_VZ_Harmony_Audit_Paired_Metrics.csv", "Xenium_VZ_Harmony_Audit_Provenance.csv"))
if (dry_run) { compact_dry_run("VZ pre-/post-Harmony audit", inputs = input, outputs = outputs); quit(save="no", status=if(file.exists(input)) 0L else 1L) }
suppressPackageStartupMessages(library(Seurat))
run_harmony_audit("VZ", input, out, overwrite, seed = config$runtime$random_seed)
