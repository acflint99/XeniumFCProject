#!/usr/bin/env Rscript
# Residual entry point for the same clustree/stability review used by VZ/RL.
# Expected 5-60 min for ~10k cells / five repeats; low confidence, HPC only.
# --list/--dry-run inspect paths and small configuration tables (<1 min).
suppressPackageStartupMessages(library(here))
here::i_am("scripts/xenium_residual_05c_resolution_stability.R")
source(here("scripts", "xenium_regional_05c_resolution_stability.R"))
main(c("--region=Residual", commandArgs(trailingOnly = TRUE)))
