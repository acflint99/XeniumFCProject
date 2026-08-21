#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
check_files <- "--check-files" %in% args

config_utils <- c(
  here::here("R", "config.R"),
  here::here("scripts", "R", "config.R")
)
config_utils <- config_utils[file.exists(config_utils)]
if (length(config_utils) == 0L) stop("Could not find R/config.R or scripts/R/config.R.")

source(config_utils[[1]])
validate_pipeline_config(check_files = check_files)

