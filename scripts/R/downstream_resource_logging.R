# Observational logging only: R-managed memory is not enforced job memory or peak RSS.
log_downstream_memory <- function(stage, object = NULL) {
  usage <- gc(full = TRUE)
  object_gib <- if (is.null(object)) NA_real_ else as.numeric(object.size(object)) / 1024^3
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", stage,
          " | R-managed MiB=", round(sum(usage[, 2]), 1),
          " | object GiB=", round(object_gib, 2),
          " | job=", Sys.getenv("SLURM_JOB_ID", unset = "local"),
          " | allocated memory MiB=", Sys.getenv("SLURM_MEM_PER_NODE", unset = "unknown"))
  invisible(NULL)
}
