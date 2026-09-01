#!/usr/bin/env Rscript

# Merge Xenium cluster-label comparison tables from the Aldinger, Sepp, and
# Science reference runs. Production mode is retained for the 34-sample
# workflow; --pilot-res2 through --pilot-res5 use only the three
# pilot samples and write beneath their respective isolated output trees.

suppressPackageStartupMessages(library(here))

source(here("scripts", "R", "config.R"))

config <- load_pipeline_config()
pilot_manifest <- load_resolution2_pilot_manifest(config)
sample_manifest <- load_sample_manifest(config)
references <- c("Aldinger", "Sepp", "Science")

args <- commandArgs(trailingOnly = TRUE)
valid_options <- c(
  "--pilot-res2", "--pilot-res3", "--pilot-res4", "--pilot-res5", "--dry-run",
  "--overwrite", "--weighted-2of3", "--all-samples-res4", "--all-samples-res5",
  "--selected-sample",
  "--list"
)
unknown_options <- args[startsWith(args, "--") & !args %in% valid_options]
if (length(unknown_options)) {
  stop("Unknown option(s): ", paste(unknown_options, collapse = ", "))
}

pilot_res2 <- "--pilot-res2" %in% args
pilot_res3 <- "--pilot-res3" %in% args
pilot_res4 <- "--pilot-res4" %in% args
pilot_res5 <- "--pilot-res5" %in% args
all_samples_res4 <- "--all-samples-res4" %in% args
all_samples_res5 <- "--all-samples-res5" %in% args
all_samples_mode <- all_samples_res4 || all_samples_res5
selected_sample_mode <- "--selected-sample" %in% args
pilot_flags <- c(pilot_res2, pilot_res3, pilot_res4, pilot_res5)
if (sum(pilot_flags) > 1L) {
  stop(
    "Choose only one of --pilot-res2, --pilot-res3, --pilot-res4, or --pilot-res5."
  )
}
pilot_mode <- any(pilot_flags)
if (all_samples_res4 && all_samples_res5) {
  stop("Choose only one of --all-samples-res4 or --all-samples-res5.")
}
if (all_samples_mode && pilot_mode) {
  stop("An all-samples mode cannot be combined with a pilot-resolution option.")
}
if (selected_sample_mode && (!pilot_res5 || all_samples_mode)) {
  stop("--selected-sample currently requires --pilot-res5.")
}
selected_sample_id <- if (selected_sample_mode) {
  trimws(Sys.getenv("SELECTED_SAMPLE_ID", unset = ""))
} else {
  ""
}
if (selected_sample_mode) {
  selected_matches <- which(sample_manifest$sample_id == selected_sample_id)
  if (!nzchar(selected_sample_id) || length(selected_matches) != 1L) {
    stop(
      "SELECTED_SAMPLE_ID must match exactly one config/samples.csv row; found ",
      length(selected_matches), " match(es) for '", selected_sample_id, "'."
    )
  }
}
resolution_mode <- pilot_mode || all_samples_mode
pilot_resolution_tag <- if (all_samples_res5) "5.0" else if (all_samples_res4) "4.0" else if (pilot_res5) "5.0" else if (pilot_res4) "4.0" else if (pilot_res3) "3.0" else "2.0"
pilot_stage <- if (selected_sample_mode) {
  "03h_resolution5_selected_sample"
} else if (all_samples_res4) {
  "03g_resolution4_all_samples"
} else if (all_samples_res5) {
  "03i_resolution5_all_samples"
} else if (pilot_res5) {
  "03e_resolution5_pilot"
} else if (pilot_res4) {
  "03d_resolution4_pilot"
} else if (pilot_res3) {
  "03c_resolution3_pilot"
} else {
  "03b_resolution2_pilot"
}
weighted_2of3 <- "--weighted-2of3" %in% args
consensus_method <- if (weighted_2of3) "weighted_2of3" else "legacy_six_vote"
dry_run <- "--dry-run" %in% args
overwrite <- "--overwrite" %in% args
list_requested <- "--list" %in% args
if (dry_run && overwrite) stop("--dry-run and --overwrite cannot be combined.")

positional_args <- args[!args %in% valid_options]
if (length(positional_args) > 1L) {
  stop(
    "Usage: Rscript scripts/xenium_annotate_02_build_consensus.R ",
    "[--pilot-res2|--pilot-res3|--pilot-res4|--pilot-res5|--all-samples-res4|--all-samples-res5] ",
    "[--selected-sample] ",
    "[--weighted-2of3] ",
    "[--list|--dry-run|--overwrite] [OUTPUT_ROOT]"
  )
}
output_root <- if (length(positional_args)) positional_args[[1]] else "outputs"

annotation_root <- if (all_samples_res4) {
  file.path(output_root, "xenium", "annotation", "resolution4_all_samples")
} else if (all_samples_res5) {
  file.path(output_root, "xenium", "annotation", "resolution5_all_samples")
} else if (pilot_mode) {
  file.path(
    output_root, "xenium", "preprocess", pilot_stage, "annotation"
  )
} else {
  file.path(output_root, "xenium", "annotation")
}
table_dirs <- setNames(
  file.path(annotation_root, "01_label_transfer", tolower(references), "tables"),
  references
)
consensus_stage <- if (weighted_2of3) {
  "02_consensus_weighted_2of3"
} else {
  "02_consensus"
}
merged_dir <- file.path(annotation_root, consensus_stage, "tables")

file_suffix <- function(reference) {
  paste0("_", reference, "_majority_vs_weighted.csv")
}

pilot_samples <- as.character(pilot_manifest$sample_id)
if (anyDuplicated(pilot_samples) || any(!nzchar(pilot_samples))) {
  stop("Clustering pilot manifest contains duplicate or blank sample IDs.")
}
all_samples <- as.character(sample_manifest$sample_id)
if (anyDuplicated(all_samples) || any(!nzchar(all_samples))) {
  stop("config/samples.csv contains duplicate or blank sample IDs.")
}
mode_description <- if (selected_sample_mode) {
  paste0("resolution-5.0 selected-sample pilot: ", selected_sample_id)
} else if (all_samples_res4) {
  "resolution-4.0 all-sample analysis"
} else if (all_samples_res5) {
  "resolution-5.0 all-sample analysis"
} else if (pilot_mode) {
  paste0("resolution-", pilot_resolution_tag, " pilot")
} else {
  "production"
}

if (list_requested) {
  if (length(positional_args)) stop("--list does not accept OUTPUT_ROOT.")
  sample_names <- if (selected_sample_mode) {
    selected_sample_id
  } else if (pilot_mode) {
    pilot_samples
  } else {
    all_samples
  }
  cat(
    "Mode:", mode_description, "\n"
  )
  cat("Consensus method:", consensus_method, "\n")
  write.table(
    data.frame(task_id = seq_along(sample_names), sample_id = sample_names),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
  quit(save = "no", status = 0L)
}

comparison_files <- lapply(references, function(reference) {
  list.files(
    table_dirs[[reference]],
    pattern = paste0("_", reference, "_majority_vs_weighted\\.csv$"),
    full.names = TRUE
  )
})
names(comparison_files) <- references

if (selected_sample_mode) {
  sample_names <- selected_sample_id
} else if (pilot_mode) {
  sample_names <- pilot_samples
} else if (all_samples_mode) {
  sample_names <- all_samples
} else {
  sample_names <- sort(unique(unlist(Map(
    function(files, reference) {
      sub(file_suffix(reference), "", basename(files), fixed = TRUE)
    },
    comparison_files,
    references
  ))))
}
if (!length(sample_names)) {
  stop("No '*_majority_vs_weighted.csv' files were found.")
}

input_grid <- do.call(rbind, lapply(sample_names, function(sample_name) {
  data.frame(
    sample_id = sample_name,
    reference = references,
    input = vapply(
      references,
      function(reference) {
        file.path(
          table_dirs[[reference]],
          paste0(sample_name, file_suffix(reference))
        )
      },
      character(1)
    ),
    stringsAsFactors = FALSE
  )
}))
output_paths <- file.path(merged_dir, paste0(sample_names, "_comparison_merged.csv"))

if (dry_run) {
  compact_dry_run(
    paste0("Consensus build [", consensus_method, "]"),
    inputs = input_grid$input,
    outputs = output_paths
  )
  inputs_ready <- all(file.exists(input_grid$input))
  quit(save = "no", status = if (inputs_ready) 0L else 1L)
}

missing_dirs <- table_dirs[!dir.exists(table_dirs)]
if (length(missing_dirs)) {
  stop(
    "These reference table directories do not exist:\n",
    paste(missing_dirs, collapse = "\n")
  )
}

if (resolution_mode && any(!file.exists(input_grid$input))) {
  stop(
    "The resolution-", pilot_resolution_tag,
    if (all_samples_mode) {
      " all-sample analysis requires all 102 reference/sample tables. Missing:\n- "
    } else if (selected_sample_mode) {
      " selected-sample pilot requires all three reference tables. Missing:\n- "
    } else {
      " pilot requires all nine reference/sample tables. Missing:\n- "
    },
    paste(input_grid$input[!file.exists(input_grid$input)], collapse = "\n- ")
  )
}

existing_outputs <- output_paths[file.exists(output_paths)]
if (length(existing_outputs) && !overwrite) {
  stop(
    "Refusing to overwrite existing consensus tables:\n- ",
    paste(existing_outputs, collapse = "\n- "),
    "\nRerun with --overwrite only after reviewing these files."
  )
}
if (length(existing_outputs)) {
  warning("Overwriting ", length(existing_outputs), " existing consensus table(s).")
}

read_reference_table <- function(sample_name, reference) {
  path <- file.path(
    table_dirs[[reference]],
    paste0(sample_name, file_suffix(reference))
  )
  if (!file.exists(path)) return(NULL)

  tab <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("seurat_clusters", "cluster_majority", "cluster_weighted")
  if (!all(required %in% names(tab))) {
    stop("Unexpected columns in ", path, ". Expected: ", paste(required, collapse = ", "))
  }
  tab <- tab[required]
  tab$seurat_clusters <- as.character(tab$seurat_clusters)
  if (!nrow(tab)) stop("Comparison table has no cluster rows: ", path)
  if (anyNA(tab$seurat_clusters) || any(!nzchar(tab$seurat_clusters))) {
    stop("Comparison table has blank cluster IDs: ", path)
  }
  if (anyDuplicated(tab$seurat_clusters)) {
    stop("Comparison table has duplicate cluster IDs: ", path)
  }

  names(tab) <- c(
    "seurat_clusters",
    paste0(tolower(reference), "_cluster_majority"),
    paste0(tolower(reference), "_cluster_weighted")
  )
  tab
}

merged_tables <- vector("list", length(sample_names))
names(merged_tables) <- sample_names
for (sample_name in sample_names) {
  tables <- lapply(references, function(reference) {
    read_reference_table(sample_name, reference)
  })
  names(tables) <- references

  absent <- names(tables)[vapply(tables, is.null, logical(1))]
  if (length(absent)) {
    warning(
      sample_name, ": no comparison table found for ", paste(absent, collapse = ", "),
      "; corresponding columns will be NA."
    )
    for (reference in absent) {
      empty_table <- data.frame(seurat_clusters = character(), check.names = FALSE)
      empty_table[[paste0(tolower(reference), "_cluster_majority")]] <- character()
      empty_table[[paste0(tolower(reference), "_cluster_weighted")]] <- character()
      tables[[reference]] <- empty_table
    }
  }

  if (resolution_mode) {
    cluster_sets <- lapply(tables, function(tab) sort(tab$seurat_clusters))
    if (!all(vapply(cluster_sets[-1], identical, logical(1), cluster_sets[[1]]))) {
      stop(
        "Reference tables do not contain identical resolution-", pilot_resolution_tag,
        " cluster IDs for ",
        sample_name, "."
      )
    }
  }

  # Each input was checked for unique cluster IDs, so these are one-to-one joins.
  merged <- Reduce(
    function(x, y) merge(x, y, by = "seurat_clusters", all = TRUE, sort = FALSE),
    tables
  )
  if (anyDuplicated(merged$seurat_clusters)) {
    stop("Consensus join unexpectedly duplicated cluster IDs for ", sample_name, ".")
  }

  if (weighted_2of3) {
    majority_columns <- paste0(tolower(references), "_cluster_majority")
    weighted_columns <- paste0(tolower(references), "_cluster_weighted")
    if (!all(c(majority_columns, weighted_columns) %in% names(merged))) {
      stop("Weighted 2-of-3 consensus lacks required reference columns for ", sample_name, ".")
    }

    # Each biological reference contributes exactly one vote: its weighted
    # cluster label. Majority/weighted agreement is retained only as an audit
    # measure and is not counted as another independent vote.
    for (i in seq_along(references)) {
      reference_key <- tolower(references[[i]])
      majority_values <- trimws(as.character(merged[[majority_columns[[i]]]]))
      weighted_values <- trimws(as.character(merged[[weighted_columns[[i]]]]))
      merged[[paste0(reference_key, "_reference_vote")]] <- weighted_values
      merged[[paste0(reference_key, "_majority_weighted_agree")]] <- (
        !is.na(majority_values) & nzchar(majority_values) &
          !is.na(weighted_values) & nzchar(weighted_values) &
          majority_values == weighted_values
      )
    }

    vote_columns <- paste0(tolower(references), "_reference_vote")
    consensus_result <- lapply(seq_len(nrow(merged)), function(i) {
      votes <- trimws(as.character(unlist(merged[i, vote_columns, drop = FALSE])))
      valid_votes <- votes[!is.na(votes) & nzchar(votes) & tolower(votes) != "unknown"]
      if (!length(valid_votes)) {
        return(list(label = "Unknown", support = 0L, nonunknown = 0L,
                    status = "insufficient_support"))
      }

      counts <- table(valid_votes)
      support <- max(counts)
      winners <- names(counts)[counts == support]
      if (support < 2L || length(winners) != 1L) {
        return(list(label = "Unknown", support = as.integer(support),
                    nonunknown = length(valid_votes), status = "insufficient_support"))
      }

      status <- if (support == 3L) "unanimous_3of3" else "majority_2of3"
      list(label = winners[[1]], support = as.integer(support),
           nonunknown = length(valid_votes), status = status)
    })
    merged$consensus_label <- vapply(consensus_result, `[[`, character(1), "label")
    merged$consensus_support_n <- vapply(consensus_result, `[[`, integer(1), "support")
    merged$consensus_nonunknown_n <- vapply(consensus_result, `[[`, integer(1), "nonunknown")
    merged$consensus_status <- vapply(consensus_result, `[[`, character(1), "status")
    merged$consensus_method <- consensus_method
  } else {
    # Preserve the established legacy rule exactly: ignore Unknown/missing
    # labels, then choose the most common value across all six columns.
    label_columns <- setdiff(names(merged), "seurat_clusters")
    merged$consensus_label <- vapply(seq_len(nrow(merged)), function(i) {
      labels <- trimws(as.character(unlist(merged[i, label_columns, drop = FALSE])))
      labels <- labels[!is.na(labels) & nzchar(labels) & tolower(labels) != "unknown"]
      if (!length(labels)) return("Unknown")

      counts <- table(labels)
      winners <- names(counts)[counts == max(counts)]
      labels[match(winners, labels)][1]
    }, character(1))
  }

  cluster_number <- suppressWarnings(as.numeric(as.character(merged$seurat_clusters)))
  merged <- merged[
    order(is.na(cluster_number), cluster_number, merged$seurat_clusters),
    ,
    drop = FALSE
  ]
  unknown_rows <- which(tolower(trimws(as.character(merged$consensus_label))) == "unknown")
  if (length(unknown_rows)) {
    merged$consensus_label[unknown_rows] <- sprintf(
      "Unknown-%02d", seq_along(unknown_rows)
    )
  }
  merged_tables[[sample_name]] <- merged
}

# Do not create or partially update the destination until every input and join
# has passed validation.
dir.create(merged_dir, recursive = TRUE, showWarnings = FALSE)
for (sample_name in sample_names) {
  out_path <- file.path(merged_dir, paste0(sample_name, "_comparison_merged.csv"))
  write.csv(merged_tables[[sample_name]], out_path, row.names = FALSE, na = "")
  message("Wrote: ", out_path)
}

unknown_consensus <- do.call(rbind, lapply(sample_names, function(sample_name) {
  merged <- merged_tables[[sample_name]]
  unknown <- merged[
    grepl("^Unknown-[0-9]{2,}$", merged$consensus_label),
    c("seurat_clusters", "consensus_label"),
    drop = FALSE
  ]
  if (!nrow(unknown)) return(NULL)
  data.frame(
    sample = sample_name,
    seurat_clusters = unknown$seurat_clusters,
    consensus_label = unknown$consensus_label
  )
}))

if (is.null(unknown_consensus) || !nrow(unknown_consensus)) {
  message("No samples had Unknown consensus clusters.")
} else {
  print(unknown_consensus)
}
