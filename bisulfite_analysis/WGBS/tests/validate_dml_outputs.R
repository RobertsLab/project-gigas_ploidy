arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L) {
  stop("Usage: validate_dml_outputs.R RUN_DIRECTORY", call. = FALSE)
}
run_directory <- normalizePath(arguments[[1]], mustWork = TRUE)

coordinate_key <- function(data) paste(data$chr, data$start, data$end, sep = "\t")

all_paths <- list.files(
  run_directory,
  pattern = "^dml_stats_.*_all[.]csv$",
  full.names = TRUE
)
if (length(all_paths) != 1L) {
  stop("Expected exactly one complete, unfiltered DML statistics CSV.", call. = FALSE)
}
all_stats <- read.csv(all_paths, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("chr", "start", "end", "pvalue", "qvalue", "meth.diff")
if (!all(required %in% names(all_stats))) {
  stop("The complete DML statistics CSV lacks required columns.", call. = FALSE)
}
all_keys <- coordinate_key(all_stats)
if (anyDuplicated(all_keys) || any(all_stats$start < 0) || any(all_stats$end <= all_stats$start)) {
  stop("The complete DML statistics CSV contains invalid or duplicate coordinates.", call. = FALSE)
}
if (any(!is.finite(all_stats$meth.diff)) || any(abs(all_stats$meth.diff) > 100)) {
  stop("The complete DML statistics CSV contains invalid methylation differences.", call. = FALSE)
}

subset_paths <- list.files(
  run_directory,
  pattern = "^dml_.*_diff[0-9.]+_q[0-9.]+[.]csv$",
  full.names = TRUE
)
if (!length(subset_paths)) {
  stop("No thresholded DML CSV files were found.", call. = FALSE)
}
observed_counts <- data.frame(
  subset = "all_tested_cpgs",
  rows = nrow(all_stats),
  stringsAsFactors = FALSE
)
for (csv_path in subset_paths) {
  stem <- sub("[.]csv$", "", basename(csv_path))
  matched <- regexec("_diff([0-9.]+)_q([0-9.]+)$", stem)
  values <- regmatches(stem, matched)[[1]]
  if (length(values) != 3L) {
    stop(sprintf("Cannot parse thresholds from %s.", basename(csv_path)), call. = FALSE)
  }
  difference <- as.numeric(values[[2]])
  max_q <- as.numeric(values[[3]])
  subset <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(required %in% names(subset))) {
    stop(sprintf("Thresholded DML CSV lacks required columns: %s", basename(csv_path)), call. = FALSE)
  }
  subset_keys <- coordinate_key(subset)
  if (anyDuplicated(subset_keys) || any(!subset_keys %in% all_keys)) {
    stop(sprintf("Thresholded DML coordinates are invalid: %s", basename(csv_path)), call. = FALSE)
  }
  if (nrow(subset) &&
      (any(is.na(subset$qvalue)) || any(subset$qvalue > max_q) ||
       any(abs(subset$meth.diff) < difference))) {
    stop(sprintf("Threshold violations found in %s.", basename(csv_path)), call. = FALSE)
  }

  bed_path <- file.path(run_directory, paste0(stem, ".bed"))
  if (!file.exists(bed_path)) {
    stop(sprintf("Missing paired BED file for %s.", basename(csv_path)), call. = FALSE)
  }
  if (file.info(bed_path)$size == 0) {
    bed <- data.frame(chr = character(), start = integer(), end = integer())
  } else {
    bed <- read.delim(
      bed_path,
      header = FALSE,
      col.names = c("chr", "start", "end"),
      stringsAsFactors = FALSE
    )
  }
  if (!identical(coordinate_key(bed), subset_keys)) {
    stop(sprintf("CSV/BED coordinate disagreement for %s.", stem), call. = FALSE)
  }
  observed_counts <- rbind(
    observed_counts,
    data.frame(subset = stem, rows = nrow(subset), stringsAsFactors = FALSE)
  )
}

count_path <- file.path(run_directory, "output_counts.csv")
if (!file.exists(count_path)) {
  stop("Missing output_counts.csv.", call. = FALSE)
}
declared_counts <- read.csv(count_path, stringsAsFactors = FALSE)
observed_counts <- observed_counts[order(observed_counts$subset), , drop = FALSE]
declared_counts <- declared_counts[order(declared_counts$subset), , drop = FALSE]
rownames(observed_counts) <- NULL
rownames(declared_counts) <- NULL
if (!identical(observed_counts, declared_counts)) {
  stop("output_counts.csv does not match the generated DML tables.", call. = FALSE)
}

message(sprintf("DML output validation passed for %d tested CpGs and %d threshold pairs.", nrow(all_stats), length(subset_paths)))
