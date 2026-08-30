#!/usr/bin/env Rscript

options(scipen = 999)

parse_arguments <- function(arguments) {
  values <- list()
  for (argument in arguments) {
    if (!grepl("^--[^=]+=.*$", argument)) {
      stop(sprintf("Invalid argument '%s'; expected --name=value.", argument), call. = FALSE)
    }
    parts <- strsplit(sub("^--", "", argument), "=", fixed = TRUE)[[1]]
    key <- parts[[1]]
    value <- paste(parts[-1], collapse = "=")
    if (!is.null(values[[key]])) {
      stop(sprintf("Argument '--%s' was supplied more than once.", key), call. = FALSE)
    }
    values[[key]] <- value
  }

  required <- c(
    "metadata", "cov-dir", "cov-template", "output-dir", "assembly",
    "min-read-coverage", "min-analysis-coverage", "high-coverage-percentile",
    "normalize-coverage", "require-all-samples", "overdispersion", "test",
    "adjust", "effect", "max-q", "differences", "r-version",
    "bioconductor-version", "methylkit-version", "chunk-size", "mc-cores",
    "preflight-only"
  )
  missing <- setdiff(required, names(values))
  unexpected <- setdiff(names(values), required)
  if (length(missing)) {
    stop(sprintf("Missing arguments: %s", paste(paste0("--", missing), collapse = ", ")), call. = FALSE)
  }
  if (length(unexpected)) {
    stop(sprintf("Unknown arguments: %s", paste(paste0("--", unexpected), collapse = ", ")), call. = FALSE)
  }
  values
}

parse_flag <- function(value, name) {
  if (!value %in% c("true", "false")) {
    stop(sprintf("--%s must be true or false.", name), call. = FALSE)
  }
  identical(value, "true")
}

parse_number <- function(value, name, integer = FALSE, minimum = -Inf, maximum = Inf) {
  result <- suppressWarnings(as.numeric(value))
  if (length(result) != 1L || is.na(result) || result < minimum || result > maximum) {
    stop(sprintf("--%s is outside its valid range.", name), call. = FALSE)
  }
  if (integer && result != floor(result)) {
    stop(sprintf("--%s must be an integer.", name), call. = FALSE)
  }
  if (integer) as.integer(result) else result
}

validate_cov_preview <- function(path, seq_id, preview_lines = 100L) {
  lines <- readLines(path, n = preview_lines, warn = FALSE)
  if (!length(lines)) {
    stop(sprintf("COV file is empty for %s: %s", seq_id, path), call. = FALSE)
  }
  fields <- strsplit(lines, "\t", fixed = TRUE)
  field_counts <- lengths(fields)
  if (any(field_counts != 6L)) {
    stop(
      sprintf("COV preview for %s must have exactly six tab-separated columns.", seq_id),
      call. = FALSE
    )
  }
  values <- do.call(rbind, fields)
  if (any(!nzchar(values[, 1])) || any(grepl("[[:space:]]", values[, 1]))) {
    stop(sprintf("COV preview contains an invalid chromosome for %s.", seq_id), call. = FALSE)
  }
  numeric_values <- suppressWarnings(matrix(
    as.numeric(values[, 2:6, drop = FALSE]),
    nrow = nrow(values),
    ncol = 5L
  ))
  if (anyNA(numeric_values) || any(!is.finite(numeric_values))) {
    stop(sprintf("COV preview contains a non-numeric value for %s.", seq_id), call. = FALSE)
  }

  start <- numeric_values[, 1]
  end <- numeric_values[, 2]
  percent <- numeric_values[, 3]
  methylated <- numeric_values[, 4]
  unmethylated <- numeric_values[, 5]
  integer_fields <- cbind(start, end, methylated, unmethylated)
  if (any(integer_fields != floor(integer_fields)) || any(integer_fields < 0)) {
    stop(sprintf("COV coordinates and counts must be non-negative integers for %s.", seq_id), call. = FALSE)
  }
  if (any(end - start != 2)) {
    stop(sprintf("Merged CpG intervals must be zero-based, half-open, and two bases wide for %s.", seq_id), call. = FALSE)
  }
  if (any(percent < 0 | percent > 100) || any(methylated + unmethylated <= 0)) {
    stop(sprintf("COV percentages or coverage counts are invalid for %s.", seq_id), call. = FALSE)
  }
  recalculated <- 100 * methylated / (methylated + unmethylated)
  if (any(abs(percent - recalculated) > 0.00001)) {
    stop(sprintf("COV percentages disagree with methylated/unmethylated counts for %s.", seq_id), call. = FALSE)
  }
  invisible(TRUE)
}

write_bed <- function(data, path) {
  if (!nrow(data)) {
    file.create(path)
    return(invisible(path))
  }
  if (any(data$start < 0) || any(data$end <= data$start)) {
    stop(sprintf("Cannot write invalid BED coordinates to %s.", path), call. = FALSE)
  }
  bed <- data.frame(chr = data$chr, start = data$start, end = data$end)
  data.table::fwrite(bed, path, sep = "\t", col.names = FALSE, quote = FALSE)
  invisible(path)
}

main <- function() {
  args <- parse_arguments(commandArgs(trailingOnly = TRUE))
  metadata_path <- normalizePath(args[["metadata"]], mustWork = TRUE)
  cov_directory <- normalizePath(args[["cov-dir"]], mustWork = TRUE)
  cov_template <- args[["cov-template"]]
  output_directory <- args[["output-dir"]]
  preflight_only <- parse_flag(args[["preflight-only"]], "preflight-only")
  normalize_coverage <- parse_flag(args[["normalize-coverage"]], "normalize-coverage")
  require_all_samples <- parse_flag(args[["require-all-samples"]], "require-all-samples")
  min_read_coverage <- parse_number(args[["min-read-coverage"]], "min-read-coverage", TRUE, 1)
  min_analysis_coverage <- parse_number(args[["min-analysis-coverage"]], "min-analysis-coverage", TRUE, 1)
  high_coverage_percentile <- parse_number(
    args[["high-coverage-percentile"]], "high-coverage-percentile", FALSE, 0, 100
  )
  max_q <- parse_number(args[["max-q"]], "max-q", FALSE, 0, 1)
  differences <- suppressWarnings(as.numeric(strsplit(args[["differences"]], ",", fixed = TRUE)[[1]]))
  if (!length(differences) || anyNA(differences) || any(differences <= 0 | differences > 100)) {
    stop("--differences must be a comma-separated list of values in (0, 100].", call. = FALSE)
  }
  if (anyDuplicated(differences)) {
    stop("--differences must not contain duplicates.", call. = FALSE)
  }
  chunk_size <- parse_number(args[["chunk-size"]], "chunk-size", TRUE, 1)
  mc_cores <- parse_number(args[["mc-cores"]], "mc-cores", TRUE, 1)

  if (!grepl("{seq_id}", cov_template, fixed = TRUE)) {
    stop("--cov-template must contain the literal placeholder {seq_id}.", call. = FALSE)
  }
  if (dir.exists(output_directory) && length(list.files(output_directory, all.files = TRUE, no.. = TRUE))) {
    stop(sprintf("Output directory already exists and is not empty: %s", output_directory), call. = FALSE)
  }
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  output_directory <- normalizePath(output_directory, mustWork = TRUE)
  run_id <- basename(output_directory)
  if (!grepl("^[A-Za-z0-9._-]+$", run_id)) {
    stop("The output-directory basename must be a safe run identifier.", call. = FALSE)
  }

  metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
  required_metadata <- c(
    "seq_id", "library_name", "tissue", "ploidy", "desiccation", "heat_shock",
    "library_kit", "sra_bioproject", "sra_accession"
  )
  if (!identical(names(metadata), required_metadata)) {
    stop("Canonical WGBS metadata columns or column order are invalid.", call. = FALSE)
  }
  if (nrow(metadata) != 10L || anyDuplicated(metadata$seq_id) || anyNA(metadata)) {
    stop("Canonical metadata must contain ten complete, unique WGBS samples.", call. = FALSE)
  }
  if (!all(metadata$ploidy %in% c("diploid", "triploid")) ||
      !all(metadata$heat_shock %in% c("no", "yes"))) {
    stop("Ploidy and heat-shock metadata contain unsupported values.", call. = FALSE)
  }

  treatment <- ifelse(metadata$ploidy == "triploid", 1L, 0L)
  heat_shock <- ifelse(metadata$heat_shock == "yes", 1L, 0L)
  design <- stats::model.matrix(~ treatment + heat_shock)
  if (qr(design)$rank != ncol(design)) {
    stop("The primary treatment/covariate design matrix is rank deficient.", call. = FALSE)
  }
  if (any(table(treatment) < 2L)) {
    stop("Both ploidy groups must retain biological replication.", call. = FALSE)
  }

  cov_names <- vapply(
    metadata$seq_id,
    function(seq_id) gsub("{seq_id}", seq_id, cov_template, fixed = TRUE),
    character(1)
  )
  if (anyDuplicated(cov_names)) {
    stop("The COV filename template does not produce unique sample files.", call. = FALSE)
  }
  cov_paths <- file.path(cov_directory, cov_names)
  missing_files <- cov_names[!file.exists(cov_paths)]
  if (length(missing_files)) {
    stop(
      sprintf("Missing required COV files: %s", paste(missing_files, collapse = ", ")),
      call. = FALSE
    )
  }
  file_details <- file.info(cov_paths)
  if (any(file_details$isdir) || anyNA(file_details$size) || any(file_details$size <= 0)) {
    stop("Every required COV input must be a non-empty regular file.", call. = FALSE)
  }
  invisible(Map(validate_cov_preview, cov_paths, metadata$seq_id))

  input_manifest <- data.frame(
    seq_id = metadata$seq_id,
    library_name = metadata$library_name,
    ploidy = metadata$ploidy,
    heat_shock = metadata$heat_shock,
    cov_file = cov_names,
    bytes = file_details$size,
    stringsAsFactors = FALSE
  )
  write.csv(input_manifest, file.path(output_directory, "cov_input_manifest.csv"), row.names = FALSE)
  write.csv(metadata, file.path(output_directory, "metadata_snapshot.csv"), row.names = FALSE)

  parameter_manifest <- data.frame(
    parameter = c(
      "run_id", "assembly", "input_coordinates", "bed_coordinates", "treatment_0", "treatment_1",
      "covariate", "meth_diff_direction", "min_read_coverage", "min_analysis_coverage",
      "high_coverage_percentile", "normalize_coverage", "require_all_samples",
      "overdispersion", "test", "adjust", "effect", "maximum_q_value",
      "minimum_absolute_differences", "chunk_size", "mc_cores"
    ),
    value = c(
      run_id, args[["assembly"]], "zero-based half-open", "zero-based half-open", "diploid", "triploid",
      "heat_shock (no=0, yes=1)", "positive means triploid > diploid, adjusted for heat shock",
      min_read_coverage, min_analysis_coverage, high_coverage_percentile,
      normalize_coverage, require_all_samples, args[["overdispersion"]], args[["test"]],
      args[["adjust"]], args[["effect"]], max_q, paste(differences, collapse = ","),
      chunk_size, mc_cores
    ),
    stringsAsFactors = FALSE
  )
  write.csv(parameter_manifest, file.path(output_directory, "analysis_parameters.csv"), row.names = FALSE)

  message(sprintf("COV preflight passed for %d canonical samples.", nrow(metadata)))
  if (preflight_only) {
    message(sprintf("Preflight outputs written to %s", output_directory))
    return(invisible(0L))
  }

  if (as.character(getRversion()) != args[["r-version"]]) {
    stop(
      sprintf("Expected R %s; observed R %s.", args[["r-version"]], getRversion()),
      call. = FALSE
    )
  }
  if (!requireNamespace("BiocManager", quietly = TRUE) ||
      as.character(BiocManager::version()) != args[["bioconductor-version"]]) {
    stop(sprintf("Expected Bioconductor %s.", args[["bioconductor-version"]]), call. = FALSE)
  }
  if (!requireNamespace("methylKit", quietly = TRUE)) {
    stop("The pinned methylKit package is not installed.", call. = FALSE)
  }
  if (as.character(utils::packageVersion("methylKit")) != args[["methylkit-version"]]) {
    stop(
      sprintf(
        "Expected methylKit %s; observed %s.",
        args[["methylkit-version"]], utils::packageVersion("methylKit")
      ),
      call. = FALSE
    )
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("methylKit's data.table dependency is unavailable.", call. = FALSE)
  }

  database_directory <- file.path(output_directory, "methylkit_db")
  database_subdirectories <- file.path(
    database_directory,
    c("raw", "filtered", "normalized", "united", "differential")
  )
  invisible(lapply(
    database_subdirectories,
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  ))
  processed <- methylKit::methRead(
    as.list(cov_paths),
    sample.id = as.list(metadata$seq_id),
    assembly = args[["assembly"]],
    treatment = treatment,
    pipeline = "bismarkCoverage",
    header = FALSE,
    context = "CpG",
    resolution = "base",
    mincov = min_read_coverage,
    dbtype = "tabix",
    dbdir = file.path(database_directory, "raw")
  )

  filtered <- methylKit::filterByCoverage(
    processed,
    lo.count = min_analysis_coverage,
    lo.perc = NULL,
    hi.count = NULL,
    hi.perc = high_coverage_percentile,
    chunk.size = chunk_size,
    save.db = TRUE,
    suffix = "cov_filtered",
    dbdir = file.path(database_directory, "filtered")
  )
  if (normalize_coverage) {
    filtered <- methylKit::normalizeCoverage(
      filtered,
      method = "median",
      chunk.size = chunk_size,
      save.db = TRUE,
      suffix = "normalized",
      dbdir = file.path(database_directory, "normalized")
    )
  }

  united <- methylKit::unite(
    filtered,
    destrand = FALSE,
    min.per.group = if (require_all_samples) NULL else 2L,
    chunk.size = chunk_size,
    mc.cores = mc_cores,
    save.db = TRUE,
    suffix = "all_samples",
    dbdir = file.path(database_directory, "united")
  )
  if (!identical(methylKit::getSampleID(united), metadata$seq_id)) {
    stop("methylKit sample order changed before model fitting.", call. = FALSE)
  }

  differential <- methylKit::calculateDiffMeth(
    united,
    covariates = data.frame(heat_shock = heat_shock),
    overdispersion = args[["overdispersion"]],
    adjust = args[["adjust"]],
    effect = args[["effect"]],
    test = args[["test"]],
    mc.cores = mc_cores,
    chunk.size = chunk_size,
    save.db = TRUE,
    suffix = "ploidy_adjusted_heat",
    dbdir = file.path(database_directory, "differential")
  )

  all_stats <- data.table::as.data.table(methylKit::getData(differential))
  required_stats <- c("chr", "start", "end", "pvalue", "qvalue", "meth.diff")
  if (!all(required_stats %in% names(all_stats))) {
    stop("The complete methylKit statistics table lacks required fields.", call. = FALSE)
  }
  if (anyDuplicated(all_stats[, c("chr", "start", "end"), with = FALSE])) {
    stop("The complete statistics table contains duplicate genomic coordinates.", call. = FALSE)
  }
  if (any(!is.finite(all_stats$meth.diff)) || any(abs(all_stats$meth.diff) > 100)) {
    stop("The complete statistics table contains invalid methylation differences.", call. = FALSE)
  }
  data.table::setorder(all_stats, chr, start, end)
  model_tag <- paste0(
    "ploidy-adjusted-heat_",
    tolower(args[["overdispersion"]]), "-", tolower(args[["test"]]), "_",
    tolower(args[["adjust"]]), "_cov", min_analysis_coverage
  )
  run_tag <- paste(model_tag, run_id, sep = "_")
  all_stats_path <- file.path(output_directory, paste0("dml_stats_", run_tag, "_all.csv"))
  data.table::fwrite(all_stats, all_stats_path)

  output_counts <- data.frame(
    subset = "all_tested_cpgs",
    rows = nrow(all_stats),
    stringsAsFactors = FALSE
  )
  for (difference in differences) {
    selected <- !is.na(all_stats$qvalue) & all_stats$qvalue <= max_q &
      abs(all_stats$meth.diff) >= difference
    subset <- all_stats[selected]
    if (nrow(subset) &&
        (any(subset$qvalue > max_q) || any(abs(subset$meth.diff) < difference))) {
      stop(sprintf("The %s%% DML subset violates its declared thresholds.", difference), call. = FALSE)
    }
    stem <- sprintf(
      "dml_%s_diff%s_q%s",
      run_tag,
      format(difference, trim = TRUE, scientific = FALSE),
      format(max_q, trim = TRUE, scientific = FALSE)
    )
    csv_path <- file.path(output_directory, paste0(stem, ".csv"))
    bed_path <- file.path(output_directory, paste0(stem, ".bed"))
    data.table::fwrite(subset, csv_path)
    write_bed(subset, bed_path)
    output_counts <- rbind(
      output_counts,
      data.frame(subset = stem, rows = nrow(subset), stringsAsFactors = FALSE)
    )
  }
  write.csv(output_counts, file.path(output_directory, "output_counts.csv"), row.names = FALSE)
  writeLines(capture.output(utils::sessionInfo()), file.path(output_directory, "session_info.txt"))
  message(sprintf("Clean MethylKit DML run completed in %s", output_directory))
  invisible(0L)
}

main()
