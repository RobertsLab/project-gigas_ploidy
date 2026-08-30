#!/usr/bin/env Rscript

options(scipen = 999)

parse_arguments <- function(arguments) {
  values <- list()
  for (argument in arguments) {
    if (!grepl("^--[^=]+=.+$", argument)) {
      stop(sprintf("Invalid argument '%s'; expected --name=value.", argument), call. = FALSE)
    }
    parts <- strsplit(sub("^--", "", argument), "=", fixed = TRUE)[[1]]
    key <- parts[[1]]
    value <- paste(parts[-1], collapse = "=")
    if (!key %in% c("metadata", "reports", "output-dir")) {
      stop(sprintf("Unknown argument '--%s'.", key), call. = FALSE)
    }
    if (!is.null(values[[key]])) {
      stop(sprintf("Argument '--%s' was supplied more than once.", key), call. = FALSE)
    }
    values[[key]] <- value
  }

  missing <- setdiff(c("metadata", "reports", "output-dir"), names(values))
  if (length(missing)) {
    stop(sprintf("Missing arguments: %s", paste(paste0("--", missing), collapse = ", ")), call. = FALSE)
  }
  values
}

extract_number <- function(lines, label, report_file, percent = FALSE) {
  matches <- lines[startsWith(lines, paste0(label, ":"))]
  if (length(matches) != 1L) {
    stop(
      sprintf("Expected one '%s' line in %s; found %d.", label, report_file, length(matches)),
      call. = FALSE
    )
  }
  value <- trimws(sub("^[^:]+:", "", matches))
  if (percent) {
    value <- sub("%$", "", value)
  }
  result <- suppressWarnings(as.numeric(value))
  if (length(result) != 1L || is.na(result)) {
    stop(sprintf("Could not parse '%s' in %s.", label, report_file), call. = FALSE)
  }
  result
}

parse_report <- function(report_path) {
  report_file <- basename(report_path)
  lines <- readLines(report_path, warn = FALSE)
  seq_id <- sub("_R1.*$", "", report_file)

  header <- lines[startsWith(lines, "Bismark report for:")]
  if (length(header) != 1L || !grepl("\\(version: [^)]+\\)$", header)) {
    stop(sprintf("Could not parse the Bismark header in %s.", report_file), call. = FALSE)
  }
  bismark_version <- sub(".*\\(version: ([^)]+)\\)$", "\\1", header)

  labels <- c(
    sequence_pairs_total = "Sequence pairs analysed in total",
    unique_best_hit = "Number of paired-end alignments with a unique best hit",
    mapping_efficiency_percent = "Mapping efficiency",
    total_cytosines = "Total number of C's analysed",
    methylated_cpg = "Total methylated C's in CpG context",
    methylated_chg = "Total methylated C's in CHG context",
    methylated_chh = "Total methylated C's in CHH context",
    methylated_unknown = "Total methylated C's in Unknown context",
    unmethylated_cpg = "Total unmethylated C's in CpG context",
    unmethylated_chg = "Total unmethylated C's in CHG context",
    unmethylated_chh = "Total unmethylated C's in CHH context",
    unmethylated_unknown = "Total unmethylated C's in Unknown context",
    cpg_percent = "C methylated in CpG context",
    chg_percent = "C methylated in CHG context",
    chh_percent = "C methylated in CHH context",
    unknown_percent = "C methylated in unknown context (CN or CHN)"
  )
  percent_fields <- c(
    "mapping_efficiency_percent", "cpg_percent", "chg_percent", "chh_percent", "unknown_percent"
  )
  values <- vapply(
    names(labels),
    function(field) extract_number(lines, labels[[field]], report_file, field %in% percent_fields),
    numeric(1)
  )

  data.frame(
    seq_id = seq_id,
    report_file = report_file,
    bismark_version = bismark_version,
    as.list(values),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

assert_report_percentages <- function(report_row) {
  known_context_total <- sum(unlist(report_row[c(
    "methylated_cpg", "methylated_chg", "methylated_chh",
    "unmethylated_cpg", "unmethylated_chg", "unmethylated_chh"
  )]))
  if (!identical(report_row$total_cytosines, known_context_total)) {
    stop(
      sprintf(
        "%s reports %s total cytosines but its CpG/CHG/CHH counts sum to %s.",
        report_row$report_file, report_row$total_cytosines, known_context_total
      ),
      call. = FALSE
    )
  }

  checks <- list(
    mapping_efficiency_percent = c("unique_best_hit", "sequence_pairs_total"),
    cpg_percent = c("methylated_cpg", "unmethylated_cpg"),
    chg_percent = c("methylated_chg", "unmethylated_chg"),
    chh_percent = c("methylated_chh", "unmethylated_chh"),
    unknown_percent = c("methylated_unknown", "unmethylated_unknown")
  )
  for (reported_field in names(checks)) {
    numerator <- report_row[[checks[[reported_field]][[1]]]]
    denominator_component <- report_row[[checks[[reported_field]][[2]]]]
    denominator <- if (reported_field == "mapping_efficiency_percent") {
      denominator_component
    } else {
      numerator + denominator_component
    }
    recalculated <- round(100 * numerator / denominator, 1)
    if (!isTRUE(all.equal(report_row[[reported_field]], recalculated, tolerance = 1e-12))) {
      stop(
        sprintf(
          "%s reports %s=%0.1f but its counts give %0.1f.",
          report_row$report_file, reported_field, report_row[[reported_field]], recalculated
        ),
        call. = FALSE
      )
    }
  }
}

arguments <- parse_arguments(commandArgs(trailingOnly = TRUE))
metadata_path <- normalizePath(arguments[["metadata"]], mustWork = TRUE)
report_directory <- normalizePath(arguments[["reports"]], mustWork = TRUE)
output_directory <- arguments[["output-dir"]]

metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
required_metadata <- c(
  "seq_id", "library_name", "tissue", "ploidy", "desiccation", "heat_shock",
  "library_kit", "sra_bioproject", "sra_accession"
)
if (!identical(names(metadata), required_metadata)) {
  stop("Canonical WGBS metadata columns or column order are invalid.", call. = FALSE)
}
if (anyDuplicated(metadata$seq_id) || anyNA(metadata$seq_id) || any(!nzchar(metadata$seq_id))) {
  stop("Canonical metadata seq_id values must be complete and unique.", call. = FALSE)
}

report_paths <- sort(list.files(
  report_directory,
  pattern = "_bismark_bt2_PE_report[.]txt$",
  full.names = TRUE
))
if (!length(report_paths)) {
  stop(sprintf("No Bismark PE reports found in %s.", report_directory), call. = FALSE)
}
parsed <- do.call(rbind, lapply(report_paths, parse_report))
if (anyDuplicated(parsed$seq_id)) {
  stop("Each sequence ID must have exactly one Bismark PE report.", call. = FALSE)
}

missing_reports <- setdiff(metadata$seq_id, parsed$seq_id)
unexpected_reports <- setdiff(parsed$seq_id, metadata$seq_id)
if (length(missing_reports) || length(unexpected_reports)) {
  stop(
    sprintf(
      "Report/metadata mismatch. Missing reports: [%s]; unexpected reports: [%s].",
      paste(missing_reports, collapse = ", "), paste(unexpected_reports, collapse = ", ")
    ),
    call. = FALSE
  )
}

invisible(lapply(seq_len(nrow(parsed)), function(index) assert_report_percentages(parsed[index, ])))
parsed <- parsed[match(metadata$seq_id, parsed$seq_id), , drop = FALSE]
per_sample <- cbind(metadata, parsed[setdiff(names(parsed), "seq_id")])

group_keys <- unique(per_sample[c("ploidy", "heat_shock")])
group_summary <- do.call(rbind, lapply(seq_len(nrow(group_keys)), function(index) {
  key <- group_keys[index, , drop = FALSE]
  selected <- per_sample$ploidy == key$ploidy & per_sample$heat_shock == key$heat_shock
  values <- per_sample$cpg_percent[selected]
  standard_deviation <- stats::sd(values)
  data.frame(
    ploidy = key$ploidy,
    heat_shock = key$heat_shock,
    n = length(values),
    mean_cpg_percent = mean(values),
    sd_cpg_percent = standard_deviation,
    se_cpg_percent = standard_deviation / sqrt(length(values)),
    stringsAsFactors = FALSE
  )
}))
rownames(group_summary) <- NULL

if (dir.exists(output_directory) && length(list.files(output_directory, all.files = TRUE, no.. = TRUE))) {
  stop(sprintf("Output directory already exists and is not empty: %s", output_directory), call. = FALSE)
}
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
output_directory <- normalizePath(output_directory, mustWork = TRUE)
write.csv(per_sample, file.path(output_directory, "per_sample.csv"), row.names = FALSE, na = "")
write.csv(group_summary, file.path(output_directory, "group_summary.csv"), row.names = FALSE, na = "")

message(sprintf("Parsed %d Bismark PE reports.", nrow(per_sample)))
message(sprintf("Wrote source-derived summaries to %s", output_directory))
