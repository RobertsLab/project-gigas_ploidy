arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) {
  stop("Usage: validate_pe_summary.R OUTPUT_DIRECTORY REPORT_DIRECTORY", call. = FALSE)
}

output_directory <- normalizePath(arguments[[1]], mustWork = TRUE)
report_directory <- normalizePath(arguments[[2]], mustWork = TRUE)
per_sample <- read.csv(file.path(output_directory, "per_sample.csv"), stringsAsFactors = FALSE)
group_summary <- read.csv(file.path(output_directory, "group_summary.csv"), stringsAsFactors = FALSE)

expected_ids <- paste0("zr3534_", seq_len(10))
if (nrow(per_sample) != 10L || anyDuplicated(per_sample$seq_id) || !identical(per_sample$seq_id, expected_ids)) {
  stop("Generated per-sample output does not contain each expected sequence ID exactly once.", call. = FALSE)
}

expected_cpg <- c(10.5, 10.4, 11.8, 11.2, 11.2, 10.3, 10.1, 10.3, 10.0, 10.0)
if (!isTRUE(all.equal(per_sample$cpg_percent, expected_cpg, tolerance = 1e-12))) {
  stop("Generated CpG percentages differ from the independently expected values.", call. = FALSE)
}
if (!all(per_sample$bismark_version == "v0.21.0")) {
  stop("The committed reports must identify Bismark v0.21.0.", call. = FALSE)
}

report_files <- list.files(
  report_directory,
  pattern = "_bismark_bt2_PE_report[.]txt$",
  full.names = TRUE
)
direct_values <- vapply(report_files, function(report_file) {
  line <- grep("^C methylated in CpG context:", readLines(report_file, warn = FALSE), value = TRUE)
  if (length(line) != 1L) {
    stop(sprintf("Expected one CpG percentage in %s.", basename(report_file)), call. = FALSE)
  }
  as.numeric(sub("%$", "", trimws(sub("^[^:]+:", "", line))))
}, numeric(1))
names(direct_values) <- sub("_R1.*$", "", basename(report_files))
if (!isTRUE(all.equal(per_sample$cpg_percent, unname(direct_values[per_sample$seq_id]), tolerance = 1e-12))) {
  stop("Generated CpG percentages do not match direct extraction from reports.", call. = FALSE)
}

if (nrow(group_summary) != 4L || anyDuplicated(group_summary[c("ploidy", "heat_shock")])) {
  stop("Generated group summary must contain the four unique ploidy/heat-shock groups.", call. = FALSE)
}
for (index in seq_len(nrow(group_summary))) {
  row <- group_summary[index, ]
  selected <- per_sample$ploidy == row$ploidy & per_sample$heat_shock == row$heat_shock
  values <- per_sample$cpg_percent[selected]
  expected <- c(
    n = length(values),
    mean_cpg_percent = mean(values),
    sd_cpg_percent = stats::sd(values),
    se_cpg_percent = stats::sd(values) / sqrt(length(values))
  )
  observed <- unlist(row[names(expected)], use.names = TRUE)
  if (!isTRUE(all.equal(observed, expected, tolerance = 1e-12, check.attributes = FALSE))) {
    stop(sprintf("Group statistics are invalid for %s/%s.", row$ploidy, row$heat_shock), call. = FALSE)
  }
}

message("PE summary validation passed for 10 samples and 4 treatment groups.")
