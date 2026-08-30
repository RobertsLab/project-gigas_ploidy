args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
if (length(script_arg) != 1L) {
  stop("Run this validation with Rscript.", call. = FALSE)
}

script_path <- normalizePath(sub("^--file=", "", script_arg))
repo_root <- normalizePath(file.path(dirname(script_path), "../../.."))

metadata_relative <- Sys.getenv("SAMPLE_METADATA_PATH")
if (!nzchar(metadata_relative)) {
  stop("SAMPLE_METADATA_PATH was not supplied by run_phase1_checks.sh.", call. = FALSE)
}

metadata_path <- file.path(repo_root, metadata_relative)
samples <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)

required_columns <- c(
  "seq_id", "library_name", "tissue", "ploidy", "desiccation",
  "heat_shock", "library_kit", "sra_bioproject", "sra_accession"
)
if (!identical(names(samples), required_columns)) {
  stop("Canonical metadata columns or column order are invalid.", call. = FALSE)
}

expected_ids <- paste0("zr3534_", seq_len(10))
if (!setequal(samples$seq_id, expected_ids) || nrow(samples) != 10L) {
  stop("Canonical metadata must contain exactly zr3534_1 through zr3534_10.", call. = FALSE)
}

unique_fields <- c("seq_id", "library_name", "sra_accession")
for (field in unique_fields) {
  if (anyDuplicated(samples[[field]]) || anyNA(samples[[field]]) || any(!nzchar(samples[[field]]))) {
    stop(sprintf("Field '%s' must be complete and unique.", field), call. = FALSE)
  }
}

expected_libraries <- c(
  zr3534_1 = "D11-C", zr3534_2 = "D12-C", zr3534_3 = "D13-C",
  zr3534_4 = "D19-C", zr3534_5 = "D20-C", zr3534_6 = "T11-C",
  zr3534_7 = "T12-C", zr3534_8 = "T13-C", zr3534_9 = "T19-C",
  zr3534_10 = "T20-C"
)
observed_libraries <- setNames(samples$library_name, samples$seq_id)
if (!identical(observed_libraries[names(expected_libraries)], expected_libraries)) {
  stop("The seq_id-to-library_name mapping is invalid.", call. = FALSE)
}

if (!all(samples$tissue == "ctenidia") || !all(samples$desiccation == "yes")) {
  stop("All WGBS samples must be desiccated ctenidia in this experiment.", call. = FALSE)
}

design <- with(samples, table(ploidy, heat_shock))
expected_design <- matrix(
  c(3L, 2L, 3L, 2L),
  nrow = 2L,
  byrow = TRUE,
  dimnames = list(c("diploid", "triploid"), c("no", "yes"))
)
if (!isTRUE(all.equal(unclass(design), expected_design, check.attributes = FALSE))) {
  stop("Expected three non-heat-shocked and two heat-shocked samples per ploidy.", call. = FALSE)
}

report_directory <- file.path(repo_root, "bisulfite_analysis/WGBS/data/PE_reports")
report_files <- list.files(
  report_directory,
  pattern = "_bismark_bt2_PE_report[.]txt$",
  full.names = FALSE
)
report_ids <- sub("_R1.*$", "", report_files)
if (anyDuplicated(report_ids)) {
  stop("More than one Bismark PE report was found for a sequence ID.", call. = FALSE)
}

allowed_missing_raw <- Sys.getenv("ALLOWED_MISSING_PE_REPORTS")
allowed_missing <- if (nzchar(allowed_missing_raw)) {
  strsplit(allowed_missing_raw, ",", fixed = TRUE)[[1]]
} else {
  character()
}
observed_missing <- setdiff(samples$seq_id, report_ids)
unexpected_reports <- setdiff(report_ids, samples$seq_id)
if (!setequal(observed_missing, allowed_missing)) {
  stop(
    sprintf(
      "PE report exceptions differ from configuration. Observed missing: %s",
      paste(observed_missing, collapse = ", ")
    ),
    call. = FALSE
  )
}
if (length(unexpected_reports)) {
  stop(
    sprintf("PE reports lack canonical metadata: %s", paste(unexpected_reports, collapse = ", ")),
    call. = FALSE
  )
}

message("Metadata validation passed for 10 WGBS samples.")
message("Configured missing PE reports: ", if (length(allowed_missing)) paste(allowed_missing, collapse = ", ") else "none")
