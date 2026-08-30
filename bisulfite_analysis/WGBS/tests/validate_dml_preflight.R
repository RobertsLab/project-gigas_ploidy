arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L) {
  stop("Usage: validate_dml_preflight.R OUTPUT_DIRECTORY", call. = FALSE)
}

output_directory <- normalizePath(arguments[[1]], mustWork = TRUE)
manifest <- read.csv(file.path(output_directory, "cov_input_manifest.csv"), stringsAsFactors = FALSE)
metadata <- read.csv(file.path(output_directory, "metadata_snapshot.csv"), stringsAsFactors = FALSE)
parameters <- read.csv(file.path(output_directory, "analysis_parameters.csv"), stringsAsFactors = FALSE)

expected_ids <- paste0("zr3534_", seq_len(10))
if (!identical(manifest$seq_id, expected_ids) || anyDuplicated(manifest$cov_file)) {
  stop("DML preflight did not resolve each canonical sample exactly once.", call. = FALSE)
}
if (!identical(metadata$seq_id, expected_ids) || nrow(metadata) != 10L) {
  stop("DML preflight metadata snapshot is invalid.", call. = FALSE)
}
expected_parameters <- c(
  input_coordinates = "zero-based half-open",
  bed_coordinates = "zero-based half-open",
  treatment_0 = "diploid",
  treatment_1 = "triploid",
  covariate = "heat_shock (no=0, yes=1)",
  overdispersion = "MN",
  test = "F",
  adjust = "SLIM",
  effect = "predicted"
)
observed_parameters <- setNames(parameters$value, parameters$parameter)
if (!identical(unname(observed_parameters[names(expected_parameters)]), unname(expected_parameters))) {
  stop("DML preflight parameters do not match the declared primary model.", call. = FALSE)
}
if (!all(manifest$ploidy == metadata$ploidy) || !all(manifest$heat_shock == metadata$heat_shock)) {
  stop("DML preflight treatment assignments differ from canonical metadata.", call. = FALSE)
}

message("DML preflight validation passed for 10 COV fixtures and the declared primary model.")
