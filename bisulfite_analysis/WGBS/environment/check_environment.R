expected <- c(
  R = "4.6.0",
  Bioconductor = "3.23",
  methylKit = "1.38.0"
)

observed <- c(
  R = as.character(getRversion()),
  Bioconductor = if (requireNamespace("BiocManager", quietly = TRUE)) {
    as.character(BiocManager::version())
  } else {
    NA_character_
  },
  methylKit = if (requireNamespace("methylKit", quietly = TRUE)) {
    as.character(utils::packageVersion("methylKit"))
  } else {
    NA_character_
  }
)

invalid <- is.na(observed) | observed != expected
if (any(invalid)) {
  details <- paste0(names(expected), ": expected ", expected, ", observed ", observed)
  stop(
    paste(c("The WGBS analysis environment is invalid:", details[invalid]), collapse = "\n"),
    call. = FALSE
  )
}

message(
  sprintf(
    "Environment validated: R %s, Bioconductor %s, methylKit %s.",
    observed[["R"]], observed[["Bioconductor"]], observed[["methylKit"]]
  )
)
