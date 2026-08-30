#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

# Phase 2 Acceptance Testing Script
# Validates that WGBS preprocessing and report summaries meet all acceptance criteria

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../../.." && pwd)"

echo "=========================================="
echo "Phase 2: WGBS Preprocessing Reproducibility"
echo "=========================================="
echo ""

# Test 1: Run PE report summary generation
echo "Test 1: Summary table generation from raw reports"
echo "---"
test_output_dir="/tmp/phase2_test_output_$$"
mkdir -p "${test_output_dir}"

cd "${repo_root}"
bash bisulfite_analysis/WGBS/code/run_pe_report_summary.sh "${test_output_dir}" > /dev/null 2>&1

if [[ -f "${test_output_dir}/per_sample.csv" && -f "${test_output_dir}/group_summary.csv" ]]; then
  echo "✅ Summary tables generated successfully"
  echo "   - per_sample.csv: $(wc -l < "${test_output_dir}/per_sample.csv" | tr -d ' ') rows"
  echo "   - group_summary.csv: $(wc -l < "${test_output_dir}/group_summary.csv" | tr -d ' ') rows"
else
  echo "❌ FAILED: Summary tables not generated"
  exit 1
fi
echo ""

# Test 2: Verify acceptance criteria
echo "Test 2: Acceptance criteria validation"
echo "---"

# Run R validation script
Rscript - "${test_output_dir}" << 'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
test_output_dir <- args[[1]]
repo_root <- normalizePath(getwd())

per_sample <- read.csv(file.path(test_output_dir, "per_sample.csv"), stringsAsFactors = FALSE)
group_summary <- read.csv(file.path(test_output_dir, "group_summary.csv"), stringsAsFactors = FALSE)
metadata <- read.csv(file.path(repo_root, "metadata/wgbs_samples.csv"), stringsAsFactors = FALSE)

# Criterion 1: Summary generated from raw reports
if (nrow(per_sample) == 10 && nrow(group_summary) == 4) {
  cat("✅ Criterion 1: Summary tables contain expected rows\n")
} else {
  cat("❌ FAILED Criterion 1\n")
  quit(status = 1)
}

# Criterion 2: Group summaries equal recalculation
all_match <- TRUE
for (i in seq_len(nrow(group_summary))) {
  row <- group_summary[i, ]
  sel <- per_sample$ploidy == row$ploidy & per_sample$heat_shock == row$heat_shock
  values <- per_sample$cpg_percent[sel]
  
  expected_n <- length(values)
  expected_mean <- mean(values)
  expected_sd <- stats::sd(values)
  expected_se <- expected_sd / sqrt(expected_n)
  
  n_match <- row$n == expected_n
  mean_match <- isTRUE(all.equal(row$mean_cpg_percent, expected_mean, tolerance = 1e-10))
  sd_match <- isTRUE(all.equal(row$sd_cpg_percent, expected_sd, tolerance = 1e-10))
  se_match <- isTRUE(all.equal(row$se_cpg_percent, expected_se, tolerance = 1e-10))
  
  if (!n_match || !mean_match || !sd_match || !se_match) {
    all_match <- FALSE
  }
}

if (all_match) {
  cat("✅ Criterion 2: Group summaries match recalculation from per-sample\n")
} else {
  cat("❌ FAILED Criterion 2\n")
  quit(status = 1)
}

# Criterion 3: Sample names and treatment assignments match metadata
if (identical(per_sample$seq_id, metadata$seq_id) &&
    identical(per_sample$library_name, metadata$library_name) &&
    identical(per_sample$ploidy, metadata$ploidy) &&
    identical(per_sample$heat_shock, metadata$heat_shock)) {
  cat("✅ Criterion 3: Sample names and treatment assignments match metadata\n")
} else {
  cat("❌ FAILED Criterion 3\n")
  quit(status = 1)
}
RSCRIPT

if [[ $? -ne 0 ]]; then
  exit 1
fi
echo ""

# Test 3: Documentation completeness
echo "Test 3: Documentation completeness"
echo "---"

if [[ -f "${repo_root}/bisulfite_analysis/WGBS/PREPROCESSING.md" ]]; then
  echo "✅ Preprocessing documentation exists (PREPROCESSING.md)"
  if grep -q "Bismark version" "${repo_root}/bisulfite_analysis/WGBS/PREPROCESSING.md" && \
     grep -q "Bowtie 2" "${repo_root}/bisulfite_analysis/WGBS/PREPROCESSING.md" && \
     grep -q "zero-based\|half-open" "${repo_root}/bisulfite_analysis/WGBS/PREPROCESSING.md"; then
    echo "✅ Documentation includes Bismark versions, tools, and coordinate system"
  else
    echo "❌ Documentation missing required content"
    exit 1
  fi
else
  echo "❌ PREPROCESSING.md not found"
  exit 1
fi
echo ""

# Test 4: PE report documentation
echo "Test 4: PE report documentation"
echo "---"

if [[ -f "${repo_root}/bisulfite_analysis/WGBS/data/PE_reports/README.md" ]]; then
  echo "✅ PE report documentation exists"
  if grep -q "zr3534_7\|recovered\|Recovered" "${repo_root}/bisulfite_analysis/WGBS/data/PE_reports/README.md" && \
     grep -q "e858d1b34f890323f1d8370214512c10ff5a595d326b497838a3cf264cbec866" "${repo_root}/bisulfite_analysis/WGBS/data/PE_reports/README.md"; then
    echo "✅ Documentation includes zr3534_7 recovery information and checksum"
  else
    echo "❌ Documentation missing zr3534_7 recovery details"
    exit 1
  fi
else
  echo "❌ PE reports README not found"
  exit 1
fi
echo ""

# Test 5: Verify Phase 1 still passes
echo "Test 5: Phase 1 validation remains intact"
echo "---"

if bash "${repo_root}/bisulfite_analysis/WGBS/tests/run_phase1_checks.sh" > /dev/null 2>&1; then
  echo "✅ Phase 1 checks still pass"
else
  echo "❌ Phase 1 checks failed"
  exit 1
fi
echo ""

# Cleanup
rm -rf "${test_output_dir}"

echo "=========================================="
echo "✅ Phase 2 acceptance tests PASSED"
echo "=========================================="
echo ""
echo "Summary:"
echo "  ✅ PE report summaries generated reproducibly"
echo "  ✅ Group statistics match recalculation"
echo "  ✅ Sample metadata correctly joined"
echo "  ✅ Preprocessing documented with tool versions"
echo "  ✅ Coordinate systems documented"
echo "  ✅ Missing zr3534_7 report recovered and checksummed"
echo "  ✅ Phase 1 validation preserved"
