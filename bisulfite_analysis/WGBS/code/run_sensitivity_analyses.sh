#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

# Phase 3 Sensitivity Analysis Framework
# Tests DML robustness across coverage thresholds and p-value adjustment methods

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../../.." && pwd)"
config_path="${repo_root}/bisulfite_analysis/WGBS/config/analysis_config.json"

cov_environment_variable="$(jq -er '.inputs.cov_files.local_directory_environment_variable' "${config_path}")"
cov_directory="${!cov_environment_variable:-}"
if [[ -z "${cov_directory}" ]]; then
  echo "Error: Set ${cov_environment_variable} to the COV directory." >&2
  exit 1
fi

baseline_results_dir="${repo_root}/bisulfite_analysis/WGBS/results/dml"
sensitivity_dir="${repo_root}/bisulfite_analysis/WGBS/results/sensitivity_analyses"
mkdir -p "${sensitivity_dir}"

# Load baseline results directory from primary run
baseline_run_dir=$(ls -dt "${baseline_results_dir}"/dml/* 2>/dev/null | head -1)
if [[ -z "${baseline_run_dir}" ]]; then
  echo "Error: No baseline DML results found. Run baseline analysis first." >&2
  exit 1
fi

baseline_dml="${baseline_run_dir}/filtered_dml_20pct_effect.txt"
if [[ ! -f "${baseline_dml}" ]]; then
  echo "Error: Baseline DML file not found: ${baseline_dml}" >&2
  exit 1
fi

echo "Baseline run directory: ${baseline_run_dir}"
echo "Baseline DML rows: $(wc -l < "${baseline_dml}")"
echo ""

# Test 1: Coverage threshold sensitivity
echo "=== Testing Coverage Threshold Sensitivity ==="
coverage_results_dir="${sensitivity_dir}/coverage_threshold"
mkdir -p "${coverage_results_dir}"

coverage_thresholds=(5 10 15 20)
for cov_thresh in "${coverage_thresholds[@]}"; do
  echo "Testing min-analysis-coverage=${cov_thresh}..."
  output_dir="${repo_root}/bisulfite_analysis/WGBS/results/dml/sensitivity_coverage_${cov_thresh}"
  
  Rscript --vanilla "${script_dir}/2_WGBS_Methylkit.R" \
    "--metadata=${repo_root}/$(jq -er '.sample_metadata' "${config_path}")" \
    "--cov-dir=${cov_directory}" \
    "--cov-template=$(jq -er '.inputs.cov_files.file_name_template' "${config_path}")" \
    "--output-dir=${output_dir}" \
    "--assembly=$(jq -er '.reference.assembly' "${config_path}")" \
    "--min-read-coverage=$(jq -er '.filtering.methread_minimum_coverage' "${config_path}")" \
    "--min-analysis-coverage=${cov_thresh}" \
    "--high-coverage-percentile=$(jq -er '.filtering.high_coverage_percentile' "${config_path}")" \
    "--normalize-coverage=$(jq -er '.filtering.normalize_coverage' "${config_path}")" \
    "--require-all-samples=$(jq -er '.filtering.require_cpg_in_all_samples' "${config_path}")" \
    "--overdispersion=$(jq -er '.primary_model.overdispersion' "${config_path}")" \
    "--test=$(jq -er '.primary_model.test' "${config_path}")" \
    "--adjust=$(jq -er '.primary_model.p_value_adjustment' "${config_path}")" \
    "--effect=$(jq -er '.primary_model.effect' "${config_path}")" \
    "--max-q=$(jq -er '.primary_model.maximum_q_value' "${config_path}")" \
    "--differences=$(jq -er '.primary_model.minimum_absolute_methylation_differences | join(",")' "${config_path}")" \
    "--r-version=$(jq -er '.runtime.r_version' "${config_path}")" \
    "--bioconductor-version=$(jq -er '.runtime.bioconductor_version' "${config_path}")" \
    "--methylkit-version=$(jq -er '.runtime.methylkit_version' "${config_path}")" \
    "--chunk-size=$(jq -er '.runtime.chunk_size' "${config_path}")" \
    "--mc-cores=$(jq -er '.runtime.mc_cores' "${config_path}")" \
    "--preflight-only=false" 2>&1 | tee "${coverage_results_dir}/cov_${cov_thresh}.log"
  
  # Count results
  result_file="${output_dir}/filtered_dml_20pct_effect.txt"
  if [[ -f "${result_file}" ]]; then
    row_count=$(wc -l < "${result_file}")
    echo "  cov=${cov_thresh}: ${row_count} DML sites" >> "${coverage_results_dir}/summary.txt"
  fi
done

# Test 2: P-value adjustment method sensitivity
echo ""
echo "=== Testing Multiple-Testing Correction Sensitivity ==="
correction_results_dir="${sensitivity_dir}/correction_methods"
mkdir -p "${correction_results_dir}"

correction_methods=("SLIM" "BH" "bonferroni")
for method in "${correction_methods[@]}"; do
  echo "Testing adjustment method=${method}..."
  output_dir="${repo_root}/bisulfite_analysis/WGBS/results/dml/sensitivity_correction_${method}"
  
  Rscript --vanilla "${script_dir}/2_WGBS_Methylkit.R" \
    "--metadata=${repo_root}/$(jq -er '.sample_metadata' "${config_path}")" \
    "--cov-dir=${cov_directory}" \
    "--cov-template=$(jq -er '.inputs.cov_files.file_name_template' "${config_path}")" \
    "--output-dir=${output_dir}" \
    "--assembly=$(jq -er '.reference.assembly' "${config_path}")" \
    "--min-read-coverage=$(jq -er '.filtering.methread_minimum_coverage' "${config_path}")" \
    "--min-analysis-coverage=$(jq -er '.filtering.analysis_minimum_coverage' "${config_path}")" \
    "--high-coverage-percentile=$(jq -er '.filtering.high_coverage_percentile' "${config_path}")" \
    "--normalize-coverage=$(jq -er '.filtering.normalize_coverage' "${config_path}")" \
    "--require-all-samples=$(jq -er '.filtering.require_cpg_in_all_samples' "${config_path}")" \
    "--overdispersion=$(jq -er '.primary_model.overdispersion' "${config_path}")" \
    "--test=$(jq -er '.primary_model.test' "${config_path}")" \
    "--adjust=${method}" \
    "--effect=$(jq -er '.primary_model.effect' "${config_path}")" \
    "--max-q=$(jq -er '.primary_model.maximum_q_value' "${config_path}")" \
    "--differences=$(jq -er '.primary_model.minimum_absolute_methylation_differences | join(",")' "${config_path}")" \
    "--r-version=$(jq -er '.runtime.r_version' "${config_path}")" \
    "--bioconductor-version=$(jq -er '.runtime.bioconductor_version' "${config_path}")" \
    "--methylkit-version=$(jq -er '.runtime.methylkit_version' "${config_path}")" \
    "--chunk-size=$(jq -er '.runtime.chunk_size' "${config_path}")" \
    "--mc-cores=$(jq -er '.runtime.mc_cores' "${config_path}")" \
    "--preflight-only=false" 2>&1 | tee "${correction_results_dir}/${method}.log"
  
  # Count results
  result_file="${output_dir}/filtered_dml_20pct_effect.txt"
  if [[ -f "${result_file}" ]]; then
    row_count=$(wc -l < "${result_file}")
    echo "  adjust=${method}: ${row_count} DML sites" >> "${correction_results_dir}/summary.txt"
  fi
done

echo ""
echo "=== Sensitivity Analysis Complete ==="
echo "Results: ${sensitivity_dir}/"
