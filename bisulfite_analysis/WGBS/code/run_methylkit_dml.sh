#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../../.." && pwd)"
config_path="${repo_root}/bisulfite_analysis/WGBS/config/analysis_config.json"

preflight_only=false
if [[ "${1:-}" = "--preflight-only" ]]; then
  preflight_only=true
  shift
fi
if [[ "$#" -gt 1 ]]; then
  echo "Usage: $0 [--preflight-only] [output-directory]" >&2
  exit 2
fi

cov_environment_variable="$(jq -er '.inputs.cov_files.local_directory_environment_variable' "${config_path}")"
cov_directory="${!cov_environment_variable:-}"
if [[ -z "${cov_directory}" ]]; then
  echo "Set ${cov_environment_variable} to the directory containing the ten external COV files." >&2
  exit 1
fi
if [[ ! -d "${cov_directory}" ]]; then
  echo "COV directory does not exist: ${cov_directory}" >&2
  exit 1
fi

run_root_relative="$(jq -er '.outputs.run_directory' "${config_path}")"
if [[ "$#" -eq 1 ]]; then
  if [[ "$1" = /* ]]; then
    output_directory="$1"
  else
    output_directory="${repo_root}/$1"
  fi
else
  run_id="$(date -u +%Y%m%dT%H%M%SZ)"
  mode="dml"
  if [[ "${preflight_only}" = true ]]; then
    mode="dml_preflight"
  fi
  output_directory="${repo_root}/${run_root_relative}/${mode}/${run_id}"
fi

Rscript --vanilla "${script_dir}/2_WGBS_Methylkit.R" \
  "--metadata=${repo_root}/$(jq -er '.sample_metadata' "${config_path}")" \
  "--cov-dir=${cov_directory}" \
  "--cov-template=$(jq -er '.inputs.cov_files.file_name_template' "${config_path}")" \
  "--output-dir=${output_directory}" \
  "--assembly=$(jq -er '.reference.assembly' "${config_path}")" \
  "--min-read-coverage=$(jq -er '.filtering.methread_minimum_coverage' "${config_path}")" \
  "--min-analysis-coverage=$(jq -er '.filtering.analysis_minimum_coverage' "${config_path}")" \
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
  "--preflight-only=${preflight_only}"

if [[ "${preflight_only}" = false ]]; then
  Rscript --vanilla "${repo_root}/bisulfite_analysis/WGBS/tests/validate_dml_outputs.R" \
    "${output_directory}"
fi
