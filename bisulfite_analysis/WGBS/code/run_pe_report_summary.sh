#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../../.." && pwd)"
config_path="${repo_root}/bisulfite_analysis/WGBS/config/analysis_config.json"

if [[ "$#" -gt 1 ]]; then
  echo "Usage: $0 [output-directory]" >&2
  exit 2
fi

metadata_relative="$(jq -er '.sample_metadata' "${config_path}")"
reports_relative="$(jq -er '.inputs.pe_report_directory' "${config_path}")"
run_root_relative="$(jq -er '.outputs.run_directory' "${config_path}")"

if [[ "$#" -eq 1 ]]; then
  if [[ "$1" = /* ]]; then
    output_directory="$1"
  else
    output_directory="${repo_root}/$1"
  fi
else
  run_id="$(date -u +%Y%m%dT%H%M%SZ)"
  output_directory="${repo_root}/${run_root_relative}/pe_report_summary/${run_id}"
fi

Rscript "${script_dir}/1_WGBS_PE_report_analysis.R" \
  "--metadata=${repo_root}/${metadata_relative}" \
  "--reports=${repo_root}/${reports_relative}" \
  "--output-dir=${output_directory}"
