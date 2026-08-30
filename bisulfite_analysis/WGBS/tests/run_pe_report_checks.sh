#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../../.." && pwd)"
config_path="${repo_root}/bisulfite_analysis/WGBS/config/analysis_config.json"
report_directory="${repo_root}/$(jq -er '.inputs.pe_report_directory' "${config_path}")"

if [[ "$#" -gt 1 ]]; then
  echo "Usage: $0 [output-directory]" >&2
  exit 2
fi

cleanup_output=false
if [[ "$#" -eq 1 ]]; then
  output_directory="$1"
else
  test_root="$(mktemp -d /tmp/gigas-wgbs-pe-test.XXXXXX)"
  output_directory="${test_root}/result"
  cleanup_output=true
fi

"${repo_root}/bisulfite_analysis/WGBS/code/run_pe_report_summary.sh" "${output_directory}"
Rscript "${script_dir}/validate_pe_summary.R" "${output_directory}" "${report_directory}"

cd "${repo_root}"
shasum -a 256 -c bisulfite_analysis/WGBS/provenance/pe-reports.sha256 >/dev/null
git -c core.fsmonitor=false diff --check

if [[ "${cleanup_output}" = true ]]; then
  rm -rf "${test_root}"
fi

echo "PE report checks passed: parsing, source percentages, group statistics, and input checksums."
