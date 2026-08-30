#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../../.." && pwd)"
config_path="${repo_root}/bisulfite_analysis/WGBS/config/analysis_config.json"
manifest_path="${repo_root}/bisulfite_analysis/WGBS/provenance/legacy-results-a910a18.sha256"

cd "${repo_root}"

jq -e . "${config_path}" >/dev/null
export SAMPLE_METADATA_PATH
SAMPLE_METADATA_PATH="$(jq -r '.sample_metadata' "${config_path}")"
export ALLOWED_MISSING_PE_REPORTS
ALLOWED_MISSING_PE_REPORTS="$(jq -r '.inputs.allowed_missing_pe_reports | join(",")' "${config_path}")"

Rscript "${script_dir}/validate_phase1.R"

manifest_lines="$(wc -l < "${manifest_path}" | tr -d ' ')"
if [[ "${manifest_lines}" != "59" ]]; then
  echo "Expected 59 legacy result checksums, found ${manifest_lines}." >&2
  exit 1
fi
shasum -a 256 -c "${manifest_path}" >/dev/null

git -c core.fsmonitor=false diff --check

echo "Phase 1 checks passed: configuration, metadata, PE report exceptions, and 59 legacy checksums."
