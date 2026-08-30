#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../../.." && pwd)"
fixture_directory="${script_dir}/fixtures/cov"
test_root="$(mktemp -d /tmp/gigas-wgbs-dml-preflight.XXXXXX)"
output_directory="${test_root}/result"

export GIGAS_WGBS_COV_DIR="${fixture_directory}"
"${repo_root}/bisulfite_analysis/WGBS/code/run_methylkit_dml.sh" \
  --preflight-only "${output_directory}"
Rscript --vanilla "${script_dir}/validate_dml_preflight.R" "${output_directory}"
Rscript --vanilla "${script_dir}/validate_dml_outputs.R" "${script_dir}/fixtures/dml_output"

if "${repo_root}/bisulfite_analysis/WGBS/code/run_methylkit_dml.sh" \
  --preflight-only "${output_directory}" >/dev/null 2>&1; then
  echo "DML workflow unexpectedly accepted a non-empty output directory." >&2
  exit 1
fi

missing_directory="${test_root}/missing"
mkdir "${missing_directory}"
if GIGAS_WGBS_COV_DIR="${missing_directory}" \
  "${repo_root}/bisulfite_analysis/WGBS/code/run_methylkit_dml.sh" \
  --preflight-only "${test_root}/missing-result" >/dev/null 2>&1; then
  echo "DML workflow unexpectedly accepted missing COV inputs." >&2
  exit 1
fi

rm -rf "${test_root}"
git -C "${repo_root}" -c core.fsmonitor=false diff --check
echo "DML preflight checks passed: filenames, COV schema, coordinates, counts, metadata, and model parameters."
