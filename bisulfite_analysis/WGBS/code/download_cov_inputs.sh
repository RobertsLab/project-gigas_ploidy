#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../../.." && pwd)"
config_path="${repo_root}/bisulfite_analysis/WGBS/config/analysis_config.json"
metadata_path="${repo_root}/$(jq -er '.sample_metadata' "${config_path}")"

if [[ "$#" -ne 1 ]]; then
  echo "Usage: $0 target-directory" >&2
  exit 2
fi
target_directory="$1"
mkdir -p "${target_directory}"
target_directory="$(cd "${target_directory}" && pwd)"

source_url="$(jq -er '.inputs.cov_files.source_url' "${config_path}")"
template="$(jq -er '.inputs.cov_files.file_name_template' "${config_path}")"
sequence_ids=()
while IFS= read -r seq_id; do
  sequence_ids+=("${seq_id}")
done < <(awk -F, 'NR > 1 {print $1}' "${metadata_path}")

for seq_id in "${sequence_ids[@]}"; do
  cov_file="${template//\{seq_id\}/${seq_id}}"
  curl --fail --location --retry 3 --continue-at - \
    --output "${target_directory}/${cov_file}" \
    "${source_url}${cov_file}"
done

(
  cd "${target_directory}"
  shasum -a 256 ./*.cov > wgbs_cov_inputs.sha256
)

echo "Downloaded ${#sequence_ids[@]} COV files and wrote ${target_directory}/wgbs_cov_inputs.sha256"
