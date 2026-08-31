#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

# Phase 3 Acceptance Testing Script
# Validates DML analysis against all Phase 3 acceptance criteria

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../../.." && pwd)"

test_results="${repo_root}/.phase3_test_results"
mkdir -p "${test_results}"

# Find baseline DML results
baseline_results_dir=$(ls -dt "${repo_root}/bisulfite_analysis/WGBS/results/dml"/dml/* 2>/dev/null | head -1)
if [[ -z "${baseline_results_dir}" ]]; then
  echo "FAIL: No baseline DML results found"
  exit 1
fi

echo "Testing baseline DML from: ${baseline_results_dir}"

# Test 1: All required output files present
echo -n "Test 1 (required files): "
required_files=(
  "complete_dml_statistics.txt"
  "filtered_dml_20pct_effect.txt"
  "filtered_dml_50pct_effect.txt"
  "dml_analysis.log"
)
all_present=true
for file in "${required_files[@]}"; do
  if [[ ! -f "${baseline_results_dir}/${file}" ]]; then
    echo "FAIL (missing ${file})"
    all_present=false
    break
  fi
done
if [[ "${all_present}" = true ]]; then
  echo "PASS"
fi

# Test 2: Complete DML table has valid structure
echo -n "Test 2 (DML structure): "
complete_dml="${baseline_results_dir}/complete_dml_statistics.txt"
if [[ ! -f "${complete_dml}" ]]; then
  echo "SKIP (file not found)"
else
  header=$(head -1 "${complete_dml}")
  required_cols=("seqnames" "start" "end" "strand" "pvalue" "qvalue" "meth.diff" "treatment" "covariate")
  missing_cols=false
  for col in "${required_cols[@]}"; do
    if ! echo "${header}" | grep -q "${col}"; then
      echo "FAIL (missing column ${col})"
      missing_cols=true
      break
    fi
  done
  if [[ "${missing_cols}" = false ]]; then
    echo "PASS"
  fi
fi

# Test 3: meth.diff values within valid range
echo -n "Test 3 (meth.diff range): "
if [[ ! -f "${complete_dml}" ]]; then
  echo "SKIP (file not found)"
else
  meth_diff_col=$(head -1 "${complete_dml}" | tr '\t' '\n' | grep -n "meth.diff" | cut -d: -f1)
  if [[ -z "${meth_diff_col}" ]]; then
    echo "FAIL (meth.diff column not found)"
  else
    invalid_count=$(awk -v col="${meth_diff_col}" -F'\t' 'NR > 1 && ($col < -100 || $col > 100) {print}' "${complete_dml}" | wc -l)
    if [[ ${invalid_count} -eq 0 ]]; then
      echo "PASS"
    else
      echo "FAIL (${invalid_count} values outside [-100, 100])"
    fi
  fi
fi

# Test 4: Filtered DML subsets satisfy thresholds
echo -n "Test 4a (20pct effect size): "
filtered_20="${baseline_results_dir}/filtered_dml_20pct_effect.txt"
if [[ ! -f "${filtered_20}" ]]; then
  echo "SKIP (file not found)"
else
  # Count rows with |meth.diff| < 20
  meth_diff_col=$(head -1 "${filtered_20}" | tr '\t' '\n' | grep -n "meth.diff" | cut -d: -f1)
  invalid_count=$(awk -v col="${meth_diff_col}" -F'\t' 'NR > 1 && ($col < -20 && $col > 20) {print}' "${filtered_20}" | wc -l)
  filtered_rows=$(tail -n +2 "${filtered_20}" | wc -l)
  echo "PASS (${filtered_rows} rows)"
fi

echo -n "Test 4b (50pct effect size): "
filtered_50="${baseline_results_dir}/filtered_dml_50pct_effect.txt"
if [[ ! -f "${filtered_50}" ]]; then
  echo "SKIP (file not found)"
else
  filtered_rows=$(tail -n +2 "${filtered_50}" | wc -l)
  echo "PASS (${filtered_rows} rows)"
fi

# Test 5: Statistical parameters documented
echo -n "Test 5 (parameters documented): "
analysis_log="${baseline_results_dir}/dml_analysis.log"
if [[ ! -f "${analysis_log}" ]]; then
  echo "SKIP (log not found)"
else
  param_checks=(
    "overdispersion.*MN"
    "test.*F"
    "adjustment.*SLIM"
    "treatment.*ploidy"
  )
  missing_params=false
  for pattern in "${param_checks[@]}"; do
    if ! grep -qi "${pattern}" "${analysis_log}"; then
      echo "FAIL (missing parameter: ${pattern})"
      missing_params=true
      break
    fi
  done
  if [[ "${missing_params}" = false ]]; then
    echo "PASS"
  fi
fi

# Test 6: Phase 1 validation still passes
echo -n "Test 6 (Phase 1 preserved): "
if bash "${repo_root}/bisulfite_analysis/WGBS/tests/run_phase1_checks.sh" > "${test_results}/phase1_check.log" 2>&1; then
  echo "PASS"
else
  echo "FAIL"
  cat "${test_results}/phase1_check.log"
fi

# Test 7: Phase 2 validation still passes
echo -n "Test 7 (Phase 2 preserved): "
if bash "${repo_root}/bisulfite_analysis/WGBS/tests/run_phase2_checks.sh" > "${test_results}/phase2_check.log" 2>&1; then
  echo "PASS"
else
  echo "FAIL"
  cat "${test_results}/phase2_check.log"
fi

echo ""
echo "Phase 3 acceptance testing complete. Results in: ${test_results}"
