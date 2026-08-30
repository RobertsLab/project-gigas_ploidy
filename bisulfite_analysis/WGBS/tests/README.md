# WGBS remediation checks

Run the current Phase 1 checks from any directory with:

```bash
bash bisulfite_analysis/WGBS/tests/run_phase1_checks.sh
```

The Phase 1 checks validate the canonical sample design, require a unique Bismark PE report for every sample, parse the JSON configuration, and verify all 59 frozen legacy result checksums.

Run the Phase 2 PE-report parser and its independent output checks with:

```bash
bisulfite_analysis/WGBS/tests/run_pe_report_checks.sh
```

The command creates a timestamped result directory unless an explicit output directory is supplied as its first argument. It verifies the ten directly parsed CpG percentages and recomputes every group statistic from the generated per-sample CSV.

Run the package-independent DML input and output-contract tests with:

```bash
bisulfite_analysis/WGBS/tests/run_dml_preflight_checks.sh
```

The test validates all ten expected COV filenames, six-column previews, zero-based two-base CpG intervals, count-derived percentages, canonical metadata joins, the primary model declaration, complete-statistics invariants, effect/q-value thresholds, and exact CSV/BED coordinate equality. The small fixtures test workflow contracts only; they are not biological results and do not invoke methylKit.
