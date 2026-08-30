# Remediated WGBS results

This directory contains outputs created by remediated workflows. Each execution writes to a new run-specific subdirectory and refuses to overwrite a non-empty run directory.

`pe_report_summary/20260830T133401Z/` is the first verified source-derived PE-report summary. Its two CSV files were generated from all ten committed Bismark reports and `metadata/wgbs_samples.csv` with:

```bash
bisulfite_analysis/WGBS/code/run_pe_report_summary.sh
```

The generated per-sample percentages, cytosine-count consistency, group statistics, source checksums, and Bismark v0.21.0 declarations pass `bisulfite_analysis/WGBS/tests/run_pe_report_checks.sh`.
