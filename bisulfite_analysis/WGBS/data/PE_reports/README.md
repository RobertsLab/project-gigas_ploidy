# Bismark Paired-End Reports

## Overview

The ten `*_bismark_bt2_PE_report.txt` files are the source inputs for `../../code/1_WGBS_PE_report_analysis.R`. The script parses these reports directly and joins `seq_id` to treatment information from `metadata/wgbs_samples.csv`.

## Files

| Sequence ID | Library Name | File Name | CpG % | Status |
|-------------|--------------|-----------|-------|--------|
| zr3534_1 | D11-C | zr3534_1_R1.fastp-trim.20201202_bismark_bt2_PE_report.txt | 10.5 | ✅ Local |
| zr3534_2 | D12-C | zr3534_2_R1.fastp-trim.20201202_bismark_bt2_PE_report.txt | 10.4 | ✅ Local |
| zr3534_3 | D13-C | zr3534_3_R1.fastp-trim.20201202_bismark_bt2_PE_report.txt | 11.8 | ✅ Local |
| zr3534_4 | D19-C | zr3534_4_R1.fastp-trim.20201202_bismark_bt2_PE_report.txt | 11.2 | ✅ Local |
| zr3534_5 | D20-C | zr3534_5_R1.fastp-trim.20201202_bismark_bt2_PE_report.txt | 11.2 | ✅ Local |
| zr3534_6 | T11-C | zr3534_6_R1.fastp-trim.20201202_bismark_bt2_PE_report.txt | 10.3 | ✅ Local |
| zr3534_7 | T12-C | zr3534_7_R1.fastp-trim.20201202_bismark_bt2_PE_report.txt | 10.1 | ✅ Recovered |
| zr3534_8 | T13-C | zr3534_8_R1.fastp-trim.20201202_bismark_bt2_PE_report.txt | 10.3 | ✅ Local |
| zr3534_9 | T19-C | zr3534_9_R1.fastp-trim.20201202_bismark_bt2_PE_report.txt | 10.0 | ✅ Local |
| zr3534_10 | T20-C | zr3534_10_R1.fastp-trim.20201202_bismark_bt2_PE_report.txt | 10.0 | ✅ Local |

## Recovered File

The previously absent `zr3534_7` report was recovered on 2026-08-30 from:

```
https://gannet.fish.washington.edu/panopea/030521-ronrosM/zr3534_7_R1.fastp-trim.20201202_bismark_bt2_PE_report.txt
```

**SHA-256 checksum:**
```
e858d1b34f890323f1d8370214512c10ff5a595d326b497838a3cf264cbec866
```

## Report Format

Each report contains:
- Bismark version and genome information
- Alignment summary (total pairs, unique hits, mapping efficiency)
- Strand-specific alignment counts (OT, OB, CTOT, CTOB)
- Cytosine methylation summary by context (CpG, CHG, CHH)
- CpG methylation percentage (methylated CpG / all CpG × 100)

Example header:
```
Bismark report for: /gscratch/srlab/sr320/data/cg/zr3534_1_R1.fastp-trim.20201202.fq.gz
                   and /gscratch/srlab/sr320/data/cg/zr3534_1_R2.fastp-trim.20201202.fq.gz
                   (version: v0.21.0)
```

## Validation

All ten reports are validated by `1_WGBS_PE_report_analysis.R`:
- ✅ Exactly one report per sequence ID
- ✅ All expected sequence IDs (zr3534_1 through zr3534_10) are present
- ✅ CpG percentages match direct extraction from reports
- ✅ Reported percentages are mathematically consistent with their underlying counts
- ✅ All reports generated with Bismark v0.21.0

**Validation script**: `../../tests/validate_pe_summary.R`

## Provenance

- **Checksums**: `../../provenance/pe-reports.sha256`
- **Canonical metadata**: `../../../../metadata/wgbs_samples.csv`
- **Legacy metadata** (superseded): `percent_methylation_summary.xlsx` (retained for history only)

## Downstream Use

These reports are parsed by the PE report analysis script:
```bash
bisulfite_analysis/WGBS/code/run_pe_report_summary.sh
```

This generates:
- `per_sample.csv`: Per-sample Bismark statistics and CpG methylation joined with canonical metadata
- `group_summary.csv`: Group statistics (mean, SD, SE) by ploidy and heat-shock treatment

See `../../code/README.md` for the full analysis workflow.
