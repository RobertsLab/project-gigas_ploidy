# Bismark paired-end reports

The ten `*_bismark_bt2_PE_report.txt` files are the source inputs for `../../code/1_WGBS_PE_report_analysis.R`. The script parses these reports directly and joins `seq_id` to treatment information from `metadata/wgbs_samples.csv`.

`percent_methylation_summary.xlsx` is retained only as legacy analysis history. It contains manually maintained sample labels and summary values and is not an input to the remediated workflow.

The previously absent `zr3534_7` report was recovered on 2026-08-30 from:

`https://gannet.fish.washington.edu/panopea/030521-ronrosM/zr3534_7_R1.fastp-trim.20201202_bismark_bt2_PE_report.txt`

Its SHA-256 checksum is:

`e858d1b34f890323f1d8370214512c10ff5a595d326b497838a3cf264cbec866`

Checksums for the complete ten-report input set are recorded in `../../provenance/pe-reports.sha256`.
