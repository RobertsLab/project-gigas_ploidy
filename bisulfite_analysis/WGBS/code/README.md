# WGBS Analysis Code

> **Remediation status:** PE-report summarization has been replaced with a clean, source-derived workflow. DML calling and feature annotation still contain known correctness and reproducibility issues and must not be treated as canonical until the remaining phases in `../../../METHYLATION_REMEDIATION_PLAN.md` are completed. Draft clean-run settings are recorded in `../config/analysis_config.json`, and provenance checksums are under `../provenance/`.

Run scripts in order:

## 1. PE Report Analysis
**`1_WGBS_PE_report_analysis.R`**
Parses all ten Bismark paired-end reports, validates reported percentages against their underlying counts, joins canonical sample metadata, and calculates group statistics without reading the legacy workbook. It uses base R and accepts explicit command-line paths.

Run it through the configuration-aware wrapper:

```bash
bisulfite_analysis/WGBS/code/run_pe_report_summary.sh
```

This writes `per_sample.csv` and `group_summary.csv` to a new timestamped directory under `results/pe_report_summary/`. Run the corresponding verification with `bisulfite_analysis/WGBS/tests/run_pe_report_checks.sh`.

## 2. MethylKit DML Calling
**`2_WGBS_Methylkit.R`**
Runs the primary ploidy effect adjusted for heat shock from explicit external COV inputs. The script uses methylKit's tabix-backed objects, preserves all biological replicates, writes the complete statistics table before threshold subsets, and generates paired CSV/BED files from the same rows.

Use `download_cov_inputs.sh` to retrieve the ten external inputs, then run through the configuration-aware wrapper:

```bash
export GIGAS_WGBS_COV_DIR=/path/to/wgbs-cov
bisulfite_analysis/WGBS/code/run_methylkit_dml.sh --preflight-only
bisulfite_analysis/WGBS/code/run_methylkit_dml.sh
```

The full run requires the pinned container in `../environment/`. See `../DML_ANALYSIS.md` for the input contract, model, filtering rationale, output definitions, and container command. PCA and clustering are deliberately excluded from this script until the unbiased all-CpG plotting workflow is implemented.

## 3. DML Genomic Feature Location
**`3_WGBS_DML_feature_location.ipynb`**
Jupyter notebook. Intersects DML BED files with GFF genome feature annotations (from `GFF/`) using bedtools to determine the genomic distribution of DMLs (exons, introns, intergenic, TE, etc.). Output summary in `analyses/DML_locations_genome_feature.xlsx`.

## Reference Genome
All analyses use **cgigas_uk_roslin_v1** (GFF files in `../GFF/`).

## Notes
- `Yaamini-virginica-2018-10-25-MethylKit.Rmd` is a reference template from a related project (C. virginica); not used directly for this analysis.
