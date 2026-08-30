# DNA Methylation Analysis Remediation Plan

This plan addresses the correctness and reproducibility issues identified in the WGBS and global 5mC ELISA analyses. The existing files are retained as analysis history, but affected results should not be treated as final until the acceptance criteria below pass.

## Current result status

| Result | Status | Reason |
|---|---|---|
| WGBS DML CSV and BED files | Invalid pending rerun | The Cov10 50% CSV and BED contain different locus sets, and the R script can reuse objects loaded from an external R session. |
| WGBS DML PCA and heatmap | Invalid pending rerun | DMLs are selected from stale objects and matched to methylation rows using start position alone. The PCA also reuses the data that selected the DMLs. |
| DML genomic-feature summary | Provisional | The checked-in counts reproduce from the DML CSVs and GFF/BED files, but the notebook does not successfully generate them and uses a conflicting BED file. |
| WGBS global methylation summary plot | Invalid pending recalculation | The workbook's summary values do not equal summaries calculated from its raw sheet. |
| ELISA 5mC calculations and ANOVA | Provisional | The standard-curve equation and reported tests reproduce, assuming every sample well contained 25 ng DNA. |
| ELISA treatment labels and final figure | Requires correction | Desiccation samples are labeled as heat stress in code, and the 0-3% y-axis limit removes nine observations before boxplot statistics are calculated. |

## Phase 1: Freeze provenance and establish one configuration

- [ ] Record the current commit hash and checksums for all files under `bisulfite_analysis/WGBS/DML/`, `bisulfite_analysis/WGBS/analyses/`, and `bisulfite_analysis/WGBS/plots/` in a legacy-results manifest.
- [ ] Mark current affected outputs as legacy in their README files; do not silently overwrite them during development.
- [ ] Create one machine-readable WGBS sample metadata file and use it everywhere. It must contain `seq_id`, `library_name`, `ploidy`, `desiccation`, `heat_shock`, and source accession.
- [ ] Reconcile the incorrect library names in `percent_methylation_summary.xlsx` with `gene_expression/data/zr3534_wgbs_info.csv` (D19/D20 and T19/T20 for the heat-shocked WGBS libraries).
- [ ] Add a project configuration file for input and output directories, genome assembly, coverage thresholds, methylation-difference thresholds, and statistical settings.
- [ ] Capture R, Bioconductor, methylKit, Bismark, bedtools, and other package/tool versions in a reproducible environment definition.

Acceptance criteria:

- A single metadata row controls every sample-to-treatment assignment.
- No analysis script contains user-specific `E:/`, `/Volumes/`, Dropbox, or RStudio document paths.
- No script loads a saved `.RData` workspace or depends on objects created outside that script.

## Phase 2: Make WGBS preprocessing and report summaries reproducible

- [ ] Separate cluster preprocessing from downstream analysis and document the exact Bismark reference files and command versions.
- [ ] Retain deduplication and merged-CpG generation, with zero-based half-open coordinates documented at the file boundary.
- [ ] Replace the workbook-only behavior in `1_WGBS_PE_report_analysis.R` with parsing of the committed Bismark PE reports plus the canonical metadata file.
- [ ] Add or retrieve the missing `zr3534_7` PE report, or explicitly record its external checksum and source URL.
- [ ] Calculate group means, standard deviations, sample counts, and standard errors directly from the parsed per-sample table. Do not maintain a manually edited summary sheet.
- [ ] Test that all expected sequence IDs occur exactly once and that parsed CpG methylation percentages match the Bismark reports.

Acceptance criteria:

- The summary table is generated from raw reports in one command.
- Group summaries exactly equal recalculation from the generated per-sample table.
- Sample names and heat-shock assignments match the canonical metadata.

## Phase 3: Rebuild the DML analysis from a clean session

- [ ] Rewrite `2_WGBS_Methylkit.R` as a non-interactive command-line script using project-relative paths and explicit inputs.
- [ ] Remove `load(...)`, `rstudioapi`, mutable `setwd(...)` chains, and all references to undefined objects such as `fileName`, `differentialMethylationStatsTreatment`, and `diffMethStatsTreatment50`.
- [ ] Fail immediately when a package, input file, sample, output directory, or metadata field is missing.
- [ ] Preserve sample-level biological replicates; do not pool animals.
- [ ] Use an overdispersion-aware model as the primary analysis, such as methylKit with `overdispersion="MN"` and its corresponding test, or a documented beta-binomial DSS analysis.
- [ ] Estimate the ploidy effect while accounting for heat shock. Define the primary estimand before running the model:
  - primary: ploidy effect adjusted for heat-shock status;
  - secondary: ploidy-by-heat interaction or stratified ploidy effects, using a model that supports that design.
- [ ] Document the rationale for coverage normalization, high-coverage filtering, the 10x threshold, and the requirement that a CpG be observed in all samples.
- [ ] Run sensitivity analyses for the coverage/occupancy rule and multiple-testing method.
- [ ] Write a complete, unfiltered DML statistics table before producing 20% and 50% subsets.

Acceptance criteria:

- A clean-session rerun produces all DML tables without warnings or pre-existing objects.
- Every `meth.diff` is within -100 to 100.
- Every filtered DML satisfies both its declared adjusted-p-value and effect-size thresholds.
- Treatment coefficients and the sign convention for `meth.diff` are documented.
- The primary result is stable enough under the predefined sensitivity analyses to support the stated conclusion, or the conclusion is revised.

## Phase 4: Generate consistent DML files and exploratory plots

- [ ] Generate CSV and BED outputs from the same in-memory DML object in the same run.
- [ ] Use unambiguous names containing the comparison, coverage threshold, effect threshold, model, and run identifier.
- [ ] Add automated checks that paired CSV/BED files have identical row counts and identical chromosome/start/end keys.
- [ ] Match methylation data by chromosome, start, and end; never by start position alone.
- [ ] Base the primary unsupervised PCA and clustering figures on the predefined all-CpG matrix, not on loci selected for group differences.
- [ ] If DML-only PCA or heatmaps are retained, label them explicitly as descriptive and circular. Use held-out samples or cross-validation for any claim of predictive separation.
- [ ] Derive axis variance labels from the fitted PCA object rather than hard-coding percentages.

Acceptance criteria:

- Each CSV/BED pair passes exact coordinate-set equality.
- No genomic join has one-to-many matches unless explicitly expected and tested.
- PCA labels equal the variance explained in the current run.
- Plot legends and annotations come from canonical metadata rather than positional assumptions.

## Phase 5: Replace the feature-location notebook with an executable workflow

- [ ] Convert `3_WGBS_DML_feature_location.ipynb` into a script or fully executable notebook with project-relative paths and a documented bedtools dependency.
- [ ] Remove saved error output and Windows/macOS user-specific paths.
- [ ] Use only the validated BED generated in Phase 4.
- [ ] Encode coordinate conventions explicitly and add boundary tests for BED versus GFF coordinates.
- [ ] Generate every reported feature class in code, including the currently undocumented exon and CDS rows.
- [ ] Report both overlap counts and denominators; document that feature categories may overlap and therefore should not necessarily sum to the DML total.
- [ ] Generate `DML_locations_genome_feature.xlsx` directly from the computed results and include run metadata.

Acceptance criteria:

- The workflow runs from start to finish without manual shell cells.
- Recomputed counts match direct interval-overlap checks.
- The workbook total equals the validated DML CSV/BED total for each threshold.

## Phase 6: Correct and harden the ELISA analysis

- [ ] Rename the D11-D18 and T11-T18 contrast from `heat_stress` to `desiccation` throughout code, prose, tables, and plot labels.
- [ ] Keep acute heat-shock terminology only for D19-D20 and T19-T20 samples where applicable.
- [ ] Confirm from the laboratory record that every analyzed sample well received exactly 25 ng DNA. If inputs varied, calculate 5mC using the actual per-well DNA amount.
- [ ] Preserve raw plate-reader values, plate maps, and technical-replicate flags without editing them in place.
- [ ] Flag D02-C and T08-C as having only one usable technical well and decide prospectively whether they remain in the primary analysis.
- [ ] Compute standard concentrations from their sample labels or an explicit lookup table rather than relying on row order.
- [ ] Check model residuals and variance assumptions, and document the chosen model and any sensitivity analysis.
- [ ] Run Tukey comparisons from the same final model used for the reported ANOVA.
- [ ] Remove `scale_y_continuous(limits=c(0,3))`. Display the full data range, or use `coord_cartesian()` only for a clearly labeled visual zoom.
- [ ] Overlay individual biological observations on boxplots and report group sample sizes.

Acceptance criteria:

- The standard curve, negative-control subtraction, DNA-mass adjustment, and all sample 5mC values reproduce from raw wells.
- No plotted observation is removed by an axis scale.
- Statistical tables, prose conclusions, and figures use the same treatment labels and model.
- The reported ANOVA and post-hoc values reproduce in a clean session.

## Phase 7: Validation, documentation, and release

- [ ] Add automated tests for metadata completeness, coordinate uniqueness, DML thresholds, CSV/BED equality, summary statistics, and plot input ranges.
- [ ] Add a single documented command or workflow target for each of: WGBS report summary, DML calling, feature annotation, and ELISA analysis.
- [ ] Generate outputs in a new run-specific directory before promoting them to the canonical result paths.
- [ ] Perform an independent review of sample metadata, model design, and regenerated figures.
- [ ] Update `WORKFLOW.md`, the WGBS README, and the project README only after all acceptance checks pass.
- [ ] Publish a result manifest containing input checksums, parameters, software versions, output checksums, and the commit used for the final run.

Final release criteria:

- All analyses run from a fresh checkout without RStudio or a saved workspace.
- All canonical results are traceable to raw or externally checksummed inputs.
- No unresolved disagreement exists among metadata, CSV, BED, workbook, and figure outputs.
- Scientific conclusions are rewritten to match the validated model and treatment definitions.

