# qPCR Scripts

## Canonical Script

**`dct_analysis_NA-to-45.R`** — current version. Run this for all analyses.

### What it does
1. Reads `data/qpcr_ct_values/qpcr_data_consolidated.csv`
2. Parses `Sample` column into Ploidy, Desiccation, HeatShock, SampleNum
3. Replaces `NA` Ct values with 45 (convention: undetected = no amplification)
4. Calculates normalized expression: `2^-(gene_Ct - Actin_Ct)` for all targets
5. Log-transforms for ANOVA normality assumption
6. Runs two-way ANOVA (`Ploidy + Desiccation + Ploidy:Desiccation`) + Tukey HSD for each gene
7. Produces boxplots for all genes

### Genes analyzed
HSC70, DNMT1, MBD2, MeCP2, HIF1A, HATHaP2, HAT, HSP90, SOD, ATPsynthetase, COX1

## Deprecated

**`dct_analysis_DEPRECATED.R`** — do not use. Earlier version with an unresolved merge conflict and missing genes (SOD, ATPsynthetase, COX1).
