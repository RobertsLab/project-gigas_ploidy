# Analysis Workflow

End-to-end pipeline for the gigas ploidy desiccation WGBS + qPCR study.

---

## Samples

10 ctenidia tissue samples sequenced by WGBS (Zymo-Seq kit):

Canonical machine-readable metadata: `metadata/wgbs_samples.csv`.

| SeqID | Library | Ploidy | Desiccation | Heat Shock | SRA |
|-------|---------|--------|-------------|------------|-----|
| zr3534_1 | D11-C | Diploid | Yes | No | SRX9508698 |
| zr3534_2 | D12-C | Diploid | Yes | No | SRX9508699 |
| zr3534_3 | D13-C | Diploid | Yes | No | SRX9508700 |
| zr3534_4 | D19-C | Diploid | Yes | Yes | SRX9508701 |
| zr3534_5 | D20-C | Diploid | Yes | Yes | SRX9508702 |
| zr3534_6 | T11-C | Triploid | Yes | No | SRX9508703 |
| zr3534_7 | T12-C | Triploid | Yes | No | SRX9508704 |
| zr3534_8 | T13-C | Triploid | Yes | No | SRX9508705 |
| zr3534_9 | T19-C | Triploid | Yes | Yes | SRX9508706 |
| zr3534_10 | T20-C | Triploid | Yes | Yes | SRX9508707 |

SRA BioProject: [PRJNA678408](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA678408)

All ten Bismark PE reports are committed. The previously absent `zr3534_7` report was recovered from the project archive; see `bisulfite_analysis/WGBS/data/PE_reports/README.md` for provenance.

---

## WGBS Pipeline

### Step 1 — Adapter trimming
**Script:** `bisulfite_analysis/WGBS/sbatch_scripts/20210316_cgig_fastp_ronit-ploidy-wgbs.sh`  
**Tool:** fastp  
**Input:** Raw FASTQ (paired-end)  
**Output:** Trimmed FASTQ

### Step 2 — Alignment
**Tool:** Bismark + bowtie2  
**Reference genome:** *C. gigas* Roslin v1 (`cgigas_uk_roslin_v1`)  
**Output:** BAM files at https://gannet.fish.washington.edu/panopea/030521-ronrosM/  
Alignment statistics summarized in PE report `.txt` files (`bisulfite_analysis/WGBS/data/PE_reports/`).

### Step 3 — PE Report Analysis
**Script:** `bisulfite_analysis/WGBS/code/1_WGBS_PE_report_analysis.R`  
**Input:** Bismark PE report `.txt` files  
**Command:** `bisulfite_analysis/WGBS/code/run_pe_report_summary.sh`
**Output:** a run-specific directory under `bisulfite_analysis/WGBS/results/pe_report_summary/` containing source-derived per-sample and group CSV files. The old `data/PE_reports/percent_methylation_summary.xlsx` is retained as legacy history only.

### Step 4 — DML Calling
**Script:** `bisulfite_analysis/WGBS/code/2_WGBS_Methylkit.R`  
**Tool:** MethylKit (R)  
**Input:** COV files at https://gannet.fish.washington.edu/panopea/WGBS-gigas-ploidy-desiccation/cov_files/  
**Primary model:** ploidy effect adjusted for heat shock; MN overdispersion; F test; SLIM q-values; predicted effect; q ≤0.01
**Command:** `bisulfite_analysis/WGBS/code/run_methylkit_dml.sh`
**Output:** complete and thresholded tables in a new directory under `bisulfite_analysis/WGBS/results/dml/`. See `bisulfite_analysis/WGBS/DML_ANALYSIS.md`.

The historical files in `bisulfite_analysis/WGBS/DML/` and `plots/` remain legacy pending a successful full rerun in the pinned environment.

### Step 5 — DML Genomic Feature Location
**Script:** `bisulfite_analysis/WGBS/code/3_WGBS_DML_feature_location.ipynb`  
**Tool:** bedtools (via Jupyter notebook)  
**Input:** DML BED files + GFF annotations in `bisulfite_analysis/WGBS/GFF/`  
**Genome features available:** gene, mRNA, exonUTR, intron, intergenic, lncRNA, upstream, downstream, flanks, TE  
**Output:** `bisulfite_analysis/WGBS/analyses/DML_locations_genome_feature.xlsx`

---

## qPCR Pipeline

### Step 6 — dCt Analysis
**Script:** `gene_expression/scripts/dct_analysis_NA-to-45.R`  
**Input:** `gene_expression/data/qpcr_ct_values/qpcr_data_consolidated.csv`  
**Reference gene:** Actin  
**Genes analyzed:** HSC70, DNMT1, MBD2, MeCP2, HIF1A, HATHaP2, HAT, HSP90, SOD, ATPsynthetase, COX1  
**NA handling:** Undetermined Ct values set to 45  
**Statistics:** Two-way ANOVA (`Ploidy + Desiccation + Ploidy:Desiccation`) + Tukey HSD  
**Output:** Boxplots in `gene_expression/analyses/`

---

## Global DNA Methylation (ELISA)

### Step 7 — MethylFlash ELISA Analysis
**Script:** `bisulfite_analysis/ELISA/GlobalDNAMeth_Polyploids.Rmd`  
**Input:** Raw plate reader data in `bisulfite_analysis/ELISA/docs/`  
**Kit:** MethylFlash Global DNA Methylation ELISA (EpiGentek P-1030)  
**Output:** `bisulfite_analysis/ELISA/figures/5mC_figure.png` — 5mC % by ploidy/treatment

---

## Key Results Files

| File | Description |
|------|-------------|
| `bisulfite_analysis/WGBS/DML/DML-getMethylDiff-ploidy-Cov10-20.csv` | All DMLs (≥10x cov, ≥20% diff) |
| `bisulfite_analysis/WGBS/DML/DML-getMethylDiff-ploidy-Cov10-50.bed` | High-confidence DMLs (≥50% diff) |
| `bisulfite_analysis/WGBS/analyses/DML_locations_genome_feature.xlsx` | DML genomic distribution |
| `bisulfite_analysis/WGBS/analyses/coverage_summary.xlsx` | Per-sample CpG coverage stats |
| `gene_expression/analyses/` | Per-gene qPCR expression boxplots |
| `bisulfite_analysis/ELISA/figures/5mC_figure.png` | Global 5mC % |
