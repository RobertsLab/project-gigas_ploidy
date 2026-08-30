# project-gigas_ploidy

Project investigating the effects of desiccation and elevated temperature exposure on **triploid and diploid Pacific oysters** (*Crassostrea gigas*).

> [!IMPORTANT]
> The current methylation analyses contain known correctness and reproducibility issues. See the [DNA Methylation Analysis Remediation Plan](METHYLATION_REMEDIATION_PLAN.md) before using the WGBS or ELISA results.

## Manuscripts & Presentations

- [Manuscript (Google Doc)](https://docs.google.com/document/d/17mcGDI-TWmU4vgBXmiXmeofe4qEuFH5inBKBHhG9tzg/edit)
- [Ronit's manuscript draft (Google Doc)](https://docs.google.com/document/d/1nwY9I3pVzF5Xlfdzb7SdO1QKaNp5qyj0DuPCCa_YXi8/edit?usp=sharing)
- [Poster (presentations/)](presentations/Triploid-Diploid-Oyster-Poster_RonitJain.pptx)

---

## Repository Structure

```
project-gigas_ploidy/
├── bisulfite_analysis/
│   ├── ELISA/          # Global DNA methylation assay (R scripts, raw data, figures)
│   └── WGBS/           # Whole-genome bisulfite sequencing analysis
│       ├── code/       # R scripts and Jupyter notebook (run in order 1→3)
│       ├── data/       # Bismark PE alignment reports
│       ├── DML/        # Differentially methylated loci (CSV, BED)
│       ├── GFF/        # Genome feature annotations (Roslin v1)
│       ├── analyses/   # Summary tables
│       ├── plots/      # PCA, heatmap, methylation/coverage figures
│       ├── results/    # New run-specific remediated outputs
│       └── sbatch_scripts/  # Cluster job scripts
├── gene_expression/
│   ├── scripts/        # qPCR dCt analysis R scripts
│   ├── data/           # qPCR Ct values (per gene CSVs) and WGBS sample info
│   └── analyses/       # Per-gene expression plots
├── manuscript/         # Local .docx backup of manuscript
└── presentations/      # Conference posters and slides
```

---

## Remote Data

Large files are hosted externally:

| Type | Location |
|------|----------|
| BED files | https://gannet.fish.washington.edu/panopea/WGBS-gigas-ploidy-desiccation/bed_files/ |
| COV files | https://gannet.fish.washington.edu/panopea/WGBS-gigas-ploidy-desiccation/cov_files/ |
| BAM files (Bismark output) | https://gannet.fish.washington.edu/panopea/030521-ronrosM/ |

---

## Analysis Workflow

1. **Trim reads** — `bisulfite_analysis/WGBS/sbatch_scripts/20210316_cgig_fastp_ronit-ploidy-wgbs.sh`
2. **Align** — Bismark (run on cluster; BAMs at gannet link above)
3. **PE report analysis** — `bisulfite_analysis/WGBS/code/run_pe_report_summary.sh`
4. **DML calling** — `bisulfite_analysis/WGBS/code/run_methylkit_dml.sh` (pinned MethylKit environment; see `bisulfite_analysis/WGBS/DML_ANALYSIS.md`)
5. **DML genomic feature location** — `bisulfite_analysis/WGBS/code/3_WGBS_DML_feature_location.ipynb`
6. **qPCR dCt analysis** — `gene_expression/scripts/dct_analysis_NA-to-45.R`
7. **Global DNA methylation (ELISA)** — `bisulfite_analysis/ELISA/GlobalDNAMeth_Polyploids.Rmd`

---

## Sample Guide

| ID Range | Ploidy | Treatment |
|----------|--------|-----------|
| D01–D08 | Diploid | Control (aquarium water) |
| D09–D10 | Diploid | Control + 1 hr heat shock (45°C) |
| D11–D18 | Diploid | Desiccation + 27°C for 24 hrs |
| D19–D20 | Diploid | Desiccation + 27°C for 24 hrs + 1 hr heat shock (45°C) |
| T01–T08 | Triploid | Control (aquarium water) |
| T09–T10 | Triploid | Control + 1 hr heat shock (45°C) |
| T11–T18 | Triploid | Desiccation + 27°C for 24 hrs |
| T19–T20 | Triploid | Desiccation + 27°C for 24 hrs + 1 hr heat shock (45°C) |

WGBS sequencing IDs: `zr3534_1` through `zr3534_10` (see the canonical metadata in `metadata/wgbs_samples.csv`).
