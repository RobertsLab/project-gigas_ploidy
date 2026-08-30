# WGBS Analysis Code

Run scripts in order:

## 1. PE Report Analysis
**`1_WGBS_PE_report_analysis.R`**
Parses Bismark paired-end alignment reports from `data/PE_reports/` and summarizes mapping rates and percent methylation across all samples.

## 2. MethylKit DML Calling
**`2_WGBS_Methylkit.R`**
Loads COV files (from gannet server), filters by coverage (≥10x), and uses MethylKit to identify differentially methylated loci (DMLs) between diploid and triploid samples. Outputs:
- `DML/DML-getMethylDiff-ploidy-Cov10-*.csv`
- `DML/DML-getMethylDiff-ploidy-Cov10-*.bed`
- PCA and clustering plots in `plots/`

## 3. DML Genomic Feature Location
**`3_WGBS_DML_feature_location.ipynb`**
Jupyter notebook. Intersects DML BED files with GFF genome feature annotations (from `GFF/`) using bedtools to determine the genomic distribution of DMLs (exons, introns, intergenic, TE, etc.). Output summary in `analyses/DML_locations_genome_feature.xlsx`.

## Reference Genome
All analyses use **cgigas_uk_roslin_v1** (GFF files in `../GFF/`).

## Notes
- `Yaamini-virginica-2018-10-25-MethylKit.Rmd` is a reference template from a related project (C. virginica); not used directly for this analysis.
