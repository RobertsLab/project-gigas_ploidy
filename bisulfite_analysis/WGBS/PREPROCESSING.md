# WGBS Preprocessing and Quality Control

## Overview

This document records the Bismark preprocessing pipeline used to generate the paired-end reports that serve as input to downstream analysis. The preprocessing steps include read trimming, genome-guided alignment, deduplication, and methylation extraction, with outputs ultimately used for DML calling via methylKit.

## Cluster Preprocessing Pipeline

The complete preprocessing was performed on the University of Washington's Hyak computing cluster using the sbatch script: `sbatch_scripts/20210316_cgig_fastp_ronit-ploidy-wgbs.sh`

### Environment
- **Cluster**: Hyak (coenv partition)
- **Working directory**: `/gscratch/scrubbed/sr320/030521-ronrosM`
- **Date executed**: March 2021
- **Resource allocation**: 28 CPU threads per job, 100 GB RAM, 15-day wall time

### Reference Genome
- **Organism**: *Crassostrea gigas* (Pacific oyster)
- **Assembly**: Roslin M (also labeled as roslin_M or cgigas_uk_roslin_v1)
- **Location**: `/gscratch/srlab/sr320/data/Cgig-genome/roslin_M/`
- **Preparation**: Performed with Bismark `bismark_genome_preparation`

### Software Versions

| Tool | Version | Command |
|------|---------|---------|
| Bismark | 0.21.0 | `/gscratch/srlab/programs/Bismark-0.21.0/bismark` |
| Bowtie 2 | 2.3.4.1 | `/gscratch/srlab/programs/bowtie2-2.3.4.1-linux-x86_64/` |
| Samtools | 1.9 | `/gscratch/srlab/programs/samtools-1.9/samtools` |

### Input Files
- **Raw sequences**: Trimmed with fastp (run externally; input files: `{seq_id}_R1/R2.fastp-trim.20201202.fq.gz`)
- **Location**: `/gscratch/srlab/sr320/data/cg/`
- **Naming convention**: `zr3534_{1-10}_R{1,2}.fastp-trim.20201202.fq.gz`

### Processing Steps

#### 1. Genome Preparation
```bash
${bismark_dir}/bismark_genome_preparation \
  --verbose \
  --parallel 28 \
  --path_to_aligner ${bowtie2_dir} \
  ${genome_folder}
```
Prepares the bisulfite-converted reference genome (top and bottom strands) for Bismark alignment.

#### 2. Paired-End Alignment
```bash
${bismark_dir}/bismark \
  --path_to_bowtie ${bowtie2_dir} \
  -genome ${genome_folder} \
  -p 8 \
  -score_min L,0,-0.6 \
  --non_directional \
  -1 ${reads_dir}{}_R1.fastp-trim.20201202.fq.gz \
  -2 ${reads_dir}{}_R2.fastp-trim.20201202.fq.gz
```

**Key parameters:**
- **Non-directional**: Alignments performed to all four DNA strands (OT, OB, CTOT, CTOB)
- **Score threshold**: L,0,-0.6 (length-dependent)
- **Threads**: 8 per sample
- **Output**: `{seq_id}_R1.fastp-trim.20201202_bismark_bt2_pe.bam`

**Outputs for each sample:**
- Paired-end alignment BAM file
- Bismark paired-end report (text) containing Bismark version, alignment statistics, and final cytosine methylation summaries
  - **Location**: Stored in `data/PE_reports/`
  - **Naming**: `{seq_id}_R1.fastp-trim.20201202_bismark_bt2_PE_report.txt`

#### 3. Deduplication
```bash
${bismark_dir}/deduplicate_bismark \
  --bam \
  --paired \
  {}.bam
```
Removes PCR duplicates while retaining all biological replicate information.

**Output**: `{seq_id}_R1.fastp-trim.20201202_bismark_bt2_pe.deduplicated.bam`

#### 4. Methylation Extraction
```bash
${bismark_dir}/bismark_methylation_extractor \
  --bedGraph --counts --scaffolds \
  --multicore 28 \
  --buffer_size 75% \
  *deduplicated.bam
```
Extracts methylation calls and generates intermediate coverage files.

**Output files:**
- `{seq_id}_R1.fastp-trim.20201202_bismark_bt2_pe.deduplicated.bismark.cov.gz` (unmerged coverage)

#### 5. Coordinate Conversion and CpG Merging
```bash
${bismark_dir}/coverage2cytosine \
  --genome_folder ${genome_folder} \
  -o {} \
  --merge_CpG \
  --zero_based \
  {}_bismark_bt2_pe.deduplicated.bismark.cov.gz
```

**Key parameters:**
- **--zero_based**: Outputs coordinates as zero-based, half-open intervals (CHR, START, END)
- **--merge_CpG**: Merges CpG dinucleotides at overlapping positions on opposite strands
- **Genome folder**: Used for reference sequence extraction

**Output files:**
- `{seq_id}_R1.CpG_report.merged_CpG_evidence.cov` (final COV file for downstream analysis)
  - **Format**: Tab-separated, 6 columns
  - **Columns**: Chromosome, Start (0-based), End (0-based, half-open), Methylation%, Methylated count, Unmethylated count
  - **Location**: Externally hosted at `https://gannet.fish.washington.edu/panopea/WGBS-gigas-ploidy-desiccation/cov_files/`

### Coordinate System

All genomic coordinates in the analysis use the following conventions:

| Stage | Coordinate System | Description |
|-------|-------------------|-------------|
| **Merged COV files** (output of this pipeline) | Zero-based, half-open | CHR START END format; compatible with BED files |
| **GFF annotations** | One-based, closed | CHR START END format (standard GFF) |
| **DML BED output** | Zero-based, half-open | CHR START END format; generated from DML calling |

**Coordinate boundary test**: All merged CpG intervals span exactly 2 bases (END - START = 2) and satisfy END > START.

### Quality Checks

Quality metrics reported in the Bismark PE reports (`data/PE_reports/`) include:

- **Mapping efficiency**: Percentage of sequence pairs aligning uniquely (range: 60-62%)
- **Methylation context distribution**: Counts of methylated and unmethylated cytosines in CpG, CHG, and CHH contexts
- **CpG methylation percentage**: Calculated as methylated CpG / (methylated + unmethylated CpG) × 100

### Outputs Used Downstream

The following outputs from this preprocessing pipeline serve as inputs to Phase 3 (DML calling):

1. **Bismark PE reports** → `1_WGBS_PE_report_analysis.R`
   - Used to verify sample identity, mapping efficiency, and CpG methylation summaries
   - Parsed by command-line script that validates percentage calculations and joins canonical metadata

2. **Merged COV files** → `2_WGBS_Methylkit.R` (DML calling)
   - External files downloaded via `download_cov_inputs.sh`
   - Format: One per sample (10 total)
   - Used as tabix-backed input to methylKit

### Provenance and Reproducibility

- **Legacy results checksums**: `provenance/legacy-results-a910a18.sha256`
- **PE report checksums**: `provenance/pe-reports.sha256`
- **Canonical metadata**: `metadata/wgbs_samples.csv` (maps seq_id to treatment)
- **Configuration file**: `config/analysis_config.json` (documents filtering parameters, model settings, and input/output paths)

All analysis code that uses these preprocessing outputs must:
1. Accept preprocessing file paths or download locations as explicit arguments or configuration
2. Validate checksums or checksums recorded externally when available
3. Never hardcode assumptions about file location, sample order, or treatment assignment

### Missing Samples

**zr3534_7** (T12-C triploid sample) PE report was initially absent from the repository. It was recovered on 2026-08-30 from:

- **URL**: `https://gannet.fish.washington.edu/panopea/030521-ronrosM/zr3534_7_R1.fastp-trim.20201202_bismark_bt2_PE_report.txt`
- **SHA-256**: `e858d1b34f890323f1d8370214512c10ff5a595d326b497838a3cf264cbec866`

All ten samples are now present and validated.

## Downstream Analysis

Preprocessing outputs feed into:

1. **Phase 2** (current): PE report summary and validation
2. **Phase 3**: DML calling with methylKit
3. **Phase 4**: DML output files and exploratory figures
4. **Phase 5**: Genomic feature annotation

See `code/README.md` for a complete analysis workflow.
