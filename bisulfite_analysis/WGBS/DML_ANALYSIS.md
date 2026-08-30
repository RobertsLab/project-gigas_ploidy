# Remediated MethylKit DML analysis

## Scope and status

`code/2_WGBS_Methylkit.R` is a clean-session, command-line replacement for the legacy interactive script. It contains no RStudio integration, saved-workspace loading, user-specific path, mutable working directory, or dependence on an object from a previous session.

The workflow and its preflight checks are implemented, but the canonical full result remains pending until the ten external COV files are downloaded and the analysis is run in the pinned environment. Existing files under `DML/` remain legacy and must not be promoted based on this code-only remediation.

## External input contract

The source directory is:

`https://gannet.fish.washington.edu/panopea/WGBS-gigas-ploidy-desiccation/cov_files/`

The ten required files use this exact template:

`{seq_id}_R1.CpG_report.merged_CpG_evidence.cov`

The files published for `zr3534_1` through `zr3534_10` are approximately 428–446 MB each. This naming differs from the `.fastp-trim.20201202` filenames hard-coded in the legacy R script. Download or resume all ten files and create a local SHA-256 manifest with:

```bash
bisulfite_analysis/WGBS/code/download_cov_inputs.sh /path/to/wgbs-cov
```

Each COV row has six tab-separated fields: chromosome, start, end, percent methylation, methylated count, and unmethylated count. The upstream `coverage2cytosine --merge_CpG --zero_based` command produced zero-based, half-open two-base CpG intervals. The remediated workflow validates that boundary convention and preserves it in BED outputs without coordinate conversion.

## Primary model

The predefined primary estimand is the ploidy effect adjusted for heat-shock status:

- treatment 0: diploid;
- treatment 1: triploid;
- covariate: heat shock (`no=0`, `yes=1`);
- overdispersion: McCullagh-Nelder (`MN`);
- statistical test: F test;
- multiple-testing adjustment: SLIM;
- effect calculation: predicted group means;
- significance threshold: q ≤ 0.01;
- reported effect subsets: absolute methylation difference ≥20 and ≥50 percentage points.

Positive `meth.diff` means higher predicted methylation in triploids than diploids after accounting for heat shock. The ten animals remain independent biological replicates; samples are never pooled.

The secondary ploidy-by-heat interaction is intentionally not implemented. It requires a model/interface selected specifically for an interaction coefficient and must not be inferred by comparing primary-model DML lists.

## Filtering rationale

- `methRead` minimum coverage 2 removes essentially unsupported calls while retaining data for explicit downstream filtering.
- The primary analysis requires coverage ≥10 in every sample so every tested CpG has adequate counts and a complete design matrix.
- The upper 99.9th coverage percentile is excluded per sample to reduce influence from abnormally high coverage that may reflect repeats, mapping artifacts, or residual duplication.
- Median coverage normalization is enabled to reduce sample-level sequencing-depth imbalance before count modeling.
- Requiring a CpG in all ten samples avoids missing-value-dependent changes to the fitted design across loci.

These choices define the primary analysis; they are not proof of robustness. Planned sensitivity analyses must compare the 10x/all-sample result with predefined alternatives such as 5x coverage, relaxed occupancy, and BH adjustment before biological conclusions are finalized.

## Commands

Validate local COV filenames, schema, coordinates, counts, sample mappings, and model parameters without loading methylKit:

```bash
export GIGAS_WGBS_COV_DIR=/path/to/wgbs-cov
bisulfite_analysis/WGBS/code/run_methylkit_dml.sh --preflight-only
```

Build the pinned environment, then run the full workflow:

```bash
docker build -t gigas-wgbs:bioc-3.23 bisulfite_analysis/WGBS/environment
docker run --rm \
  -v "$PWD:/project-gigas_ploidy" \
  -v "$GIGAS_WGBS_COV_DIR:/wgbs-cov:ro" \
  -e GIGAS_WGBS_COV_DIR=/wgbs-cov \
  -w /project-gigas_ploidy \
  gigas-wgbs:bioc-3.23 \
  bisulfite_analysis/WGBS/code/run_methylkit_dml.sh
```

Every execution uses a new timestamped directory. The script writes input and metadata snapshots, parameters, a complete unfiltered statistics CSV, paired CSV/BED threshold subsets, output counts, methylKit flat-file databases, and `sessionInfo()`. It refuses to write into a non-empty output directory.

The post-run validator checks coordinate uniqueness, `meth.diff` bounds, q/effect thresholds, exact CSV/BED coordinate equality, subset membership in the full table, and declared row counts.
