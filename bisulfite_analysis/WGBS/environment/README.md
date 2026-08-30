# Reproducible WGBS environment

The remediated downstream analysis is pinned to R 4.6.0, Bioconductor 3.23, and methylKit 1.38.0. Bioconductor 3.23 is the official release for R 4.6, and methylKit 1.38.0 is the package version in that release:

- https://bioconductor.org/install/
- https://bioconductor.org/packages/3.23/bioc/html/methylKit.html

Build the environment from the repository root:

```bash
docker build -t gigas-wgbs:bioc-3.23 bisulfite_analysis/WGBS/environment
```

The build fails unless the expected R, Bioconductor, and methylKit versions are present. The multi-platform base image is pinned to index digest `sha256:b10002b39efa30c3779ad839549806ebdbb29b3266f0d2428478b04426e55929`, in addition to the explicit Bioconductor/R tag. The exact methylKit source tarball is verified against SHA-256 `008990fdf453e7ec66781ea6736b8f373dfae12d891c76ef3bf3456e97ff264e` before installation. Record the final locally built image digest with each released result manifest.

Run the environment check inside the image with:

```bash
docker run --rm gigas-wgbs:bioc-3.23 Rscript --vanilla /opt/gigas/check_environment.R
```

The local development machine currently has R 4.3.2/Bioconductor 3.18 and cannot execute the canonical full DML analysis. Package-independent input preflight tests remain runnable locally; the full analysis must use the pinned container or an independently verified equivalent environment.

`software_versions.tsv` separates historical upstream preprocessing versions from the remediated downstream environment. The original cluster script records Bismark 0.21.0, Bowtie 2 2.3.4.1, and samtools 1.9. `bedtools` is not part of DML calling and will be pinned when the feature-annotation workflow is remediated.
