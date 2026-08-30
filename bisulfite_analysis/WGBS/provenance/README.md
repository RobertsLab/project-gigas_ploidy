# WGBS result provenance

`legacy-results-a910a18.sha256` freezes the pre-remediation result files as they existed at commit `a910a18a91b6733b52f3ba7337e35a889d491136` on 2026-08-30.

The manifest covers 59 files under:

- `bisulfite_analysis/WGBS/DML/`
- `bisulfite_analysis/WGBS/analyses/`
- `bisulfite_analysis/WGBS/plots/`

These checksums preserve provenance; they do not validate the scientific correctness of the files. See `METHYLATION_REMEDIATION_PLAN.md` and the README in each result directory for known limitations.

From the repository root, verify the frozen files with:

```bash
shasum -a 256 -c bisulfite_analysis/WGBS/provenance/legacy-results-a910a18.sha256
```

Remediated outputs must be generated in a new run-specific directory and must not overwrite these files during development.

`pe-reports.sha256` records the complete ten-report source set used by the remediated PE-report summary. It includes the recovered `zr3534_7` report. Verify it from the repository root with:

```bash
shasum -a 256 -c bisulfite_analysis/WGBS/provenance/pe-reports.sha256
```
