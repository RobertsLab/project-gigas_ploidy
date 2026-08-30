# WGBS configuration

`analysis_config.json` is the draft configuration for the clean WGBS rerun. Paths are relative to the repository root. Large COV inputs remain external; set `GIGAS_WGBS_COV_DIR` to the local directory containing them.

The primary model is explicitly identified as a draft until its design receives scientific review. The secondary interaction model is intentionally unset so that code cannot imply an unsupported interaction analysis.

`allowed_missing_pe_reports` is empty because all ten Bismark PE reports are now committed. The previously absent `zr3534_7` report was recovered from the project archive; its source and checksum are recorded in `../data/PE_reports/README.md`.

The COV filename template matches the files currently published at the configured gannet source. Those files omit the `.fastp-trim.20201202` segment used in the legacy script's hard-coded local paths.

The primary model treats diploid as treatment 0 and triploid as treatment 1, adjusts for heat-shock status, uses methylKit's predicted effect, McCullagh-Nelder overdispersion correction, F test, and SLIM q-values. Positive `meth.diff` therefore means higher predicted methylation in triploids relative to diploids after adjustment for heat shock.
