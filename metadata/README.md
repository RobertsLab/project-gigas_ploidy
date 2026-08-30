# Project metadata

`wgbs_samples.csv` is the canonical sample metadata source for the WGBS analysis. Analysis scripts must join metadata by `seq_id`; they must not infer treatment or ploidy from file order or sample-name prefixes.

The file supersedes the historical `gene_expression/data/zr3534_wgbs_info.csv`, whose field names included a misspelling of desiccation. The historical values were migrated without changing the sample-to-SRA mapping.

Required fields:

- `seq_id`: sequencing identifier used in Bismark and coverage filenames.
- `library_name`: Roberts Lab biological sample identifier.
- `tissue`: sampled tissue.
- `ploidy`: `diploid` or `triploid`.
- `desiccation`: whether the animal received the 24-hour desiccation treatment.
- `heat_shock`: whether the animal received the subsequent acute heat shock.
- `library_kit`: WGBS library preparation kit.
- `sra_bioproject`: SRA BioProject accession.
- `sra_accession`: sample SRA accession.

Changes to sample assignments require a documented source and review because they affect every downstream contrast.
