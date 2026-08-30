# Legacy DML outputs

The data files in this directory predate the methylation remediation effort and are frozen in `../provenance/legacy-results-a910a18.sha256`.

Do not treat the current CSV and BED files as interchangeable canonical results. In particular, the Cov10 50% CSV contains 1,083 loci while the same-named BED contains 3,804 loci, and their coordinate sets differ. New DML files will be generated together from one clean analysis object in a run-specific results directory.
