# TaxSEA News

# TaxSEA 1.5.4

## Breaking changes

- `TaxSEA()` no longer includes BugSigDB by default; use `bugsigdb = TRUE`
  to include it. Previously it was included automatically whenever
  `bugsigdbr` was installed.

  Why this matters: TaxSEA tests each taxon set by comparing its members
  against all the other taxa covered by the sets being analysed. BugSigDB
  adds thousands of signatures covering many extra taxa, so switching it on
  changes that comparison for every set, not just the BugSigDB ones. The
  same metabolite-producer set can get a different p-value depending on
  whether BugSigDB was included, and on which BugSigDB release was
  downloaded; on the bundled test data some p-values moved by as much as
  0.26. Leaving it off by default keeps results reproducible and lets
  `TaxSEA()` run offline.

## Bug fixes

- If BugSigDB cannot be downloaded (for example during a Zenodo outage),
  `TaxSEA(bugsigdb = TRUE)` now warns and continues without it instead of
  failing. The vignette and tests likewise skip their BugSigDB-dependent
  parts, so an upstream outage no longer breaks the package build.
- The main vignette and README accessed BugSigDB results as
  `taxsea_results$BugSigdB`, which silently returns `NULL`; corrected to
  `$BugSigDB`. The vignette also documented a nonexistent `database`
  argument.
- `get_taxon_sets()` had a default argument referring to a nonexistent object,
  so calling it with no argument failed with `object 'taxon' not found`. It
  now reports the missing argument properly. It also no longer runs
  unreachable code after `return()`, and correctly handles a lookup that
  resolves to more than one NCBI ID.

## Other

- Raised the R dependency to 4.6.0 to match Bioconductor 3.24.
- `get_ncbi_taxon_ids()` uses `vapply()` instead of `sapply()`, and qualifies
  its `utils` calls.
- `taxon_rank_sets()` examples use `\donttest` rather than `\dontrun`.

# TaxSEA 1.5.3

## New

- `TaxSEA()` gains a `bugsigdb` argument controlling whether BugSigDB
  signatures are downloaded and included at run time.

## Testing

- Added unit tests for `ssTaxSEA()`, for ORA mode and for `TaxSEA()` input
  validation, none of which were previously covered.
- Added regression snapshots of the results produced from the bundled test
  data, so that database updates surface as a reviewable diff.
- Added database integrity checks covering set naming, set sizes and
  `NCBI_ids` coverage. These document two pre-existing issues: 84 set members
  have no `NCBI_ids` entry and are silently dropped (worst affected are
  `Siderophore_producers`, 73 members of which 25 are usable, and the
  Valles-Colomer2019 Gut-Brain Modules), and
  `GutMGene_producers_of_Phenylalanine` appears twice, so two of its taxa
  are unreachable.

## Documentation

- Corrected the documented default for `max_set_size` in `TaxSEA()`, which
  said 100 while the signature has been 300.
- Documented that including BugSigDB changes the p-values of all taxon sets,
  not only the BugSigDB ones.

# TaxSEA 1.5.2

## Breaking changes

- `ssTaxSEA()` now scores each taxon set as the **mean centered log-ratio
  (CLR) of the set's members within a sample**, replacing the previous
  ranked, cohort z-scored, ssGSEA-style running-sum statistic.
- `ssTaxSEA()` returns a **numeric matrix with taxon sets as rows and samples
  as columns**, instead of a list of two matrices oriented samples x sets.
  Code written as `res$scores[sample, set]` becomes `res[set, sample]`.
- `ssTaxSEA()` no longer returns p-values. The score is a descriptive
  statistic; test it across samples with a test appropriate to your design.

## New

- `ssTaxSEA()` gains a `pseudocount` argument (default 0.5), added to every
  value before the log transform rather than only to zeros.
- `ssTaxSEA()` now works on a **single sample**, because the score no longer
  references the rest of the cohort.
- `ssTaxSEA()` rejects proportion-like input (columns summing to 1), negative
  values, non-finite values and empty samples, instead of silently producing
  meaningless scores.

## Bug fixes

- The CLR is now computed across all supplied taxa before subsetting to set
  members. Previously the matrix was filtered to set members first, which made
  the geometric mean the CLR divides by depend on which sets were being
  tested and could manufacture apparent signal in sets that had none.

# TaxSEA 1.5.1

- Removed `R/TaxSEA_export.R`, an outdated duplicate of an internal
  function that shadowed the current version.
- Declared the `methods` import and the remaining `stats`/`utils` imports,
  clearing the Bioconductor devel check warning.
- Recompressed the bundled data (`NCBI_ids` 136 Kb to 42 Kb, `TaxSEA_db`
  48 Kb to 35 Kb), clearing the data compression warning.
- Added a GitHub Actions workflow running `R CMD check` and `BiocCheck`
  against Bioconductor devel.

# TaxSEA 1.3.3

- Add `ssTaxSEA()` function for single sample enrichment testing and documentation

# TaxSEA 1.3.2

- Add `taxon_rank_sets()` helper with documentation and tests
- Add functionality to interact with tse objects

## TaxSEA 1.3.1

- Internal refactor to modular analysis pipeline 
- Added infrastructure for over-representation analysis (ORA).
- Improved internal input preparation and validation.

## Version 1.1.8 (2025-08-12)
- Added output for all taxon sets 
- Fixed spelling in GutMGene_producers_of_Phenyalanine

## Version 1.1.7 (2025-08-11)  
- Added Gut-Brain modules database from Valles-Colomer et al.
Nature Microbiology. 2019. 

## Version 1.1.6 (2025-08-08)  
- Bug fixes

## Version 1.1.5 (2025-08-06)
- Fixed bug in BacDive Database

## Version 1.1.4 (2025-08-05)
- Added BacDive Database
- Added Blossum and VANISH taxon sets in disease associations category

## Version 0.99.2 (2025-02-20)
- Added ability to test custom taxon sets using `custom_db` parameter.
- Removed deprecated `plotting.R` file.
- TaxSEA now reports the KS test statistic for each taxon set.

## Version 0.99.1 (2025-02-11)
- Bringing into line with requirements for BioConductor submission
