# TaxSEA News

# TaxSEA 1.5.3

## New

- `TaxSEA()` gains a `bugsigdb` argument (default TRUE). Set it to FALSE to
  skip the run-time BugSigDB download, for offline or reproducible analyses.

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
- Documented that including BugSigDB changes p-values in every other output
  category, because its signatures enlarge the background the competitive KS
  test runs against.

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
