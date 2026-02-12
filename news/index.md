# Changelog

## TaxSEA 1.3.3

- Add
  [`ssTaxSEA()`](https://feargalr.github.io/TaxSEA/reference/ssTaxSEA.md)
  function for single sample enrichment testing and documentation

## TaxSEA 1.3.2

- Add
  [`taxon_rank_sets()`](https://feargalr.github.io/TaxSEA/reference/taxon_rank_sets.md)
  helper with documentation and tests
- Add functionality to interact with tse objects

### TaxSEA 1.3.1

- Internal refactor to modular analysis pipeline
- Added infrastructure for over-representation analysis (ORA).
- Improved internal input preparation and validation.

### Version 1.1.8 (2025-08-12)

- Added output for all taxon sets
- Fixed spelling in GutMGene_producers_of_Phenyalanine

### Version 1.1.7 (2025-08-11)

- Added Gut-Brain modules database from Valles-Colomer et al. Nature
  Microbiology. 2019.

### Version 1.1.6 (2025-08-08)

- Bug fixes

### Version 1.1.5 (2025-08-06)

- Fixed bug in BacDive Database

### Version 1.1.4 (2025-08-05)

- Added BacDive Database
- Added Blossum and VANISH taxon sets in disease associations category

### Version 0.99.2 (2025-02-20)

- Added ability to test custom taxon sets using `custom_db` parameter.
- Removed deprecated `plotting.R` file.
- TaxSEA now reports the KS test statistic for each taxon set.

### Version 0.99.1 (2025-02-11)

- Bringing into line with requirements for BioConductor submission
