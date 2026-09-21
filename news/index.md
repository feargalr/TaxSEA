# Changelog

## TaxSEA 1.5.6

- Rewrote the package description: it now leads with what TaxSEA is for
  and describes both approaches,
  [`TaxSEA()`](https://feargalr.github.io/TaxSEA/reference/TaxSEA.md) on
  ranked taxa and
  [`ssTaxSEA()`](https://feargalr.github.io/TaxSEA/reference/ssTaxSEA.md)
  on a single sample. Also fixed two long-standing typos in it
  (“alternations”, and a missing “can be applied”).

## TaxSEA 1.5.5

### Breaking changes

- The BacDive sets in `TaxSEA_db` have been rebuilt from a fresh,
  complete harvest of BacDive: **1,776 sets covering 9,045 taxa**, up
  from 267 sets covering 716. Results for `BacDive_bacterial_physiology`
  will differ from earlier releases, and code that refers to BacDive
  sets by name needs updating.

  - **Set names changed.** Names no longer contain spaces and follow
    `BacDive_<Family>_<value>`. For example:

    | Before | Now |
    |----|----|
    | `BacDive_facultative anaerobe` | `BacDive_Oxygen_facultative_anaerobe` |
    | `BacDive_anaerobe`, `BacDive_aerobe`, `BacDive_microaerophile` | `BacDive_Oxygen_anaerobe`, `BacDive_Oxygen_aerobe`, `BacDive_Oxygen_microaerophile` |
    | `BacDive_rod-shaped`, `BacDive_coccus-shaped` | `BacDive_Shape_rod`, `BacDive_Shape_coccus` (and `BacDive_Shape_coccoid`, coccus or ovoid) |
    | `BacDive_Utilizes_<substrate>` | `BacDive_Uses_<substrate>`, split by assay (below) |
    | `BacDive_Enzyme_alkaline phosphatase` | `BacDive_Enzyme_alkaline_phosphatase` |

  - **Memberships changed.** Sets now contain measured BacDive data
    only, aggregated over every strain of each species. Membership that
    came from the earlier LLM-derived table is gone. For example,
    *Faecalibacterium prausnitzii* is no longer in the anaerobe set,
    because BacDive holds no oxygen tolerance record for it.

  - **Substrate use is split by assay.** Growing on a substrate,
    fermenting it, producing acid or gas from it, reducing it, respiring
    it and using it as a nitrogen source are different physiology and
    are now separate sets. The former `BacDive_Utilizes_nitrate` mixed
    several of these; the nearest equivalent is
    `BacDive_Uses_nitrate_reduction`.

  - Because TaxSEA compares each set against all other taxa covered by
    the sets being analysed, the larger BacDive collection and the
    extended `NCBI_ids` also shift p-values slightly for other sources.
    On the bundled test data, the top-ranked metabolite-producer and
    health-association sets are unchanged, with small movements in
    p-value and order.

### New features

- New BacDive families: Gram stain, motility, spore formation, cell
  shape, growth temperature class, salt tolerance, pH range, GC content,
  enzyme activities, metabolite production, respiratory quinones,
  peptidoglycan type and more. Sets of species tested and found
  **negative** are included (suffix `_negative`), which earlier releases
  could not express.
- `NCBI_ids` extended with the species names of every BacDive set
  member, with spaces and with underscores, including former names of
  reclassified species (e.g. `Bacteroides vulgatus` and
  `Phocaeicola vulgatus` both map to 821), and each member’s taxid as
  its own name. No existing entry was changed. Names that were present
  with no taxid, and so never mapped, were filled where BacDive resolves
  them. This also recovers members of other sources: taxa in `TaxSEA_db`
  with no `NCBI_ids` entry fell from 84 to 57.
- Provenance for the BacDive sets ships in `inst/extdata`:
  `BacDive_set_provenance.tsv` (per set: source family, analyte, number
  of species and strains) and `BacDive_build_info.tsv` (harvester
  release, BacDive client version, harvest dates). How the sets were
  built is documented in `inst/scripts/make_BacDive_sets.R`; the
  pipeline is <https://github.com/feargalr/bacdive_harvester>.

## TaxSEA 1.5.4

### Breaking changes

- [`TaxSEA()`](https://feargalr.github.io/TaxSEA/reference/TaxSEA.md) no
  longer includes BugSigDB by default; use `bugsigdb = TRUE` to include
  it. Previously it was included automatically whenever `bugsigdbr` was
  installed.

  Why this matters: TaxSEA tests each taxon set by comparing its members
  against all the other taxa covered by the sets being analysed.
  BugSigDB adds thousands of signatures covering many extra taxa, so
  switching it on changes that comparison for every set, not just the
  BugSigDB ones. The same metabolite-producer set can get a different
  p-value depending on whether BugSigDB was included, and on which
  BugSigDB release was downloaded; on the bundled test data some
  p-values moved by as much as 0.26. Leaving it off by default keeps
  results reproducible and lets
  [`TaxSEA()`](https://feargalr.github.io/TaxSEA/reference/TaxSEA.md)
  run offline.

### Bug fixes

- If BugSigDB cannot be downloaded (for example during a Zenodo outage),
  `TaxSEA(bugsigdb = TRUE)` now warns and continues without it instead
  of failing. The vignette and tests likewise skip their
  BugSigDB-dependent parts, so an upstream outage no longer breaks the
  package build.
- The main vignette and README accessed BugSigDB results as
  `taxsea_results$BugSigdB`, which silently returns `NULL`; corrected to
  `$BugSigDB`. The vignette also documented a nonexistent `database`
  argument.
- [`get_taxon_sets()`](https://feargalr.github.io/TaxSEA/reference/get_taxon_sets.md)
  had a default argument referring to a nonexistent object, so calling
  it with no argument failed with `object 'taxon' not found`. It now
  reports the missing argument properly. It also no longer runs
  unreachable code after
  [`return()`](https://rdrr.io/r/base/function.html), and correctly
  handles a lookup that resolves to more than one NCBI ID.

### Other

- Raised the R dependency to 4.6.0 to match Bioconductor 3.24.
- [`get_ncbi_taxon_ids()`](https://feargalr.github.io/TaxSEA/reference/get_ncbi_taxon_ids.md)
  uses [`vapply()`](https://rdrr.io/r/base/lapply.html) instead of
  [`sapply()`](https://rdrr.io/r/base/lapply.html), and qualifies its
  `utils` calls.
- [`taxon_rank_sets()`](https://feargalr.github.io/TaxSEA/reference/taxon_rank_sets.md)
  examples use `\donttest` rather than `\dontrun`.

## TaxSEA 1.5.3

### New

- [`TaxSEA()`](https://feargalr.github.io/TaxSEA/reference/TaxSEA.md)
  gains a `bugsigdb` argument controlling whether BugSigDB signatures
  are downloaded and included at run time.

### Testing

- Added unit tests for
  [`ssTaxSEA()`](https://feargalr.github.io/TaxSEA/reference/ssTaxSEA.md),
  for ORA mode and for
  [`TaxSEA()`](https://feargalr.github.io/TaxSEA/reference/TaxSEA.md)
  input validation, none of which were previously covered.
- Added regression snapshots of the results produced from the bundled
  test data, so that database updates surface as a reviewable diff.
- Added database integrity checks covering set naming, set sizes and
  `NCBI_ids` coverage. These document two pre-existing issues: 84 set
  members have no `NCBI_ids` entry and are silently dropped (worst
  affected are `Siderophore_producers`, 73 members of which 25 are
  usable, and the Valles-Colomer2019 Gut-Brain Modules), and
  `GutMGene_producers_of_Phenylalanine` appears twice, so two of its
  taxa are unreachable.

### Documentation

- Corrected the documented default for `max_set_size` in
  [`TaxSEA()`](https://feargalr.github.io/TaxSEA/reference/TaxSEA.md),
  which said 100 while the signature has been 300.
- Documented that including BugSigDB changes the p-values of all taxon
  sets, not only the BugSigDB ones.

## TaxSEA 1.5.2

### Breaking changes

- [`ssTaxSEA()`](https://feargalr.github.io/TaxSEA/reference/ssTaxSEA.md)
  now scores each taxon set as the **mean centered log-ratio (CLR) of
  the set’s members within a sample**, replacing the previous ranked,
  cohort z-scored, ssGSEA-style running-sum statistic.
- [`ssTaxSEA()`](https://feargalr.github.io/TaxSEA/reference/ssTaxSEA.md)
  returns a **numeric matrix with taxon sets as rows and samples as
  columns**, instead of a list of two matrices oriented samples x sets.
  Code written as `res$scores[sample, set]` becomes `res[set, sample]`.
- [`ssTaxSEA()`](https://feargalr.github.io/TaxSEA/reference/ssTaxSEA.md)
  no longer returns p-values. The score is a descriptive statistic; test
  it across samples with a test appropriate to your design.

### New

- [`ssTaxSEA()`](https://feargalr.github.io/TaxSEA/reference/ssTaxSEA.md)
  gains a `pseudocount` argument (default 0.5), added to every value
  before the log transform rather than only to zeros.
- [`ssTaxSEA()`](https://feargalr.github.io/TaxSEA/reference/ssTaxSEA.md)
  now works on a **single sample**, because the score no longer
  references the rest of the cohort.
- [`ssTaxSEA()`](https://feargalr.github.io/TaxSEA/reference/ssTaxSEA.md)
  rejects proportion-like input (columns summing to 1), negative values,
  non-finite values and empty samples, instead of silently producing
  meaningless scores.

### Bug fixes

- The CLR is now computed across all supplied taxa before subsetting to
  set members. Previously the matrix was filtered to set members first,
  which made the geometric mean the CLR divides by depend on which sets
  were being tested and could manufacture apparent signal in sets that
  had none.

## TaxSEA 1.5.1

- Removed `R/TaxSEA_export.R`, an outdated duplicate of an internal
  function that shadowed the current version.
- Declared the `methods` import and the remaining `stats`/`utils`
  imports, clearing the Bioconductor devel check warning.
- Recompressed the bundled data (`NCBI_ids` 136 Kb to 42 Kb, `TaxSEA_db`
  48 Kb to 35 Kb), clearing the data compression warning.
- Added a GitHub Actions workflow running `R CMD check` and `BiocCheck`
  against Bioconductor devel.

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
