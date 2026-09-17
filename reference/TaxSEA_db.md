# TaxSEA Database A dataset containing taxon sets. Each item in the list is a taxon set, and each member within a taxon set is a taxon.

TaxSEA Database A dataset containing taxon sets. Each item in the list
is a taxon set, and each member within a taxon set is a taxon.

## Usage

``` r
data(TaxSEA_db)
```

## Format

A list of vectors. Each vector contains character strings representing
taxa as NCBI taxonomy IDs.

## Source

See the README. BacDive sets: Schober et al. (2025) BacDive in 2025: the
core database for prokaryotic strain data, via
<https://bacdive.dsmz.de>.

## Details

Set names start with their source: `MiMeDB_`, `GutMGene_`, `GMRepoV2_`,
`mBodyMap_`, `BacDive_` and so on.

The `BacDive_` sets describe measured phenotypes of cultured strains
(oxygen tolerance, Gram stain, substrate use, enzyme activities, growth
temperature, salt and pH ranges and others), aggregated to species over
every strain BacDive holds. They are built by the bacdive_harvester
pipeline (<https://github.com/feargalr/bacdive_harvester>); see
`inst/scripts/make_BacDive_sets.R`. Only sets with at least 3 members
are included. Names read `BacDive_<Family>_<value>`, optionally followed
by an assay context (e.g. `_fermentation`, `_reduction`) and `_negative`
for species tested and found negative, e.g.
`BacDive_Oxygen_facultative_anaerobe`, `BacDive_Uses_nitrate_reduction`.

Per-set provenance (source family, number of species and strains behind
each set) and a record of the build (harvester release, BacDive client
version, harvest dates) are shipped with the package:


    system.file("extdata", "BacDive_set_provenance.tsv", package = "TaxSEA")
    system.file("extdata", "BacDive_build_info.tsv", package = "TaxSEA")

## Examples

``` r
data(TaxSEA_db)
all_sets <- names(TaxSEA_db)
GABA_producers<-TaxSEA_db[["MiMeDB_producers_of_GABA"]]
```
