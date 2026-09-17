# NCBI IDs Dataset

A lookup from taxon names to NCBI taxonomy IDs, used to translate the
names supplied to
[`TaxSEA()`](https://feargalr.github.io/TaxSEA/reference/TaxSEA.md) into
the IDs that `TaxSEA_db` is built on. Names are present both with spaces
and with underscores, and each ID is also present as its own name so
that IDs can be supplied directly. Species names, former names and IDs
of the `BacDive_` set members were added from the BacDive harvest
without changing any existing entry.

## Usage

``` r
data(NCBI_ids)
```

## Format

A named list where:

- names:

  Taxon names (e.g. "Bifidobacterium breve", "Bifidobacterium_breve") or
  NCBI IDs

- values:

  NCBI taxonomy IDs, as character

## Source

NCBI Taxonomy; BacDive (<https://bacdive.dsmz.de>).

## Examples

``` r
data(NCBI_ids)
# Can look up either with or without spaces
NCBI_ids["Bifidobacterium_breve"]
#> $Bifidobacterium_breve
#> [1] "1685"
#> 
NCBI_ids["Bifidobacterium breve"]
#> $`Bifidobacterium breve`
#> [1] "1685"
#> 
```
