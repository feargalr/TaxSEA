# TaxSEA: Taxon Set Enrichment Analysis

TaxSEA enables rapid annotation of changes by testing for enrichment of
pre-defined taxon sets.

## Usage

``` r
TaxSEA(
  taxon_ranks,
  lookup_missing = FALSE,
  min_set_size = 5,
  max_set_size = 100,
  custom_db = NULL
)
```

## Arguments

- taxon_ranks:

  A named vector of log2 fold changes between control and test groups.

- lookup_missing:

  Logical indicating whether to fetch missing NCBI IDs. Default is
  FALSE.

- min_set_size:

  Minimum size of taxon sets to include in the analysis. Default is 5.

- max_set_size:

  Maximum size of taxon sets to include in the analysis. Default is 100.

- custom_db:

  A user-provided list of taxon sets. If NULL (default), the built-in
  database is used.

## Value

A list of data frames with taxon set enrichment results.

## See also

- <https://doi.org/10.1093/nar/gkac868> for MiMeDB

- <https://doi.org/10.1093/nar/gkab1019> for GMrepo

- <https://doi.org/10.1093/nar/gkab786> for gutMGene

- <https://doi.org/10.1038/s41587-023-01872-y> for BugSigDB

## Examples

``` r
data("TaxSEA_test_data")
taxsea_results <- TaxSEA(TaxSEA_test_data)
#> 
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
```
