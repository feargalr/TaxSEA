# TaxSEA Test Data

A dataset containing taxon ranks and taxon IDs.

## Usage

``` r
TaxSEA_test_data
```

## Format

A data frame with two columns:

- rank:

  Character vector representing taxon ranks

- id:

  Character vector representing taxon IDs

## Source

See READ ME.

## Examples

``` r
data(TaxSEA_test_data)
test_results <- TaxSEA(TaxSEA_test_data)
#> Using cached version from 2026-01-13 06:16:43
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
