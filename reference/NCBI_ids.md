# NCBI IDs Dataset

A dataset for mapping NCBI IDs to species/genus names. This named vector
allows for lookup of NCBI IDs associated with species or genus names.

## Usage

``` r
NCBI_ids
```

## Format

A named vector where:

- names:

  NCBI IDs

- values:

  Species or genus names

## Source

NCBI

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
