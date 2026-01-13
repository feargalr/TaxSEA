# TaxSEA Database A dataset containing taxon sets. Each item in the list is a taxon set, and each member within a taxon set is a taxon.

TaxSEA Database A dataset containing taxon sets. Each item in the list
is a taxon set, and each member within a taxon set is a taxon.

## Usage

``` r
TaxSEA_db
```

## Format

A list of vectors. Each vector contains character strings representing
taxa.

## Source

See READ ME.

## Examples

``` r
data(TaxSEA_db)
all_sets <- names(TaxSEA_db)
GABA_producers<-TaxSEA_db[["MiMeDB_producers_of_GABA"]]
```
