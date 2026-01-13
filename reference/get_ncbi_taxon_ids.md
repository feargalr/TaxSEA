# Retrieve NCBI Taxonomy IDs for a list of taxon names

This function takes a vector of taxon names and returns a vector of NCBI
taxonomy IDs by querying the NCBI Entrez API.

## Usage

``` r
get_ncbi_taxon_ids(taxon_names)
```

## Arguments

- taxon_names:

  A character vector of taxon names

## Value

A character vector of NCBI taxonomy IDs corresponding to the input taxon
names

## Examples

``` r
taxon_names <- c("Escherichia coli", "Staphylococcus aureus")
taxon_ids <- get_ncbi_taxon_ids(taxon_names)
```
