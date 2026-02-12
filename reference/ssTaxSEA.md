# Single-Sample Taxon Set Enrichment Analysis

Computes per-sample enrichment scores for taxon sets using an
ssGSEA-style approach. Counts are CLR-transformed, then each taxon is
z-scored across samples to capture between-sample variation. Per-sample
enrichment is then computed by ranking each sample's z-scores and
applying a weighted running-sum statistic against taxon sets from the
TaxSEA database.

## Usage

``` r
ssTaxSEA(
  counts,
  lookup_missing = FALSE,
  min_set_size = 5,
  max_set_size = 300,
  custom_db = NULL
)
```

## Arguments

- counts:

  A numeric matrix, data.frame, or
  `SummarizedExperiment`/`TreeSummarizedExperiment` object. For
  matrix/data.frame input: rows are taxa, columns are samples, and row
  names must be taxon names (e.g. species names or NCBI IDs). For
  SummarizedExperiment input: the first assay is used and
  [`rownames()`](https://rdrr.io/r/base/colnames.html) provide taxon
  identifiers.

- lookup_missing:

  Logical indicating whether to fetch missing NCBI IDs via the NCBI API.
  Default is FALSE.

- min_set_size:

  Minimum size of taxon sets to include. Default is 5.

- max_set_size:

  Maximum size of taxon sets to include. Default is 300.

- custom_db:

  A user-provided list of taxon sets. If NULL (default), the built-in
  TaxSEA database is used (excluding BugSigDB).

## Value

A list with two elements:

- scores:

  A matrix (samples x taxon sets) of enrichment scores. Positive scores
  indicate the set taxa tend to have higher abundance in that sample
  relative to the cohort.

- pvalues:

  A matrix (samples x taxon sets) of KS test p-values for each
  sample-set combination.

## Details

The approach works as follows:

1.  Raw counts are CLR-transformed (centered log-ratio with pseudocount
    of 0.5 for zeros).

2.  Each taxon is then z-scored across all samples, so that values
    represent how much higher or lower a taxon is in a given sample
    relative to the cohort mean.

3.  For each sample, the z-scores are ranked and an ssGSEA-style
    weighted running-sum enrichment score is computed for each taxon
    set.

4.  A KS test p-value is also computed per sample per set.

This cohort-relative approach ensures that taxon sets which are
consistently elevated in a subset of samples (e.g. disease samples) will
produce high enrichment scores in those samples, even if the taxa are
not the most abundant within any single sample.

## Examples

``` r
if (FALSE) { # \dontrun{
# From a count matrix (taxa x samples)
counts <- matrix(rpois(500, lambda = 10), nrow = 50, ncol = 10)
rownames(counts) <- paste0("Taxon_", seq_len(50))
colnames(counts) <- paste0("Sample_", seq_len(10))
res <- ssTaxSEA(counts, custom_db = list(
  set1 = paste0("Taxon_", 1:10),
  set2 = paste0("Taxon_", 20:30)
), min_set_size = 2)
head(res$scores)
head(res$pvalues)
} # }
```
