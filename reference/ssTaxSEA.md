# Single-Sample Taxon Set Enrichment Analysis

Computes a per-sample score for each taxon set as the mean centered
log-ratio (CLR) of the set's members within that sample.

## Usage

``` r
ssTaxSEA(
  counts,
  lookup_missing = FALSE,
  min_set_size = 5,
  max_set_size = 300,
  custom_db = NULL,
  pseudocount = 0.5
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

- pseudocount:

  Numeric value added to every count before the log transform, to handle
  zeros. Default is 0.5. This choice affects the scores, so report it
  alongside your results.

## Value

A numeric matrix of scores with taxon sets as rows and samples as
columns. Higher values indicate that the set's members are more abundant
in that sample, relative to the average taxon in that same sample.

## Details

For each sample the counts are centered log-ratio transformed, and the
score for a taxon set is the mean of the CLR values of the set's
members:

\$\$score\_{S,j} = \frac{1}{\|S\|} \sum\_{i \in S} clr\_{ij}\$\$

The CLR is computed across *all* taxa supplied in `counts`, before any
subsetting to set members. This matters: restricting the matrix to set
members first would make the geometric mean that the CLR divides by
depend on the sets being tested, which manufactures apparent signal in
sets that have none.

Because the score is a within-sample quantity, `ssTaxSEA` does not need
a cohort and is well defined for a single sample.

## Interpreting the scores

Scores are comparable **across samples within a taxon set**: a higher
score in sample A than sample B means the set's members make up more of
sample A's community.

Scores are **not** comparable across taxon sets. A set of abundant taxa
will score higher than a set of rare taxa in every sample, regardless of
biology, simply because its members are more abundant. When visualising
the matrix, center or scale the rows first (for example
`t(scale(t(scores)))`), and when comparing groups, compare a single
set's scores between groups rather than comparing different sets to each
other.

No p-values are returned. The score is a descriptive statistic; to test
a hypothesis, compare scores across samples using a test appropriate to
your design (for example
[`wilcox.test`](https://rdrr.io/r/stats/wilcox.test.html) or a linear
model on one row of the returned matrix).

## See also

[`TaxSEA`](https://feargalr.github.io/TaxSEA/reference/TaxSEA.md) for
group-level enrichment.

## Examples

``` r
# Toy count matrix: 30 taxa x 8 samples
set.seed(42)
counts <- matrix(rpois(240, lambda = 10), nrow = 30, ncol = 8)
rownames(counts) <- paste0("Taxon_", seq_len(30))
colnames(counts) <- paste0("Sample_", seq_len(8))

# Raise one set's members in the last four samples
counts[1:6, 5:8] <- counts[1:6, 5:8] * 5

sets <- list(
  elevated_set = paste0("Taxon_", 1:6),
  control_set  = paste0("Taxon_", 20:27)
)

scores <- ssTaxSEA(counts, custom_db = sets, min_set_size = 3)
dim(scores)        # sets x samples
#> [1] 2 8
round(scores, 2)
#>              Sample_1 Sample_2 Sample_3 Sample_4 Sample_5 Sample_6 Sample_7
#> elevated_set     0.02     0.20     0.04    -0.13     1.15     1.37     1.31
#> control_set     -0.09     0.02     0.01    -0.01    -0.29    -0.34    -0.52
#>              Sample_8
#> elevated_set     1.07
#> control_set     -0.14

# The planted signal shows up as higher scores in samples 5-8
rowMeans(scores[, 1:4]) - rowMeans(scores[, 5:8])
#> elevated_set  control_set 
#>   -1.1948884    0.3033688 
```
