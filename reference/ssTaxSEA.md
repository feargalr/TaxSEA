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
  pseudocount = NULL,
  scale = c("counts", "proportion"),
  reference = c("all", "top_k"),
  top_k = 20
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

  Numeric value added to every value before the log transform, to handle
  zeros. It is on the scale set by `scale`: if NULL (default), 0.5 for
  `scale = "counts"` and 1e-6 for `scale = "proportion"`. This choice
  affects the scores, so report it alongside your results.

- scale:

  Either `"counts"` (default) or `"proportion"`. With `"counts"` the
  input is used as supplied and must be counts or another count-like
  measure. With `"proportion"` each sample is first divided by its
  total, so counts, proportions and percentages (for example MetaPhlAn
  output) are all accepted, and `pseudocount` is a proportion. See
  Details.

- reference:

  Either `"all"` (default) or `"top_k"`. What each sample's log values
  are centered on: the mean over all taxa (the standard CLR), or the
  mean over that sample's `top_k` most abundant detected taxa. See
  Details.

- top_k:

  Number of most abundant detected taxa used as the reference when
  `reference = "top_k"`. Default is 20. A sample with fewer detected
  taxa uses all of its detected taxa.

## Value

A numeric matrix of scores with taxon sets as rows and samples as
columns. Higher values indicate that the set's members are more abundant
in that sample, relative to the average taxon in that same sample (or,
with `reference = "top_k"`, relative to that sample's most abundant
taxa).

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

## Choosing a reference

The standard CLR (`reference = "all"`) centers each sample on the mean
log value of every taxon in the table, zeros included. That reference
moves with richness: when many rare taxa drop out, as after antibiotics
or chemotherapy, the reference falls and every surviving taxon's score
rises, even if its abundance relative to the other survivors has not
changed. Many taxon sets share the same surviving members, so hundreds
of scores can rise together. Adjusting for richness in a downstream
model does not remove this. The reference also depends on which rows are
in the table, so filtering taxa differently changes every score.

`reference = "top_k"` centers each sample on the mean log value of its
`top_k` most abundant detected taxa. Rare taxa are still scored and
still count towards their sets; they just no longer move everyone else's
score. Losing rare taxa, or adding or removing rare and all-zero rows,
leaves the reference unchanged. Zeros never enter the reference: a
sample with fewer than `top_k` detected taxa uses all of them. The score
then reads as abundance relative to the sample's dominant taxa.

A pseudocount set from each sample's own data would move with the rare
tail in the same way, which is why `pseudocount` is a fixed value. On
the count scale a pseudocount of 0.5 reads is a different proportion in
samples of different depth; `scale = "proportion"` with a small fixed
pseudocount avoids that. The default of 1e-6 is near the detection limit
of deep shotgun profiles; raise it if your smallest detected proportions
are larger. For relative abundance input, such as MetaPhlAn percentages,
use `scale = "proportion"`.

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

# A reference that is not moved by changes in rare taxa, for
# relative abundance input
pct <- sweep(counts, 2, colSums(counts), "/") * 100
ssTaxSEA(pct, custom_db = sets, min_set_size = 3,
         scale = "proportion", reference = "top_k", top_k = 10)
#>                Sample_1    Sample_2   Sample_3   Sample_4   Sample_5   Sample_6
#> elevated_set -0.2824106 -0.09203292 -0.2866433 -0.5232409  0.4322547  0.5063751
#> control_set  -0.3933116 -0.28047262 -0.3170458 -0.3956802 -1.0473972 -1.2533668
#>                Sample_7   Sample_8
#> elevated_set  0.4450376  0.3641358
#> control_set  -1.4454834 -0.8702312
```
