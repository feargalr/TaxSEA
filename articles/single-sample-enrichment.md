# Single-sample enrichment with ssTaxSEA

## Why per-person enrichment?

Standard TaxSEA works at the group level. You run a differential
abundance analysis, get a fold change for each taxon, and ask whether
biologically related taxa shift together. That is powerful for
characterising what differs between groups.

But sometimes the more interesting question is about individuals:

- Which patients in my IBD cohort show an oral bacteria signature?
- Does this person’s microbiome look more “disease-like” for a
  particular trait?
- Can I stratify my cohort by enrichment patterns rather than by
  diagnosis alone?

Group-level enrichment cannot answer these questions because it
collapses all individuals into a single summary statistic.
[`ssTaxSEA()`](https://feargalr.github.io/TaxSEA/reference/ssTaxSEA.md)
produces a score for every sample, for every taxon set.

## How it works

For each sample, the counts are centered log-ratio (CLR) transformed,
and the score for a taxon set is simply the **mean of the CLR values of
the set’s members**:

``` math
score_{S,j} = \frac{1}{|S|} \sum_{i \in S} clr_{ij}
```

That is the whole method. There is no ranking, no running-sum statistic
and no cohort-relative standardisation.

Two details are worth stating explicitly.

**The CLR is computed across all taxa you supply, before any subsetting
to set members.** The CLR divides each taxon by the geometric mean of
the sample, so that denominator has to represent the whole community. If
you were to subset the matrix to set members first and transform
afterwards, the denominator would depend on which sets you happened to
be testing — and sets with no real signal would pick up an apparent one.
[`ssTaxSEA()`](https://feargalr.github.io/TaxSEA/reference/ssTaxSEA.md)
does this in the right order; the point matters if you are ever tempted
to hand-roll the calculation.

**The score is a within-sample quantity.** It does not reference the
rest of the cohort at all, which means
[`ssTaxSEA()`](https://feargalr.github.io/TaxSEA/reference/ssTaxSEA.md)
is well defined for a single sample.

## Usage

### From a count matrix

Rows are taxa, columns are samples, and row names are taxon names in the
same format standard TaxSEA expects.

``` r

set.seed(42)
counts <- matrix(rpois(30 * 8, lambda = 10), nrow = 30, ncol = 8)
rownames(counts) <- paste0("Taxon_", seq_len(30))
colnames(counts) <- paste0("Sample_", seq_len(8))

# Plant a signal: raise one set's members in the last four samples
counts[1:6, 5:8] <- counts[1:6, 5:8] * 5

my_sets <- list(
  elevated_set = paste0("Taxon_", 1:6),
  control_set  = paste0("Taxon_", 20:27)
)

scores <- ssTaxSEA(counts, custom_db = my_sets, min_set_size = 3)
dim(scores)
#> [1] 2 8
round(scores, 2)
#>              Sample_1 Sample_2 Sample_3 Sample_4 Sample_5 Sample_6 Sample_7
#> elevated_set     0.02     0.20     0.04    -0.13     1.15     1.37     1.31
#> control_set     -0.09     0.02     0.01    -0.01    -0.29    -0.34    -0.52
#>              Sample_8
#> elevated_set     1.07
#> control_set     -0.14
```

The returned object is a plain numeric matrix with **taxon sets as rows
and samples as columns**. The planted signal is visible as higher scores
in the last four samples of `elevated_set`:

``` r

rowMeans(scores[, 5:8]) - rowMeans(scores[, 1:4])
#> elevated_set  control_set 
#>    1.1948884   -0.3033688
```

### A single sample

Because the score never references the cohort, one sample is enough:

``` r

ssTaxSEA(counts[, 1, drop = FALSE], custom_db = my_sets, min_set_size = 3)
#>                 Sample_1
#> elevated_set  0.01662292
#> control_set  -0.08585010
```

### With the built-in database

Leave `custom_db` unset to use the built-in TaxSEA sets. Row names are
mapped to NCBI taxonomy IDs, exactly as in standard
[`TaxSEA()`](https://feargalr.github.io/TaxSEA/reference/TaxSEA.md).

``` r

scores <- ssTaxSEA(counts)
dim(scores)
```

### From a TreeSummarizedExperiment

[`ssTaxSEA()`](https://feargalr.github.io/TaxSEA/reference/ssTaxSEA.md)
accepts `SummarizedExperiment` and `TreeSummarizedExperiment` objects
directly, extracting the first assay and using
[`rownames()`](https://rdrr.io/r/base/colnames.html) as taxon
identifiers.

``` r

library(SummarizedExperiment)

se <- SummarizedExperiment(assays = list(counts = counts))
scores_se <- ssTaxSEA(se, custom_db = my_sets, min_set_size = 3)

# Identical to the matrix result
all.equal(scores, scores_se)
#> [1] TRUE
```

## Interpreting the scores

**Compare across samples, within a set.** A higher score for
`elevated_set` in sample 6 than in sample 2 means the set’s members make
up more of sample 6’s community. This is the comparison the score is
built for.

**Do not compare across sets.** A set of abundant taxa scores higher
than a set of rare taxa in every sample, regardless of biology, because
its members are simply more abundant. In the example above the two sets
sit at different baselines for exactly this reason:

``` r

rowMeans(scores)
#> elevated_set  control_set 
#>    0.6292661   -0.1682476
```

For heatmaps and clustering, center or scale the rows first so that you
are looking at variation between samples rather than differences in
baseline abundance:

``` r

scaled <- t(scale(t(scores)))
round(scaled, 2)
#>              Sample_1 Sample_2 Sample_3 Sample_4 Sample_5 Sample_6 Sample_7
#> elevated_set    -0.94    -0.66    -0.90    -1.17     0.80     1.14     1.05
#> control_set      0.42     0.96     0.94     0.80    -0.61    -0.89    -1.79
#>              Sample_8
#> elevated_set     0.68
#> control_set      0.17
#> attr(,"scaled:center")
#> elevated_set  control_set 
#>    0.6292661   -0.1682476 
#> attr(,"scaled:scale")
#> elevated_set  control_set 
#>    0.6513343    0.1946639
```

**No p-values are returned.** The score is a descriptive statistic, and
a p-value computed within a single sample would mostly be telling you
whether the set’s members are more abundant than the average taxon — a
property of the set, not of the sample. To test a hypothesis, take one
row of the matrix and compare it across samples with a test suited to
your design:

``` r

group <- rep(c("control", "case"), each = 4)
wilcox.test(scores["elevated_set", ] ~ group)$p.value
#> [1] 0.02857143
```

**Remember that the data are compositional.** Because CLR values within
a sample are constrained, raising one set’s members slightly lowers
everything else. Small shifts in sets you believe to be null are
expected, and with enough samples they can reach statistical
significance while being negligible in size. Report effect sizes
alongside p-values rather than relying on p-values alone.

## Practical considerations

**The pseudocount matters.** Zeros cannot be log-transformed, so
`pseudocount` (default 0.5) is added to every value before the
transform. Adding it to all values rather than only to zeros keeps the
transform monotonic in the counts. The choice does affect the scores, so
report the value you used.

``` r

p1 <- ssTaxSEA(counts, custom_db = my_sets, min_set_size = 3, pseudocount = 1)
round(p1["elevated_set", ] - scores["elevated_set", ], 3)
#> Sample_1 Sample_2 Sample_3 Sample_4 Sample_5 Sample_6 Sample_7 Sample_8 
#>   -0.002   -0.008   -0.002    0.007   -0.032   -0.037   -0.037   -0.027
```

**Supply counts, not proportions.** Adding a pseudocount of 0.5 to
values that sum to 1 would overwhelm the data completely, so
[`ssTaxSEA()`](https://feargalr.github.io/TaxSEA/reference/ssTaxSEA.md)
refuses input whose columns sum to 1. If you only have relative
abundances, rescale them to a count-like scale first.

**BugSigDB is excluded.**
[`ssTaxSEA()`](https://feargalr.github.io/TaxSEA/reference/ssTaxSEA.md)
uses the built-in TaxSEA database (GMRepoV2, MiMeDB, gutMGene, mBodyMap,
BacDive) but does not load BugSigDB. This keeps the function fast and
avoids the `bugsigdbr` dependency.

## When to use ssTaxSEA vs TaxSEA

|  | **TaxSEA** | **ssTaxSEA** |
|----|----|----|
| **Input** | Named vector of fold changes or correlations | Count matrix (taxa x samples) or TreeSummarizedExperiment |
| **Output** | One enrichment result per taxon set | One score per sample per taxon set |
| **Question** | What trait-based shifts characterise the difference between groups? | Which individual samples show a given enrichment pattern? |
| **Use case** | Interpreting DA results, contextualising group differences | Patient stratification, individual-level phenotyping |
| **Significance** | p-values and FDR per set | none; test scores across samples yourself |

In practice you might use both: run
[`TaxSEA()`](https://feargalr.github.io/TaxSEA/reference/TaxSEA.md) to
identify which taxon sets are enriched at the group level, then
[`ssTaxSEA()`](https://feargalr.github.io/TaxSEA/reference/ssTaxSEA.md)
to see how those enrichments distribute across individuals.

## Downstream analysis

The score matrix is designed to slot into standard workflows:

- Heatmaps of per-sample enrichment patterns (scale the rows first)
- Correlation with clinical metadata
- Group comparisons on a single set’s scores
- Clustering or ordination to identify patient subgroups
- Regression or survival analysis using scores as predictors

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.5 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] SummarizedExperiment_1.42.0 Biobase_2.72.0             
#>  [3] GenomicRanges_1.64.0        Seqinfo_1.2.0              
#>  [5] IRanges_2.46.0              S4Vectors_0.50.3           
#>  [7] BiocGenerics_0.58.1         generics_0.1.4             
#>  [9] MatrixGenerics_1.24.0       matrixStats_1.5.0          
#> [11] TaxSEA_1.5.6               
#> 
#> loaded via a namespace (and not attached):
#>  [1] Matrix_1.7-5        jsonlite_2.0.0      compiler_4.6.1     
#>  [4] jquerylib_0.1.4     systemfonts_1.3.2   textshaping_1.0.5  
#>  [7] yaml_2.3.12         fastmap_1.2.0       lattice_0.22-9     
#> [10] R6_2.6.1            XVector_0.52.0      S4Arrays_1.12.0    
#> [13] knitr_1.52          DelayedArray_0.38.2 desc_1.4.3         
#> [16] bslib_0.12.0        rlang_1.3.0         cachem_1.1.0       
#> [19] xfun_0.61           fs_2.1.0            sass_0.4.10        
#> [22] otel_0.2.0          SparseArray_1.12.2  cli_3.6.6          
#> [25] pkgdown_2.2.1       digest_0.6.39       grid_4.6.1         
#> [28] lifecycle_1.0.5     evaluate_1.0.5      ragg_1.5.2         
#> [31] abind_1.4-8         rmarkdown_2.32      tools_4.6.1        
#> [34] htmltools_0.5.9
```
