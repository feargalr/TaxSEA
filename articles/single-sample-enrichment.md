# Single-sample enrichment with ssTaxSEA

## Why per-person enrichment?

Standard TaxSEA works at the group level. You run a differential
abundance analysis, get a fold change for each taxon, and ask whether
biologically related taxa shift together. That’s powerful for
characterising what differs between groups.

But sometimes the more interesting question is about individuals:

- Which patients in my IBD cohort show an oral bacteria signature?
- Does this person’s microbiome look more “disease-like” or
  “health-like” for a particular trait?
- Can I stratify my cohort by enrichment patterns rather than by
  diagnosis alone?

Group-level enrichment can’t answer these questions because it collapses
all individuals into a single summary statistic. `ssTaxSEA` extends
TaxSEA to produce an enrichment score for every sample, for every taxon
set.

## How it works

The core idea is borrowed from single-sample Gene Set Enrichment
Analysis (ssGSEA), adapted for microbiome count data.

The approach has four steps:

1.  **CLR transform.** Raw counts are transformed using the centered
    log-ratio, which accounts for the compositional nature of microbiome
    data. Zeros are replaced with a small pseudocount (0.5) before
    transformation.

2.  **Z-score each taxon across the cohort.** After CLR transformation,
    each taxon is z-scored across all samples. This is the critical
    step. It means that values no longer reflect within-sample
    abundance, but instead reflect how much higher or lower a taxon is
    in a given sample *relative to the cohort mean*. Without this step,
    you would only be asking “are these taxa the most abundant in this
    person?” rather than “does this person have unusually high levels of
    these taxa compared to everyone else?”.

3.  **Rank and score.** For each sample, the z-scores are ranked and an
    ssGSEA-style weighted running-sum enrichment score is computed for
    each taxon set.

4.  **KS test.** A Kolmogorov-Smirnov p-value is computed per sample per
    set, testing whether the set members are shifted relative to all
    taxa in that sample.

The output is a matrix of enrichment scores (samples x taxon sets) and a
corresponding matrix of p-values.

## Why z-scoring matters

This is worth emphasising because it’s the key difference between a
naive per-sample approach and one that actually works.

Consider oral bacteria in an IBD cohort. At the group level, these taxa
are consistently elevated in IBD relative to healthy controls, and
TaxSEA picks this up easily from fold changes. But in any individual
sample, oral bacteria are unlikely to be the *most abundant* taxa. They
might sit at moderate abundance in IBD samples and low abundance in
controls.

If you simply rank each sample’s taxa by their own abundance and run
enrichment, oral bacteria won’t stand out in either group because both
groups have the same dominant gut commensals. The within-sample ranking
buries the between-group signal.

Z-scoring across the cohort solves this. A taxon gets a high z-score in
a sample only if that sample has more of it than average. So in IBD
samples, oral bacteria will have positive z-scores (elevated relative to
the cohort), and in healthy samples, they will have negative z-scores.
The enrichment test then picks up exactly the pattern you expect.

## When to use ssTaxSEA vs TaxSEA

|              | **TaxSEA**                                                          | **ssTaxSEA**                                                                 |
|--------------|---------------------------------------------------------------------|------------------------------------------------------------------------------|
| **Input**    | Named vector of fold changes or correlations                        | Count matrix (taxa x samples) or TreeSummarizedExperiment                    |
| **Output**   | One enrichment result per taxon set                                 | One enrichment score per sample per taxon set                                |
| **Question** | What trait-based shifts characterise the difference between groups? | Which individual samples show a given enrichment pattern?                    |
| **Use case** | Interpreting DA results, contextualising group differences          | Patient stratification, individual-level phenotyping, heterogeneity analysis |

In practice, you might use both: run TaxSEA first to identify which
taxon sets are enriched at the group level, then use ssTaxSEA to see how
those enrichments distribute across individuals.

## Usage

### From a count matrix

The input should be a matrix with taxa as rows and samples as columns.
Row names must be taxon names (species or genus) in the same format as
standard TaxSEA.

``` r
library(TaxSEA)

# counts: taxa (rows) x samples (columns)
res <- ssTaxSEA(counts)

# Enrichment scores: samples x taxon sets
dim(res$scores)
res$scores[1:5, 1:3]

# P-values: samples x taxon sets
res$pvalues[1:5, 1:3]
```

### From a TreeSummarizedExperiment

`ssTaxSEA` accepts TSE objects directly. It will extract the count assay
and use [`rownames()`](https://rdrr.io/r/base/colnames.html) as taxon
identifiers.

``` r
library(mia)
library(TaxSEA)

data(GlobalPatterns, package = "mia")
tse <- GlobalPatterns

# Subset to two sample types
tse <- tse[, colData(tse)$SampleType %in% c("Feces", "Soil")]
tse <- tse[rowSums(assay(tse, "counts") > 0) >= 2, ]

res <- ssTaxSEA(tse, min_set_size = 5)
dim(res$scores)
```

### With custom taxon sets

As with standard TaxSEA, you can provide your own taxon sets. When using
a custom database, names in the sets should match your row names
directly (no NCBI ID conversion is performed).

``` r
my_sets <- list(
  oral_taxa = c("Fusobacterium_nucleatum", "Streptococcus_mutans",
                "Porphyromonas_gingivalis", "Prevotella_intermedia",
                "Treponema_denticola", "Aggregatibacter_actinomycetemcomitans"),
  scfa_producers = c("Faecalibacterium_prausnitzii", "Roseburia_intestinalis",
                     "Eubacterium_rectale", "Coprococcus_catus",
                     "Butyricicoccus_pullicaecorum", "Anaerostipes_caccae")
)

res <- ssTaxSEA(counts, custom_db = my_sets, min_set_size = 3)
```

## Interpreting the output

**Scores.** A positive enrichment score for a given sample and taxon set
means that the taxa in that set tend to have higher abundance in that
sample relative to the rest of the cohort. A negative score means they
tend to be lower. The magnitude reflects how strong and coordinated the
shift is.

**P-values.** The KS test p-value tells you whether the distribution of
z-scores for set members differs from the background in that sample.
Note that these are unadjusted p-values. Since you are testing many sets
across many samples, apply your own multiple testing correction as
appropriate for your analysis.

## Practical considerations

**Cohort composition matters.** Because z-scoring is relative to the
cohort mean, the composition of your cohort affects the scores. If your
cohort is 90% disease samples, the “average” will look disease-like and
scores will be less extreme. A balanced cohort (or at least a
well-defined reference group) generally gives more interpretable
results.

**Sample size.** Z-scoring requires a reasonable number of samples to
estimate the mean and standard deviation for each taxon. With very few
samples (fewer than ~10), estimates will be noisy.

**BugSigDB is excluded.** `ssTaxSEA` uses the built-in TaxSEA database
(GMRepoV2, MiMeDB, gutMGene, mBodyMap, BacDive) but does not load
BugSigDB. This keeps the function fast and avoids the need for the
`bugsigdbr` dependency.

**Downstream analysis.** The score matrix is designed to slot into
standard statistical workflows. You can use it for:

- Heatmaps to visualise per-sample enrichment patterns
- Correlation with clinical metadata
- Group comparisons (e.g. Wilcoxon test on scores between disease and
  control)
- Clustering or ordination to identify patient subgroups
- Survival analysis or regression using scores as predictors
