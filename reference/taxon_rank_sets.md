# Taxon Rank Set Enrichment Analysis

Groups species by taxonomic ranks and performs TaxSEA enrichment
analysis at each rank level. Returns a named list of data frames, one
per taxonomic rank (excluding species).

## Usage

``` r
taxon_rank_sets(taxon_ranks, lineage_df, min_set_size = 5, max_set_size = 100)
```

## Arguments

- taxon_ranks:

  A named numeric vector of log2 fold changes. Names should be feature
  identifiers matching the `species` column in `lineage_df`, or matching
  [`rownames()`](https://rdrr.io/r/base/colnames.html) of a
  SummarizedExperiment/TreeSummarizedExperiment.

- lineage_df:

  Either a data frame or a
  `SummarizedExperiment`/`TreeSummarizedExperiment` object.

  **Data frame input:** Must include a `species` column and one or more
  taxonomic rank columns (e.g., kingdom, phylum, class, order, family,
  genus). The `species` column is used for matching and is excluded from
  the enrichment analysis.

  **SummarizedExperiment input:** Taxonomy is extracted from
  [`rowData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  and feature identifiers from
  [`rownames()`](https://rdrr.io/r/base/colnames.html). If the `mia`
  package is installed,
  [`taxonomyRanks()`](https://microbiome.github.io/mia/reference/taxonomy-methods.html)
  is used to identify taxonomic rank columns; otherwise all
  [`rowData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  columns are used. The Species rank column is automatically excluded
  from the analysis. Requires the `SummarizedExperiment` package.

- min_set_size:

  Minimum number of species in a set to include in the analysis. Default
  is 5.

- max_set_size:

  Maximum number of species in a set to include in the analysis. Default
  is 100.

## Value

A named list of data frames, one per taxonomic rank. Each data frame
contains columns: taxonSetName, median_rank_of_set_members, PValue,
Test_statistic, and FDR. Ranks that produce no valid sets (e.g., due to
size filtering) are included as empty data frames with a message.

## Examples

``` r
# --- Example 1: Data frame input ---
# Create a lineage data frame (e.g., parsed from curatedMetagenomicData)
# The 'species' column must match the names in taxon_ranks.

lineage_df <- data.frame(
  species = c("Cutibacterium_acnes", "Klebsiella_pneumoniae",
              "Propionibacterium_humerusii", "Moraxella_osloensis",
              "Enhydrobacter_aerosaccus", "Staphylococcus_capitis",
              "Staphylococcus_epidermidis", "Staphylococcus_aureus",
              "Escherichia_coli", "Enterobacter_cloacae",
              "Pseudomonas_aeruginosa", "Acinetobacter_baumannii",
              "Lactobacillus_rhamnosus", "Lactobacillus_acidophilus",
              "Bifidobacterium_longum", "Bifidobacterium_breve"),
  kingdom = rep("Bacteria", 16),
  phylum = c("Actinobacteria", "Proteobacteria",
             "Actinobacteria", "Proteobacteria",
             "Proteobacteria", "Firmicutes",
             "Firmicutes", "Firmicutes",
             "Proteobacteria", "Proteobacteria",
             "Proteobacteria", "Proteobacteria",
             "Firmicutes", "Firmicutes",
             "Actinobacteria", "Actinobacteria"),
  class = c("Actinobacteria", "Gammaproteobacteria",
            "Actinobacteria", "Gammaproteobacteria",
            "Alphaproteobacteria", "Bacilli",
            "Bacilli", "Bacilli",
            "Gammaproteobacteria", "Gammaproteobacteria",
            "Gammaproteobacteria", "Gammaproteobacteria",
            "Bacilli", "Bacilli",
            "Actinobacteria", "Actinobacteria"),
  order = c("Propionibacteriales", "Enterobacterales",
            "Propionibacteriales", "Pseudomonadales",
            "Rhodospirillales", "Bacillales",
            "Bacillales", "Bacillales",
            "Enterobacterales", "Enterobacterales",
            "Pseudomonadales", "Pseudomonadales",
            "Lactobacillales", "Lactobacillales",
            "Bifidobacteriales", "Bifidobacteriales"),
  family = c("Propionibacteriaceae", "Enterobacteriaceae",
             "Propionibacteriaceae", "Moraxellaceae",
             "Rhodospirillaceae", "Staphylococcaceae",
             "Staphylococcaceae", "Staphylococcaceae",
             "Enterobacteriaceae", "Enterobacteriaceae",
             "Pseudomonadaceae", "Moraxellaceae",
             "Lactobacillaceae", "Lactobacillaceae",
             "Bifidobacteriaceae", "Bifidobacteriaceae"),
  genus = c("Cutibacterium", "Klebsiella",
            "Cutibacterium", "Moraxella",
            "Enhydrobacter", "Staphylococcus",
            "Staphylococcus", "Staphylococcus",
            "Escherichia", "Enterobacter",
            "Pseudomonas", "Acinetobacter",
            "Lactobacillus", "Lactobacillus",
            "Bifidobacterium", "Bifidobacterium"),
  stringsAsFactors = FALSE
)

set.seed(42)
fc <- setNames(rnorm(16), lineage_df$species)
results <- taxon_rank_sets(fc, lineage_df, min_set_size = 2)
names(results)
#> [1] "kingdom" "phylum"  "class"   "order"   "family"  "genus"  
results$family
#>           taxonSetName median_rank_of_set_members    PValue Test_statistic FDR
#> 1     Lactobacillaceae                -0.83382473 0.1666667      0.7857143   1
#> 2        Moraxellaceae                 1.45975400 0.4666667      0.5714286   1
#> 3 Propionibacteriaceae                 0.86704343 0.7000000      0.5000000   1
#> 4    Staphylococcaceae                -0.09465904 0.9470588      0.2857143   1
#> 5   Bifidobacteriaceae                 0.25131453 1.0000000      0.2857143   1
#> 6   Enterobacteriaceae                -0.06271410 1.0000000      0.1904762   1

# --- Example 2: SummarizedExperiment / TreeSummarizedExperiment input ---
if (FALSE) { # \dontrun{
library(mia)
data(GlobalPatterns, package = "mia")
tse <- GlobalPatterns

# Run differential abundance (e.g., ALDEx2) to get fold changes
# aldex_out <- ... (your DA analysis)
# fc <- aldex_out$effect
# names(fc) <- rownames(aldex_out)

# Run taxon rank set enrichment directly from the TSE
results <- taxon_rank_sets(fc, tse, min_set_size = 5)
names(results)  # Kingdom, Phylum, Class, Order, Family, Genus
results$Family
} # }
```
