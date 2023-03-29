# TaxSEA: Taxon Set Enrichment Analysis

TaxSEA is an R package designed to enable rapid annotation of changes observed in a microbiome association study by testing for enrichment for producers of particular metabolites or associations with marker taxa for specific diseases. It focuses specifically on human gut microbiome associations and uses a Kolmogorov-Smirnov test to test if a particular set of taxa is changed relative to a control group.

## Functions

- `get_ncbi_taxon_ids(taxon_names)`: Retrieves NCBI Taxonomy IDs for a list of taxon names.
- `TaxSEA(taxon_ranks, database = "All")`: Performs taxon set enrichment analysis.
- `TaxSEA_barplot(taxsea_results, threshold = 0.2, custom_colors = NULL)`: Creates a bar plot of TaxSEA results.
- `TaxSEA_barcode(taxsea_results, taxon_ranks, taxon_sets, axis_limits = c(-7,7), n_to_plot = 10, boxplot = FALSE)`: Creates a barcode plot of TaxSEA results.

## Usage

```R
# Load required libraries
library(TaxSEA)

# Run TaxSEA
enriched_taxon_sets <- TaxSEA(taxon_ranks)

# Create a bar plot of the results
TaxSEA_barplot(enriched_taxon_sets)


# Create a barcode plot of the results
TaxSEA_barcode(enriched_taxon_sets, taxon_ranks, taxon_sets)
