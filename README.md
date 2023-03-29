# TaxSEA: Taxon Set Enrichment Analysis

TaxSEA is an R package designed to enable rapid annotation of changes observed in a microbiome association study by testing for enrichment for producers of particular metabolites or associations with marker taxa for specific diseases. It focuses specifically on human gut microbiome associations and uses a Kolmogorov-Smirnov test to test if a particular set of taxa is changed relative to a control group.

**_Note - TaxSEA has not yet undergone extensive validation. Please treat results with appropriate caution_**


TaxSEA utilizes taxon sets generated from three reference databases (**gutMGene**, **GMrepo v2**, **MiMeDB**). All of gutMgene is utilized, wheras only marker taxa from GMRepo v2 and a manually curated subset from MiMeDB are included. Please cite the appropriate database if using:
- Cheng et al. gutMGene: a comprehensive database for target genes of gut microbes and microbial metabolites_ Nucleic Acids Res. 2022 Jan 7;50(D1):D795-D800.
- Dai et al. GMrepo v2: a curated human gut microbiome database with special focus on disease markers and cross-dataset comparison_ Nucleic Acids Res. 2022 Jan 7;50(D1):D777-D784.
- Wishart et. al. MiMeDB: the Human Microbial Metabolome Database_ Nucleic Acids Res. 2023 Jan 6;51(D1):D611-D620. 

The test data provided with TaxSEA consists of log2 fold changes comparing between healthy and inflammatory bowel disease. The count data for this was downloaded from curatedMetagenomeSeq and fold changes generated with LinDA.
- Hall et al. A novel Ruminococcus gnavus clade enriched in inflammatory bowel disease patients** Genome Med. 2017 Nov 28;9(1):103.
- Pasolli et al. Accessible, curated metagenomic data through ExperimentHub. Nat Methods. 2017 Oct 31;14(11):1023-1024. doi: 10.1038/nmeth.4468.
- Zhou et al. LinDA: linear models for differential abundance analysis of microbiome compositional data. Genome Biol. 2022 Apr 14;23(1):95.
# Installation


```{r example}
library(devtools)
install_github("feargalr/TaxSEA")
library(TaxSEA)
```


## Functions

- `get_ncbi_taxon_ids(taxon_names)`: Retrieves NCBI Taxonomy IDs for a list of taxon names.
- `TaxSEA(taxon_ranks, database = "All")`: Performs taxon set enrichment analysis.
- `TaxSEA_barplot(taxsea_results, threshold = 0.2, custom_colors = NULL)`: Creates a bar plot of TaxSEA results.
- `TaxSEA_barcode(taxsea_results, taxon_ranks, taxon_sets, axis_limits = c(-7,7), n_to_plot = 10, boxplot = FALSE)`: Creates a barcode plot of TaxSEA results.

## Usage
#### Quick start
```{r example}
# Load required libraries
library(TaxSEA)

# Run TaxSEA with test data provided
head(TaxSEA_test_data)
enriched_taxon_sets <- TaxSEA(taxon_ranks=TaxSEA_test_data)
```

#### Input 
All that is required for TaxSEA is a named vector of log2 fold changes between groups. See example below. TaxSEA will lookup and convert taxon names to NCBI taxonomic identifiers. TaxSEA stores a commonly used identifiers internally and so will only look up whatever is not covered to save time. 
```{r output}
> head(sample(TaxSEA_test_data),4)
Bacteroides_thetaiotaomicron           Blautia_sp_CAG_257          Ruminococcus bromii       Clostridium_disporicum 
                       1.908                        3.650                       -5.038                        0.300 
```

#### Output
The output is a dataframe with 5 columns
- taxonSetName - The name of the taxon set tested
- NES - Normalized enrichment score. This is simply the median fold change across the entire set
- P value - Kolmogorov-Smirnov test P value.
- FDR - P value adjusted for multiple testing. 
- TaxonSet - Returns list of taxa in the set to show what is driving the signal

#### Visualisation 
TaxSEA contains two functions which uses ggplot2 to plot results. The test data is a comparison between healthy and IBD. Higher fold change values (labelled as upregulated) indicate sets of taxa increased in IBD relative to control. 

```{r example}
# Create a bar plot of the results
TaxSEA_barplot(enriched_taxon_sets)
```
![TaxSEA Barplot 1](https://user-images.githubusercontent.com/7561275/228441264-f233b7ac-6030-4208-a48a-a43a92163b33.png)

Barcode plots are routinely used for plotting enrichment signatures. Here all taxa are ordered by the their fold change along the X axis, and a single black line indicates where the taxa in that set are positioned along that order. 
```{r example}
# Create a barcode plot of the results
TaxSEA_barcode(enriched_taxon_sets[grepl("GutMGene",enriched_taxon_sets$taxonSetName),], taxon_ranks=TaxSEA_test_data)

```
![Barcode plot2](https://user-images.githubusercontent.com/7561275/228441385-cf9ffd3c-ab51-4e2b-beda-7398148447df.png)


