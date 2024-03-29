# TaxSEA: Taxon Set Enrichment Analysis

TaxSEA is an R package designed to enable rapid annotation of changes observed in a microbiome association study by testing for enrichment for groups of species associated with either a particular phenotype, or known to producer a particular metabolite.It focuses specifically on human gut microbiome associations and uses a Kolmogorov-Smirnov test to test if a particular set of taxa is changed relative to a control group. A key goal of the TaxSEA project is to move researchers away from describing microbiome alterations simply as 'dysbiotic' and towards more specific descriptions that allow easier interpretation (e.g. depleted for taxa known to produce GABA and enriched for taxa asscociated with diabetes in other studies)

## Taxon set database
TaxSEA utilizes taxon sets generated from five reference databases (**gutMGene**, **GMrepo v2**, **MiMeDB**, **mBodyMap**, **BugSigDB**). 

Please cite the appropriate database if using:
- Cheng et al. gutMGene: a comprehensive database for target genes of gut microbes and microbial metabolites Nucleic Acids Res. 2022.
- Dai et al. GMrepo v2: a curated human gut microbiome database with special focus on disease markers and cross-dataset comparison Nucleic Acids Res. 2022.
- Wishart et. al. MiMeDB: the Human Microbial Metabolome Database Nucleic Acids Res. 2023.
- Jin et al. mBodyMap: a curated database for microbes across human body and their associations with health and diseases. Nucleic Acids Res. 2022.
- Geistlinger et al. BugSigDB captures patterns of differential abundance across a broad range of host-associated microbial signatures. Nature Biotech. 2023. 

## Interpretation of results

There is an enourmous amount we still do not know about the human microbiome. TaxSEA, or any similar analysis, will ultimately be limited by the quality of the databases available, which in turn are limited by our current knowledge. Thus, our reccomendation is the treat results as a starting point for further interpretation and hypothesis generation. 


## Test data
The test data provided with TaxSEA consists of log2 fold changes comparing between healthy and IBD. The count data for this was downloaded from curatedMetagenomeData and fold changes generated with LinDA.
- Hall et al. A novel Ruminococcus gnavus clade enriched in inflammatory bowel disease patients** Genome Med. 2017 Nov 28;9(1):103.
- Pasolli et al. Accessible, curated metagenomic data through ExperimentHub. Nat Methods. 2017 Oct 31;14(11):1023-1024. doi: 10.1038/nmeth.4468.
- Zhou et al. LinDA: linear models for differential abundance analysis of microbiome compositional data. Genome Biol. 2022 Apr 14;23(1):95.

## Installation
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
taxsea_results <- TaxSEA(taxon_ranks=TaxSEA_test_data)

#Enrichments among metabolite producers from gutMgene and MiMeDB
metabolites.df = taxsea_results$Metabolite_producers

#Enrichments among health and disease signatures from GMRepoV2 and mBodyMap
disease.df = taxsea_results$Health_associations

#Enrichments amongh published associations from BugSigDB
bsdb.df = taxsea_results$BugSigdB

```

#### Input 
All that is required for TaxSEA is a named vector of log2 fold changes between groups. This should be for all taxa tested, and not limited to only a pre-defined set (e.g. do not use a threshold for significance or remove any taxa). See example below for format. TaxSEA will lookup and convert taxon names to NCBI taxonomic identifiers. TaxSEA stores a commonly used identifiers internally and so will only look up whatever is not covered to save time. 
```{r output}
> head(sample(TaxSEA_test_data),4)
Bacteroides_thetaiotaomicron           Blautia_sp_CAG_257          Ruminococcus bromii       Clostridium_disporicum 
                       1.908                        3.650                       -5.038                        0.300 
```

#### Output
The output is a list of three dataframes providing enrichment results for metabolite produers, health/disease associations, and published signatures from BugSigDB.
Each dataframe has 5 columns
- taxonSetName - The name of the taxon set tested
- NES - Normalized enrichment score. This is simply the median fold change across the entire set
- P value - Kolmogorov-Smirnov test P value.
- FDR - P value adjusted for multiple testing. 
- TaxonSet - Returns list of taxa in the set to show what is driving the signal


#### TaxSEA database with other enrichment tools
The TaxSEA function by default uses the Kolmogorov Smirnov test. This was the original test used for gene set enrichment analysis, however subsequent methods have been developed based on alternative approaches. One particularly popular package is fast gene set enrichment analysis (fgsea) which is based upon a permutation type approach. The TaxSEA database is compatible with fgsea as long as the input identifiers are NCBI IDs. If they are not currently in that format we provide a function for converting. See example code below. 

```{r output}
library(fgsea) #This package is installable via Bioconductor

#Convert input names to NCBI taxon ids
names(TaxSEA_test_data) = get_ncbi_taxon_ids(names(TaxSEA_test_data))
TaxSEA_test_data = TaxSEA_test_data[!is.na(names(TaxSEA_test_data))]

#Run fgsea
fgsea_results <- fgsea(TaxSEA_db, TaxSEA_test_data, minSize=5, maxSize=500)
```

#### Visualisation 
TaxSEA contains two functions which uses ggplot2 to plot results. 
```{r example}
# Create a bar plot of the results
TaxSEA_barplot(enriched_taxon_sets)
```
![TaxSEA Barplot 1](https://user-images.githubusercontent.com/7561275/228441264-f233b7ac-6030-4208-a48a-a43a92163b33.png)

Barcode plots are routinely used for plotting enrichment signatures. Here all taxa are ordered by the their fold change along the X axis, and a single black line indicates where the taxa in that set are positioned along that order. In the example below you can see all butyrate producers from the GutMGene taxa set are decreased. 
```{r example}
# Create a barcode plot of the results
TaxSEA_barcode(enriched_taxon_sets[grepl("GutMGene",enriched_taxon_sets$taxonSetName),], taxon_ranks=TaxSEA_test_data)

```
![Barcode plot2](https://user-images.githubusercontent.com/7561275/228441385-cf9ffd3c-ab51-4e2b-beda-7398148447df.png)


