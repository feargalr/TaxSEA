# TaxSEA: Taxon Set Enrichment Analysis

[![Bioconductor](https://www.bioconductor.org/shields/years-in-bioc/TaxSEA.svg)](https://bioconductor.org/packages/TaxSEA)  
📦 Available on [**Bioconductor**](https://bioconductor.org/packages/TaxSEA)  
📝 [Read the preprint](https://www.biorxiv.org/content/10.1101/2024.11.20.624438v1)

TaxSEA helps microbiome researchers test for enrichment in known microbial signatures, including:

- Metabolite producers  
- Disease associations  
- Previously published microbiome signatures

TaxSEA takes as input a vector of genus or species names and a rank. 
For example log2 fold changes or Spearman's rho. TaxSEA then uses
a Kolmogorov-Smirnov test to identify if a particular group of species
or genera (i.e. a set of taxa such as butyrate producers) are skewed to 
one end of the distribution . 

Note: Although TaxSEA in principle can be applied to microbiome data from
any source, the databases utilized largely cover human associated microbiomes
and the human gut microbiome in particular. As such TaxSEA will likely perform
best on human gut microbiome data. 

### Taxon set database
By default TaxSEA utilizes taxon sets generated from five reference databases 
(**gutMGene**, **GMrepo v2**, **MiMeDB**, **mBodyMap**, **BugSigDB**). See below for 
examples of using custom databases or taxonomically defined taxon sets. 

Please cite the appropriate database if using:

- Cheng et al. gutMGene: a comprehensive database for target genes of gut microbes and
microbial metabolites Nucleic Acids Res. 2022.
- Dai et al. GMrepo v2: a curated human gut microbiome database with special focus on
disease markers and cross-dataset comparison Nucleic Acids Res. 2022.
- Wishart et. al. MiMeDB: the Human Microbial Metabolome Database Nucleic Acids Res. 2023.
- Jin et al. mBodyMap: a curated database for microbes across human body and their
associations with health and diseases. Nucleic Acids Res. 2022.
- Geistlinger et al. BugSigDB captures patterns of differential abundance across a broad
range of host-associated microbial signatures. Nature Biotech. 2023. 

### Installation
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("TaxSEA")

```


### Usage
#### Quick start
```r
library(TaxSEA)

# Retrieve taxon sets containing Bifidobacterium longum.
blong.sets <- get_taxon_sets(taxon="Bifidobacterium_longum")

# Run TaxSEA with test data provided
data(TaxSEA_test_data)
taxsea_results <- TaxSEA(taxon_ranks=TaxSEA_test_data)

#Enrichments among metabolite producers from gutMgene and MiMeDB
metabolites.df <- taxsea_results$Metabolite_producers

#Enrichments among health and disease signatures from GMRepoV2 and mBodyMap
disease.df <- taxsea_results$Health_associations

#Enrichments amongh published associations from BugSigDB
bsdb.df <- taxsea_results$BugSigdB

```
