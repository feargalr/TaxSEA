# TaxSEA: Taxon Set Enrichment Analysis

[Check out the preprint by clicking here](https://www.biorxiv.org/content/10.1101/2024.11.20.624438v1) 

TaxSEA is designed for microbiome researchers seeking deeper functional
insights into their datasets.A common question in human microbiome 
research is whether a detected pattern of microbial changes has been 
previously reported in other diseases or contexts. However, current 
methods to answer this question are highly limited, often requiring 
researchers to manually search literature, cross-reference multiple 
databases, or repurpose tools not originally designed for microbiome 
data. TaxSEA directly addresses this gap by integrating multiple public
microbiome databases, allowing researchers to systematically test for 
enrichment of known disease signatures, metabolite producers, or other
biologically relevant taxon sets. TaxSEA is available as R package.

TaxSEA takes as input a vector of genus or species names and a rank. 
For example log2 fold changes or Spearman's rho. TaxSEA then uses
a Kolmogorov-Smirnov test to identify if a particular group of species
or genera (i.e. a set of taxa such as butyrate producers) are skewed to 
one end of the distribution . 

Simply put, TaxSEA allows users to rapidly convert species/genus level 
changes to alterations in 
- Metabolite producers
- Disease signatures
- Previously published associations

Note: Although TaxSEA in principle can be applied to microbiome data from
any source, the databases utilized largely cover human associated microbiomes
and the human gut microbiome in particular. As such TaxSEA will likely perform
best on human gut microbiome data. 


## Installation
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("TaxSEA")
```

## Taxon set database
TaxSEA utilizes taxon sets generated from five reference databases 
(**gutMGene**, **GMrepo v2**, **MiMeDB**, **mBodyMap**, **BugSigDB**). 

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


## Usage
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

