# TaxSEA: Taxon Set Enrichment Analysis

TaxSEA is an R package designed to enable rapid interpretation of 
differential analysis output or correlation analysis output. 
TaxSEA takes as input a vector of genus or species names and a rank. 
For example log2 fold changes, or Spearman's rho. TaxSEA then uses
a Kolmogorov-Smirnov test to identify if a particular group of species
or genera (i.e. a taxon set such as producers of butyrate) skewed to 
one end of the distribution of fold changes or correlation coefficients. 

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

## Installation
```{r example}
library(devtools)
install_github("feargalr/TaxSEA")
library(TaxSEA)
```


## Usage
#### Quick start
```{r example}
library(TaxSEA)

# Retrieve taxon sets containing Bifidobacterium longum.
blong.sets <- get_taxon_sets(taxon="Bifidobacterium_longum")

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
All that is required for TaxSEA is a named vector of log2 fold changes 
between groups. This should be for all taxa tested, and not limited to 
only a pre-defined set (e.g. do not use a threshold for significance or 
remove any taxa). See example below for format. TaxSEA will lookup and 
convert taxon names to NCBI taxonomic identifiers. TaxSEA stores a 
commonly observed identifiers internally and so will only look up whatever 
is not covered to save time. 

#### Test data
The test data provided with TaxSEA consists of log2 fold changes comparing between healthy 
and IBD. The count data for this was downloaded from curatedMetagenomeData and 
fold changes generated with LinDA.
- Hall et al. A novel Ruminococcus gnavus clade enriched in inflammatory 
bowel disease patients** Genome Med. 2017 Nov 28;9(1):103.
- Pasolli et al. Accessible, curated metagenomic data through 
ExperimentHub. Nat Methods. 2017 Oct 31;14(11):1023-1024. doi: 10.1038/nmeth.4468.
- Zhou et al. LinDA: linear models for differential abundance analysis of microbiome compositional data.
- Genome Biol. 2022 Apr 14;23(1):95.

```{r output}
> head(sample(TaxSEA_test_data),4)
Bacteroides_thetaiotaomicron           Blautia_sp_CAG_257          Ruminococcus bromii       Clostridium_disporicum 
                       1.908                        3.650                       -5.038                        0.300 
```

#### Output
The output is a list of three data frames providing enrichment results for metabolite produers, 
health/disease associations, and published signatures from BugSigDB.
Each dataframe has 5 columns
- taxonSetName - The name of the taxon set tested
- NES - Normalized enrichment score. This is simply the median rank across 
all detected members in the set. 
- P value - Kolmogorov-Smirnov test P value.
- FDR - P value adjusted for multiple testing. 
- TaxonSet - Returns list of taxa in the set to show what is driving the signal


#### BugSigDB
The format of BugSigDB is that each publication is entered as a "Study", and within this there is different 
experiments and signatures. 
Should users wish to find out more information about the signatures. For 
example different signatures may be the taxa that are increased or 
decreased. 
from the BugSigDB, they can do so by querying that database.  

```{r output}
library(bugsigdbr) #This package is installable via Bioconductor
bsdb <- importBugSigDB() #Import database 

#E.g. if the BugSigDB identifier you found enriched was #bsdb:74/1/2_obesity:obese_vs_non-obese_DOWN
#This is Study 74, Experiment 1, Signature 2
bsdb[bsdb$Study=="Study 74" & 
     bsdb$Experiment=="Experiment 1" & 
     bsdb$Signature=="Signature 2",]



```


#### TaxSEA database with other enrichment tools
The TaxSEA function by default uses the Kolmogorov Smirnov test and the 
original idea was inspired by gene set enrichment analyses from RNASeq.
Should users wish to use an alternative gene set enrichment analysis tool
the database is formatted in such a way that should be possible. See below
for an example with fast gene set enrichment analysis (fgsea). 

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
