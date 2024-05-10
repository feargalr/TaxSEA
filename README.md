# TaxSEA: Taxon Set Enrichment Analysis

TaxSEA is an R package designed to enable rapid interpretation of 
differential abundance analysis  or correlation analysis output for 
microbiota data. 
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
data(TaxSEA_test_data)
taxsea_results <- TaxSEA(taxon_ranks=TaxSEA_test_data)

#Enrichments among metabolite producers from gutMgene and MiMeDB
metabolites.df <- taxsea_results$Metabolite_producers

#Enrichments among health and disease signatures from GMRepoV2 and mBodyMap
disease.df <- taxsea_results$Health_associations

#Enrichments amongh published associations from BugSigDB
bsdb.df <- taxsea_results$BugSigdB

```

#### Input 
All that is required for TaxSEA is a vector in R containing ranks (e.g.
log2 fold changes) and names (E.g. species/genus). TaxSEA will not work for 
ranks higher than species or genus. The input should be for all taxa tested, and not 
limited to only a pre-defined set (e.g. do not use a threshold for 
significance or remove any taxa). See example below for format. TaxSEA will
lookup and convert taxon names to NCBI taxonomic identifiers. TaxSEA stores
a commonly observed identifiers internally and so will only look up whatever 
is not covered to save time. 

Input IDs should be in the format of like one of the following
- Species name. E.g. "Bifidobacterium longum", "Bifidobacterium_longum"
- Genus name. E.g. "Bifidobacterium"
- NCBI ID E.g. 216816


```{r input_data}
#Input IDs with the full taxonomic lineage should be split up. E.g.
x <- "d__Bacteria.p__Actinobacteriota.c__Actinomycetes.o__Bifidobacteriales.f__Bifidobacteriaceae.g__Bifidobacterium"
x <- strsplit(x,split="\\.")[[1]][6]
x <- gsub("g__","",x)

## Running this through a vector of IDs may look something like the following
#new_ids <- sapply(as.character(old_ids),function(y) {strsplit(x = y,split="\\.")[[1]][6]})
#new_ids <- gsub("g__","",new_ids)

## Example test data
library(TaxSEA)
data("TaxSEA_test_data")
head(sample(TaxSEA_test_data),4)

```



### Run TaxSEA with test data
```{r example2}
data("TaxSEA_test_data")
taxsea_results <- TaxSEA(taxon_ranks=TaxSEA_test_data)

#Enrichments among metabolite producers from gutMgene and MiMeDB
metabolites.df <- taxsea_results$Metabolite_producers

#Enrichments among health and disease signatures from GMRepoV2 and mBodyMap
disease.df <- taxsea_results$Health_associations

#Enrichments among published associations from BugSigDB
bsdb.df <- taxsea_results$BugSigdB

```

#### Test data
The test data provided with TaxSEA consists of log2 fold changes comparing between healthy 
and IBD. The count data for this was downloaded from curatedMetagenomeData and 
fold changes generated with LinDA.
- Hall et al. A novel Ruminococcus gnavus clade enriched in inflammatory 
bowel disease patients** Genome Med. 2017 Nov 28;9(1):103.
- Pasolli et al. Accessible, curated metagenomic data through 
ExperimentHub. Nat Methods. 2017 Oct 31;14(11):1023-1024. doi: 10.1038/nmeth.4468.
- Zhou et al. LinDA: linear models for differential abundance analysis of microbiome compositional data. Genome Biol. 2022 Apr 14;23(1):95.

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


#### Visualisation of TaxSEA output. 
In other contexts, enrichment scores are visualised by barplots of the
-log10 transformed adjusted P value. TaxSEA contains a function which 
uses ggplot2 to create a barplot results. This splits the enriched and
depleted signatures to plot them going in opposite directions. The Y axis
is ordered by the -log10 FDR also. 

```{r example}
library(TaxSEA)

data("TaxSEA_test_data")
taxsea_results <- TaxSEA(taxon_ranks=TaxSEA_test_data)

#Enrichments among health and disease signatures from GMRepoV2 and mBodyMap
disease.df <- taxsea_results$Health_associations

# Create a bar plot of the results
TaxSEA_barplot(disease.df,threshold=0.05)
```



#### BugSigDB
The format of BugSigDB is that each publication is entered as a "Study", and within this there is different 
experiments and signatures. For example signature 1 may be taxa increased in an experiment, and signature 2
taxa that are decreased. Users can find out more by querying the BugSigDB. See below for an example. 

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

