# TaxSEA

## TaxSEA: Taxon Set Enrichment Analysis

TaxSEA is an R package designed to enable rapid interpretation of
differential abundance analysis or correlation analysis output for
microbiota data. TaxSEA takes as input a vector of genus or species
names and a rank. For example log2 fold changes, or Spearman’s rho.
TaxSEA then uses a Kolmogorov-Smirnov test to identify if a particular
group of species or genera (i.e. a set of taxa such as butyrate
producers) are skewed to one end of the distribution .

Simply put, TaxSEA allows users to rapidly go from a list of species to
which metabolite producers are altered, and if the findings are similar
to a previously published study.

### Installation

To install the latest version of TaxSEA from Bioconductor:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("TaxSEA")
```

### Taxon set database

TaxSEA utilizes taxon sets generated from five reference databases
(**gutMGene**, **GMrepo v2**, **MiMeDB**, **mBodyMap**, **BugSigDB**).

Please cite the appropriate database if using:

- Cheng et al. gutMGene: a comprehensive database for target genes of
  gut microbes and microbial metabolites Nucleic Acids Res. 2022.

- Dai et al. GMrepo v2: a curated human gut microbiome database with
  special focus on disease markers and cross-dataset comparison Nucleic
  Acids Res. 2022.

- Wishart et. al. MiMeDB: the Human Microbial Metabolome Database
  Nucleic Acids Res. 2023.

- Jin et al. mBodyMap: a curated database for microbes across human body
  and their associations with health and diseases. Nucleic Acids Res.
  2022.

- Geistlinger et al. BugSigDB captures patterns of differential
  abundance across a broad range of host-associated microbial
  signatures. Nature Biotech. 2023.

### Test data

The test data provided with TaxSEA consists of log2 fold changes
comparing between healthy and IBD. The count data for this was
downloaded from curatedMetagenomeData and fold changes generated with
LinDA.

- Hall et al. A novel Ruminococcus gnavus clade enriched in inflammatory
  bowel disease patients Genome Med. 2017 Nov 28;9(1):103.
- Pasolli et al. Accessible, curated metagenomic data through
  ExperimentHub. Nat Methods. 2017 Oct 31;14(11):1023-1024. doi:
  10.1038/nmeth.4468.
- Zhou et al. LinDA: linear models for differential abundance analysis
  of microbiome compositional data. Genome Biol. 2022 Apr 14;23(1):95.

### Functions

- `get_taxon_sets(taxon)`: Retrieves taxon sets which contain a
  particular taxon for a list of taxon names.
- `get_ncbi_taxon_ids(taxon_names)`: Retrieves NCBI Taxonomy IDs for a
  list of taxon names.
- `TaxSEA(taxon_ranks, database = "All")`: Taxon set enrichment
  analysis.

### Usage

#### Retrieve sets containing a particular taxon

``` r
library(TaxSEA)

# Retrieve taxon sets containing Bifidobacterium longum.
blong.sets <- get_taxon_sets(taxon="Bifidobacterium_longum")
```

#### Running an enrichment analysis

All that is required for TaxSEA is a named vector of log2 fold changes
between groups for species or genera. TaxSEA will not work for ranks
higher than species or genus. The input should be for all taxa tested,
and not limited to only a pre-defined set (e.g. do not use a threshold
for significance or remove any taxa). See example below for format.
TaxSEA will lookup and convert taxon names to NCBI taxonomic
identifiers. TaxSEA stores commonly observed identifiers internally to
save time.

TaxSEA can also utilise custom databases which should be a named list of
taxon sets. In this case the ID conversion is disabled and it is
expected that the input names and database names will be in the same
format

Input IDs should be in the format of like one of the following - Species
name. E.g. “Bifidobacterium longum”, “Bifidobacterium_longum” - Genus
name. E.g. “Bifidobacterium” - NCBI ID E.g. 216816

``` r
#Input IDs with the full taxonomic lineage should be split up. E.g.
x <- paste0(
  "d__Bacteria.p__Actinobacteriota.c__Actinomycetes.",
  "o__Bifidobacteriales.f__Bifidobacteriaceae.g__Bifidobacterium")
x = strsplit(x,split="\\.")[[1]][6]
x = gsub("g__","",x)
print(x)
```

    ## [1] "Bifidobacterium"

``` r
## Example test data
library(TaxSEA)
data(TaxSEA_test_data)
head(sample(TaxSEA_test_data),4)
```

    ##   Holdemania_filiformis  Clostridium_sp_CAG_242 Bifidobacterium_bifidum 
    ##                  -1.609                   0.275                   3.146 
    ##        Actinomyces_oris 
    ##                   1.336

#### Run TaxSEA with test data

``` r
data("TaxSEA_test_data")
taxsea_results <- TaxSEA(taxon_ranks=TaxSEA_test_data)
```

    ## Using cached version from 2026-01-13 05:56:45

    ## Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
    ## approximate in the presence of ties
    ## Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
    ## approximate in the presence of ties
    ## Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
    ## approximate in the presence of ties
    ## Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
    ## approximate in the presence of ties
    ## Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
    ## approximate in the presence of ties
    ## Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
    ## approximate in the presence of ties
    ## Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
    ## approximate in the presence of ties
    ## Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
    ## approximate in the presence of ties
    ## Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
    ## approximate in the presence of ties
    ## Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
    ## approximate in the presence of ties

``` r
#Enrichments among metabolite producers from gutMgene and MiMeDB
metabolites.df = taxsea_results$Metabolite_producers

#Enrichments among health and disease signatures from GMRepoV2 and mBodyMap
disease.df = taxsea_results$Health_associations

#Enrichments amongh published associations from BugSigDB
bsdb.df = taxsea_results$BugSigdB
```

##### Output

The output is a list of three dataframes providing enrichment results
for metabolite produers, health/disease associations, and published
signatures from BugSigDB. Each dataframe has 5 columns - taxonSetName -
The name of the taxon set tested - median_rank - The median rank of set
members - P value - Kolmogorov-Smirnov test P value. - FDR - P value
adjusted for multiple testing. - TaxonSet - Returns list of taxa in the
set to show what is driving the signal

##### BugSigDB

The format of BugSigDB is that each publication is entered as a “Study”,
and within this there is different experiments and signatures. Should
users wish to find out more information about the signatures, they can
do so by querying that database.

``` r
library(bugsigdbr) #This package is installable via Bioconductor
bsdb <- importBugSigDB() #Import database 
```

    ## Using cached version from 2026-01-13 05:56:45

``` r
#E.g. if the BugSigDB identifier you found enriched was 
#bsdb:74/1/2_obesity:obese_vs_non-obese_DOWN
#This is Study 74, Experiment 1, Signature 2
bsdb[bsdb$Study=="Study 74" & 
     bsdb$Experiment=="Experiment 1" & 
     bsdb$Signature=="Signature 2",]
```

    ##  [1] BSDB ID                    Study                     
    ##  [3] Study design               PMID                      
    ##  [5] DOI                        URL                       
    ##  [7] Authors list               Title                     
    ##  [9] Journal                    Year                      
    ## [11] Keywords                   Experiment                
    ## [13] Location of subjects       Host species              
    ## [15] Body site                  UBERON ID                 
    ## [17] Condition                  EFO ID                    
    ## [19] Group 0 name               Group 1 name              
    ## [21] Group 1 definition         Group 0 sample size       
    ## [23] Group 1 sample size        Antibiotics exclusion     
    ## [25] Sequencing type            16S variable region       
    ## [27] Sequencing platform        Statistical test          
    ## [29] Significance threshold     MHT correction            
    ## [31] LDA Score above            Matched on                
    ## [33] Confounders controlled for Pielou                    
    ## [35] Shannon                    Chao1                     
    ## [37] Simpson                    Inverse Simpson           
    ## [39] Richness                   Signature page name       
    ## [41] Source                     Curated date              
    ## [43] Curator                    Revision editor           
    ## [45] Description                Abundance in Group 1      
    ## [47] MetaPhlAn taxon names      NCBI Taxonomy IDs         
    ## [49] State                      Reviewer                  
    ## <0 rows> (or 0-length row.names)

##### TaxSEA database with other enrichment tools

The TaxSEA function by default uses the Kolmogorov Smirnov test and the
original idea was inspired by gene set enrichment analyses from RNASeq.
Should users wish to use an alternative gene set enrichment analysis
tool the database is formatted in such a way that should be possible.
See below for an example with fast gene set enrichment analysis (fgsea).

``` r
library(fgsea) #This package is installable via Bioconductor
data("TaxSEA_test_data")
data("TaxSEA_db")
#Convert input names to NCBI taxon ids
names(TaxSEA_test_data) = get_ncbi_taxon_ids(names(TaxSEA_test_data))
TaxSEA_test_data = TaxSEA_test_data[!is.na(names(TaxSEA_test_data))]

#Run fgsea
fgsea_results <- fgsea(TaxSEA_db, TaxSEA_test_data, minSize=5, maxSize=500)
```

    ## Warning in preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam, : There are ties in the preranked stats (1.83% of the list).
    ## The order of those tied genes will be arbitrary, which may produce unexpected results.

``` r
sessionInfo()
```

    ## R version 4.5.2 (2025-10-31)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] fgsea_1.36.2     bugsigdbr_1.16.2 TaxSEA_1.3.0     BiocStyle_2.38.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] fastmatch_1.1-6     gtable_0.3.6        xfun_0.55          
    ##  [4] bslib_0.9.0         ggplot2_4.0.1       httr2_1.2.2        
    ##  [7] lattice_0.22-7      vctrs_0.6.5         tools_4.5.2        
    ## [10] generics_0.1.4      curl_7.0.0          parallel_4.5.2     
    ## [13] tibble_3.3.1        RSQLite_2.4.5       blob_1.2.4         
    ## [16] pkgconfig_2.0.3     Matrix_1.7-4        data.table_1.18.0  
    ## [19] dbplyr_2.5.1        RColorBrewer_1.1-3  S7_0.2.1           
    ## [22] desc_1.4.3          lifecycle_1.0.5     compiler_4.5.2     
    ## [25] farver_2.1.2        textshaping_1.0.4   codetools_0.2-20   
    ## [28] htmltools_0.5.9     sass_0.4.10         yaml_2.3.12        
    ## [31] pillar_1.11.1       pkgdown_2.2.0       jquerylib_0.1.4    
    ## [34] BiocParallel_1.44.0 cachem_1.1.0        tidyselect_1.2.1   
    ## [37] digest_0.6.39       dplyr_1.1.4         purrr_1.2.1        
    ## [40] bookdown_0.46       cowplot_1.2.0       fastmap_1.2.0      
    ## [43] grid_4.5.2          cli_3.6.5           magrittr_2.0.4     
    ## [46] withr_3.0.2         filelock_1.0.3      scales_1.4.0       
    ## [49] rappdirs_0.3.3      bit64_4.6.0-1       rmarkdown_2.30     
    ## [52] bit_4.6.0           ragg_1.5.0          memoise_2.0.1      
    ## [55] evaluate_1.0.5      knitr_1.51          BiocFileCache_3.0.0
    ## [58] rlang_1.1.7         Rcpp_1.1.1          glue_1.8.0         
    ## [61] DBI_1.2.3           BiocManager_1.30.27 jsonlite_2.0.0     
    ## [64] R6_2.6.1            systemfonts_1.3.1   fs_1.6.6
