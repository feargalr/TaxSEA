# BacDive physiology sets

TaxSEA includes taxon sets built from
[BacDive](https://bacdive.dsmz.de), the Bacterial Diversity
Metadatabase. Each set groups species that share a measured
physiological property: they tolerate oxygen, stain Gram-negative,
ferment glucose, reduce nitrate, grow above 45 °C, and so on. Testing
these sets asks whether a shift in your data is a shift in *what the
bacteria can do*, not just *which bacteria are there*.

Since TaxSEA 1.5.5 these sets are rebuilt from a complete harvest of
BacDive and contain measured laboratory data only.

## What is in them

``` r

library(TaxSEA)
data(TaxSEA_db)

bacdive <- TaxSEA_db[startsWith(names(TaxSEA_db), "BacDive_")]
length(bacdive)                                   # number of sets
#> [1] 1776
length(unique(unlist(bacdive, use.names = FALSE))) # taxa covered
#> [1] 9045
```

Set names read `BacDive_<Family>_<value>`, optionally followed by an
assay context and `_negative`:

``` r

family <- sub("^BacDive_([^_]+)_.*", "\\1", names(bacdive))
sort(table(family), decreasing = TRUE)
#> family
#>          Uses        Enzyme          Prod        Murein      Produces 
#>          1348           240            53            50            13 
#>         Shape          NaCl          Temp          Test     Nutrition 
#>            11            10             9             8             6 
#>   Observation        Oxygen            pH            GC          Gram 
#>             5             4             4             3             3 
#>     Tolerance      Motility PathogenHuman         Spore 
#>             3             2             2             2
```

| Family | What it describes | Examples |
|----|----|----|
| `Oxygen` | oxygen tolerance | `BacDive_Oxygen_anaerobe`, `BacDive_Oxygen_facultative_anaerobe` |
| `Gram`, `Motility`, `Spore`, `Shape` | cell properties | `BacDive_Gram_negative`, `BacDive_Spore_yes`, `BacDive_Shape_coccoid` |
| `Temp`, `NaCl`, `pH`, `GC` | growth range and genome GC, from numeric measurements | `BacDive_Temp_thermophile`, `BacDive_NaCl_halotolerant` |
| `Uses` | substrate use, split by assay | `BacDive_Uses_glucose_fermentation`, `BacDive_Uses_nitrate_reduction` |
| `Enzyme` | enzyme activities, including hydrolysis | `BacDive_Enzyme_catalase`, `BacDive_Enzyme_esculin_hydrolysis` |
| `Prod`, `Produces`, `Test` | metabolites produced; classic biochemical tests | `BacDive_Prod_indole`, `BacDive_Produces_vitamin_B12` |
| `Observation`, `Murein`, `Nutrition`, `Tolerance`, `PathogenHuman` | quinones, peptidoglycan type, nutrition type and more | `BacDive_Observation_menaquinone` |

### Assay context matters

Growing on a substrate, fermenting it, producing acid or gas from it,
reducing it and respiring it are different physiology, so they are
different sets:

``` r

grep("^BacDive_Uses_nitrate", names(bacdive), value = TRUE)
#>  [1] "BacDive_Uses_nitrate"                     
#>  [2] "BacDive_Uses_nitrate_gas"                 
#>  [3] "BacDive_Uses_nitrate_gas_negative"        
#>  [4] "BacDive_Uses_nitrate_negative"            
#>  [5] "BacDive_Uses_nitrate_Nsource"             
#>  [6] "BacDive_Uses_nitrate_Nsource_negative"    
#>  [7] "BacDive_Uses_nitrate_reduction"           
#>  [8] "BacDive_Uses_nitrate_reduction_negative"  
#>  [9] "BacDive_Uses_nitrate_respiration"         
#> [10] "BacDive_Uses_nitrate_respiration_negative"
```

### Negative sets

A species that was tested and found negative is not the same as one that
was never tested. Sets ending in `_negative` hold the tested negatives,
so you can ask, for example, whether catalase-negative species shift:

``` r

lengths(bacdive[c("BacDive_Enzyme_catalase", "BacDive_Enzyme_catalase_negative")])
#>          BacDive_Enzyme_catalase BacDive_Enzyme_catalase_negative 
#>                             3828                              857
```

## Using them

BacDive results come back in their own table:

``` r

data(TaxSEA_test_data)
res <- TaxSEA(taxon_ranks = TaxSEA_test_data)
#> Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
#> approximate in the presence of ties
#> Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
#> approximate in the presence of ties
#> Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
#> approximate in the presence of ties
#> Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
#> approximate in the presence of ties
#> Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
#> approximate in the presence of ties
#> Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
#> approximate in the presence of ties
#> Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
#> approximate in the presence of ties
#> Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
#> approximate in the presence of ties
#> Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
#> approximate in the presence of ties
#> Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
#> approximate in the presence of ties
#> Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
#> approximate in the presence of ties
#> Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
#> approximate in the presence of ties
#> Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
#> approximate in the presence of ties
#> Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
#> approximate in the presence of ties
#> Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
#> approximate in the presence of ties
#> Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
#> approximate in the presence of ties
#> Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
#> approximate in the presence of ties
#> Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
#> approximate in the presence of ties
#> Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
#> approximate in the presence of ties
#> Warning in ks.test.default(taxon_set_ranks, taxon_ranks): p-value will be
#> approximate in the presence of ties
head(res$BacDive_bacterial_physiology[, c("taxonSetName",
                                         "median_rank_of_set_members",
                                         "FDR")])
#>                                        taxonSetName median_rank_of_set_members
#> 1                    BacDive_Uses_nitrate_reduction                     1.7750
#> 2         BacDive_Uses_ribose_fermentation_negative                     1.7505
#> 3 BacDive_Enzyme_gamma-glutamyltransferase_negative                     1.3675
#> 4       BacDive_Uses_D-xylose_fermentation_negative                     1.4100
#> 5      BacDive_Enzyme_lysine_decarboxylase_negative                     1.3990
#> 6               BacDive_Uses_pullulan_acid_negative                     1.7750
#>          FDR
#> 1 0.01220060
#> 2 0.02854591
#> 3 0.02854591
#> 4 0.02854591
#> 5 0.02854591
#> 6 0.02856967
```

To see which BacDive sets a species belongs to:

``` r

ecoli <- get_taxon_sets("Escherichia_coli")
head(grep("^BacDive_Oxygen|^BacDive_Gram|^BacDive_Uses_nitrate", ecoli,
          value = TRUE))
#> [1] "BacDive_Gram_negative"               "BacDive_Oxygen_facultative_anaerobe"
#> [3] "BacDive_Uses_nitrate_reduction"
```

Species names can be given with spaces or underscores, and former names
of reclassified species work too (*Bacteroides vulgatus* and
*Phocaeicola vulgatus* both resolve). NCBI taxonomy IDs can be given
directly.

## How many strains stand behind a set

Every set ships with its provenance: which BacDive field it came from
and how many species and strains support it.

``` r

prov <- read.delim(system.file("extdata", "BacDive_set_provenance.tsv",
                               package = "TaxSEA"))
head(prov[order(-prov$n_taxa),
          c("set_name", "family", "n_species", "n_strains_total")])
#>                             set_name            family n_species
#> 1             BacDive_Temp_mesophile      culture_temp      8061
#> 2                  BacDive_Shape_rod             shape      4277
#> 3       BacDive_Prod_indole_negative        production      4103
#> 4            BacDive_Enzyme_catalase            enzyme      3921
#> 5     BacDive_Enzyme_urease_negative enzyme|hydrolysis      3912
#> 6 BacDive_Enzyme_leucine_arylamidase    api_zym|enzyme      3648
#>   n_strains_total
#> 1           25750
#> 2            7983
#> 3           12252
#> 4            8262
#> 5           29314
#> 6           17868
```

`BacDive_build_info.tsv` in the same folder records when BacDive was
harvested and which harvester release built the sets.

## Things to keep in mind

- **BacDive is silent for some well-known gut taxa.** It holds no oxygen
  tolerance record for *Faecalibacterium prausnitzii* or *Roseburia
  intestinalis*, so they are not in the oxygen sets. Earlier TaxSEA
  releases listed them from a curated table; the sets now reflect
  BacDive only.
- **Sets are correlated.** A species that ferments glucose often
  ferments mannose, so several substrate sets can be significant for the
  same underlying reason.
- **Very small sets are left out.** Only sets with at least 3 members
  are included.

## Changes in TaxSEA 1.5.5

Set names changed. The most common old names map as follows:

| Before 1.5.5 | From 1.5.5 |
|----|----|
| `BacDive_facultative anaerobe` | `BacDive_Oxygen_facultative_anaerobe` |
| `BacDive_anaerobe` | `BacDive_Oxygen_anaerobe` |
| `BacDive_rod-shaped` | `BacDive_Shape_rod` |
| `BacDive_Utilizes_nitrate` | `BacDive_Uses_nitrate_reduction` (nearest equivalent) |
| `BacDive_Enzyme_alkaline phosphatase` | `BacDive_Enzyme_alkaline_phosphatase` |

The sets are built by
[bacdive_harvester](https://github.com/feargalr/bacdive_harvester),
which documents every curation decision.

Please cite BacDive if you use these sets: Schober et al. BacDive in
2025: the core database for prokaryotic strain data. *Nucleic Acids
Res.* 2025.
