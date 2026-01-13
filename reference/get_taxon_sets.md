# Retrieve Taxon Sets from TaxSEA Library

Retrieve from the TaxSEA database which taxon sets (metabolite producers
and disease signatures) contain a taxon of interest.

## Usage

``` r
get_taxon_sets(taxon_to_fetch = taxon)
```

## Arguments

- taxon_to_fetch:

  The taxon to search for in the TaxSEA database.

## Value

A character vector containing the names of taxonomic sets where the
specified taxon is present.

## Examples

``` r
# Retrieve sets for Bifidobacterium longum
get_taxon_sets(taxon="Bifidobacterium_longum")
#>  [1] "GMRepoV2_Decreased_in_Adenoma"                                        
#>  [2] "GMRepoV2_Decreased_in_COVID-19"                                       
#>  [3] "GMRepoV2_Increased_in_Spondylitis,_Ankylosing"                        
#>  [4] "GMRepoV2_Increased_in_Tuberculosis"                                   
#>  [5] "GutMGene_producers_of_4-Hydroxyphenylacetic_acid"                     
#>  [6] "GutMGene_producers_of_5-Acetamido-2-azaniumylpentanoate"              
#>  [7] "GutMGene_producers_of_Acetate"                                        
#>  [8] "GutMGene_producers_of_Arctigenin"                                     
#>  [9] "GutMGene_producers_of_Daidzein"                                       
#> [10] "GutMGene_producers_of_Dihydrocaffeic_acid"                            
#> [11] "GutMGene_producers_of_Genipin"                                        
#> [12] "GutMGene_producers_of_Genistein"                                      
#> [13] "GutMGene_producers_of_Ginsenoside-Rd"                                 
#> [14] "GutMGene_producers_of_Glycitein"                                      
#> [15] "GutMGene_producers_of_Hydroquinone"                                   
#> [16] "GutMGene_producers_of_Indole-3-acetic_acid"                           
#> [17] "GutMGene_producers_of_Indole-3-lactic_acid"                           
#> [18] "GutMGene_producers_of_Kaempferol"                                     
#> [19] "GutMGene_producers_of_Lactate"                                        
#> [20] "GutMGene_producers_of_Loganetin"                                      
#> [21] "GutMGene_producers_of_Oxyberberine"                                   
#> [22] "GutMGene_producers_of_Phenylalanine"                                  
#> [23] "MiMeDB_producers_of_3IAA"                                             
#> [24] "MiMeDB_producers_of_Acetic_acid"                                      
#> [25] "MiMeDB_producers_of_Butyric_acid"                                     
#> [26] "MiMeDB_producers_of_Dimethylamine_"                                   
#> [27] "MiMeDB_producers_of_Folic_acid"                                       
#> [28] "MiMeDB_producers_of_GABA"                                             
#> [29] "MiMeDB_producers_of_Glycocholic_acid"                                 
#> [30] "MiMeDB_producers_of_Hydrogen_sulfide"                                 
#> [31] "MiMeDB_producers_of_Indole"                                           
#> [32] "MiMeDB_producers_of_Propionic_acid"                                   
#> [33] "MiMeDB_producers_of_Taurocholic_acid"                                 
#> [34] "BacDive_anaerobe"                                                     
#> [35] "Mucin_utilisers"                                                      
#> [36] "Valles-Colomer2019_MGB005.Tryptophan.synthesis"                       
#> [37] "Valles-Colomer2019_MGB006.Glutamate.synthesis.I"                      
#> [38] "Valles-Colomer2019_MGB007.Glutamate.synthesis.II"                     
#> [39] "Valles-Colomer2019_MGB029.ClpB.(ATP-dependent.chaperone.protein)"     
#> [40] "Valles-Colomer2019_MGB032.Quinolinic.acid.synthesis"                  
#> [41] "Valles-Colomer2019_MGB033.Quinolinic.acid.degradation"                
#> [42] "Valles-Colomer2019_MGB035.Isovaleric.acid.synthesis.II.(KADC.pathway)"
#> [43] "Valles-Colomer2019_MGB036.S-Adenosylmethionine.(SAM).synthesis"       
#> [44] "Valles-Colomer2019_MGB043.Acetate.synthesis.I"                        
```
