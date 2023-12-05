library(bugsigdbr) #Installable via Bioconductor
library(TaxSEA)
bsdb <- importBugSigDB() #Import database
mp.sigs <- getSignatures(bsdb, tax.id.type = "ncbi") #get a list of signatures with NCBI IDs

TaxSEA_db = c(TaxSEA_db,mp.sigs)

my_results =  TaxSEA(taxon_ranks = TaxSEA_test_data)
