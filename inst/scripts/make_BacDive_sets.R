## make_BacDive_sets.R -- how the BacDive_* sets in data/TaxSEA_db.rda and the
## BacDive additions to data/NCBI_ids.rda were produced.
##
## The sets are built outside this package by bacdive_harvester, a staged and
## cached pipeline over the BacDive API:
##
##   https://github.com/feargalr/bacdive_harvester   release v1.0.0
##
## The exact build this package ships is recorded in
## inst/extdata/BacDive_build_info.tsv (harvester commit, BacDive R client
## version, harvest dates, counts), and every shipped set is described in
## inst/extdata/BacDive_set_provenance.tsv (source family, analyte, polarity,
## number of species and strains behind it).
##
## Requirements: BacDive API credentials (free registration at
## https://api.bacdive.dsmz.de) and the BacDive R client, which is distributed
## on R-Forge rather than CRAN:
##
##   install.packages("BacDive", repos = "https://R-Forge.R-project.org")
##
## The harvest is the only step that touches the network. BacDive is updated
## continuously, so a fresh harvest will not reproduce these sets exactly; the
## curation that turns records into sets is deterministic given a harvest.
##
## Run from a shell:
##
##   git clone --branch v1.0.0 https://github.com/feargalr/bacdive_harvester
##   cd bacdive_harvester
##   cp .Renviron.example .Renviron        # add BACDIVE_USER / BACDIVE_PASSWORD
##
##   Rscript R/01a_target_list.R            # genera to harvest, from
##                                          # curatedMetagenomicData species and
##                                          # current TaxSEA_db members
##   Rscript R/01_harvest.R                 # every strain of every target genus
##
##   ## Second pass, repeated until it reports no new genera: finds target
##   ## species reclassified into genera not yet queried (e.g. Agathobacter,
##   ## Phocaeicola), then harvests those genera.
##   Rscript R/01b_resolve_synonyms.R
##   ## merge data/target_genera_extra.csv into data/target_genera.csv, then
##   Rscript R/01_harvest.R
##
##   ## Flatten, aggregate to species, build and validate sets, then write
##   ## TaxSEA's data files. TAXSEA_DATA_DIR points at this package's data/.
##   TAXSEA_DATA_DIR=/path/to/TaxSEA/data Rscript run_all.R
##
## Stage 06 (R/06_merge_taxsea_db.R) writes to the harvester's output/:
##
##   TaxSEA_db.rda   every existing BacDive_* set removed, every other source
##                   kept unchanged, BacDive sets with >= 3 members added
##   NCBI_ids.rda    add-only: species names and former names of the shipped
##                   taxa (with spaces and with underscores, binomials only)
##                   and each shipped taxid as its own name. No existing
##                   taxid is changed; disagreements are listed in
##                   NCBI_ids_conflicts.tsv instead
##   BacDive_set_provenance.tsv, BacDive_build_info.tsv
##
## These were copied into data/ and inst/extdata/ by hand.
