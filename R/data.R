#' TaxSEA Database
#' A dataset containing taxon sets. Each item in the list is a taxon set,
#' and each member within a taxon set is a taxon.
#'
#' @details
#' Set names start with their source: \code{MiMeDB_}, \code{GutMGene_},
#' \code{GMRepoV2_}, \code{mBodyMap_}, \code{BacDive_} and so on.
#'
#' The \code{BacDive_} sets describe measured phenotypes of cultured strains
#' (oxygen tolerance, Gram stain, substrate use, enzyme activities, growth
#' temperature, salt and pH ranges and others), aggregated to species over
#' every strain BacDive holds. They are built by the bacdive_harvester
#' pipeline (\url{https://github.com/feargalr/bacdive_harvester}); see
#' \code{inst/scripts/make_BacDive_sets.R}. Only sets with at least 3 members
#' are included. Names read \code{BacDive_<Family>_<value>}, optionally
#' followed by an assay context (e.g. \code{_fermentation},
#' \code{_reduction}) and \code{_negative} for species tested and found
#' negative, e.g. \code{BacDive_Oxygen_facultative_anaerobe},
#' \code{BacDive_Uses_nitrate_reduction}.
#'
#' Per-set provenance (source family, number of species and strains behind
#' each set) and a record of the build (harvester release, BacDive client
#' version, harvest dates) are shipped with the package:
#' \preformatted{
#' system.file("extdata", "BacDive_set_provenance.tsv", package = "TaxSEA")
#' system.file("extdata", "BacDive_build_info.tsv", package = "TaxSEA")
#' }
#'
#' @format A list of vectors. Each vector contains character strings
#' representing taxa as NCBI taxonomy IDs.
#' @source See the README. BacDive sets: Schober et al. (2025) BacDive in
#' 2025: the core database for prokaryotic strain data, via
#' \url{https://bacdive.dsmz.de}.
#' @name TaxSEA_db
#' @examples
#' data(TaxSEA_db)
#' all_sets <- names(TaxSEA_db)
#' GABA_producers<-TaxSEA_db[["MiMeDB_producers_of_GABA"]]
"TaxSEA_db"
#'
#' TaxSEA Test Data
#'
#' A dataset containing taxon ranks and taxon IDs.
#'
#' @format A data frame with two columns:
#' \describe{
#'   \item{rank}{Character vector representing taxon ranks}
#'   \item{id}{Character vector representing taxon IDs}
#' }
#' @source See READ ME.
#' @format A data frame with columns 'rank' and 'id' representing taxon
#' ranks and taxon IDs, respectively.
#' @name TaxSEA_test_data
#' @examples
#' data(TaxSEA_test_data)
#' test_results <- TaxSEA(TaxSEA_test_data)
"TaxSEA_test_data"
#'
#' NCBI IDs Dataset
#'
#' A lookup from taxon names to NCBI taxonomy IDs, used to translate the
#' names supplied to \code{TaxSEA()} into the IDs that \code{TaxSEA_db} is
#' built on. Names are present both with spaces and with underscores, and
#' each ID is also present as its own name so that IDs can be supplied
#' directly. Species names, former names and IDs of the \code{BacDive_} set
#' members were added from the BacDive harvest without changing any existing
#' entry.
#'
#' @format A named list where:
#' \describe{
#'   \item{names}{Taxon names (e.g. "Bifidobacterium breve",
#'   "Bifidobacterium_breve") or NCBI IDs}
#'   \item{values}{NCBI taxonomy IDs, as character}
#' }
#' @source NCBI Taxonomy; BacDive (\url{https://bacdive.dsmz.de}).
#' @name NCBI_ids
#' @examples
#' data(NCBI_ids)
#' # Can look up either with or without spaces
#' NCBI_ids["Bifidobacterium_breve"]
#' NCBI_ids["Bifidobacterium breve"]
"NCBI_ids"
