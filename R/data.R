#' TaxSEA Database
#'
#' A dataset containing taxon sets. Each item in the list is a taxon set,
#' and each member within a taxon set is a taxon.
#'
#' @format A list of vectors. Each vector contains character strings
#' representing taxa
#' @source See READ ME.
#' @examples
#' data(TaxSEA_db)
#' all_sets = names(TaxSEA_db)
#' GABA_producers = TaxSEA_db[["MiMeDB_producers_of_GABA"]]
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
#' @source Generated internally.
#' @examples
#' data(TaxSEA_test_data)
#' test_results <- TaxSEA(TaxSEA_test_data)
#'
#' #' NCBI IDs Dataset
#'
#' A dataset for mapping NCBI IDs to species/genus names.
#'
#' @format A named vector where:
#' \describe{
#'   \item{names}{NCBI IDs}
#'   \item{values}{Species or genus names}
#' }
#' @source NCBI
#' @examples
#' data(NCBI_ids)
#' Can look up either with or without spaces
#' NCBI_ids["Bifidobacterium_breve"]
#' NCBI_ids["Bifidobacterium breve"]
