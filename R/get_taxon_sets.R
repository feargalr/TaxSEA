#' Retrieve Taxon Sets from TaxSEA Library
#'
#' Retrieve from the TaxSEA database which taxon sets (metabolite producers and disease signatures) contain a taxon of interest.
#'
#' @param taxon_to_fetch The taxon to search for in the TaxSEA database.
#'
#' @return A character vector containing the names of taxonomic sets where the specified taxon is present.
#'
#' @examples
#' # Retrieve sets for Bifidobacterium longum
#' get_taxon_sets(taxon="Bifidobacterium_longum")
#'
#' @export
get_taxon_sets <- function(taxon_to_fetch=taxon) {
  taxon = get_ncbi_taxon_ids(taxon_to_fetch)
  taxon_sets <- TaxSEA_db
  
  # Using lapply to check for taxon presence
  taxon_presence <- lapply(taxon_sets, function(taxon_set) {
    taxon %in% taxon_set
  })
  # Filter the original list based on presence of the taxon
  filtered_taxon_sets <- taxon_sets[unlist(taxon_presence)]
  return(names(filtered_taxon_sets))
}