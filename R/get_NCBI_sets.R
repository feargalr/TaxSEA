#' Retrieve NCBI Taxonomy IDs for a list of taxon names
#'
#' This function takes a vector of taxon names and returns a vector of
#' NCBI taxonomy IDs
#' by querying the NCBI Entrez API.
#' @param taxon_names A character vector of taxon names
#' @return A character vector of NCBI taxonomy IDs corresponding to the
#' input taxon names
#' @examples
#' taxon_names <- c("Escherichia coli", "Staphylococcus aureus")
#' taxon_ids <- get_ncbi_taxon_ids(taxon_names)
#' @export
get_ncbi_taxon_ids <- function(taxon_names) {
  get_ncbi_taxon_id <- function(taxon_name) {
    # Define the URL for the API call
    base_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    esearch_url <- paste0(base_url, "esearch.fcgi?db=taxonomy&term=",
                          URLencode(taxon_name, reserved = TRUE),
                          "&retmode=xml")

    # Retrieve the XML data
    xml_data <- readLines(url(esearch_url, "r"))

    # Find the line containing the taxonomy ID
    id_line <- grep("<Id>", xml_data, value = TRUE)

    # Extract the taxonomy ID
    taxon_id <- gsub("<Id>|</Id>", "", id_line)

    # Return the taxonomy ID
    if (length(taxon_id) > 0) {
      return(taxon_id)
    } else {
      warning("No taxonomy ID found for taxon:", taxon_name)
      return(NULL)
    }
  }
  data("TaxSEA_db", package = "TaxSEA", envir = environment())
data("NCBI_ids", package = "TaxSEA", envir = environment())
  ids2fetch<-taxon_names[!taxon_names %in% names(NCBI_ids)]
  taxon_names<-taxon_names[taxon_names %in% names(NCBI_ids)]
  local_ids<-unlist(NCBI_ids[taxon_names])
  # Retrieve the taxonomy IDs for each taxon name in the vector
  if (length(ids2fetch) > 0){
    message("Fetching some NCBI IDs for input taxa. Please wait")
    taxon_ids <- vapply(ids2fetch, get_ncbi_taxon_id,character(1))
    taxon_ids<-c(local_ids,unlist(taxon_ids))
  } else {
    taxon_ids<-local_ids
  }

  return(taxon_ids)
  close(xml_data)

}

