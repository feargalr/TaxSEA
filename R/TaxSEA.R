#' TaxSEA: Taxon Set Enrichment Analysis
#'
#' Modular TaxSEA implementation supporting enrichment (KS) and 
#' ORA (Fisher).
#' Provide either \code{taxon_ranks} for enrichment or 
#' \code{input_taxa} for ORA.
#'
#' @param taxon_ranks Named numeric vector of statistics 
#' (e.g., log2 fold changes).
#'   Required for enrichment.
#' @param input_taxa Character vector of taxa to treat as 
#' "hits"/selected taxa.
#'   Required for ORA.
#' @param mode Character. One of \code{"enrichment"} or \code{"ora"}.
#'   If NULL, inferred from which input is provided.
#' @param lookup_missing Logical indicating whether to fetch
#'  missing NCBI IDs.
#'   Default is FALSE.
#' @param min_set_size Minimum size of taxon sets to include
#'  in the analysis.
#'   Default is 5.
#' @param max_set_size Maximum size of taxon sets to include
#'  in the analysis.
#'   Default is 100.
#' @param custom_db A user-provided list of taxon sets. 
#' If NULL (default),
#'   the built-in database is used.
#'
#' @return A list of data frames with taxon set results.
#'
#' @examples
#' data("TaxSEA_test_data")
#' res <- TaxSEA(taxon_ranks = TaxSEA_test_data)
#' head(res$All_databases)
#'
#' # ORA example (toy): treat taxa with positive values as "hits"
#' hits <- names(TaxSEA_test_data)
#' res_ora <- TaxSEA(input_taxa = hits, mode = "ora")
#' head(res_ora$All_databases)
#'
#' @export
TaxSEA <- function(taxon_ranks = NULL,
                    input_taxa = NULL,
                    mode = NULL,
                    lookup_missing = FALSE,
                    min_set_size = 5,
                    max_set_size = 300,
                    custom_db = NULL) {
  
  # Infer mode if not provided (strict, no magic)
  if (is.null(mode)) {
    if (!is.null(taxon_ranks) && is.null(input_taxa)) {
      mode <- "enrichment"
    } else if (is.null(taxon_ranks) && !is.null(input_taxa)) {
      mode <- "ora"
    } else if (!is.null(taxon_ranks) && !is.null(input_taxa)) {
      stop("Provide either 'taxon_ranks' or 'input_taxa', not both")
    } else {
      stop("Provide either 'taxon_ranks' (enrichment) or
           'input_taxa' (ORA).")
    }
  }
  
  mode <- match.arg(mode, c("enrichment", "ora"))
  
  # Enforce consistency between mode and supplied input
  if (mode == "enrichment") {
    if (is.null(taxon_ranks)) stop("Enrichment mode requires 
                                   'taxon_ranks'.")
    if (!is.null(input_taxa)) stop("Enrichment mode does not accept 
                                   'input_taxa'.")
  } else {
    if (is.null(input_taxa)) stop("ORA mode requires 
                                  'input_taxa'.")
    if (!is.null(taxon_ranks)) stop("ORA mode does not 
                                    accept 'taxon_ranks'.")
  }
  
  prep <- taxsea_prepare(
    taxon_ranks = taxon_ranks,
    input_taxa = input_taxa,
    lookup_missing = lookup_missing,
    min_set_size = min_set_size,
    max_set_size = max_set_size,
    custom_db = custom_db
  )
  
  core <- switch(
    mode,
    enrichment = taxsea_KS(prep),   # your enrichment core
    ora = taxsea_ORA(prep)          # your ORA core (Fisher + OR)
  )
  
  taxsea_format_results(core, prep)
}