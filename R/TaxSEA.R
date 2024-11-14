#' TaxSEA: Taxon Set Enrichment Analysis
#'
#' TaxSEA is designed to enable rapid annotation of changes observed in a
#' microbiome association study
#' by testing for enrichment for producers of particular metabolites, or
#' associations with marker taxa
#' for particular diseases. It focuses specifically on human gut
#' microbiome #' associations and uses
#' a Kolmogorov-Smirnov test to test if a particular set of taxa is
#' changed
#' relative to a control group.
#' The input taxon_ranks are log2 fold changes between control and test
#' group
#' (e.g., healthy and IBD).
#'
#' @param taxon_ranks A named vector of log2 fold changes between control
#' and test group.
#' @param lookup_missing Logical value indicating whether to fetch missing
#' NCBI IDs. Default is FALSE.
#' @param min_set_size Minimum size of taxon sets to include in the
#' analysis.
#' Default is 5.
#' @param max_set_size Maximum size of taxon sets to include in the
#' analysis.
#' Default is 100.
#' @return A list of data frames with taxon set enrichment results.
#' @seealso
#' \itemize{
#'   \item \url{https://doi.org/10.1093/nar/gkac868} for MiMeDB
#'   \item \url{https://doi.org/10.1093/nar/gkab1019} for GMrepo
#'   \item \url{https://doi.org/10.1093/nar/gkab786} for gutMGene
#'   \item \url{https://doi.org/10.1038/s41587-023-01872-y} for BugSigDB

#' }
#' @examples
#' data("TaxSEA_test_data")
#' taxsea_results <- TaxSEA(TaxSEA_test_data)
#'
#' @import stats
#' @importFrom utils data
#' @export
TaxSEA <- function(taxon_ranks, lookup_missing = FALSE,
                   min_set_size = 5, max_set_size = 100) {
  # Function implementation
}
#Declaring global availabilability
utils::globalVariables(c("TaxSEA_db","Rank", "TaxonSet", "InSet",
                         "log10FDR",
                         "taxonSetName","fetched_ids","URLencode",
                         "median_rank",
                         "xml_data","NCBI_ids","TaxSEA_test_data"))

TaxSEA <- function(taxon_ranks, lookup_missing = FALSE,
                   min_set_size = 5, max_set_size = 100) {

data("TaxSEA_db", package = "TaxSEA", envir = environment())
data("NCBI_ids", package = "TaxSEA", envir = environment())

  taxon_sets <- TaxSEA_db

  if (lookup_missing == TRUE) {

  #Relabel taxa with NCBI IDs
  ids2fetch <- names(taxon_ranks[!(names(taxon_ranks) %in%
                                     names(NCBI_ids))])
  if (length(ids2fetch) > 0){
    message("Fetching some NCBI IDs for input taxa. Please wait")
    fetched_ids <- (get_ncbi_taxon_ids(ids2fetch))
    if (length(unlist(fetched_ids)) > 0){
      NCBI_ids <- c(NCBI_ids,unlist(fetched_ids))
    }
  }
}

  taxon_ranks <- taxon_ranks[names(taxon_ranks) %in% names(NCBI_ids)]
  original_ranks <- taxon_ranks
  names(taxon_ranks) <- NCBI_ids[names(taxon_ranks)]
  NCBI_ids <- NCBI_ids[names(NCBI_ids) %in% names(original_ranks)]
  NCBI_ids <- NCBI_ids[!duplicated(NCBI_ids)]
  legible_names <- names(NCBI_ids)
  names(legible_names)  <- NCBI_ids[legible_names]

  if(length(taxon_ranks) < 3) {
    warning("Error: Very few taxa provided")
  }

  # Filter taxon_ranks and taxon_sets
  taxon_ranks <- taxon_ranks[names(taxon_ranks) %in%
                               unique(unlist(taxon_sets))]
  taxon_sets <- lapply(taxon_sets, function(taxon_set) {
    return(unique(taxon_set[taxon_set %in% names(taxon_ranks)]))})
  set_sizes <- vapply(taxon_sets, function(taxon_set)
    { length(taxon_set) },numeric(1))
  taxon_sets <- taxon_sets[set_sizes>= min_set_size &
                             set_sizes <= max_set_size]
  taxon_ranks <- taxon_ranks[names(taxon_ranks) %in%
                               unique(unlist(taxon_sets))]

  if (!is.vector(taxon_ranks) || !is.list(taxon_sets)) {
    warning("Input doesn't look correct")
  }

  taxon_sets <- lapply(taxon_sets, function(taxon_set) {
    intersect(names(taxon_ranks), taxon_set)
  })

  # Perform KS tests
  ks_results <- lapply(taxon_sets, function(taxon_set) {
    taxon_set_ranks <- taxon_ranks[taxon_set]


    nes <- median(taxon_ranks[taxon_set])

    ks_result <- (ks.test(taxon_set_ranks, taxon_ranks))

    return(list(nes = nes, ks_result = ks_result))
  })
  taxon_sets <- lapply(taxon_sets,function(X){legible_names[X]})
  result_df <- data.frame(
    taxonSetName <- names(taxon_sets),
    median_rank<-vapply(ks_results, function(res) res$nes, numeric(1)),
    PValue<-vapply(ks_results, function(res) res$ks_result$p.value,
                    numeric(1) ),
    FDR<-p.adjust(vapply(ks_results, function(res) res$ks_result$p.value,
                          numeric(1)),
                   method = "fdr"),
    TaxonSet<-vapply(taxon_sets, function(taxon_set)
      paste(taxon_set,collapse = ", "),
                      character(1))
  )
  result_df <- result_df[order(result_df$PValue,decreasing = FALSE),]
  colnames(result_df) = c("taxonSetName","median_rank_of_set_members",
                          "PValue","FDR","TaxonSet")
  metabolites_df <-
    result_df[grepl("producers_of",result_df$taxonSetName),]
  bsdb_df <- result_df[grepl("bsdb",result_df$taxonSetName),]

  disease_df <- result_df[!grepl("producers_of",result_df$taxonSetName),]
  disease_df <- disease_df[!grepl("bsdb",disease_df$taxonSetName),]

  metabolites_df$FDR <- p.adjust(metabolites_df$PValue,method = "fdr")
  disease_df$FDR <- p.adjust(disease_df$PValue,method = "fdr")
  bsdb_df$FDR <- p.adjust(bsdb_df$PValue,method = "fdr")

  res_list <- list(Metabolite_producers = metabolites_df,
                  Health_associations = disease_df,
                  BugSigdB = bsdb_df)
  return(res_list)
  rm(TaxSEA_db)
  rm(NCBI_ids)

}
