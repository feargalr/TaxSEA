#' TaxSEA: Taxon Set Enrichment Analysis
#'
#' TaxSEA enables rapid annotation of changes by testing for enrichment 
#' of pre-defined taxon sets. 
#'
#' @param taxon_ranks A named vector of log2 fold changes between control 
#' and test groups.
#' @param lookup_missing Logical indicating whether to fetch missing 
#' NCBI IDs. Default is FALSE.
#' @param min_set_size Minimum size of taxon sets to include in the 
#' analysis. Default is 5.
#' @param max_set_size Maximum size of taxon sets to include in the 
#' analysis. Default is 100.
#' @param custom_db A user-provided list of taxon sets. 
#' If NULL (default), the built-in database is used.
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
                   min_set_size = 5, max_set_size = 100,
                   custom_db = NULL) {
  if (any(grepl("\\[|\\]", names(taxon_ranks)))) {
    stop("Taxon names contain square brackets [ ]. Please remove or rename these entries before running TaxSEA.")
  }
  
  
  # Load built-in database unless a custom one is provided
  if (is.null(custom_db)) {
    data("TaxSEA_db", package = "TaxSEA", envir = environment())
    taxon_sets <- TaxSEA_db
    
    # Append BugSigDB dynamically
    if (requireNamespace("bugsigdbr", quietly = TRUE)) {
      bsdb <- bugsigdbr::importBugSigDB()
      mp.sigs <- bugsigdbr::getSignatures(bsdb, tax.id.type = "ncbi")
      bugsigdb_list <- stack(mp.sigs)
      names(bugsigdb_list) <- c("Species", "MSID")
      bugsigdb_list <- split(bugsigdb_list$Species, bugsigdb_list$MSID)
      names(bugsigdb_list) <- paste0("bsdb_", names(bugsigdb_list))
      taxon_sets <- c(taxon_sets, bugsigdb_list)
    } else {
      warning("bugsigdbr not installed; skipping BugSigDB integration.")
    }
  } else {
    taxon_sets <- custom_db
    if (!is.list(taxon_sets)) {
      stop("Custom database must be a list of taxon sets.")
    }
  }
  
  # Load NCBI ID mappings only if using the default database
  if (is.null(custom_db)) {
    data("NCBI_ids", package = "TaxSEA", envir = environment())
    
    if (lookup_missing) {
      ids2fetch <- names(taxon_ranks[!(names(taxon_ranks) 
                                       %in% names(NCBI_ids))])
      if (length(ids2fetch) > 0) {
        fetched_ids <- get_ncbi_taxon_ids(ids2fetch)
        if (length(unlist(fetched_ids)) > 0) {
          NCBI_ids <- c(NCBI_ids, unlist(fetched_ids))
        }
      }
    }
    
    # Convert taxon names to NCBI IDs
    taxon_ranks <- taxon_ranks[names(taxon_ranks) %in% names(NCBI_ids)]
    original_ranks <- taxon_ranks
    names(taxon_ranks) <- NCBI_ids[names(taxon_ranks)]
    NCBI_ids <- NCBI_ids[names(NCBI_ids) %in% names(original_ranks)]
    NCBI_ids <- NCBI_ids[!duplicated(NCBI_ids)]
    legible_names <- names(NCBI_ids)
    names(legible_names) <- NCBI_ids[legible_names]
  } else {
    # If custom database is used, retain original names
    legible_names <- names(taxon_ranks)
  }
  
  # Ensure sufficient data for testing
  if (length(taxon_ranks) < 3) {
    warning("Error: Very few taxa provided; results may be unreliable.")
  }
  
  # Filter taxon_ranks and taxon_sets
  taxon_ranks <- 
    taxon_ranks[names(taxon_ranks) %in% unique(unlist(taxon_sets))]
  taxon_sets <- 
    lapply(taxon_sets, function(set) 
      unique(set[set %in% names(taxon_ranks)]))
  set_sizes <- vapply(taxon_sets, length, numeric(1))
  taxon_sets <- 
    taxon_sets[set_sizes >= min_set_size & set_sizes <= max_set_size]
  taxon_ranks <- 
    taxon_ranks[names(taxon_ranks) %in% unique(unlist(taxon_sets))]
  
  taxon_sets <- 
    lapply(taxon_sets, function(set) intersect(names(taxon_ranks), set))
  
  # Perform KS tests
  ks_results <- lapply(taxon_sets, function(set) {
    taxon_set_ranks <- taxon_ranks[set]
    nes <- median(taxon_ranks[set])
    ks_result <- ks.test(taxon_set_ranks, taxon_ranks)
    list(nes = nes, ks_result = ks_result)
  })
  
  taxon_sets <- lapply(taxon_sets, function(X) legible_names[X])
  
  result_df <- data.frame(
    taxonSetName = 
      names(taxon_sets),
    median_rank = 
      vapply(ks_results, function(res) 
        res$nes, numeric(1)),
    PValue = 
      vapply(ks_results, function(res) 
        res$ks_result$p.value, numeric(1)),
    Test_statistic = 
      vapply(ks_results, function(res) 
        res$ks_result$statistic, numeric(1)),
    FDR = 
      p.adjust(vapply(ks_results, function(res) 
        res$ks_result$p.value, numeric(1)), method = "fdr"),
    TaxonSet = vapply(taxon_sets, function(set) 
      paste(set, collapse = ", "), character(1))
  )
  result_df <- result_df[order(result_df$PValue), ]
  colnames(result_df) <- c("taxonSetName", 
                           "median_rank_of_set_members", 
                           "PValue","Test_statistic", 
                           "FDR", "TaxonSet")
  
  if (is.null(custom_db)) {
    
  # Subset results into categories
  metabolites_df <- 
    result_df[grepl("producers_of", result_df$taxonSetName), ]
  bsdb_df <- 
    result_df[grepl("bsdb", result_df$taxonSetName), ]
  disease_df <- 
    result_df[!grepl("producers_of|bsdb", result_df$taxonSetName), ]
  
  # Adjust FDR separately for each subset
  
  metabolites_df$FDR <- p.adjust(metabolites_df$PValue, method = "fdr")
  disease_df$FDR <- p.adjust(disease_df$PValue, method = "fdr")
  bsdb_df$FDR <- p.adjust(bsdb_df$PValue, method = "fdr")
  
  res_list <- list(
    Metabolite_producers = metabolites_df,
    Health_associations = disease_df,
    BugSigDB = bsdb_df
  )
  
  return(res_list)}
  res_list <- list(custom_sets = result_df[,seq_len(5)])
  return(res_list)
}

