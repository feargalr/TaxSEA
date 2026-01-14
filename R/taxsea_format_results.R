#' Format TaxSEA results into output list
#'
#' Internal helper that reproduces the current TaxSEA() output structure.
#'
#' @param result_df Result table from a testing function.
#' @param prep A \code{"TaxSEA_prep"} object from \code{taxsea_prepare()}.
#'
#' @return A named list of data.frames matching TaxSEA() output.
#' @keywords internal
#' @noRd
taxsea_format_results <- function(result_df, prep) {
  if (!inherits(prep, "TaxSEA_prep")) {
    stop("taxsea_format_results() expects a 'TaxSEA_prep' object from taxsea_prepare().")
  }
  
  # Custom DB: mimic current TaxSEA() behaviour
  if (!prep$used_default_db) {
    # Your current TaxSEA() returns only first 5 cols for custom sets
    out <- result_df
    if (ncol(out) >= 5) out <- out[, seq_len(5), drop = FALSE]
    return(list(custom_sets = out))
  }
  
  if (!is.null(prep$mode) && prep$mode == "ora") {
    # Rename statistic + drop median rank column for ORA
    if ("Test_statistic" %in% colnames(result_df)) {
      colnames(result_df)[colnames(result_df) == "Test_statistic"] <- "Odds_ratio"
    }
    if ("median_rank_of_set_members" %in% colnames(result_df)) {
      result_df$median_rank_of_set_members <- NULL
    }
  }
  
  # Default DB: split into categories (exactly as in TaxSEA())
  metabolites_df <- result_df[grepl("producers_of", result_df$taxonSetName), , drop = FALSE]
  bsdb_df        <- result_df[grepl("bsdb",        result_df$taxonSetName), , drop = FALSE]
  disease_df     <- result_df[!grepl("producers_of|bsdb|BacDive|Colomer2019_MGB",
                                     result_df$taxonSetName), , drop = FALSE]
  bacdive_df     <- result_df[grepl("BacDive",     result_df$taxonSetName), , drop = FALSE]
  GBM_df         <- result_df[grepl("Colomer2019_MGB", result_df$taxonSetName), , drop = FALSE]
  
  # Adjust FDR separately for each subset (same as TaxSEA())
  if (nrow(metabolites_df) > 0) metabolites_df$FDR <- stats::p.adjust(metabolites_df$PValue, method = "fdr")
  if (nrow(disease_df) > 0)     disease_df$FDR     <- stats::p.adjust(disease_df$PValue, method = "fdr")
  if (nrow(bsdb_df) > 0)        bsdb_df$FDR        <- stats::p.adjust(bsdb_df$PValue, method = "fdr")
  if (nrow(bacdive_df) > 0)     bacdive_df$FDR     <- stats::p.adjust(bacdive_df$PValue, method = "fdr")
  if (nrow(GBM_df) > 0)         GBM_df$FDR         <- stats::p.adjust(GBM_df$PValue, method = "fdr")
  
  # Reformat BugSigDB output: add PubMedID column + rename (same as TaxSEA())
  if (nrow(bsdb_df) > 0) {
    pubmed_ids <- bsdb_df$taxonSetName
    pubmed_ids <- sub(".*bsdb:([0-9]+)/.*", "\\1", pubmed_ids)
    pubmed_ids <- ifelse(nchar(pubmed_ids) < 6, NA, pubmed_ids)
    bsdb_df$PubMedID <- pubmed_ids
    
    # Put PubMedID first, then the rest; rename BugSigDB_ID
    bsdb_df <- bsdb_df[, c("PubMedID", colnames(bsdb_df)), drop = FALSE]
    colnames(bsdb_df)[2] <- "BugSigDB_ID"  # taxonSetName becomes BugSigDB_ID
    rownames(bsdb_df) <- NULL
  }
  
  rownames(metabolites_df) <- NULL
  rownames(disease_df)     <- NULL
  rownames(bacdive_df)     <- NULL
  rownames(GBM_df)         <- NULL
  
  list(
    All_databases = result_df,
    Metabolite_producers = metabolites_df,
    Health_associations = disease_df,
    BacDive_bacterial_physiology = bacdive_df,
    BugSigDB = bsdb_df,
    Gut_Brain_Modules_VallesColomer2019 = GBM_df
  )
}
