#' Run TaxSEA enrichment (KS) test
#'
#' Internal helper that mirrors the current TaxSEA() enrichment behaviour.
#'
#' @param prep A \code{"TaxSEA_prep"} object from \code{taxsea_prepare()}.
#'
#' @return A data.frame of results (All_databases-style).
#' @keywords internal
#' @noRd
taxsea_KS <- function(prep) {
  if (!inherits(prep, "TaxSEA_prep")) {
    stop("taxsea_KS() expects a 'TaxSEA_prep' object from taxsea_prepare().")
  }
  
  taxon_ranks <- prep$taxon_ranks
  taxon_sets  <- prep$taxon_sets
  
  if (length(taxon_sets) == 0) {
    warning("No taxon sets remain after filtering; returning empty result.")
    return(data.frame(
      taxonSetName = character(0),
      median_rank_of_set_members = numeric(0),
      PValue = numeric(0),
      Test_statistic = numeric(0),
      FDR = numeric(0),
      TaxonSet = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  # Perform KS tests (exactly as in TaxSEA())
  ks_results <- lapply(taxon_sets, function(set) {
    taxon_set_ranks <- taxon_ranks[set]
    nes <- median(taxon_ranks[set])
    ks_result <- stats::ks.test(taxon_set_ranks, taxon_ranks)
    list(nes = nes, ks_result = ks_result)
  })
  
  # Convert set members to legible names (same behaviour as TaxSEA())
  legible_names <- prep$legible_names
  taxon_sets_legible <- lapply(taxon_sets, function(x) legible_names[x])
  
  # Build result table (same columns/order as TaxSEA() All_databases)
  pvals <- vapply(ks_results, function(res) res$ks_result$p.value, numeric(1))
  
  result_df <- data.frame(
    taxonSetName = names(taxon_sets_legible),
    median_rank_of_set_members = vapply(ks_results, function(res) res$nes, numeric(1)),
    PValue = pvals,
    Test_statistic = vapply(ks_results, function(res) res$ks_result$statistic, numeric(1)),
    FDR = stats::p.adjust(pvals, method = "fdr"),
    TaxonSet = vapply(taxon_sets_legible, function(set) paste(set, collapse = ", "),
                      character(1)),
    stringsAsFactors = FALSE
  )
  
  result_df <- result_df[order(result_df$PValue), ]
  rownames(result_df) <- NULL
  result_df
}
