#' Run TaxSEA ORA (Fisher's exact test)
#'
#' Internal helper. Tests over-representation of input taxa (hits) in each taxon set
#' using a one-sided Fisher's exact test (alternative = "greater"). The universe is
#' defined as \code{prep$background_taxa}.
#'
#' @param prep A \code{"TaxSEA_prep"} object from \code{taxsea_prepare(input_taxa=...)}.
#'
#' @return A data.frame with columns consistent with enrichment core results.
#'   For ORA, \code{Test_statistic} is the odds ratio from Fisher's test.
#' @keywords internal
#' @noRd
taxsea_ORA <- function(prep) {
  if (!inherits(prep, "TaxSEA_prep")) {
    stop("taxsea_ORA() expects a 'TaxSEA_prep' object.")
  }
  if (is.null(prep$input_taxa)) {
    stop("ORA prep missing input_taxa. Call taxsea_prepare(input_taxa = ...).")
  }
  if (is.null(prep$background_taxa)) {
    stop("ORA prep missing background_taxa.")
  }
  if (is.null(prep$taxon_sets) || !is.list(prep$taxon_sets)) {
    stop("ORA prep missing taxon_sets.")
  }
  
  hits <- unique(prep$input_taxa)
  universe <- unique(prep$background_taxa)
  
  universe <- universe[!is.na(universe)]
  hits <- hits[hits %in% universe]
  
  if (length(universe) == 0) stop("ORA: universe is empty.")
  if (length(hits) == 0) {
    warning("ORA: none of the input taxa are in the universe; returning empty results.")
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
  
  sets <- prep$taxon_sets
  if (length(sets) == 0) {
    warning("ORA: no taxon sets available after filtering; returning empty results.")
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
  
  N <- length(universe)
  K <- length(hits)
  
  pvals <- numeric(length(sets))
  odds  <- numeric(length(sets))
  names(pvals) <- names(sets)
  
  for (i in seq_along(sets)) {
    set <- unique(sets[[i]])
    set <- set[set %in% universe]
    n <- length(set)
    
    a <- sum(set %in% hits)   # hits in set
    b <- n - a                # non-hits in set
    c <- K - a                # hits not in set
    d <- (N - K) - b          # non-hits not in set
    
    if (any(c(a, b, c, d) < 0)) {
      pvals[i] <- NA_real_
      odds[i]  <- NA_real_
      next
    }
    
    ft <- stats::fisher.test(
      matrix(c(a, b, c, d), nrow = 2, byrow = TRUE),
      alternative = "greater"
    )
    
    pvals[i] <- ft$p.value
    odds[i]  <- unname(ft$estimate)  # odds ratio (may be Inf)
  }
  
  fdr <- stats::p.adjust(pvals, method = "fdr")
  
  taxon_sets_legible <- lapply(sets, function(x) prep$legible_names[x])
  
  result_df <- data.frame(
    taxonSetName = names(sets),
    median_rank_of_set_members = NA_real_,  # ORA has no ranks; keep schema stable
    PValue = pvals,
    Test_statistic = odds,                  # ORA: odds ratio
    FDR = fdr,
    TaxonSet = vapply(taxon_sets_legible, function(set) paste(set, collapse = ", "),
                      character(1)),
    stringsAsFactors = FALSE
  )
  
  result_df <- result_df[order(result_df$PValue), ]
  rownames(result_df) <- NULL
  result_df
}