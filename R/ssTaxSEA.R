#' Single-Sample Taxon Set Enrichment Analysis
#'
#' Computes per-sample enrichment scores for taxon sets using an
#' ssGSEA-style approach. Counts are CLR-transformed, then each
#' taxon is z-scored across samples to capture between-sample
#' variation. Per-sample enrichment is then computed by ranking
#' each sample's z-scores and applying a weighted running-sum
#' statistic against taxon sets from the TaxSEA database.
#'
#' @param counts A numeric matrix, data.frame, or
#'   \code{SummarizedExperiment}/\code{TreeSummarizedExperiment} object.
#'   For matrix/data.frame input: rows are taxa, columns are samples,
#'   and row names must be taxon names (e.g. species names or NCBI IDs).
#'   For SummarizedExperiment input: the first assay is used and
#'   \code{rownames()} provide taxon identifiers.
#' @param lookup_missing Logical indicating whether to fetch missing
#'   NCBI IDs via the NCBI API. Default is FALSE.
#' @param min_set_size Minimum size of taxon sets to include.
#'   Default is 5.
#' @param max_set_size Maximum size of taxon sets to include.
#'   Default is 300.
#' @param custom_db A user-provided list of taxon sets. If NULL
#'   (default), the built-in TaxSEA database is used (excluding
#'   BugSigDB).
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{scores}{A matrix (samples x taxon sets) of enrichment
#'       scores. Positive scores indicate the set taxa tend to have
#'       higher abundance in that sample relative to the cohort.}
#'     \item{pvalues}{A matrix (samples x taxon sets) of KS test
#'       p-values for each sample-set combination.}
#'   }
#'
#' @details
#' The approach works as follows:
#' \enumerate{
#'   \item Raw counts are CLR-transformed (centered log-ratio with
#'     pseudocount of 0.5 for zeros).
#'   \item Each taxon is then z-scored across all samples, so that
#'     values represent how much higher or lower a taxon is in a
#'     given sample relative to the cohort mean.
#'   \item For each sample, the z-scores are ranked and an
#'     ssGSEA-style weighted running-sum enrichment score is
#'     computed for each taxon set.
#'   \item A KS test p-value is also computed per sample per set.
#' }
#'
#' This cohort-relative approach ensures that taxon sets which are
#' consistently elevated in a subset of samples (e.g. disease
#' samples) will produce high enrichment scores in those samples,
#' even if the taxa are not the most abundant within any single
#' sample.
#'
#' @examples
#' \dontrun{
#' # From a count matrix (taxa x samples)
#' counts <- matrix(rpois(500, lambda = 10), nrow = 50, ncol = 10)
#' rownames(counts) <- paste0("Taxon_", seq_len(50))
#' colnames(counts) <- paste0("Sample_", seq_len(10))
#' res <- ssTaxSEA(counts, custom_db = list(
#'   set1 = paste0("Taxon_", 1:10),
#'   set2 = paste0("Taxon_", 20:30)
#' ), min_set_size = 2)
#' head(res$scores)
#' head(res$pvalues)
#' }
#'
#' @export
ssTaxSEA <- function(counts,
                     lookup_missing = FALSE,
                     min_set_size = 5,
                     max_set_size = 300,
                     custom_db = NULL) {

  # --- Handle SummarizedExperiment / TreeSummarizedExperiment ---
  if (methods::is(counts, "SummarizedExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("The 'SummarizedExperiment' package is required when ",
           "providing a SummarizedExperiment or ",
           "TreeSummarizedExperiment object.")
    }
    message("Detected SummarizedExperiment input. ",
            "Extracting counts from assay().")
    counts <- as.matrix(SummarizedExperiment::assay(counts))
  }

  # --- Validate input ---
  if (is.data.frame(counts)) counts <- as.matrix(counts)
  if (!is.matrix(counts) || !is.numeric(counts)) {
    stop("'counts' must be a numeric matrix, data.frame, or ",
         "SummarizedExperiment.")
  }
  if (is.null(rownames(counts))) {
    stop("'counts' must have row names (taxon names).")
  }
  if (ncol(counts) < 1) stop("'counts' must have at least one sample.")

  if (any(grepl("\\[|\\]", rownames(counts)))) {
    stop("Taxon names contain square brackets [ ]. ",
         "Please remove or rename these entries before running ssTaxSEA.")
  }

  # --- Prepare taxon sets and ID mapping ---
  prep <- ss_prepare(
    taxon_names = rownames(counts),
    lookup_missing = lookup_missing,
    min_set_size = min_set_size,
    max_set_size = max_set_size,
    custom_db = custom_db
  )

  taxon_sets <- prep$taxon_sets
  id_map <- prep$id_map

  if (length(taxon_sets) == 0) {
    warning("No taxon sets remain after filtering. ",
            "Returning empty result.")
    return(list(
      scores = matrix(nrow = ncol(counts), ncol = 0,
                      dimnames = list(colnames(counts), character(0))),
      pvalues = matrix(nrow = ncol(counts), ncol = 0,
                       dimnames = list(colnames(counts), character(0)))
    ))
  }

  # --- Map row names to NCBI IDs (if using default DB) ---
  if (is.null(custom_db)) {
    mapped_rows <- rownames(counts) %in% names(id_map)
    counts <- counts[mapped_rows, , drop = FALSE]
    rownames(counts) <- id_map[rownames(counts)]
  }

  # Keep only taxa that appear in at least one set
  all_set_taxa <- unique(unlist(taxon_sets))
  counts <- counts[rownames(counts) %in% all_set_taxa, , drop = FALSE]

  if (nrow(counts) < 3) {
    warning("Very few taxa overlap with taxon sets (", nrow(counts),
            "). Results may be unreliable.")
  }

  # --- CLR transform ---
  clr_mat <- ss_clr(counts)

  # --- Z-score each taxon across samples (cohort-relative) ---
  if (ncol(clr_mat) < 3) {
    warning("Fewer than 3 samples: z-scoring across samples may ",
            "be unreliable.")
  }
  row_means <- rowMeans(clr_mat)
  row_sds <- apply(clr_mat, 1, stats::sd)

  # Handle zero-variance taxa (constant across all samples)
  zero_sd <- row_sds == 0
  if (any(zero_sd)) {
    message(sum(zero_sd), " taxa have zero variance across samples ",
            "and will be assigned a z-score of 0.")
  }
  row_sds[zero_sd] <- 1  # avoid division by zero; result will be 0

  z_mat <- (clr_mat - row_means) / row_sds

  # --- Compute enrichment scores and p-values per sample ---
  n_samples <- ncol(z_mat)
  n_sets <- length(taxon_sets)
  set_names <- names(taxon_sets)

  score_mat <- matrix(NA_real_, nrow = n_samples, ncol = n_sets,
                      dimnames = list(colnames(z_mat), set_names))
  pval_mat <- matrix(NA_real_, nrow = n_samples, ncol = n_sets,
                     dimnames = list(colnames(z_mat), set_names))

  taxa_ids <- rownames(z_mat)

  for (j in seq_len(n_samples)) {
    sample_vals <- z_mat[, j]
    ranked_vals <- rank(sample_vals, ties.method = "average")
    names(ranked_vals) <- taxa_ids

    for (k in seq_len(n_sets)) {
      set_members <- intersect(taxon_sets[[k]], taxa_ids)
      if (length(set_members) < 2) next

      score_mat[j, k] <- ss_enrichment_score(ranked_vals, set_members)

      set_vals <- sample_vals[set_members]
      ks <- suppressWarnings(stats::ks.test(set_vals, sample_vals))
      pval_mat[j, k] <- ks$p.value
    }
  }

  list(scores = score_mat, pvalues = pval_mat)
}


#' Prepare taxon sets and ID mapping for ssTaxSEA
#'
#' @param taxon_names Character vector of taxon names (row names
#'   from count matrix).
#' @param lookup_missing Logical; fetch missing NCBI IDs.
#' @param min_set_size Minimum set size.
#' @param max_set_size Maximum set size.
#' @param custom_db Optional custom database.
#'
#' @return A list with \code{taxon_sets} (filtered) and \code{id_map}
#'   (named character vector mapping taxon names to NCBI IDs).
#' @keywords internal
#' @noRd
ss_prepare <- function(taxon_names,
                       lookup_missing = FALSE,
                       min_set_size = 5,
                       max_set_size = 300,
                       custom_db = NULL) {

  # --- Load database ---
  if (is.null(custom_db)) {
    utils::data("TaxSEA_db", package = "TaxSEA", envir = environment())
    taxon_sets <- TaxSEA_db

    # Load NCBI ID mapping
    utils::data("NCBI_ids", package = "TaxSEA", envir = environment())

    if (lookup_missing) {
      ids2fetch <- taxon_names[!(taxon_names %in% names(NCBI_ids))]
      if (length(ids2fetch) > 0) {
        fetched_ids <- get_ncbi_taxon_ids(ids2fetch)
        if (length(unlist(fetched_ids)) > 0) {
          NCBI_ids <- c(NCBI_ids, unlist(fetched_ids))
        }
      }
    }

    # Map taxon names to NCBI IDs
    matched <- taxon_names[taxon_names %in% names(NCBI_ids)]
    id_map <- NCBI_ids[matched]
    id_map <- id_map[!duplicated(id_map)]

    # Filter sets to taxa present in our data
    mapped_ncbi <- unname(id_map)
    taxon_sets <- lapply(taxon_sets, function(set) {
      unique(set[set %in% mapped_ncbi])
    })
  } else {
    taxon_sets <- custom_db
    if (!is.list(taxon_sets)) {
      stop("Custom database must be a list of taxon sets.")
    }
    # For custom DB, taxon names are used directly
    taxon_sets <- lapply(taxon_sets, function(set) {
      intersect(set, taxon_names)
    })
    id_map <- stats::setNames(taxon_names, taxon_names)
  }

  # Filter by set size
  set_sizes <- vapply(taxon_sets, length, numeric(1))
  taxon_sets <- taxon_sets[set_sizes >= min_set_size &
                             set_sizes <= max_set_size]

  list(taxon_sets = taxon_sets, id_map = id_map)
}


#' CLR transform a count matrix
#'
#' Applies centered log-ratio transformation per sample (column).
#' Zeros are replaced with a pseudocount of 0.5 before transformation.
#'
#' @param mat Numeric matrix (taxa x samples).
#' @return CLR-transformed matrix of the same dimensions.
#' @keywords internal
#' @noRd
ss_clr <- function(mat) {
  mat[mat == 0] <- 0.5
  log_mat <- log(mat)
  geom_means <- colMeans(log_mat)
  sweep(log_mat, 2, geom_means, FUN = "-")
}


#' Compute ssGSEA enrichment score
#'
#' Calculates a single-sample enrichment score using the weighted
#' running-sum statistic (Barbie et al., 2009).
#'
#' @param ranked_vals Named numeric vector of ranks for one sample.
#' @param set_members Character vector of taxon IDs in the set.
#' @param alpha Weighting exponent. Default 0.25.
#' @return Numeric enrichment score.
#' @keywords internal
#' @noRd
ss_enrichment_score <- function(ranked_vals, set_members, alpha = 0.25) {
  N <- length(ranked_vals)
  n_set <- length(set_members)

  sorted_idx <- order(ranked_vals, decreasing = TRUE)
  sorted_names <- names(ranked_vals)[sorted_idx]
  sorted_ranks <- ranked_vals[sorted_idx]

  is_in_set <- sorted_names %in% set_members

  # Weighted step-up for hits
  hit_weights <- abs(sorted_ranks[is_in_set])^alpha
  norm_hit <- sum(hit_weights)

  # Step-down for misses
  miss_penalty <- 1 / (N - n_set)

  # Running sum
  running_sum <- numeric(N)
  cumulative <- 0
  for (i in seq_len(N)) {
    if (is_in_set[i]) {
      cumulative <- cumulative + abs(sorted_ranks[i])^alpha / norm_hit
    } else {
      cumulative <- cumulative - miss_penalty
    }
    running_sum[i] <- cumulative
  }

  sum(running_sum)
}
