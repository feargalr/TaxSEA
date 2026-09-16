#' Single-Sample Taxon Set Enrichment Analysis
#'
#' Computes a per-sample score for each taxon set as the mean
#' centered log-ratio (CLR) of the set's members within that sample.
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
#' @param pseudocount Numeric value added to every count before the
#'   log transform, to handle zeros. Default is 0.5. This choice
#'   affects the scores, so report it alongside your results.
#'
#' @return A numeric matrix of scores with taxon sets as rows and
#'   samples as columns. Higher values indicate that the set's members
#'   are more abundant in that sample, relative to the average taxon
#'   in that same sample.
#'
#' @details
#' For each sample the counts are centered log-ratio transformed, and
#' the score for a taxon set is the mean of the CLR values of the set's
#' members:
#'
#' \deqn{score_{S,j} = \frac{1}{|S|} \sum_{i \in S} clr_{ij}}
#'
#' The CLR is computed across \emph{all} taxa supplied in \code{counts},
#' before any subsetting to set members. This matters: restricting the
#' matrix to set members first would make the geometric mean that the
#' CLR divides by depend on the sets being tested, which manufactures
#' apparent signal in sets that have none.
#'
#' Because the score is a within-sample quantity, \code{ssTaxSEA} does
#' not need a cohort and is well defined for a single sample.
#'
#' @section Interpreting the scores:
#' Scores are comparable \strong{across samples within a taxon set}:
#' a higher score in sample A than sample B means the set's members
#' make up more of sample A's community.
#'
#' Scores are \strong{not} comparable across taxon sets. A set of
#' abundant taxa will score higher than a set of rare taxa in every
#' sample, regardless of biology, simply because its members are more
#' abundant. When visualising the matrix, center or scale the rows
#' first (for example \code{t(scale(t(scores)))}), and when comparing
#' groups, compare a single set's scores between groups rather than
#' comparing different sets to each other.
#'
#' No p-values are returned. The score is a descriptive statistic;
#' to test a hypothesis, compare scores across samples using a test
#' appropriate to your design (for example \code{\link[stats]{wilcox.test}}
#' or a linear model on one row of the returned matrix).
#'
#' @examples
#' # Toy count matrix: 30 taxa x 8 samples
#' set.seed(42)
#' counts <- matrix(rpois(240, lambda = 10), nrow = 30, ncol = 8)
#' rownames(counts) <- paste0("Taxon_", seq_len(30))
#' colnames(counts) <- paste0("Sample_", seq_len(8))
#'
#' # Raise one set's members in the last four samples
#' counts[1:6, 5:8] <- counts[1:6, 5:8] * 5
#'
#' sets <- list(
#'   elevated_set = paste0("Taxon_", 1:6),
#'   control_set  = paste0("Taxon_", 20:27)
#' )
#'
#' scores <- ssTaxSEA(counts, custom_db = sets, min_set_size = 3)
#' dim(scores)        # sets x samples
#' round(scores, 2)
#'
#' # The planted signal shows up as higher scores in samples 5-8
#' rowMeans(scores[, 1:4]) - rowMeans(scores[, 5:8])
#'
#' @seealso \code{\link{TaxSEA}} for group-level enrichment.
#' @export
ssTaxSEA <- function(counts,
                     lookup_missing = FALSE,
                     min_set_size = 5,
                     max_set_size = 300,
                     custom_db = NULL,
                     pseudocount = 0.5) {

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
  if (nrow(counts) < 1) stop("'counts' must have at least one taxon.")

  if (!is.numeric(pseudocount) || length(pseudocount) != 1 ||
      is.na(pseudocount) || pseudocount <= 0) {
    stop("'pseudocount' must be a single positive number.")
  }

  if (any(!is.finite(counts))) {
    stop("'counts' contains missing or non-finite values. ",
         "Replace or remove them before running ssTaxSEA.")
  }
  if (any(counts < 0)) {
    stop("'counts' contains negative values. ssTaxSEA expects counts ",
         "or another non-negative abundance measure.")
  }

  # Proportions are not safe here: adding a pseudocount of 0.5 to values
  # that sum to 1 would overwhelm the data entirely.
  col_totals <- colSums(counts)
  if (all(abs(col_totals - 1) < 1e-6)) {
    stop("'counts' appears to contain proportions (columns sum to 1). ",
         "Supply counts instead, or rescale to a count-like scale: ",
         "adding a pseudocount to proportions would dominate the ",
         "transform.")
  }
  if (any(col_totals == 0)) {
    stop("One or more samples have zero total abundance. ",
         "Remove empty samples before running ssTaxSEA.")
  }

  if (any(grepl("\\[|\\]", rownames(counts)))) {
    stop("Taxon names contain square brackets [ ]. ",
         "Please remove or rename these entries before running ssTaxSEA.")
  }

  # --- CLR transform on the FULL table ---
  # This must happen before any subsetting to set members, so that the
  # geometric mean each value is divided by reflects the whole community
  # rather than whichever taxa happen to be in the tested sets.
  clr_mat <- ss_clr(counts, pseudocount = pseudocount)

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
    return(matrix(numeric(0), nrow = 0, ncol = ncol(counts),
                  dimnames = list(character(0), colnames(counts))))
  }

  # --- Map row names to NCBI IDs (if using default DB) ---
  if (is.null(custom_db)) {
    mapped_rows <- rownames(clr_mat) %in% names(id_map)
    clr_mat <- clr_mat[mapped_rows, , drop = FALSE]
    rownames(clr_mat) <- id_map[rownames(clr_mat)]
  }

  taxa_ids <- rownames(clr_mat)
  if (length(taxa_ids) < 3) {
    warning("Very few taxa overlap with taxon sets (", length(taxa_ids),
            "). Results may be unreliable.")
  }

  # --- Score: mean CLR of each set's members, per sample ---
  # rbind over a named list keeps the sets-as-rows orientation for any
  # number of samples, including one. vapply would simplify to a vector
  # in the single-sample case and silently transpose the result.
  score_mat <- do.call(rbind, lapply(taxon_sets, function(set) {
    members <- intersect(set, taxa_ids)
    if (length(members) == 0) {
      return(stats::setNames(rep(NA_real_, ncol(clr_mat)),
                             colnames(clr_mat)))
    }
    colMeans(clr_mat[members, , drop = FALSE])
  }))

  dimnames(score_mat) <- list(names(taxon_sets), colnames(clr_mat))
  score_mat
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
#' Applies the centered log-ratio transformation per sample (column).
#' The pseudocount is added to every value, not only to zeros, so that
#' the transform stays monotonic in the counts.
#'
#' @param mat Numeric matrix (taxa x samples).
#' @param pseudocount Value added to every count before the log.
#' @return CLR-transformed matrix of the same dimensions.
#' @keywords internal
#' @noRd
ss_clr <- function(mat, pseudocount = 0.5) {
  log_mat <- log(mat + pseudocount)
  sweep(log_mat, 2, colMeans(log_mat), FUN = "-")
}
