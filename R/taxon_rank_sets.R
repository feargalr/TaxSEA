#' Taxon Rank Set Enrichment Analysis
#'
#' Groups species by taxonomic ranks and performs TaxSEA enrichment
#' analysis at each rank level. Returns a named list of data frames,
#' one per taxonomic rank (excluding species).
#'
#' @param taxon_ranks A named numeric vector of log2 fold changes.
#'   Names should be feature identifiers matching the
#'   \code{species} column in \code{lineage_df}, or matching
#'   \code{rownames()} of a SummarizedExperiment/TreeSummarizedExperiment.
#' @param lineage_df Either a data frame or a
#'   \code{SummarizedExperiment}/\code{TreeSummarizedExperiment} object.
#'
#'   \strong{Data frame input:} Must include a \code{species} column
#'   and one or more taxonomic rank columns (e.g., kingdom, phylum,
#'   class, order, family, genus). The \code{species} column is used
#'   for matching and is excluded from the enrichment analysis.
#'
#'   \strong{SummarizedExperiment input:} Taxonomy is extracted from
#'   \code{rowData()} and feature identifiers from \code{rownames()}.
#'   If the \code{mia} package is installed, \code{taxonomyRanks()} is
#'   used to identify taxonomic rank columns; otherwise all
#'   \code{rowData()} columns are used. The Species rank column is
#'   automatically excluded from the analysis. Requires the
#'   \code{SummarizedExperiment} package.
#' @param min_set_size Minimum number of species in a set to include
#'   in the analysis. Default is 5.
#' @param max_set_size Maximum number of species in a set to include
#'   in the analysis. Default is 100.
#' @return A named list of data frames, one per taxonomic rank. Each
#'   data frame contains columns: taxonSetName,
#'   median_rank_of_set_members, PValue, Test_statistic, and FDR.
#'   Ranks that produce no valid sets (e.g., due to size filtering)
#'   are included as empty data frames with a message.
#' @examples
#' # --- Example 1: Data frame input ---
#' # Create a lineage data frame (e.g., parsed from curatedMetagenomicData)
#' # The 'species' column must match the names in taxon_ranks.
#'
#' lineage_df <- data.frame(
#'   species = c("Cutibacterium_acnes", "Klebsiella_pneumoniae",
#'               "Propionibacterium_humerusii", "Moraxella_osloensis",
#'               "Enhydrobacter_aerosaccus", "Staphylococcus_capitis",
#'               "Staphylococcus_epidermidis", "Staphylococcus_aureus",
#'               "Escherichia_coli", "Enterobacter_cloacae",
#'               "Pseudomonas_aeruginosa", "Acinetobacter_baumannii",
#'               "Lactobacillus_rhamnosus", "Lactobacillus_acidophilus",
#'               "Bifidobacterium_longum", "Bifidobacterium_breve"),
#'   kingdom = rep("Bacteria", 16),
#'   phylum = c("Actinobacteria", "Proteobacteria",
#'              "Actinobacteria", "Proteobacteria",
#'              "Proteobacteria", "Firmicutes",
#'              "Firmicutes", "Firmicutes",
#'              "Proteobacteria", "Proteobacteria",
#'              "Proteobacteria", "Proteobacteria",
#'              "Firmicutes", "Firmicutes",
#'              "Actinobacteria", "Actinobacteria"),
#'   class = c("Actinobacteria", "Gammaproteobacteria",
#'             "Actinobacteria", "Gammaproteobacteria",
#'             "Alphaproteobacteria", "Bacilli",
#'             "Bacilli", "Bacilli",
#'             "Gammaproteobacteria", "Gammaproteobacteria",
#'             "Gammaproteobacteria", "Gammaproteobacteria",
#'             "Bacilli", "Bacilli",
#'             "Actinobacteria", "Actinobacteria"),
#'   order = c("Propionibacteriales", "Enterobacterales",
#'             "Propionibacteriales", "Pseudomonadales",
#'             "Rhodospirillales", "Bacillales",
#'             "Bacillales", "Bacillales",
#'             "Enterobacterales", "Enterobacterales",
#'             "Pseudomonadales", "Pseudomonadales",
#'             "Lactobacillales", "Lactobacillales",
#'             "Bifidobacteriales", "Bifidobacteriales"),
#'   family = c("Propionibacteriaceae", "Enterobacteriaceae",
#'              "Propionibacteriaceae", "Moraxellaceae",
#'              "Rhodospirillaceae", "Staphylococcaceae",
#'              "Staphylococcaceae", "Staphylococcaceae",
#'              "Enterobacteriaceae", "Enterobacteriaceae",
#'              "Pseudomonadaceae", "Moraxellaceae",
#'              "Lactobacillaceae", "Lactobacillaceae",
#'              "Bifidobacteriaceae", "Bifidobacteriaceae"),
#'   genus = c("Cutibacterium", "Klebsiella",
#'             "Cutibacterium", "Moraxella",
#'             "Enhydrobacter", "Staphylococcus",
#'             "Staphylococcus", "Staphylococcus",
#'             "Escherichia", "Enterobacter",
#'             "Pseudomonas", "Acinetobacter",
#'             "Lactobacillus", "Lactobacillus",
#'             "Bifidobacterium", "Bifidobacterium"),
#'   stringsAsFactors = FALSE
#' )
#'
#' set.seed(42)
#' fc <- setNames(rnorm(16), lineage_df$species)
#' results <- taxon_rank_sets(fc, lineage_df, min_set_size = 2)
#' names(results)
#' results$family
#'
#' # --- Example 2: SummarizedExperiment / TreeSummarizedExperiment input ---
#' \dontrun{
#' library(mia)
#' data(GlobalPatterns, package = "mia")
#' tse <- GlobalPatterns
#'
#' # Run differential abundance (e.g., ALDEx2) to get fold changes
#' # aldex_out <- ... (your DA analysis)
#' # fc <- aldex_out$effect
#' # names(fc) <- rownames(aldex_out)
#'
#' # Run taxon rank set enrichment directly from the TSE
#' results <- taxon_rank_sets(fc, tse, min_set_size = 5)
#' names(results)  # Kingdom, Phylum, Class, Order, Family, Genus
#' results$Family
#' }
#'
#' @export
taxon_rank_sets <- function(taxon_ranks, lineage_df,
                            min_set_size = 5, max_set_size = 100) {

  if (!is.numeric(taxon_ranks) || is.null(names(taxon_ranks))) {
    stop("'taxon_ranks' must be a named numeric vector.")
  }

  # --- Handle SummarizedExperiment / TreeSummarizedExperiment input ---
  if (is(lineage_df, "SummarizedExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("The 'SummarizedExperiment' package is required when ",
           "providing a SummarizedExperiment or TreeSummarizedExperiment ",
           "object. Install it with: BiocManager::install(",
           "'SummarizedExperiment')")
    }

    message("Detected SummarizedExperiment input. Extracting taxonomy ",
            "from rowData().")

    # Extract taxonomy from rowData
    rd <- as.data.frame(SummarizedExperiment::rowData(lineage_df))

    if (ncol(rd) == 0) {
      stop("rowData() of the SummarizedExperiment is empty. ",
           "Taxonomic rank columns are required.")
    }

    # Determine which columns are taxonomic ranks
    if (requireNamespace("mia", quietly = TRUE)) {
      rank_cols <- mia::taxonomyRanks(lineage_df)
      message("Using mia::taxonomyRanks() to identify rank columns: ",
              paste(rank_cols, collapse = ", "))
    } else {
      rank_cols <- colnames(rd)
      message("mia package not installed. Using all rowData columns ",
              "as taxonomic ranks: ", paste(rank_cols, collapse = ", "))
    }

    # Use rownames as feature identifiers
    feature_ids <- rownames(SummarizedExperiment::rowData(lineage_df))
    rd$species <- feature_ids

    # Exclude Species column from rank columns (case-insensitive)
    rank_cols <- rank_cols[tolower(rank_cols) != "species"]

    if (length(rank_cols) == 0) {
      stop("No taxonomic rank columns found after excluding Species.")
    }

    # Convert to data.frame format expected by downstream code
    lineage_df <- rd[, c("species", rank_cols), drop = FALSE]
  }

  # --- Validate data frame input ---
  if (!is.data.frame(lineage_df)) {
    stop("'lineage_df' must be a data frame or a ",
         "SummarizedExperiment/TreeSummarizedExperiment object.")
  }

  if (!"species" %in% colnames(lineage_df)) {
    stop("'lineage_df' must contain a 'species' column.")
  }

  # Identify taxonomic rank columns (all columns except 'species')
  rank_cols <- setdiff(colnames(lineage_df), "species")

  if (length(rank_cols) == 0) {
    stop("'lineage_df' must contain at least one taxonomic rank column ",
         "in addition to 'species'.")
  }

  # Filter lineage_df to species present in taxon_ranks
  shared_species <- intersect(names(taxon_ranks), lineage_df$species)

  if (length(shared_species) == 0) {
    stop("No overlap between names in 'taxon_ranks' and ",
         "'lineage_df$species'. Ensure that the names of your ",
         "taxon_ranks vector match the feature identifiers in ",
         "your lineage data.")
  }

  if (length(shared_species) < length(names(taxon_ranks))) {
    n_missing <- length(names(taxon_ranks)) - length(shared_species)
    message(n_missing, " species in 'taxon_ranks' not found in ",
            "'lineage_df' and will be excluded.")
  }

  lineage_sub <- lineage_df[lineage_df$species %in% shared_species, ]
  taxon_ranks_sub <- taxon_ranks[shared_species]

  # Run TaxSEA for each taxonomic rank
  results <- list()

  for (rank in rank_cols) {
    # Build custom taxon sets: named list of species vectors grouped by rank
    rank_groups <- split(lineage_sub$species, lineage_sub[[rank]])

    # Remove groups with NA or empty rank names
    rank_groups <- rank_groups[!is.na(names(rank_groups)) &
                                 names(rank_groups) != ""]

    if (length(rank_groups) == 0) {
      message("No valid groups found for rank '", rank,
              "'. Returning empty data frame.")
      results[[rank]] <- data.frame(
        taxonSetName = character(0),
        median_rank_of_set_members = numeric(0),
        PValue = numeric(0),
        Test_statistic = numeric(0),
        FDR = numeric(0),
        stringsAsFactors = FALSE
      )
      next
    }

    res <- tryCatch({
      taxsea_res <- TaxSEA(taxon_ranks = taxon_ranks_sub,
                            custom_db = rank_groups,
                            min_set_size = min_set_size,
                            max_set_size = max_set_size)
      taxsea_res$custom_sets
    }, error = function(e) {
      message("TaxSEA returned no results for rank '", rank,
              "': ", conditionMessage(e))
      data.frame(
        taxonSetName = character(0),
        median_rank_of_set_members = numeric(0),
        PValue = numeric(0),
        Test_statistic = numeric(0),
        FDR = numeric(0),
        stringsAsFactors = FALSE
      )
    })

    results[[rank]] <- res
  }

  return(results)
}
