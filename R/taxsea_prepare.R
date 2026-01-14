#' Prepare TaxSEA inputs
#'
#' Internal helper that loads databases, maps taxa IDs, and filters sets/taxa.
#' Supports enrichment (ranked vector) or ORA (input taxa vector).
#'
#' @param taxon_ranks Named numeric vector (e.g., log2 fold changes). For enrichment.
#' @param input_taxa Character vector of taxa (hits). For ORA.
#' @param lookup_missing Logical; fetch missing NCBI IDs (default FALSE).
#' @param min_set_size Minimum set size to include.
#' @param max_set_size Maximum set size to include.
#' @param custom_db Optional list of taxon sets; if NULL uses built-in DB.
#'
#' @return An object of class \code{"TaxSEA_prep"}.
#' @keywords internal
#' @noRd
taxsea_prepare <- function(taxon_ranks = NULL,
                           input_taxa = NULL,
                           lookup_missing = FALSE,
                           min_set_size = 5,
                           max_set_size = 100,
                           custom_db = NULL) {
  
  has_ranks <- !is.null(taxon_ranks)
  has_taxa  <- !is.null(input_taxa)
  
  if (has_ranks && has_taxa) {
    stop("Provide either 'taxon_ranks' (enrichment) or 'input_taxa' (ORA), but not both.")
  }
  if (!has_ranks && !has_taxa) {
    stop("Provide either 'taxon_ranks' (enrichment) or 'input_taxa' (ORA).")
  }
  
  mode <- if (has_ranks) "enrichment" else "ora"
  
  # Shared bracket check
  if (mode == "enrichment") {
    if (any(grepl("\\[|\\]", names(taxon_ranks)))) {
      stop("Taxon names contain square brackets [ ]. 
           Please remove or rename these entries before running TaxSEA.")
    }
  } else {
    if (!is.character(input_taxa)) stop("'input_taxa' must be a character vector.")
    input_taxa <- unique(input_taxa)
    if (any(grepl("\\[|\\]", input_taxa))) {
      stop("Input taxa contain square brackets [ ]. 
           Please remove or rename these entries before running TaxSEA.")
    }
  }
  
  # ---- Load DB (same logic as your TaxSEA()) ----
  if (is.null(custom_db)) {
    utils::data("TaxSEA_db", package = "TaxSEA", envir = environment())
    taxon_sets <- TaxSEA_db
    
    if (requireNamespace("bugsigdbr", quietly = TRUE)) {
      bsdb <- bugsigdbr::importBugSigDB()
      mp.sigs <- bugsigdbr::getSignatures(bsdb, tax.id.type = "ncbi")
      bugsigdb_list <- utils::stack(mp.sigs)
      names(bugsigdb_list) <- c("Species", "MSID")
      bugsigdb_list <- split(bugsigdb_list$Species, bugsigdb_list$MSID)
      names(bugsigdb_list) <- paste0("bsdb_", names(bugsigdb_list))
      taxon_sets <- c(taxon_sets, bugsigdb_list)
    } else {
      warning("bugsigdbr not installed; skipping BugSigDB integration.")
    }
  } else {
    taxon_sets <- custom_db
    if (!is.list(taxon_sets)) stop("Custom database must be a list of taxon sets.")
  }
  
  # ---- ENRICHMENT PATH ----
  if (mode == "enrichment") {
    
    if (is.null(custom_db)) {
      utils::data("NCBI_ids", package = "TaxSEA", envir = environment())
      
      if (lookup_missing) {
        ids2fetch <- names(taxon_ranks[!(names(taxon_ranks) %in% names(NCBI_ids))])
        if (length(ids2fetch) > 0) {
          fetched_ids <- get_ncbi_taxon_ids(ids2fetch)
          if (length(unlist(fetched_ids)) > 0) {
            NCBI_ids <- c(NCBI_ids, unlist(fetched_ids))
          }
        }
      }
      
      taxon_ranks <- taxon_ranks[names(taxon_ranks) %in% names(NCBI_ids)]
      original_ranks <- taxon_ranks
      names(taxon_ranks) <- NCBI_ids[names(taxon_ranks)]
      NCBI_ids <- NCBI_ids[names(NCBI_ids) %in% names(original_ranks)]
      NCBI_ids <- NCBI_ids[!duplicated(NCBI_ids)]
      legible_names <- names(NCBI_ids)
      names(legible_names) <- NCBI_ids[legible_names]
    } else {
      legible_names <- names(taxon_ranks)
    }
    
    if (length(taxon_ranks) < 3) {
      warning("Error: Very few taxa provided; results may be unreliable.")
    }
    
    taxon_ranks <- taxon_ranks[names(taxon_ranks) %in% unique(unlist(taxon_sets))]
    taxon_sets <- lapply(taxon_sets, function(set) unique(set[set %in% names(taxon_ranks)]))
    set_sizes <- vapply(taxon_sets, length, numeric(1))
    taxon_sets <- taxon_sets[set_sizes >= min_set_size & set_sizes <= max_set_size]
    taxon_ranks <- taxon_ranks[names(taxon_ranks) %in% unique(unlist(taxon_sets))]
    taxon_sets <- lapply(taxon_sets, function(set) intersect(names(taxon_ranks), set))
    
    background_taxa <- names(taxon_ranks)
    
    return(structure(
      list(
        mode = "enrichment",
        taxon_ranks = taxon_ranks,
        input_taxa = NULL,
        background_taxa = background_taxa,
        taxon_sets = taxon_sets,
        legible_names = legible_names,
        custom_db = custom_db,
        used_default_db = is.null(custom_db),
        min_set_size = min_set_size,
        max_set_size = max_set_size,
        lookup_missing = lookup_missing
      ),
      class = "TaxSEA_prep"
    ))
  }
  
  # ---- ORA PATH ----
  if (is.null(custom_db)) {
    utils::data("NCBI_ids", package = "TaxSEA", envir = environment())
    
    if (lookup_missing) {
      ids2fetch <- input_taxa[!(input_taxa %in% names(NCBI_ids))]
      if (length(ids2fetch) > 0) {
        fetched_ids <- get_ncbi_taxon_ids(ids2fetch)
        if (length(unlist(fetched_ids)) > 0) {
          NCBI_ids <- c(NCBI_ids, unlist(fetched_ids))
        }
      }
    }
    NCBI_ids <- NCBI_ids[!duplicated(NCBI_ids)]
    legible_names <- names(NCBI_ids)
    names(legible_names) <- NCBI_ids[legible_names]
    input_taxa <- input_taxa[input_taxa %in% names(NCBI_ids)]
    input_taxa <- unname(NCBI_ids[input_taxa])

  } else {
    # custom db: no mapping
    all_taxa <- unique(unlist(taxon_sets, use.names = FALSE))
    legible_names <- all_taxa
    names(legible_names) <- all_taxa
  }
  
  # Filter sets by size *as stored* (no rank shrinking)
  taxon_sets <- lapply(taxon_sets, unique)
  set_sizes <- vapply(taxon_sets, length, numeric(1))
  taxon_sets <- taxon_sets[set_sizes >= min_set_size & set_sizes <= max_set_size]
  
  background_taxa <- unique(unlist(taxon_sets, use.names = FALSE))
  
  # Keep only hits in background
  input_taxa <- unique(input_taxa[input_taxa %in% background_taxa])
  
  # Intersect sets with background (mostly redundant)
  taxon_sets <- lapply(taxon_sets, function(set) intersect(background_taxa, set))
  taxon_sets <- taxon_sets[vapply(taxon_sets, length, numeric(1)) > 0]
  
  structure(
    list(
      mode = "ora",
      taxon_ranks = NULL,
      input_taxa = input_taxa,
      background_taxa = background_taxa,
      taxon_sets = taxon_sets,
      legible_names = legible_names,
      custom_db = custom_db,
      used_default_db = is.null(custom_db),
      min_set_size = min_set_size,
      max_set_size = max_set_size,
      lookup_missing = lookup_missing
    ),
    class = "TaxSEA_prep"
  )
}