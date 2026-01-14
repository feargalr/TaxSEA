# Internal helper: prepare TaxSEA inputs (db + id mapping + filtering)
# Not exported; used to prototype new testing modes (e.g., ORA).
taxsea_prepare <- function(taxon_ranks,
                           lookup_missing = FALSE,
                           min_set_size = 5,
                           max_set_size = 100,
                           custom_db = NULL) {
  
  if (any(grepl("\\[|\\]", names(taxon_ranks)))) {
    stop("Taxon names contain square brackets [ ]. 
         Please remove or rename these entries before running TaxSEA.")
  }
  
  ## 1) Load database (built-in unless custom_db is provided)
  if (is.null(custom_db)) {
    data("TaxSEA_db", package = "TaxSEA", envir = environment())
    taxon_sets <- TaxSEA_db
    
    # Append BugSigDB dynamically (same as TaxSEA())
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
  
  ## 2) Map taxon names -> NCBI IDs (only for default db)
  if (is.null(custom_db)) {
    data("NCBI_ids", package = "TaxSEA", envir = environment())
    
    if (lookup_missing) {
      ids2fetch <- names(taxon_ranks[!(names(taxon_ranks) %in% names(NCBI_ids))])
      if (length(ids2fetch) > 0) {
        fetched_ids <- get_ncbi_taxon_ids(ids2fetch)
        if (length(unlist(fetched_ids)) > 0) {
          NCBI_ids <- c(NCBI_ids, unlist(fetched_ids))
        }
      }
    }
    
    # Convert taxon names to NCBI IDs (exactly as in TaxSEA())
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
  
  ## 3) Sanity warning (same behaviour)
  if (length(taxon_ranks) < 3) {
    warning("Error: Very few taxa provided; results may be unreliable.")
  }
  
  ## 4) Filter taxon_ranks and taxon_sets (exactly as in TaxSEA())
  taxon_ranks <- taxon_ranks[names(taxon_ranks) %in% unique(unlist(taxon_sets))]
  
  taxon_sets <- lapply(taxon_sets, function(set) {
    unique(set[set %in% names(taxon_ranks)])
  })
  
  set_sizes <- vapply(taxon_sets, length, numeric(1))
  taxon_sets <- taxon_sets[set_sizes >= min_set_size & set_sizes <= max_set_size]
  
  taxon_ranks <- taxon_ranks[names(taxon_ranks) %in% unique(unlist(taxon_sets))]
  
  taxon_sets <- lapply(taxon_sets, function(set) {
    intersect(names(taxon_ranks), set)
  })
  
  ## Return a standardized prep object for downstream testing
  structure(
    list(
      taxon_ranks = taxon_ranks,     # named numeric (IDs if default db)
      taxon_sets  = taxon_sets,      # list of IDs (or original names if custom_db)
      legible_names = legible_names, # id -> readable label (default db) OR original names
      custom_db = custom_db,
      used_default_db = is.null(custom_db),
      min_set_size = min_set_size,
      max_set_size = max_set_size,
      lookup_missing = lookup_missing
    ),
    class = "TaxSEA_prep"
  )
}