## Regression guard for database updates ---------------------------------
#
# These snapshots pin the results TaxSEA produces from the bundled test
# data against the bundled database. They are expected to change when
# TaxSEA_db is updated -- that is the point. A failure here means "the
# database changed the answers", and the diff should be reviewed and
# accepted deliberately with testthat::snapshot_accept() rather than
# regenerated on autopilot.
#
# BugSigDB is downloaded at run time and drifts independently of this
# package, so these tests pass bugsigdb = FALSE to stay hermetic.

skip_if_no_db <- function() {
  skip_if_not_installed("testthat")
}

# Deterministic, network-free summary of one result table.
summarise_category <- function(df, n = 20) {
  if (is.null(df) || nrow(df) == 0) {
    return(data.frame(
      taxonSetName = character(0),
      PValue = character(0),
      FDR = character(0),
      stringsAsFactors = FALSE
    ))
  }
  name_col <- if ("taxonSetName" %in% colnames(df)) {
    "taxonSetName"
  } else {
    colnames(df)[1]
  }
  df <- df[order(df$PValue, df[[name_col]]), , drop = FALSE]
  df <- utils::head(df, n)
  data.frame(
    taxonSetName = as.character(df[[name_col]]),
    PValue = formatC(df$PValue, format = "e", digits = 6),
    FDR = formatC(df$FDR, format = "e", digits = 6),
    stringsAsFactors = FALSE
  )
}

test_that("enrichment results on the bundled test data are stable", {
  skip_if_no_db()
  data("TaxSEA_test_data", package = "TaxSEA")

  res <- suppressWarnings(
    TaxSEA(taxon_ranks = TaxSEA_test_data, bugsigdb = FALSE)
  )

  for (nm in c("Metabolite_producers", "Health_associations",
               "BacDive_bacterial_physiology",
               "Gut_Brain_Modules_VallesColomer2019")) {
    expect_snapshot_value(
      summarise_category(res[[nm]]),
      style = "json2",
      tolerance = 1e-8
    )
  }
})

test_that("set counts per category are stable", {
  skip_if_no_db()
  data("TaxSEA_test_data", package = "TaxSEA")

  res <- suppressWarnings(TaxSEA(taxon_ranks = TaxSEA_test_data, bugsigdb = FALSE))
  counts <- vapply(
    res[c("Metabolite_producers", "Health_associations",
          "BacDive_bacterial_physiology",
          "Gut_Brain_Modules_VallesColomer2019")],
    nrow, integer(1)
  )

  expect_snapshot_value(counts, style = "json2")
})

test_that("the bundled database has the expected shape", {
  data("TaxSEA_db", package = "TaxSEA")

  # Guards against an update that accidentally truncates or explodes the
  # database, or that introduces empty or unnamed sets.
  expect_true(is.list(TaxSEA_db))
  expect_true(all(nzchar(names(TaxSEA_db))))
  expect_true(all(vapply(TaxSEA_db, is.character, logical(1))))
  expect_true(all(lengths(TaxSEA_db) > 0))

  # Known, pre-existing duplicate. The v1.1.8 spelling fix renamed
  # "Phenyalanine" onto the correctly spelled name rather than merging the
  # two sets, leaving two entries of 14 and 2 members. List lookup returns
  # only the first, so those 2 taxa never reach any analysis, and the
  # 2-member entry is in any case below min_set_size. Merge the two sets in
  # TaxSEA_db to resolve; do not extend this list to cover a new duplicate.
  expect_identical(
    unique(names(TaxSEA_db)[duplicated(names(TaxSEA_db))]),
    "GutMGene_producers_of_Phenylalanine"
  )

  expect_snapshot_value(
    c(n_sets = length(TaxSEA_db),
      n_members = length(unlist(TaxSEA_db, use.names = FALSE)),
      n_unique_members = length(unique(unlist(TaxSEA_db,
                                              use.names = FALSE)))),
    style = "json2"
  )
})

# Known, pre-existing gap: 57 members of TaxSEA_db have no entry in
# NCBI_ids and are therefore silently dropped by taxsea_prepare(), which
# shrinks the affected sets. The worst affected are Siderophore_producers
# (73 members -> 33 effective) and the Valles-Colomer2019 Gut-Brain
# Modules (87 member-slots lost across 26 sets). No BacDive member is
# affected: the 1.5.5 rebuild added every BacDive taxid to NCBI_ids, which
# also lowered this from 84.
#
# This ceiling stops the gap growing. Lower it as NCBI_ids is extended;
# do not raise it to accommodate a database update -- extend NCBI_ids
# instead.
KNOWN_UNRESOLVED_MEMBERS <- 57L

unresolved_members <- function() {
  db <- get("TaxSEA_db")
  ids <- get("NCBI_ids")
  setdiff(unique(unlist(db, use.names = FALSE)), unique(unname(ids)))
}

test_that("database members missing from NCBI_ids do not increase", {
  data("TaxSEA_db", package = "TaxSEA")
  data("NCBI_ids", package = "TaxSEA")

  unresolved <- unresolved_members()

  expect_lte(
    length(unresolved), KNOWN_UNRESOLVED_MEMBERS,
    label = paste0(
      length(unresolved), " set member(s) are absent from NCBI_ids ",
      "(known baseline ", KNOWN_UNRESOLVED_MEMBERS, "), e.g. ",
      paste(utils::head(unresolved, 5), collapse = ", "),
      ". These are dropped silently and shrink their sets; extend ",
      "NCBI_ids rather than raising this ceiling."
    )
  )
})

test_that("the per-source breakdown of unresolved members is stable", {
  data("TaxSEA_db", package = "TaxSEA")
  data("NCBI_ids", package = "TaxSEA")

  unresolved <- unresolved_members()
  affected <- vapply(TaxSEA_db, function(s) sum(s %in% unresolved),
                     integer(1))
  by_source <- tapply(affected, sub("_.*", "", names(TaxSEA_db)), sum)

  expect_snapshot_value(by_source[by_source > 0], style = "json2")
})

test_that("no set is emptied below min_set_size by unresolved members", {
  data("TaxSEA_db", package = "TaxSEA")
  data("NCBI_ids", package = "TaxSEA")

  unresolved <- unresolved_members()
  sizes <- lengths(TaxSEA_db)
  lost <- vapply(TaxSEA_db, function(s) sum(s %in% unresolved), integer(1))
  vanishing <- names(TaxSEA_db)[sizes >= 5 & (sizes - lost) < 5]

  # Sets here look usable in TaxSEA_db but never reach the output.
  expect_snapshot_value(sort(vanishing), style = "json2")
})
