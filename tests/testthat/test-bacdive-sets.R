## BacDive sets ------------------------------------------------------------
#
# Invariants and biological spot checks for the BacDive_* sets built by
# bacdive_harvester (see inst/scripts/make_BacDive_sets.R). Taxa are NCBI
# taxonomy IDs.

bacdive_sets <- function() {
  data("TaxSEA_db", package = "TaxSEA", envir = environment())
  TaxSEA_db[startsWith(names(TaxSEA_db), "BacDive_")]
}

test_that("BacDive set names are well formed and unique", {
  sets <- bacdive_sets()
  expect_gt(length(sets), 1000)
  expect_false(any(grepl("[[:space:]]", names(sets))))
  expect_false(anyDuplicated(names(sets)) > 0)
  # TaxSEA intersects sets with observed taxa before size filtering, so a
  # set smaller than this could never be tested.
  expect_gte(min(lengths(sets)), 3L)
  expect_false(any(vapply(sets, anyDuplicated, integer(1)) > 0))
})

test_that("every BacDive set member can be looked up", {
  data("NCBI_ids", package = "TaxSEA")
  members <- unique(unlist(bacdive_sets(), use.names = FALSE))
  # Present as a name, so users supplying taxids reach it...
  expect_true(all(members %in% names(NCBI_ids)))
  # ...and as a value, so it survives the name -> id translation.
  expect_true(all(members %in% unlist(NCBI_ids, use.names = FALSE)))
})

test_that("BacDive provenance ships and matches the sets", {
  f <- system.file("extdata", "BacDive_set_provenance.tsv", package = "TaxSEA")
  expect_true(nzchar(f))
  prov <- utils::read.delim(f, stringsAsFactors = FALSE)
  sets <- bacdive_sets()
  expect_setequal(prov$set_name, names(sets))
  expect_equal(prov$n_taxa, unname(lengths(sets)[prov$set_name]))

  info <- utils::read.delim(
    system.file("extdata", "BacDive_build_info.tsv", package = "TaxSEA"),
    stringsAsFactors = FALSE)
  expect_true(all(c("harvester_commit", "harvest_last_date", "sets_shipped")
                  %in% info$key))
  expect_equal(as.integer(info$value[info$key == "sets_shipped"]),
               length(sets))
})

test_that("oxygen sets agree with textbook physiology", {
  sets <- bacdive_sets()
  has <- function(set, taxid) taxid %in% sets[[set]]

  # Escherichia coli: facultative anaerobe, reduces nitrate.
  expect_true(has("BacDive_Oxygen_facultative_anaerobe", "562"))
  expect_false(has("BacDive_Oxygen_aerobe", "562"))
  expect_true(has("BacDive_Uses_nitrate_reduction", "562"))
  # Bacteroides uniformis: anaerobe, not microaerophile.
  expect_true(has("BacDive_Oxygen_anaerobe", "820"))
  expect_false(has("BacDive_Oxygen_microaerophile", "820"))
  # Helicobacter pylori: microaerophile.
  expect_true(has("BacDive_Oxygen_microaerophile", "210"))
  # Reclassified taxa, reached only through their current genus:
  # Agathobacter rectalis and Phocaeicola vulgatus.
  expect_true(has("BacDive_Oxygen_anaerobe", "39491"))
  expect_true(has("BacDive_Oxygen_anaerobe", "821"))
  # Faecalibacterium prausnitzii has no oxygen record in BacDive. Earlier
  # releases listed it as an anaerobe from an LLM-derived table; the sets
  # now contain measured data only.
  expect_false(has("BacDive_Oxygen_anaerobe", "853"))
})

test_that("current and former species names map to the same taxid", {
  data("NCBI_ids", package = "TaxSEA")
  id <- function(n) as.character(unlist(NCBI_ids[[n]]))
  expect_identical(id("Phocaeicola vulgatus"), "821")
  expect_identical(id("Bacteroides vulgatus"), "821")
  expect_identical(id("Phocaeicola_vulgatus"), "821")
  expect_identical(id("Agathobacter rectalis"), "39491")
  expect_identical(id("Eubacterium rectale"), "39491")
  expect_identical(id("821"), "821")
})
