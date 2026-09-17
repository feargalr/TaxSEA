# skip_if_offline() only checks DNS. BugSigDB is served from Zenodo, which
# can be unreachable while the network is fine, so probe the download itself.
skip_if_no_bugsigdb <- function() {
  skip_if_not_installed("bugsigdbr")
  ok <- tryCatch({
    bugsigdbr::importBugSigDB()
    TRUE
  }, error = function(e) FALSE)
  if (!ok) skip("BugSigDB could not be downloaded")
}

test_that("TaxSEA returns the expected output structure", {
  data("TaxSEA_test_data", package = "TaxSEA")

  # bugsigdb = FALSE keeps this offline and fast. The BugSigDB element is
  # still present in the output, just empty.
  res <- suppressWarnings(
    TaxSEA(TaxSEA_test_data, bugsigdb = FALSE)
  )

  expect_type(res, "list")
  expect_named(res, c("All_databases", "Metabolite_producers",
                      "Health_associations",
                      "BacDive_bacterial_physiology", "BugSigDB",
                      "Gut_Brain_Modules_VallesColomer2019"))

  for (df in res) {
    expect_s3_class(df, "data.frame")
    expect_true(all(c("median_rank_of_set_members", "PValue",
                      "Test_statistic", "FDR") %in% colnames(df)))
  }

  # Every category is a subset of All_databases
  total <- sum(vapply(res[setdiff(names(res), "All_databases")],
                      nrow, integer(1)))
  expect_equal(total, nrow(res$All_databases))
})

test_that("BugSigDB signatures are included when requested", {
  skip_if_no_bugsigdb()
  data("TaxSEA_test_data", package = "TaxSEA")

  res <- suppressWarnings(TaxSEA(TaxSEA_test_data, bugsigdb = TRUE))

  expect_gt(nrow(res$BugSigDB), 0)
  expect_true("PubMedID" %in% colnames(res$BugSigDB))
  expect_true("BugSigDB_ID" %in% colnames(res$BugSigDB))
})

test_that("BugSigDB is excluded by default", {
  data("TaxSEA_test_data", package = "TaxSEA")

  expect_false(formals(TaxSEA)$bugsigdb)

  by_default <- suppressWarnings(TaxSEA(TaxSEA_test_data))
  explicit <- suppressWarnings(TaxSEA(TaxSEA_test_data, bugsigdb = FALSE))

  expect_equal(nrow(by_default$BugSigDB), 0)
  expect_false(any(grepl("bsdb", by_default$All_databases$taxonSetName)))
  expect_equal(by_default, explicit)
})

test_that("including BugSigDB shifts p-values in the other categories", {
  # Documents a real coupling rather than endorsing it: BugSigDB sets
  # enlarge the union of set members, which enlarges the background the
  # competitive KS test runs against, so every other category's p-values
  # move. Results are therefore not comparable between runs made with and
  # without BugSigDB, nor across BugSigDB releases.
  skip_if_no_bugsigdb()
  data("TaxSEA_test_data", package = "TaxSEA")

  with_bsdb <- suppressWarnings(TaxSEA(TaxSEA_test_data, bugsigdb = TRUE))
  without <- suppressWarnings(TaxSEA(TaxSEA_test_data, bugsigdb = FALSE))

  a <- with_bsdb$Metabolite_producers
  b <- without$Metabolite_producers

  # Same sets survive filtering either way ...
  expect_setequal(a$taxonSetName, b$taxonSetName)
  # ... but the p-values differ
  merged <- merge(a[, c("taxonSetName", "PValue")],
                  b[, c("taxonSetName", "PValue")],
                  by = "taxonSetName")
  expect_gt(max(abs(merged$PValue.x - merged$PValue.y)), 0)
})
