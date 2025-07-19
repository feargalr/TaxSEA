test_that("TaxSEA returns valid output with test data", {
  data("TaxSEA_test_data", package = "TaxSEA")
  res <- TaxSEA(TaxSEA_test_data)

  expect_type(res, "list")
  expect_named(res, c("Metabolite_producers", "Health_associations", "BugSigDB"))

  for (df in res) {
    expect_s3_class(df, "data.frame")
    expect_true(all(c("median_rank_of_set_members", "PValue", "Test_statistic", "FDR") %in% colnames(df)))
  }
})
