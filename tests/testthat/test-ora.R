## Over-representation analysis (Fisher's exact test) --------------------

# A custom database keeps these tests offline and independent of the
# built-in database's contents.
ora_db <- function() {
  list(
    hit_rich = paste0("tx", 1:10),   # 8 of 10 are hits
    hit_poor = paste0("tx", 21:30),  # 0 of 10 are hits
    mixed    = paste0("tx", 5:14)    # spans both
  )
}

ora_hits <- function() paste0("tx", 1:8)

test_that("ORA mode is inferred from input_taxa alone", {
  res <- TaxSEA(input_taxa = ora_hits(), custom_db = ora_db())

  expect_type(res, "list")
  expect_named(res, "custom_sets")
  expect_s3_class(res$custom_sets, "data.frame")
})

test_that("ORA p-value and odds ratio match a hand-built Fisher test", {
  db <- ora_db()
  hits <- ora_hits()
  res <- TaxSEA(input_taxa = hits, custom_db = db, mode = "ora")$custom_sets

  universe <- unique(unlist(db))
  in_universe <- hits[hits %in% universe]

  set_members <- db$hit_rich
  a <- sum(set_members %in% in_universe)          # hits in set
  b <- length(set_members) - a                    # non-hits in set
  cc <- length(in_universe) - a                   # hits outside set
  d <- length(universe) - length(in_universe) - b # non-hits outside set

  expected <- stats::fisher.test(
    matrix(c(a, b, cc, d), nrow = 2, byrow = TRUE),
    alternative = "greater"
  )

  row <- res[res$taxonSetName == "hit_rich", ]
  expect_equal(row$PValue, expected$p.value)
  expect_equal(row$Test_statistic, unname(expected$estimate))
})

test_that("a set with no hits is less significant than a set full of them", {
  res <- TaxSEA(input_taxa = ora_hits(), custom_db = ora_db(),
                mode = "ora")$custom_sets

  p_rich <- res$PValue[res$taxonSetName == "hit_rich"]
  p_poor <- res$PValue[res$taxonSetName == "hit_poor"]

  expect_lt(p_rich, 0.05)
  expect_gt(p_poor, p_rich)
})

test_that("ORA results are sorted by p-value and carry an FDR column", {
  res <- TaxSEA(input_taxa = ora_hits(), custom_db = ora_db(),
                mode = "ora")$custom_sets

  expect_false(is.unsorted(res$PValue))
  expect_true(all(c("PValue", "FDR") %in% colnames(res)))
  expect_true(all(res$FDR >= res$PValue - 1e-12))
})

test_that("ORA against the built-in database splits results by category", {
  data("TaxSEA_test_data", package = "TaxSEA")
  hits <- names(TaxSEA_test_data)[TaxSEA_test_data > 0]

  res <- suppressWarnings(TaxSEA(input_taxa = hits, mode = "ora", bugsigdb = FALSE))

  expect_true(all(c("All_databases", "Metabolite_producers",
                    "Health_associations") %in% names(res)))
  # ORA reports an odds ratio, and has no ranks to take a median of
  expect_true("Odds_ratio" %in% colnames(res$All_databases))
  expect_false("median_rank_of_set_members" %in%
                 colnames(res$All_databases))
})

## Mode selection and input validation -----------------------------------

test_that("TaxSEA rejects ambiguous or missing input", {
  data("TaxSEA_test_data", package = "TaxSEA")

  expect_error(
    TaxSEA(taxon_ranks = TaxSEA_test_data, input_taxa = c("a", "b")),
    "not both"
  )
  expect_error(TaxSEA(), "Provide either")
})

test_that("mode must agree with the input that was supplied", {
  data("TaxSEA_test_data", package = "TaxSEA")

  expect_error(
    TaxSEA(taxon_ranks = TaxSEA_test_data, mode = "ora"),
    "requires"
  )
  expect_error(
    TaxSEA(input_taxa = names(TaxSEA_test_data), mode = "enrichment"),
    "requires"
  )
})

test_that("square-bracketed taxon names are rejected in both modes", {
  ranks <- c(`[Ruminococcus] gnavus` = 1.5, Escherichia_coli = -0.5)

  expect_error(TaxSEA(taxon_ranks = ranks), "square brackets")
  expect_error(TaxSEA(input_taxa = names(ranks), mode = "ora"),
               "square brackets")
})

test_that("input_taxa must be a character vector", {
  expect_error(TaxSEA(input_taxa = c(1, 2, 3), mode = "ora"),
               "character vector")
})

## Enrichment core -------------------------------------------------------

test_that("KS statistic and p-value match a direct ks.test call", {
  set.seed(11)
  taxa <- paste0("tx", seq_len(60))
  ranks <- stats::setNames(rnorm(60), taxa)
  db <- list(set_a = taxa[1:12], set_b = taxa[20:34])

  res <- TaxSEA(taxon_ranks = ranks, custom_db = db,
                min_set_size = 5)$custom_sets

  # The test is competitive against the taxa covered by the surviving sets,
  # not against every taxon supplied: taxsea_prepare() drops ranks for taxa
  # that appear in no set before the KS test runs.
  background <- ranks[names(ranks) %in% unique(unlist(db))]

  for (nm in names(db)) {
    expected <- suppressWarnings(
      stats::ks.test(ranks[db[[nm]]], background)
    )
    row <- res[res$taxonSetName == nm, ]
    expect_equal(row$PValue, expected$p.value)
    expect_equal(row$Test_statistic, unname(expected$statistic))
    expect_equal(row$median_rank_of_set_members,
                 stats::median(ranks[db[[nm]]]))
  }
})

test_that("taxa in no taxon set are excluded from the background", {
  set.seed(13)
  taxa <- paste0("tx", seq_len(40))
  ranks <- stats::setNames(rnorm(40), taxa)
  db <- list(set_a = taxa[1:10], set_b = taxa[11:20])

  covered <- TaxSEA(taxon_ranks = ranks, custom_db = db,
                    min_set_size = 5)$custom_sets

  # Adding taxa that belong to no set must not change any result
  extra <- c(ranks, stats::setNames(rnorm(20), paste0("extra", 1:20)))
  with_extra <- TaxSEA(taxon_ranks = extra, custom_db = db,
                       min_set_size = 5)$custom_sets

  expect_equal(covered$PValue, with_extra$PValue)
  expect_equal(covered$Test_statistic, with_extra$Test_statistic)
})

test_that("custom_db results are returned as a five-column table", {
  set.seed(12)
  taxa <- paste0("tx", seq_len(40))
  ranks <- stats::setNames(rnorm(40), taxa)
  res <- TaxSEA(taxon_ranks = ranks,
                custom_db = list(s1 = taxa[1:10], s2 = taxa[11:25]),
                min_set_size = 5)

  expect_named(res, "custom_sets")
  expect_identical(ncol(res$custom_sets), 5L)
  expect_identical(
    colnames(res$custom_sets),
    c("taxonSetName", "median_rank_of_set_members", "PValue",
      "Test_statistic", "FDR")
  )
})

## FDR handling ----------------------------------------------------------

test_that("FDR is recomputed within each output category", {
  data("TaxSEA_test_data", package = "TaxSEA")
  res <- suppressWarnings(TaxSEA(taxon_ranks = TaxSEA_test_data, bugsigdb = FALSE))

  for (nm in c("Metabolite_producers", "Health_associations",
               "BacDive_bacterial_physiology")) {
    df <- res[[nm]]
    if (nrow(df) == 0) next
    # Each category is corrected over its own rows, not over All_databases
    expect_equal(df$FDR, stats::p.adjust(df$PValue, method = "fdr"),
                 info = nm)
  }
})
