## Fixtures -------------------------------------------------------------

# Small deterministic matrix, used where an exact value is checked.
toy_counts <- function() {
  m <- matrix(c(10, 20, 30, 40,
                 1,  2,  3,  4,
                 5,  5,  5,  5,
               100, 50, 25, 10),
              nrow = 4, byrow = TRUE)
  rownames(m) <- paste0("t", seq_len(4))
  colnames(m) <- paste0("s", seq_len(4))
  m
}

toy_sets <- list(A = c("t1", "t2"), B = c("t3", "t4"))

# Cohort with a planted signal in a set of rare taxa, plus a null set of
# abundant taxa. Used for the CLR-ordering regression test.
planted_cohort <- function() {
  set.seed(1)
  n_taxa <- 300
  n_samples <- 40
  base <- exp(rnorm(n_taxa, 4, 2))
  m <- vapply(seq_len(n_samples), function(j) rpois(n_taxa, base),
              numeric(n_taxa))
  rownames(m) <- paste0("t", seq_len(n_taxa))
  colnames(m) <- paste0("s", seq_len(n_samples))
  group <- rep(c("ctrl", "dis"), each = 20)
  signal <- rownames(m)[order(base)][seq_len(15)]
  null_set <- rownames(m)[order(-base)][seq_len(15)]
  m[signal, group == "dis"] <- round(m[signal, group == "dis"] * 4)
  list(counts = m, group = group,
       sets = list(signal = signal, null_set = null_set))
}

## Core statistic --------------------------------------------------------

test_that("score is the mean CLR of set members", {
  m <- toy_counts()
  res <- ssTaxSEA(m, custom_db = toy_sets, min_set_size = 2)

  log_m <- log(m + 0.5)
  clr <- sweep(log_m, 2, colMeans(log_m), FUN = "-")
  expected <- rbind(
    A = colMeans(clr[c("t1", "t2"), ]),
    B = colMeans(clr[c("t3", "t4"), ])
  )

  expect_equal(res, expected)
})

test_that("result is a matrix with sets as rows and samples as columns", {
  res <- ssTaxSEA(toy_counts(), custom_db = toy_sets, min_set_size = 2)

  expect_true(is.matrix(res))
  expect_type(res, "double")
  expect_identical(dim(res), c(2L, 4L))
  expect_identical(rownames(res), c("A", "B"))
  expect_identical(colnames(res), paste0("s", seq_len(4)))
  # Not a list: the pre-1.5.2 API returned list(scores, pvalues)
  expect_false(is.list(res))
})

test_that("pseudocount is applied to all values and changes the scores", {
  m <- toy_counts()
  default <- ssTaxSEA(m, custom_db = toy_sets, min_set_size = 2)
  bigger <- ssTaxSEA(m, custom_db = toy_sets, min_set_size = 2,
                     pseudocount = 1)

  expect_false(isTRUE(all.equal(default, bigger)))

  log_m <- log(m + 1)
  clr <- sweep(log_m, 2, colMeans(log_m), FUN = "-")
  expect_equal(bigger["A", ], colMeans(clr[c("t1", "t2"), ]))
})

## The CLR must be computed before subsetting -----------------------------

test_that("CLR uses the full taxon table, not just set members", {
  # Subsetting to set members before transforming makes the geometric mean
  # depend on the sets being tested, which drives a null set's scores in the
  # mirror image of a set that does have signal. Guard against a regression.
  fx <- planted_cohort()
  res <- ssTaxSEA(fx$counts, custom_db = fx$sets, min_set_size = 5)

  diff_of <- function(set) {
    x <- res[set, ]
    mean(x[fx$group == "dis"]) - mean(x[fx$group == "ctrl"])
  }

  # Planted signal is recovered
  expect_gt(diff_of("signal"), 0.5)
  # Null set stays put. With the subset-first ordering it moved to about
  # -0.435, a near-perfect mirror of the signal set.
  expect_lt(abs(diff_of("null_set")), 0.1)
})

test_that("scores for one set do not depend on which other sets are tested", {
  fx <- planted_cohort()
  both <- ssTaxSEA(fx$counts, custom_db = fx$sets, min_set_size = 5)
  alone <- ssTaxSEA(fx$counts, custom_db = fx$sets["null_set"],
                    min_set_size = 5)

  expect_equal(both["null_set", ], alone["null_set", ])
})

## Single sample ---------------------------------------------------------

test_that("a single sample works and agrees with the full matrix", {
  m <- toy_counts()
  full <- ssTaxSEA(m, custom_db = toy_sets, min_set_size = 2)
  one <- ssTaxSEA(m[, 1, drop = FALSE], custom_db = toy_sets,
                  min_set_size = 2)

  expect_identical(dim(one), c(2L, 1L))
  expect_identical(rownames(one), c("A", "B"))
  expect_equal(as.vector(one), as.vector(full[, 1]))
})

## Input validation ------------------------------------------------------

test_that("invalid input is rejected with an informative message", {
  m <- toy_counts()

  expect_error(ssTaxSEA(unname(m), custom_db = toy_sets), "row names")

  neg <- m; neg[1, 1] <- -5
  expect_error(ssTaxSEA(neg, custom_db = toy_sets), "negative")

  na_mat <- m; na_mat[1, 1] <- NA
  expect_error(ssTaxSEA(na_mat, custom_db = toy_sets), "non-finite")

  empty <- m; empty[, 1] <- 0
  expect_error(ssTaxSEA(empty, custom_db = toy_sets), "zero total")

  brackets <- m; rownames(brackets)[1] <- "[Ruminococcus] gnavus"
  expect_error(ssTaxSEA(brackets, custom_db = toy_sets), "square brackets")

  expect_error(ssTaxSEA(m, custom_db = toy_sets, pseudocount = 0),
               "positive")
  expect_error(ssTaxSEA(m, custom_db = toy_sets, pseudocount = -1),
               "positive")

  expect_error(ssTaxSEA(m, custom_db = "not a list"), "list of taxon sets")
})

test_that("proportions are rejected rather than silently transformed", {
  # A pseudocount of 0.5 on values summing to 1 would swamp the data.
  props <- sweep(toy_counts(), 2, colSums(toy_counts()), FUN = "/")
  expect_error(ssTaxSEA(props, custom_db = toy_sets, min_set_size = 2),
               "proportions")
})

test_that("counts scaled to a count-like total are accepted", {
  scaled <- sweep(toy_counts(), 2, colSums(toy_counts()), FUN = "/") * 1e4
  expect_silent(
    res <- ssTaxSEA(scaled, custom_db = toy_sets, min_set_size = 2)
  )
  expect_identical(dim(res), c(2L, 4L))
})

test_that("percentages are rejected on the count scale", {
  pct <- sweep(toy_counts(), 2, colSums(toy_counts()), FUN = "/") * 100
  expect_error(ssTaxSEA(pct, custom_db = toy_sets, min_set_size = 2),
               "percentages")

  # MetaPhlAn species rows without UNCLASSIFIED sum to less than 100
  expect_error(ssTaxSEA(pct * 0.85, custom_db = toy_sets, min_set_size = 2),
               "percentages")

  # MetaPhlAn rounding can leave columns summing to just over 100
  rounded <- pct
  rounded[1, ] <- rounded[1, ] + 1e-5
  expect_error(ssTaxSEA(rounded, custom_db = toy_sets, min_set_size = 2),
               "percentages")
})

test_that("small integer count tables are not mistaken for percentages", {
  m <- toy_counts()
  m[] <- pmin(m, 20)
  expect_true(all(colSums(m) <= 100))
  expect_silent(ssTaxSEA(m, custom_db = toy_sets, min_set_size = 2))
})

## Proportion scale ------------------------------------------------------

test_that("proportion scale gives the same scores for counts, proportions and percentages", {
  m <- toy_counts()
  props <- sweep(m, 2, colSums(m), FUN = "/")

  from_counts <- ssTaxSEA(m, custom_db = toy_sets, min_set_size = 2,
                          scale = "proportion")
  from_props <- ssTaxSEA(props, custom_db = toy_sets, min_set_size = 2,
                         scale = "proportion")
  from_pct <- ssTaxSEA(props * 100, custom_db = toy_sets, min_set_size = 2,
                       scale = "proportion")

  expect_equal(from_props, from_counts)
  expect_equal(from_pct, from_counts)

  # Default pseudocount on this scale is 1e-6
  log_p <- log(props + 1e-6)
  clr <- sweep(log_p, 2, colMeans(log_p), FUN = "-")
  expect_equal(from_counts["A", ], colMeans(clr[c("t1", "t2"), ]))
})

test_that("a count-sized pseudocount on the proportion scale warns", {
  expect_warning(ssTaxSEA(toy_counts(), custom_db = toy_sets,
                          min_set_size = 2, scale = "proportion",
                          pseudocount = 0.5),
                 "dominate")
})

## Top-k reference -------------------------------------------------------

test_that("top_k reference centers on the most abundant detected taxa", {
  m <- toy_counts()
  res <- ssTaxSEA(m, custom_db = toy_sets, min_set_size = 2,
                  reference = "top_k", top_k = 2)

  log_m <- log(m + 0.5)
  ref <- apply(log_m, 2, function(x) mean(sort(x, decreasing = TRUE)[1:2]))
  clr <- sweep(log_m, 2, ref, FUN = "-")
  expect_equal(res["A", ], colMeans(clr[c("t1", "t2"), ]))
  expect_equal(res["B", ], colMeans(clr[c("t3", "t4"), ]))
})

test_that("top_k scores of surviving taxa do not move when rare taxa are lost", {
  # The motivating case: rare taxa drop out (antibiotics, chemotherapy)
  # while the survivors are unchanged. The all-taxa reference falls and
  # inflates every survivor's score; the top_k reference does not.
  set.seed(3)
  before <- sort(rpois(80, exp(rnorm(80, 3, 1.5))) + 1, decreasing = TRUE)
  after <- before
  after[41:80] <- 0
  m <- cbind(before = before, after = after)
  rownames(m) <- paste0("t", seq_len(80))
  sets <- list(survivors = paste0("t", 21:40))

  topk <- ssTaxSEA(m, custom_db = sets, reference = "top_k", top_k = 20)
  all_taxa <- ssTaxSEA(m, custom_db = sets)

  expect_equal(topk[, "after"], topk[, "before"])
  expect_gt(all_taxa[, "after"] - all_taxa[, "before"], 0.5)
})

test_that("top_k scores do not depend on extra rare or empty rows", {
  m <- toy_counts()
  padded <- rbind(m, matrix(c(0, 1, 0, 0,
                              0, 0, 0, 0), nrow = 2, byrow = TRUE,
                            dimnames = list(c("rare", "empty"), NULL)))

  topk <- ssTaxSEA(m, custom_db = toy_sets, min_set_size = 2,
                   reference = "top_k", top_k = 2)
  topk_padded <- ssTaxSEA(padded, custom_db = toy_sets, min_set_size = 2,
                          reference = "top_k", top_k = 2)
  expect_equal(topk_padded, topk)

  all_taxa <- ssTaxSEA(m, custom_db = toy_sets, min_set_size = 2)
  all_padded <- ssTaxSEA(padded, custom_db = toy_sets, min_set_size = 2)
  expect_false(isTRUE(all.equal(all_padded, all_taxa)))
})

test_that("top_k is capped at the number of detected taxa", {
  m <- toy_counts()
  m[c("t2", "t3"), "s1"] <- 0

  expect_message(
    res <- ssTaxSEA(m, custom_db = toy_sets, min_set_size = 2,
                    reference = "top_k", top_k = 3),
    "1 of 4 samples"
  )

  # s1 has two detected taxa, so its reference is their mean, not one
  # that includes a zero
  log_s1 <- log(m[, "s1"] + 0.5)
  ref_s1 <- mean(log_s1[c("t1", "t4")])
  expect_equal(res["A", "s1"], mean(log_s1[c("t1", "t2")] - ref_s1))
})

test_that("top_k is validated", {
  m <- toy_counts()
  expect_error(ssTaxSEA(m, custom_db = toy_sets, reference = "top_k",
                        top_k = 0), "top_k")
  expect_error(ssTaxSEA(m, custom_db = toy_sets, reference = "top_k",
                        top_k = 2.5), "top_k")
  expect_error(ssTaxSEA(m, custom_db = toy_sets, reference = "median"),
               "should be one of")
})

## Set filtering ---------------------------------------------------------

test_that("sets outside the size bounds are dropped", {
  m <- toy_counts()

  # min_set_size = 3 excludes both two-member sets
  expect_warning(res <- ssTaxSEA(m, custom_db = toy_sets, min_set_size = 3),
                 "No taxon sets remain")
  expect_identical(dim(res), c(0L, 4L))
  expect_identical(colnames(res), colnames(m))

  # max_set_size = 1 excludes them too
  expect_warning(ssTaxSEA(m, custom_db = toy_sets, min_set_size = 1,
                          max_set_size = 1),
                 "No taxon sets remain")
})

test_that("set members absent from the data are ignored", {
  m <- toy_counts()
  padded <- list(A = c("t1", "t2", "not_present_1", "not_present_2"))
  res <- ssTaxSEA(m, custom_db = padded, min_set_size = 2)

  expected <- ssTaxSEA(m, custom_db = list(A = c("t1", "t2")),
                       min_set_size = 2)
  expect_equal(res, expected)
})

## SummarizedExperiment input --------------------------------------------

test_that("SummarizedExperiment input matches matrix input", {
  skip_if_not_installed("SummarizedExperiment")

  m <- toy_counts()
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = m)
  )

  from_matrix <- ssTaxSEA(m, custom_db = toy_sets, min_set_size = 2)
  expect_message(
    from_se <- ssTaxSEA(se, custom_db = toy_sets, min_set_size = 2),
    "SummarizedExperiment"
  )

  expect_equal(from_se, from_matrix)
})

## Built-in database -----------------------------------------------------

test_that("the built-in database maps taxon names to NCBI IDs", {
  data("TaxSEA_test_data", package = "TaxSEA")
  data("NCBI_ids", package = "TaxSEA")

  taxa <- names(TaxSEA_test_data)
  taxa <- taxa[taxa %in% names(NCBI_ids)]

  set.seed(7)
  m <- matrix(rpois(length(taxa) * 6, lambda = 20),
              nrow = length(taxa), ncol = 6)
  rownames(m) <- taxa
  colnames(m) <- paste0("s", seq_len(6))

  res <- ssTaxSEA(m)

  expect_true(is.matrix(res))
  expect_identical(colnames(res), colnames(m))
  expect_gt(nrow(res), 0)
  expect_true(all(is.finite(res)))
  # Set names come from the database, not from the input
  expect_true(any(grepl("BacDive|MiMeDB|GMRepoV2", rownames(res))))
})
