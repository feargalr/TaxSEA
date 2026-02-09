test_that("taxon_rank_sets returns a named list of data frames for all ranks", {
  # Create mock lineage data with multiple taxonomic levels
  lineage_df <- data.frame(
    species = paste0("sp", seq_len(12)),
    genus = rep(c("genusA", "genusB"), each = 6),
    family = rep(c("famA", "famB"), each = 6),
    order = rep(c("ordA", "ordB"), each = 6),
    stringsAsFactors = FALSE
  )

  # Create fold changes: first 6 species positive, last 6 negative
  set.seed(42)
  fc <- c(rnorm(6, mean = 2), rnorm(6, mean = -2))
  names(fc) <- lineage_df$species

  res <- taxon_rank_sets(fc, lineage_df, min_set_size = 3)

  # Should return a list with one element per rank (genus, family, order)
  expect_type(res, "list")
  expect_named(res, c("genus", "family", "order"))

  # Each element should be a data frame with the expected columns
  expected_cols <- c("taxonSetName", "median_rank_of_set_members",
                     "PValue", "Test_statistic", "FDR")
  for (rank_name in names(res)) {
    expect_s3_class(res[[rank_name]], "data.frame")
    expect_true(all(expected_cols %in% colnames(res[[rank_name]])))
  }

  # With min_set_size = 3, each rank should have results
  expect_true(nrow(res$genus) > 0)
  expect_true(nrow(res$family) > 0)
  expect_true(nrow(res$order) > 0)
})

test_that("taxon_rank_sets validates inputs correctly", {
  lineage_df <- data.frame(
    species = c("sp1", "sp2"),
    genus = c("g1", "g2"),
    stringsAsFactors = FALSE
  )
  fc <- c(sp1 = 1.0, sp2 = -1.0)

  # Missing species column
  bad_df <- data.frame(taxon = c("sp1", "sp2"), genus = c("g1", "g2"))
  expect_error(taxon_rank_sets(fc, bad_df),
               "species")

  # No overlap
  fc_no_match <- c(spX = 1.0, spY = -1.0)
  expect_error(taxon_rank_sets(fc_no_match, lineage_df),
               "No overlap")

  # No rank columns besides species
  species_only <- data.frame(species = c("sp1", "sp2"))
  expect_error(taxon_rank_sets(fc, species_only),
               "at least one taxonomic rank column")
})

test_that("taxon_rank_sets handles ranks with no valid groups gracefully", {
  lineage_df <- data.frame(
    species = paste0("sp", seq_len(12)),
    genus = rep(c("genusA", "genusB"), each = 6),
    family = rep(c("famA", "famB"), each = 6),
    stringsAsFactors = FALSE
  )

  set.seed(42)
  fc <- c(rnorm(6, mean = 2), rnorm(6, mean = -2))
  names(fc) <- lineage_df$species

  # With high min_set_size, some ranks may return empty results
  res <- taxon_rank_sets(fc, lineage_df, min_set_size = 50)

  expect_type(res, "list")
  expect_named(res, c("genus", "family"))

  # Results should still be data frames even if empty
  for (rank_name in names(res)) {
    expect_s3_class(res[[rank_name]], "data.frame")
  }
})

test_that("taxon_rank_sets works with SummarizedExperiment input", {
  skip_if_not_installed("SummarizedExperiment")

  # Create a SummarizedExperiment with taxonomic rowData
  feature_names <- paste0("OTU", seq_len(12))

  row_data <- S4Vectors::DataFrame(
    Kingdom = rep("Bacteria", 12),
    Phylum = rep(c("Firmicutes", "Proteobacteria"), each = 6),
    Class = rep(c("Bacilli", "Gammaproteobacteria"), each = 6),
    Order = rep(c("Bacillales", "Enterobacterales"), each = 6),
    Family = rep(c("Staphylococcaceae", "Enterobacteriaceae"), each = 6),
    Genus = rep(c("Staphylococcus", "Escherichia"), each = 6),
    Species = paste0("Species_", seq_len(12))
  )
  rownames(row_data) <- feature_names

  # Create a minimal count matrix
  counts <- matrix(rpois(12 * 4, lambda = 100), nrow = 12, ncol = 4)
  rownames(counts) <- feature_names
  colnames(counts) <- paste0("Sample", seq_len(4))

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowData = row_data
  )

  # Create fold changes with OTU IDs as names (matching rownames)
  set.seed(42)
  fc <- c(rnorm(6, mean = 2), rnorm(6, mean = -2))
  names(fc) <- feature_names

  res <- taxon_rank_sets(fc, se, min_set_size = 3)

  # Should return a list with rank names (Species excluded)
  expect_type(res, "list")
  expect_true(!"Species" %in% names(res))
  expect_true(all(c("Kingdom", "Phylum", "Class", "Order",
                     "Family", "Genus") %in% names(res)))

  # Each element should be a data frame
  expected_cols <- c("taxonSetName", "median_rank_of_set_members",
                     "PValue", "Test_statistic", "FDR")
  for (rank_name in names(res)) {
    expect_s3_class(res[[rank_name]], "data.frame")
    expect_true(all(expected_cols %in% colnames(res[[rank_name]])))
  }
})

test_that("taxon_rank_sets SE results match equivalent data frame results", {
  skip_if_not_installed("SummarizedExperiment")

  # Create matching data frame and SE inputs
  feature_names <- paste0("feat", seq_len(12))

  lineage_df <- data.frame(
    species = feature_names,
    Phylum = rep(c("PhylumA", "PhylumB"), each = 6),
    Genus = rep(c("GenusA", "GenusB"), each = 6),
    stringsAsFactors = FALSE
  )

  row_data <- S4Vectors::DataFrame(
    Phylum = rep(c("PhylumA", "PhylumB"), each = 6),
    Genus = rep(c("GenusA", "GenusB"), each = 6),
    Species = feature_names
  )
  rownames(row_data) <- feature_names

  counts <- matrix(1L, nrow = 12, ncol = 2)
  rownames(counts) <- feature_names
  colnames(counts) <- c("S1", "S2")

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowData = row_data
  )

  set.seed(42)
  fc <- c(rnorm(6, mean = 2), rnorm(6, mean = -2))
  names(fc) <- feature_names

  res_df <- taxon_rank_sets(fc, lineage_df, min_set_size = 3)
  res_se <- taxon_rank_sets(fc, se, min_set_size = 3)

  # Both should have the same ranks (Phylum and Genus)
  expect_true("Phylum" %in% names(res_df))
  expect_true("Genus" %in% names(res_df))
  expect_true("Phylum" %in% names(res_se))
  expect_true("Genus" %in% names(res_se))

  # P-values should match
  expect_equal(res_df$Phylum$PValue, res_se$Phylum$PValue)
  expect_equal(res_df$Genus$PValue, res_se$Genus$PValue)
})
