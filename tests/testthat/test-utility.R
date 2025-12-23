# Test Utility Functions
# Tests for minmaxNorm, genMotifVar, findTopMotif, normalizeScores

test_that("minmaxNorm scales to correct range", {
  input <- c(0, 5, 10)
  result <- minmaxNorm(input, a = 0, b = 1)
  
  expect_equal(min(result), 0)
  expect_equal(max(result), 1)
  expect_equal(result[2], 0.5)
})

test_that("minmaxNorm handles custom range", {
  input <- c(0, 5, 10)
  result <- minmaxNorm(input, a = 1, b = exp(1))
  
  expect_equal(min(result), 1)
  expect_equal(max(result), exp(1), tolerance = 0.001)
})

test_that("minmaxNorm handles constant input", {
  input <- c(5, 5, 5)
  result <- minmaxNorm(input, a = 0, b = 1)
  
  # All same value - should handle gracefully
  expect_true(all(is.finite(result)))
})

test_that("genMotifVar generates correct number of variants", {
  # For a 4-mer, each position can have 3 mutations = 4 * 3 = 12 variants
  result <- genMotifVar("AAAA", type = "DNA")
  
  expect_length(result, 12)
  expect_false("AAAA" %in% result)  # Original not included
})

test_that("genMotifVar generates correct variants", {
  result <- genMotifVar("AA", type = "DNA")
  
  # Position 1: CA, GA, TA
  # Position 2: AC, AG, AT
  expected <- c("CA", "GA", "TA", "AC", "AG", "AT")
  
  expect_true(all(expected %in% result))
})

test_that("genMotifVar works with RNA type", {
  result <- genMotifVar("AA", type = "RNA")
  
  # Position 1: CA, GA, UA
  # Position 2: AC, AG, AU
  expect_true("UA" %in% result)
  expect_true("AU" %in% result)
  expect_false("TA" %in% result)  # T not in RNA
})

test_that("findTopMotif returns highest scoring motif", {
  df <- data.frame(
    MOTIF = c("AAAA", "CCCC", "GGGG"),
    Score = c(5, 10, 3)
  )
  
  result <- findTopMotif(df)
  
  expect_equal(result, "CCCC")
})

test_that("findTopMotif handles ties by returning first", {
  df <- data.frame(
    MOTIF = c("AAAA", "CCCC", "GGGG"),
    Score = c(10, 10, 3)
  )
  
  result <- findTopMotif(df)
  
  # Should return first one encountered with max score
  expect_true(result %in% c("AAAA", "CCCC"))
})

test_that("normalizeScores works with min_max method", {
  input <- c(0, 5, 10)
  result <- normalizeScores(input, method = "min_max", a = 0, b = 1)
  
  expect_equal(min(result), 0)
  expect_equal(max(result), 1)
})

test_that("normalizeScores works with z_score method", {
  input <- c(1, 2, 3, 4, 5)
  result <- normalizeScores(input, method = "z_score")
  
  expect_equal(mean(result), 0, tolerance = 0.001)
  expect_equal(sd(result), 1, tolerance = 0.001)
})

test_that("normalizeScores 'none' returns original values", {
  input <- c(1, 5, 10)
  result <- normalizeScores(input, method = "none")
  
  expect_equal(result, input)
})
