# Test Mutational Sensitivity (MS) Functions
# Tests for returnMS()

test_that("returnMS returns numeric value with output_type='number'", {
  df <- data.frame(
    MOTIF = c("AAAA", "CAAA", "GAAA", "UAAA", "ACAA", "AGAA", "AUAA"),
    Score = c(10, 5, 4, 3, 6, 5, 4)
  )
  
  result <- returnMS(df, output_type = "number")
  
  expect_type(result, "double")
  expect_true(is.finite(result))
})

test_that("returnMS returns matrix with output_type='matrix'", {
  df <- data.frame(
    MOTIF = c("AAAA", "CAAA", "GAAA", "UAAA", "ACAA", "AGAA", "AUAA"),
    Score = c(10, 5, 4, 3, 6, 5, 4)
  )
  
  result <- returnMS(df, output_type = "matrix")
  
  expect_true(is.matrix(result) || is.data.frame(result))
})

test_that("returnMS works with specific motif", {
  df <- data.frame(
    MOTIF = c("CCCC", "ACCC", "GCCC", "UCCC", "CACC", "CGCC", "CUCC"),
    Score = c(10, 8, 7, 6, 5, 4, 3)
  )
  
  result <- returnMS(df, motif = "CCCC", output_type = "number")
  
  expect_type(result, "double")
})

test_that("returnMS returns all MS values with return_type='all'", {
  df <- data.frame(
    MOTIF = c("AAAA", "CAAA", "GAAA", "UAAA"),
    Score = c(10, 5, 4, 3)
  )
  
  result <- returnMS(df, return_type = "all", output_type = "number")
  
  expect_s3_class(result, "data.frame")
  expect_true("MOTIF" %in% names(result))
  expect_true("MS" %in% names(result))
})

test_that("returnMS errors on invalid input", {
  # Missing required columns
  df_bad <- data.frame(Motif = c("AAAA"), Value = c(10))
  
  expect_error(returnMS(df_bad))
})

test_that("returnMS errors on empty data frame", {
  df_empty <- data.frame(MOTIF = character(), Score = numeric())
  
  expect_error(returnMS(df_empty))
})

test_that("returnMS handles motifs with no variants gracefully", {
  # Only one motif - no variants to compare
  df <- data.frame(
    MOTIF = c("AAAA"),
    Score = c(10)
  )
  
  # Should return NA or handle gracefully
  result <- returnMS(df, output_type = "number")
  expect_true(is.na(result) || is.finite(result))
})
