# Test Variation Sensitivity (VS) Functions
# Tests for returnVS()

test_that("returnVS returns numeric value with output_type='number'", {
  df <- data.frame(
    MOTIF = c("AAAA", "CAAA", "GAAA", "UAAA", "ACAA", "AGAA", "AUAA"),
    Score = c(10, 5, 4, 3, 6, 5, 4)
  )
  
  result <- returnVS(df, output_type = "number")
  
  expect_type(result, "double")
  expect_true(is.finite(result))
})

test_that("returnVS returns matrix with output_type='matrix'", {
  df <- data.frame(
    MOTIF = c("AAAA", "CAAA", "GAAA", "UAAA", "ACAA", "AGAA", "AUAA"),
    Score = c(10, 5, 4, 3, 6, 5, 4)
  )
  
  result <- returnVS(df, output_type = "matrix")
  
  expect_true(is.matrix(result) || is.data.frame(result))
})

test_that("returnVS works with specific motif", {
  df <- data.frame(
    MOTIF = c("CCCC", "ACCC", "GCCC", "UCCC", "CACC", "CGCC", "CUCC"),
    Score = c(10, 8, 7, 6, 5, 4, 3)
  )
  
  result <- returnVS(df, motif = "CCCC", output_type = "number")
  
  expect_type(result, "double")
})

test_that("returnVS returns all VS values with return_type='all'", {
  df <- data.frame(
    MOTIF = c("AAAA", "CAAA", "GAAA", "UAAA"),
    Score = c(10, 5, 4, 3)
  )
  
  result <- returnVS(df, return_type = "all", output_type = "number")
  
  expect_s3_class(result, "data.frame")
  expect_true("MOTIF" %in% names(result))
  expect_true("VS" %in% names(result))
})

test_that("returnVS errors on invalid input", {
  # Missing required columns
  df_bad <- data.frame(Motif = c("AAAA"), Value = c(10))
  
  expect_error(returnVS(df_bad))
})

test_that("returnVS errors on empty data frame", {
  df_empty <- data.frame(MOTIF = character(), Score = numeric())
  
  expect_error(returnVS(df_empty))
})

test_that("returnVS handles motifs with no variants gracefully", {
  # Only one motif - no variants to compare
  df <- data.frame(
    MOTIF = c("AAAA"),
    Score = c(10)
  )
  
  # Should return NA or handle gracefully
  result <- returnVS(df, output_type = "number")
  expect_true(is.na(result) || is.finite(result))
})

