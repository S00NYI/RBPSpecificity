# Test Inherent Specificity (IS) Functions
# Tests for returnIS()

test_that("returnIS returns numeric value for valid input", {
  df <- data.frame(
    MOTIF = c("AAAA", "CCCC", "GGGG", "UUUU"),
    Score = c(10, 5, 3, 2)
  )
  
  result <- returnIS(df)
  
  expect_type(result, "double")
  expect_true(is.finite(result))
})

test_that("returnIS calculates IS correctly", {
  # IS = top score / median score
  df <- data.frame(
    MOTIF = c("AAAA", "CCCC", "GGGG", "UUUU"),
    Score = c(10, 4, 3, 2)  # median of 4,3,2,10 = 3.5
  )
  
  result <- returnIS(df)
  expected <- 10 / median(c(10, 4, 3, 2))  # 10 / 3.5 = 2.857
  
  expect_equal(result, expected, tolerance = 0.01)
})

test_that("returnIS returns all IS values with return_type='all'", {
  df <- data.frame(
    MOTIF = c("AAAA", "CCCC", "GGGG"),
    Score = c(10, 5, 2)
  )
  
  result <- returnIS(df, return_type = "all")
  
  expect_s3_class(result, "data.frame")
  expect_true("MOTIF" %in% names(result))
  expect_true("IS" %in% names(result))
  expect_equal(nrow(result), 3)
})

test_that("returnIS works with specific motif", {
  df <- data.frame(
    MOTIF = c("AAAA", "CCCC", "GGGG"),
    Score = c(10, 5, 2)
  )
  
  result <- returnIS(df, motif = "CCCC")
  expected <- 5 / median(c(10, 5, 2))  # 5 / 5 = 1
  
  expect_equal(result, expected, tolerance = 0.01)
})

test_that("returnIS errors on invalid input", {
  # Missing required columns
  df_bad <- data.frame(Motif = c("AAAA"), Value = c(10))
  
  expect_error(returnIS(df_bad))
})

test_that("returnIS errors on empty data frame", {
  df_empty <- data.frame(MOTIF = character(), Score = numeric())
  
  expect_error(returnIS(df_empty))
})
