# Test Visualization Functions
# Tests for plotIS() and plotVS()

test_that("plotIS returns a ggplot object", {
  df <- data.frame(
    MOTIF = c("AAAA", "CCCC", "GGGG", "UUUU"),
    Score = c(10, 5, 3, 2)
  )
  
  result <- plotIS(df)
  
  expect_s3_class(result, "ggplot")
})

test_that("plotIS works with specific motif", {
  df <- data.frame(
    MOTIF = c("AAAA", "CCCC", "GGGG", "UUUU"),
    Score = c(10, 5, 3, 2)
  )
  
  result <- plotIS(df, motif = "CCCC")
  
  expect_s3_class(result, "ggplot")
})

test_that("plotIS respects bins parameter", {
  df <- data.frame(
    MOTIF = c("AAAA", "CCCC", "GGGG", "UUUU"),
    Score = c(10, 5, 3, 2)
  )
  
  result <- plotIS(df, bins = 20)
  
  expect_s3_class(result, "ggplot")
})

test_that("plotIS errors on invalid input", {
  df_bad <- data.frame(Motif = c("AAAA"), Value = c(10))
  
  expect_error(plotIS(df_bad))
})

test_that("plotVS returns a ggplot object", {
  df <- data.frame(
    MOTIF = c("AAAA", "CAAA", "GAAA", "UAAA", "ACAA", "AGAA", "AUAA"),
    Score = c(10, 5, 4, 3, 6, 5, 4)
  )
  
  result <- plotVS(df)
  
  expect_s3_class(result, "ggplot")
})

test_that("plotVS works with specific motif", {
  df <- data.frame(
    MOTIF = c("CCCC", "ACCC", "GCCC", "UCCC", "CACC", "CGCC", "CUCC"),
    Score = c(10, 8, 7, 6, 5, 4, 3)
  )
  
  result <- plotVS(df, motif = "CCCC")
  
  expect_s3_class(result, "ggplot")
})

test_that("plotVS handles incomplete variant data", {
  # Motif with missing variants
  df <- data.frame(
    MOTIF = c("AAAA", "CAAA"),  # Missing GAAA, UAAA, etc.
    Score = c(10, 5)
  )
  
  result <- plotVS(df)
  
  expect_s3_class(result, "ggplot")
})
