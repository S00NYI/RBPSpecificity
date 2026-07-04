# Test Visualization Functions
# Tests for plotSpecificity() and plotSensitivity()

test_that("plotSpecificity returns a ggplot object", {
    df <- data.frame(
        MOTIF = c("AAAA", "CCCC", "GGGG", "UUUU"),
        Score = c(10, 5, 3, 2)
    )

    result <- plotSpecificity(df)

    expect_s3_class(result, "ggplot")
})

test_that("plotSpecificity works with specific motif", {
    df <- data.frame(
        MOTIF = c("AAAA", "CCCC", "GGGG", "UUUU"),
        Score = c(10, 5, 3, 2)
    )

    result <- plotSpecificity(df, motif = "CCCC")

    expect_s3_class(result, "ggplot")
})

test_that("plotSpecificity respects bins parameter", {
    df <- data.frame(
        MOTIF = c("AAAA", "CCCC", "GGGG", "UUUU"),
        Score = c(10, 5, 3, 2)
    )

    result <- plotSpecificity(df, bins = 20)

    expect_s3_class(result, "ggplot")
})

test_that("plotSpecificity errors on invalid input", {
    df_bad <- data.frame(Motif = c("AAAA"), Value = c(10))

    expect_error(plotSpecificity(df_bad))
})

test_that("plotSensitivity returns a ggplot object", {
    df <- data.frame(
        MOTIF = c("AAAA", "CAAA", "GAAA", "UAAA", "ACAA", "AGAA", "AUAA"),
        Score = c(10, 5, 4, 3, 6, 5, 4)
    )

    result <- plotSensitivity(df)

    expect_s3_class(result, "ggplot")
})

test_that("plotSensitivity works with specific motif", {
    df <- data.frame(
        MOTIF = c("CCCC", "ACCC", "GCCC", "UCCC", "CACC", "CGCC", "CUCC"),
        Score = c(10, 8, 7, 6, 5, 4, 3)
    )

    result <- plotSensitivity(df, motif = "CCCC")

    expect_s3_class(result, "ggplot")
})

test_that("plotSensitivity handles incomplete variant data", {
    # Motif with missing variants
    df <- data.frame(
        MOTIF = c("AAAA", "CAAA"), # Missing GAAA, UAAA, etc.
        Score = c(10, 5)
    )

    result <- plotSensitivity(df)

    expect_s3_class(result, "ggplot")
})
