# Test Sensitivity Functions
# Tests for returnSensitivity()

test_that("returnSensitivity returns numeric value with output_type='number'", {
    df <- data.frame(
        MOTIF = c("AAAA", "CAAA", "GAAA", "UAAA", "ACAA", "AGAA", "AUAA"),
        Score = c(10, 5, 4, 3, 6, 5, 4)
    )

    result <- returnSensitivity(df, output_type = "number")

    expect_type(result, "double")
    expect_true(is.finite(result))
})

test_that("returnSensitivity returns matrix with output_type='matrix'", {
    df <- data.frame(
        MOTIF = c("AAAA", "CAAA", "GAAA", "UAAA", "ACAA", "AGAA", "AUAA"),
        Score = c(10, 5, 4, 3, 6, 5, 4)
    )

    result <- returnSensitivity(df, output_type = "matrix")

    expect_true(is.matrix(result) || is.data.frame(result))
})

test_that("returnSensitivity works with specific motif", {
    df <- data.frame(
        MOTIF = c("CCCC", "ACCC", "GCCC", "UCCC", "CACC", "CGCC", "CUCC"),
        Score = c(10, 8, 7, 6, 5, 4, 3)
    )

    result <- returnSensitivity(df, motif = "CCCC", output_type = "number")

    expect_type(result, "double")
})

test_that("returnSensitivity returns all sensitivity values with return_type='all'", {
    df <- data.frame(
        MOTIF = c("AAAA", "CAAA", "GAAA", "UAAA"),
        Score = c(10, 5, 4, 3)
    )

    result <- returnSensitivity(df, return_type = "all")

    expect_s3_class(result, "data.frame")
    expect_true("MOTIF" %in% names(result))
    expect_true("Sensitivity" %in% names(result))
})

test_that("returnSensitivity errors on invalid input", {
    # Missing required columns
    df_bad <- data.frame(Motif = c("AAAA"), Value = c(10))

    expect_error(returnSensitivity(df_bad))
})

test_that("returnSensitivity errors on empty data frame", {
    df_empty <- data.frame(MOTIF = character(), Score = numeric())

    expect_error(returnSensitivity(df_empty))
})

test_that("returnSensitivity handles motifs with no variants gracefully", {
    # Only one motif - no variants to compare
    df <- data.frame(
        MOTIF = c("AAAA"),
        Score = c(10)
    )

    # Should return NA or handle gracefully
    # Since genMotifVar will generate variants but none are in df, it should warning/return NA
    expect_warning(result <- returnSensitivity(df, output_type = "number"))
    expect_true(is.na(result))
})
