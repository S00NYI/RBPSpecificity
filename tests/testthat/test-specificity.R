# Test Specificity Functions
# Tests for returnSpecificity()

test_that("returnSpecificity returns numeric value for valid input", {
    df <- data.frame(
        MOTIF = c("AAAA", "CCCC", "GGGG", "UUUU"),
        Score = c(10, 5, 3, 2)
    )

    result <- returnSpecificity(df)

    expect_type(result, "double")
    expect_true(is.finite(result))
})

test_that("returnSpecificity calculates specificity correctly", {
    # specificity = top score / median score
    df <- data.frame(
        MOTIF = c("AAAA", "CCCC", "GGGG", "UUUU"),
        Score = c(10, 4, 3, 2) # median of 4,3,2,10 = 3.5
    )

    result <- returnSpecificity(df)
    expected <- 10 / median(c(10, 4, 3, 2)) # 10 / 3.5 = 2.857

    expect_equal(result, expected, tolerance = 0.01)
})

test_that("returnSpecificity returns all specificity values with return_type='all'", {
    df <- data.frame(
        MOTIF = c("AAAA", "CCCC", "GGGG"),
        Score = c(10, 5, 2)
    )

    result <- returnSpecificity(df, return_type = "all")

    expect_s3_class(result, "data.frame")
    expect_true("MOTIF" %in% names(result))
    expect_true("Specificity" %in% names(result))
    expect_equal(nrow(result), 3)
})

test_that("returnSpecificity works with specific motif", {
    df <- data.frame(
        MOTIF = c("AAAA", "CCCC", "GGGG"),
        Score = c(10, 5, 2)
    )

    result <- returnSpecificity(df, motif = "CCCC")
    expected <- 5 / median(c(10, 5, 2)) # 5 / 5 = 1

    expect_equal(result, expected, tolerance = 0.01)
})

test_that("returnSpecificity errors on invalid input", {
    # Missing required columns
    df_bad <- data.frame(Motif = c("AAAA"), Value = c(10))

    expect_error(returnSpecificity(df_bad))
})

test_that("returnSpecificity errors on empty data frame", {
    df_empty <- data.frame(MOTIF = character(), Score = numeric())

    expect_error(returnSpecificity(df_empty))
})
