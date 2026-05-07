# --- Test fixtures ---
mock_peak_counts <- data.frame(
    MOTIF = c("AAA", "AAC", "AAG", "AAT"),
    COUNT = c(100L, 50L, 30L, 10L),
    SEQ_WITH_MOTIF = c(80L, 40L, 25L, 8L),
    SEQ_TOTAL = c(200L, 200L, 200L, 200L),
    stringsAsFactors = FALSE
)

mock_bkg_data <- list(
    MOTIF = c("AAA", "AAC", "AAG", "AAT"),
    bkg_total_count = c(5000L, 3000L, 2000L, 1000L),
    bkg_presence_count = c(3000L, 2000L, 1500L, 800L),
    bkg_total_seqs = 20000L
)


# --- calEnrichment tests ---
test_that("calEnrichment returns 10-column data.frame", {
    result <- calEnrichment(
        mock_peak_counts, mock_bkg_data, method = "zoops"
    )
    expect_true(is.data.frame(result))
    expect_equal(ncol(result), 10)
    expect_named(result, c(
        "MOTIF", "Score", "log2FC", "pvalue", "padj",
        "target_fraction", "bkg_fraction",
        "target_occurence", "bkg_occurence", "multiplicity"
    ))
})

test_that("calEnrichment ZOOPS produces valid output", {
    result <- calEnrichment(
        mock_peak_counts, mock_bkg_data, method = "zoops"
    )

    # Score should be [0, 1]
    expect_true(all(result$Score >= 0 & result$Score <= 1))

    # Fractions should be in [0, 1]
    expect_true(all(
        result$target_fraction >= 0 &
        result$target_fraction <= 1
    ))
    expect_true(all(
        result$bkg_fraction >= 0 &
        result$bkg_fraction <= 1
    ))

    # p-values should be in [0, 1]
    expect_true(all(result$pvalue >= 0 & result$pvalue <= 1))
    expect_true(all(result$padj >= 0 & result$padj <= 1))

    # padj >= pvalue (BH only inflates)
    expect_true(all(result$padj >= result$pvalue))
})

test_that("calEnrichment ANR produces valid output", {
    result <- calEnrichment(
        mock_peak_counts, mock_bkg_data, method = "anr"
    )

    expect_equal(ncol(result), 10)
    expect_true(all(result$Score >= 0 & result$Score <= 1))
    expect_true(all(result$pvalue >= 0 & result$pvalue <= 1))
})

test_that("calEnrichment ZOOPS and ANR give different pvalues", {
    zoops <- calEnrichment(
        mock_peak_counts, mock_bkg_data, method = "zoops"
    )
    anr <- calEnrichment(
        mock_peak_counts, mock_bkg_data, method = "anr"
    )

    # Different statistical tests -> different p-values
    expect_false(all(zoops$pvalue == anr$pvalue))
    # But same raw data columns
    expect_equal(zoops$target_fraction, anr$target_fraction)
    expect_equal(zoops$target_occurence, anr$target_occurence)
    expect_equal(zoops$multiplicity, anr$multiplicity)
})

test_that("calEnrichment ZOOPS phyper matches manual", {
    result <- calEnrichment(
        mock_peak_counts, mock_bkg_data, method = "zoops"
    )

    # Manual calculation for AAA
    eps <- 1e-4
    expected_log2fc <- log2((80/200 + eps) / (3000/20000 + eps))

    expect_equal(
        result$log2FC[1], expected_log2fc,
        tolerance = 1e-6
    )

    expected_pval <- phyper(
        q = 80L - 1L,
        m = 80L + 3000L,
        n = (200L - 80L) + (20000L - 3000L),
        k = 200L,
        lower.tail = FALSE
    )
    expect_equal(
        result$pvalue[1], expected_pval,
        tolerance = 1e-10
    )
})

test_that("calEnrichment ANR pbinom matches manual", {
    result <- calEnrichment(
        mock_peak_counts, mock_bkg_data, method = "anr"
    )

    eps <- 1e-4
    expected_log2fc <- log2((100/200 + eps) / (5000/20000 + eps))

    expect_equal(
        result$log2FC[1], expected_log2fc,
        tolerance = 1e-6
    )

    p_null <- 200 / (200 + 20000)
    expected_pval <- pbinom(
        q = 100L - 1L,
        size = 100L + 5000L,
        prob = p_null,
        lower.tail = FALSE
    )
    expect_equal(
        result$pvalue[1], expected_pval,
        tolerance = 1e-10
    )
})

test_that("calEnrichment multiplicity is correct", {
    result <- calEnrichment(
        mock_peak_counts, mock_bkg_data, method = "zoops"
    )

    # AAA: 100 total / 80 sequences with motif = 1.25
    expect_equal(result$multiplicity[1], 1.25)
    # AAC: 50 / 40 = 1.25
    expect_equal(result$multiplicity[2], 1.25)
})

test_that("calEnrichment ZOOPS Score uses fraction difference", {
    result <- calEnrichment(
        mock_peak_counts, mock_bkg_data, method = "zoops"
    )

    # AAA: target_fraction=80/200=0.4, bkg_fraction=3000/20000=0.15
    # AAT: target_fraction=8/200=0.04, bkg_fraction=800/20000=0.04
    # AAA has highest difference -> should get Score = 1
    # AAT has lowest difference (0) -> should get Score = 0
    expect_equal(result$Score[1], 1.0)  # AAA
    expect_equal(result$Score[4], 0.0)  # AAT
})

test_that("calEnrichment ANR Score uses rate difference", {
    result <- calEnrichment(
        mock_peak_counts, mock_bkg_data, method = "anr"
    )

    # AAA: target_rate=100/200=0.5, bkg_rate=5000/20000=0.25
    # Highest rate difference -> Score = 1
    expect_equal(result$Score[1], 1.0)
})

test_that("calEnrichment rejects invalid method", {
    expect_error(
        calEnrichment(
            mock_peak_counts, mock_bkg_data,
            method = "subtract"
        ),
        "arg"
    )
})
