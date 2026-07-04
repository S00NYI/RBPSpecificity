# Test motifEnrichment Function
# Tests for motifEnrichment() and related helpers

test_that("motifEnrichment validates input arguments correctly", {
    # Missing arguments
    expect_error(motifEnrichment())
    expect_error(motifEnrichment(K = 5))

    # Invalid coordinates
    expect_error(motifEnrichment(coordinates = data.frame(a = 1), species_or_build = "hg38", K = 5))
})

test_that("motifEnrichment runs successfully with hg38 genome", {
    skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")

    # Create a small dummy peak coordinate set on chr9 (matching ENCODE HNRNPC coordinates)
    peaks <- data.frame(
        chr = c("chr9", "chr9"),
        start = c(120760581, 120760600),
        end = c(120760630, 120760650),
        strand = c("-", "-")
    )

    # Run motifEnrichment with low iterations for speed
    result <- motifEnrichment(
        coordinates = peaks,
        species_or_build = "hg38",
        K = 5,
        method = "anr",
        bkg_iter = 5,
        bkg_min_dist = 100,
        bkg_max_dist = 200,
        scramble_bkg = FALSE
    )

    expect_s3_class(result, "data.frame")
    expect_equal(ncol(result), 10)
    expect_true(all(c("MOTIF", "Score", "pvalue", "padj") %in% colnames(result)))
    expect_true(all(result$Score >= 0 & result$Score <= 1))
})

test_that("motifEnrichment works in both ZOOPS and ANR modes", {
    skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")

    peaks <- data.frame(
        chr = c("chr9"),
        start = c(120760581),
        end = c(120760645),
        strand = c("-")
    )

    # ZOOPS mode
    res_zoops <- motifEnrichment(
        coordinates = peaks,
        species_or_build = "hg38",
        K = 5,
        method = "zoops",
        bkg_iter = 2,
        bkg_min_dist = 50,
        bkg_max_dist = 100
    )

    # ANR mode
    res_anr <- motifEnrichment(
        coordinates = peaks,
        species_or_build = "hg38",
        K = 5,
        method = "anr",
        bkg_iter = 2,
        bkg_min_dist = 50,
        bkg_max_dist = 100
    )

    expect_s3_class(res_zoops, "data.frame")
    expect_s3_class(res_anr, "data.frame")
    expect_equal(nrow(res_zoops), 1024) # 4^5 possible 5-mers
    expect_equal(nrow(res_anr), 1024)
})
