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

# ------------------------------------------------------------------------------
# Tests for getSequence (Internal Function)
# ------------------------------------------------------------------------------

# Note: getSequence is internal, so we access it if needed, but testthat usually runs in the pkg env.
# We need a BSgenome object. We'll use BSgenome.Hsapiens.UCSC.hg38 as it is available.
if (requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
  hg38 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  
  test_that("getSequence handles bidirectional extension (c(X, Y)) correctly on + strand", {
    # Create a 10bp range on chr1 + strand
    # chr1: 100-109 (width 10)
    gr <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges::IRanges(start = 100, end = 109), strand = "+")
    
    # Extend 5' by 5, 3' by 10
    # Expected: 5' (start) - 5 = 95
    # Expected: 3' (end) + 10 = 119
    # Final width: 119 - 95 + 1 = 25
    seqs <- getSequence(gr, hg38, extension = c(5, 10), min_length = 1)
    
    expect_equal(Biostrings::width(seqs), 25)
  })
  
  test_that("getSequence handles bidirectional extension correctly on - strand", {
    # Create a 10bp range on chr1 - strand
    # chr1: 100-109 (width 10)
    # 5' end is at 109, 3' end is at 100
    gr <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges::IRanges(start = 100, end = 109), strand = "-")
    
    # Extend 5' by 5, 3' by 10
    # Expected: 5' (end) + 5 = 114 (new end)
    # Expected: 3' (start) - 10 = 90 (new start)
    # Final width: 114 - 90 + 1 = 25
    seqs <- getSequence(gr, hg38, extension = c(5, 10), min_length = 1)
    
    expect_equal(Biostrings::width(seqs), 25)
  })
  
  test_that("getSequence handles trimming (negative extension)", {
    # Create a 20bp range
    gr <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges::IRanges(start = 100, end = 119), strand = "+")
    
    # Trim 5' by 2, 3' by 3
    # Expected: start + 2 = 102
    # Expected: end - 3 = 116
    # Final width: 116 - 102 + 1 = 15
    seqs <- getSequence(gr, hg38, extension = c(-2, -3), min_length = 1)
    
    expect_equal(Biostrings::width(seqs), 15)
  })
  
  test_that("getSequence filters ranges shorter than min_length", {
    # Create two ranges
    # gr1: 10bp
    # gr2: 10bp
    gr <- GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges = IRanges::IRanges(start = c(100, 200), width = 10),
      strand = "+"
    )
    
    # Trim both by 8bp total (4 from each side) -> width 2
    # Set min_length = 5.
    # Result: both should be removed.
    
    expect_message(
      seqs <- getSequence(gr, hg38, extension = c(-4, -4), min_length = 5),
      "ranges removed due to final length < 5"
    )
    
    expect_equal(length(seqs), 0)
    
    # Try one valid, one invalid
    # gr1: 10bp -> becomes 2bp (invalid)
    # gr2: 20bp -> becomes 12bp (valid)
    gr2 <- GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges = IRanges::IRanges(start = c(100, 200), width = c(10, 20)),
      strand = "+"
    )
    
    expect_message(
      seqs2 <- getSequence(gr2, hg38, extension = c(-4, -4), min_length = 5),
      "1 ranges removed due to final length < 5"
    )
    
    expect_equal(length(seqs2), 1)
    # If the second one is kept, its width: 20 -> 20 - 4 - 4 = 12
    expect_equal(Biostrings::width(seqs2), 12)
  })
  
  test_that("getSequence handles c(0,0) (no change)", {
    gr <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges::IRanges(start = 100, width = 10), strand = "+")
    seqs <- getSequence(gr, hg38, extension = c(0, 0))
    expect_equal(Biostrings::width(seqs), 10)
  })
  
  test_that("getSequence errors on invalid extension input", {
     gr <- GenomicRanges::GRanges(seqnames = "chr1", ranges = IRanges::IRanges(start = 100, width = 10), strand = "+")
     expect_error(getSequence(gr, hg38, extension = 5), "'extension' must be a numeric vector of length 2")
  })
}
