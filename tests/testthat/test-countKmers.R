test_that("countKmers returns correct columns", {
    seqs <- Biostrings::DNAStringSet(c("AAACCC", "AAAGGG"))
    result <- countKmers(seqs, K = 3, type = "DNA")

    expect_true(is.data.frame(result))
    expect_named(
        result,
        c("MOTIF", "COUNT", "SEQ_WITH_MOTIF", "SEQ_TOTAL")
    )
    expect_equal(nrow(result), 4^3)
})

test_that("countKmers SEQ_TOTAL matches input length", {
    seqs <- Biostrings::DNAStringSet(c("AAACCC", "AAAGGG", "TTTAAA"))
    result <- countKmers(seqs, K = 3, type = "DNA")

    expect_true(all(result$SEQ_TOTAL == 3L))
})

test_that("countKmers SEQ_WITH_MOTIF counts presence correctly", {
    # AAA appears in seq1 and seq2, but NOT in seq3
    seqs <- Biostrings::DNAStringSet(c("AAACCC", "AAAGGG", "TTTCCC"))
    result <- countKmers(seqs, K = 3, type = "DNA")

    aaa_row <- result[result$MOTIF == "AAA", ]
    expect_equal(aaa_row$SEQ_WITH_MOTIF, 2L)
    # AAACCC -> AAA at pos1, AAC at pos2, ACC at pos3, CCC at pos4
    # So AAA appears once per sequence, twice total
    expect_equal(aaa_row$COUNT, 2L)
})

test_that("countKmers distinguishes COUNT from SEQ_WITH_MOTIF", {
    # AAAAAA has 4 overlapping AAA occurrences in one seq
    seqs <- Biostrings::DNAStringSet(c("AAAAAA"))
    result <- countKmers(seqs, K = 3, type = "DNA")

    aaa_row <- result[result$MOTIF == "AAA", ]
    # ANR: 4 overlapping AAA in AAAAAA
    expect_equal(aaa_row$COUNT, 4L)
    # ZOOPS: only 1 sequence contains it
    expect_equal(aaa_row$SEQ_WITH_MOTIF, 1L)
})

test_that("countKmers handles empty input", {
    seqs <- Biostrings::DNAStringSet()
    result <- countKmers(seqs, K = 3, type = "DNA")

    expect_equal(nrow(result), 4^3)
    expect_true(all(result$COUNT == 0L))
    expect_true(all(result$SEQ_WITH_MOTIF == 0L))
    expect_true(all(result$SEQ_TOTAL == 0L))
})

test_that("countKmers RNA mode converts T to U", {
    seqs <- Biostrings::DNAStringSet(c("AAATTT"))
    result <- countKmers(seqs, K = 3, type = "RNA")

    expect_true("UUU" %in% result$MOTIF)
    expect_false("TTT" %in% result$MOTIF)
    uuu_row <- result[result$MOTIF == "UUU", ]
    expect_equal(uuu_row$COUNT, 1L)
    expect_equal(uuu_row$SEQ_WITH_MOTIF, 1L)
})
