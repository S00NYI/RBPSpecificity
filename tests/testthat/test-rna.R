
# Tests for RNA support in countKmers

test_that("countKmers handles RNA type (DNA input converted to RNA)", {
  # Input DNA: "TTT"
  # Type: RNA
  # Converted to: "UUU"
  # K=3
  # Expected: "UUU" count 1
  
  seqs <- Biostrings::DNAStringSet("TTT")
  res <- countKmers(seqs, K = 3, type = "RNA")
  
  expect_true("UUU" %in% res$MOTIF)
  expect_equal(res$COUNT[res$MOTIF == "UUU"], 1)
  
  # Check that TTT is NOT in the motifs (since it's RNA mode)
  expect_false("TTT" %in% res$MOTIF)
})
