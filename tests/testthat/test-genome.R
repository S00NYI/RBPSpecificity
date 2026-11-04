# Tests for genome selection and genome-aware chromosome filtering

test_that("selectGenome validates input", {
    expect_error(selectGenome(NA_character_), "must be a BSgenome object")
    expect_error(selectGenome(c("hg38", "mm10")), "must be a BSgenome object")
    expect_error(selectGenome(5), "must be a BSgenome object")
    expect_error(selectGenome("not_a_real_genome_xyz"), "Could not load a genome")
})

test_that("selectGenome accepts objects, shortcuts and package names", {
    skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")
    hg38 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

    expect_identical(selectGenome(hg38), hg38)
    expect_equal(genomeLabel(selectGenome("HG38")), "hg38")
    expect_equal(genomeLabel(selectGenome("human")), "hg38")
    expect_equal(
        genomeLabel(selectGenome("BSgenome.Hsapiens.UCSC.hg38")), "hg38"
    )
})

test_that("mapToGenomeSeqlevels matches names across naming styles", {
    skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")
    hg38 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

    expect_equal(
        mapToGenomeSeqlevels(c("chr1", "1", "MT", "X", "scaffold_1"), hg38),
        c("chr1", "chr1", "chrM", "chrX", NA)
    )
    # Nothing matches: still returns one NA per input
    expect_equal(mapToGenomeSeqlevels(c("foo", "bar"), hg38), c(NA_character_, NA_character_))
})

test_that("resolveChromosomes handles keywords and name vectors", {
    skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")
    hg38 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

    std <- resolveChromosomes("standard", hg38)
    expect_true(all(c("chr1", "chrX", "chrY", "chrM") %in% std))
    expect_false(any(grepl("_", std)))

    expect_equal(resolveChromosomes("all", hg38), GenomeInfoDb::seqlevels(hg38))
    expect_equal(resolveChromosomes(c("1", "chrX"), hg38), c("chr1", "chrX"))

    expect_warning(
        res <- resolveChromosomes(c("chr1", "chr99"), hg38),
        "not found in genome"
    )
    expect_equal(res, "chr1")
    expect_error(resolveChromosomes("chr99", hg38), "None of the requested chromosomes")
    expect_error(resolveChromosomes(character(0), hg38), "'chromosomes' must be")
    expect_error(resolveChromosomes(NA_character_, hg38), "'chromosomes' must be")
})

test_that("filterPeakChromosomes renames, filters and summarizes", {
    skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")
    hg38 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

    gr <- GenomicRanges::GRanges(
        c("1", "MT", "chr2", "scaffold_1", "chrUn_KI270742v1"),
        IRanges::IRanges(start = 1000, width = 50),
        strand = "+",
        score = 1:5
    )
    res <- filterPeakChromosomes(gr, hg38, "standard")

    expect_equal(as.character(GenomicRanges::seqnames(res$gr)), c("chr1", "chrM", "chr2"))
    expect_equal(res$gr$score, 1:3)

    s <- res$summary
    expect_equal(colnames(s), c("step", "removed", "remaining", "details"))
    expect_equal(s$remaining, c(5L, 5L, 4L, 3L))
    expect_equal(s$removed[3:4], c(1L, 1L))
    expect_match(s$details[2], "1 -> chr1")
    expect_equal(s$details[3], "scaffold_1")
    expect_equal(s$details[4], "chrUn_KI270742v1")

    # "all" keeps the unplaced contig
    expect_length(filterPeakChromosomes(gr, hg38, "all")$gr, 4)

    # No renaming needed: details stay empty
    plain <- filterPeakChromosomes(gr[3], hg38, "standard")
    expect_equal(plain$summary$details[2], "")
})

test_that("formatFilterSummary returns an aligned table", {
    s <- data.frame(
        step = c("Input peaks", "Chromosome not in genome"),
        removed = c(NA, 1200L),
        remaining = c(12480L, 11280L),
        details = c("", "a, b, c, d, e"),
        stringsAsFactors = FALSE
    )
    txt <- formatFilterSummary(s, "Title")
    lines <- strsplit(txt, "\n")[[1]]

    expect_equal(lines[1], "Title")
    expect_length(lines, 4)
    expect_match(txt, "12,480")
    expect_match(txt, "a, b, c \\(\\+2 more\\)")
    expect_equal(addFilterStep(s, "Next", 80L)$remaining[3], 11200L)
})

test_that("motifEnrichment accepts Ensembl names and returns filter summary", {
    skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")
    hg38 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

    peaks <- data.frame(
        chr = c("9", "9", "scaffold_1"),
        start = c(120760581, 120760600, 1),
        end = c(120760630, 120760650, 100),
        strand = c("-", "-", "+")
    )
    result <- suppressMessages(motifEnrichment(
        peaks, hg38,
        K = 5, bkg_iter = 2, bkg_min_dist = 100, bkg_max_dist = 200
    ))
    s <- attr(result, "filter_summary")

    expect_s3_class(result, "data.frame")
    expect_s3_class(s, "data.frame")
    expect_equal(s$remaining[nrow(s)], 2L)
    expect_equal(s$details[s$step == "Chromosome not in genome"], "scaffold_1")
})

test_that("motifEnrichment explains why no peaks remain", {
    skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")

    expect_error(
        suppressMessages(motifEnrichment(
            data.frame(chr = "scaffold_1", start = 1, end = 100),
            "hg38",
            K = 5
        )),
        "do not match genome 'hg38'"
    )
    expect_error(
        suppressMessages(motifEnrichment(
            data.frame(chr = "chr9", start = 120760581, end = 120760630),
            "hg38",
            K = 5, chromosomes = "chr1"
        )),
        "No peaks are on the selected chromosomes"
    )
})
