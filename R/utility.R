# Utility and helper functions used across multiple modules

#' Min-Max Normalization
#'
#' @description Scale a numeric vector to a specified range.
#'   Handles cases where min and max are equal.
#'
#' @param x Numeric vector.
#' @param a Numeric, the minimum value of the target range (default: 0).
#' @param b Numeric, the maximum value of the target range (default: 1).
#'
#' @return Numeric vector scaled to the range [a, b]. NA values in input
#'   will result in NA values in output.
#' @keywords internal
minmaxNorm <- function(x, a = 0, b = 1) {
  # Ensure x is numeric
  if (!is.numeric(x)) {
    stop("Input 'x' must be a numeric vector.")
  }
  # Get min and max, removing NA values for calculation
  min_x <- min(x, na.rm = TRUE)
  max_x <- max(x, na.rm = TRUE)

  # Check for edge cases: all NA, or constant value
  if (all(is.na(x))) {
    return(x) # Return vector of NAs
  }
  if (min_x == max_x) {
    # If all non-NA values are the same, return a vector where non-NA values
    # are the midpoint of the target range [a, b]
    # Or could return 'a' or 'b' - midpoint seems reasonable.
    mid_point <- (a + b) / 2
    result <- rep(NA_real_, length(x))
    result[!is.na(x)] <- mid_point
    return(result)
  }

  # Perform scaling
  scaled_x <- a + (x - min_x) * (b - a) / (max_x - min_x)
  return(scaled_x)
}


#' Scramble DNA Sequence
#'
#' @description Randomly shuffles nucleotides within a sequence. Input can be
#'   character string or DNAString. Output is always DNAString.
#'
#' @param seq A DNAString or character string representing a sequence.
#'
#' @return A DNAString object with shuffled nucleotides. Returns empty DNAString
#'   if input is empty or NA.
#' @importFrom Biostrings DNAString BString subject
#' @importFrom methods is
#' @keywords internal
scrambleDNA <- function(seq) {
  # Handle NULL, NA, or empty input
  if (is.null(seq) || is.na(seq) || nchar(as.character(seq)) == 0) {
    return(Biostrings::DNAString())
  }

  # Ensure input is treated as character for splitting and sampling
  if (methods::is(seq, "BString")) {
    # Extract sequence if it's a Biostrings object (like DNAString)
    char_seq <- as.character(Biostrings::subject(seq))
  } else {
    char_seq <- as.character(seq)
  }

  # Split, sample, paste, and return as DNAString
  nucleotides <- unlist(strsplit(char_seq, split = ""))
  scrambled_nucleotides <- sample(nucleotides)
  scrambled_seq_char <- paste(scrambled_nucleotides, collapse = "")

  return(Biostrings::DNAString(scrambled_seq_char))
}


#' Generate Random Genomic Distances
#'
#' @description Generates N random distances, ensuring a minimum separation from zero.
#' Used for shifting background regions. Uses runif for sampling.
#'
#' @param N Integer, the number of distances to generate. Must be positive.
#' @param X Integer, the minimum absolute distance from zero. Must be non-negative (default: 500).
#' @param max_dist Integer, the maximum absolute distance for sampling (default: 1000).
#'
#' @return A numeric vector of N random distances.
#' @importFrom stats runif
#' @keywords internal
genRNDist <- function(N, X = 500, max_dist = 1000) {
  if (!is.numeric(N) || N <= 0 || N %% 1 != 0) {
    stop("'N' must be a positive integer.")
  }
  if (!is.numeric(X) || X < 0) {
    stop("'X' must be a non-negative number.")
  }
  if (!is.numeric(max_dist) || max_dist <= X) {
    stop("'max_dist' must be a number greater than 'X'.")
  }

  # Ensure N is even for splitting, adjust if necessary (e.g., floor/ceiling)
  n_neg <- floor(N / 2)
  n_pos <- ceiling(N / 2)

  distances_negative <- stats::runif(n_neg, min = -max_dist, max = -X)
  distances_positive <- stats::runif(n_pos, min = X, max = max_dist)

  distances <- sample(c(distances_negative, distances_positive)) # Shuffle results
  distances <- round(distances, 0) # Round to integer
  return(distances)
}


#' Generate Single Nucleotide Variants for a Motif
#'
#' @description Create all possible single nucleotide variants.
#'
#' @param motif Character string, the reference motif. Should only contain valid nucleotides.
#' @param type Character string specifying the nucleic acid type for the inputmotifand the generated variants. Accepted values are "DNA" (A,C,G,T) or "RNA" (A,C,G,U) (default: 'DNA').
#'
#' @return A character vector of all unique SNV motifs (excluding the original motif).
#' @keywords internal
genMotifVar <- function(motif, type = 'DNA') {
  # Validate 'type' input
  if (!is.character(type) || length(type) != 1) {
    stop("'type' must be a single character string ('DNA' or 'RNA').")
  }

  type_upper <- toupper(type) # Normalize for case-insensitive comparison

  nucleotides <- switch(type_upper,
                        "DNA" = c('A', 'C', 'G', 'T'),
                        "RNA" = c('A', 'C', 'G', 'U'),
                        # Default case if type_upper matches neither:
                        stop("Invalid 'type': must be 'DNA' or 'RNA'. Received: '", type, "'")
  )

  # Validate 'motif' input (basic checks)
  if (!is.character(motif) || length(motif) != 1 || nchar(motif) == 0) {
    stop("'motif' must be a non-empty character string.")
  }

  # Normalize motif to uppercase for consistent character checking and variant generation
  motif_upper <- toupper(motif)
  original_chars <- strsplit(motif_upper, "")[[1]]

  # Check if all characters in motif are valid nucleotides
  are_chars_valid <- all(original_chars %in% nucleotides)

  if (!are_chars_valid) {
    invalid_chars <- unique(original_chars[!(original_chars %in% nucleotides)])
    stop(paste0("Motif '", motif, "' contains invalid characters: '",
                paste(invalid_chars, collapse = "', '"),
                "'. Expected characters for type '", type_upper, "' are: '",
                paste(nucleotides, collapse = "', '"), "'."))
  }

  variants <- character() # Initialize empty vector
  motif_len <- nchar(motif)
  original_chars <- strsplit(motif, "")[[1]]

  # Loop through each position
  for (i in 1:motif_len) {
    # Loop through each possible nucleotide
    for (nuc in nucleotides) {
      # If the nucleotide is different from the original at this position
      if (original_chars[i] != nuc) {
        # Construct the new variant
        variant_chars <- original_chars
        variant_chars[i] <- nuc
        variants <- c(variants, paste0(variant_chars, collapse = ""))
      }
    }
  }
  # Return only unique variants (though the loop structure avoids duplicates here)
  # Also explicitly remove the original motif if it somehow got generated (shouldn't happen with if check)
  variants <- unique(variants)
  variants <- variants[variants != motif]
  return(variants)
}


#' Find Top Motif by Score
#'
#' @description Identifies the motif with the highest score in an enrichment data frame.
#'   Handles ties by returning the first one encountered.
#'
#' @param motif_enrichment Data frame with 'MOTIF' and 'Score' (or 'Enrichment') columns.
#'                         Assumes higher scores are better.
#'
#' @return Character string, the top-scoring motif.
#' @keywords internal
findTopMotif <- function(motif_enrichment) {
  # Basic input validation
  if (!is.data.frame(motif_enrichment)) {
    stop("'motif_enrichment' must be a data frame.")
  }
  required_cols <- c("MOTIF")
  score_col <- NULL
  if ("Score" %in% colnames(motif_enrichment)) {
    score_col <- "Score"
  } else if ("Enrichment" %in% colnames(motif_enrichment)) {
    score_col <- "Enrichment"
  } else {
    stop("Data frame must contain a 'Score' or 'Enrichment' column.")
  }
  required_cols <- c(required_cols, score_col)

  if (!all(required_cols %in% colnames(motif_enrichment))) {
    stop("Data frame must contain 'MOTIF' and either 'Score' or 'Enrichment' columns.")
  }
  if(nrow(motif_enrichment) == 0) {
    stop("Input data frame 'motif_enrichment' is empty.")
  }

  # Find the index of the maximum score
  # which.max returns the index of the first maximum in case of ties
  top_index <- which.max(motif_enrichment[[score_col]])

  # Return the corresponding MOTIF
  # Use top_index[1] just in case which.max somehow returns multiple (shouldn't)
  top_motif <- motif_enrichment$MOTIF[top_index[1]]

  return(top_motif)
}


#' Validate Motif Input
#'
#' @description Checks if a user-provided motif is valid (correct length, exists
#'   within a provided set of available motifs). Stops execution if invalid.
#'
#' @param motif Character string, the motif to validate.
#' @param kmer_size Integer, the expected K-mer length.
#' @param available_motifs Character vector of all motifs present in the data.
#'
#' @return Returns `TRUE` invisibly if valid, otherwise stops with an error.
#' @keywords internal
valInputMotif <- function(motif, kmer_size, available_motifs) {
  # Check motif is a single string
  if (!is.character(motif) || length(motif) != 1) {
    stop("'motif' must be a single character string.")
  }
  # Check length
  if (nchar(motif) != kmer_size) {
    stop("Provided motif '", motif, "' has length ", nchar(motif),
         ", but expected length is ", kmer_size, ".")
  }
  # Check if motif exists in the available set
  if (!(motif %in% available_motifs)) {
    stop("Provided motif '", motif, "' not found in the available motifs dataset.")
  }
  # If all checks pass:
  invisible(TRUE)
}


#' Select and Load BSgenome Object
#'
#' @description Loads the appropriate BSgenome package based on a species or
#' build identifier. Checks if the package is installed.
#'
#' @param species_or_build Character string (e.g., "hg38", "human", "mm10", "mouse").
#'   Currently supports hg19, hg38, mm9, mm10.
#'
#' @return A BSgenome object.
#' @importFrom BSgenome getBSgenome
#' @importFrom methods is
#' @keywords internal
selectGenome <- function(species_or_build) {
  # Normalize input to lower case
  lookup <- tolower(species_or_build)

  # Simple mapping (extend this as needed)
  genome_map <- list(
    "hg19" = "BSgenome.Hsapiens.UCSC.hg19",
    "grch37" = "BSgenome.Hsapiens.UCSC.hg19",
    "human" = "BSgenome.Hsapiens.UCSC.hg38", # Default human to hg38
    "hg38" = "BSgenome.Hsapiens.UCSC.hg38",
    "grch38" = "BSgenome.Hsapiens.UCSC.hg38",
    "mm9" = "BSgenome.Mmusculus.UCSC.mm9",
    "mouse" = "BSgenome.Mmusculus.UCSC.mm10", # Default mouse to mm10
    "mm10" = "BSgenome.Mmusculus.UCSC.mm10",
    "grcm38" = "BSgenome.Mmusculus.UCSC.mm10"
    # Add more mappings here if needed
  )

  pkg_name <- genome_map[[lookup]]

  if (is.null(pkg_name)) {
    stop("Unknown species or build identifier: '", species_or_build,
         "'. Supported identifiers include: ", paste(names(genome_map), collapse=", "))
  }

  # Check if the package is installed
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    stop("The required genome package '", pkg_name, "' is not installed.\n",
         "Please install it using BiocManager, e.g.:\n",
         "if (!requireNamespace(\"BiocManager\", quietly = TRUE))\n",
         "    install.packages(\"BiocManager\")\n",
         "BiocManager::install(\"", pkg_name, "\")")
  }

  # Load and return the genome object
  genome <- BSgenome::getBSgenome(pkg_name)
  return(genome)
}


#' Parse Peak Input Data
#'
#' @description Validates and converts peak input (currently data frame or GRanges)
#' into a standardized GRanges object. Basic column checking included.
#'
#' @param input Can be a data frame (with chr, start, end columns), or a GRanges object.
#'   File path input is not yet supported by this helper.
#' @param required_cols Character vector of essential column names expected if input is a data frame
#'   (default: c("chr", "start", "end")). Common alternatives like "seqnames" are handled.
#' @param keep_extra Logical, if TRUE, keep additional metadata columns when converting
#'   data frame to GRanges (default: TRUE).
#'
#' @return A GRanges object representing the peaks.
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame seqnames start end strand
#' @importFrom methods is
#' @keywords internal
peakParse <- function(input, required_cols = c("chr", "start", "end"), keep_extra = TRUE) {

  if (methods::is(input, "GRanges")) {
    # If already GRanges, maybe perform minimal validation if needed?
    # For now, assume it's valid if it's the right class.
    return(input)

  } else if (is.data.frame(input)) {
    df <- input
    col_names <- colnames(df)

    # Check for required columns - flexible naming (chr/seqnames, start, end)
    chr_col <- NULL
    start_col <- NULL
    end_col <- NULL
    strand_col <- NULL # Optional

    if ("seqnames" %in% col_names) chr_col <- "seqnames"
    else if ("chr" %in% col_names) chr_col <- "chr"
    else if ("chromosome" %in% col_names) chr_col <- "chromosome"
    else stop("Input data frame must contain a chromosome column (e.g., 'chr', 'seqnames').")

    if ("start" %in% col_names) start_col <- "start"
    else stop("Input data frame must contain a 'start' column.")

    if ("end" %in% col_names) end_col <- "end"
    else stop("Input data frame must contain an 'end' column.")

    # Check for optional strand column
    if ("strand" %in% col_names) strand_col <- "strand"


    # Create GRanges
    message("Attempting to create GRanges from data frame...")
    gr <- tryCatch({
      GenomicRanges::makeGRangesFromDataFrame(
        df = df,
        keep.extra.columns = keep_extra,
        seqnames.field = chr_col,
        start.field = start_col,
        end.field = end_col,
        strand.field = if(!is.null(strand_col)) strand_col else character(0) # Handle optional strand
        # ignore.strand = is.null(strand_col) # Handled by strand.field=character(0)
      )
    }, error = function(e) {
      stop("Failed to create GRanges object from data frame.\n",
           "Please ensure columns '", chr_col, "', '", start_col, "', '", end_col,
           "' (and optionally '", strand_col, "') are present and contain valid genomic coordinates.\n",
           "Original error: ", e$message)
    })
    message("GRanges object created successfully.")
    return(gr)

  } else {
    stop("Input must be a data frame or a GRanges object. File path input not yet supported by this helper.")
  }
}


#' Get Sequences from GRanges
#'
#' @description Extracts DNA sequences for a set of genomic ranges from a BSgenome object.
#'
#' @param granges_obj A GRanges object.
#' @param genome_obj A BSgenome object.
#'
#' @return A DNAStringSet object containing the sequences.
#' @importFrom Biostrings getSeq DNAStringSet
#' @importFrom GenomicRanges GRanges
#' @keywords internal
getSequence <- function(granges_obj, genome_obj) {
  # 1. Use Biostrings::getSeq(genome_obj, granges_obj)
  # 2. Return DNAStringSet
  {}
}


#' Count K-mers in Sequences
#'
#' @description Counts occurrences of all possible K-mers within a DNAStringSet.
#'
#' @param sequences A DNAStringSet object.
#' @param K Integer, the K-mer size.
#' @param nucleotides Character vector of nucleotides (default: DNA).
#'
#' @return A data frame with 'MOTIF' and 'COUNT' columns.
#' @importFrom Biostrings DNAStringSet PDict vcountPDict oligonucleotideFrequency
#' @keywords internal
countKmers <- function(sequences, K, nucleotides = c("A", "C", "G", "T")) {
  # 1. Generate all possible K-mer strings
  # 2. Create a PDict object from the K-mer strings
  # 3. Use Biostrings::vcountPDict(pdict, sequences, collapse=1) to count
  # 4. Alternatively, use Biostrings::oligonucleotideFrequency(sequences, K, step=1, as.prob=FALSE)
  #    and sum across sequences if needed. vcountPDict is usually better here.
  # 5. Format results into a data frame ('MOTIF', 'COUNT')
  # 6. Return data frame
  {}
}


#-------------------------------------------------------------------------------
