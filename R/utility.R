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
#' @importFrom Biostrings DNAString BString
#' @importFrom methods is
#' @keywords internal
scrambleDNA <- function(seq) {
  # Robust check for NULL, NA, or empty string content
  # Use any(is.na(seq)) to handle cases where is.na(seq) might return a vector,
  # ensuring the condition for `||` is scalar.
  if (is.null(seq) || any(is.na(seq)) || nchar(as.character(seq)) == 0) {
    return(Biostrings::DNAString())
  }

  # Ensure input is treated as character for splitting and sampling
  # For DNAString (which inherits XString), as.character() is well-defined.
  if (methods::is(seq, "XString")) {
    char_seq <- as.character(seq)
  } else {
    char_seq <- as.character(seq) # For plain character strings
  }

  # The nchar check above should already catch char_seq == ""
  # but as an additional safeguard if seq was e.g. c("", "") and got through.
  # However, the function is documented for a single sequence.
  if (nchar(char_seq) == 0) { # Redundant if the first check worked, but safe.
    return(Biostrings::DNAString())
  }

  nucleotides <- unlist(strsplit(char_seq, split = ""))

  # Ensure there's something to sample
  if (length(nucleotides) == 0) {
    return(Biostrings::DNAString())
  }

  scrambled_nucleotides <- sample(nucleotides)
  scrambled_seq_char <- paste(scrambled_nucleotides, collapse = "")

  return(Biostrings::DNAString(scrambled_seq_char))
}


#' Scramble a Set of DNA Sequences
#'
#' @description Vectorized version of \code{\link{scrambleDNA}} that operates on
#'   an entire DNAStringSet. More efficient than applying scrambleDNA element-wise
#'   because it avoids N individual DNAString constructor calls, constructing one
#'   DNAStringSet from a character vector at the end.
#'
#' @param dna_set A DNAStringSet object.
#'
#' @return A DNAStringSet object with each sequence independently shuffled.
#'   Returns an empty DNAStringSet if input is empty.
#' @importFrom Biostrings DNAStringSet
#' @keywords internal
scrambleDNASet <- function(dna_set) {
  if (length(dna_set) == 0) return(Biostrings::DNAStringSet())

  char_seqs <- as.character(dna_set)
  scrambled <- vapply(char_seqs, function(s) {
    if (is.na(s) || nchar(s) == 0) return("")
    chars <- strsplit(s, "")[[1]]
    paste0(sample(chars), collapse = "")
  }, character(1), USE.NAMES = FALSE)

  Biostrings::DNAStringSet(scrambled)
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
genMotifVar <- function(motif, type = "DNA") {
  # Validate 'type' input
  if (!is.character(type) || length(type) != 1) {
    stop("'type' must be a single character string ('DNA' or 'RNA').")
  }

  type_upper <- toupper(type) # Normalize for case-insensitive comparison

  nucleotides <- switch(type_upper,
    "DNA" = c("A", "C", "G", "T"),
    "RNA" = c("A", "C", "G", "U"),
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
    stop(paste0(
      "Motif '", motif, "' contains invalid characters: '",
      paste(invalid_chars, collapse = "', '"),
      "'. Expected characters for type '", type_upper, "' are: '",
      paste(nucleotides, collapse = "', '"), "'."
    ))
  }

  variants <- character() # Initialize empty vector
  motif_len <- nchar(motif)
  original_chars <- strsplit(motif, "")[[1]]

  # Loop through each position
  for (i in seq_len(motif_len)) {
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
  if (nrow(motif_enrichment) == 0) {
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
    stop(
      "Provided motif '", motif, "' has length ", nchar(motif),
      ", but expected length is ", kmer_size, "."
    )
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
    stop(
      "Unknown species or build identifier: '", species_or_build,
      "'. Supported identifiers include: ", paste(names(genome_map), collapse = ", ")
    )
  }

  # Check if the package is installed
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    stop(
      "The required genome package '", pkg_name, "' is not installed.\n",
      "Please install it using BiocManager, e.g.:\n",
      "if (!requireNamespace(\"BiocManager\", quietly = TRUE))\n",
      "    install.packages(\"BiocManager\")\n",
      "BiocManager::install(\"", pkg_name, "\")"
    )
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

#' @param keep_extra Logical, if TRUE, keep additional metadata columns when converting
#'   data frame to GRanges (default: TRUE).
#'
#' @return A GRanges object representing the peaks.
#' @importFrom GenomicRanges GRanges makeGRangesFromDataFrame seqnames start end strand
#' @importFrom methods is
#' @keywords internal
peakParse <- function(input, standard_chroms_only = TRUE, keep_extra = TRUE) {
  initial_gr <- NULL

  if (methods::is(input, "GRanges")) {
    initial_gr <- input
  } else if (is.data.frame(input)) {
    df <- input
    col_names <- colnames(df)

    chr_col <- NULL
    start_col <- NULL
    end_col <- NULL
    strand_col <- NULL # Optional

    # Try to find chromosome column
    if ("seqnames" %in% col_names) {
      chr_col <- "seqnames"
    } else if ("chr" %in% col_names) {
      chr_col <- "chr"
    } else if ("chromosome" %in% col_names) {
      chr_col <- "chromosome"
    } else {
      stop("Input data frame must contain a chromosome column (e.g., 'chr', 'seqnames').")
    }

    if ("start" %in% col_names) {
      start_col <- "start"
    } else {
      stop("Input data frame must contain a 'start' column.")
    }

    if ("end" %in% col_names) {
      end_col <- "end"
    } else {
      stop("Input data frame must contain an 'end' column.")
    }

    if ("strand" %in% col_names) strand_col <- "strand"

    # message("Attempting to create GRanges from data frame in peakParse...")
    initial_gr <- tryCatch(
      {
        GenomicRanges::makeGRangesFromDataFrame(
          df = df,
          keep.extra.columns = keep_extra,
          seqnames.field = chr_col,
          start.field = start_col,
          end.field = end_col,
          strand.field = if (!is.null(strand_col)) strand_col else character(0),
          starts.in.df.are.0based = FALSE # Assume 1-based starts unless specified
        )
      },
      error = function(e) {
        stop("Failed to create GRanges object from data frame in peakParse.\n",
          "Ensure columns '", chr_col, "', '", start_col, "', '", end_col,
          "' (and optionally '", strand_col, "') are present and valid.\n",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
    # message("GRanges object created successfully from data frame in peakParse.")
  } else {
    stop("Input to peakParse must be a data frame or a GRanges object.")
  }

  if (length(initial_gr) == 0) {
    # message("peakParse received or resulted in an empty GRanges object before standard chromosome filtering.")
    return(initial_gr) # Return empty GRanges if it's already empty
  }

  # Filter for standard chromosomes if requested
  if (standard_chroms_only) {
    original_len <- length(initial_gr)

    # Note: keepStandardChromosomes works best if seqlevelsStyle is already UCSC (e.g. 'chr1', 'chrM').
    # If your input uses 'MT' and BSgenome uses 'chrM', 'MT' might be dropped.
    # You may need to ensure consistent chromosome naming for 'chrM' vs 'MT' *before* this step
    # or add logic here to standardize mitochondrial chromosome name if needed.
    # For example, if "chrMT" exists and "chrM" does not, and style is UCSC-like:
    # current_sl <- GenomeInfoDb::seqlevels(initial_gr)
    # if ("chrMT" %in% current_sl && !("chrM" %in% current_sl) && grepl("chr", current_sl[1])) {
    #   message("Renaming 'chrMT' to 'chrM' in GRanges for standard chromosome filtering.")
    #   GenomeInfoDb::seqlevels(initial_gr)[GenomeInfoDb::seqlevels(initial_gr) == "chrMT"] <- "chrM"
    # }

    filtered_gr <- GenomeInfoDb::keepStandardChromosomes(initial_gr, pruning.mode = "tidy")
    num_removed <- original_len - length(filtered_gr)
    if (num_removed > 0) {
      message(num_removed, " ranges removed by filtering for standard chromosomes in peakParse().")
    }
    initial_gr <- filtered_gr
  }

  # message("peakParse returning GRanges object with ", length(initial_gr), " ranges.")
  return(initial_gr)
}


#' Get Sequences from GRanges
#'
#' @description Extracts DNA sequences for a set of genomic ranges from a
#'   BSgenome object. Optionally extends or trims ranges from their 5'-end
#'   and/or 3'-end in a strand-aware manner, filters by minimum length, and
#'   ensures ranges are valid within chromosome boundaries. Original metadata
#'   columns (`mcols`) and names from the input `granges_obj` are preserved
#'   for the ranges that successfully yield sequences.
#'
#' @param granges_obj A GRanges object.
#' @param genome_obj A BSgenome object.
#' @param extension A numeric vector of length 2: `c(five_prime, three_prime)`.
#'   Positive values extend the range, negative values trim the range.
#'   Extension/trimming is strand-aware:
#'   - For + strand (and * strand): 5'-end is start, 3'-end is end
#'   - For - strand: 5'-end is end, 3'-end is start
#'   Default: `c(0, 0)` (no change).
#' @param min_length Integer, the minimum acceptable length of a range after
#'   any extension/trimming. Ranges shorter than this will be removed.
#'   Default: 1.
#'
#' @return A DNAStringSet object containing the extracted sequences. Names of the
#'   sequences in the set will correspond to the names of the processed GRanges
#'   object (if names were present on the input and preserved). Original metadata
#'   columns are preserved on the intermediate GRanges object used for sequence
#'   extraction. Returns an empty DNAStringSet if input `granges_obj` is empty
#'   or if all ranges are filtered out.
#' @importFrom Biostrings getSeq DNAStringSet
#' @importFrom GenomicRanges GRanges width strand seqnames trim start end start<- end<-
#' @importFrom GenomeInfoDb seqlengths<- seqlevels seqlevelsInUse keepSeqlevels seqlengths
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom methods is
#' @keywords internal
getSequence <- function(granges_obj, genome_obj, extension = c(0, 0), min_length = 1) {
  # Input Validations
  if (!methods::is(granges_obj, "GRanges")) {
    stop("'granges_obj' must be a GRanges object.")
  }
  if (!methods::is(genome_obj, "BSgenome")) {
    stop("'genome_obj' must be a BSgenome object.")
  }
  if (!is.numeric(extension) || length(extension) != 2) {
    stop("'extension' must be a numeric vector of length 2: c(five_prime, three_prime).")
  }
  extension <- as.integer(extension) # Ensure integers
  ext_5prime <- extension[1]
  ext_3prime <- extension[2]

  if (!is.numeric(min_length) || length(min_length) != 1 || min_length < 1 || min_length %% 1 != 0) {
    stop("'min_length' must be a single positive integer.")
  }

  # Handle empty granges_obj
  if (length(granges_obj) == 0) {
    return(Biostrings::DNAStringSet())
  }

  # Store original metadata and add an index for robust re-attachment
  original_mcols <- S4Vectors::mcols(granges_obj)
  original_names <- names(granges_obj)
  original_count <- length(granges_obj)

  # Work with a copy to modify coordinates
  current_gr <- granges_obj
  current_gr$..original_index.. <- seq_along(current_gr) # Temporary index column

  # Get strand information
  strand_info <- as.character(GenomicRanges::strand(current_gr))

  # Identify strand types
  plus_strand_idx <- strand_info %in% c("+", "*")
  minus_strand_idx <- strand_info == "-"

  # Get current coordinates
  new_start <- GenomicRanges::start(current_gr)
  new_end <- GenomicRanges::end(current_gr)

  # Apply strand-aware extension/trimming

  # For + strand (and * strand): 5' = start, 3' = end
  #   5' extension: start = start - ext_5prime (extend upstream)
  #   3' extension: end = end + ext_3prime (extend downstream)
  new_start[plus_strand_idx] <- new_start[plus_strand_idx] - ext_5prime
  new_end[plus_strand_idx] <- new_end[plus_strand_idx] + ext_3prime

  # For - strand: 5' = end, 3' = start (reversed)
  #   5' extension: end = end + ext_5prime (extend upstream on - strand)
  #   3' extension: start = start - ext_3prime (extend downstream on - strand)
  new_end[minus_strand_idx] <- new_end[minus_strand_idx] + ext_5prime
  new_start[minus_strand_idx] <- new_start[minus_strand_idx] - ext_3prime

  # Calculate final widths (before applying to GRanges, to check validity)
  final_widths <- new_end - new_start + 1

  # Filter ranges that would become too short (< min_length)
  valid_idx <- final_widths >= min_length

  num_removed_short <- sum(!valid_idx)
  if (num_removed_short > 0) {
    message(
      num_removed_short, " ranges removed due to final length < ", min_length,
      " after extension/trimming."
    )
  }

  if (!any(valid_idx)) {
    message("No ranges remaining after filtering for minimum length.")
    return(Biostrings::DNAStringSet())
  }

  # Apply changes only to valid ranges
  current_gr <- current_gr[valid_idx]
  new_start <- new_start[valid_idx]
  new_end <- new_end[valid_idx]

  # Update coordinates
  GenomicRanges::start(current_gr) <- new_start
  GenomicRanges::end(current_gr) <- new_end

  # Filter for valid seqlevels & trim to chromosome boundaries
  valid_seqlevels <- intersect(GenomeInfoDb::seqlevels(current_gr), GenomeInfoDb::seqlevels(genome_obj))

  if (length(valid_seqlevels) < length(GenomeInfoDb::seqlevelsInUse(current_gr))) {
    original_gr_count_sl <- length(current_gr)
    current_gr <- GenomeInfoDb::keepSeqlevels(current_gr, valid_seqlevels, pruning.mode = "tidy")
    num_removed_sl <- original_gr_count_sl - length(current_gr)
    if (num_removed_sl > 0) {
      message(num_removed_sl, " ranges removed due to seqlevels not present in the genome object.")
    }
  }

  if (length(current_gr) == 0) {
    message("No ranges remaining after filtering for seqlevels present in genome object.")
    return(Biostrings::DNAStringSet())
  }

  genome_seqlengths_subset <- GenomeInfoDb::seqlengths(genome_obj)[GenomeInfoDb::seqlevels(current_gr)]
  genome_seqlengths_subset <- genome_seqlengths_subset[!is.na(genome_seqlengths_subset)]
  common_levels_for_seqlengths <- intersect(names(genome_seqlengths_subset), GenomeInfoDb::seqlevels(current_gr))
  GenomeInfoDb::seqlengths(current_gr)[common_levels_for_seqlengths] <- genome_seqlengths_subset[common_levels_for_seqlengths]

  # Trim to chromosome boundaries (handles negative starts and exceeding ends)
  trimmed_gr <- GenomicRanges::trim(current_gr)

  # After trimming, check min_length again (coordinates may have been clipped)
  widths_after_trim <- GenomicRanges::width(trimmed_gr)
  keep_after_trim <- widths_after_trim >= min_length

  if (any(!keep_after_trim)) {
    num_removed_after_trim <- sum(!keep_after_trim)
    message(
      num_removed_after_trim, " ranges removed after chromosome boundary trimming ",
      "resulted in length < ", min_length, "."
    )
    trimmed_gr <- trimmed_gr[keep_after_trim]
  }

  if (length(trimmed_gr) == 0) {
    message("No valid ranges remaining after trimming to chromosome boundaries.")
    return(Biostrings::DNAStringSet())
  }

  # Metadata re-attachment using original index
  final_original_indices <- trimmed_gr$..original_index..

  if (!is.null(final_original_indices)) {
    if (!is.null(original_mcols) && ncol(original_mcols) > 0) {
      S4Vectors::mcols(trimmed_gr) <- original_mcols[final_original_indices, , drop = FALSE]
    } else {
      S4Vectors::mcols(trimmed_gr) <- NULL
    }
    if (!is.null(original_names)) {
      names(trimmed_gr) <- original_names[final_original_indices]
    } else {
      names(trimmed_gr) <- NULL
    }
    trimmed_gr$..original_index.. <- NULL # Remove temporary index column
  } else if (length(trimmed_gr) > 0) {
    warning("Could not re-attach original metadata due to missing internal index. Metadata might be lost.")
    S4Vectors::mcols(trimmed_gr) <- NULL
    names(trimmed_gr) <- NULL
  }

  # Get sequences
  sequences <- Biostrings::getSeq(genome_obj, trimmed_gr)

  return(sequences)
}


#' Count K-mers in Sequences
#'
#' @description Counts occurrences of all possible K-mers of a given length
#'   within a DNAStringSet. This function is case-sensitive for motifs once
#'   generated (all K-mers are generated in uppercase).
#'
#' @param sequences A DNAStringSet object, typically the output from
#'   \code{\link{getSequence}}.
#' @param K Integer, the K-mer size (length of motifs). Must be a single
#'   positive integer.
#' @param type Character string specifying the nucleic acid type for K-mer
#'   generation. Accepted values are "DNA" (using A,C,G,T) or "RNA"
#'   (using A,C,G,U). The input is case-insensitive. Default: "DNA".
#'
#' @return A data frame with two columns: 'MOTIF' (character, all possible
#'   K-mers in uppercase) and 'COUNT' (integer). 'COUNT' is the total
#'   count of each K-mer across all input sequences. Returns a data frame
#'   with all possible K-mers and 0 counts if input sequences are empty or
#'   if all sequences are shorter than K.
#'
#' @importFrom Biostrings DNAStringSet PDict vcountPDict
#' @importFrom S4Vectors elementNROWS
#' @importFrom methods is
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Assuming 'seqs' is a DNAStringSet object and 'getSequence' is available
#' # library(Biostrings)
#' # seqs_dna <- DNAStringSet(c("ATGCGATGC", "GGCCTTAA", "TTT"))
#' # countKmers(seqs_dna, K = 3, type = "DNA")
#'
#' # seqs_rna <- DNAStringSet(c("AUGCGAUGC", "GGCCUUAA", "UUU"))
#' # countKmers(seqs_rna, K = 3, type = "RNA")
#'
#' # Empty input
#' # countKmers(DNAStringSet(), K = 3)
#'
#' # Sequences shorter than K
#' # countKmers(DNAStringSet(c("A", "CG")), K = 3)
#' }
countKmers <- function(sequences, K, type = "DNA") {
  # Input Validation for 'sequences'
  if (!methods::is(sequences, "DNAStringSet")) {
    stop("'sequences' must be a DNAStringSet object.")
  }

  # Input Validation for 'K'
  if (!is.numeric(K) || length(K) != 1 || K <= 0 || K %% 1 != 0) {
    stop("'K' must be a single positive integer.")
  }
  if (K > 12) { # Practical warning for very large K
    warning("'K' is ", K, ", generating all K-mer combinations can be memory/time intensive.")
  }

  # Input Validation and setup for 'type'
  if (!is.character(type) || length(type) != 1) {
    stop("'type' must be a single character string ('DNA' or 'RNA').")
  }
  type_upper <- toupper(type)
  if (!(type_upper %in% c("DNA", "RNA"))) {
    stop("Invalid 'type': must be 'DNA' or 'RNA'. Received: '", type, "'")
  }

  # oligonucleotideFrequency returns a matrix with length(sequences) rows.
  # colSums sums counts over all sequences in the DNAStringSet.
  # This returns counts for all possible K-mers of length K in the alphabet.
  kmer_counts <- colSums(Biostrings::oligonucleotideFrequency(sequences, width = K))

  # Create and return results data frame
  results_df <- data.frame(
    MOTIF = names(kmer_counts),
    COUNT = as.integer(kmer_counts),
    stringsAsFactors = FALSE
  )

  # If type is RNA, replace T with U in motifs
  if (type_upper == "RNA") {
    results_df$MOTIF <- gsub("T", "U", results_df$MOTIF)
  }

  return(results_df)
}


#' Generate a Single Set of Background Sequences
#'
#' @description Creates one set of background sequences by randomly
#'   shifting original peak locations, removing overlaps with original peaks,
#'   and extracting sequences. Optionally scrambles the sequences.
#'
#' @param peak_gr A GRanges object representing the original genomic peaks.
#' @param genome_obj A BSgenome object from which sequences will be extracted.
#' @param K Integer, the K-mer length. This is used to set a minimum length
#'   for sequences extracted by the internal call to `getSequence`.
#' @param bkg_min_dist Integer, the minimum absolute distance (in base pairs) by which
#'   peaks should be shifted. Passed as `X` to `genRNDist()`.
#' @param bkg_max_dist Integer, the maximum absolute distance (in base pairs)
#'   for shifting peaks. Passed as `max_dist` to `genRNDist()`.
#' @param scramble Logical, if TRUE (default), scrambles the background sequences
#'   after extraction. Set to FALSE for faster execution when scrambling is not needed.
#'
#' @return A DNAStringSet object of background sequences (scrambled or not).
#'   Returns an empty DNAStringSet if no valid background regions/sequences
#'   could be generated.
#'
#' @importFrom GenomicRanges GRanges shift findOverlaps width
#' @importFrom Biostrings DNAStringSet BString
#' @importFrom S4Vectors queryHits
#' @importFrom methods is
#'
#' @keywords internal
generateBkgSet <- function(peak_gr, genome_obj, K,
                           bkg_min_dist, bkg_max_dist, scramble = TRUE) {
  # Input Validation
  if (!methods::is(peak_gr, "GRanges") || length(peak_gr) == 0) {
    stop("'peak_gr' must be a non-empty GRanges object in generateBkgSet.")
  }
  if (!methods::is(genome_obj, "BSgenome")) {
    stop("'genome_obj' must be a BSgenome object in generateBkgSet.")
  }
  if (!is.numeric(K) || length(K) != 1 || K <= 0 || K %% 1 != 0) {
    stop("'K' must be a single positive integer in generateBkgSet.")
  }
  # bkg_min_dist and bkg_max_dist validation is handled by genRNDist

  num_peaks <- length(peak_gr)

  # Generate random shift distances
  rand_distances <- genRNDist(N = num_peaks, X = bkg_min_dist, max_dist = bkg_max_dist)

  # Shift original peaks
  # Suppress warnings for out-of-bounds shifts; getSequence will trim.
  shifted_gr <- suppressWarnings(GenomicRanges::shift(peak_gr, shift = rand_distances))

  # Remove shifted regions that overlap with ANY original peak
  # type = "any" means any kind of overlap.
  # select = "all" gives all pairs. We need queryHits to index shifted_gr.
  # ignore.strand = TRUE because overlaps are purely positional for background exclusion.
  overlaps_with_original <- GenomicRanges::findOverlaps(
    shifted_gr,
    peak_gr,
    type = "any",
    select = "all",
    ignore.strand = TRUE
  )

  indices_to_remove <- S4Vectors::queryHits(overlaps_with_original)
  if (length(indices_to_remove) > 0) {
    shifted_gr_filtered <- shifted_gr[-indices_to_remove]
  } else {
    shifted_gr_filtered <- shifted_gr
  }

  if (length(shifted_gr_filtered) == 0) {
    # No valid non-overlapping regions generated
    return(Biostrings::DNAStringSet())
  }

  # Get sequences for the valid background regions
  # No additional extension (extension = c(0, 0)).
  # Sequences must be at least K long for countKmers to work.
  # Suppress messages to avoid redundant per-iteration output.
  background_seqs_gr <- suppressMessages(getSequence(
    granges_obj = shifted_gr_filtered,
    genome_obj = genome_obj,
    extension = c(0, 0),
    min_length = K # Ensure sequences are viable for K-mer counting
  ))

  if (length(background_seqs_gr) == 0) {
    # No valid sequences obtained (e.g., all trimmed to < K length)
    return(Biostrings::DNAStringSet())
  }

  # Conditionally scramble sequences
  if (scramble) {
    # scrambleDNA works on a single sequence. Use lapply for DNAStringSet.
    # as.list converts DNAStringSet to a list of individual DNAString objects.
    scrambled_seq_list <- lapply(as.list(background_seqs_gr), scrambleDNA)
    final_seqs <- Biostrings::DNAStringSet(scrambled_seq_list)
  } else {
    # Use unscrambled shifted sequences
    final_seqs <- background_seqs_gr
  }

  # Filter out any truly empty sequences that might result if scrambleDNA returned empty
  # (Our scrambleDNA returns empty DNAString for empty/NA input, which then results in 0-width here)
  final_seqs <- final_seqs[Biostrings::width(final_seqs) > 0]

  if (length(final_seqs) == 0) {
    return(Biostrings::DNAStringSet())
  }

  return(final_seqs)
}


#' Generate Batched Background Sequences Across Multiple Iterations
#'
#' @description Batched version of \code{\link{generateBkgSet}} that processes
#'   multiple iterations at once, reducing per-call overhead for GRanges shifting,
#'   overlap removal, and sequence extraction.
#'
#' @param peak_gr A GRanges object representing the original genomic peaks.
#' @param genome_obj A BSgenome object.
#' @param min_seq_length Integer, minimum length for extracted sequences.
#' @param bkg_min_dist Integer, minimum absolute shift distance.
#' @param bkg_max_dist Integer, maximum absolute shift distance.
#' @param scramble Logical, whether to scramble sequences.
#' @param n_iter Integer, number of iterations to batch.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{sequences}{A DNAStringSet of all background sequences.}
#'     \item{iter_tags}{Integer vector of iteration indices (1 to n_iter),
#'       same length as sequences, indicating which iteration each sequence
#'       belongs to.}
#'   }
#' @importFrom GenomicRanges shift findOverlaps trim width
#' @importFrom GenomeInfoDb seqlevels seqlevelsInUse seqlengths seqlengths<- seqnames keepSeqlevels
#' @importFrom Biostrings getSeq DNAStringSet width
#' @importFrom S4Vectors queryHits
#' @importFrom methods is
#' @keywords internal
generateBkgSetBatched <- function(peak_gr, genome_obj, min_seq_length,
                                   bkg_min_dist, bkg_max_dist,
                                   scramble, n_iter) {
  empty_result <- list(sequences = Biostrings::DNAStringSet(), iter_tags = integer(0))
  num_peaks <- length(peak_gr)
  total_regions <- num_peaks * n_iter

  # 1. Replicate peaks and tag with iteration index
  batched_gr <- rep(peak_gr, times = n_iter)
  iter_tags <- rep(seq_len(n_iter), each = num_peaks)

  # 2. Generate all random distances and shift at once
  rand_distances <- genRNDist(N = total_regions, X = bkg_min_dist, max_dist = bkg_max_dist)
  shifted_gr <- suppressWarnings(GenomicRanges::shift(batched_gr, shift = rand_distances))

  # 3. Remove shifted regions overlapping with original peaks
  overlaps <- GenomicRanges::findOverlaps(
    shifted_gr, peak_gr,
    type = "any", select = "all", ignore.strand = TRUE
  )
  indices_to_remove <- unique(S4Vectors::queryHits(overlaps))
  if (length(indices_to_remove) > 0) {
    keep_idx <- setdiff(seq_along(shifted_gr), indices_to_remove)
    shifted_gr <- shifted_gr[keep_idx]
    iter_tags <- iter_tags[keep_idx]
  }
  if (length(shifted_gr) == 0) return(empty_result)

  # 4. Filter to valid seqlevels present in genome
  valid_seqlevels <- intersect(
    GenomeInfoDb::seqlevelsInUse(shifted_gr),
    GenomeInfoDb::seqlevels(genome_obj)
  )
  if (length(valid_seqlevels) < length(GenomeInfoDb::seqlevelsInUse(shifted_gr))) {
    keep_sl <- as.character(GenomeInfoDb::seqnames(shifted_gr)) %in% valid_seqlevels
    shifted_gr <- shifted_gr[keep_sl]
    iter_tags <- iter_tags[keep_sl]
  }
  if (length(shifted_gr) == 0) return(empty_result)

  # Set seqlengths for trim
  shifted_gr <- GenomeInfoDb::keepSeqlevels(shifted_gr, valid_seqlevels, pruning.mode = "tidy")
  genome_sl <- GenomeInfoDb::seqlengths(genome_obj)[GenomeInfoDb::seqlevels(shifted_gr)]
  genome_sl <- genome_sl[!is.na(genome_sl)]
  common_sl <- intersect(names(genome_sl), GenomeInfoDb::seqlevels(shifted_gr))
  if (length(common_sl) > 0) {
    GenomeInfoDb::seqlengths(shifted_gr)[common_sl] <- genome_sl[common_sl]
  }

  # 5. Trim to chromosome boundaries and filter by min length
  shifted_gr <- GenomicRanges::trim(shifted_gr)
  valid_width <- GenomicRanges::width(shifted_gr) >= min_seq_length
  if (any(!valid_width)) {
    shifted_gr <- shifted_gr[valid_width]
    iter_tags <- iter_tags[valid_width]
  }
  if (length(shifted_gr) == 0) return(empty_result)

  # 6. Extract sequences directly (no metadata handling needed for background)
  sequences <- Biostrings::getSeq(genome_obj, shifted_gr)

  # 7. Scramble if requested
  if (scramble) {
    sequences <- scrambleDNASet(sequences)
    valid_after <- Biostrings::width(sequences) > 0
    if (any(!valid_after)) {
      sequences <- sequences[valid_after]
      iter_tags <- iter_tags[valid_after]
    }
  }

  list(sequences = sequences, iter_tags = iter_tags)
}


#' Normalize a Vector of Scores
#'
#' @description Applies a specified normalization method to a numeric vector of scores.
#'
#' @param scores A numeric vector.
#' @param method Character string (case-insensitive), the normalization method to apply.
#'   Supported methods: "min_max" (default), "z_score", "log2", "none".
#' @param pseudocount Numeric, pseudocount added before log2 transformation
#'   (relevant only if `method = "log2"`). Default: 1.
#' @param ... Additional arguments passed to specific normalization methods.
#'   For "min_max", these are `a` and `b` for the target range (defaults to 0 and 1
#'   respectively if not provided, via `minmaxNorm` defaults).
#'
#' @return A numeric vector of normalized scores.
#'
#' @importFrom stats sd
#' @keywords internal # Or could be exported if generally useful
#'
#' @examples
#' \dontrun{
#' test_scores <- c(10, 20, 30, 40, 50, NA)
#' normalizeScores(test_scores, method = "min_max")
#' normalizeScores(test_scores, method = "min_max", a = 0, b = 10)
#' normalizeScores(test_scores, method = "z_score")
#' normalizeScores(c(0, 1, 3, 7, 15, NA), method = "log2") # Uses pseudocount = 1
#' normalizeScores(c(0, 1, 3, 7, 15), method = "log2", pseudocount = 0.1)
#' normalizeScores(test_scores, method = "none") # Returns original scores
#' }
normalizeScores <- function(scores, method = "min_max", pseudocount = 1, ...) {
  # Validate scores input
  if (!is.numeric(scores)) {
    stop("'scores' must be a numeric vector.")
  }
  if (length(scores) == 0) {
    return(numeric(0)) # Return empty numeric if input is empty
  }

  # Validate method input
  if (!is.character(method) || length(method) != 1) {
    stop("'method' must be a single character string.")
  }

  method_upper <- toupper(method)
  additional_args <- list(...)

  normalized_scores <- switch(method_upper,
    "MIN_MAX" = {
      # Pass additional_args (like a, b) to minmaxNorm
      # minmaxNorm has defaults a=0, b=1 if not provided in ...
      do.call(minmaxNorm, c(list(x = scores), additional_args))
    },
    "Z_SCORE" = {
      if (length(additional_args) > 0) {
        warning("Additional arguments in '...' are ignored for Z-score normalization.")
      }
      mean_s <- mean(scores, na.rm = TRUE)
      sd_s <- stats::sd(scores, na.rm = TRUE)
      if (is.na(sd_s)) { # Happens if too few non-NA values
        warning("Cannot compute Z-scores (standard deviation is NA). Returning NAs.")
        return(rep(NA_real_, length(scores)))
      }
      if (sd_s == 0) {
        # If all non-NA values are the same, Z-score is typically 0
        # or could be NA/NaN depending on convention. Returning 0 for non-NA.
        warning("Standard deviation is zero. Z-scores for non-NA values will be 0.")
        result <- rep(NA_real_, length(scores))
        result[!is.na(scores)] <- 0
        return(result)
      }
      (scores - mean_s) / sd_s
    },
    "LOG2" = {
      if (length(additional_args) > 0) {
        warning("Additional arguments in '...' are ignored for log2 normalization (except 'pseudocount').")
      }
      if (!is.numeric(pseudocount) || length(pseudocount) != 1) {
        warning("'pseudocount' is not a single number, using default of 1.")
        pseudocount <- 1
      }

      scores_plus_pseudo <- scores + pseudocount
      if (any(scores_plus_pseudo <= 0, na.rm = TRUE)) {
        warning("Some scores + pseudocount are <= 0. log2 will produce -Inf or NaN for these values.")
      }
      log2(scores_plus_pseudo)
    },
    "NONE" = {
      if (length(additional_args) > 0) {
        warning("Additional arguments in '...' are ignored when method is 'none'.")
      }
      scores
    },
    # Default case for switch if method_upper is not matched
    stop(
      "Unsupported normalization method: '", method,
      "'. Supported methods are 'min_max', 'z_score', 'log2', 'none'."
    )
  )
  return(normalized_scores)
}
#-------------------------------------------------------------------------------
