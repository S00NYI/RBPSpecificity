# Functions related to calculating Mutational Sensitivity (MS)

#' Calculate Mutational Sensitivity (MS)
#'
#' @description Calculates the Mutational Sensitivity (MS) for a target motif.
#' MS reflects how much the enrichment score changes upon single nucleotide
#' variations (SNVs) in the motif.
#'
#' @param motif_enrichment A data frame containing motif enrichment scores
#'   (output from `motifEnrichment()`). Requires 'MOTIF' and 'Score' columns.
#' @param motif Character string, the reference motif for which to calculate MS.
#'   If NULL or "top", the motif with the highest score is used (default: NULL).
#' @param output_type Character string, specifying the output format:
#'   "matrix" returns a matrix of sensitivity scores for each SNV (Nucleotide x Position),
#'   "number" returns a single average sensitivity score across all SNVs (default: "matrix").
#' @param sensitivity_method Character string defining how sensitivity is calculated
#'   (e.g., "1_minus_norm_score", "score_diff", "log2_ratio") (default: "1_minus_norm_score").
#'
#' @return Depending on `output_type`:
#'   - A matrix with nucleotides as rows and positions as columns, containing sensitivity scores.
#'   - A single numeric value representing the average sensitivity score.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming enrichment_results is output from motifEnrichment()
#' ms_matrix <- returnMS(enrichment_results, motif = "UGUGU", output_type = "matrix")
#' ms_average <- returnMS(enrichment_results, motif = "UGUGU", output_type = "number")
#' }
returnMS <- function(motif_enrichment, motif = NULL, output_type = "matrix", sensitivity_method = "1_minus_norm_score") {
  # 1. Validate inputs
  if (!is.data.frame(motif_enrichment) || !all(c("MOTIF", "Score") %in% colnames(motif_enrichment))) {
    stop("'motif_enrichment' must be a data frame with 'MOTIF' and 'Score' columns.")
  }
  if (nrow(motif_enrichment) == 0) stop("Input 'motif_enrichment' data frame is empty.")
  if (!output_type %in% c("matrix", "number")) stop("'output_type' must be 'matrix' or 'number'.")

  # 2. Find and validate the target motif
  if (is.null(motif) || motif == "top") {
    target_motif <- findTopMotif(motif_enrichment)
    message("No motif provided, using top-scoring motif: ", target_motif)
  } else {
    target_motif <- motif
  }
  kmer_size <- nchar(motif_enrichment$MOTIF[1])
  valInputMotif(motif = target_motif, kmer_size = kmer_size, available_motifs = motif_enrichment$MOTIF)

  # 3. Generate all SNVs for the reference motif
  # For now, assume DNA. Could be made flexible later by checking motif content.
  nucleotides_used <- if (grepl("U", target_motif)) c('A', 'C', 'G', 'U') else c('A', 'C', 'G', 'T')
  all_variants <- genMotifVar(motif = target_motif, type = if('U' %in% nucleotides_used) "RNA" else "DNA")

  # 4. Look up scores for the reference and all SNVs
  reference_score <- motif_enrichment$Score[motif_enrichment$MOTIF == target_motif]

  # Match all variants at once for efficiency
  variant_indices <- match(all_variants, motif_enrichment$MOTIF)
  valid_indices <- !is.na(variant_indices)

  variant_scores <- motif_enrichment$Score[variant_indices[valid_indices]]
  names(variant_scores) <- all_variants[valid_indices]

  if (length(variant_scores) == 0) {
    warning("No variants of the target motif were found in the enrichment data. Cannot calculate MS.")
    return(if (output_type == "matrix") NA else NA_real_)
  }

  # 5. Calculate sensitivity for each SNV
  snv_sensitivities <- calSNVSens(reference_score, variant_scores, method = sensitivity_method)

  # 6. Return result based on output_type
  if (output_type == "matrix") {
    # Format the matrix
    ms_matrix <- formatMS(snv_sensitivities, target_motif, kmer_size, nucleotides = nucleotides_used)
    # NEW LINE: Attach the motif name as an attribute
    attr(ms_matrix, "motif_name") <- target_motif
    return(ms_matrix)
  } else { # "number"
    return(mean(snv_sensitivities, na.rm = TRUE))
  }
}


#' Calculate Sensitivity Score for SNVs
#'
#' @description Internal helper to calculate sensitivity based on reference
#' and variant scores using a specified method.
#'
#' @param reference_score Numeric score of the reference motif.
#' @param variant_scores Named numeric vector or list, where names are SNV motifs
#'   and values are their scores.
#' @param method Character string defining calculation method (e.g.,
#'   "1_minus_norm_score", "score_diff", "log2_ratio").
#'
#' @return A named numeric vector/list similar to `variant_scores` but containing
#'   the calculated sensitivity scores.
#' @keywords internal
calSNVSens <- function(reference_score, variant_scores, method = "1_minus_norm_score") {
  # Ensure variant_scores is a simple numeric vector for calculations
  scores_vec <- unlist(variant_scores)

  sensitivity_scores <- switch(
    toupper(method),
    "1_MINUS_NORM_SCORE" = {
      # Assumes scores were normalized from 0-1, where 1 is the best.
      # Sensitivity is high if the variant score is low.
      1 - scores_vec
    },
    "SCORE_DIFF" = {
      # Sensitivity is the difference between the reference and the variant.
      reference_score - scores_vec
    },
    "LOG2_RATIO" = {
      # Use a pseudocount to avoid division by zero or log of zero.
      pseudocount <- 1e-6
      log2((reference_score + pseudocount) / (scores_vec + pseudocount))
    },
    stop("Unsupported sensitivity 'method': ", method)
  )

  # Restore original names to the calculated scores
  names(sensitivity_scores) <- names(variant_scores)
  return(sensitivity_scores)
}


#' Format Mutational Sensitivity Scores into a Matrix
#'
#' @description Internal helper to arrange SNV sensitivity scores into the
#' standard Nucleotide x Position matrix format.
#'
#' @param snv_sensitivity_scores Named numeric vector/list of sensitivity scores
#'   (output from `calSNVSens`).
#' @param reference_motif Character string, the original motif.
#' @param K Integer, the length of the motif.
#' @param nucleotides Character vector of possible nucleotides.
#'
#' @return A matrix (rows=nucleotides, cols=positions 1 to K) with sensitivity scores.
#'   The cell corresponding to the original nucleotide at a position might be NA or 0.
#' @keywords internal
formatMS <- function(snv_sensitivity_scores, reference_motif, K, nucleotides = c('A', 'C', 'G', 'T')) {
  # Initialize matrix (Nucs x Pos) with NAs
  ms_matrix <- matrix(NA_real_,
                      nrow = length(nucleotides),
                      ncol = K,
                      dimnames = list(nucleotides, paste0("Pos", 1:K)))

  ref_chars <- strsplit(reference_motif, "")[[1]]

  # Loop through each position and nucleotide to fill the matrix
  for (pos in 1:K) {
    for (nuc in nucleotides) {
      # If this is the original nucleotide at this position, score is NA (no change)
      if (nuc == ref_chars[pos]) {
        ms_matrix[nuc, pos] <- NA_real_
        next
      }

      # Construct the variant motif string to look up its score
      variant_chars <- ref_chars
      variant_chars[pos] <- nuc
      variant_motif <- paste0(variant_chars, collapse = "")

      # Find the sensitivity score for this variant
      if (variant_motif %in% names(snv_sensitivity_scores)) {
        ms_matrix[nuc, pos] <- snv_sensitivity_scores[[variant_motif]]
      }
    }
  }
  return(ms_matrix)
}

#-------------------------------------------------------------------------------
