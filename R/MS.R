#' Calculate Mutational Sensitivity (MS)
#'
#' @description Calculates the Mutational Sensitivity (MS) for a target motif or
#'   for all motifs. MS reflects how much the enrichment score changes upon
#'   single nucleotide variations (SNVs).
#'
#' @param motif_enrichment A data frame with 'MOTIF' and 'Score' columns.
#' @param motif Character string, the reference motif. Only used if
#'   `return_type = "specific"`. If NULL or "top", the top-scoring motif is used.
#' @param return_type Character string, either "specific" (default) or "all".
#'   If "specific", calculates MS for one motif. If "all", calculates an
#'   average MS score for every motif.
#' @param output_type Character string, specifying the output format if
#'   `return_type = "specific"`: "matrix" (default) or "number". Ignored if
#'   `return_type = "all"`.
#' @param sensitivity_method Character string defining how sensitivity is calculated.
#'   Default: "1_minus_norm_score".
#'
#' @return Depends on inputs:
#'   - If `return_type = "specific"` and `output_type = "matrix"`: An MS matrix.
#'   - If `return_type = "specific"` and `output_type = "number"`: A single numeric MS score.
#'   - If `return_type = "all"`: A data frame with 'MOTIF' and 'MS' columns.
#' @export
returnMS <- function(motif_enrichment, motif = NULL, return_type = "specific",
                     output_type = "matrix", sensitivity_method = "1_minus_norm_score") {
  # --- Input Validation ---
  if (!is.data.frame(motif_enrichment) || !all(c("MOTIF", "Score") %in% colnames(motif_enrichment))) {
    stop("'motif_enrichment' must be a data frame with 'MOTIF' and 'Score' columns.")
  }
  if (nrow(motif_enrichment) == 0) stop("Input 'motif_enrichment' data frame is empty.")
  if (!return_type %in% c("specific", "all")) stop("'return_type' must be 'specific' or 'all'.")

  # --- Logic based on return_type ---
  if (return_type == "specific") {
    # ... (The code to find and validate the target motif remains the same) ...
    if (!output_type %in% c("matrix", "number")) stop("'output_type' must be 'matrix' or 'number' for specific return.")
    if (is.null(motif) || motif == "top") {
      target_motif <- findTopMotif(motif_enrichment)
      message("No motif provided, using top-scoring motif: ", target_motif)
    } else {
      target_motif <- motif
    }
    kmer_size <- nchar(motif_enrichment$MOTIF[1])
    valInputMotif(motif = target_motif, kmer_size = kmer_size, available_motifs = motif_enrichment$MOTIF)

    nucleotides_used <- if (grepl("U", target_motif)) c('A', 'C', 'G', 'U') else c('A', 'C', 'G', 'T')
    all_variants <- genMotifVar(motif = target_motif, type = if('U' %in% nucleotides_used) "RNA" else "DNA")

    reference_score <- motif_enrichment$Score[motif_enrichment$MOTIF == target_motif]
    variant_indices <- match(all_variants, motif_enrichment$MOTIF)
    valid_indices <- !is.na(variant_indices)
    variant_scores <- motif_enrichment$Score[variant_indices[valid_indices]]
    names(variant_scores) <- all_variants[valid_indices]

    if (length(variant_scores) == 0) {
      warning("No variants of the target motif found. Cannot calculate MS.")
      return(if (output_type == "matrix") NA else NA_real_)
    }

    snv_sensitivities <- calSNVSens(reference_score, variant_scores, method = sensitivity_method)

    if (output_type == "matrix") {
      ms_matrix <- formatMS(snv_sensitivities, target_motif, reference_score, kmer_size,
                            nucleotides = nucleotides_used, sensitivity_method = sensitivity_method)
      attr(ms_matrix, "motif_name") <- target_motif
      return(ms_matrix)
    } else { # "number"
      score_diff_sensitivities <- calSNVSens(reference_score, variant_scores, method = "score_diff")
      non_zero_sensitivities <- score_diff_sensitivities[score_diff_sensitivities != 0]

      return(mean(non_zero_sensitivities, na.rm = TRUE))
    }

  } else { # return_type == "all"
    # --- Calculate for all motifs ---
    message("Calculating MS for all ", nrow(motif_enrichment), " motifs. This may take a while...")
    all_motifs <- motif_enrichment$MOTIF
    ms_values <- numeric(length(all_motifs))
    pb <- utils::txtProgressBar(min = 0, max = length(all_motifs), style = 3)

    for (i in seq_along(all_motifs)) {
      target_motif <- all_motifs[i]
      nucleotides_used <- if (grepl("U", target_motif)) c('A', 'C', 'G', 'U') else c('A', 'C', 'G', 'T')
      all_variants <- genMotifVar(motif = target_motif, type = if('U' %in% nucleotides_used) "RNA" else "DNA")

      reference_score <- motif_enrichment$Score[i]
      variant_indices <- match(all_variants, all_motifs)
      valid_indices <- !is.na(variant_indices)
      variant_scores <- motif_enrichment$Score[variant_indices[valid_indices]]

      if (length(variant_scores) > 0) {
        score_diffs <- reference_score - variant_scores
        non_zero_diffs <- score_diffs[score_diffs != 0]
        ms_values[i] <- mean(non_zero_diffs, na.rm = TRUE)
      } else {
        ms_values[i] <- NA_real_
      }
      utils::setTxtProgressBar(pb, i)
    }
    close(pb)

    results_df <- data.frame(MOTIF = all_motifs, MS = ms_values)
    return(results_df)
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
formatMS <- function(snv_sensitivity_scores, reference_motif, reference_score, K,
                     nucleotides = c('A', 'C', 'G', 'T'), sensitivity_method = "1_minus_norm_score") {
  # Initialize matrix (Nucs x Pos) with NAs
  ms_matrix <- matrix(NA_real_,
                      nrow = length(nucleotides),
                      ncol = K,
                      dimnames = list(nucleotides, paste0("Pos", 1:K)))

  ref_chars <- strsplit(reference_motif, "")[[1]]

  # Loop through each position and nucleotide to fill the matrix
  for (pos in 1:K) {
    for (nuc in nucleotides) {
      # If this is the original nucleotide at this position
      if (nuc == ref_chars[pos]) {
        # The "mutation" is to itself. Calculate sensitivity based on the chosen method.
        # This makes the row for the reference nucleotide consistent.
        self_sens <- calSNVSens(reference_score, reference_score, method = sensitivity_method)
        ms_matrix[nuc, pos] <- self_sens
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
