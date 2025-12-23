# Functions related to calculating Inherent Specificity (IS)

#' Calculate Inherent Specificity (IS)
#'
#' @description Calculates the Inherent Specificity (IS) for a given motif based
#' on its enrichment score relative to the median enrichment score of all motifs.
#'
#' @param motif_enrichment A data frame containing motif enrichment scores,
#'   typically output from `motifEnrichment()`. Must have columns 'MOTIF'
#'   and 'Score'.
#' @param motif Character string, the specific motif for which to calculate IS.
#'   Only used if `return_type = "specific"`. If NULL or "top", the top-scoring
#'   motif is used. Default: NULL.
#' @param return_type Character string, either "specific" (default) or "all".
#'   If "specific", returns a single IS value for the target motif. If "all",
#'   returns a data frame with IS values for every motif.
#'
#' @return A single numeric IS value if `return_type = "specific"`, or a data
#'   frame with 'MOTIF' and 'IS' columns if `return_type = "all"`.
#' @export
#'
#' @examples
#' # Create a dummy enrichment dataframe
#' df <- data.frame(MOTIF = c("AAAAA", "CCCCC", "GGGGG"), Score = c(10, 2, 5))
#' 
#' # Calculate IS for a specific motif
#' returnIS(df, motif = "AAAAA")
#' 
#' # Calculate IS for all motifs
#' returnIS(df, return_type = "all")
returnIS <- function(motif_enrichment, motif = NULL, return_type = "specific") {
  # --- Input Validation ---
  if (!is.data.frame(motif_enrichment) || !all(c("MOTIF", "Score") %in% colnames(motif_enrichment))) {
    stop("'motif_enrichment' must be a data frame with 'MOTIF' and 'Score' columns.")
  }
  if (nrow(motif_enrichment) == 0) stop("Input 'motif_enrichment' data frame is empty.")
  if (!return_type %in% c("specific", "all")) stop("'return_type' must be 'specific' or 'all'.")

  # --- Logic based on return_type ---
  if (return_type == "specific") {
    # --- Calculate for a single motif (original logic) ---
    if (is.null(motif) || motif == "top") {
      target_motif <- findTopMotif(motif_enrichment)
      message("No motif provided, using top-scoring motif: ", target_motif)
    } else {
      target_motif <- motif
    }
    kmer_size <- nchar(motif_enrichment$MOTIF[1])
    valInputMotif(motif = target_motif, kmer_size = kmer_size, available_motifs = motif_enrichment$MOTIF)
    target_score <- motif_enrichment$Score[motif_enrichment$MOTIF == target_motif]
    all_scores <- motif_enrichment$Score
    is_value <- calIS(scores = all_scores, target_score = target_score)
    return(is_value)

  } else { # return_type == "all"
    # --- Calculate for all motifs (vectorized) ---
    message("Calculating IS for all motifs...")
    all_scores <- motif_enrichment$Score
    median_score <- stats::median(all_scores, na.rm = TRUE)

    if (is.na(median_score) || median_score == 0) {
      is_values <- rep(NA_real_, nrow(motif_enrichment))
    } else {
      is_values <- all_scores / median_score
    }

    results_df <- data.frame(
      MOTIF = motif_enrichment$MOTIF,
      IS = is_values
    )
    return(results_df)
  }
}


#' Calculate Inherent Specificity Core Logic
#'
#' @description Internal helper to calculate the ratio of a target score
#' to the median score.
#'
#' @param scores Numeric vector of all scores.
#' @param target_score Numeric value, the score of the specific item of interest.
#'
#' @return Single numeric value (target_score / median(scores)).
#' @importFrom stats median
#' @keywords internal
calIS <- function(scores, target_score) {
  # 1. Calculate median of scores (handle NAs)
  median_score <- stats::median(scores, na.rm = TRUE)

  # 2. Handle division by zero or NA median if necessary
  if (is.na(median_score) || median_score == 0) {
    return(NA_real_)
  }

  # 3. Calculate and return ratio
  return(target_score / median_score)
}

#-------------------------------------------------------------------------------
