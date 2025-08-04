# Functions related to calculating Inherent Specificity (IS)

#' Calculate Inherent Specificity (IS)
#'
#' @description Calculates the Inherent Specificity (IS) for a given motif based
#' on its enrichment score relative to the median enrichment score of all motifs.
#'
#' @param motif_enrichment A data frame containing motif enrichment scores,
#'   typically output from `motifEnrichment()`. Must have columns like 'MOTIF'
#'   and 'Score' (or 'Enrichment').
#' @param motif Character string, the specific motif for which to calculate IS.
#'   If NULL or "top", the motif with the highest score is used (default: NULL).
#'
#' @return A single numeric value representing the IS score for the motif.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming enrichment_results is output from motifEnrichment()
#' is_value_top <- returnIS(enrichment_results)
#' is_value_specific <- returnIS(enrichment_results, motif = "UGUGU")
#' }
returnIS <- function(motif_enrichment, motif = NULL) {
  # Validate motif_enrichment input
  if (!is.data.frame(motif_enrichment) || !all(c("MOTIF", "Score") %in% colnames(motif_enrichment))) {
    stop("'motif_enrichment' must be a data frame with 'MOTIF' and 'Score' columns.")
  }
  if (nrow(motif_enrichment) == 0) {
    stop("Input 'motif_enrichment' data frame is empty.")
  }

  # 1. Find target motif
  if (is.null(motif) || motif == "top") {
    target_motif <- findTopMotif(motif_enrichment)
    message("No motif provided, using top-scoring motif: ", target_motif)
  } else {
    target_motif <- motif
  }

  # 2. Validate the target motif
  kmer_size <- nchar(motif_enrichment$MOTIF[1])
  available_motifs <- motif_enrichment$MOTIF
  valInputMotif(motif = target_motif, kmer_size = kmer_size, available_motifs = available_motifs)

  # 3. Get score for target motif and all scores
  target_score <- motif_enrichment$Score[motif_enrichment$MOTIF == target_motif]
  all_scores <- motif_enrichment$Score

  # 4. Calculate IS by calling the helper
  is_value <- calIS(scores = all_scores, target_score = target_score)

  # 5. Return IS value
  return(is_value)
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
