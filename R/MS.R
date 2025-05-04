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
  # 1. Find target motif: findTopMotif() if motif is NULL/"top"
  # 2. Validate target motif: valInputMotif()
  # 3. Get score for the target reference motif
  # 4. Generate all SNVs for the reference motif: genMotifVar()
  # 5. Look up scores for all SNVs from motif_enrichment (handle missing variants)
  # 6. Calculate sensitivity for each SNV: calSNVSens()
  # 7. If output_type == "matrix":
  #    a. Format results into matrix: formatMS()
  #    b. Return matrix
  # 8. If output_type == "number":
  #    a. Calculate mean of SNV sensitivity scores (handle NAs)
  #    b. Return average number
  {}
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
  # 1. Implement logic based on `method`
  #    - e.g., for "1_minus_norm_score": `1 - variant_scores` (assumes scores were normalized 0-1)
  #    - e.g., for "score_diff": `reference_score - variant_scores`
  #    - e.g., for "log2_ratio": `log2(reference_score / variant_scores)` (handle zeros/negatives)
  # 2. Return calculated sensitivity scores
  {}
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
  # (Based on original return_MS logic)
  # 1. Initialize matrix (Nucs x Pos) with NAs
  # 2. Loop through positions 1 to K
  # 3. Loop through nucleotides
  # 4. Construct the corresponding variant motif string
  # 5. Find the sensitivity score for this variant in snv_sensitivity_scores
  # 6. Place score in the matrix[nucleotide, position]
  # 7. Handle the cell for the original nucleotide (e.g., set to NA or 0)
  # 8. Return the formatted matrix
  {}
}

#-------------------------------------------------------------------------------
