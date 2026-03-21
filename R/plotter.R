# Plotting functions for visualization

#' Plot Inherent Specificity (IS) Distribution
#'
#' @description Generates a histogram of motif enrichment scores, annotated with
#' the median score, the score of a specific motif, and its IS value.
#'
#' @param motif_enrichment A data frame containing motif enrichment scores
#'   (output from `motifEnrichment()`). Requires 'MOTIF' and 'Score' columns.
#' @param motif Character string, the specific motif to highlight and calculate
#'   IS for. If NULL or "top", the top-scoring motif is used (default: NULL).
#' @param bins Integer, number of bins for the histogram (default: 50).
#' @param ... Additional arguments passed to ggplot theme layers or geoms.
#'
#' @return A ggplot object representing the annotated histogram.
#' @export
#' @import ggplot2
#'
#' @examples
#' # Dummy data
#' df <- data.frame(MOTIF = c("AAAA", "CCCC", "GGGG", "UUUU"), Score = c(10, 2, 5, 1))
#' 
#' # Plot specific motif
#' plotIS(df, motif = "AAAA")
#' 
#' # Plot top motif (AAAA)
#' plotIS(df)
#' 
#' # Change number of bins
#' plotIS(df, bins = 20)
plotIS <- function(motif_enrichment, motif = NULL, bins = 50, ...) {
  # 1. Validate input and find target motif
  if (!is.data.frame(motif_enrichment) || !all(c("MOTIF", "Score") %in% colnames(motif_enrichment))) {
    stop("'motif_enrichment' must be a data frame with 'MOTIF' and 'Score' columns.")
  }
  if (is.null(motif) || motif == "top") {
    target_motif <- findTopMotif(motif_enrichment)
  } else {
    target_motif <- motif
  }
  valInputMotif(target_motif, nchar(motif_enrichment$MOTIF[1]), motif_enrichment$MOTIF)

  # 2. Get scores and calculate IS value
  target_score <- motif_enrichment$Score[motif_enrichment$MOTIF == target_motif]
  median_score <- stats::median(motif_enrichment$Score, na.rm = TRUE)
  is_value <- returnIS(motif_enrichment, motif = target_motif)

  # 3. Create the histogram using ggplot2
  plot <- ggplot2::ggplot(motif_enrichment, ggplot2::aes(x = Score)) +
    ggplot2::geom_histogram(bins = bins, fill = "grey50", alpha = 0.7) +
    ggplot2::geom_vline(xintercept = median_score, color = "skyblue", linetype = "dashed", linewidth = 1) +
    ggplot2::geom_vline(xintercept = target_score, color = "salmon", linetype = "solid", linewidth = 1) +
    ggplot2::annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
                      label = paste("IS for", target_motif, "=", round(is_value, 2)), size = 5) +
    ggplot2::labs(
      title = "Distribution of K-mer Enrichment Scores",
      subtitle = paste("Target Motif:", target_motif),
      x = "Normalized Enrichment Score",
      y = "Frequency"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text = ggplot2::element_text(size=14),
                   axis.title = ggplot2::element_text(size=14, face = 'bold'),
                   legend.text = ggplot2::element_text(size=14))

  return(plot)
}

#' Plot Variation Sensitivity (VS) Matrix
#'
#' @description Visualizes the variation sensitivity scores for a motif,
#' typically showing sensitivity to SNVs at each position.
#'
#' @param motif_enrichment A data frame containing motif enrichment scores
#'   (output from `motifEnrichment()`). Requires 'MOTIF' and 'Score' columns.
#' @param motif Character string, the reference motif for which VS was calculated.
#'   If NULL or "top", the top-scoring motif is used (default: NULL).
#' @param ... Additional arguments passed to ggplot theme layers or geoms.
#'
#' @return A ggplot object visualizing the VS matrix (e.g., using points sized
#' or colored by sensitivity).
#' @export
#' @import ggplot2
#' @importFrom reshape2 melt
#'
#' @examples
#' # Mock data with variants
#' motifs <- c("AAAA", "CAAA", "GAAA", "TAAA")
#' scores <- c(10,      5,      4,      3)
#' df <- data.frame(MOTIF = motifs, Score = scores)
#' 
#' # Plot VS for "AAAA"
#' plotVS(df, motif = "AAAA")
plotVS <- function(motif_enrichment, motif = NULL, ...) {
  # 1. Get the VS matrix by calling returnVS
  vs_matrix <- returnVS(motif_enrichment, motif = motif, output_type = "matrix")

  if(all(is.na(vs_matrix))) {
    warning("VS matrix contains all NAs, cannot generate plot.")
    return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = "VS data not available"))
  }

  # 2. Melt the matrix for ggplot
  vs_df <- as.data.frame(vs_matrix)
  vs_df$Nucleotide <- rownames(vs_df)
  vs_long <- reshape2::melt(vs_df, id.vars = "Nucleotide", variable.name = "Position", value.name = "Sensitivity")
  vs_long$Position <- as.numeric(gsub("Pos", "", vs_long$Position))

  # 3. Determine the subtitle text safely
  ref_motif_name <- attr(vs_matrix, "motif_name")
  subtitle_text <- if (!is.null(ref_motif_name)) {
    paste("Reference Motif:", ref_motif_name)
  } else {
    # This is a fallback in case the attribute is missing
    "Reference Motif: Top Motif"
  }

  # 4. Create the plot
  plot <- ggplot2::ggplot(vs_long, ggplot2::aes(x = Position, y = Sensitivity, fill = Nucleotide)) +
    ggplot2::geom_point(shape = 21, size = 8, stroke = 1, na.rm = TRUE) +
    ggplot2::scale_fill_manual(values = c("G" = "#F5C714", "A" = "#70BF52", "C" = "#3D94D1", "U" = "#E0546C", "T" = "#E0546C")) +
    ggplot2::scale_x_continuous(breaks = seq_len(ncol(vs_matrix))) +
    ggplot2::scale_y_continuous(limits = c(1.05, -0.05), breaks = seq(0, 1, by = 0.2), trans = 'reverse') +
    ggplot2::labs(
      title = "Variation Sensitivity Profile",
      subtitle = subtitle_text,
      x = "Position in Motif",
      y = "Variation Sensitivity (1 - Variant Score)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text = ggplot2::element_text(size=14),
                   axis.title = ggplot2::element_text(size=14, face = 'bold'),
                   legend.text = ggplot2::element_text(size=14))

  return(plot)
}

#-------------------------------------------------------------------------------

