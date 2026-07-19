# Plotting functions for visualization

#' Plot Specificity Distribution
#'
#' @description Generates a histogram of motif enrichment scores, annotated with
#' the median score, the score of a specific motif, and its specificity value.
#'
#' @param motif_enrichment A data frame containing motif enrichment scores
#'   (output from `motifEnrichment()`). Requires 'MOTIF' and 'Score' columns.
#' @param motif Character string, the specific motif to highlight and calculate
#'   specificity for. If NULL or "top", the top-scoring motif is used (default: NULL).
#' @param bins Integer, number of bins for the histogram (default: 50).
#' @param ... Additional arguments passed to ggplot theme layers or geoms.
#'
#' @return A ggplot object representing the annotated histogram.
#' @export
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline annotate labs theme_bw theme element_text
#'
#' @examples
#' # Dummy data
#' df <- data.frame(MOTIF = c("AAAA", "CCCC", "GGGG", "UUUU"), Score = c(10, 2, 5, 1))
#'
#' # Plot specific motif
#' plotSpecificity(df, motif = "AAAA")
#'
#' # Plot top motif (AAAA)
#' plotSpecificity(df)
#'
#' # Change number of bins
#' plotSpecificity(df, bins = 20)
plotSpecificity <- function(motif_enrichment, motif = NULL, bins = 50, ...) {
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

    # 2. Get scores and calculate specificity value
    target_score <- motif_enrichment$Score[motif_enrichment$MOTIF == target_motif]
    median_score <- stats::median(motif_enrichment$Score, na.rm = TRUE)
    spec_value <- returnSpecificity(motif_enrichment, motif = target_motif)

    # 3. Create the histogram using ggplot2
    plot <- ggplot2::ggplot(motif_enrichment, ggplot2::aes(x = Score)) +
        ggplot2::geom_histogram(bins = bins, fill = "grey50", alpha = 0.7) +
        ggplot2::geom_vline(xintercept = median_score, color = "skyblue", linetype = "dashed", linewidth = 1) +
        ggplot2::geom_vline(xintercept = target_score, color = "salmon", linetype = "solid", linewidth = 1) +
        ggplot2::annotate("text",
            x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
            label = paste("Specificity for", target_motif, "=", round(spec_value, 2)), size = 5
        ) +
        ggplot2::labs(
            title = "Distribution of K-mer Enrichment Scores",
            subtitle = paste("Target Motif:", target_motif),
            x = "Enrichment Score",
            y = "Frequency"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.text = ggplot2::element_text(size = 14),
            axis.title = ggplot2::element_text(size = 14, face = "bold"),
            legend.text = ggplot2::element_text(size = 14)
        )

    return(plot)
}

#' Plot Sensitivity Matrix
#'
#' @description Visualizes the sensitivity scores for a motif,
#' typically showing sensitivity to SNVs at each position.
#'
#' @param motif_enrichment A data frame containing motif enrichment scores
#'   (output from `motifEnrichment()`). Requires 'MOTIF' and 'Score' columns.
#' @param motif Character string, the reference motif for which sensitivity was calculated.
#'   If NULL or "top", the top-scoring motif is used (default: NULL).
#' @param ... Additional arguments passed to ggplot theme layers or geoms.
#'
#' @return A ggplot object visualizing the sensitivity matrix (e.g., using points sized
#' or colored by sensitivity).
#' @export
#' @importFrom ggplot2 ggplot aes geom_point scale_fill_manual scale_x_continuous scale_y_continuous labs theme_bw theme element_text theme_void
#' @importFrom reshape2 melt
#'
#' @examples
#' # Mock data with variants
#' motifs <- c("AAAA", "CAAA", "GAAA", "TAAA")
#' scores <- c(10, 5, 4, 3)
#' df <- data.frame(MOTIF = motifs, Score = scores)
#'
#' # Plot sensitivity for "AAAA"
#' plotSensitivity(df, motif = "AAAA")
plotSensitivity <- function(motif_enrichment, motif = NULL, ...) {
    # 1. Get the sensitivity matrix by calling returnSensitivity
    sens_matrix <- returnSensitivity(motif_enrichment, motif = motif, output_type = "matrix")

    if (all(is.na(sens_matrix))) {
        warning("Sensitivity matrix contains all NAs, cannot generate plot.")
        return(ggplot2::ggplot() +
            ggplot2::theme_void() +
            ggplot2::labs(title = "Sensitivity data not available"))
    }

    # 2. Melt the matrix for ggplot
    sens_df <- as.data.frame(sens_matrix)
    sens_df$Nucleotide <- rownames(sens_df)
    sens_long <- reshape2::melt(sens_df, id.vars = "Nucleotide", variable.name = "Position", value.name = "Sensitivity")
    sens_long$Position <- as.numeric(gsub("Pos", "", sens_long$Position))

    # 3. Determine the subtitle text safely
    ref_motif_name <- attr(sens_matrix, "motif_name")
    subtitle_text <- if (!is.null(ref_motif_name)) {
        paste("Reference Motif:", ref_motif_name)
    } else {
        "Reference Motif: Top Motif"
    }

    # 4. Create the plot
    plot <- ggplot2::ggplot(sens_long, ggplot2::aes(x = Position, y = Sensitivity, fill = Nucleotide)) +
        ggplot2::geom_point(shape = 21, size = 8, stroke = 1, na.rm = TRUE) +
        ggplot2::scale_fill_manual(values = c("G" = "#F5C714", "A" = "#70BF52", "C" = "#3D94D1", "U" = "#E0546C", "T" = "#E0546C")) +
        ggplot2::scale_x_continuous(breaks = seq_len(ncol(sens_matrix))) +
        ggplot2::scale_y_continuous(limits = c(1.05, -0.05), breaks = seq(0, 1, by = 0.2), trans = "reverse") +
        ggplot2::labs(
            title = "Sensitivity Profile",
            subtitle = subtitle_text,
            x = "Position in Motif",
            y = "Sensitivity"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.text = ggplot2::element_text(size = 14),
            axis.title = ggplot2::element_text(size = 14, face = "bold"),
            legend.text = ggplot2::element_text(size = 14)
        )

    return(plot)
}
