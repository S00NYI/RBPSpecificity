#' RBPSpecificity: Analyze RBP Binding Specificity and Mutational Sensitivity
#'
#' @description
#' The `RBPSpecificity` package provides a suite of tools for researchers
#' working with RNA Binding Proteins (RBPs) to analyze binding specificity from
#' high-throughput sequencing data such as eCLIP (enhanced CrossLinking and
#' ImmunoPrecipitation) and RBNS (RNA Bind-n-Seq).
#'
#' It allows users to quantify K-mer enrichment within RBP binding sites,
#' assess the Inherent Specificity (IS) of binding motifs, and evaluate the
#' Mutational Sensitivity (MS) of these motifs to single nucleotide variations.
#' The package aims to streamline these common bioinformatic analyses and
#' provide clear visualizations for interpreting RBP binding characteristics.
#'
#' **Note:** This package is currently under active development.
#'
#' @section Key Functions:
#' The main workflow of the package revolves around the following exported functions:
#' \itemize{
#'   \item \code{\link{motifEnrichment}}: Calculates K-mer enrichment from peak data (e.g., BED files or GRanges representing eCLIP peaks), with options for background correction.
#'   \item \code{\link{returnIS}}: Determines the Inherent Specificity (IS) score for a given motif based on its enrichment relative to other motifs.
#'   \item \code{\link{plotIS}}: Visualizes the distribution of motif enrichment scores and annotates the IS for a selected motif.
#'   \item \code{\link{returnMS}}: Quantifies the Mutational Sensitivity (MS) of a motif to all possible single nucleotide variations, providing results as a matrix or a single average score.
#'   \item \code{\link{plotMS}}: Generates plots to visualize Mutational Sensitivity profiles, highlighting the impact of nucleotide changes at each position within a motif.
#' }
#'
#' @docType package
#' @name RBPSpecificity-package
#' @aliases RBPSpecificity
#' @author Soon Yi \email{cu.soonyi@gmail.com}
#'
#' @seealso
#' Useful links:
#' \itemize{
#'   \item \url{https://github.com/S00NYI/RBPSpecificity} (Development repository)
#'   \item Report bugs at \url{https://github.com/S00NYI/RBPSpecificity/issues}
#' }
NULL # This NULL is important; it signifies the end of the roxygen block's association.

#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
