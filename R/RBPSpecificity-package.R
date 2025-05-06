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
#'   \item \code{\link{motifEnrichment}}: Calculates K-mer enrichment from peak data...
#'   \item \code{\link{returnIS}}: Quantifies the Inherent Specificity (IS) of RBP...
#'   \item \code{\link{plotIS}}: Visualizes the distribution of motif enrichment scores with annotation...
#'   \item \code{\link{returnMS}}: Quantifies the Mutational Sensitivity (MS) of RBP...
#'   \item \code{\link{plotMS}}: Generates plots to visualize Mutational Sensitivity
#' }
#'
#' @author Soon Yi \email{cu.soonyi@gmail.com}
#' @aliases RBPSpecificity RBPSpecificity-package # Ensures ?RBPSpecificity and ?RBPSpecificity-package work
#'
#' @seealso
#' Useful links:
#' \itemize{
#'   \item \url{https://github.com/S00NYI/RBPSpecificity} (Development repository)
#'   \item Report bugs at \url{https://github.com/S00NYI/RBPSpecificity/issues}
#' }
"_PACKAGE" # Document this special string directly

## usethis namespace: start
## usethis namespace: end
NULL
