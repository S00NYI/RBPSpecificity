#' RBPSpecificity: Analyze RBP Binding Specificity and Variation Sensitivity
#'
#' @description
#' The `RBPSpecificity` package provides a suite of tools for researchers
#' working with RNA Binding Proteins (RBPs) to analyze binding specificity from
#' high-throughput sequencing data such as CLIP (CrossLinking and
#' ImmunoPrecipitation) and RBNS (RNA Bind-n-Seq).
#'
#' It allows users to quantify K-mer enrichment within RBP binding sites,
#' assess the inherent specificity (IS) and the variation sensitivity (VS) of
#' RBP towards the target motifs. The package aims to streamline bioinformatic
#' analyses and provide clear visualizations for interpreting RBP specificity.
#'
#' **Note:** This package is under active development.
#'
#' @section Key Functions:
#' The main workflow of the package revolves around the following exported functions:
#' \itemize{
#'   \item \code{\link{motifEnrichment}}: Calculates K-mer enrichment from peak data...
#'   \item \code{\link{returnSpecificity}}: Quantifies the specificity of RBP...
#'   \item \code{\link{plotSpecificity}}: Visualizes the distribution of motif enrichment scores with annotation...
#'   \item \code{\link{returnSensitivity}}: Quantifies the sensitivity of RBP...
#'   \item \code{\link{plotSensitivity}}: Generates plots to visualize sensitivity
#' }
#'
#' @author Soon Yi \email{cu.soonyi@gmail.com}
#' @aliases RBPSpecificity RBPSpecificity-package
#'
#' @seealso
#' Useful links:
#' \itemize{
#'   \item \url{https://github.com/S00NYI/RBPSpecificity}
#'         (Development repository)
#'   \item Report bugs at \url{https://github.com/S00NYI/RBPSpecificity/issues}
#' }
"_PACKAGE" # Document this special string directly

# Suppress R CMD check notes for NSE variables
utils::globalVariables(c(
    "COUNT", "AVG_BKG_COUNT", "MOTIF", "EnrichmentScore", "Score",
    "Sensitivity", "Nucleotide", "."
))

## usethis namespace: start
## usethis namespace: end
NULL
