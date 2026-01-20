<p align="center">
  <img src="man/figures/logo.png" alt="RBPSpecificity" height="200" />
</p>

<!-- badges: start -->
[![R-CMD-check](https://img.shields.io/badge/R%20CMD%20check-passing-brightgreen)](https://github.com/S00NYI/RBPSpecificity)
[![License: GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R Version](https://img.shields.io/badge/R-%3E%3D%204.3-blue)](https://cran.r-project.org/)
<!-- badges: end -->

**RBPSpecificity** provides tools to calculate Inherent Specificity (IS) and Mutational Sensitivity (MS) metrics for RNA-binding proteins from enrichment data.

## Installation

### From GitHub (Development)

```r
# Install devtools if needed
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("S00NYI/RBPSpecificity")
```

<!-- ### From Bioconductor (After Acceptance)

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("RBPSpecificity")
``` -->

## Quick Start

```r
library(RBPSpecificity)

# Prepare enrichment data (requires MOTIF and Score columns)
enrichment <- data.frame(
  MOTIF = c("UUUUU", "UUUUC", "UUUUG", ...),
  Score = c(10.5, 5.2, 4.8, ...)
)

# Calculate Inherent Specificity
is_value <- returnIS(enrichment)

# Calculate Mutational Sensitivity
ms_value <- returnMS(enrichment, output_type = "number")
ms_matrix <- returnMS(enrichment, output_type = "matrix")

# Visualize
plotIS(enrichment)
plotMS(enrichment)
```

## Features

- **Inherent Specificity (IS)** - Quantify binding specificity from enrichment distributions
- **Mutational Sensitivity (MS)** - Measure sensitivity to single nucleotide changes
- **De novo enrichment** - Calculate k-mer enrichment directly from CLIP peaks
- **Visualization** - Distribution plots and sensitivity profiles

## Available Functions

### `motifEnrichment`
Main workflow function to calculate K-mer enrichment.
- **coordinates**: data frame or GRanges containing genomic coordinates.
- **species_or_build**: text string for genome build (e.g. "hg38").
- **K**: integer K-mer length (e.g. 5).
- **extension**: numeric vector `c(5', 3')` to shift 5'/3' ends. Positive extends, negative trims. (default: `c(0,0)`).
- **enrichment_method**: "subtract" (default), "fold_change", or "log2_fold_change".
- **normalization_method**: "min_max" (default), "z_score", "log2", or "none".
- **log_transform**: logical, apply log transformation (default: TRUE).
- **bkg_iter**: number of background iterations (default: 100).
- **bkg_min_dist**: min shift for background (default: 500).
- **bkg_max_dist**: max shift for background (default: 1000).
- **nucleic_acid_type**: "DNA" or "RNA" (default: "DNA").

### `returnIS`
Calculate Inherent Specificity (IS).
- **motif_enrichment**: result from `motifEnrichment`.
- **motif**: target motif (default: top scoring).
- **return_type**: "specific" (default) or "all".

### `returnMS`
Calculate Mutational Sensitivity (MS).
- **motif_enrichment**: result from `motifEnrichment`.
- **motif**: reference motif (default: top scoring).
- **return_type**: "specific" (default) or "all".
- **output_type**: "matrix" (default) or "number".
- **sensitivity_method**: calculation method (default: "1_minus_norm_score").

### `plotIS`
Visualize IS value on score distribution.
- **motif_enrichment**: result from `motifEnrichment`.
- **motif**: target motif (default: top scoring).
- **bins**: histogram bins (default: 50).

### `plotMS`
Visualize Mutational Sensitivity matrix.
- **motif_enrichment**: result from `motifEnrichment`.
- **motif**: reference motif (default: top scoring).


## Documentation

See the [package vignette](vignettes/RBPSpecificity.Rmd) for detailed examples covering:

1. Pre-computed enrichment data (RBNS)
2. De novo enrichment from peaks (eCLIP)
3. Multi-RBP comparison
4. Visualization options

## Related Packages

**RBPSpecificity** works alongside [RBPBind](https://github.com/S00NYI/RBPBind), which simulates competitive RBP binding using equilibrium kinetics. Together, these packages enable:

- Calculating specificity metrics from experimental data (RBPSpecificity)
- Simulating competitive binding scenarios (RBPBind)

## Citation

If you use RBPSpecificity in your research, please cite:

> Yi S, Singh SS, Ye X, Krishna R, Jankowsky E, Luna JM. (2025). 
> *Inherent Specificity and Mutational Sensitivity as Quantitative Metrics for RBP Binding.* 
> bioRxiv. https://www.biorxiv.org/content/10.1101/2025.03.28.646018v2

## Acknowledgments

This package was developed in part with **Google Antigravity**.

For the original R implementation of RBP specificity analysis, see:
https://github.com/S00NYI/BITS_Specificity
