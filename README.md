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

### From Bioconductor (After Acceptance)

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("RBPSpecificity")
```

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
