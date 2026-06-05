<p align="center">
  <img src="man/figures/logo.png" alt="RBPSpecificity" height="200" />
</p>

<!-- badges: start -->
[![R-CMD-check](https://img.shields.io/badge/R%20CMD%20check-passing-brightgreen)](https://github.com/S00NYI/RBPSpecificity)
[![License: GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R Version](https://img.shields.io/badge/R-%3E%3D%204.3-blue)](https://cran.r-project.org/)
<!-- badges: end -->

**RBPSpecificity** provides tools to calculate Inherent Specificity (IS) and Variation Sensitivity (VS) metrics for RNA-binding proteins from enrichment data.

## Installation

### From GitHub (Development)

```r
# Install pak if needed
if (!requireNamespace("pak", quietly = TRUE))
    install.packages("pak")

pak::pkg_install("S00NYI/RBPSpecificity")
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

# Prepare enrichment data (requires MOTIF and Score columns in [0, 1] range)
enrichment <- data.frame(
  MOTIF = c("UUUUU", "UUUUC", "UUUUG"),
  Score = c(0.95, 0.45, 0.40)
)

# Calculate Inherent Specificity
is_value <- returnIS(enrichment)

# Calculate Variation Sensitivity
vs_value <- returnVS(enrichment, output_type = "number")
vs_matrix <- returnVS(enrichment, output_type = "matrix")

# Visualize
plotIS(enrichment)
plotVS(enrichment)
```

## Features

- **Inherent Specificity (IS)** - Quantify binding specificity from enrichment distributions
- **Variation Sensitivity (VS)** - Measure sensitivity to single nucleotide changes
- **De novo enrichment** - Calculate k-mer enrichment directly from CLIP peaks
- **Visualization** - Distribution plots and sensitivity profiles

## Available Functions

### `motifEnrichment`
Main workflow function to calculate K-mer enrichment.
- **coordinates**: data frame or GRanges containing genomic coordinates.
- **species_or_build**: text string for genome build (e.g. "hg38").
- **K**: integer K-mer length (e.g. 5).
- **extension**: numeric vector `c(five_prime, three_prime)` indicating shifting of 5' and 3' ends. Positive extends, negative trims (default: `c(0,0)`).
- **method**: enrichment model: "anr" (default, any number of repetitions, conditional Binomial test) or "zoops" (zero-or-one occurrence per sequence, Hypergeometric test).
- **bkg_iter**: number of background iterations for local shift (default: 100).
- **bkg_min_dist**: min shift distance for background (default: 500).
- **bkg_max_dist**: max shift distance for background (default: 1000).
- **scramble_bkg**: logical, scramble background sequences to control for nucleotide composition (default: FALSE).
- **nucleic_acid_type**: "DNA" or "RNA" (default: "DNA").

### `returnIS`
Calculate Inherent Specificity (IS).
- **motif_enrichment**: result from `motifEnrichment`.
- **motif**: target motif (default: top scoring).
- **return_type**: "specific" (default) or "all".

### `returnVS`
Calculate Variation Sensitivity (VS).
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

### `plotVS`
Visualize Variation Sensitivity matrix.
- **motif_enrichment**: result from `motifEnrichment`.
- **motif**: reference motif (default: top scoring).

> **Note**: `returnIS()`, `returnVS()`, `plotIS()`, and `plotVS()` accept any data frame with `MOTIF` and `Score` columns. Pre-computed enrichment data from RBNS, RNACompete, or other sources can be used directly without running `motifEnrichment()`, as long as the input is formatted with the correct column names and `Score` values are normalized to a [0, 1] scale.

## Enrichment Calculation

`motifEnrichment()` calculates k-mer enrichment from sequencing results such as peak files from CLIP experiments.

### Motif Occurrence Models

From peak sequences, k-mer enrichment is calculated using one of two motif occurrence models via the `method` parameter:

- **ANR (Any Number of Repetitions)** — Default. Counts every occurrence of a motif within each sequence, retaining multiplicity information. A sequence with 5× UGUGU contributes 5 to the UGUGU count. Enrichment significance is evaluated using a **conditional binomial test** on aggregate occurrence counts, under the assumption that occurrences are Poisson-distributed.

- **ZOOPS (Zero-or-One Occurrence Per Sequence)** — Treats each sequence as a binary outcome: motif present (≥1 occurrence) or absent. A sequence with 5× UGUGU contributes the same as one with 1× UGUGU. Enrichment significance is evaluated using the **cumulative hypergeometric distribution**.

ANR is sensitive to motif multiplicity and cooperative binding hotspots, whereas ZOOPS is robust to repetitive elements and compositional noise. Due to these differing statistical assumptions, ANR and ZOOPS rankings and significance scores are expected to diverge, particularly for motifs with high multiplicity or repetitive backgrounds.

### Background Correction

To correct for motif, sequence, and nucleotide bias in the genomic background, `motifEnrichment()` also calculates background enrichment:

1. Background regions are generated by locally shifting peak coordinates within a configurable distance range (default: ±500–1000 bp, set via `bkg_min_dist` and `bkg_max_dist`). This preserves local genomic and transcript context — if a peak is in a 3'UTR, the background is also in/near that 3'UTR, controlling for regional composition biases.

2. The shifting process is repeated `bkg_iter` times (default: 100) and background k-mer counts are averaged across iterations for robust estimation.

3. Background k-mer enrichment is calculated using the same model (ANR or ZOOPS) as the peak sequences.

4. Background enrichment is subtracted from peak enrichment to produce the corrected enrichment score.

5. Optionally, background sequences can be scrambled (`scramble_bkg = TRUE`, default: FALSE) to control for nucleotide composition while removing sequence-specific motif signals.

### Output

The enrichment difference (target − background) is computed for each k-mer, then min-max normalized across all k-mers to [0, 1] to produce the `Score` column. This `Score` is the primary metric consumed by `returnIS()` and `returnVS()`. The output also includes `log2FC` and `padj` (BH-adjusted p-values) as supplementary columns.

## Documentation

See the [package vignette](vignettes/RBPSpecificity.Rmd) for detailed examples covering:

1. Pre-computed enrichment data (e.g. RBNS or RNACompete)
2. De novo enrichment from CLIP peaks
3. Multi-RBP comparison
4. Visualization options

For detailed usage, see **Figure_Scripts** folder in:
https://github.com/S00NYI/BITS_Specificity

## Related Packages

**RBPSpecificity** works alongside [RBPBind](https://github.com/S00NYI/RBPBind), which simulates competitive RBP binding using equilibrium kinetics. Together, these packages enable:

- Calculating specificity metrics from experimental data (RBPSpecificity)
- Simulating competitive binding scenarios (RBPBind)

## Citation

If you use RBPSpecificity in your research, please cite:

> Yi S, Singh SS, Ye X, Krishna R, Jankowsky E, Luna JM. (2025). 
> *Inherent Specificity and Mutational Sensitivity as Quantitative Metrics for RBP Binding.* 
> bioRxiv. https://www.biorxiv.org/content/10.1101/2025.03.28.646018v2

This package was developed in part with **Google Antigravity**.

For the original R implementation of RBP specificity analysis, see **Deprecated** folder in:
https://github.com/S00NYI/BITS_Specificity

## Further Reading

1. Kjosmoen, T. & Ryen, T. (2009). [Exploring the Combinatorics of Motif Alignments](https://www.semanticscholar.org/paper/Exploring-the-Combinatorics-of-Motif-Alignments-Kjosmoen-Ryen/8a9546c437dad20983fefdeb7222d7ec65a04176). *Semantic Scholar*.

2. Jayaram, N., Bhowmick, P., & Martin, A.C.R. (2021). [Computational Methods for Identification of Regulatory Elements in Genomic Sequences](https://link.springer.com/protocol/10.1007/978-1-0716-1307-8_3). *Methods in Molecular Biology*.

3. Oregon State University. [Chapter 2: Sequence Motifs — Applied Bioinformatics](https://open.oregonstate.education/appliedbioinformatics/chapter/chapter-2-sequence-motifs/).

4. University of Wisconsin–Madison, BMI 776. [Motif Modeling (Lecture Notes)](https://www.biostat.wisc.edu/bmi776/lectures/motif-modeling.pdf).
