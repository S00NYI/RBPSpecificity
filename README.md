# RBPSpecificity

**Note:** This package is currently under active development. Function names, arguments, and return values may change.

## Overview

`RBPSpecificity` is an R package designed for analyzing RNA Binding Protein (RBP) binding specificity and characteristics using data from high-throughput sequencing experiments like RBNS (RNA Bind-n-Seq) and CLIP (CrossLinking and ImmunoPrecipitation).

The package provides tools to:

* For CLIP peak data, calculate K-mer enrichment within binding peak regions, with options for background correction.
* Quantify the Inherent Specificity (IS) and the Mutational Sensitivity (MS) of RBP based on the K-mer enrichment.
* Visualize these specificity and sensitivity metrics.

For detailed information of these metrics, please check our [manuscript](https://www.biorxiv.org/content/10.1101/2025.03.28.646018v2).

## Core Functionality

RBPSpecificity provides a step-by-step approach to analyze RBP binding data:

### 1. K-mer Motif Enrichment (`motifEnrichment()`)

* **Purpose:** This initial step processes your RBP binding peak data (e.g., from CLIP experiments) to calculate the enrichment scores for all possible K-mers of a specified length.
* **Key Features:**
    * Accepts peak data as data frames or GRanges objects.
    * Supports various genome builds (e.g., human, mouse).
    * Includes a robust background correction method using randomly shifted and scrambled genomic sequences to provide more accurate enrichment values.
* **Output:** A data frame containing K-mers and their corresponding normalized enrichment scores, which serves as the input for subsequent specificity and sensitivity analyses.

### 2. Inherent Specificity (IS) Analysis (`returnIS()` & `plotIS()`)

* **Purpose:** After calculating K-mer enrichment, this module allows you to quantify the Inherent Specificity (IS) of an RBP for a chosen K-mer. The IS metric reflects how preferentially an RBP binds to that K-mer compared to all other K-mers.
* **Key Features:**
    * Calculates the IS score based on the relative enrichment of a target motif.
    * `plotIS()` visualizes the overall distribution of K-mer enrichment scores, highlighting the position of the selected motif and its IS value.
* **Output:** A numeric IS score and a ggplot object for visualization.

### 3. Mutational Sensitivity (MS) Analysis (`returnMS()` & `plotMS()`)

* **Purpose:** This module assesses how single nucleotide variations (SNVs) within a specific K-mer affect its binding enrichment, a measure known as Mutational Sensitivity (MS). This analysis helps identify critical nucleotide positions within the motif that are essential for RBP recognition.
* **Key Features:**
    * Calculates the impact of every possible single base change across all positions of a selected motif.
    * `returnMS()` can output a detailed matrix of sensitivity scores for each SNV or a single, averaged MS score.
    * `plotMS()` generates a visual profile of the mutational sensitivity across the motif positions.
* **Output:** An MS matrix or a summary MS score, and a ggplot object for visualization.

## Installation

You can install the development version of `RBPSpecificity` from [GitHub](https://github.com/your_github_username/RBPSpecificity) using the `devtools` package (or the `remotes` package).

First, ensure you have `devtools` installed:

```r
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}
