################################################################################
## RBPSpecificity Alpha Build Test Script
## Written by Soon Yi
## Last Updated: 2025-12-23
##
## Purpose:
##   This script tests the core functionalities of the RBPSpecificity package.
##   It serves as the basis for the package vignette.
##
## Data Files (in inst/extdata/):
##   - RBNS_normalized_5mer.csv: Pre-computed RBNS enrichment scores for EIF4G2 & HNRNPC
##   - eCLIP_Peaks_HNRNPC.bed: 516 eCLIP peaks for HNRNPC (K562)
##   - eCLIP_Peaks_EIF4G2.bed: 4628 eCLIP peaks for EIF4G2 (K562)
################################################################################

# ==============================================================================
# Section 0: Setup
# ==============================================================================
# For development: load from source instead of installed package
# library(RBPSpecificity)  # Use this after package is installed
devtools::load_all(".")    # Use this during development

library(ggplot2)

# Define K-mer size for all analyses
K <- 5

# Path to data files - uses inst/extdata/ during development
# After installation, use system.file("extdata", ..., package = "RBPSpecificity")
data_dir <- "inst/extdata/"

# ==============================================================================
# Section 1: Workflow A - Pre-computed Enrichment Data (RBNS)
# ==============================================================================
# Goal: Demonstrate returnIS(), returnMS(), plotIS(), and plotMS() using
#       pre-computed RBNS scores without requiring genome access.

message("\n========== Section 1: RBNS Analysis ==========\n")

# --- Load RBNS data ---
rbns_file <- file.path(data_dir, "RBNS_normalized_5mer.csv")
rbns_data <- read.csv(rbns_file, stringsAsFactors = FALSE)

# Check loaded data
message("Loaded RBNS data with ", nrow(rbns_data), " motifs and ", ncol(rbns_data) - 1, " RBPs.")
message("RBPs: ", paste(colnames(rbns_data)[-1], collapse = ", "))

# --- Prepare data for a single RBP: HNRNPC ---
# The package expects a data.frame with 'MOTIF' and 'Score' columns.
hnrnpc_rbns <- data.frame(
  MOTIF = rbns_data$Motif,
  Score = rbns_data$HNRNPC
)

# --- Calculate Inherent Specificity (IS) ---
message("\n--- HNRNPC RBNS: Inherent Specificity ---")
is_value_hnrnpc <- returnIS(hnrnpc_rbns)
message("IS for top motif: ", round(is_value_hnrnpc, 2))

# Get IS for ALL motifs
is_all_hnrnpc <- returnIS(hnrnpc_rbns, return_type = "all")
message("IS calculated for all ", nrow(is_all_hnrnpc), " motifs.")
print(head(is_all_hnrnpc[order(-is_all_hnrnpc$IS), ], 10)) # Top 10

# --- Calculate Mutational Sensitivity (MS) ---
message("\n--- HNRNPC RBNS: Mutational Sensitivity ---")
ms_matrix_hnrnpc <- returnMS(hnrnpc_rbns, output_type = "matrix")
message("MS matrix for top motif:")
print(ms_matrix_hnrnpc)

ms_value_hnrnpc <- returnMS(hnrnpc_rbns, output_type = "number")
message("Average MS score: ", round(ms_value_hnrnpc, 4))

# --- Visualization ---
message("\n--- HNRNPC RBNS: Visualization ---")
p1 <- plotIS(hnrnpc_rbns)
print(p1)

p2 <- plotMS(hnrnpc_rbns)
print(p2)

# --- Repeat for EIF4G2 ---
message("\n--- EIF4G2 RBNS: Metrics ---")
eif4g2_rbns <- data.frame(
  MOTIF = rbns_data$Motif,
  Score = rbns_data$EIF4G2
)
is_value_eif4g2 <- returnIS(eif4g2_rbns)
ms_value_eif4g2 <- returnMS(eif4g2_rbns, output_type = "number")
message("EIF4G2 IS: ", round(is_value_eif4g2, 2), ", MS: ", round(ms_value_eif4g2, 4))

p3 <- plotIS(eif4g2_rbns)
print(p3)
p4 <- plotMS(eif4g2_rbns)
print(p4)

# ==============================================================================
# Section 2: Workflow B - De Novo Enrichment from Peaks (eCLIP)
# ==============================================================================
# Goal: Demonstrate the full motifEnrichment() pipeline using eCLIP peaks.
# NOTE: This requires BSgenome.Hsapiens.UCSC.hg38 to be installed.

message("\n========== Section 2: eCLIP Peak Analysis ==========\n")

# --- Load eCLIP peak data ---
# The BED files are in narrowPeak format (10 columns)
hnrnpc_peaks_file <- file.path(data_dir, "eCLIP_Peaks_HNRNPC.bed")
hnrnpc_peaks <- read.table(hnrnpc_peaks_file, header = FALSE)
colnames(hnrnpc_peaks) <- c("chr", "start", "end", "name", "score", "strand",
                            "signalValue", "pValue", "qValue", "peak")

message("Loaded ", nrow(hnrnpc_peaks), " HNRNPC eCLIP peaks.")

# --- Run motifEnrichment ---
# NOTE: Bkg_number is set low (10) for speed; use 100+ for real analysis.
message("\nRunning motifEnrichment() for HNRNPC peaks...")
message("(Using Bkg_number = 10 for demo; increase for production runs)")

hnrnpc_enrichment <- motifEnrichment(
  peak_data = hnrnpc_peaks,
  species_or_build = "hg38",
  K = K,
  extension = 25,
  Bkg_number = 10,   # Low for demo speed
  Bkg_dist = 500,
  max_shift_dist = 1000,
  log_transform = FALSE
)

message("Enrichment calculation complete.")
print(head(hnrnpc_enrichment[order(-hnrnpc_enrichment$Score), ], 10))

# --- Calculate IS/MS from enrichment ---
message("\n--- HNRNPC eCLIP: Metrics ---")
is_eclip_hnrnpc <- returnIS(hnrnpc_enrichment)
ms_eclip_hnrnpc <- returnMS(hnrnpc_enrichment, output_type = "number")
message("eCLIP IS: ", round(is_eclip_hnrnpc, 2), ", MS: ", round(ms_eclip_hnrnpc, 4))

# --- Visualize ---
p5 <- plotIS(hnrnpc_enrichment)
print(p5)
p6 <- plotMS(hnrnpc_enrichment)
print(p6)

# --- Repeat for EIF4G2 (using fewer peaks for speed) ---
message("\n--- EIF4G2 eCLIP: Metrics (subsampled) ---")
eif4g2_peaks_file <- file.path(data_dir, "eCLIP_Peaks_EIF4G2.bed")
eif4g2_peaks <- read.table(eif4g2_peaks_file, header = FALSE)
colnames(eif4g2_peaks) <- c("chr", "start", "end", "name", "score", "strand",
                            "signalValue", "pValue", "qValue", "peak")
message("Loaded ", nrow(eif4g2_peaks), " EIF4G2 eCLIP peaks.")

# Subsample to 500 peaks for speed
set.seed(42)
eif4g2_peaks_sub <- eif4g2_peaks[sample(nrow(eif4g2_peaks), 500), ]
message("Subsampled to 500 peaks for demo.")

eif4g2_enrichment <- motifEnrichment(
  peak_data = eif4g2_peaks_sub,
  species_or_build = "hg38",
  K = K,
  extension = 25,
  Bkg_number = 10,
  Bkg_dist = 500,
  max_shift_dist = 1000,
  log_transform = FALSE
)

is_eclip_eif4g2 <- returnIS(eif4g2_enrichment)
ms_eclip_eif4g2 <- returnMS(eif4g2_enrichment, output_type = "number")
message("EIF4G2 eCLIP IS: ", round(is_eclip_eif4g2, 2), ", MS: ", round(ms_eclip_eif4g2, 4))

# ==============================================================================
# Section 3: Comparing Metrics Across Multiple RBPs
# ==============================================================================
# Goal: Show how to loop through RBPs in a dataset and build a summary table.

message("\n========== Section 3: Multi-RBP Comparison ==========\n")

# Get list of RBPs from RBNS data column names
rbps <- colnames(rbns_data)[-1] # Exclude 'Motif'
message("Analyzing RBPs: ", paste(rbps, collapse = ", "))

# Prepare results data.frame
results <- data.frame(
  RBP = rbps,
  IS = NA_real_,
  MS = NA_real_
)

for (i in seq_along(rbps)) {
  rbp <- rbps[i]
  temp_df <- data.frame(
    MOTIF = rbns_data$Motif,
    Score = rbns_data[[rbp]]
  )
  results$IS[i] <- returnIS(temp_df)
  results$MS[i] <- returnMS(temp_df, output_type = "number")
}

message("Summary Table:")
print(results)

# ==============================================================================
# Section 4: RBNS vs eCLIP Comparison
# ==============================================================================
# Goal: Compare inherent (RBNS) vs. apparent (eCLIP) specificity for HNRNPC.

message("\n========== Section 4: RBNS vs eCLIP Comparison ==========\n")

comparison_df <- data.frame(
  RBP = "HNRNPC",
  RBNS_IS = is_value_hnrnpc,
  eCLIP_IS = is_eclip_hnrnpc,
  RBNS_MS = ms_value_hnrnpc,
  eCLIP_MS = ms_eclip_hnrnpc
)

message("Comparison of RBNS (Inherent) vs eCLIP (Apparent):")
print(comparison_df)

# ==============================================================================
# Done!
# ==============================================================================
message("\n========== Alpha Build Test Complete ==========")
message("All core functions executed successfully.")
