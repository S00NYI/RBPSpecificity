# R script to document and retrieve the raw external datasets stored in inst/extdata

# This script documents how the raw data files under inst/extdata/ were obtained
# and processed, satisfying the Bioconductor guidelines.

# ==============================================================================
# 1. ENCODE eCLIP Peak Datasets
# ==============================================================================
# The package includes raw eCLIP narrowPeak BED files for four RNA-binding
# proteins (EIF4G2, HNRNPC, PCBP2, and RBFOX2) retrieved from the ENCODE project
# portal (https://www.encodeproject.org/).
#
# Reference:
# Luo, Y., Hitz, B.C., Gabdank, I., Hilton, J.A., Kagda, M.S., Lam, B.,
# Myers, Z., Sud, P., Jou, J., Lin, K., et al. (2020). New developments
# on the Encyclopedia of DNA Elements (ENCODE) data portal. Nucleic Acids Res.
# 48, D882–D889.
#
# Files and Accessions:
# - ENCODE_eCLIP_EIF4G2_K562_narrowPeak.bed
#   * Experiment Accession: ENCSR307YIW
#   * Source URL: https://www.encodeproject.org/experiments/ENCSR307YIW/
# - ENCODE_eCLIP_HNRNPC_K562_narrowPeak.bed
#   * File Accession: ENCFF167CDB
#   * Source URL: https://www.encodeproject.org/files/ENCFF167CDB/
# - ENCODE_eCLIP_PCBP2_HepG2_narrowPeak.bed
#   * File Accession: ENCFF642GNE
#   * Source URL: https://www.encodeproject.org/files/ENCFF642GNE/
# - ENCODE_eCLIP_RBFOX2_K562_narrowPeak.bed
#   * File Accession: ENCFF206RIM
#   * Source URL: https://www.encodeproject.org/files/ENCFF206RIM/
#
# Preprocessing:
# The BED files were downloaded and trimmed to include only peak coordinates.
#
# Example R code to download HNRNPC peaks:
# download.file("https://www.encodeproject.org/files/ENCFF167CDB/@@download/ENCFF167CDB.bed.gz",
#               destfile = "ENCODE_eCLIP_HNRNPC_K562_narrowPeak.bed.gz")

# ==============================================================================
# 2. RNA Bind-n-Seq (RBNS) Dataset
# ==============================================================================
# File: RBNS_normalized_5mer.csv
#
# A comma-separated values (CSV) file containing normalized 5-mer enrichment
# scores (using U instead of T) for four RBPs: EIF4G2, HNRNPC, RBFOX2, and PCBP2.
#
# Reference:
# Lambert, N., Robertson, A., Jangi, M., McGeary, S., Sharp, P.A., and Burge,
# C.B. (2014). RNA Bind-n-Seq: quantitative assessment of the sequence and
# structural binding specificity of RNA binding proteins. Mol. Cell 54, 887–900.
#
# Source:
# The raw counts were obtained from the ENCODE project portal (https://www.encodeproject.org/).
# For each K-mer enrichment dataset, the RBP concentration generating the
# highest R-value was selected. The collected RBNS R scores were feature-scaled
# to [1, e] followed by natural log transformation to normalize the enrichment
# scores to a [0, 1] range.
