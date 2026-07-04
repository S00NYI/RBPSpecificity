library(RBPSpecificity)

message("=== Starting Manual RBP Specificity and Sensitivity Checks ===")

# 1. Check genome dependency
if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
    stop(
        "BSgenome.Hsapiens.UCSC.hg38 is required to run this script. Install it using:\n",
        "BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')"
    )
}

# 2. Load Peak files from package data (10-column format)
hnrnpc_bed <- system.file("extdata", "ENCODE_eCLIP_HNRNPC_K562_narrowPeak.bed", package = "RBPSpecificity")
eif4g2_bed <- system.file("extdata", "ENCODE_eCLIP_EIF4G2_K562_narrowPeak.bed", package = "RBPSpecificity")
rbfox2_bed <- system.file("extdata", "ENCODE_eCLIP_RBFOX2_K562_narrowPeak.bed", package = "RBPSpecificity")
pcbp2_bed <- system.file("extdata", "ENCODE_eCLIP_PCBP2_HepG2_narrowPeak.bed", package = "RBPSpecificity")

bed_cols <- c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")

hnrnpc_peaks <- read.table(hnrnpc_bed, header = FALSE, sep = "\t")
colnames(hnrnpc_peaks) <- bed_cols

eif4g2_peaks <- read.table(eif4g2_bed, header = FALSE, sep = "\t")
colnames(eif4g2_peaks) <- bed_cols

rbfox2_peaks <- read.table(rbfox2_bed, header = FALSE, sep = "\t")
colnames(rbfox2_peaks) <- bed_cols

pcbp2_peaks <- read.table(pcbp2_bed, header = FALSE, sep = "\t")
colnames(pcbp2_peaks) <- bed_cols

message("Loaded peak datasets:")
message("  HNRNPC: ", nrow(hnrnpc_peaks), " peaks")
message("  EIF4G2: ", nrow(eif4g2_peaks), " peaks")
message("  RBFOX2: ", nrow(rbfox2_peaks), " peaks")
message("  PCBP2 : ", nrow(pcbp2_peaks), " peaks")

# 3. Parameters for motifEnrichment
extension5 <- 25
extension3 <- 0
bootstrap <- 100

message("\n=== Running motifEnrichment (bootstrap = ", bootstrap, ") ===")

# ZOOPS mode
message("ZOOPS Mode:")
message("  HNRNPC...")
hnrnpc_zoops <- motifEnrichment(hnrnpc_peaks, "hg38", K = 5, extension = c(extension5, extension3), method = "zoops", bkg_iter = bootstrap, scramble_bkg = FALSE)
message("  EIF4G2...")
eif4g2_zoops <- motifEnrichment(eif4g2_peaks, "hg38", K = 5, extension = c(extension5, extension3), method = "zoops", bkg_iter = bootstrap, scramble_bkg = FALSE)
message("  RBFOX2...")
rbfox2_zoops <- motifEnrichment(rbfox2_peaks, "hg38", K = 5, extension = c(extension5, extension3), method = "zoops", bkg_iter = bootstrap, scramble_bkg = FALSE)
message("  PCBP2...")
pcbp2_zoops <- motifEnrichment(pcbp2_peaks, "hg38", K = 5, extension = c(extension5, extension3), method = "zoops", bkg_iter = bootstrap, scramble_bkg = FALSE)

# ANR mode
message("ANR Mode:")
message("  HNRNPC...")
hnrnpc_anr <- motifEnrichment(hnrnpc_peaks, "hg38", K = 5, extension = c(extension5, extension3), method = "anr", bkg_iter = bootstrap, scramble_bkg = FALSE)
message("  EIF4G2...")
eif4g2_anr <- motifEnrichment(eif4g2_peaks, "hg38", K = 5, extension = c(extension5, extension3), method = "anr", bkg_iter = bootstrap, scramble_bkg = FALSE)
message("  RBFOX2...")
rbfox2_anr <- motifEnrichment(rbfox2_peaks, "hg38", K = 5, extension = c(extension5, extension3), method = "anr", bkg_iter = bootstrap, scramble_bkg = FALSE)
message("  PCBP2...")
pcbp2_anr <- motifEnrichment(pcbp2_peaks, "hg38", K = 5, extension = c(extension5, extension3), method = "anr", bkg_iter = bootstrap, scramble_bkg = FALSE)

# Calculate Specificity and Sensitivity metrics
cs_hnrnpc_zoops <- returnSpecificity(hnrnpc_zoops)
cvs_hnrnpc_zoops <- returnSensitivity(hnrnpc_zoops, output_type = "number")
cs_eif4g2_zoops <- returnSpecificity(eif4g2_zoops)
cvs_eif4g2_zoops <- returnSensitivity(eif4g2_zoops, output_type = "number")
cs_rbfox2_zoops <- returnSpecificity(rbfox2_zoops)
cvs_rbfox2_zoops <- returnSensitivity(rbfox2_zoops, output_type = "number")
cs_pcbp2_zoops <- returnSpecificity(pcbp2_zoops)
cvs_pcbp2_zoops <- returnSensitivity(pcbp2_zoops, output_type = "number")

cs_hnrnpc_anr <- returnSpecificity(hnrnpc_anr)
cvs_hnrnpc_anr <- returnSensitivity(hnrnpc_anr, output_type = "number")
cs_eif4g2_anr <- returnSpecificity(eif4g2_anr)
cvs_eif4g2_anr <- returnSensitivity(eif4g2_anr, output_type = "number")
cs_rbfox2_anr <- returnSpecificity(rbfox2_anr)
cvs_rbfox2_anr <- returnSensitivity(rbfox2_anr, output_type = "number")
cs_pcbp2_anr <- returnSpecificity(pcbp2_anr)
cvs_pcbp2_anr <- returnSensitivity(pcbp2_anr, output_type = "number")

# 4. Load in vitro RBNS Data and compute metrics
message("\n=== Calculating In Vitro Specificity/Sensitivity from RBNS ===")
rbns_file <- system.file("extdata", "RBNS_normalized_5mer.csv", package = "RBPSpecificity")
rbns_data <- read.csv(rbns_file, header = TRUE)

rbns_metrics <- list()
for (rbp in c("EIF4G2", "HNRNPC", "RBFOX2", "PCBP2")) {
    rbns_enrichment <- data.frame(
        MOTIF = rbns_data$Motif,
        Score = rbns_data[[rbp]]
    )
    is_val <- returnSpecificity(rbns_enrichment)
    vs_val <- returnSensitivity(rbns_enrichment, output_type = "number")
    rbns_metrics[[rbp]] <- list(Specificity = is_val, Sensitivity = vs_val)
}

# 5. Display Comparative Results
message("\n=== COMPARISON SUMMARY ===")
message("Structural Context Independent (Direct motif binding, should show strong correlation/matching trends):")
message("  HNRNPC:")
message("    In Vitro (RBNS)   - Spec: ", round(rbns_metrics[["HNRNPC"]]$Specificity, 2), " | Sens: ", round(rbns_metrics[["HNRNPC"]]$Sensitivity, 4))
message("    Cellular (ZOOPS)  - Spec: ", round(cs_hnrnpc_zoops, 2), " | Sens: ", round(cvs_hnrnpc_zoops, 4))
message("    Cellular (ANR)    - Spec: ", round(cs_hnrnpc_anr, 2), " | Sens: ", round(cvs_hnrnpc_anr, 4))
message("  PCBP2:")
message("    In Vitro (RBNS)   - Spec: ", round(rbns_metrics[["PCBP2"]]$Specificity, 2), " | Sens: ", round(rbns_metrics[["PCBP2"]]$Sensitivity, 4))
message("    Cellular (ZOOPS)  - Spec: ", round(cs_pcbp2_zoops, 2), " | Sens: ", round(cvs_pcbp2_zoops, 4))
message("    Cellular (ANR)    - Spec: ", round(cs_pcbp2_anr, 2), " | Sens: ", round(cvs_pcbp2_anr, 4))

message("\nStructural Context Dependent (Divergence or low correlation expected between in vitro and cellular contexts):")
message("  RBFOX2:")
message("    In Vitro (RBNS)   - Spec: ", round(rbns_metrics[["RBFOX2"]]$Specificity, 2), " | Sens: ", round(rbns_metrics[["RBFOX2"]]$Sensitivity, 4))
message("    Cellular (ZOOPS)  - Spec: ", round(cs_rbfox2_zoops, 2), " | Sens: ", round(cvs_rbfox2_zoops, 4))
message("    Cellular (ANR)    - Spec: ", round(cs_rbfox2_anr, 2), " | Sens: ", round(cvs_rbfox2_anr, 4))
message("  EIF4G2:")
message("    In Vitro (RBNS)   - Spec: ", round(rbns_metrics[["EIF4G2"]]$Specificity, 2), " | Sens: ", round(rbns_metrics[["EIF4G2"]]$Sensitivity, 4))
message("    Cellular (ZOOPS)  - Spec: ", round(cs_eif4g2_zoops, 2), " | Sens: ", round(cvs_eif4g2_zoops, 4))
message("    Cellular (ANR)    - Spec: ", round(cs_eif4g2_anr, 2), " | Sens: ", round(cvs_eif4g2_anr, 4))

message("\nManual check complete.")
