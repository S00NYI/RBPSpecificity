# RBPSpecificity News

## Version 0.99.0

### New Features
- Initial Bioconductor submission
- **Inherent Specificity (IS)** calculation from k-mer enrichment data
- **Variation Sensitivity (VS)** quantification for motif variants
- **De novo enrichment** from CLIP peak data via `motifEnrichment()`

### Core Functions
- `returnIS()` - Calculate Inherent Specificity for motifs
- `returnVS()` - Calculate Variation Sensitivity matrices or scores
- `motifEnrichment()` - Generate k-mer enrichment from peak data

### Visualization
- `plotIS()` - Score distribution histograms with IS annotation
- `plotVS()` - Variation sensitivity profile visualization

### Data Support
- RBNS enrichment data input
- eCLIP/iCLIP peak data processing
- Multiple genome support (hg38, mm10, etc.)

### Future Implementations
- PWM exploration

