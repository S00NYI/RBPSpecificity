# RBPSpecificity News

## Version 0.99.0

### New Features
- Initial Bioconductor submission
- **Inherent Specificity (IS)** calculation from k-mer enrichment data
- **Mutational Sensitivity (MS)** quantification for motif variants
- **De novo enrichment** from CLIP peak data via `motifEnrichment()`

### Core Functions
- `returnIS()` - Calculate Inherent Specificity for motifs
- `returnMS()` - Calculate Mutational Sensitivity matrices or scores
- `motifEnrichment()` - Generate k-mer enrichment from peak data

### Visualization
- `plotIS()` - Score distribution histograms with IS annotation
- `plotMS()` - Mutational sensitivity profile visualization

### Data Support
- RBNS enrichment data input
- eCLIP/iCLIP peak data processing
- Multiple genome support (hg38, mm10, etc.)

### Future Implementations
- PWM exploration
