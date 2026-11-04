# RBPSpecificity News

## Version 1.1.1 (development)

### New Features
- `motifEnrichment()` now works with any genome that has a BSgenome package.
  `species_or_build` accepts a `BSgenome` object, or the name of any installed
  BSgenome build (e.g. `"danRer11"`) or package
  (e.g. `"BSgenome.Drerio.UCSC.danRer11"`). The existing human and mouse
  shortcuts (`"hg38"`, `"mm39"`, etc.) still work.
- New `chromosomes` argument in `motifEnrichment()`: `"standard"` (default),
  `"all"`, or a vector of chromosome names.
- `motifEnrichment()` prints a table of how many peaks each filtering step
  removed, and attaches it to the result as `attr(result, "filter_summary")`.

### Changes
- Peak chromosome names are matched to the genome's naming style before
  filtering, so Ensembl-style peaks (`"1"`, `"MT"`) work with UCSC genomes.
  Previously these peaks were dropped.
- Standard chromosomes are now taken from the selected genome's organism,
  instead of being guessed before the genome was known.
- If no peaks remain after filtering, the error now shows the summary table
  and a hint about the likely cause.

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

