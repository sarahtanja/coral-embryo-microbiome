# Example: Mantel Test Workflow for Coral Embryo Project

## Quick Start Example

This document provides a practical example of running a Mantel test using the actual file paths and data structure in this repository.

## Prerequisites

### Required Data Files

**Currently available in this repository:**
- ✅ Microbiome ASV table: `salipante/241121_StonyCoral/270x200/1st_collapse/feature-table.tsv`
- ✅ Metadata: `metadata/metadata.csv`
- ❌ RNA-seq data: **NOT YET AVAILABLE**

### What You Need to Add

To run the complete Mantel test, you need to add RNA-seq count data from the same samples. The RNA-seq file should be:
- A CSV or TSV file with genes/transcripts as columns and samples as rows (or vice versa)
- Sample IDs matching those in `metadata/metadata.csv`
- Raw or normalized count data

Example RNA-seq file structure:
```
sample_id,Gene1,Gene2,Gene3,...
101112C14,1245,892,456,...
101112C4,1189,934,512,...
...
```

## Step-by-Step Example

### 1. Setup R Environment

```r
# Set working directory to the repository root
setwd("/path/to/coral-embryo-microbiome")

# Load required packages
library(vegan)
library(tidyverse)
library(compositions)
library(DESeq2)

# Set seed for reproducibility
set.seed(123)
```

### 2. Import Microbiome Data

```r
# Read the collapsed feature table from QIIME2
asv_table <- read_tsv(
  "salipante/241121_StonyCoral/270x200/1st_collapse/feature-table.tsv",
  skip = 1,
  show_col_types = FALSE
)

# Transpose so samples are rows
asv_matrix <- asv_table %>%
  column_to_rownames("#OTU ID") %>%
  t() %>%
  as.data.frame()

cat("Microbiome data:", nrow(asv_matrix), "samples x", 
    ncol(asv_matrix), "features\n")

# Preview
asv_matrix[1:5, 1:5]
```

Expected output:
```
Microbiome data: 63 samples x 127 features

# A sample preview showing counts
```

### 3. Filter and Transform Microbiome Data

```r
# Filter rare features (present in <10% of samples)
prevalence_threshold <- 0.1 * nrow(asv_matrix)
asv_filtered <- asv_matrix[, colSums(asv_matrix > 0) >= prevalence_threshold]

cat("After filtering:", ncol(asv_filtered), "features retained\n")

# Apply CLR transformation for compositional data
asv_clr <- clr(asv_filtered + 0.5)  # Add pseudocount for zeros

# Calculate Aitchison distance
microbiome_dist <- dist(asv_clr, method = "euclidean")

cat("Microbiome distance matrix:", attr(microbiome_dist, "Size"), "samples\n")
```

### 4. Import RNA-seq Data (Placeholder)

```r
# PLACEHOLDER: Replace with your actual RNA-seq file path
rnaseq_file <- "data/rnaseq_counts.csv"

# Check if file exists
if (!file.exists(rnaseq_file)) {
  stop("RNA-seq data not found. Please add your RNA-seq count matrix to: ", rnaseq_file)
}

# Import RNA-seq counts
rnaseq_counts <- read.csv(rnaseq_file, row.names = 1)

# Ensure samples are rows (transpose if needed)
if (nrow(rnaseq_counts) > ncol(rnaseq_counts)) {
  rnaseq_counts <- t(rnaseq_counts) %>% as.data.frame()
}

cat("RNA-seq data:", nrow(rnaseq_counts), "samples x",
    ncol(rnaseq_counts), "genes\n")
```

### 5. Transform RNA-seq Data

```r
# Use DESeq2 variance stabilizing transformation
dds <- DESeqDataSetFromMatrix(
  countData = t(rnaseq_counts),  # DESeq2 expects genes as rows
  colData = data.frame(sample = rownames(rnaseq_counts)),
  design = ~ 1
)

# Apply VST
vst_data <- vst(dds, blind = TRUE)
rnaseq_transformed <- t(assay(vst_data)) %>% as.data.frame()

# Calculate Euclidean distance
rnaseq_dist <- dist(rnaseq_transformed, method = "euclidean")

cat("RNA-seq distance matrix:", attr(rnaseq_dist, "Size"), "samples\n")
```

### 6. Match Samples Between Datasets

```r
# Find common samples
common_samples <- intersect(labels(microbiome_dist), labels(rnaseq_dist))

if (length(common_samples) == 0) {
  stop("No common samples found between datasets!")
}

cat("Found", length(common_samples), "matched samples\n")

# Subset distance matrices
microbiome_dist_matched <- as.dist(
  as.matrix(microbiome_dist)[common_samples, common_samples]
)

rnaseq_dist_matched <- as.dist(
  as.matrix(rnaseq_dist)[common_samples, common_samples]
)

# Verify sample order
stopifnot(all(labels(microbiome_dist_matched) == labels(rnaseq_dist_matched)))
```

### 7. Run Mantel Test

```r
# Perform Mantel test
mantel_result <- mantel(
  xdis = microbiome_dist_matched,
  ydis = rnaseq_dist_matched,
  method = "spearman",
  permutations = 9999
)

# Print results
print(mantel_result)

# Extract statistics
cat("\n=== Results Summary ===\n")
cat("Mantel statistic (r):", round(mantel_result$statistic, 4), "\n")
cat("P-value:", round(mantel_result$signif, 4), "\n")
cat("Significance:", ifelse(mantel_result$signif < 0.05, "YES", "NO"), "\n")
```

### 8. Visualize Results

```r
# Create scatterplot
library(ggplot2)

plot_data <- data.frame(
  Microbiome = as.vector(microbiome_dist_matched),
  RNAseq = as.vector(rnaseq_dist_matched)
)

mantel_plot <- ggplot(plot_data, aes(x = Microbiome, y = RNAseq)) +
  geom_point(alpha = 0.3, size = 2, color = "steelblue") +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(
    title = "Mantel Test: Microbiome vs. RNA-seq",
    subtitle = paste0("Spearman r = ", round(mantel_result$statistic, 3),
                     ", p = ", round(mantel_result$signif, 4)),
    x = "Microbiome Distance (Aitchison)",
    y = "RNA-seq Distance (Euclidean)",
    caption = paste(length(common_samples), "samples")
  ) +
  theme_bw()

# Save plot
ggsave("output/mantel_test_result.png", 
       mantel_plot, 
       width = 8, 
       height = 6, 
       dpi = 300)

print(mantel_plot)
```

### 9. Partial Mantel Test (Control for Developmental Stage)

```r
# Import metadata
metadata <- read.csv("metadata/metadata.csv", row.names = 1)

# Match to common samples
metadata_matched <- metadata[common_samples, ]

# Create distance matrix for developmental timepoint
timepoint_matrix <- model.matrix(~ hpf - 1, data = metadata_matched)
timepoint_dist <- dist(timepoint_matrix, method = "euclidean")

# Run partial Mantel test
partial_result <- mantel.partial(
  xdis = microbiome_dist_matched,
  ydis = rnaseq_dist_matched,
  zdis = timepoint_dist,
  method = "spearman",
  permutations = 9999
)

print(partial_result)

cat("\n=== Partial Mantel (controlling for timepoint) ===\n")
cat("Mantel statistic:", round(partial_result$statistic, 4), "\n")
cat("P-value:", round(partial_result$signif, 4), "\n")
```

### 10. Stratified Analysis by Treatment

```r
# Analyze each leachate treatment separately
treatments <- c("control", "low", "mid", "high")
stratified_results <- list()

for (trt in treatments) {
  # Get samples for this treatment
  trt_samples <- rownames(metadata_matched)[metadata_matched$leachate == trt]
  trt_samples <- intersect(trt_samples, common_samples)
  
  if (length(trt_samples) < 4) {
    cat("Skipping", trt, "- insufficient samples\n")
    next
  }
  
  # Subset distance matrices
  mb_trt <- as.dist(
    as.matrix(microbiome_dist_matched)[trt_samples, trt_samples]
  )
  rna_trt <- as.dist(
    as.matrix(rnaseq_dist_matched)[trt_samples, trt_samples]
  )
  
  # Run Mantel test
  result_trt <- mantel(mb_trt, rna_trt, 
                       method = "spearman", 
                       permutations = 999)
  
  stratified_results[[trt]] <- list(
    n = length(trt_samples),
    r = result_trt$statistic,
    p = result_trt$signif
  )
  
  cat(trt, "treatment: n =", length(trt_samples), 
      ", r =", round(result_trt$statistic, 3),
      ", p =", round(result_trt$signif, 4), "\n")
}
```

### 11. Save Results

```r
# Compile results
results_df <- data.frame(
  Analysis = "Mantel Test",
  N_Samples = length(common_samples),
  Mantel_r = mantel_result$statistic,
  P_value = mantel_result$signif,
  Method = "Spearman",
  Permutations = 9999,
  Date = Sys.Date()
)

# Save to CSV
write.csv(results_df, 
          "output/mantel_test_results.csv", 
          row.names = FALSE)

cat("\nResults saved to: output/mantel_test_results.csv\n")
```

## Expected Sample Sizes

Based on the metadata in this repository:

| Category | Levels | Samples per level |
|----------|--------|------------------|
| Timepoint | 4, 9, 14 hpf | ~21 per timepoint |
| Leachate | Control, Low, Mid, High | ~15-16 per treatment |
| Total | - | 63 samples |

**Note:** For stratified analyses by treatment, you'll have ~15 samples per group, which is adequate for Mantel tests.

## Troubleshooting

### Issue: "RNA-seq data not found"
**Solution:** Create your RNA-seq count matrix and save it to `data/rnaseq_counts.csv` with sample IDs matching those in `metadata/metadata.csv`.

### Issue: "No common samples found"
**Solution:** Check that your RNA-seq sample IDs exactly match the microbiome sample IDs (e.g., "101112C14", "101112C4", etc.).

### Issue: Sample order mismatch
**Solution:** The script includes checks (`stopifnot()`) to ensure sample order matches. If this fails, verify your data import.

### Issue: Non-significant result
**Possible reasons:**
1. True lack of correlation (microbiome and transcriptome not related)
2. Insufficient sample size (need more matched samples)
3. High variability within groups
4. Confounding variables (try partial Mantel test)

## Next Steps

After running the Mantel test:

1. **If significant:**
   - Identify which microbial taxa correlate with specific genes
   - Use Procrustes analysis to compare ordinations
   - Explore functional pathways (PICRUSt2 for microbiome, GO/KEGG for transcriptome)

2. **If not significant:**
   - Try partial Mantel controlling for timepoint or treatment
   - Analyze subsets (e.g., separate analyses per timepoint)
   - Consider other multi-omics methods (MOFA, DIABLO)

## Complete Script

A complete executable version of this workflow is available in:
- **[mantel_test_analysis.qmd](mantel_test_analysis.qmd)**

For detailed explanations and theory:
- **[MANTEL_TEST_GUIDE.md](MANTEL_TEST_GUIDE.md)**

For quick reference:
- **[MANTEL_QuickRef.md](MANTEL_QuickRef.md)**
