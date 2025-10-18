# Mantel Test for Microbiome-Transcriptome Correlation Analysis

## Overview

The **Mantel test** is a statistical method used to assess the correlation between two distance or dissimilarity matrices. In the context of coral embryo research, this test can be used to determine whether patterns in the microbiome composition (ASV/16S rRNA data) are correlated with patterns in host gene expression (RNA-seq data).

This guide explains how to perform a Mantel test to compare your microbiome ASV matrix with RNA-seq expression data, following the vegan package implementation in R.

## What is a Mantel Test?

The Mantel test evaluates the relationship between two distance/dissimilarity matrices by:

1. Converting both datasets into distance matrices
2. Calculating a correlation coefficient (typically Pearson or Spearman) between corresponding elements of the two matrices
3. Assessing statistical significance through permutation testing

### When to Use a Mantel Test

Use a Mantel test when you want to answer questions like:
- Is the similarity in microbial community composition between samples correlated with similarity in host gene expression?
- Do samples that are similar in their microbiome also have similar transcriptome profiles?
- Are changes in the microbiome associated with changes in host gene expression patterns?

### Key Considerations

**Strengths:**
- Non-parametric (no distributional assumptions)
- Handles complex multivariate relationships
- Appropriate for compositional and count data (after transformation)
- Well-suited for ecological and microbiome studies

**Limitations:**
- Tests correlation, not causation
- Assumes independence of observations (violations can inflate Type I error)
- Statistical power depends on sample size
- May have inflated Type I error with spatial autocorrelation

## Required Data

### 1. Microbiome Data (ASV/16S rRNA)

You need a feature table (ASV table) where:
- **Rows** = samples
- **Columns** = ASVs/taxa
- **Values** = abundance counts or relative abundances

**Example structure:**
```
          ASV_001  ASV_002  ASV_003  ...
Sample_1    142      58       23     ...
Sample_2     89      72       41     ...
Sample_3    156      45       18     ...
```

**Location in this repository:**
- QIIME2 feature tables: `salipante/241121_StonyCoral/270x200/feature-table.tsv`
- Or collapsed taxonomy tables: `salipante/241121_StonyCoral/270x200/1st_collapse/feature-table.tsv`

### 2. RNA-seq Data (Gene Expression)

You need a gene expression matrix where:
- **Rows** = samples (matching the microbiome samples)
- **Columns** = genes/transcripts
- **Values** = expression counts (raw or normalized)

**Example structure:**
```
          Gene_1  Gene_2  Gene_3  ...
Sample_1   1245    892    456    ...
Sample_2   1189    934    512    ...
Sample_3   1302    867    423    ...
```

**Important:** Sample IDs must match between the microbiome and RNA-seq datasets!

## Step-by-Step Analysis

### Step 1: Install and Load Required Packages

```r
# Install packages if needed
if (!requireNamespace("vegan", quietly = TRUE)) {
  install.packages("vegan")
}
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}
if (!requireNamespace("compositions", quietly = TRUE)) {
  install.packages("compositions")
}
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}

# Load libraries
library(vegan)
library(tidyverse)
library(compositions)
library(DESeq2)
```

### Step 2: Import and Prepare Microbiome Data

```r
# Import ASV table from QIIME2 export
asv_table <- read_tsv("salipante/241121_StonyCoral/270x200/1st_collapse/feature-table.tsv",
                      skip = 1) # Skip the first line if it's a comment

# Set sample IDs as rownames and remove metadata columns
asv_matrix <- asv_table %>%
  column_to_rownames("#OTU ID") %>%
  t() %>%  # Transpose so samples are rows
  as.data.frame()

# Remove rare features (optional but recommended)
# Keep ASVs present in at least 10% of samples
prevalence_threshold <- 0.1 * nrow(asv_matrix)
asv_matrix_filtered <- asv_matrix[, colSums(asv_matrix > 0) >= prevalence_threshold]

# Apply CLR transformation for compositional data
# Add a small pseudocount to handle zeros
pseudocount <- 0.5
asv_clr <- clr(asv_matrix_filtered + pseudocount)

# Calculate Aitchison distance (Euclidean distance on CLR-transformed data)
microbiome_dist <- dist(asv_clr, method = "euclidean")

# Alternatively, use Bray-Curtis dissimilarity on relative abundances
# asv_rel_abund <- decostand(asv_matrix_filtered, method = "total")
# microbiome_dist <- vegdist(asv_rel_abund, method = "bray")
```

### Step 3: Import and Prepare RNA-seq Data

```r
# Import RNA-seq count matrix
# Assuming you have a file with genes as rows and samples as columns
rnaseq_counts <- read.csv("path/to/rnaseq_counts.csv", row.names = 1)

# Transpose if needed so samples are rows
if (nrow(rnaseq_counts) > ncol(rnaseq_counts)) {
  rnaseq_counts <- t(rnaseq_counts)
}

# Normalize RNA-seq data using DESeq2 (recommended)
# This accounts for library size and compositional effects
dds <- DESeqDataSetFromMatrix(
  countData = t(rnaseq_counts),  # DESeq2 expects genes as rows
  colData = data.frame(sample = rownames(rnaseq_counts)),
  design = ~ 1  # No design for simple normalization
)

# Apply variance stabilizing transformation
vst_data <- vst(dds, blind = TRUE)
rnaseq_transformed <- t(assay(vst_data))  # Transpose back to samples as rows

# Alternative: Use log2(CPM + 1) transformation
# cpm <- sweep(rnaseq_counts, 1, rowSums(rnaseq_counts), FUN = "/") * 1e6
# rnaseq_transformed <- log2(cpm + 1)

# Calculate Euclidean distance on transformed expression data
rnaseq_dist <- dist(rnaseq_transformed, method = "euclidean")

# Alternatively, use correlation-based distance
# cor_matrix <- cor(t(rnaseq_transformed), method = "spearman")
# rnaseq_dist <- as.dist(1 - cor_matrix)
```

### Step 4: Match Samples Between Datasets

```r
# Ensure both distance matrices have the same samples in the same order
common_samples <- intersect(labels(microbiome_dist), labels(rnaseq_dist))

# Subset distance matrices to common samples
microbiome_dist_matched <- as.dist(as.matrix(microbiome_dist)[common_samples, common_samples])
rnaseq_dist_matched <- as.dist(as.matrix(rnaseq_dist)[common_samples, common_samples])

# Verify sample order matches
stopifnot(all(labels(microbiome_dist_matched) == labels(rnaseq_dist_matched)))

cat("Number of matched samples:", length(common_samples), "\n")
```

### Step 5: Run the Mantel Test

```r
# Perform Mantel test
# method: "pearson" for linear relationships, "spearman" for monotonic relationships
# permutations: number of permutations for significance testing (typically 999 or 9999)

mantel_result <- mantel(microbiome_dist_matched, 
                        rnaseq_dist_matched,
                        method = "spearman",
                        permutations = 9999)

# Print results
print(mantel_result)

# Extract key statistics
cat("\nMantel statistic (r):", mantel_result$statistic, "\n")
cat("Significance (p-value):", mantel_result$signif, "\n")
```

### Step 6: Interpret Results

The Mantel test returns:

- **Mantel statistic (r)**: Correlation coefficient between the two distance matrices
  - Range: -1 to +1
  - Positive values: samples similar in microbiome tend to be similar in transcriptome
  - Negative values: samples similar in microbiome tend to be dissimilar in transcriptome
  - Values near 0: no correlation
  
- **p-value**: Statistical significance based on permutation test
  - p < 0.05: significant correlation
  - p ≥ 0.05: no significant correlation

**Example interpretation:**
- r = 0.35, p = 0.001: "There is a significant positive correlation between microbiome composition and host gene expression patterns (Mantel r = 0.35, p < 0.01). Samples with similar microbial communities tend to have similar transcriptome profiles."

### Step 7: Visualize the Relationship

```r
# Create a Mantel scatterplot
# Compare the pairwise distances from both matrices

# Convert distance matrices to vectors
microbiome_vec <- as.vector(microbiome_dist_matched)
rnaseq_vec <- as.vector(rnaseq_dist_matched)

# Create dataframe for plotting
plot_data <- data.frame(
  Microbiome_Distance = microbiome_vec,
  RNAseq_Distance = rnaseq_vec
)

# Create scatterplot
library(ggplot2)
mantel_plot <- ggplot(plot_data, aes(x = Microbiome_Distance, y = RNAseq_Distance)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(
    title = "Mantel Test: Microbiome vs. RNA-seq Distances",
    subtitle = paste0("Mantel r = ", round(mantel_result$statistic, 3), 
                     ", p = ", round(mantel_result$signif, 4)),
    x = "Microbiome Distance (Aitchison)",
    y = "RNA-seq Distance (Euclidean)",
    caption = paste("Based on", length(common_samples), "samples")
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11)
  )

# Save plot
ggsave("../output/mantel_test_plot.png", 
       mantel_plot, 
       width = 8, 
       height = 6, 
       dpi = 300)

print(mantel_plot)
```

## Partial Mantel Test

If you want to control for confounding variables (e.g., developmental stage, treatment), use a **partial Mantel test**:

```r
# Create a distance matrix for the confounding variable
# For example, if you want to control for developmental stage:

# Import metadata
metadata <- read.csv("../metadata/metadata.csv", row.names = 1)

# Match to common samples
metadata_matched <- metadata[common_samples, ]

# Create a distance matrix based on the confounding variable
# For categorical variables (e.g., timepoint):
confound_matrix <- model.matrix(~ timepoint - 1, data = metadata_matched)
confound_dist <- dist(confound_matrix, method = "euclidean")

# Or for continuous variables (e.g., leachate concentration):
confound_dist <- dist(metadata_matched$leachate_mgL, method = "euclidean")

# Run partial Mantel test
partial_mantel_result <- mantel.partial(microbiome_dist_matched,
                                       rnaseq_dist_matched,
                                       confound_dist,
                                       method = "spearman",
                                       permutations = 9999)

print(partial_mantel_result)
```

## Complete Analysis Script Template

```r
#!/usr/bin/env Rscript
# Mantel Test: Microbiome-Transcriptome Correlation Analysis
# Author: Your Name
# Date: YYYY-MM-DD

# ============================================================================
# Setup
# ============================================================================

library(vegan)
library(tidyverse)
library(compositions)
library(DESeq2)
library(ggplot2)

set.seed(123)  # For reproducibility

# ============================================================================
# Import Data
# ============================================================================

# Microbiome data
asv_table <- read_tsv("salipante/241121_StonyCoral/270x200/1st_collapse/feature-table.tsv",
                      skip = 1)
asv_matrix <- asv_table %>%
  column_to_rownames("#OTU ID") %>%
  t() %>%
  as.data.frame()

# RNA-seq data
rnaseq_counts <- read.csv("path/to/rnaseq_counts.csv", row.names = 1)
if (nrow(rnaseq_counts) > ncol(rnaseq_counts)) {
  rnaseq_counts <- t(rnaseq_counts)
}

# Metadata
metadata <- read.csv("../metadata/metadata.csv", row.names = 1)

# ============================================================================
# Data Preparation
# ============================================================================

# Filter rare ASVs
prevalence_threshold <- 0.1 * nrow(asv_matrix)
asv_matrix_filtered <- asv_matrix[, colSums(asv_matrix > 0) >= prevalence_threshold]

# CLR transformation for microbiome
pseudocount <- 0.5
asv_clr <- clr(asv_matrix_filtered + pseudocount)
microbiome_dist <- dist(asv_clr, method = "euclidean")

# VST transformation for RNA-seq
dds <- DESeqDataSetFromMatrix(
  countData = t(rnaseq_counts),
  colData = data.frame(sample = rownames(rnaseq_counts)),
  design = ~ 1
)
vst_data <- vst(dds, blind = TRUE)
rnaseq_transformed <- t(assay(vst_data))
rnaseq_dist <- dist(rnaseq_transformed, method = "euclidean")

# ============================================================================
# Match Samples
# ============================================================================

common_samples <- intersect(labels(microbiome_dist), labels(rnaseq_dist))
microbiome_dist_matched <- as.dist(as.matrix(microbiome_dist)[common_samples, common_samples])
rnaseq_dist_matched <- as.dist(as.matrix(rnaseq_dist)[common_samples, common_samples])

cat("Analyzing", length(common_samples), "matched samples\n")

# ============================================================================
# Mantel Test
# ============================================================================

mantel_result <- mantel(microbiome_dist_matched,
                        rnaseq_dist_matched,
                        method = "spearman",
                        permutations = 9999)

print(mantel_result)

# ============================================================================
# Visualization
# ============================================================================

plot_data <- data.frame(
  Microbiome_Distance = as.vector(microbiome_dist_matched),
  RNAseq_Distance = as.vector(rnaseq_dist_matched)
)

mantel_plot <- ggplot(plot_data, aes(x = Microbiome_Distance, y = RNAseq_Distance)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(
    title = "Mantel Test: Microbiome vs. RNA-seq",
    subtitle = paste0("r = ", round(mantel_result$statistic, 3),
                     ", p = ", round(mantel_result$signif, 4)),
    x = "Microbiome Distance",
    y = "RNA-seq Distance"
  ) +
  theme_bw()

ggsave("../output/mantel_test_plot.png", mantel_plot, width = 8, height = 6, dpi = 300)

# ============================================================================
# Save Results
# ============================================================================

results_summary <- data.frame(
  Test = "Mantel Test",
  N_Samples = length(common_samples),
  Mantel_r = mantel_result$statistic,
  P_value = mantel_result$signif,
  Permutations = 9999,
  Method = "Spearman",
  Microbiome_Metric = "Aitchison (CLR + Euclidean)",
  RNAseq_Metric = "Euclidean (VST)"
)

write.csv(results_summary, "../output/mantel_test_results.csv", row.names = FALSE)

cat("\nAnalysis complete. Results saved to ../output/\n")
```

## Troubleshooting

### Problem: No RNA-seq data available

**Solution:** If you don't have matched RNA-seq data yet:
1. Focus on within-microbiome analyses (PERMANOVA, differential abundance)
2. Plan RNA-seq experiment to match your microbiome samples
3. Ensure sample IDs, timepoints, and treatments match between datasets

### Problem: Different numbers of samples

**Solution:** Only include samples that exist in both datasets. The script above handles this with the `intersect()` function.

### Problem: Zeros in the data

**Solution:** 
- Microbiome: Add a small pseudocount (0.5 or 1) before CLR transformation
- RNA-seq: Use VST or rlog transformation which handle zeros gracefully

### Problem: Non-significant result

**Possible explanations:**
1. True lack of correlation between microbiome and transcriptome
2. Insufficient sample size (Mantel tests require adequate power)
3. Strong confounding variables obscuring the relationship (try partial Mantel test)
4. Inappropriate distance metrics for your data
5. High within-group variability relative to between-group differences

## Alternative Approaches

If the Mantel test is not appropriate for your data, consider:

1. **Procrustes Analysis**: Compares the overall configuration of samples in ordination space
2. **Co-inertia Analysis**: Identifies common patterns between two datasets
3. **Canonical Correspondence Analysis (CCA)**: Tests if one dataset explains variation in another
4. **Multi-omics integration methods**: MOFA, DIABLO, mixOmics

## References

1. Mantel, N. (1967). The detection of disease clustering and a generalized regression approach. *Cancer Research*, 27(2), 209-220.

2. Legendre, P., & Legendre, L. (2012). *Numerical Ecology* (3rd ed.). Elsevier.

3. Oksanen, J., et al. (2020). vegan: Community Ecology Package. R package version 2.5-7. https://CRAN.R-project.org/package=vegan

4. Gloor, G. B., et al. (2017). Microbiome Datasets Are Compositional: And This Is Not Optional. *Frontiers in Microbiology*, 8, 2224.

5. Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12), 550.

## Additional Resources

- **vegan package documentation**: https://rdrr.io/rforge/vegan/man/mantel.html
- **Mantel test tutorial**: https://archetypalecology.wordpress.com/2018/02/21/distance-matrices-and-the-mantel-test/
- **Compositional data analysis**: Gloor et al. 2017 paper in repository (`Gloor et al. 2017 - Front. Microbiol_.pdf`)
- **QIIME2 to R workflow**: See `compositional_analysis.qmd` in this repository

## Citation

If you use this analysis in your research, please cite:

```
Oksanen, J., Simpson, G.L., Blanchet, F.G., Kindt, R., Legendre, P., Minchin, P.R., 
O'Hara, R.B., Solymos, P., Stevens, M.H.H., Szoecs, E., Wagner, H., Barbour, M., 
Bedward, M., Bolker, B., Borcard, D., Carvalho, G., Chirico, M., De Caceres, M., 
Durand, S., Evangelista, H.B.A., FitzJohn, R., Friendly, M., Furneaux, B., 
Hannigan, G., Hill, M.O., Lahti, L., McGlinn, D., Ouellette, M., Ribeiro Cunha, E., 
Smith, T., Stier, A., Ter Braak, C.J.F., Weedon, J. (2022). vegan: Community 
Ecology Package. R package version 2.6-4, https://CRAN.R-project.org/package=vegan.
```
