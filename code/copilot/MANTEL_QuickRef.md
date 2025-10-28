# Mantel Test Quick Reference

## What is it?
Tests correlation between two distance matrices using permutation.

**Use case:** Assess if patterns in microbiome composition correlate with patterns in host gene expression.

## When to use
- ✅ Comparing two multivariate datasets
- ✅ Testing microbiome-transcriptome relationships
- ✅ Non-parametric correlation needed
- ✅ Compositional and count data (after transformation)

## Basic workflow

```r
library(vegan)

# 1. Create distance matrices
microbiome_dist <- dist(microbiome_transformed)
rnaseq_dist <- dist(rnaseq_transformed)

# 2. Run Mantel test
result <- mantel(microbiome_dist, rnaseq_dist, 
                 method = "spearman", 
                 permutations = 9999)

# 3. Check results
print(result)
```

## Key parameters

| Parameter | Options | Notes |
|-----------|---------|-------|
| `method` | `"pearson"`, `"spearman"`, `"kendall"` | Spearman recommended for most cases |
| `permutations` | Integer (typically 999 or 9999) | Higher = more accurate p-value |

## Interpreting results

**Mantel statistic (r):**
- Range: -1 to +1
- `r > 0`: Positive correlation (similar samples in both datasets)
- `r < 0`: Negative correlation (similar samples in one, dissimilar in other)
- `|r| > 0.4`: Moderate to strong correlation

**P-value:**
- `p < 0.05`: Significant correlation
- `p ≥ 0.05`: No significant correlation

## Data preparation tips

### Microbiome (16S/ASV data)
```r
# CLR transformation for compositional data
library(compositions)
asv_clr <- clr(asv_counts + 0.5)
microbiome_dist <- dist(asv_clr, method = "euclidean")  # Aitchison distance
```

### RNA-seq data
```r
# VST transformation for count data
library(DESeq2)
dds <- DESeqDataSetFromMatrix(counts, colData, design = ~1)
vst_data <- vst(dds, blind = TRUE)
rnaseq_dist <- dist(t(assay(vst_data)), method = "euclidean")
```

## Common issues

| Problem | Solution |
|---------|----------|
| Different sample sizes | Use `intersect()` to find common samples |
| Sample order mismatch | Subset both matrices to same samples in same order |
| Non-significant result | Try partial Mantel (control confounders) or subset analysis |
| Zeros in data | Add pseudocount (microbiome) or use VST (RNA-seq) |

## Controlling for confounders

Use **partial Mantel test** to control for a third variable:

```r
# Create distance matrix for confounding variable
confound_dist <- dist(metadata$timepoint)

# Run partial Mantel test
partial_result <- mantel.partial(microbiome_dist, 
                                 rnaseq_dist, 
                                 confound_dist,
                                 method = "spearman",
                                 permutations = 9999)
```

## Quick checklist

- [ ] Both datasets from same samples?
- [ ] Sample IDs match exactly?
- [ ] Data appropriately transformed?
- [ ] Rare features filtered?
- [ ] Distance metrics appropriate for data type?
- [ ] Sufficient sample size (n ≥ 15 recommended)?
- [ ] Checked for confounding variables?

## Related analyses

- **Procrustes analysis**: Compare ordination configurations
- **Co-inertia analysis**: Find common patterns between datasets
- **CCA/RDA**: Test if one dataset explains variation in another
- **MOFA/DIABLO**: Multi-omics integration methods

## References

- Mantel, N. (1967). *Cancer Research*, 27(2), 209-220.
- Legendre & Legendre (2012). *Numerical Ecology* (3rd ed.)
- `vegan::mantel()` documentation: https://rdrr.io/rforge/vegan/man/mantel.html

## Full documentation

See **[MANTEL_TEST_GUIDE.md](MANTEL_TEST_GUIDE.md)** for comprehensive tutorial with examples.

See **[mantel_test_analysis.qmd](mantel_test_analysis.qmd)** for executable workflow.
