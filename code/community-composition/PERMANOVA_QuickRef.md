# PERMANOVA Quick Reference Card

## What is PERMANOVA?

**PERMANOVA** = Permutational Multivariate Analysis of Variance

- Non-parametric test for multivariate data
- Tests whether groups differ in community composition
- Uses permutations instead of parametric assumptions
- Works with distance/dissimilarity matrices

## When to Use PERMANOVA in Microbiome Analysis

✅ **Use PERMANOVA when you want to:**
- Test if microbial communities differ between treatment groups
- Assess effects of experimental factors on beta diversity
- Determine if time/development affects community composition
- Statistically validate patterns seen in PCoA plots

❌ **Don't use PERMANOVA for:**
- Testing individual taxa abundance (use DESeq2, ANCOMBC, or MaAsLin3)
- Alpha diversity comparisons (use ANOVA, Kruskal-Wallis)
- Testing ordination coordinates directly (use distance matrix instead!)

## Quick Workflow

```r
library(vegan)

# 1. Calculate distance matrix
dist_matrix <- vegdist(feature_table, method = "bray")

# 2. Check dispersions (required!)
disp_test <- betadisper(dist_matrix, metadata$treatment)
permutest(disp_test)  # Should be non-significant

# 3. Run PERMANOVA
result <- adonis2(
  dist_matrix ~ treatment * timepoint,
  data = metadata,
  permutations = 999
)

# 4. Visualize with PCoA (for visualization only!)
pcoa <- cmdscale(dist_matrix, k = 2)
```

## Understanding the Output

```
                  Df  SumOfSqs    R2      F    Pr(>F)    
treatment          3   2.584   0.242  4.862   0.001 ***
timepoint          2   1.762   0.165  4.975   0.001 ***
treatment:time     6   0.923   0.086  0.869   0.663    
Residuals         48   8.530   0.517                  
Total             59  16.800   1.000
```

### Key Metrics

| Metric | Meaning | Interpretation |
|--------|---------|----------------|
| **F** | Pseudo-F statistic | Higher = stronger group separation |
| **R²** | Effect size | Proportion of variation explained (0-1) |
| **Pr(>F)** | p-value | p < 0.05 = significant difference |
| **Df** | Degrees of freedom | Based on number of groups |

### Interpreting R²

- **R² < 0.10**: Small effect (< 10% variation)
- **R² = 0.10-0.30**: Moderate effect
- **R² > 0.30**: Large effect

**Example**: R² = 0.24 means treatment explains 24% of variation in community composition.

## Distance Metrics for Microbiome Data

| Metric | Type | When to Use |
|--------|------|-------------|
| **Aitchison** | Compositional | **Recommended** for count data (CLR-transformed) |
| **Bray-Curtis** | Abundance | Common, intuitive, ignores double-zeros |
| **Jaccard** | Presence/Absence | Focus on rare taxa |
| **Weighted UniFrac** | Phylogenetic + Abundance | When phylogeny matters |
| **Unweighted UniFrac** | Phylogenetic + P/A | Phylogeny, focus on rare taxa |

## Common Mistakes ⚠️

### ❌ WRONG: Running PERMANOVA on PCoA coordinates

```r
# This is INCORRECT!
pcoa_coords <- cmdscale(dist_matrix, k = 2)
wrong_result <- adonis2(pcoa_coords ~ treatment, data = metadata)
```

### ✅ CORRECT: Running PERMANOVA on distance matrix

```r
# This is CORRECT!
correct_result <- adonis2(dist_matrix ~ treatment, data = metadata)
```

### ❌ WRONG: Ignoring dispersion violations

```r
# Dispersions differ significantly - PERMANOVA may be unreliable!
# But still interpreting results...
```

### ✅ CORRECT: Check and report dispersions

```r
disp_test <- betadisper(dist_matrix, metadata$treatment)
disp_p <- permutest(disp_test)$tab$`Pr(>F)`[1]

if (disp_p < 0.05) {
  message("Warning: Groups differ in dispersion. PERMANOVA results may be unreliable.")
}
```

## Reporting Results

### In Methods

> "We performed PERMANOVA using the adonis2 function in vegan (v2.6) in R. We calculated Aitchison distance from CLR-transformed data and tested effects of treatment and timepoint with 999 permutations. We verified homogeneity of dispersions using PERMDISP."

### In Results

> "PERMANOVA revealed a significant effect of treatment on community composition (F₃,₄₈ = 4.86, R² = 0.24, p = 0.001), explaining 24% of variation. The treatment × timepoint interaction was not significant (p = 0.663)."

## Checklist Before Running PERMANOVA

- [ ] Data is appropriately transformed (CLR for compositional data)
- [ ] Distance metric is appropriate for your data type
- [ ] Sample sizes are adequate (≥5 per group recommended)
- [ ] Samples are independent (no pseudoreplication)
- [ ] You will check dispersion assumptions
- [ ] You will report both R² and p-values
- [ ] You will perform pairwise tests if >2 groups

## Need More Detail?

**Full Documentation:**
- [`PERMANOVA_README.md`](PERMANOVA_README.md) - Complete theoretical guide
- [`permanova_analysis.qmd`](permanova_analysis.qmd) - Working code examples

## Key Citations

- Anderson (2001) - Original PERMANOVA paper
- Anderson & Walsh (2013) - Dispersion assumptions
- Gloor et al. (2017) - Compositional data for microbiomes

---

**Remember**: PERMANOVA tells you IF groups differ in composition. To identify WHICH taxa drive differences, use differential abundance methods (ANCOMBC, DESeq2, MaAsLin3).
