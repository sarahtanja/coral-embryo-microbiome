# MaAsLin3 Analysis Documentation

## Quick Summary

The MaAsLin3 analysis in this repository can use various model specifications including **additive models** that test main effects of leachate and timepoint separately, or **interaction models** that test whether effects depend on each other. MaAsLin3 improves on MaAsLin2 by better accounting for compositionality.

### Key Files

- **[maaslin2_interpretation.md](maaslin2_interpretation.md)** - Comprehensive guide to understanding the model design and interpreting results (historical MaAsLin2 documentation)
- **[maaslin3.qmd](maaslin3.qmd)** - Current analysis code using MaAsLin3 with multiple model specifications
- **[explore.qmd](explore.qmd)** - Exploration of significant results
- **`../salipante/Sarah_StonyCoral/Level7_filtered_organism_output/significant_results.tsv`** - Results file

## The Issue: Loss of Granularity

### What the Current Model Tests

The model `fixed_effects = c("leachate", "timepoint")` creates an additive formula:

```
Abundance ~ leachate + timepoint + (1|cross)
```

This tests:
- **Leachate effect**: averaged across all timepoints (4, 9, 14 hpf)
- **Timepoint effect**: averaged across all leachate levels (0, 0.01, 0.1, 1)

### What We CANNOT Determine

From the current results alone:
- ❌ Which specific leachate concentration(s) drive the difference
- ❌ Which specific timepoint(s) drive the difference  
- ❌ Whether leachate effects differ at different timepoints
- ❌ Whether timepoint effects differ at different leachate levels

### What the Results Tell Us

From the 732 significant results:
- **607 bacteria**: Significant timepoint effects (developmental changes consistent across leachate treatments)
- **125 bacteria**: Significant leachate effects (responses consistent across developmental stages)

## Solutions

### Option 1: Use Interaction Model (MaAsLin3)

Test whether effects depend on each other using MaAsLin3:

```r
formula = ~ leachate * hpf + reads + (1|spawn_night)
```

This identifies bacteria where:
- Leachate effects differ by developmental stage
- Developmental changes differ by leachate treatment
- See [maaslin3.qmd](maaslin3.qmd) for detailed examples

### Option 2: Categorical Analysis (MaAsLin3)

Get specific pairwise comparisons using ordered factors:

```r
# MaAsLin3 can use ordered factors for categorical analysis
formula = ~ ordered(leachate) + stage + reads + (1|spawn_night)
```

This provides level contrasts and comparisons:
- Tests for linear, quadratic, and cubic trends
- See [maaslin3.qmd](maaslin3.qmd) for level contrast models

### Option 3: Post-hoc Visualization

For specific bacteria of interest:
- Plot abundance profiles across conditions
- Conduct targeted pairwise comparisons
- Examine raw data patterns

## Interpreting Current Results

### Example: Timepoint Effect

```
feature: d__Bacteria.p__Proteobacteria...g__Algicola.s__uncultured_bacterium
metadata: timepoint
coef: 1.66
```

**Means:** For each 1-hour increase in developmental time, log-abundance increases by 1.66 (averaged across all leachate levels)

### Example: Leachate Effect

```
feature: d__Bacteria.p__Proteobacteria...g__Cognatishimia.s__uncultured_bacterium
metadata: leachate
coef: -0.59
```

**Means:** For each 1-unit increase in leachate concentration, log-abundance decreases by 0.59 (averaged across all timepoints)

## Experimental Design Context

### Metadata Structure

- **Leachate concentrations**: 0, 0.01, 0.1, 1 (treated as numeric)
- **Timepoints**: 4, 9, 14 hours post fertilization (treated as numeric)
- **Cross**: Parental crosses (random effect)
- **Development stages**: cleavage (4 hpf), prawn chip (9 hpf), early gastrula (14 hpf)

### Sample Size

- 63 total samples
- 4 leachate levels × 3 timepoints = 12 treatment combinations
- Multiple crosses per combination

## Model Comparison

| Approach | Answers | Advantages | Limitations |
|----------|---------|------------|-------------|
| **Current (Additive)** | Are there overall trends? | Simple, clear main effects | Loses specificity, assumes no interaction |
| **With Interaction** | Do effects depend on each other? | Detects condition-specific responses | Complex interpretation, needs more power |
| **Categorical** | Which levels differ from reference? | Specific comparisons | Many tests, doesn't assess trends |

## Recommendations

1. **Current results are valid** - they show bacteria with consistent trends across conditions
2. **For more specificity**: Choose Option 1 (interaction) or Option 2 (categorical) depending on your research question
3. **For individual bacteria**: Use post-hoc visualization and targeted analysis

## Further Reading

- Full interpretation guide: [maaslin2_interpretation.md](maaslin2_interpretation.md) (historical MaAsLin2 documentation)
- MaAsLin3 analysis: [maaslin3.qmd](maaslin3.qmd) (current implementation)
- MaAsLin3 documentation: https://huttenhower.sph.harvard.edu/MaAsLin3
- MaAsLin3 tutorial: https://github.com/biobakery/biobakery/wiki/maaslin3
- Original MaAsLin2 paper: Mallick et al. (2021) PLOS Computational Biology https://doi.org/10.1371/journal.pcbi.1009442
