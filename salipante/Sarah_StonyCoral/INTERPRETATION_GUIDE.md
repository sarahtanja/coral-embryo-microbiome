# How to Interpret Maaslin2 Results

This directory contains Maaslin2 differential abundance analysis results. **Important:** The model design affects how we interpret these results.

## Quick Answer to the Issue

**Q: Why do the results only show 'leachate' or 'timepoint' without telling us which specific leachate concentration or timepoint drives the differences?**

**A:** The model uses an **additive design** without interaction terms. This means:
- Leachate effects are averaged across all timepoints
- Timepoint effects are averaged across all leachate concentrations
- We test for overall trends, not specific condition comparisons

## Understanding Your Results

### What's in `Level7_filtered_organism_output/significant_results.tsv`

- **607 bacteria** show significant **timepoint effects**: Their abundance changes with development (4, 9, or 14 hpf), averaged across all leachate treatments
- **125 bacteria** show significant **leachate effects**: Their abundance changes with leachate concentration (0, 0.01, 0.1, or 1), averaged across all developmental stages

### What the coefficients mean

**Example timepoint result:**
```
coef = 1.66 for timepoint
```
→ For each 1-hour increase in developmental time, log-abundance increases by 1.66 (averaged across all leachate levels)

**Example leachate result:**
```
coef = -0.59 for leachate
```
→ For each 1-unit increase in leachate concentration, log-abundance decreases by 0.59 (averaged across all timepoints)

## The Current Model

```r
fixed_effects = c("leachate", "timepoint")
```

Creates the formula: `Abundance ~ leachate + timepoint + (1|cross)`

**Tests:**
- Main effect of leachate ✓
- Main effect of timepoint ✓
- Interaction between leachate and timepoint ✗

## How to Get More Specific Information

### Option 1: Add Interaction Term

Test whether leachate effects differ by timepoint:

```r
fixed_effects = c("leachate", "timepoint", "leachate:timepoint")
```

This would identify bacteria where the effect of leachate depends on developmental stage.

### Option 2: Categorical Analysis

Get specific pairwise comparisons (e.g., 0.01 vs 0, 0.1 vs 0, 1 vs 0):

```r
# Convert to factors first
df_Level7_filtered_metadata$leachate <- as.factor(df_Level7_filtered_metadata$leachate)
df_Level7_filtered_metadata$timepoint <- as.factor(df_Level7_filtered_metadata$timepoint)

# Then run Maaslin2 with reference levels
fit = Maaslin2(..., reference = c("leachate,0", "timepoint,4"))
```

### Option 3: Post-hoc Visualization

For bacteria of interest:
- Plot abundance profiles for each leachate level across timepoints
- Conduct targeted pairwise statistical tests
- Examine patterns in raw data

## Detailed Documentation

For comprehensive explanation of the model, interpretation, and recommendations:

📖 **[../../code/maaslin2_interpretation.md](../../code/maaslin2_interpretation.md)** - Full interpretation guide  
📖 **[../../code/MAASLIN2_README.md](../../code/MAASLIN2_README.md)** - Quick reference  
📖 **[../../code/maaslin2.qmd](../../code/maaslin2.qmd)** - Updated analysis code with alternative models  

## Key Takeaway

Your current results **are valid** and show which bacteria have:
- Overall developmental trends (consistent across leachate treatments)
- Overall leachate responses (consistent across developmental stages)

To answer questions about **specific conditions** or **whether effects depend on each other**, you'll need to run one of the alternative analyses described above.

---

*Generated: 2025-10-15*
