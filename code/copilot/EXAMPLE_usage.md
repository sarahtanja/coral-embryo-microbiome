# Example Usage: Taxa Barplot Visualization

## Quick Start

To generate taxa barplots from your QIIME2 visualization file, follow these steps:

### 1. Install Required Packages (First Time Only)

Open R or RStudio and run:

```r
# Install packages if not already installed
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("scales")) install.packages("scales")
if (!require("viridis")) install.packages("viridis")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
```

### 2. Render the Document

#### From Command Line (with Quarto installed):

```bash
cd /path/to/coral-embryo-microbiome
quarto render code/taxa_barplot_visualization.qmd
```

#### From RStudio:

1. Open `code/taxa_barplot_visualization.qmd`
2. Click the "Render" button in the toolbar
3. Wait for processing to complete

### 3. View Results

After rendering, you'll find:

- **HTML Report**: `code/taxa_barplot_visualization.html` - Interactive report with all visualizations
- **Figures**: `output/figures/` directory containing:
  - `taxa_barplot_all_samples.png`
  - `taxa_barplot_by_treatment.png`
  - `taxa_barplot_phylum.png`
  - `diversity_by_treatment.png`

## Understanding the Output

### Figure 1: All Samples Barplot
Shows the taxonomic composition (top 15 genera) for each individual sample. Useful for:
- Identifying sample-to-sample variation
- Spotting outliers or contaminated samples
- Visualizing overall community composition

### Figure 2: Treatment-Faceted Barplot
Shows mean relative abundance by treatment group (leachate × timepoint). Useful for:
- Comparing communities across experimental conditions
- Identifying treatment effects
- Manuscript main figures

### Figure 3: Phylum-Level Barplot
Shows composition at phylum level (higher taxonomic rank). Useful for:
- Broad overview of community composition
- Initial screening of major community shifts
- Supplementary figures

### Figure 4: Diversity by Treatment
Line plot showing taxonomic richness across conditions. Useful for:
- Tracking diversity changes over time
- Comparing diversity between treatments
- Statistical analysis preparation

## Common Modifications

### Change Number of Taxa Shown

In the "Select top taxa for visualization" section:

```r
# Show top 20 instead of 15
top_genera <- taxa_long %>%
  group_by(genus) %>%
  summarize(mean_abundance = mean(rel_abundance)) %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = 20)  # Change this number
```

### Use Different Taxonomic Level

To use family level instead of genus:

```r
# In "Transform data for visualization" section, change:
taxa_long <- taxa_l5 %>%  # Changed from taxa_l6 to taxa_l5
  pivot_longer(...)

# And extract family instead of genus:
mutate(
  family = str_extract(taxon, "f__([^;_]+)", group = 1)
)
```

### Adjust Figure Size

In chunk options:

```r
#| fig-width: 16    # Make wider
#| fig-height: 10   # Make taller
#| fig-dpi: 600     # Increase resolution
```

### Change Color Scheme

```r
# Use different ColorBrewer palette
colors <- brewer.pal(n_colors, "Set1")  # Instead of "Set3"

# Or use custom colors
colors <- c("#FF6B6B", "#4ECDC4", "#45B7D1", "#FFA07A", ...)
```

## Troubleshooting

### "Cannot find file" Error

Make sure you're running from the repository root:

```bash
# Check current directory
pwd

# Should end with: coral-embryo-microbiome
# If not, navigate there first
cd /path/to/coral-embryo-microbiome
```

### Plots Look Cluttered

If you have many samples, try:

1. Increase figure width: `fig-width: 18`
2. Reduce text size: `theme(axis.text.x = element_text(size = 6))`
3. Create separate plots for different timepoints

### Memory Issues

If R runs out of memory:

```r
# At the start of your session
memory.limit(size = 8000)  # Windows only

# Or load only one taxonomic level at a time
```

## Integration with Other Analyses

The processed data can be exported for other analyses:

```r
# Save processed data
write_csv(taxa_long, "../output/taxa_long_genus.csv")
write_csv(top_taxa_summary, "../output/top_taxa_summary.csv")

# Save for STAMP or other tools
taxa_wide <- taxa_long %>%
  select(sample_id, genus, rel_abundance) %>%
  pivot_wider(names_from = genus, values_from = rel_abundance)
write_csv(taxa_wide, "../output/taxa_wide_for_stamp.csv")
```

## Next Steps

After generating your taxa barplots:

1. **Statistical Analysis**: Use maaslin2 or ANCOM for differential abundance testing
2. **Diversity Metrics**: Calculate alpha and beta diversity using phyloseq or qiime2
3. **Functional Prediction**: Use PICRUSt2 for functional profiling
4. **Network Analysis**: Explore co-occurrence patterns

## Tips for Publication

1. **High Resolution**: Use `fig-dpi: 600` for publication quality
2. **File Format**: Save as both PNG (for drafts) and PDF (for final submission)
3. **Color Blind Friendly**: Use `viridis` or ColorBrewer palettes
4. **Figure Legends**: Add detailed legends describing samples, treatments, and methods
5. **Consistency**: Use the same color scheme across all figures in your manuscript

## Questions?

If you encounter issues or need help customizing the visualizations, please:
- Check the README_taxa_barplot.md for detailed documentation
- Review the comments in the .qmd file
- Open an issue on GitHub with your specific question
