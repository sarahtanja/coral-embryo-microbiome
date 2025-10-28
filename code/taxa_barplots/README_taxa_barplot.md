# Taxa Barplot Visualization Guide

## Overview

The `taxa_barplot_visualization.qmd` file creates publication-ready taxa barplot visualizations from QIIME2 visualization files (.qzv). This script extracts taxonomic composition data and visualizes it using ggplot2 in R.

## Prerequisites

### Required R Packages

Install the following R packages before running the script:

```r
install.packages(c("tidyverse", "scales", "viridis", "RColorBrewer"))
```

### Required Files

1. **QIIME2 visualization file**: `salipante/Sarah_StonyCoral/250414_StonyCoral_270x200_filtered_taxa-bar-plots.qzv`
2. **Metadata file**: `metadata/metadata.tsv`

## How to Run

### Option 1: Render with Quarto (Recommended)

From the command line in the repository root:

```bash
quarto render code/taxa_barplot_visualization.qmd
```

This will generate both HTML and GitHub Markdown outputs.

### Option 2: Run in RStudio

1. Open `code/taxa_barplot_visualization.qmd` in RStudio
2. Click the "Render" button
3. The output will be generated in the same directory

### Option 3: Run Individual Code Chunks

Open the file in RStudio and run code chunks interactively using Ctrl+Enter (Windows/Linux) or Cmd+Enter (Mac).

## Output Files

The script generates the following outputs:

### Figures (saved to `output/figures/`)

1. **taxa_barplot_all_samples.png** - Bar plot showing taxonomic composition for all individual samples
2. **taxa_barplot_by_treatment.png** - Faceted bar plot showing mean composition by leachate treatment and timepoint
3. **taxa_barplot_phylum.png** - Bar plot at phylum level
4. **diversity_by_treatment.png** - Line plot showing taxonomic richness across treatments and timepoints

### Reports

- HTML report with interactive elements
- GitHub Markdown (.md) file for viewing on GitHub

## Customization

### Changing Taxonomic Level

The script is set to use Level 6 (genus level) by default. To change this:

1. Locate the "Transform data for visualization" section
2. Change `taxa_l6` to `taxa_l5` (family) or `taxa_l7` (species)

### Adjusting Number of Top Taxa

To show more or fewer taxa in the plots:

1. Find the line: `slice_head(n = 15)`
2. Change `15` to your desired number

### Color Schemes

The script uses:
- `brewer.pal("Set3")` for genus-level plots
- `brewer.pal("Paired")` for phylum-level plots
- `viridis` for diversity plots

These can be changed to other ColorBrewer palettes or custom color vectors.

### Plot Dimensions

Adjust figure dimensions in the chunk options:
```r
#| fig-width: 12
#| fig-height: 8
```

### Plot DPI

For higher quality images, increase the DPI:
```r
#| fig-dpi: 300  # Change to 600 for publication quality
```

## Troubleshooting

### Error: "Cannot open connection"

If you get file path errors:
- Check that the .qzv file exists at the specified location
- Verify the metadata file path
- Ensure you're running from the repository root or adjust paths accordingly

### Error: Package not installed

Install missing packages:
```r
install.packages("package_name")
```

### Memory Issues with Large Datasets

If you encounter memory issues:
- Reduce the number of taxonomic levels loaded
- Process one level at a time
- Increase R memory limit: `memory.limit(size = 8000)` (Windows)

## Data Structure

### Input Data Format

The .qzv file contains multiple CSV files:
- `level-1.csv` through `level-7.csv` representing different taxonomic ranks
- Each CSV has samples as rows and taxa as columns
- Taxonomic strings follow QIIME2 format: `d__Domain;p__Phylum;c__Class;...`

### Metadata Format

The metadata file should be tab-separated with:
- First line: column names
- Second line: `#q2:types` declarations
- Subsequent lines: sample data

Required columns:
- `sample-id`: Sample identifier
- `development-stage`: Developmental stage
- `timepoint`: Time in hours post fertilization
- `leachate`: Leachate concentration

## Citation

If you use this script, please cite:
- QIIME2 (Bolyen et al., 2019)
- ggplot2 (Wickham, 2016)
- tidyverse packages

## Questions?

For issues or questions about this script, please open an issue on the GitHub repository or contact the author.
