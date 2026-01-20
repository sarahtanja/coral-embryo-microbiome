# Salipante Lab Analysis Directory

This directory contains microbiome analysis results from the Salipante Lab (mim_c) for the Stony Coral embryo microbiome project.

## Directory Structure

### Sarah_StonyCoral/

Primary analysis directory containing filtered QIIME2 results and historical MaAsLin2 differential abundance analysis outputs. For current MaAsLin3 analysis, see [../code/differentially-abundant-taxa/maaslin3.qmd](../code/differentially-abundant-taxa/maaslin3.qmd).

**Contents:**
- `250405_QIIME2_filtered_StonyCloral_final.txt` - Final QIIME2 analysis commands for filtered data
- `250405_RScript_Maaslin2_StonyCloral_final.txt` - R script for historical MaAsLin2 differential abundance analysis
- `241121_StonyCoral_readcounts.xlsx` - Sequencing read count summaries
- `250425_metadata_only_maaslin2.txt` - Metadata file formatted for MaAsLin2 input
- `250425_shannon_alpha_diversity.xlsx` - Shannon diversity metrics
- `250414_StonyCoral_imported_reads_trimmed.qzv` - QIIME2 visualization of imported and trimmed reads
- `250414_StonyCoral_270x200_filtered_taxa-bar-plots.qzv` - Taxa bar plots for filtered data
- `250416_StonyCoral_270x200_alpha_rarefaction_curves_400k.qzv` - Alpha rarefaction curves
- `filtered_shannon_vector.qzv` - Shannon diversity vector visualization

**Subdirectory:**
- `Level7_filtered_organism_output/` - Historical MaAsLin2 results for species-level (Level 7) filtered data
  - `significant_results.tsv` - Table of statistically significant differential abundance results
  - `all_results.tsv` - Complete results from differential abundance analysis
  - `leachate.pdf` & `timepoint.pdf` - Summary plots for treatment effects
  - `heatmap.pdf` - Heatmap of significant features
  - `figures/` - Individual feature plots (leachate_1-10.png, timepoint_1-10.png)
  - `fitted.rds`, `ranef.rds`, `residuals.rds` - Model objects
  - `maaslin2.log` - Analysis log file

### 241121_StonyCoral/

QIIME2 pipeline outputs and diversity analyses for the StonyCoral sequencing run from November 21, 2024.

**Contents:**
- `250414_StonyCoral_read_manifest.tsv` - Manifest file linking sample IDs to raw sequencing files
- `250416_StonyCoral_metadata.tsv` - Sample metadata (categorical and continuous variables)
- `250416_StonyCoral_metadata_categoric.tsv` - Metadata with categorical variables only

**Subdirectories:**

- `270x200/` - DADA2 denoising parameters (270 forward, 200 reverse trim lengths)
  - Feature tables (filtered and unfiltered) `.qza` and `.qzv` files
  - Representative sequences and taxonomy assignments
  - Denoising statistics
  - Taxa bar plots (filtered and unfiltered)
  - Alpha rarefaction curves (200K and 400K depth)
  - Collapsed taxonomy tables at levels 1-7 (`collapsed-l1.qza` through `collapsed-l7.qza`)
  - `taxa_collapse_filtered.sh` - Shell script for taxonomy collapsing
  - `filtered_hierarchical.tsv` - Hierarchical taxonomy table
  - `1st_collapse/` - Initial taxonomy collapse outputs
  - `hierarchical-export/` - Exported hierarchical data

- `SEPP/` - Phylogenetic placement results using SEPP algorithm
  - `insertion-placements.qza` - Placement positions on reference tree
  - `insertion-tree.qza` - Phylogenetic tree with inserted sequences

- `core-metrics-results/` - Core diversity metrics (unfiltered data)
  - Alpha diversity vectors (Shannon, Faith's PD, Evenness, Observed Features)
  - Beta diversity distance matrices (Bray-Curtis, Jaccard, UniFrac weighted/unweighted)
  - PCoA results and Emperor visualizations
  - Statistical significance tests for metadata categories

- `core-metrics-results-filtered/` - Core diversity metrics (filtered data, excluding chloroplast/mitochondria)
  - Same structure as core-metrics-results but with filtered dataset

### 250520_picrust_Maaslin2/

PICRUSt2 functional predictions analyzed with Maaslin2 for differential abundance (dated May 20, 2025).

**Subdirectories:**

- `filtered_ec_TSS_output/` - Enzyme Commission (EC) number predictions, TSS normalized
  - `significant_results_TSS.tsv` - Significant differential abundant enzyme functions
  - `all_results.tsv` - Complete EC prediction results
  - Maaslin2 outputs (PDFs, model objects, figures)

- `filtered_ko_TSS_output/` - KEGG Orthology (KO) predictions, TSS normalized
  - `ko_significant_results_TSS.tsv` - Significant differential abundant KEGG functions
  - `all_results.tsv` - Complete KO prediction results
  - Maaslin2 outputs (PDFs, model objects, figures)

- `filtered_path_TSS_output/` - MetaCyc pathway predictions, TSS normalized
  - `significant_results.tsv` - Significant differential abundant metabolic pathways
  - `all_results.tsv` - Complete pathway prediction results
  - Maaslin2 outputs (PDFs, model objects, figures)

## Analysis Notes

- Data was processed using QIIME2 with DADA2 denoising (270x200 trim parameters)
- Chloroplast and mitochondrial sequences were filtered from the dataset
- Differential abundance was assessed using Maaslin2 with fixed effects for leachate treatment and timepoint, and random effects for cross
- Functional predictions were generated using PICRUSt2 and analyzed with TSS normalization
- 732 significantly different bacteria were identified across leachate levels and developmental timepoints
- After filtering, Shannon diversity increased from 4.6 to 5.7 on average
