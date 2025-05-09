# Running PICRUSt2
Sarah Tanja

- [Goals](#goals)
- [Load packages](#load-packages)
- [Activate qiime2 environment in R](#activate-qiime2-environment-in-r)
- [Install picrust2](#install-picrust2)
- [Download & install q2-picrust2
  plugin](#download--install-q2-picrust2-plugin)
- [Feature table (BIOM)](#feature-table-biom)
- [Sequence file](#sequence-file)
- [PICRUSt2 qiime pipeline](#picrust2-qiime-pipeline)
- [](#section)

# Goals

Run

> Yes it is reasonable to look at picrust output and see if there are
> functional shifts, but it is not clear how adequate picrust libraries
> are to allow metabolic inference of these particular populations.  I
> would recommend that the actual analysis of metabolic shifts be
> conducted in Maaslin, as before, to account for covariates.
>
> it is possible to give sarah the input files she would need for the
> analysis. \[see
> [https://forum.qiime2.org/t/picrust2-input-files-for-bacteria/27968](https://urldefense.com/v3/__https://forum.qiime2.org/t/picrust2-input-files-for-bacteria/27968__;!!K-Hz7m0Vt54!iZIXAiB2i2QtliEgrFMztEKZHQZuAILNAY54WW_klnDOn0rCNia6tnE-GrFq0efQ8IeJkkAwEJWtIM0$)\].
> For the record, you could isolate the ASV table from DADA2 in qiime
> (these are stored in rep-seqs-dada2.qza, according to
> [https://forum.qiime2.org/t/how-to-generate-an-asv-table-from-dada2-output/9973](https://urldefense.com/v3/__https://forum.qiime2.org/t/how-to-generate-an-asv-table-from-dada2-output/9973__;!!K-Hz7m0Vt54!iZIXAiB2i2QtliEgrFMztEKZHQZuAILNAY54WW_klnDOn0rCNia6tnE-GrFq0efQ8IeJkkAwpWyxhvY$)),
> and convert it to biom format
> \[[https://forum.qiime2.org/t/how-to-generate-an-asv-table-from-dada2-output/9973/2](https://urldefense.com/v3/__https://forum.qiime2.org/t/how-to-generate-an-asv-table-from-dada2-output/9973/2__;!!K-Hz7m0Vt54!iZIXAiB2i2QtliEgrFMztEKZHQZuAILNAY54WW_klnDOn0rCNia6tnE-GrFq0efQ8IeJkkAwc-L0gtk$)\]. 
>
> having said this it may just be easier to use the picrust2 plugin for
> qiime and hand over those outputs, since we already have the inputs
> ready in qiime
>
> [https://github.com/picrust/picrust2/wiki/q2-picrust2-Tutorial](https://urldefense.com/v3/__https://github.com/picrust/picrust2/wiki/q2-picrust2-Tutorial__;!!K-Hz7m0Vt54!iZIXAiB2i2QtliEgrFMztEKZHQZuAILNAY54WW_klnDOn0rCNia6tnE-GrFq0efQ8IeJkkAwVcNVfvQ$)
>
> [https://github.com/picrust/picrust2/wiki/Manually-install-QIIME-2-plugin](https://urldefense.com/v3/__https://github.com/picrust/picrust2/wiki/Manually-install-QIIME-2-plugin__;!!K-Hz7m0Vt54!iZIXAiB2i2QtliEgrFMztEKZHQZuAILNAY54WW_klnDOn0rCNia6tnE-GrFq0efQ8IeJkkAwAif4ef0$)  
>   
> id say generate those outputs and then run them through maaslin to
> look for significant feature abundances. – S.S.

# Load packages

``` r
library(knitr)
library(dplyr)
library(reticulate)
```

# Activate qiime2 environment in R

``` r
use_condaenv(condaenv = "qiime2-2023.5", conda = "/home/shared/8TB_HDD_02/stanja/miniconda3/condabin/conda")

# Check successful env loading
py_config()
```

If this is successful, the first line of output should show that the
Python environment being used is the one in your conda environment path.

Your conda and qiime environments should be listed in your \$PATH:

``` bash
echo $PATH
```

And the following qiime command should run:

``` bash
qiime --help
```

# Install picrust2

Follow
[this](https://github.com/picrust/picrust2/wiki/Manually-install-QIIME-2-plugin)
picrust2 instruction page. Run the following code in the Terminal and
type ‘y’ when prompted.

    conda install -c bioconda -c conda-forge picrust2

# Download & install q2-picrust2 plugin

``` bash
cd /home/shared/8TB_HDD_02/stanja
wget https://github.com/picrust/q2-picrust2/archive/refs/tags/2023.2_0.tar.gz
tar -xzvf 2023.2_0.tar.gz
```

``` bash
cd /home/shared/8TB_HDD_02/stanja/q2-picrust2-2023.2_0
pip install -e .
qiime dev refresh-cache
```

``` bash
qiime picrust2 --help
```

`Usage: qiime picrust2 [OPTIONS] COMMAND [ARGS]…`

`Description: This QIIME 2 plugin wraps the default 16S PICRUSt2 pipeline to run metagenome inference based on marker gene data. Currently only unstratified output is supported.`

`Plugin website: https://github.com/gavinmdouglas/q2-picrust2`

`Getting user support: Please post to the QIIME 2 forum for help with this plugin: https://forum.qiime2.org`

`Options: –version Show the version and exit. –example-data PATH Write example data and exit. –citations Show citations and exit. –help Show this message and exit.`

`Commands: custom-tree-pipeline 16S PICRUSt2 pipeline with custom tree full-pipeline Default 16S PICRUSt2 Pipeline`

> The required inputs are –i-table and –i-seq, which need to correspond
> to QIIME 2 artifacts of types FeatureTable\[Frequency\] and
> FeatureData\[Sequence\], respectively. The Feature Table needs to
> contain the abundances of ASVs (i.e. a BIOM table) and the sequence
> file needs to be a FASTA file containing the sequences for each ASV.

> These files correspond to the ASV count table, ASV sequences, and
> metadata for the samples.

# Feature table (BIOM)

located at path
`../input/241121_StonyCoral/270x200/250414_StonyCoral_270x200_featureTable_filtered.qza`

``` bash
qiime tools peek ../input/241121_StonyCoral/270x200/250414_StonyCoral_270x200_featureTable_filtered.qza
```

# Sequence file

located at path
`../input/241121_StonyCoral/270x200/250414_StonyCoral_270x200_representative-sequences.qza`

``` bash
qiime tools peek ../input/241121_StonyCoral/270x200/250414_StonyCoral_270x200_representative-sequences.qza
```

# PICRUSt2 qiime pipeline

``` bash
qiime picrust2 full-pipeline --help
```

``` bash
qiime picrust2 full-pipeline \
   --i-table ../input/241121_StonyCoral/270x200/250414_StonyCoral_270x200_featureTable.qza \
   --i-seq ../input/241121_StonyCoral/270x200/250414_StonyCoral_270x200_representative-sequences.qza \
   --output-dir ../q2-picrust2_output \
   --p-placement-tool sepp \
   --p-threads 4 \
   --p-hsp-method pic \
   --p-max-nsti 2 \
   --verbose
```

# 
