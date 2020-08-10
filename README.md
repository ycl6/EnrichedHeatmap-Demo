# EnrichedHeatmap Demo

### To visualize enrichment of genomic signals on targets or regions using EnrichedHeatmap

-----

**EnrichedHeatmap-Demo web version:** https://ycl6.github.io/EnrichedHeatmap-Demo

**Package URL:** https://github.com/jokergoo/EnrichedHeatmap

**Package Doc:** http://jokergoo.github.io/EnrichedHeatmap/articles/EnrichedHeatmap.html

**Demo Dataset:** [GSE80779](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80779) and [GSE80774](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80774) from *[ENL links histone acetylation to oncogenic gene expression in AML.](https://pubmed.ncbi.nlm.nih.gov/28241141/) Nature. 2017;543(7644):265-269*

**License:** GPL-3.0

## Introduction

In this tutorial, we will use some of the processed data from [GSE80779 ChIP-seq](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80779) and [GSE80774 RNA-seq](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80774) to demonstrate how to use EnrichedHeatmap to show enrichment profiles. In order to build the vignette fast, the data only includes chromosome 22 of human genome (hg19).

## Prerequisites

### Software & packages

- [R](https://www.r-project.org/) version 3.0 and above

This demo requires the following R packages

- [data.table](https://cran.r-project.org/package=data.table) from CRAN
- [circlize](https://cran.r-project.org/package=circlize) from CRAN
- [GenomicRanges](https://bioconductor.org/packages/GenomicRanges/) from Bioconductor
- [EnrichedHeatmap](https://www.bioconductor.org/packages/EnrichedHeatmap/) from Bioconductor

Optionally, [knitr](https://cran.r-project.org/package=knitr) from CRAN can be used to convert the R Markdown document into an HTML, PDF, or Word document

### Demo data used

| Data | Type | Description |
| --- | --- | --- |
| GSM2136938_MOLM13.H3K27ac | ChIP-seq | Enriched at active enhancers |
| GSM2136939_MOLM13.H3K4me1 | ChIP-seq | Enriched at enhancers (both active and poised) |
| GSM2136940_MOLM13.H3K4me3 | ChIP-seq | Enriched at actively transcribed promoters |
| GSM2136941_MOLM13.H3K9ac | ChIP-seq | Enriched at active enhancers |
| GSM2136954_MOLM13.POL2_sgGFP | ChIP-seq | Pol II, Control |
| GSM2136952_MOLM13.POL2-S2P_sgGFP | ChIP-seq | Pol II Ser2P, Control |
| GSM2136953_MOLM13.POL2_sgENL5 | ChIP-seq | Pol II - ENL KD |
| GSM2136933_MOLM13.CDK9_Control | ChIP-seq | CDK9, Control |
| GSM2136934_MOLM13.CDK9_ENLKO | ChIP-seq | CDK9 - ENL KD |
| GSE80774_MOLM13-sgGFP-DMSO | RNA-seq | Control |

### Other data used

| Data | Description |
| --- | --- |
| Gene annotation | From GENCODE v34lift37 |

## Download demo data

Download the result table from GitHub

```S
cd /ngs/EnrichedHeatmap-Demo

git clone https://github.com/ycl6/EnrichedHeatmap-Demo.git
```

# Go to online tutorial

https://ycl6.github.io/EnrichedHeatmap-Demo/
