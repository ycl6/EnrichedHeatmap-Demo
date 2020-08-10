#' ---
#' title: "Using the EnrichedHeatmap Bioconductor Package"
#' description: |
#'   Demonstration of using EnrichedHeatmap Bioconductor package
#' author:
#'   - name: "I-Hsuan Lin"
#'     url: https://github.com/ycl6
#'     affiliation: University of Manchester
#'     affiliation_url: https://www.manchester.ac.uk/
#' date: '`r format(Sys.Date(), "%B %d, %Y")`'
#' output:
#'     rmarkdown::html_document:
#'         theme: united
#'         highlight: tango
#'         self_contained: true
#'         toc: true
#'         toc_depth: 1
#'         toc_float:
#'             collapsed: false
#'             smooth_scroll: true
#' ---
#' 

#' 
#' -----
#' 
#' **Package URL:** https://github.com/jokergoo/EnrichedHeatmap
#' 
#' **Package Doc:** http://jokergoo.github.io/EnrichedHeatmap/articles/EnrichedHeatmap.html
#' 
#' **Demo Dataset:** [GSE80779](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80779) and [GSE80774](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80774) from *[ENL links histone acetylation to oncogenic gene expression in AML.](https://pubmed.ncbi.nlm.nih.gov/28241141/) Nature. 2017;543(7644):265-269*
#' 
#' **License:** GPL-3.0
#' 
#' # Introduction
#' 
#' In this tutorial, we will use some of the processed data from [GSE80779 ChIP-seq](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80779) and [GSE80774 RNA-seq](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80774) to demonstrate how to use EnrichedHeatmap to show enrichment profiles. In order to build the vignette fast, the data only includes chromosome 22 of human genome (hg19).
#' 
#' # Prerequisites
#' 
#' ## Software & packages
#' 
#' - [R](https://www.r-project.org/) version 3.0 and above
#' 
#' This demo requires the following R packages
#' 
#' - [data.table](https://cran.r-project.org/package=data.table)
#' - [R.utils](https://CRAN.R-project.org/package=R.utils)
#' - [circlize](https://cran.r-project.org/package=circlize)
#' - [edgeR](https://bioconductor.org/packages/edgeR/)
#' - [GenomicRanges](https://bioconductor.org/packages/GenomicRanges/)
#' - [EnrichedHeatmap](https://www.bioconductor.org/packages/EnrichedHeatmap/)
#' 
## ---- eval = FALSE------------------------------------------------------------
## # From CRAN
## install.packages(c("data.table", "R.utils", "circlize"))
## 
## # From Bioconductor
## if (!requireNamespace("BiocManager", quietly = TRUE))
## install.packages("BiocManager")
## BiocManager::install(c("edgeR", "GenomicRanges", "EnrichedHeatmap"))

#' 
#' ## Optional packages
#' 
#' - [rtracklayer](https://bioconductor.org/packages/rtracklayer/) from Bioconductor can be used to import annotation from GFF and GTF
#' - [knitr](https://cran.r-project.org/package=knitr) from CRAN can be used to convert the R Markdown document into an HTML, PDF, or Word document
#' 
#' ## Demo data used
#' 
#' | Data | Type | Description |
#' | ------ | -- | --------- |
#' | GSM2136938_MOLM13.H3K27ac | ChIP-seq | Enriched at active enhancers |
#' | GSM2136939_MOLM13.H3K4me1 | ChIP-seq | Enriched at enhancers (both active and poised) |
#' | GSM2136940_MOLM13.H3K4me3 | ChIP-seq | Enriched at actively transcribed promoters |
#' | GSM2136941_MOLM13.H3K9ac | ChIP-seq | Enriched at active enhancers |
#' | GSM2136954_MOLM13.POL2_sgGFP | ChIP-seq | Pol II, Control |
#' | GSM2136952_MOLM13.POL2-S2P_sgGFP | ChIP-seq | Pol II Ser2P, Control |
#' | GSM2136953_MOLM13.POL2_sgENL5 | ChIP-seq | Pol II - ENL KD |
#' | GSM2136933_MOLM13.CDK9_Control | ChIP-seq | CDK9, Control (super elongation complex (SEC) component) |
#' | GSM2136934_MOLM13.CDK9_ENLKO | ChIP-seq | CDK9 - ENL KD |
#' | GSE80774_MOLM13-sgGFP-DMSO | RNA-seq | Control |
#' 
#' ## Other data used
#' 
#' | Data | Description |
#' | ---- | ----------- |
#' | Gene annotation | From GENCODE v34lift37 |
#' 
#' # Download demo data
#' 
#' ```
#' git clone https://github.com/ycl6/EnrichedHeatmap-Demo.git
#' ```
#' 
#' -----
#' 
#' # Load package and set path
#' 
## ----load-libraries-----------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)   # Also requires R.utils to read gz and bz2 files
  library(GenomicRanges)
  library(EnrichedHeatmap)
  library(circlize)
})

#' 
#' Change the path to correspond to the location of the data files.
#' 
## ----set-path-----------------------------------------------------------------
data <- "data/"
list.files(data)

#' 
#' # Load data files
#' 
## ----load-files---------------------------------------------------------------
DT_H3K27ac <- fread(paste0(data, "GSM2136938_MOLM13.H3K27ac.chr22.bdg.gz"))
DT_H3K4me1 <- fread(paste0(data, "GSM2136939_MOLM13.H3K4me1.chr22.bdg.gz"))
DT_H3K4me3 <- fread(paste0(data, "GSM2136940_MOLM13.H3K4me3.chr22.bdg.gz"))
DT_H3K9ac <- fread(paste0(data, "GSM2136941_MOLM13.H3K9ac.chr22.bdg.gz"))

DT_POL2 <- fread(paste0(data, "GSM2136954_MOLM13.POL2_sgGFP.chr22.bdg.gz"))
DT_POL2ENLKD <- fread(paste0(data, "GSM2136953_MOLM13.POL2_sgENL5.chr22.bdg.gz"))
DT_POL2S2P <- fread(paste0(data, "GSM2136952_MOLM13.POL2-S2P_sgGFP.chr22.bdg.gz"))
DT_CDK9 <- fread(paste0(data, "GSM2136933_MOLM13.CDK9_Control.chr22.bdg.gz"))
DT_CDK9ENLKD <- fread(paste0(data, "GSM2136934_MOLM13.CDK9_ENLKO.chr22.bdg.gz"))

DT_genes <- fread(paste0(data, "gencode.chr22.txt.gz"))

rnaseq <- data.frame(fread(paste0(data, "GSE80774_MOLM13-sgGFP-DMSO.RNAseq.counttab.txt.gz")))
rownames(rnaseq) <- rnaseq$gene_id
colnames(rnaseq)[2:3] <- c("rep1","rep2")
rnaseq$mean <- apply(rnaseq[,2:3], 1, mean)

#' 
#' # Prepare data
#' 
## ----prepare-data-------------------------------------------------------------
# Calculate CPM using edgeR
rnaseq$cpm <- edgeR::cpm(rnaseq$mean)

# Subset rnaseq and genes to include genes present in both Data Frames
rnaseq <- rnaseq[rnaseq$gene_id %in% DT_genes$V4,]
DT_genes <- DT_genes[DT_genes$V4 %in% rnaseq$gene_id,]

# Reorder rows according to DT_genes
rnaseq <- rnaseq[DT_genes$V4,]

#' 
## -----------------------------------------------------------------------------
nrow(rnaseq)

#' 
## -----------------------------------------------------------------------------
nrow(DT_genes)

#' 
#' # *\*Import GFF & GTF annotation*
#' 
#' > **Note:** You can also import annotation from GFF or GTF using the `import.gff3` or `import` from the `rtracklayer` package.
#' 
#' ## Import GFF
#' 
## ----import-gff---------------------------------------------------------------
gff_file <- "gencode.v34lift37.annotation.gff3.gz"
gff_url <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/gencode.v34lift37.annotation.gff3.gz"

# Download GFF from GENCODE
download.file(gff_url, paste0(data, "/", gff_file))

# Create GRanges object from GFF
ggff <- rtracklayer::import.gff3("data/gencode.v34lift37.annotation.gff3.gz")

# To subset gene entries from chr22 from ggff
genes <- ggff[seqnames(ggff) == "chr22" & ggff$type == "gene"]
names(genes) <- substr(genes$gene_id, 1, 15)

genes

#' 
#' ## Import GTF
#' 
## ----import-gtf---------------------------------------------------------------
gtf_file <- "gencode.v34lift37.annotation.gtf.gz"
gtf_url <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/gencode.v34lift37.annotation.gtf.gz"

# Download GTF from GENCODE
download.file(gtf_url, paste0(data, "/", gtf_file))

# Create GRanges object from GTF
ggtf <- rtracklayer::import("data/gencode.v34lift37.annotation.gtf.gz")

# or from ggtf
genes <- ggtf[seqnames(ggtf) == "chr22" & ggtf$type == "gene"]
names(genes) <- substr(genes$gene_id, 1, 15)

genes

#' 
#' # Creating `GRanges` objects
#' 
## ----create-gRanges-----------------------------------------------------------
H3K27ac <- GRanges(seqnames = DT_H3K27ac$V1, 
        ranges = IRanges(start = DT_H3K27ac$V2+1, end = DT_H3K27ac$V3), 
        coverage = DT_H3K27ac$V4)

H3K4me1 <- GRanges(seqnames = DT_H3K4me1$V1, 
        ranges = IRanges(start = DT_H3K4me1$V2+1, end = DT_H3K4me1$V3), 
        coverage = DT_H3K4me1$V4)

H3K4me3 <- GRanges(seqnames = DT_H3K4me3$V1, 
        ranges = IRanges(start = DT_H3K4me3$V2+1, end = DT_H3K4me3$V3), 
        coverage = DT_H3K4me3$V4)

H3K9ac <- GRanges(seqnames = DT_H3K9ac$V1, 
        ranges = IRanges(start = DT_H3K9ac$V2+1, end = DT_H3K9ac$V3), 
        coverage = DT_H3K9ac$V4)

POL2 <- GRanges(seqnames = DT_POL2$V1, 
        ranges = IRanges(start = DT_POL2$V2+1, end = DT_POL2$V3), 
        coverage = DT_POL2$V4)

POL2ENLKD <- GRanges(seqnames = DT_POL2ENLKD$V1, 
        ranges = IRanges(start = DT_POL2ENLKD$V2+1, end = DT_POL2ENLKD$V3), 
        coverage = DT_POL2ENLKD$V4)

POL2S2P <- GRanges(seqnames = DT_POL2S2P$V1, 
        ranges = IRanges(start = DT_POL2S2P$V2+1, end = DT_POL2S2P$V3), 
        coverage = DT_POL2S2P$V4)

CDK9 <- GRanges(seqnames = DT_CDK9$V1, 
        ranges = IRanges(start = DT_CDK9$V2+1, end = DT_CDK9$V3), 
        coverage = DT_CDK9$V4)

CDK9ENLKD <- GRanges(seqnames = DT_CDK9ENLKD$V1, 
        ranges = IRanges(start = DT_CDK9ENLKD$V2+1, end = DT_CDK9ENLKD$V3), 
        coverage = DT_CDK9ENLKD$V4)

genes <- GRanges(seqnames = DT_genes$V1, 
        ranges = IRanges(start = DT_genes$V2+1, end = DT_genes$V3), 
        names = DT_genes$V4, strand = DT_genes$V5)
names(genes) <- DT_genes$V4

#' 
#' Let's have a look at `H3K27ac`.
#' 
## -----------------------------------------------------------------------------
H3K27ac

#' 
#' # Create TSS `GRanges` object
#' 
## ----create-tss---------------------------------------------------------------
tss <- promoters(genes, upstream = 0, downstream = 1)

#' 
## -----------------------------------------------------------------------------
genes

#' 
## -----------------------------------------------------------------------------
tss

#' 
#' # Create `normalizedMatrix` objects
#' 
#' Here we obtain the association between ChIP-seq peaks and targets (i.e. TSS) by normalizing into a matrix. We specify the targets regions to be extended to 5Kb upstream and downstream of TSS (`extend = 5000`). 
#' 
#' For each region, it splits into a list of 50 small windows (`w = 50`), then overlaps genomic signals to these small windows and calculates the value for every small window which is the mean value of genomic signals that intersects with the window. The value used for the calculation is controlled by `value_column` and how to calcualte the mean value is controlled by `mean_mode`. There are several modes for mean_mode according to different types of genomic signals, details can be found [here](https://www.bioconductor.org/packages/release/bioc/vignettes/EnrichedHeatmap/inst/doc/EnrichedHeatmap.html#toc_7).
#' 
#' ```
#'        40      50     20     values in signal regions (e.g. Coverage)
#'      ++++++   +++    +++++   signal regions (e.g. ChIP-seq peaks)
#'             30               values in signal region
#'           ++++++             signal region
#'        =================     a 17-bp window, and 4 bases did not overlap to any signal regions
#'        ****   ***    ***
#'          4     3      3      overlap bases
#'           ******
#'             6                overlap bases
#' ```
#' 
#' | Mode | Description | Calculation |
#' | -- | ------- | ------ |
#' | **absolute** | The mean of all signal regions regardless of their width | (40 + 30 + 50 + 20)/4 |
#' | **weighted** | The mean of all signal regions weighted by the width of their intersections | (40\*4 + 30\*6 + 50\*3 + 20\*3)/(4 + 6 + 3 + 3) |
#' | **w0** | The weighted mean between the intersected parts and un-intersected parts | (40\*4 + 30\*6 + 50\*3 + 20\*3)/(4 + 6 + 3 + 3 + 4) |
#' | **coverage** | The mean signal averged by the length of the window | (40\*4 + 30\*6 + 50\*3 + 20\*3)/17 |
#' 
## ----create-tss-matrix--------------------------------------------------------
tss_H3K27ac <-  normalizeToMatrix(H3K27ac, tss, value_column = "coverage", 
        extend = 5000, mean_mode = "w0", w = 50)

tss_H3K4me1 <- normalizeToMatrix(H3K4me1, tss, value_column = "coverage", 
        extend = 5000, mean_mode = "w0", w = 50)

tss_H3K4me3 <- normalizeToMatrix(H3K4me3, tss, value_column = "coverage", 
        extend = 5000, mean_mode = "w0", w = 50)

tss_H3K9ac <- normalizeToMatrix(H3K9ac, tss, value_column = "coverage", 
        extend = 5000, mean_mode = "w0", w = 50)

tss_POL2 <- normalizeToMatrix(POL2, tss, value_column = "coverage", 
        extend = 5000, mean_mode = "w0", w = 50)

#' 
#' Let's have a look at `tss_H3K27ac`.
#' 
## -----------------------------------------------------------------------------
tss_H3K27ac

#' 
#' How about the windows and values around TSS?
#' 
## -----------------------------------------------------------------------------
# Examine coverage in windows around TSS
tss_H3K27ac[95:100,91:110]

#' 
#' # Enrichment profile at TSS
#' 
#' Visualize the ChIP-seq profile around TSS by heatmap
#' 
## ----plot-tss-heatmap, fig.width = 4, fig.height = 6, fig.align = "center", dpi = 100, fig.cap = "H3K27ac profile (before color optimisation)"----
EnrichedHeatmap(tss_H3K27ac, col = c("white", "darkgreen"), name = "H3K27ac", 
        column_title = "H3K27ac coverage around TSS")

#' 
#' Output multiple images to PDF.
#' 
## ----plot-tss-heatmap-to-pdf--------------------------------------------------
pdf("EnrichedHeatmap_1.pdf", width = 4, height = 6, pointsize=12)
EnrichedHeatmap(tss_H3K27ac, col = c("white", "darkgreen"), name = "H3K27ac", 
        column_title = "H3K27ac coverage around TSS")
EnrichedHeatmap(tss_H3K4me1, col = c("white", "darkgreen"), name = "H3K4me1", 
        column_title = "H3K4me1 coverage around TSS")
EnrichedHeatmap(tss_H3K4me3, col = c("white", "darkgreen"), name = "H3K4me3", 
        column_title = "H3K4me3 coverage around TSS")
EnrichedHeatmap(tss_H3K9ac, col = c("white", "darkgreen"), name = "H3K9ac", 
        column_title = "H3K9ac coverage around TSS")
EnrichedHeatmap(tss_POL2, col = c("white", "blue"), name = "POL2", 
        column_title = "POL2 coverage around TSS")
invisible(dev.off())

#' 
#' # Tackle extreme values
#' 
#' If a vector of colors is specified, sequential values from minimal to maximal are mapped to the colors, and other values are linearly interpolated. Hence, extreme values can affect the color gradient.
#' 
#' ## Check distribution
#' 
## ----quantile-----------------------------------------------------------------
quantile(H3K27ac$coverage, c(0, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99, 1))

quantile(H3K4me1$coverage, c(0, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99, 1))

quantile(H3K4me3$coverage, c(0, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99, 1))

quantile(H3K9ac$coverage, c(0, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99, 1))

quantile(POL2$coverage, c(0, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99, 1))

#' 
#' ## Color mapping function
#' 
#' To get around of such extreme values, we define a color mapping function which only maps colors to values less than Nth percentile and the value larger than the Nth percentile uses same color as the 99th percentile.
#' 
## ----col-fun-1----------------------------------------------------------------
col1_fun1 <- colorRamp2(quantile(tss_H3K27ac, c(0, 0.99)), c("white", "darkgreen"))
col1_fun2 <- colorRamp2(quantile(tss_H3K4me1, c(0, 0.99)), c("white", "darkgreen"))
col1_fun3 <- colorRamp2(quantile(tss_H3K4me3, c(0, 0.95)), c("white", "darkgreen"))
col1_fun4 <- colorRamp2(quantile(tss_H3K9ac, c(0, 0.99)), c("white", "darkgreen"))
col1_fun5 <- colorRamp2(quantile(tss_POL2, c(0, 0.99)), c("white", "blue"))

#' 
## ----plot-tss-heatmap-col-fun, fig.width = 4, fig.height = 6, fig.align = "center", dpi = 100, fig.cap = "H3K27ac profile (after color optimisation)"----
EnrichedHeatmap(tss_H3K27ac, col = col1_fun1, name = "H3K27ac", 
        column_title = "H3K27ac coverage around TSS")

#' 
#' Output multiple images to PDF.
#' 
## ----plot-tss-heatmap-to-pdf-col-fun------------------------------------------
pdf("EnrichedHeatmap_2.pdf", width = 4, height = 6, pointsize=12)
EnrichedHeatmap(tss_H3K27ac, col = col1_fun1, name = "H3K27ac", 
                column_title = "H3K27ac coverage around TSS")
EnrichedHeatmap(tss_H3K4me1, col = col1_fun2, name = "H3K4me1", 
                column_title = "H3K4me1 coverage around TSS")
EnrichedHeatmap(tss_H3K4me3, col = col1_fun3, name = "H3K4me3", 
                column_title = "H3K4me3 coverage around TSS")
EnrichedHeatmap(tss_H3K9ac, col = col1_fun4, name = "H3K9ac", 
                column_title = "H3K9ac coverage around TSS")
EnrichedHeatmap(tss_POL2, col = col1_fun5, name = "POL2", 
                column_title = "POL2 coverage around TSS")
invisible(dev.off())

#' 
#' # Visualize multiple heatmaps
#' 
#' We can place multiple heatmaps into a single image.
#' 
## ----plot-multiple-heatmaps---------------------------------------------------
ht_list <- EnrichedHeatmap(tss_H3K27ac, col = col1_fun1, name = "H3K27ac", 
                          column_title = "H3K27ac") +
EnrichedHeatmap(tss_H3K4me1, col = col1_fun2, name = "H3K4me1", column_title = "H3K4me1") +
EnrichedHeatmap(tss_H3K4me3, col = col1_fun3, name = "H3K4me3", column_title = "H3K4me3") +
EnrichedHeatmap(tss_H3K9ac, col = col1_fun4, name = "H3K9ac", column_title = "H3K9ac") +
EnrichedHeatmap(tss_POL2, col = col1_fun5, name = "POL2", column_title = "POL2") +
Heatmap(log2(rnaseq$cpm+1), col = c("white", "black"), name = "log2(CPM+1)", 
        show_row_names = FALSE, show_column_names = FALSE, width = unit(10, "mm"))

png("MultipleHeatmaps_1.png", width = 10, height = 10, units = "in", res = 300)
draw(ht_list, ht_gap = unit(c(7, 7, 7, 7, 7), "mm"))
invisible(dev.off())

#' 
## ---- fig.width = 10, fig.height = 10, fig.align = "center", dpi = 100, fig.cap = "Multiple heatmaps 1"----
include_graphics("MultipleHeatmaps_1.png")

#' 
#' # Enrichment profile at gene-body
#' 
#' We can also view enrichment profile over a genomic region. 
#' 
#' Here, we will create normalized matrix of POL2 and CDK9 ChIP-seq signals over gene body. The Width of the target regions shown on heatmap can be controlled by `target_ratio` which is relative to the width of the complete heatmap.
#' 
## ----create-gene-body-matrix--------------------------------------------------
gb_POL2 <- normalizeToMatrix(POL2, genes, value_column = "coverage", 
        extend = 5000, mean_mode = "w0", w = 50, target_ratio = 0.5)

gb_POL2ENLKD <- normalizeToMatrix(POL2ENLKD, genes, value_column = "coverage", 
        extend = 5000, mean_mode = "w0", w = 50, target_ratio = 0.5)

gb_POL2S2P <- normalizeToMatrix(POL2S2P, genes, value_column = "coverage", 
        extend = 5000, mean_mode = "w0", w = 50, target_ratio = 0.5)

gb_CDK9 <- normalizeToMatrix(CDK9, genes, value_column = "coverage", 
        extend = 5000, mean_mode = "w0", w = 50, target_ratio = 0.5)

gb_CDK9ENLKD <- normalizeToMatrix(CDK9ENLKD, genes, value_column = "coverage", 
        extend = 5000, mean_mode = "w0", w = 50, target_ratio = 0.5)

#' 
#' Check distribution.
#' 
## ----quantile2----------------------------------------------------------------
quantile(POL2$coverage, c(0, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99, 1))

quantile(POL2ENLKD$coverage, c(0, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99, 1))

quantile(POL2S2P$coverage, c(0, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99, 1))

quantile(CDK9$coverage, c(0, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99, 1))

quantile(CDK9ENLKD$coverage, c(0, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99, 1))

#' 
#' Define color mapping functions.
#' 
## ----col-fun-2----------------------------------------------------------------
col2_fun1 = colorRamp2(quantile(gb_POL2, c(0, 0.99)), c("white", "blue"))
col2_fun2 = colorRamp2(quantile(gb_POL2ENLKD, c(0, 0.99)), c("white", "blue"))
col2_fun3 = colorRamp2(quantile(gb_POL2S2P, c(0, 0.99)), c("white", "blue"))
col2_fun4 = colorRamp2(quantile(gb_CDK9, c(0, 0.99)), c("white", "red"))
col2_fun5 = colorRamp2(quantile(gb_CDK9ENLKD, c(0, 0.99)), c("white", "red"))

#' 
#' We will define 3 partitions according to `kmeans` clustering of the expression profile, and split the heatmap accordingly.
#' 
## ----plot-multiple-heatmaps-with-partitions-----------------------------------
set.seed(123)
partition <- paste0("cluster", kmeans(log2(rnaseq$cpm+1), centers = 3)$cluster)
lgd = Legend(at = c("cluster1", "cluster2", "cluster3"), title = "Clusters", 
        type = "lines", legend_gp = gpar(col = 3:5))

axis_name <- c("-5kb", "TSS", "TTS", "+5kb")

ht_list <- Heatmap(partition, col = structure(3:5, names = paste0("cluster", 1:3)), 
        name = "Partition", show_row_names = FALSE, show_column_names = FALSE, 
        width = unit(3, "mm")) +
Heatmap(log2(rnaseq$cpm+1), col = c("white", "black"), name = "log2(CPM+1)", 
        show_row_names = FALSE, show_column_names = FALSE, width = unit(15, "mm"), 
        top_annotation = HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 3:5), 
        outline = FALSE, axis_param = list(side = "left")))) +
EnrichedHeatmap(gb_POL2, col = col2_fun1, name = "POL2", column_title = "POL2", 
        axis_name = axis_name, axis_name_gp = gpar(fontsize = 8), 
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 3:5), 
        axis_param = list(side = "left")))) +
EnrichedHeatmap(gb_POL2ENLKD, col = col2_fun2, name = "POL2 ENLKD", column_title = "POL2 ENL-KD", 
        axis_name = axis_name, axis_name_gp = gpar(fontsize = 8), 
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 3:5), 
        axis_param = list(side = "left")))) +
EnrichedHeatmap(gb_POL2S2P, col = col2_fun3, name = "POL2 S2P", column_title = "POL2 S2P", 
        axis_name = axis_name, axis_name_gp = gpar(fontsize = 8), 
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 3:5), 
        axis_param = list(side = "left")))) +
EnrichedHeatmap(gb_CDK9, col = col2_fun4, name = "CDK9", column_title = "CDK9", 
        axis_name = axis_name, axis_name_gp = gpar(fontsize = 8), 
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 3:5), 
        axis_param = list(side = "left")))) +
EnrichedHeatmap(gb_CDK9ENLKD, col = col2_fun5, name = "CDK9 ENLKD", column_title = "CDK9 ENL-KD", 
        axis_name = axis_name, axis_name_gp = gpar(fontsize = 8), 
        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 3:5), 
        axis_param = list(side = "left"))))

png("MultipleHeatmaps_2.png", width = 12, height = 10, units = "in", res = 300)
draw(ht_list, split = partition, ht_gap = unit(c(2, 7, 7, 7, 7, 7), "mm"))
invisible(dev.off())

#' 
## ---- fig.width = 12, fig.height = 10, fig.align = "center", dpi = 100, fig.cap = "Multiple heatmaps 2"----
include_graphics("MultipleHeatmaps_2.png")

#' 
#' # Session information
#' 
## ----session-info-------------------------------------------------------------
sessionInfo()

