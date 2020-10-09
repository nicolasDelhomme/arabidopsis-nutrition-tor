#' ---
#' title: "Biological QA after resequencing of all the batches"
#' author: "Thomas Dobrenel and Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' * Set the working dir
setwd("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/all-batches")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/all-batches")
#' ```

#' * Load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(DESeq2)
    library(gplots)
    library(here)
    library(hyperSpec)
    library(limma)
    library(LSD)
    library(magrittr)
    library(matrixStats)
    library(parallel)
    library(pander)
    library(plotly)
    library(RColorBrewer)
    library(scatterplot3d)
    library(tidyverse)
    library(tximport)
    library(VennDiagram)
    library(vsn)
})

#' * Source some helper functions
suppressPackageStartupMessages({
    source("~/Git/UPSCb/UPSCb-common/src/R/featureSelection.R")
    source("~/Git/UPSCb/UPSCb-common/src/R/gopher.R")
    source("~/Git/UPSCb/UPSCb-common/src/R/plot.multidensity.R")
    source("~/Git/UPSCb/UPSCb-common/src/R/volcanoPlot.R")
})

#' * Cutoff
#' As of Schurch et al., RNA, 2016
lfc <- 0.5
FDR <- 0.01

#' * Create palettes
pal <- c(brewer.pal(8,"Dark2"),1)
pal2 <- brewer.pal(9,"Paired") #require package RColorBrewer
cols <- rainbow(17)
hpal <- colorRampPalette(c("blue","white","red"))(100)

#' * Register the default plot margin
mar <- par("mar")

#' # Raw data
#' ## Loading of sample information
#' * Read the sample information
samples <- read.csv("~/Git/UPSCb/projects/arabidopsis-nutrition-TOR/doc/samples3.csv")

#' * Remove unnecessary samples
samples %<>% filter(!grepl("P11554_1",SciLifeID)) %>% 
    filter(! SciLifeID %in% c("P13406_101",
                        "P14066_128",
                        "P14066_133",
                        "P13406_102",
                        "P14066_131")) %>%
    mutate(Nutrition,Nutrition=relevel(Nutrition,"NPS")) %>% 
    mutate(AZD,AZD=relevel(AZD,"DMSO"))

samples <- samples[order(samples$Timepoint, samples$Nutrition, samples$AZD),]

samples %<>% mutate(Timepoint,Timepoint=factor(paste0("T",Timepoint)))

#' ### Load the data
#' * Call the data
filenames <- list.files("../Salmon", 
                    recursive = TRUE, 
                    pattern = "quant.sf",
                    full.names = TRUE)

#' * Name the data
names(filenames) <- sub("_S.*","",sapply(strsplit(filenames, "/"), .subset, 3))

#' * Match  data <=> sample list
filenames <- filenames[match(samples$SciLifeID,names(filenames))]
filenames <-filenames[!is.na(names(filenames))]
samples <- samples[match(names(filenames),samples$SciLifeID),]

#' * Annotate the samples
samples$Conditions <- factor(paste(samples$Timepoint,samples$Nutrition,samples$AZD,sep="_"))
samples$Batch <- factor(substr(samples$SciLifeID,1,8))

#' # Expression data
#' Read the expression at the transcript level
tx <- suppressMessages(tximport(files = filenames, 
                                type = "salmon", 
                                txOut = TRUE))

#' summarise to genes
tx2gene <- data.frame(TXID=rownames(tx$counts),
                      GENEID=sub("\\.[0-9]+","",rownames(tx$counts)))
gx <- summarizeToGene(tx,tx2gene=tx2gene)

kg <- round(gx$counts) 

#' Sanity check
stopifnot(all(colnames(kg) == samples$SciLifeID))

#' ## Raw data export
#dir.create(file.path("analysis_Tom","salmon"),showWarnings=FALSE,recursive=TRUE)
#write.table(kg,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/analysis_Tom/nutrition-unormalised-gene-expression_data.csv")
#save(kg, samples, file = "/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/analysis_Tom/counts.rda")

#' ## Preliminary validations
#' ### Check for the genes that are never expressed
sel <- rowSums(kg) == 0 
sprintf("%s%% (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(kg),digits=1),
        sum(sel),
        nrow(kg))

#' # Data normalisation 
#'  For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample type 
dds <- DESeqDataSetFromMatrix(
    countData = kg,
    colData = samples,
    design = ~ Conditions)
dds <- estimateSizeFactors(dds)

#' ## Perform a Variance Stabilizing Transformation for plotting
vst <- varianceStabilizingTransformation(dds,blind=FALSE)
vsd <- assay(vst)
vsd <- vsd - min(vsd)

#' Write out
#write.csv(vst,"/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/analysis_Tom/library-size-normalized_variance-stabilized_data_nutrition.csv")


#' # Multivariate analysis
#' ## PCA
#' Principal Component Analysis on the normalized data
#' * Establishment of the PCA
conditions1 <- factor(paste(samples$AZD,samples$Nutrition,sep="_"))
pc <- prcomp(t(vsd))
percent <- round(summary(pc)$importance[2,]*100);percent

#' * Graphical representation of PC1 x PC2
conds <- droplevels(conditions1)
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(conds)],
     pch=c(15,16,17)[as.factor(samples$Timepoint)],
     main="All timepoints")

legend("bottomright",pch=19,
       col=pal[1:nlevels(conds)],
       legend=levels(conds))

legend("topleft",pch=c(15,16,17),
       col="black",
       legend=c("T0","T6","T24"))

#' * Graphical representation of PC2 x PC3
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal[as.integer(conds)],
     pch=c(15,16,17)[as.factor(samples$Timepoint)],
     main="All timepoints")

legend("topleft",pch=19,
       col=pal[1:nlevels(conds)],
       legend=levels(conds))

legend("bottomleft",pch=c(15,16,17),
       col="black",
       legend=c("T0","T6","T24"))

#' * Graphical representation of PC1 x PC3
plot(pc$x[,1],
     pc$x[,3],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal[as.integer(conds)],
     pch=c(15,16,17)[as.factor(samples$Timepoint)],
     main="All timepoints")

legend("bottomleft",pch=19,
       col=pal[1:nlevels(conds)],
       legend=levels(conds))

legend("bottomright",pch=c(15,16,17),
       col="black",
       legend=c("T0","T6","T24"))

#' # Differential expression based on all conditions compared to the T0
#' ## Preparation of the design
sel <- samples$Timepoint %in% c("T0","T6","T24")
suppressMessages(dds <- DESeqDataSetFromMatrix(
    countData = kg[,sel],
    colData = samples[sel,],
    design = ~ Conditions))

#' ## Differential expression analysis
dds <- DESeq(dds)

#' ## Variance Stabilising Transformation
#' ### Perform a Variance Stabilizing Transformation for plotting
vst <- varianceStabilizingTransformation(dds,blind=FALSE)
vsd <- assay(vst)
vsd <- vsd - min(vsd)

#' ## Contrasts to T0
#' The contrast by default is the first one (not Intercept)
resultsNames(dds)

#' ### Different effects at T6
#' #### T6_NPS_DMSO
res <- results(dds,name = "Conditions_T6_NPS_DMSO_vs_T0_T0_0")

cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
T6NPSDMSO <- rownames(res[cutoffs,])
T6NPSDMSO_Low <- rownames(res[cutoff2,])
T6NPSDMSO_High <- rownames(res[cutoff1,])

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))

#' #### MA plot
DESeq2::plotMA(res)

#' #### Volcano plot 
volcanoPlot(res)

#' #### Heatmap
hm <- heatmap.2(t(scale(t(vsd[cutoffs,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))

#' #### T6_NPS_AZD
res <- results(dds,name = "Conditions_T6_NPS_AZD_vs_T0_T0_0")

cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
T6NPSAZD <- rownames(res[cutoffs,])
T6NPSAZD_Low <- rownames(res[cutoff2,])
T6NPSAZD_High <- rownames(res[cutoff1,])

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))

#' #### MA plot
DESeq2::plotMA(res)

#' #### Volcano plot 
volcanoPlot(res)

#' #### Heatmap
hm <- heatmap.2(t(scale(t(vsd[cutoffs,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))

#' #### T6_PKS_DMSO
res <- results(dds,name = "Conditions_T6_PKS_DMSO_vs_T0_T0_0")

cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
T6PKSDMSO <- rownames(res[cutoffs,])
T6PKSDMSO_Low <- rownames(res[cutoff2,])
T6PKSDMSO_High <- rownames(res[cutoff1,])

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))

#' #### MA plot
DESeq2::plotMA(res)

#' #### Volcano plot 
volcanoPlot(res)

#' #### Heatmap
hm <- heatmap.2(t(scale(t(vsd[cutoffs,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))

#' #### T6_PKS_AZD
res <- results(dds,name = "Conditions_T6_PKS_AZD_vs_T0_T0_0")

cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
T6PKSAZD <- rownames(res[cutoffs,])
T6PKSAZD_Low <- rownames(res[cutoff2,])
T6PKSAZD_High <- rownames(res[cutoff1,])

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))

#' #### MA plot
DESeq2::plotMA(res)

#' #### Volcano plot 
volcanoPlot(res)

#' #### Heatmap
hm <- heatmap.2(t(scale(t(vsd[cutoffs,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))

#' #### T6_NS_DMSO
res <- results(dds,name = "Conditions_T6_NS_DMSO_vs_T0_T0_0")

cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
T6NSDMSO <- rownames(res[cutoffs,])
T6NSDMSO_Low <- rownames(res[cutoff2,])
T6NSDMSO_High <- rownames(res[cutoff1,])

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))

#' #### MA plot
DESeq2::plotMA(res)

#' #### Volcano plot 
volcanoPlot(res)

#' #### Heatmap
hm <- heatmap.2(t(scale(t(vsd[cutoffs,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))

#' #### T6_NS_AZD
res <- results(dds,name = "Conditions_T6_NS_AZD_vs_T0_T0_0")

cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
T6NSAZD <- rownames(res[cutoffs,])
T6NSAZD_Low <- rownames(res[cutoff2,])
T6NSAZD_High <- rownames(res[cutoff1,])

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))

#' #### MA plot
DESeq2::plotMA(res)

#' #### Volcano plot 
volcanoPlot(res)

#' #### Heatmap
hm <- heatmap.2(t(scale(t(vsd[cutoffs,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))

#' #### T6_NP_DMSO
res <- results(dds,name = "Conditions_T6_NP_DMSO_vs_T0_T0_0")

cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
T6NPDMSO <- rownames(res[cutoffs,])
T6NPDMSO_Low <- rownames(res[cutoff2,])
T6NPDMSO_High <- rownames(res[cutoff1,])

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))

#' #### MA plot
DESeq2::plotMA(res)

#' #### Volcano plot 
volcanoPlot(res)

#' #### Heatmap
hm <- heatmap.2(t(scale(t(vsd[cutoffs,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))

#' #### T6_NP_AZD
res <- results(dds,name = "Conditions_T6_NP_AZD_vs_T0_T0_0")

cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
T6NPAZD <- rownames(res[cutoffs,])
T6NPAZD_Low <- rownames(res[cutoff2,])
T6NPAZD_High <- rownames(res[cutoff1,])

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))

#' #### MA plot
DESeq2::plotMA(res)

#' #### Volcano plot 
volcanoPlot(res)

#' #### Heatmap
hm <- heatmap.2(t(scale(t(vsd[cutoffs,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))






#' ### Different effects at T24
#' #### T24_NPS_DMSO
res <- results(dds,name = "Conditions_T24_NPS_DMSO_vs_T0_T0_0")

cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
T24NPSDMSO <- rownames(res[cutoffs,])
T24NPSDMSO_Low <- rownames(res[cutoff2,])
T24NPSDMSO_High <- rownames(res[cutoff1,])

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))

#' #### MA plot
DESeq2::plotMA(res)

#' #### Volcano plot 
volcanoPlot(res)

#' #### Heatmap
hm <- heatmap.2(t(scale(t(vsd[cutoffs,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))

#' #### T24_NPS_AZD
res <- results(dds,name = "Conditions_T24_NPS_AZD_vs_T0_T0_0")

cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
T24NPSAZD <- rownames(res[cutoffs,])
T24NPSAZD_Low <- rownames(res[cutoff2,])
T24NPSAZD_High <- rownames(res[cutoff1,])

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))

#' #### MA plot
DESeq2::plotMA(res)

#' #### Volcano plot 
volcanoPlot(res)

#' #### Heatmap
hm <- heatmap.2(t(scale(t(vsd[cutoffs,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))

#' #### T24_PKS_DMSO
res <- results(dds,name = "Conditions_T24_PKS_DMSO_vs_T0_T0_0")

cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
T24PKSDMSO <- rownames(res[cutoffs,])
T24PKSDMSO_Low <- rownames(res[cutoff2,])
T24PKSDMSO_High <- rownames(res[cutoff1,])

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))

#' #### MA plot
DESeq2::plotMA(res)

#' #### Volcano plot 
volcanoPlot(res)

#' #### Heatmap
hm <- heatmap.2(t(scale(t(vsd[cutoffs,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))

#' #### T24_PKS_AZD
res <- results(dds,name = "Conditions_T24_PKS_AZD_vs_T0_T0_0")

cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
T24PKSAZD <- rownames(res[cutoffs,])
T24PKSAZD_Low <- rownames(res[cutoff2,])
T24PKSAZD_High <- rownames(res[cutoff1,])

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))

#' #### MA plot
DESeq2::plotMA(res)

#' #### Volcano plot 
volcanoPlot(res)

#' #### Heatmap
hm <- heatmap.2(t(scale(t(vsd[cutoffs,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))

#' #### T24_NS_DMSO
res <- results(dds,name = "Conditions_T24_NS_DMSO_vs_T0_T0_0")

cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
T24NSDMSO <- rownames(res[cutoffs,])
T24NSDMSO_Low <- rownames(res[cutoff2,])
T24NSDMSO_High <- rownames(res[cutoff1,])

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))

#' #### MA plot
DESeq2::plotMA(res)

#' #### Volcano plot 
volcanoPlot(res)

#' #### Heatmap
hm <- heatmap.2(t(scale(t(vsd[cutoffs,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))

#' #### T24_NS_AZD
res <- results(dds,name = "Conditions_T24_NS_AZD_vs_T0_T0_0")

cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
T24NSAZD <- rownames(res[cutoffs,])
T24NSAZD_Low <- rownames(res[cutoff2,])
T24NSAZD_High <- rownames(res[cutoff1,])

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))

#' #### MA plot
DESeq2::plotMA(res)

#' #### Volcano plot 
volcanoPlot(res)

#' #### Heatmap
hm <- heatmap.2(t(scale(t(vsd[cutoffs,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))

#' #### T24_NP_DMSO
res <- results(dds,name = "Conditions_T24_NP_DMSO_vs_T0_T0_0")

cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
T24NPDMSO <- rownames(res[cutoffs,])
T24NPDMSO_Low <- rownames(res[cutoff2,])
T24NPDMSO_High <- rownames(res[cutoff1,])

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))

#' #### MA plot
DESeq2::plotMA(res)

#' #### Volcano plot 
volcanoPlot(res)

#' #### Heatmap
hm <- heatmap.2(t(scale(t(vsd[cutoffs,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))

#' #### T24_NP_AZD
res <- results(dds,name = "Conditions_T24_NP_AZD_vs_T0_T0_0")

cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
T24NPAZD <- rownames(res[cutoffs,])
T24NPAZD_Low <- rownames(res[cutoff2,])
T24NPAZD_High <- rownames(res[cutoff1,])

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))

#' #### MA plot
DESeq2::plotMA(res)

#' #### Volcano plot 
volcanoPlot(res)

#' #### Heatmap
hm <- heatmap.2(t(scale(t(vsd[cutoffs,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))










#' # Comparisons of GOI lists (Venn diagrams)
#' ## GOI of nutrition at T6
#' * All GOI
grid.draw(venn.diagram(list(T6NPSDMSO, T6PKSDMSO, T6NSDMSO, T6NPDMSO),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("Full","No N","No P", "No C")))

#' * Downregulated genes
grid.draw(venn.diagram(list(T6NPSDMSO_Low, T6PKSDMSO_Low, T6NSDMSO_Low, T6NPDMSO_Low),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("Full","No N","No P", "No C")))

#' * Upregulated genes
grid.draw(venn.diagram(list(T6NPSDMSO_High, T6PKSDMSO_High, T6NSDMSO_High, T6NPDMSO_High),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("Full","No N","No P", "No C")))

#' ## GOI of AZD at T6
#' * All GOI
grid.draw(venn.diagram(list(T6NPSAZD, T6PKSAZD, T6NSAZD, T6NPAZD),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("Full","No N","No P", "No C")))

#' * Downregulated genes
grid.draw(venn.diagram(list(T6NPSAZD_Low, T6PKSAZD_Low, T6NSAZD_Low, T6NPAZD_Low),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("Full","No N","No P", "No C")))

#' * Upregulated genes
grid.draw(venn.diagram(list(T6NPSAZD_High, T6PKSAZD_High, T6NSAZD_High, T6NPAZD_High),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("Full","No N","No P", "No C")))

#' ## Nutrition and AZD interaction at T6
#' ### For Carbon starvation
#' * Downregulated genes
grid.draw(venn.diagram(list(T6NPSAZD_Low, T6NPSDMSO_Low, T6NPAZD_Low, T6NPDMSO_Low),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("AZD","Ctrl","-Suc x AZD", "-Suc")))

#' * Upregulated genes
grid.draw(venn.diagram(list(T6NPSAZD_High, T6NPSDMSO_High, T6NPAZD_High, T6NPDMSO_High),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("AZD","Ctrl","-Suc x AZD", "-Suc")))

#' ### For Phosphorus starvation
#' * Downregulated genes
grid.draw(venn.diagram(list(T6NPSAZD_Low, T6NPSDMSO_Low, T6NSAZD_Low, T6NSDMSO_Low),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("AZD","Ctrl","-Phos x AZD", "-Phos")))

#' * Upregulated genes
grid.draw(venn.diagram(list(T6NPSAZD_High, T6NPSDMSO_High, T6NSAZD_High, T6NSDMSO_High),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("AZD","Ctrl","-Phos x AZD", "-Phos")))

#' ### For Nitrogen starvation
#' * Downregulated genes
grid.draw(venn.diagram(list(T6NPSAZD_Low, T6NPSDMSO_Low, T6PKSAZD_Low, T6PKSDMSO_Low),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("AZD","Ctrl","-Nit x AZD", "-Nit")))

#' * Upregulated genes
grid.draw(venn.diagram(list(T6NPSAZD_High, T6NPSDMSO_High, T6PKSAZD_High, T6PKSDMSO_High),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("AZD","Ctrl","-Nit x AZD", "-Nit")))

#' ## Nutrition and AZD interaction at T24
#' ### For Carbon starvation
#' * Downregulated genes
grid.draw(venn.diagram(list(T24NPSAZD_Low, T24NPSDMSO_Low, T24NPAZD_Low, T24NPDMSO_Low),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("AZD","Ctrl","-Suc x AZD", "-Suc")))

#' * Upregulated genes
grid.draw(venn.diagram(list(T24NPSAZD_High, T24NPSDMSO_High, T24NPAZD_High, T24NPDMSO_High),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("AZD","Ctrl","-Suc x AZD", "-Suc")))

#' ### For Phosphorus starvation
#' * Downregulated genes
grid.draw(venn.diagram(list(T24NPSAZD_Low, T24NPSDMSO_Low, T24NSAZD_Low, T24NSDMSO_Low),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("AZD","Ctrl","-Phos x AZD", "-Phos")))

#' * Upregulated genes
grid.draw(venn.diagram(list(T24NPSAZD_High, T24NPSDMSO_High, T24NSAZD_High, T24NSDMSO_High),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("AZD","Ctrl","-Phos x AZD", "-Phos")))

#' ### For Nitrogen starvation
#' * Downregulated genes
grid.draw(venn.diagram(list(T24NPSAZD_Low, T24NPSDMSO_Low, T24PKSAZD_Low, T24PKSDMSO_Low),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("AZD","Ctrl","-Nit x AZD", "-Nit")))

#' * Upregulated genes
grid.draw(venn.diagram(list(T24NPSAZD_High, T24NPSDMSO_High, T24PKSAZD_High, T24PKSDMSO_High),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("AZD","Ctrl","-Nit x AZD", "-Nit")))

#' # Export list of DEGs
write.table(T6NPSDMSO_Low, file = "T6NPSDMSO_Low.txt", sep = "\t", row.names = FALSE)
write.table(T6PKSDMSO_Low, file = "T6PKSDMSO_Low.txt", sep = "\t", row.names = FALSE)
write.table(T6NSDMSO_Low, file = "T6NSDMSO_Low.txt", sep = "\t", row.names = FALSE)
write.table(T6NPDMSO_Low, file = "T6NPDMSO_Low.txt", sep = "\t", row.names = FALSE)

write.table(T6NPSDMSO_High, file = "T6NPSDMSO_High.txt", sep = "\t", row.names = FALSE)
write.table(T6PKSDMSO_High, file = "T6PKSDMSO_High.txt", sep = "\t", row.names = FALSE)
write.table(T6NSDMSO_High, file = "T6NSDMSO_High.txt", sep = "\t", row.names = FALSE)
write.table(T6NPDMSO_High, file = "T6NPDMSO_High.txt", sep = "\t", row.names = FALSE)

write.table(T6NPSAZD_Low, file = "T6NPSAZD_Low.txt", sep = "\t", row.names = FALSE)
write.table(T6PKSAZD_Low, file = "T6PKSAZD_Low.txt", sep = "\t", row.names = FALSE)
write.table(T6NSAZD_Low, file = "T6NSAZD_Low.txt", sep = "\t", row.names = FALSE)
write.table(T6NPAZD_Low, file = "T6NPAZD_Low.txt", sep = "\t", row.names = FALSE)

write.table(T6NPSAZD_High, file = "T6NPSAZD_High.txt", sep = "\t", row.names = FALSE)
write.table(T6PKSAZD_High, file = "T6PKSAZD_High.txt", sep = "\t", row.names = FALSE)
write.table(T6NSAZD_High, file = "T6NSAZD_High.txt", sep = "\t", row.names = FALSE)
write.table(T6NPAZD_High, file = "T6NPAZD_High.txt", sep = "\t", row.names = FALSE)


write.table(T24NPSDMSO_Low, file = "T24NPSDMSO_Low.txt", sep = "\t", row.names = FALSE)
write.table(T24PKSDMSO_Low, file = "T24PKSDMSO_Low.txt", sep = "\t", row.names = FALSE)
write.table(T24NSDMSO_Low, file = "T24NSDMSO_Low.txt", sep = "\t", row.names = FALSE)
write.table(T24NPDMSO_Low, file = "T24NPDMSO_Low.txt", sep = "\t", row.names = FALSE)

write.table(T24NPSDMSO_High, file = "T24NPSDMSO_High.txt", sep = "\t", row.names = FALSE)
write.table(T24PKSDMSO_High, file = "T24PKSDMSO_High.txt", sep = "\t", row.names = FALSE)
write.table(T24NSDMSO_High, file = "T24NSDMSO_High.txt", sep = "\t", row.names = FALSE)
write.table(T24NPDMSO_High, file = "T24NPDMSO_High.txt", sep = "\t", row.names = FALSE)

write.table(T24NPSAZD_Low, file = "T24NPSAZD_Low.txt", sep = "\t", row.names = FALSE)
write.table(T24PKSAZD_Low, file = "T24PKSAZD_Low.txt", sep = "\t", row.names = FALSE)
write.table(T24NSAZD_Low, file = "T24NSAZD_Low.txt", sep = "\t", row.names = FALSE)
write.table(T24NPAZD_Low, file = "T24NPAZD_Low.txt", sep = "\t", row.names = FALSE)

write.table(T24NPSAZD_High, file = "T24NPSAZD_High.txt", sep = "\t", row.names = FALSE)
write.table(T24PKSAZD_High, file = "T24PKSAZD_High.txt", sep = "\t", row.names = FALSE)
write.table(T24NSAZD_High, file = "T24NSAZD_High.txt", sep = "\t", row.names = FALSE)
write.table(T24NPAZD_High, file = "T24NPAZD_High.txt", sep = "\t", row.names = FALSE)

#' # GO analysis
#' ## T6 DMSO
#' ### Induced genes
GO <- gopher(T6NPSDMSO_High,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana")
GO$go$name
GO$kegg$id
GO$pfam$name
write.table(data.frame(GO$go$id, GO$go$padj), file = "GO_T6NPSDMSO_High.txt", sep = "\t",
            row.names = FALSE)
GO_T6NPSDMSO_High <- GO

GO <- gopher(T6NSDMSO_High,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana")
GO$go$name
GO$kegg$id
GO$pfam$name
write.table(data.frame(GO$go$id, GO$go$padj), file = "GO_T6NSDMSO_High.txt", sep = "\t",
            row.names = FALSE)
GO_T6NSDMSO_High <- GO

GO <- gopher(T6PKSDMSO_High,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana")
GO$go$name
GO$kegg$id
GO$pfam$name
write.table(data.frame(GO$go$id, GO$go$padj), file = "GO_T6PKSDMSO_High.txt", sep = "\t",
            row.names = FALSE)
GO_T6PKSDMSO_High <- GO

GO <- gopher(T6NPDMSO_High,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana")
GO$go$name
GO$kegg$id
GO$pfam$name
write.table(data.frame(GO$go$id, GO$go$padj), file = "GO_T6NPDMSO_High.txt", sep = "\t",
            row.names = FALSE)
GO_T6NPDMSO_High <- GO

#' ### Repressed genes
GO <- gopher(T6NPSDMSO_Low,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana")
GO$go$name
GO$kegg$id
GO$pfam$name
write.table(data.frame(GO$go$id, GO$go$padj), file = "GO_T6NPSDMSO_Low.txt", sep = "\t",
            row.names = FALSE)
GO_T6NPSDMSO_Low <- GO

GO <- gopher(T6NSDMSO_Low,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana")
GO$go$name
GO$kegg$id
GO$pfam$name
write.table(data.frame(GO$go$id, GO$go$padj), file = "GO_T6NSDMSO_Low.txt", sep = "\t",
            row.names = FALSE)
GO_T6NSDMSO_Low <- GO

GO <- gopher(T6PKSDMSO_Low,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana")
GO$go$name
GO$kegg$id
GO$pfam$name
write.table(data.frame(GO$go$id, GO$go$padj), file = "GO_T6PKSDMSO_Low.txt", sep = "\t",
            row.names = FALSE)
GO_T6PKSDMSO_Low <- GO

GO <- gopher(T6NPDMSO_Low,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana")
GO$go$name
GO$kegg$id
GO$pfam$name
write.table(data.frame(GO$go$id, GO$go$padj), file = "GO_T6NPDMSO_Low.txt", sep = "\t",
            row.names = FALSE)
GO_T6NPDMSO_Low <- GO

#' ## T6 AZD
#' ### Induced genes
GO <- gopher(T6NPSAZD_High,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana")
GO$go$name
GO$kegg$id
GO$pfam$name
write.table(data.frame(GO$go$id, GO$go$padj), file = "GO_T6NPSAZD_High.txt", sep = "\t",
            row.names = FALSE)
GO_T6NPSAZD_High <- GO

GO <- gopher(T6NSAZD_High,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana")
GO$go$name
GO$kegg$id
GO$pfam$name
write.table(data.frame(GO$go$id, GO$go$padj), file = "GO_T6NSAZD_High.txt", sep = "\t",
            row.names = FALSE)
GO_T6NSAZD_High <- GO

GO <- gopher(T6PKSAZD_High,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana")
GO$go$name
GO$kegg$id
GO$pfam$name
write.table(data.frame(GO$go$id, GO$go$padj), file = "GO_T6PKSAZD_High.txt", sep = "\t",
            row.names = FALSE)
GO_T6PKSAZD_High <- GO

GO <- gopher(T6NPAZD_High,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana")
GO$go$name
GO$kegg$id
GO$pfam$name
write.table(data.frame(GO$go$id, GO$go$padj), file = "GO_T6NPAZD_High.txt", sep = "\t",
            row.names = FALSE)
GO_T6NPAZD_High <- GO

#' ### Repressed genes
GO <- gopher(T6NPSAZD_Low,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana")
GO$go$name
GO$kegg$id
GO$pfam$name
write.table(data.frame(GO$go$id, GO$go$padj), file = "GO_T6NPSAZD_Low.txt", sep = "\t",
            row.names = FALSE)
GO_T6NPSAZD_Low <- GO

GO <- gopher(T6NSAZD_Low,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana")
GO$go$name
GO$kegg$id
GO$pfam$name
write.table(data.frame(GO$go$id, GO$go$padj), file = "GO_T6NSAZD_Low.txt", sep = "\t",
            row.names = FALSE)
GO_T6NSAZD_Low <- GO

GO <- gopher(T6PKSAZD_Low,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana")
GO$go$name
GO$kegg$id
GO$pfam$name
write.table(data.frame(GO$go$id, GO$go$padj), file = "GO_T6PKSAZD_Low.txt", sep = "\t",
            row.names = FALSE)
GO_T6PKSAZD_Low <- GO

GO <- gopher(T6NPAZD_Low,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana")
GO$go$name
GO$kegg$id
GO$pfam$name
write.table(data.frame(GO$go$id, GO$go$padj), file = "GO_T6NPAZD_Low.txt", sep = "\t",
            row.names = FALSE)
GO_T6NPAZD_Low <- GO

#' ## T24 DMSO
#' ### Induced genes
GO <- gopher(T24NPSDMSO_High,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana")
GO$go$name
GO$kegg$id
GO$pfam$name
write.table(data.frame(GO$go$id, GO$go$padj), file = "GO_T24NPSDMSO_High.txt", sep = "\t",
            row.names = FALSE)
GO_T24NPSDMSO_High <- GO

GO <- gopher(T24NSDMSO_High,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana")
GO$go$name
GO$kegg$id
GO$pfam$name
write.table(data.frame(GO$go$id, GO$go$padj), file = "GO_T24NSDMSO_High.txt", sep = "\t",
            row.names = FALSE)
GO_T24NSDMSO_High <- GO

GO <- gopher(T24PKSDMSO_High,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana")
GO$go$name
GO$kegg$id
GO$pfam$name
write.table(data.frame(GO$go$id, GO$go$padj), file = "GO_T24PKSDMSO_High.txt", sep = "\t",
            row.names = FALSE)
GO_T24PKSDMSO_High <- GO

GO <- gopher(T24NPDMSO_High,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana")
GO$go$name
GO$kegg$id
GO$pfam$name
write.table(data.frame(GO$go$id, GO$go$padj), file = "GO_T24NPDMSO_High.txt", sep = "\t",
            row.names = FALSE)
GO_T24NPDMSO_High <- GO

#' ### Repressed genes
GO <- gopher(T24NPSDMSO_Low,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana")
GO$go$name
GO$kegg$id
GO$pfam$name
write.table(data.frame(GO$go$id, GO$go$padj), file = "GO_T24NPSDMSO_Low.txt", sep = "\t",
            row.names = FALSE)
GO_T24NPSDMSO_Low <- GO

GO <- gopher(T24NSDMSO_Low,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana")
GO$go$name
GO$kegg$id
GO$pfam$name
write.table(data.frame(GO$go$id, GO$go$padj), file = "GO_T24NSDMSO_Low.txt", sep = "\t",
            row.names = FALSE)
GO_T24NSDMSO_Low <- GO

GO <- gopher(T24PKSDMSO_Low,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana")
GO$go$name
GO$kegg$id
GO$pfam$name
write.table(data.frame(GO$go$id, GO$go$padj), file = "GO_T24PKSDMSO_Low.txt", sep = "\t",
            row.names = FALSE)
GO_T24PKSDMSO_Low <- GO

GO <- gopher(T24NPDMSO_Low,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana")
GO$go$name
GO$kegg$id
GO$pfam$name
write.table(data.frame(GO$go$id, GO$go$padj), file = "GO_T24NPDMSO_Low.txt", sep = "\t",
            row.names = FALSE)
GO_T24NPDMSO_Low <- GO

#' ## T24 AZD
#' ### Induced genes
GO <- gopher(T24NPSAZD_High,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana")
GO$go$name
GO$kegg$id
GO$pfam$name
write.table(data.frame(GO$go$id, GO$go$padj), file = "GO_T24NPSAZD_High.txt", sep = "\t",
            row.names = FALSE)
GO_T24NPSAZD_High <- GO

GO <- gopher(T24NSAZD_High,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana")
GO$go$name
GO$kegg$id
GO$pfam$name
write.table(data.frame(GO$go$id, GO$go$padj), file = "GO_T24NSAZD_High.txt", sep = "\t",
            row.names = FALSE)
GO_T24NSAZD_High <- GO

GO <- gopher(T24PKSAZD_High,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana")
GO$go$name
GO$kegg$id
GO$pfam$name
write.table(data.frame(GO$go$id, GO$go$padj), file = "GO_T24PKSAZD_High.txt", sep = "\t",
            row.names = FALSE)
GO_T24PKSAZD_High <- GO

GO <- gopher(T24NPAZD_High,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana")
GO$go$name
GO$kegg$id
GO$pfam$name
write.table(data.frame(GO$go$id, GO$go$padj), file = "GO_T24NPAZD_High.txt", sep = "\t",
            row.names = FALSE)
GO_T24NPAZD_High <- GO

#' ### Repressed genes
GO <- gopher(T24NPSAZD_Low,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana")
GO$go$name
GO$kegg$id
GO$pfam$name
write.table(data.frame(GO$go$id, GO$go$padj), file = "GO_T24NPSAZD_Low.txt", sep = "\t",
            row.names = FALSE)
GO_T24NPSAZD_Low <- GO

GO <- gopher(T24NSAZD_Low,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana")
GO$go$name
GO$kegg$id
GO$pfam$name
write.table(data.frame(GO$go$id, GO$go$padj), file = "GO_T24NSAZD_Low.txt", sep = "\t",
            row.names = FALSE)
GO_T24NSAZD_Low <- GO

GO <- gopher(T24PKSAZD_Low,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana")
GO$go$name
GO$kegg$id
GO$pfam$name
write.table(data.frame(GO$go$id, GO$go$padj), file = "GO_T24PKSAZD_Low.txt", sep = "\t",
            row.names = FALSE)
GO_T24PKSAZD_Low <- GO

GO <- gopher(T24NPAZD_Low,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana")
GO$go$name
GO$kegg$id
GO$pfam$name
write.table(data.frame(GO$go$id, GO$go$padj), file = "GO_T24NPAZD_Low.txt", sep = "\t",
            row.names = FALSE)
GO_T24NPAZD_Low <- GO





#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```

