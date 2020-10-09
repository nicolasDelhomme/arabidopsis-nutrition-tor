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
#' Different from Schurch et al., RNA, 2016 (LFC at 1 instead of 0,5)
lfc <- 1
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

#' # Differential expression based on the nutrition and treatment at T6
#' ## Filtration of samples based on timepoint
sel <- samples$Timepoint %in% c("T6")
suppressMessages(dds <- DESeqDataSetFromMatrix(
    countData = kg[,sel],
    colData = samples[sel,],
    design = ~ Nutrition * AZD))

#' ## Differential expression analysis
dds <- DESeq(dds)

#' ## Variance Stabilising Transformation
#' ### Perform a Variance Stabilizing Transformation for plotting
vst <- varianceStabilizingTransformation(dds,blind=FALSE)
vsd <- assay(vst)
vsd <- vsd - min(vsd)

#' ## Contrasts to NPS_DMSO
#' The contrast by default is the first one (not Intercept)
resultsNames(dds)

#' ### Nutrition effect of carbon starvation
#' #### Extraction of the results from the DESeq analysis
res <- results(dds,name = "Nutrition_NP_vs_NPS")

cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
SucEffect6hrs <- rownames(res[cutoffs,])
SucLow6hrs <- rownames(res[cutoff2,])
SucHigh6hrs <- rownames(res[cutoff1,])

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))

#' #### MA plot
DESeq2::plotMA(res)

#' #### Volcano plot 
volcanoPlot(res, lfc = 1)

#' #### Heatmap
hm <- heatmap.2(t(scale(t(vsd[cutoffs,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))

#' ### Nutrition effect of phosphorus starvation
#' #### Extraction of the results from the DESeq analysis
res <- results(dds,name = "Nutrition_NS_vs_NPS")

cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
PiEffect6hrs <- rownames(res[cutoffs,])
PiLow6hrs <- rownames(res[cutoff2,])
PiHigh6hrs <- rownames(res[cutoff1,])

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))

#' #### MA plot
DESeq2::plotMA(res)

#' #### Volcano plot 
volcanoPlot(res, lfc=1)

#' #### Heatmap
hm <- heatmap.2(t(scale(t(vsd[cutoffs,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))

#' ### Nutrition effect of nitrogen starvation
#' #### Extraction of the results from the DESeq analysis
res <- results(dds,name = "Nutrition_PKS_vs_NPS")

cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
NitEffect6hrs <- rownames(res[cutoffs,])
NitLow6hrs <- rownames(res[cutoff2,])
NitHigh6hrs <- rownames(res[cutoff1,])

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))

#' #### MA plot
DESeq2::plotMA(res)

#' #### Volcano plot 
volcanoPlot(res, lfc=1)

#' #### Heatmap
hm <- heatmap.2(t(scale(t(vsd[cutoffs,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))

#' ### Nutrition effect of AZD-8055
#' #### Extraction of the results from the DESeq analysis
res <- results(dds,name = "AZD_AZD_vs_DMSO")

cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
AzdEffect6hrs <- rownames(res[cutoffs,])
AzdLow6hrs <- rownames(res[cutoff2,])
AzdHigh6hrs <- rownames(res[cutoff1,])

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))

#' #### MA plot
DESeq2::plotMA(res)

#' #### Volcano plot 
volcanoPlot(res, lfc=1)

#' #### Heatmap
hm <- heatmap.2(t(scale(t(vsd[cutoffs,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))

#' # Differential expression based on the nutrition and treatment at T24
#' ## Filtration of samples based on timepoint
sel <- samples$Timepoint %in% c("T24")
suppressMessages(dds <- DESeqDataSetFromMatrix(
    countData = kg[,sel],
    colData = samples[sel,],
    design = ~ Nutrition * AZD))

#' ## Differential expression analysis
dds <- DESeq(dds)

#' ## Variance Stabilising Transformation
#' ### Perform a Variance Stabilizing Transformation for plotting
vst <- varianceStabilizingTransformation(dds,blind=FALSE)
vsd <- assay(vst)
vsd <- vsd - min(vsd)

#' ## Contrasts to NPS_DMSO
#' The contrast by default is the first one (not Intercept)
resultsNames(dds)

#' ### Nutrition effect of carbon starvation
#' #### Extraction of the results from the DESeq analysis
res <- results(dds,name = "Nutrition_NP_vs_NPS")

cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
SucEffect24hrs <- rownames(res[cutoffs,])
SucLow24hrs <- rownames(res[cutoff2,])
SucHigh24hrs <- rownames(res[cutoff1,])

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))

#' #### MA plot
DESeq2::plotMA(res)

#' #### Volcano plot 
volcanoPlot(res, lfc=1)

#' #### Heatmap
hm <- heatmap.2(t(scale(t(vsd[cutoffs,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))

#' ### Nutrition effect of phosphorus starvation
#' #### Extraction of the results from the DESeq analysis
res <- results(dds,name = "Nutrition_NS_vs_NPS")

cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
PiEffect24hrs <- rownames(res[cutoffs,])
PiLow24hrs <- rownames(res[cutoff2,])
PiHigh24hrs <- rownames(res[cutoff1,])

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))

#' #### MA plot
DESeq2::plotMA(res)

#' #### Volcano plot 
volcanoPlot(res, lfc = 1)

#' #### Heatmap
hm <- heatmap.2(t(scale(t(vsd[cutoffs,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))

#' ### Nutrition effect of nitrogen starvation
#' #### Extraction of the results from the DESeq analysis
res <- results(dds,name = "Nutrition_PKS_vs_NPS")

cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
NitEffect24hrs <- rownames(res[cutoffs,])
NitLow24hrs <- rownames(res[cutoff2,])
NitHigh24hrs <- rownames(res[cutoff1,])

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))

#' #### MA plot
DESeq2::plotMA(res)

#' #### Volcano plot 
volcanoPlot(res, lfc=1)

#' #### Heatmap
hm <- heatmap.2(t(scale(t(vsd[cutoffs,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))

#' ### Nutrition effect of AZD-8055
#' #### Extraction of the results from the DESeq analysis
res <- results(dds,name = "AZD_AZD_vs_DMSO")

cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
AzdEffect24hrs <- rownames(res[cutoffs,])
AzdLow24hrs <- rownames(res[cutoff2,])
AzdHigh24hrs <- rownames(res[cutoff1,])

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))

#' #### MA plot
DESeq2::plotMA(res)

#' #### Volcano plot 
volcanoPlot(res, lfc=1)

#' #### Heatmap
hm <- heatmap.2(t(scale(t(vsd[cutoffs,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))

#' # Comparisons of GOI lists (Venn diagrams)
#' ## GOI at T6
#' * All GOI
grid.draw(venn.diagram(list(AzdEffect6hrs, PiEffect6hrs, SucEffect6hrs, NitEffect6hrs),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("AZD","Phos. Starv.","Suc. Starv.", "Nit. Starv.")))

#' * Downregulated genes
grid.draw(venn.diagram(list(AzdLow6hrs, PiLow6hrs, SucLow6hrs, NitLow6hrs),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("AZD","Phos. Starv.","Suc. Starv.", "Nit. Starv.")))

#' * Upregulated genes
grid.draw(venn.diagram(list(AzdHigh6hrs, PiHigh6hrs, SucHigh6hrs, NitHigh6hrs),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("AZD","Phos. Starv.","Suc. Starv.", "Nit. Starv.")))

#' ## GOI at T24
#' * All GOI
grid.draw(venn.diagram(list(AzdEffect24hrs, PiEffect24hrs, SucEffect24hrs, NitEffect24hrs),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("AZD","Phos. Starv.","Suc. Starv.", "Nit. Starv.")))

#' * Downregulated genes
grid.draw(venn.diagram(list(AzdLow24hrs, PiLow24hrs, SucLow24hrs, NitLow24hrs),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("AZD","Phos. Starv.","Suc. Starv.", "Nit. Starv.")))

#' * Upregulated genes
grid.draw(venn.diagram(list(AzdHigh24hrs, PiHigh24hrs, SucHigh24hrs, NitHigh24hrs),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("AZD","Phos. Starv.","Suc. Starv.", "Nit. Starv.")))

#' ## Between timepoints
#' ### In response to AZD treatment
#' * Upregulated
grid.draw(venn.diagram(list(AzdHigh6hrs, AzdHigh24hrs),
                       filename=NULL,
                       col=pal[1:2],
                       category.names=c("6hrs","24hrs")))
#' * DOwnregulated
grid.draw(venn.diagram(list(AzdLow6hrs, AzdLow24hrs),
                       filename=NULL,
                       col=pal[1:2],
                       category.names=c("6hrs","24hrs")))

#' ### In response to Phosphate starvation
#' * Upregulated
grid.draw(venn.diagram(list(PiHigh6hrs, PiHigh24hrs),
                       filename=NULL,
                       col=pal[1:2],
                       category.names=c("6hrs","24hrs")))
#' * DOwnregulated
grid.draw(venn.diagram(list(PiLow6hrs, PiLow24hrs),
                       filename=NULL,
                       col=pal[1:2],
                       category.names=c("6hrs","24hrs")))

#' ### In response to Sucrose starvation
#' * Upregulated
grid.draw(venn.diagram(list(SucHigh6hrs, SucHigh24hrs),
                       filename=NULL,
                       col=pal[1:2],
                       category.names=c("6hrs","24hrs")))
#' * DOwnregulated
grid.draw(venn.diagram(list(SucLow6hrs, SucLow24hrs),
                       filename=NULL,
                       col=pal[1:2],
                       category.names=c("6hrs","24hrs")))

#' ### In response to Nitrogen starvation
#' * Upregulated
grid.draw(venn.diagram(list(NitHigh6hrs, NitHigh24hrs),
                       filename=NULL,
                       col=pal[1:2],
                       category.names=c("6hrs","24hrs")))
#' * DOwnregulated
grid.draw(venn.diagram(list(NitLow6hrs, NitLow24hrs),
                       filename=NULL,
                       col=pal[1:2],
                       category.names=c("6hrs","24hrs")))

#' # Export list of DEGs
write.table(AzdHigh24hrs, file = "AzdHigh24hrs.txt", sep = "\t", row.names = FALSE)
write.table(AzdHigh6hrs, file = "AzdHigh6hrs.txt", sep = "\t", row.names = FALSE)
write.table(PiHigh24hrs, file = "PiHigh24hrs.txt", sep = "\t", row.names = FALSE)
write.table(PiHigh6hrs, file = "PiHigh6hrs.txt", sep = "\t", row.names = FALSE)
write.table(SucHigh24hrs, file = "SucHigh24hrs.txt", sep = "\t", row.names = FALSE)
write.table(SucHigh6hrs, file = "SucHigh6hrs.txt", sep = "\t", row.names = FALSE)

write.table(AzdLow24hrs, file = "AzdLow24hrs.txt", sep = "\t", row.names = FALSE)
write.table(AzdLow6hrs, file = "AzdLow6hrs.txt", sep = "\t", row.names = FALSE)
write.table(PiLow24hrs, file = "PiLow24hrs.txt", sep = "\t", row.names = FALSE)
write.table(PiLow6hrs, file = "PiLow6hrs.txt", sep = "\t", row.names = FALSE)
write.table(SucLow24hrs, file = "SucLow24hrs.txt", sep = "\t", row.names = FALSE)
write.table(SucLow6hrs, file = "SucLow6hrs.txt", sep = "\t", row.names = FALSE)




#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```

