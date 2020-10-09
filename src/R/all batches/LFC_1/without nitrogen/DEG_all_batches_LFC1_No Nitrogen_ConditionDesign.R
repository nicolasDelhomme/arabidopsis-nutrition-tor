#' ---
#' title: "DEG analysis after resequencing of the 3 last batches, excluding Nitrogen, LFC=1, excluding T0, Analysis by conditions"
#' author: "Thomas Dobrenel and Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' * Set the working dir
setwd("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions")
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
    filter(! Nutrition == "PKS") %>%
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
samples$Conditions <- factor(paste(samples$Timepoint,
                                   samples$Nutrition,
                                   samples$AZD,sep="_"),
                             levels=c("T0_T0_0",
                                      "T6_NPS_DMSO", "T6_NPS_AZD", "T6_NS_DMSO", "T6_NS_AZD", "T6_NP_DMSO", "T6_NP_AZD",
                                      "T24_NPS_DMSO", "T24_NPS_AZD", "T24_NS_DMSO", "T24_NS_AZD", "T24_NP_DMSO", "T24_NP_AZD"))
samples$Batch <- factor(substr(samples$SciLifeID,1,8))
samples <- cbind(samples,
                 Id = gsub(".+-",
                           "",
                           samples$SampleName))
samples <- cbind(samples,
                 Tp_Id = factor(paste(samples$Timepoint,
                                      samples$Id,
                                      sep=".")))
write.csv(samples,"/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/analysis_Tom/samplelist.csv")


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
stopifnot(all(colnames(vsd) == samples$SciLifeID))
colnames(vsd) <- samples$Tp_Id
write.csv(vsd,"/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/analysis_Tom/library-size-normalized_variance-stabilized_data_nutrition.csv")


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


#' # Differential expression compared to T0
#' ## Filtration of samples based on timepoint
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

#' The contrast by default is the first one (not Intercept)
resultsNames(dds)

#' ### Effect of 6 hrs of treatment
res <- results(dds,name = "Conditions_T6_NPS_DMSO_vs_T0_T0_0")
data <- data.frame(rownames(res), res$log2FoldChange, res$padj)
write.table(data, file = "Conditions_T6_NPS_DMSO_vs_T0_T0_0.csv", sep = ";", row.names = FALSE, dec=",")
cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))


res <- results(dds,name = "Conditions_T6_NPS_AZD_vs_T0_T0_0")
data <- data.frame(rownames(res), res$log2FoldChange, res$padj)
write.table(data, file = "Conditions_T6_NPS_AZD_vs_T0_T0_0.csv", sep = ";", row.names = FALSE, dec=",")
cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))


res <- results(dds,name = "Conditions_T6_NS_DMSO_vs_T0_T0_0")
data <- data.frame(rownames(res), res$log2FoldChange, res$padj)
write.table(data, file = "Conditions_T6_NS_DMSO_vs_T0_T0_0.csv", sep = ";", row.names = FALSE, dec=",")
cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))


res <- results(dds,name = "Conditions_T6_NS_AZD_vs_T0_T0_0")
data <- data.frame(rownames(res), res$log2FoldChange, res$padj)
write.table(data, file = "Conditions_T6_NS_AZD_vs_T0_T0_0.csv", sep = ";", row.names = FALSE, dec=",")
cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))


res <- results(dds,name = "Conditions_T6_NP_DMSO_vs_T0_T0_0")
data <- data.frame(rownames(res), res$log2FoldChange, res$padj)
write.table(data, file = "Conditions_T6_NP_DMSO_vs_T0_T0_0.csv", sep = ";", row.names = FALSE, dec=",")
cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))


res <- results(dds,name = "Conditions_T6_NP_AZD_vs_T0_T0_0")
data <- data.frame(rownames(res), res$log2FoldChange, res$padj)
write.table(data, file = "Conditions_T6_NP_AZD_vs_T0_T0_0.csv", sep = ";", row.names = FALSE, dec=",")
cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))




#' ### Effect of 24 hrs of treatment
res <- results(dds,name = "Conditions_T24_NPS_DMSO_vs_T0_T0_0")
data <- data.frame(rownames(res), res$log2FoldChange, res$padj)
write.table(data, file = "Conditions_T24_NPS_DMSO_vs_T0_T0_0.csv", sep = ";", row.names = FALSE, dec=",")
cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))


res <- results(dds,name = "Conditions_T24_NPS_AZD_vs_T0_T0_0")
data <- data.frame(rownames(res), res$log2FoldChange, res$padj)
write.table(data, file = "Conditions_T24_NPS_AZD_vs_T0_T0_0.csv", sep = ";", row.names = FALSE, dec=",")
cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))


res <- results(dds,name = "Conditions_T24_NS_DMSO_vs_T0_T0_0")
data <- data.frame(rownames(res), res$log2FoldChange, res$padj)
write.table(data, file = "Conditions_T24_NS_DMSO_vs_T0_T0_0.csv", sep = ";", row.names = FALSE, dec=",")
cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))


res <- results(dds,name = "Conditions_T24_NS_AZD_vs_T0_T0_0")
data <- data.frame(rownames(res), res$log2FoldChange, res$padj)
write.table(data, file = "Conditions_T24_NS_AZD_vs_T0_T0_0.csv", sep = ";", row.names = FALSE, dec=",")
cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))


res <- results(dds,name = "Conditions_T24_NP_DMSO_vs_T0_T0_0")
data <- data.frame(rownames(res), res$log2FoldChange, res$padj)
write.table(data, file = "Conditions_T24_NP_DMSO_vs_T0_T0_0.csv", sep = ";", row.names = FALSE, dec=",")
cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))


res <- results(dds,name = "Conditions_T24_NP_AZD_vs_T0_T0_0")
data <- data.frame(rownames(res), res$log2FoldChange, res$padj)
write.table(data, file = "Conditions_T24_NP_AZD_vs_T0_T0_0.csv", sep = ";", row.names = FALSE, dec=",")
cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))









#' # Differential expression based on the nutrition and treatment at T6
#' ## Filtration of samples based on timepoint
samples %<>% mutate(Conditions,Conditions=relevel(Conditions,"T6_NPS_DMSO"))
sel <- samples$Timepoint %in% c("T6")
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

#' ## Contrasts to NPS_DMSO
#' The contrast by default is the first one (not Intercept)
resultsNames(dds)

#' ### Nutrition effect of carbon starvation
#' #### Extraction of the results from the DESeq analysis
res <- results(dds,name = "Conditions_T6_NP_DMSO_vs_T6_NPS_DMSO")
data <- data.frame(rownames(res), res$log2FoldChange, res$padj)
write.table(data, file = "Conditions_T6_NP_DMSO_vs_T6_NPS_DMSO.csv", sep = ";", row.names = FALSE, dec=",")
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
res <- results(dds,name = "Conditions_T6_NS_DMSO_vs_T6_NPS_DMSO")
data <- data.frame(rownames(res), res$log2FoldChange, res$padj)
write.table(data, file = "Conditions_T6_NS_DMSO_vs_T6_NPS_DMSO.csv", sep = ";", row.names = FALSE, dec=",")
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


#' ### Effect of AZD-8055
#' #### Extraction of the results from the DESeq analysis
res <- results(dds,name = "Conditions_T6_NPS_AZD_vs_T6_NPS_DMSO")
data <- data.frame(rownames(res), res$log2FoldChange, res$padj)
write.table(data, file = "Conditions_T6_NPS_AZD_vs_T6_NPS_DMSO.csv", sep = ";", row.names = FALSE, dec=",")
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

#' ### Combined effect of AZD treatment and sugar starvation
#' #### Extraction of the results from the DESeq analysis
res <- results(dds,name = "Conditions_T6_NP_AZD_vs_T6_NPS_DMSO")
data <- data.frame(rownames(res), res$log2FoldChange, res$padj)
write.table(data, file = "Conditions_T6_NP_AZD_vs_T6_NPS_DMSO.csv", sep = ";", row.names = FALSE, dec=",")
cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
AzdNP6hrs <- rownames(res[cutoffs,])
AzdNPLow6hrs <- rownames(res[cutoff2,])
AzdNPHigh6hrs <- rownames(res[cutoff1,])

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

#' ### Combined effect of AZD treatment and phosphorus starvation
#' #### Extraction of the results from the DESeq analysis
res <- results(dds,name = "Conditions_T6_NS_AZD_vs_T6_NPS_DMSO")
data <- data.frame(rownames(res), res$log2FoldChange, res$padj)
write.table(data, file = "Conditions_T6_NS_AZD_vs_T6_NPS_DMSO.csv", sep = ";", row.names = FALSE, dec=",")
cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
AzdNS6hrs <- rownames(res[cutoffs,])
AzdNSLow6hrs <- rownames(res[cutoff2,])
AzdNSHigh6hrs <- rownames(res[cutoff1,])

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
samples %<>% mutate(Conditions,Conditions=relevel(Conditions,"T24_NPS_DMSO"))
sel <- samples$Timepoint %in% c("T24")
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

#' ## Contrasts to NPS_DMSO
#' The contrast by default is the first one (not Intercept)
resultsNames(dds)

#' ### Nutrition effect of carbon starvation
#' #### Extraction of the results from the DESeq analysis
res <- results(dds,name = "Conditions_T24_NP_DMSO_vs_T24_NPS_DMSO")
data <- data.frame(rownames(res), res$log2FoldChange, res$padj)
write.table(data, file = "Conditions_T24_NP_DMSO_vs_T24_NPS_DMSO.csv", sep = ";", row.names = FALSE, dec=",")
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
res <- results(dds,name = "Conditions_T24_NS_DMSO_vs_T24_NPS_DMSO")
data <- data.frame(rownames(res), res$log2FoldChange, res$padj)
write.table(data, file = "Conditions_T24_NS_DMSO_vs_T24_NPS_DMSO.csv", sep = ";", row.names = FALSE, dec=",")
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

#' ### Effect of AZD-8055
#' #### Extraction of the results from the DESeq analysis
res <- results(dds,name = "Conditions_T24_NPS_AZD_vs_T24_NPS_DMSO")
data <- data.frame(rownames(res), res$log2FoldChange, res$padj)
write.table(data, file = "Conditions_T24_NPS_AZD_vs_T24_NPS_DMSO.csv", sep = ";", row.names = FALSE, dec=",")
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

#' ### Combined effect of AZD treatment and sugar starvation
#' #### Extraction of the results from the DESeq analysis
res <- results(dds,name = "Conditions_T24_NP_AZD_vs_T24_NPS_DMSO")
data <- data.frame(rownames(res), res$log2FoldChange, res$padj)
write.table(data, file = "Conditions_T24_NP_AZD_vs_T24_NPS_DMSO.csv", sep = ";", row.names = FALSE, dec=",")
cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
AzdNP24hrs <- rownames(res[cutoffs,])
AzdNPLow24hrs <- rownames(res[cutoff2,])
AzdNPHigh24hrs <- rownames(res[cutoff1,])

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

#' ### Combined effect of AZD treatment and phosphorus starvation
#' #### Extraction of the results from the DESeq analysis
res <- results(dds,name = "Conditions_T24_NS_AZD_vs_T24_NPS_DMSO")
data <- data.frame(rownames(res), res$log2FoldChange, res$padj)
write.table(data, file = "Conditions_T24_NS_AZD_vs_T24_NPS_DMSO.csv", sep = ";", row.names = FALSE, dec=",")
cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
AzdNS24hrs <- rownames(res[cutoffs,])
AzdNSLow24hrs <- rownames(res[cutoff2,])
AzdNSHigh24hrs <- rownames(res[cutoff1,])

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
#' ### Comparison of AZD treatment with starvations
#' * All GOI
grid.draw(venn.diagram(list(AzdEffect6hrs, PiEffect6hrs, SucEffect6hrs),
                       filename=NULL,
                       col=pal[1:3],
                       category.names=c("AZD","Phos. Starv.","Suc. Starv.")))

#' * Downregulated genes
grid.draw(venn.diagram(list(AzdLow6hrs, PiLow6hrs, SucLow6hrs),
                       filename=NULL,
                       col=pal[1:3],
                       category.names=c("AZD","Phos. Starv.","Suc. Starv.")))

#' * Upregulated genes
grid.draw(venn.diagram(list(AzdHigh6hrs, PiHigh6hrs, SucHigh6hrs),
                       filename=NULL,
                       col=pal[1:3],
                       category.names=c("AZD","Phos. Starv.","Suc. Starv.")))

#' ### Interaction between AZD treatment and starvations
#' * Induced genes for sucrose starvation
grid.draw(venn.diagram(list(AzdHigh6hrs, AzdNPHigh6hrs, SucHigh6hrs),
                       filename=NULL,
                       col=pal[1:3],
                       category.names=c("AZD","Suc. Starv. + AZD","Suc. Starv.")))

#' * Repressed genes for sucrose starvation
grid.draw(venn.diagram(list(AzdLow6hrs, AzdNPLow6hrs, SucLow6hrs),
                       filename=NULL,
                       col=pal[1:3],
                       category.names=c("AZD","Suc. Starv. + AZD","Suc. Starv.")))

#' * Induced genes for phosphorus starvation
grid.draw(venn.diagram(list(AzdHigh6hrs, AzdNSHigh6hrs, PiHigh6hrs),
                       filename=NULL,
                       col=pal[1:3],
                       category.names=c("AZD","Pi Starv. + AZD","Pi Starv.")))

#' * Repressed genes for phosphorus starvation
grid.draw(venn.diagram(list(AzdLow6hrs, AzdNSLow6hrs, PiLow6hrs),
                       filename=NULL,
                       col=pal[1:3],
                       category.names=c("AZD","Pi Starv. + AZD","Pi Starv.")))

#' ## GOI at T24
#' #' ### Comparison of AZD treatment with starvations
#' * All GOI
grid.draw(venn.diagram(list(AzdEffect24hrs, PiEffect24hrs, SucEffect24hrs),
                       filename=NULL,
                       col=pal[1:3],
                       category.names=c("AZD","Phos. Starv.","Suc. Starv.")))

#' * Downregulated genes
grid.draw(venn.diagram(list(AzdLow24hrs, PiLow24hrs, SucLow24hrs),
                       filename=NULL,
                       col=pal[1:3],
                       category.names=c("AZD","Phos. Starv.","Suc. Starv.")))

#' * Upregulated genes
grid.draw(venn.diagram(list(AzdHigh24hrs, PiHigh24hrs, SucHigh24hrs),
                       filename=NULL,
                       col=pal[1:3],
                       category.names=c("AZD","Phos. Starv.","Suc. Starv.")))

#' ### Interaction between AZD treatment and starvations
#' * Induced genes for sucrose starvation
grid.draw(venn.diagram(list(AzdHigh24hrs, AzdNPHigh24hrs, SucHigh24hrs),
                       filename=NULL,
                       col=pal[1:3],
                       category.names=c("AZD","Suc. Starv. + AZD","Suc. Starv.")))

#' * Repressed genes for sucrose starvation
grid.draw(venn.diagram(list(AzdLow24hrs, AzdNPLow24hrs, SucLow24hrs),
                       filename=NULL,
                       col=pal[1:3],
                       category.names=c("AZD","Suc. Starv. + AZD","Suc. Starv.")))

#' * Induced genes for phosphorus starvation
grid.draw(venn.diagram(list(AzdHigh24hrs, AzdNSHigh24hrs, PiHigh24hrs),
                       filename=NULL,
                       col=pal[1:3],
                       category.names=c("AZD","Pi Starv. + AZD","Pi Starv.")))

#' * Repressed genes for phosphorus starvation
grid.draw(venn.diagram(list(AzdLow24hrs, AzdNSLow24hrs, PiLow24hrs),
                       filename=NULL,
                       col=pal[1:3],
                       category.names=c("AZD","Pi Starv. + AZD","Pi Starv.")))

#' ## Between timepoints
#' ### In response to AZD treatment
#' * Upregulated
grid.draw(venn.diagram(list(AzdHigh6hrs, AzdHigh24hrs),
                       filename=NULL,
                       col=pal[1:2],
                       category.names=c("6hrs","24hrs")))
#' * Downregulated
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
#' * Downregulated
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
#' * Downregulated
grid.draw(venn.diagram(list(SucLow6hrs, SucLow24hrs),
                       filename=NULL,
                       col=pal[1:2],
                       category.names=c("6hrs","24hrs")))

#' # Export list of DEGs
write.table(AzdHigh24hrs, file = "AzdHigh24hrs.csv", sep=";", row.names = FALSE)
write.table(AzdHigh6hrs, file = "AzdHigh6hrs.csv", sep=";", row.names = FALSE)
write.table(PiHigh24hrs, file = "PiHigh24hrs.csv", sep=";", row.names = FALSE)
write.table(PiHigh6hrs, file = "PiHigh6hrs.csv", sep=";", row.names = FALSE)
write.table(SucHigh24hrs, file = "SucHigh24hrs.csv", sep=";", row.names = FALSE)
write.table(SucHigh6hrs, file = "SucHigh6hrs.csv", sep=";", row.names = FALSE)

write.table(AzdLow24hrs, file = "AzdLow24hrs.csv", sep=";", row.names = FALSE)
write.table(AzdLow6hrs, file = "AzdLow6hrs.csv", sep=";", row.names = FALSE)
write.table(PiLow24hrs, file = "PiLow24hrs.csv", sep=";", row.names = FALSE)
write.table(PiLow6hrs, file = "PiLow6hrs.csv", sep=";", row.names = FALSE)
write.table(SucLow24hrs, file = "SucLow24hrs.csv", sep=";", row.names = FALSE)
write.table(SucLow6hrs, file = "SucLow6hrs.csv", sep=";", row.names = FALSE)

write.table(AzdNPLow6hrs, file = "AzdNPLow6hrs.csv", sep=";", row.names = FALSE)
write.table(AzdNPHigh6hrs, file = "AzdNPHigh6hrs.csv", sep=";", row.names = FALSE)
write.table(AzdNPLow24hrs, file = "AzdNPLow24hrs.csv", sep=";", row.names = FALSE)
write.table(AzdNPHigh24hrs, file = "AzdNPHigh24hrs.csv", sep=";", row.names = FALSE)

write.table(AzdNSLow6hrs, file = "AzdNSLow6hrs.csv", sep=";", row.names = FALSE)
write.table(AzdNSHigh6hrs, file = "AzdNSHigh6hrs.csv", sep=";", row.names = FALSE)
write.table(AzdNSLow24hrs, file = "AzdNSLow24hrs.csv", sep=";", row.names = FALSE)
write.table(AzdNSHigh24hrs, file = "AzdNSHigh24hrs.csv", sep=";", row.names = FALSE)

#' # GO enrichment analysis
#' ## For AZD treatment
#' * Induced genes after 6hrs
GO <- gopher(AzdHigh6hrs,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana", alpha=2)
write.table(GO$go, file = "GO_AzdHigh6hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$kegg, file = "KEGG_AzdHigh6hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$pfam, file = "PFAM_AzdHigh6hrs.csv", sep = ";", row.names = FALSE)

#' * Induced genes after 24hrs
GO <- gopher(AzdHigh24hrs,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana",alpha=2)
write.table(GO$go, file = "GO_AzdHigh24hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$kegg, file = "KEGG_AzdHigh24hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$pfam, file = "PFAM_AzdHigh24hrs.csv", sep = ";", row.names = FALSE)


#' * Repressed genes after 6hrs
GO <- gopher(AzdLow6hrs,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana",alpha=2)
write.table(GO$go, file = "GO_AzdLow6hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$kegg, file = "KEGG_AzdLow6hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$pfam, file = "PFAM_AzdLow6hrs.csv", sep = ";", row.names = FALSE)

#' * Repressed genes after 24hrs
GO <- gopher(AzdLow24hrs,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana",alpha=2)
write.table(GO$go, file = "GO_AzdLow24hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$kegg, file = "KEGG_AzdLow24hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$pfam, file = "PFAM_AzdLow24hrs.csv", sep = ";", row.names = FALSE)

#' ## For Sucrose starvation
#' * Induced genes after 6hrs
GO <- gopher(SucHigh6hrs,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana", alpha=2)
write.table(GO$go, file = "GO_SucHigh6hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$kegg, file = "KEGG_SucHigh6hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$pfam, file = "PFAM_SucHigh6hrs.csv", sep = ";", row.names = FALSE)

#' * Induced genes after 24hrs
GO <- gopher(SucHigh24hrs,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana", alpha=2)
write.table(GO$go, file = "GO_SucHigh24hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$kegg, file = "KEGG_SucHigh24hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$pfam, file = "PFAM_SucHigh24hrs.csv", sep = ";", row.names = FALSE)


#' * Repressed genes after 6hrs
GO <- gopher(SucLow6hrs,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana", alpha=2)
write.table(GO$go, file = "GO_SucLow6hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$kegg, file = "KEGG_SucLow6hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$pfam, file = "PFAM_SucLow6hrs.csv", sep = ";", row.names = FALSE)

#' * Repressed genes after 24hrs
GO <- gopher(SucLow24hrs,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana", alpha=2)
write.table(GO$go, file = "GO_SucLow24hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$kegg, file = "KEGG_SucLow24hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$pfam, file = "PFAM_SucLow24hrs.csv", sep = ";", row.names = FALSE)

#' ## For Phosphorus starvation
#' * Induced genes after 6hrs
GO <- gopher(PiHigh6hrs,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana", alpha=2)
write.table(GO$go, file = "GO_PiHigh6hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$kegg, file = "KEGG_PiHigh6hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$pfam, file = "PFAM_PiHigh6hrs.csv", sep = ";", row.names = FALSE)

#' * Induced genes after 24hrs
GO <- gopher(PiHigh24hrs,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana", alpha=2)
write.table(GO$go, file = "GO_PiHigh24hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$kegg, file = "KEGG_PiHigh24hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$pfam, file = "PFAM_PiHigh24hrs.csv", sep = ";", row.names = FALSE)

#' * Repressed genes after 6hrs
GO <- gopher(PiLow6hrs,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana", alpha=2)
write.table(GO$go, file = "GO_PiLow6hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$kegg, file = "KEGG_PiLow6hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$pfam, file = "PFAM_PiLow6hrs.csv", sep = ";", row.names = FALSE)

#' * Repressed genes after 24hrs
GO <- gopher(PiLow24hrs,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana", alpha=2)
write.table(GO$go, file = "GO_PiLow24hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$kegg, file = "KEGG_PiLow24hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$pfam, file = "PFAM_PiLow24hrs.csv", sep = ";", row.names = FALSE)



#' ## For the interaction between AZD treatment and phosphorus starvation
#' * Induced genes after 6hrs
GO <- gopher(AzdNSHigh6hrs,background=rownames(vsd),url="athaliana", alpha=2)
write.table(GO$go, file = "GO_AzdNSHigh6hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$kegg, file = "KEGG_AzdNSHigh6hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$pfam, file = "PFAM_AzdNSHigh6hrs.csv", sep = ";", row.names = FALSE)

#' * Induced genes after 24hrs
GO <- gopher(AzdNSHigh24hrs,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana", alpha=2)
write.table(GO$go, file = "GO_AzdNSHigh24hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$kegg, file = "KEGG_AzdNSHigh24hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$pfam, file = "PFAM_AzdNSHigh24hrs.csv", sep = ";", row.names = FALSE)


#' * Repressed genes after 6hrs
GO <- gopher(AzdNSLow6hrs,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana", alpha=2)
write.table(GO$go, file = "GO_AzdNSLow6hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$kegg, file = "KEGG_AzdNSLow6hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$pfam, file = "PFAM_AzdNSLow6hrs.csv", sep = ";", row.names = FALSE)

#' * Repressed genes after 24hrs
GO <- gopher(AzdNSLow24hrs,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana", alpha=2)
write.table(GO$go, file = "GO_AzdNSLow24hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$kegg, file = "KEGG_AzdNSLow24hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$pfam, file = "PFAM_AzdNSLow24hrs.csv", sep = ";", row.names = FALSE)


#' ## For the interaction between AZD treatment and sugar starvation
#' * Induced genes after 6hrs
GO <- gopher(AzdNPHigh6hrs,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana", alpha=2)
write.table(GO$go, file = "GO_AzdNPHigh6hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$kegg, file = "KEGG_AzdNPHigh6hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$pfam, file = "PFAM_AzdNPHigh6hrs.csv", sep = ";", row.names = FALSE)

#' * Induced genes after 24hrs
GO <- gopher(AzdNPHigh24hrs,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana", alpha=2)
write.table(GO$go, file = "GO_AzdNPHigh24hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$kegg, file = "KEGG_AzdNPHigh24hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$pfam, file = "PFAM_AzdNPHigh24hrs.csv", sep = ";", row.names = FALSE)

#' * Repressed genes after 6hrs
GO <- gopher(AzdNPLow6hrs,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana", alpha=2)
write.table(GO$go, file = "GO_AzdNPLow6hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$kegg, file = "KEGG_AzdNPLow6hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$pfam, file = "PFAM_AzdNPLow6hrs.csv", sep = ";", row.names = FALSE)

#' * Repressed genes after 24hrs
GO <- gopher(AzdNPLow24hrs,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana", alpha=2)
write.table(GO$go, file = "GO_AzdNPLow24hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$kegg, file = "KEGG_AzdNPLow24hrs.csv", sep = ";", row.names = FALSE) 
write.table(GO$pfam, file = "PFAM_AzdNPLow24hrs.csv", sep = ";", row.names = FALSE)






#' # Differential expression of +/-AZD in -Pi samples at 24hrs
#' ## Filtration of samples based on timepoint
samples %<>% mutate(Conditions,Conditions=relevel(Conditions,"T24_NS_DMSO"))
sel <- samples$Timepoint %in% c("T24")
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

#' ## Contrasts to NPS_DMSO
#' The contrast by default is the first one (not Intercept)
resultsNames(dds)

#' ### Nutrition effect of carbon starvation
#' #### Extraction of the results from the DESeq analysis
res <- results(dds,name = "Conditions_T24_NS_AZD_vs_T24_NS_DMSO")
data <- data.frame(rownames(res), res$log2FoldChange, res$padj)
write.table(data, file = "Conditions_T24_NS_AZD_vs_T24_NS_DMSO.csv", sep = ";", row.names = FALSE, dec=",")
cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
AzdInPiEffect24hrs <- rownames(res[cutoffs,])
AzdInPiLow24hrs <- rownames(res[cutoff2,])
AzdInPiHigh24hrs <- rownames(res[cutoff1,])

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))
message(sprintf("There are %s genes that are induced",sum(cutoff1)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2)))
write.table(AzdInPiLow24hrs, file = "AzdInPiLow24hrs.csv", sep = ";", row.names = FALSE)
write.table(AzdInPiHigh24hrs, file = "AzdInPiHigh24hrs.csv", sep = ";", row.names = FALSE)

#' #### MA plot
DESeq2::plotMA(res)

#' #### Volcano plot 
volcanoPlot(res, lfc = 0.5)

#' #### Heatmap
hm <- heatmap.2(t(scale(t(vsd[cutoffs,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))

#' ## GOpher analysis
#' * Induced genes
GO <- gopher(AzdInPiHigh24hrs,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana", alpha=2)
write.table(GO$go, file = "GO_AzdInPiHigh24hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$kegg, file = "KEGG_AzdInPiHigh24hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$pfam, file = "PFAM_AzdInPiHigh24hrs.csv", sep = ";", row.names = FALSE)
#View(GO$go)
#View(GO$kegg)

#' * Repressed genes
GO <- gopher(AzdInPiLow24hrs,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana", alpha=2)
write.table(GO$go, file = "GO_AzdInPiLow24hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$kegg, file = "KEGG_AzdInPiLow24hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$pfam, file = "PFAM_AzdInPiLow24hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$pfam, file = "PFAM_AzdInPiHigh24hrs.csv", sep = ";", row.names = FALSE)
#View(GO$go)
#View(GO$kegg)


#' # Expression level analysis of specific gene lists
levels(samples$Conditions)

sel_T0 <- samples$Timepoint %in% c("T0","T6","T24")

suppressMessages(dds_T0 <- DESeqDataSetFromMatrix(
    countData = kg[,sel_T0],
    colData = samples[sel_T0,],
    design = ~ Conditions))

#' ## Differential expression analysis
dds_T0 <- DESeq(dds_T0)

#' ## Variance Stabilising Transformation
#' ### Perform a Variance Stabilizing Transformation for plotting
vst_T0 <- varianceStabilizingTransformation(dds_T0,blind=FALSE)
vsd_T0 <- assay(vst_T0)
vsd_T0 <- vsd_T0 - min(vsd_T0)

#' ## Export the complete list of averages
sel1 <- rownames(vsd_T0)
expr <- vsd_T0[match(sel1,rownames(vsd_T0)),]
Avg <- sapply(split.data.frame(t(expr),f = droplevels(samples$Conditions[sel_T0])),colMeans)
write.table(Avg, file = "Averaged_Normalized_Counts.csv", sep = ";", row.names = FALSE, dec=",")




#' ## TOR complex members
#' ### Preparation of the complex member list
sel1 <- c("AT1G50030","AT3G18140","AT2G22040","AT3G08850","AT5G01770") 


#' ### Preparation of the expression table for the shortlisted genes
expr <- vsd_T0[match(sel1,rownames(vsd_T0)),]
#expr <- vsd_T0[match(rownames(vsd_T0),sel1),]
#expr <- na.omit(expr)
AvgTOR <- sapply(split.data.frame(t(expr),f = droplevels(samples$Conditions[sel_T0])),colMeans)
SDTOR <- sapply(split.data.frame(t(expr),f = droplevels(samples$Conditions[sel_T0])),colSds)

#' ### Graphical representations
#' #### Combined barplots
barplot(AvgTOR,beside=T, cex.names=0.7,legend.text = sel1,ylim=c(0,4))
barcenters <- barplot(AvgTOR,beside=T, cex.names=0.7,legend.text = sel1,ylim=c(0,4))
arrows(x0=barcenters, x1=barcenters, y0=AvgTOR-SDTOR, y1=AvgTOR+SDTOR,
       lwd=1.5, angle=90, length=0.05,code=3)

#' #### Separated barplots
for (i in 1:length(AvgTOR[,1])){
    barcenters <- barplot(AvgTOR[i,],beside=T, cex.names=0.7,legend.text = sel1[i],ylim=c(0,4))
    arrows(x0=barcenters, x1=barcenters, y0=AvgTOR[i,]-SDTOR[i,], y1=AvgTOR[i,]+SDTOR[i,],
           lwd=1.5, angle=90, length=0.2,code=3)
}


#' #### Individual graphs
for (i in 1:length(AvgTOR[,1]))
{
    plot(AvgTOR[i,],ylim=c(min(AvgTOR[i,]-SDTOR[i,]),max(AvgTOR[i,]+SDTOR[i,])),
         xaxt="n",xlab="",ylab=sprintf("gene = %s",rownames(AvgTOR)[i]))
    axis(side = 1,at=1:13, colnames(AvgTOR),cex.axis=0.7,las=2)
    arrows(x0=1:length(AvgTOR[i,]), x1=1:length(AvgTOR[i,]),
           y0=AvgTOR[i,]-SDTOR[i,], y1=AvgTOR[i,]+SDTOR[i,],
           lwd=1.5, angle=90, length=0.05,code=3)
}

#' #### Heatmap
AvgTOR[AvgTOR == 0] <- 0.000001
AvgTOR_norm <- log2(AvgTOR[,] / AvgTOR[,1])
AvgTOR_norm[AvgTOR_norm < -0.5] <- -0.5; AvgTOR_norm[AvgTOR_norm > 0.5] <- 0.5


heatmap.2(AvgTOR_norm,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

#' ## Genes involved in the cell cycle
#' ### Preparation of a list of genes of the cell cycle
sel1 <- read.csv("~/arabidopsis-nutrition-TOR/seidr/cell_cycle_genes.csv", sep=";")[,1]

#' ## Preparation of the expression table for the shortlisted genes
expr <- vsd_T0[match(sel1,rownames(vsd_T0)),]
#expr <- vsd_T0[match(rownames(vsd_T0),sel1),]
#expr <- na.omit(expr)
AvgTOR <- sapply(split.data.frame(t(expr),f = droplevels(samples$Conditions[sel_T0])),colMeans)
SDTOR <- sapply(split.data.frame(t(expr),f = droplevels(samples$Conditions[sel_T0])),colSds)

#' ## Graphical representations
#' ### Combined barplots
barplot(AvgTOR,beside=T, cex.names=0.7,legend.text = sel1,ylim=c(0,4))
barcenters <- barplot(AvgTOR,beside=T, cex.names=0.7,legend.text = sel1,ylim=c(0,4))
arrows(x0=barcenters, x1=barcenters, y0=AvgTOR-SDTOR, y1=AvgTOR+SDTOR,
       lwd=1.5, angle=90, length=0.05,code=3)

#' ### Individual graphs without colors
for (i in 1:length(AvgTOR[,1]))
{
    plot(AvgTOR[i,],ylim=c(min(AvgTOR[i,]-SDTOR[i,]),max(AvgTOR[i,]+SDTOR[i,])),
         xaxt="n",xlab="",ylab=sprintf("gene = %s",rownames(AvgTOR)[i]))
    axis(side = 1,at=1:13, colnames(AvgTOR),cex.axis=0.7,las=2)
    arrows(x0=1:length(AvgTOR[i,]), x1=1:length(AvgTOR[i,]),
           y0=AvgTOR[i,]-SDTOR[i,], y1=AvgTOR[i,]+SDTOR[i,],
           lwd=1.5, angle=90, length=0.05,code=3)
}

#' ### Individual graphs with colors
a <- matrix(NA, nrow=6, ncol=3)
colnames(a) <- c(0,6,24)
rownames(a) <- c("NPS_DMSO","NPS_AZD","NS_DMSO","NS_AZD","NP_DMSO","NP_AZD")
b <- matrix(NA, nrow=6, ncol=3)
colnames(b) <- c(0,6,24)
rownames(b) <- c("NPS_DMSO","NPS_AZD","NS_DMSO","NS_AZD","NP_DMSO","NP_AZD")


for (i in 1:length(AvgTOR[,1]))
    {
    a[1:6,1] <- AvgTOR[i,1]
    a[1:6,2] <- AvgTOR[i,2:7]
    a[1:6,3] <- AvgTOR[i,8:13]
    
    b[1:6,1] <- SDTOR[i,1]
    b[1:6,2] <- SDTOR[i,2:7]
    b[1:6,3] <- SDTOR[i,8:13]
    
    plot(colnames(a),a[1,], type="l",ylim=c(min(a)-max(b),max(a)+max(b)),ylab=rownames(AvgTOR)[i])
    lines(colnames(a),a[2,], type="l",ylim=c(min(a)-max(b),max(a)+max(b)),lty=2)
    lines(colnames(a),a[3,], type="l",ylim=c(min(a)-max(b),max(a)+max(b)),col="turquoise")
    lines(colnames(a),a[4,], type="l",ylim=c(min(a)-max(b),max(a)+max(b)),col="turquoise",lty=2)
    lines(colnames(a),a[5,], type="l",ylim=c(min(a)-max(b),max(a)+max(b)),col="hotpink3")
    lines(colnames(a),a[6,], type="l",ylim=c(min(a)-max(b),max(a)+max(b)),col="hotpink3",lty=2)
    
    arrows(x0=c(0,6,24),
           x1=c(0,6,24),
           y0=a[1,] - b[1,],
           y1=a[1,] + b[1,],
           lwd=1, angle=90, length=0.05, code=3)
    arrows(x0=c(0,6,24),
           x1=c(0,6,24),
           y0=a[2,] - b[2,],
           y1=a[2,] + b[2,],
           lwd=1, angle=90, length=0.05, code=3)
    arrows(x0=c(0,6,24),
           x1=c(0,6,24),
           y0=a[3,] - b[3,],
           y1=a[3,] + b[3,],
           lwd=1, angle=90, length=0.05, code=3, col="turquoise")
    arrows(x0=c(0,6,24),
           x1=c(0,6,24),
           y0=a[4,] - b[4,],
           y1=a[4,] + b[4,],
           lwd=1, angle=90, length=0.05, code=3, col="turquoise")
    arrows(x0=c(0,6,24),
           x1=c(0,6,24),
           y0=a[5,] - b[5,],
           y1=a[5,] + b[5,],
           lwd=1, angle=90, length=0.05, code=3, col="hotpink3")
    arrows(x0=c(0,6,24),
           x1=c(0,6,24),
           y0=a[6,] - b[6,],
           y1=a[6,] + b[6,],
           lwd=1, angle=90, length=0.05, code=3, col="hotpink3")
}







#' ### Heatmap
AvgTOR[AvgTOR == 0] <- 0.000001
AvgTOR_norm <- log2(AvgTOR[,] / AvgTOR[,1])

AvgTOR_norm[AvgTOR_norm < -3] <- -3; AvgTOR_norm[AvgTOR_norm > 3] <- 3

heatmap.2(AvgTOR_norm,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8),
          labRow = read.csv("~/arabidopsis-nutrition-TOR/seidr/cell_cycle_genes.csv", sep=";")[,2])

heatmap.2(AvgTOR_norm,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          labRow = read.csv("~/arabidopsis-nutrition-TOR/seidr/cell_cycle_genes.csv", sep=";")[,2],
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})



#' ## Preparation of other lists
#' * Preparation of the hexokinase list
sel1 <- c("AT1G05205","AT1G47840","AT1G47845","AT2G19860","AT4G29130")
#' * Preparation of the AtHXK1
sel1 <- c("AT1G47845","AT4G29130")



#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```

