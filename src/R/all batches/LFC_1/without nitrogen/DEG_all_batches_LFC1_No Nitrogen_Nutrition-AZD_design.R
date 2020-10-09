#' ---
#' title: "DEG analysis after resequencing of the 3 last batches, excluding Nitrogen, LFC=1, excluding T0, Combined analysis of nutrition and AZD"
#' author: "Thomas Dobrenel and Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' * Set the working dir
setwd("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Nutrition_AZD")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Nutrition_AZD")
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
samples %<>% mutate(Conditions,Conditions=relevel(Conditions,"T6_NPS_DMSO"))
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
# See output of "Condition design"

#' ### Nutrition effect of phosphorus starvation
# See output of "Condition design"

#' ### Effect of AZD-8055
# See output of "Condition design"


#' ### Differential response of AZD & sugar starvation from AZD + sugar starvation responses
#' #### Extraction of the results from the DESeq analysis
res <- results(dds,name = "NutritionNP.AZDAZD")
data <- data.frame(rownames(res), res$log2FoldChange, res$padj)
write.table(data, file = "NutritionNP.AZDAZD_6h.csv", sep = ";", row.names = FALSE, dec=",")
cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
Azd_NP6hrs <- rownames(res[cutoffs,])
Azd_NPLow6hrs <- rownames(res[cutoff2,])
Azd_NPHigh6hrs <- rownames(res[cutoff1,])

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

#' ### Differential response of AZD & phosphorus starvation from AZD + phosphorus starvation responses
#' #### Extraction of the results from the DESeq analysis
res <- results(dds,name = "NutritionNS.AZDAZD")
data <- data.frame(rownames(res), res$log2FoldChange, res$padj)
write.table(data, file = "NutritionNS.AZDAZD_6h.csv", sep = ";", row.names = FALSE, dec=",")
cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
Azd_NS6hrs <- rownames(res[cutoffs,])
Azd_NSLow6hrs <- rownames(res[cutoff2,])
Azd_NSHigh6hrs <- rownames(res[cutoff1,])

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
# See output of "Condition design"

#' ### Nutrition effect of phosphorus starvation
# See output of "Condition design"

#' ### Effect of AZD-8055
# See output of "Condition design"


#' ### Differential response of AZD & sugar starvation from AZD + sugar starvation responses
#' #### Extraction of the results from the DESeq analysis
res <- results(dds,name = "NutritionNP.AZDAZD")
data <- data.frame(rownames(res), res$log2FoldChange, res$padj)
write.table(data, file = "NutritionNP.AZDAZD_24h.csv", sep = ";", row.names = FALSE, dec=",")
cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
Azd_NP24hrs <- rownames(res[cutoffs,])
Azd_NPLow24hrs <- rownames(res[cutoff2,])
Azd_NPHigh24hrs <- rownames(res[cutoff1,])

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

#' ### Differential response of AZD & phosphorus starvation from AZD + phosphorus starvation responses
#' #### Extraction of the results from the DESeq analysis
res <- results(dds,name = "NutritionNS.AZDAZD")
data <- data.frame(rownames(res), res$log2FoldChange, res$padj)
write.table(data, file = "NutritionNS.AZDAZD_24h.csv", sep = ";", row.names = FALSE, dec=",")
cutoffs <- abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff1 <- res$log2FoldChange >= lfc & ! is.na(res$padj) & res$padj <= FDR
cutoff2 <- res$log2FoldChange <= -lfc & ! is.na(res$padj) & res$padj <= FDR
Azd_NS24hrs <- rownames(res[cutoffs,])
Azd_NSLow24hrs <- rownames(res[cutoff2,])
Azd_NSHigh24hrs <- rownames(res[cutoff1,])

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





#' # Export list of DEGs
write.table(Azd_NPLow6hrs, file = "Azd_NPLow6hrs.csv", sep = ";", row.names = FALSE)
write.table(Azd_NPHigh6hrs, file = "Azd_NPHigh6hrs.csv", sep = ";", row.names = FALSE)

write.table(Azd_NSLow6hrs, file = "Azd_NSLow6hrs.csv", sep = ";", row.names = FALSE)
write.table(Azd_NSHigh6hrs, file = "Azd_NSHigh6hrs.csv", sep = ";", row.names = FALSE)

write.table(Azd_NPLow24hrs, file = "Azd_NPLow24hrs.csv", sep = ";", row.names = FALSE)
write.table(Azd_NPHigh24hrs, file = "Azd_NPHigh24hrs.csv", sep = ";", row.names = FALSE)

write.table(Azd_NSLow24hrs, file = "Azd_NSLow24hrs.csv", sep = ";", row.names = FALSE)
write.table(Azd_NSHigh24hrs, file = "Azd_NSHigh24hrs.csv", sep = ";", row.names = FALSE)





#' # GO enrichment analysis
#' ## For the interaction before AZD and sugar starvation
#' * Induced genes after 6hrs
GO <- gopher(Azd_NPHigh6hrs,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana", alpha=2)
write.table(GO$go, file = "GO_Azd_NPHigh6hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$kegg, file = "KEGG_Azd_NPHigh6hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$pfam, file = "PFAM_Azd_NPHigh6hrs.csv", sep = ";", row.names = FALSE)

#' * Induced genes after 24hrs
GO <- gopher(Azd_NPHigh24hrs,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana",alpha=2)
write.table(GO$go, file = "GO_Azd_NPHigh24hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$kegg, file = "KEGG_Azd_NPHigh24hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$pfam, file = "PFAM_Azd_NPHigh24hrs.csv", sep = ";", row.names = FALSE)


#' * Repressed genes after 6hrs
GO <- gopher(Azd_NPLow6hrs,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana",alpha=2)
write.table(GO$go, file = "GO_Azd_NPLow6hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$kegg, file = "KEGG_Azd_NPLow6hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$pfam, file = "PFAM_Azd_NPLow6hrs.csv", sep = ";", row.names = FALSE)

#' * Repressed genes after 24hrs
GO <- gopher(Azd_NPLow24hrs,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana",alpha=2)
write.table(GO$go, file = "GO_Azd_NPLow24hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$kegg, file = "KEGG_Azd_NPLow24hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$pfam, file = "PFAM_Azd_NPLow24hrs.csv", sep = ";", row.names = FALSE)

#' ## For the interaction between AZD and phosphorus starvation
#' * Induced genes after 6hrs
GO <- gopher(Azd_NSHigh6hrs,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana", alpha=2)
write.table(GO$go, file = "GO_Azd_NSHigh6hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$kegg, file = "KEGG_Azd_NSHigh6hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$pfam, file = "PFAM_Azd_NSHigh6hrs.csv", sep = ";", row.names = FALSE)

#' * Induced genes after 24hrs
GO <- gopher(Azd_NSHigh24hrs,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana", alpha=2)
write.table(GO$go, file = "GO_Azd_NSHigh24hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$kegg, file = "KEGG_Azd_NSHigh24hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$pfam, file = "PFAM_Azd_NSHigh24hrs.csv", sep = ";", row.names = FALSE)


#' * Repressed genes after 6hrs
GO <- gopher(Azd_NSLow6hrs,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana", alpha=2)
write.table(GO$go, file = "GO_Azd_NSLow6hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$kegg, file = "KEGG_Azd_NSLow6hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$pfam, file = "PFAM_Azd_NSLow6hrs.csv", sep = ";", row.names = FALSE)

#' * Repressed genes after 24hrs
GO <- gopher(Azd_NSLow24hrs,background=rownames(vsd)[rowSums(vsd)>0],url="athaliana", alpha=2)
write.table(GO$go, file = "GO_Azd_NSLow24hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$kegg, file = "KEGG_Azd_NSLow24hrs.csv", sep = ";", row.names = FALSE)
write.table(GO$pfam, file = "PFAM_Azd_NSLow24hrs.csv", sep = ";", row.names = FALSE)



#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```

