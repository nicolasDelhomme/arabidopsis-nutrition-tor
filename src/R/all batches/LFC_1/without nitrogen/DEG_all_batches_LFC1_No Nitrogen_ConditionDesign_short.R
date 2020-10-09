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

legend("bottomleft",pch=19,
       col=pal[1:nlevels(conds)],
       legend=levels(conds))

legend("topleft",pch=c(15,16,17),
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

#' ## Export the complete list of averages and SDs
sel1 <- rownames(vsd_T0)
expr <- vsd_T0[match(sel1,rownames(vsd_T0)),]
Avg <- sapply(split.data.frame(t(expr),
                               f = droplevels(samples$Conditions[sel_T0])),
              colMeans)
SD <- sapply(split.data.frame(t(expr),f = droplevels(samples$Conditions[sel_T0])),colSds)
rownames(SD) <- rownames(Avg)
write.table(Avg, file = "Averaged_Normalized_Counts.csv", sep = ";", row.names = FALSE, dec=",")
write.table(rownames(Avg), file = "Averaged_Normalized_Counts_genes.csv", sep = ";", row.names = FALSE, dec=",")




#' ## TOR complex members
#' ### Preparation of the complex member list
sel1 <- c("AT1G50030","AT3G18140","AT2G22040","AT3G08850","AT5G01770") 


#' ### Preparation of the expression table for the shortlisted genes
AvgTOR <- Avg[match(sel1,rownames(Avg)),]
SDTOR <- SD[match(sel1,rownames(SD)),]

#' ### Graphical representations
#' #### Combined barplots
barcenters <- barplot(AvgTOR,beside=T, cex.names=0.7,legend.text = sel1,ylim=c(0,4))
arrows(x0=barcenters, x1=barcenters, y0=AvgTOR-SDTOR, y1=AvgTOR+SDTOR,
       lwd=1.5, angle=90, length=0.05,code=3)

#' #### Separated barplots
for (i in 1:length(AvgTOR[,1])){
    barcenters <- barplot(AvgTOR[i,],beside=T, cex.names=0.7,legend.text = sel1[i],ylim=c(0,4))
    arrows(x0=barcenters, x1=barcenters, y0=AvgTOR[i,]-SDTOR[i,], y1=AvgTOR[i,]+SDTOR[i,],
           lwd=1.5, angle=90, length=0.2,code=3)
}


#' #### Individual graphs without colors
for (i in 1:length(AvgTOR[,1]))
{
    plot(AvgTOR[i,],ylim=c(min(AvgTOR[i,]-SDTOR[i,]),max(AvgTOR[i,]+SDTOR[i,])),
         xaxt="n",xlab="",ylab=sprintf("gene = %s",rownames(AvgTOR)[i]))
    axis(side = 1,at=1:13, colnames(AvgTOR),cex.axis=0.7,las=2)
    arrows(x0=1:length(AvgTOR[i,]), x1=1:length(AvgTOR[i,]),
           y0=AvgTOR[i,]-SDTOR[i,], y1=AvgTOR[i,]+SDTOR[i,],
           lwd=1.5, angle=90, length=0.05,code=3)
}

#' #### Individual graphs with colors
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
          margins=c(8,12))

#' ### Graphical representations of all the samples
#' #### Preparation of the sample data
sel1 <- c("AT1G50030","AT3G18140","AT3G08850","AT5G01770")
AvgTOR <- Avg[match(sel1,rownames(Avg)),]
TOR <- vsd_T0[match(sel1,rownames(vsd_T0)),]
TOR_norm <- TOR/AvgTOR[,1]
TOR_norm <- log2(TOR_norm)
colnames(TOR_norm) <- samples$Conditions


a <- match(colnames(TOR_norm),samples$Conditions)
plot(TOR_norm[4,] ~ a, ylim=c(-0.5,0.5),main="AT5G01770")
b <- c(1,9,12,21,24,15,18,27,30,39,42,33,36)
lines(b,AvgTOR_norm[5,])





#' ## Genes involved in the cell cycle
#' ### Preparation of a list of genes of the cell cycle
sel1 <- read.csv("~/arabidopsis-nutrition-TOR/seidr/cell_cycle_genes.csv", sep=";")[,1]

#' ## Preparation of the expression table for the shortlisted genes
AvgTOR <- Avg[match(sel1,rownames(Avg)),]
SDTOR <- SD[match(sel1,rownames(SD)),]

#' ### Graphical representations
#' #### Combined barplots
barcenters <- barplot(AvgTOR,beside=T, cex.names=0.7,legend.text = sel1,ylim=c(0,4))
arrows(x0=barcenters, x1=barcenters, y0=AvgTOR-SDTOR, y1=AvgTOR+SDTOR,
       lwd=1.5, angle=90, length=0.05,code=3)

#' #### Individual graphs without colors
for (i in 1:length(AvgTOR[,1]))
{
    plot(AvgTOR[i,],ylim=c(min(AvgTOR[i,]-SDTOR[i,]),max(AvgTOR[i,]+SDTOR[i,])),
         xaxt="n",xlab="",ylab=sprintf("gene = %s",rownames(AvgTOR)[i]))
    axis(side = 1,at=1:13, colnames(AvgTOR),cex.axis=0.7,las=2)
    arrows(x0=1:length(AvgTOR[i,]), x1=1:length(AvgTOR[i,]),
           y0=AvgTOR[i,]-SDTOR[i,], y1=AvgTOR[i,]+SDTOR[i,],
           lwd=1.5, angle=90, length=0.05,code=3)
}

#' #### Individual graphs with colors
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



#' #### Heatmap
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

#' ## Histone genes
#' ### List of Histone H3 family members excluding the "AT1G75610" which was not detected and the "AT1G19890" which was found detected only in the T0 sample
sel1 <- c("AT5G65360","AT1G09200","AT3G27360","AT5G10400","AT5G65350","AT5G10390","AT4G40030",
          "AT4G40040","AT5G10980","AT1G13370","AT1G75600","AT1G01370") 



#' ### Preparation of the expression table for the shortlisted genes
AvgTOR <- Avg[match(sel1,rownames(Avg)),]
SDTOR <- SD[match(sel1,rownames(SD)),]

#' ### Graphical representations
#' #### Combined barplots
barcenters <- barplot(AvgTOR,beside=T, cex.names=0.7,legend.text = sel1)
arrows(x0=barcenters, x1=barcenters, y0=AvgTOR-SDTOR, y1=AvgTOR+SDTOR,
       lwd=1.5, angle=90, length=0.05,code=3)

#' #### Separated barplots
for (i in 1:length(AvgTOR[,1])){
    barcenters <- barplot(AvgTOR[i,],beside=T, cex.names=0.7,legend.text = sel1[i], las=2)
    arrows(x0=barcenters, x1=barcenters, y0=AvgTOR[i,]-SDTOR[i,], y1=AvgTOR[i,]+SDTOR[i,],
           lwd=1.5, angle=90, length=0.2,code=3)
}


#' #### Individual graphs without colors
for (i in 1:length(AvgTOR[,1]))
{
    plot(AvgTOR[i,],ylim=c(min(AvgTOR[i,]-SDTOR[i,]),max(AvgTOR[i,]+SDTOR[i,])),
         xaxt="n",xlab="",ylab=sprintf("gene = %s",rownames(AvgTOR)[i]))
    axis(side = 1,at=1:13, colnames(AvgTOR),cex.axis=0.7,las=2)
    arrows(x0=1:length(AvgTOR[i,]), x1=1:length(AvgTOR[i,]),
           y0=AvgTOR[i,]-SDTOR[i,], y1=AvgTOR[i,]+SDTOR[i,],
           lwd=1.5, angle=90, length=0.05,code=3)
}

#' #### Individual graphs with colors
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

#' #### Heatmap
AvgTOR[AvgTOR == 0] <- 0.000001
AvgTOR_norm <- log2(AvgTOR[,] / AvgTOR[,1])
AvgTOR_norm[AvgTOR_norm < -3] <- -3; AvgTOR_norm[AvgTOR_norm > 3] <- 3


heatmap.2(AvgTOR_norm,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,12))



#' ## Genes involved in the response to phosphate starvation
#' ### Preparation of the transporter list (Poirier and Bucher, 2002) and phosphatases (Hanchi, PhD thesis) and Pi starvation responsive genes (Wang et al., 2018, Front. Plant Sci.)
sel1 <- c("AT5G43350","AT5G43370","AT5G43360","AT2G38940","AT2G32830","AT5G43340","AT3G54700",
          "AT1G20860","AT1G76430","AT3G26570","AT5G14040","AT3G48850","AT2G17270","AT5G46110",
          "AT5G33320","AT5G54800","AT1G61800","AT5G17640","AT1G73010","AT1G17710",
          "AT3G09922" ,"AT2G02990","AT3G17790") 

sel1 <- c("AT5G43350","AT5G43370","AT5G43360","AT2G38940","AT2G32830","AT5G43340","AT3G54700",
          "AT1G20860","AT1G76430","AT3G26570","AT5G14040","AT3G48850","AT2G17270","AT5G46110",
          "AT5G33320","AT5G54800","AT1G61800","AT5G17640","AT1G73010","AT1G17710") 

#' ### Preparation of the expression table for the shortlisted genes
AvgTOR <- Avg[match(sel1,rownames(Avg)),]
SDTOR <- SD[match(sel1,rownames(SD)),]

#' ### Graphical representations
#' #### Combined barplots
barcenters <- barplot(AvgTOR,beside=T, cex.names=0.7,legend.text = sel1)
arrows(x0=barcenters, x1=barcenters, y0=AvgTOR-SDTOR, y1=AvgTOR+SDTOR,
       lwd=1.5, angle=90, length=0.05,code=3)

#' #### Separated barplots
for (i in 1:length(AvgTOR[,1])){
    barcenters <- barplot(AvgTOR[i,],beside=T, cex.names=0.7,legend.text = sel1[i], las=2)
    arrows(x0=barcenters, x1=barcenters, y0=AvgTOR[i,]-SDTOR[i,], y1=AvgTOR[i,]+SDTOR[i,],
           lwd=1.5, angle=90, length=0.2,code=3)
}


#' #### Individual graphs without colors
for (i in 1:length(AvgTOR[,1]))
{
    plot(AvgTOR[i,],ylim=c(min(AvgTOR[i,]-SDTOR[i,]),max(AvgTOR[i,]+SDTOR[i,])),
         xaxt="n",xlab="",ylab=sprintf("gene = %s",rownames(AvgTOR)[i]))
    axis(side = 1,at=1:13, colnames(AvgTOR),cex.axis=0.7,las=2)
    arrows(x0=1:length(AvgTOR[i,]), x1=1:length(AvgTOR[i,]),
           y0=AvgTOR[i,]-SDTOR[i,], y1=AvgTOR[i,]+SDTOR[i,],
           lwd=1.5, angle=90, length=0.05,code=3)
}

#' #### Individual graphs with colors
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

#' #### Heatmap
AvgTOR[AvgTOR == 0] <- 0.000001
AvgTOR_norm <- log2(AvgTOR[,] / AvgTOR[,1])
AvgTOR_norm[AvgTOR_norm < -3] <- -3; AvgTOR_norm[AvgTOR_norm > 3] <- 3


heatmap.2(AvgTOR_norm,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,12))

#' ### Preparation of the systemically-regulated genes by Pi starvation (Thibaud et al., 2010)
sel1 <- c('AT5G01220','AT3G02870','AT3G17790','AT5G64000','AT3G52820','AT2G11810','AT1G73010',
          'AT3G03540','AT3G02040','AT5G20410','AT3G05630','AT4G33030','AT2G27190','AT4G00550',
          'AT2G45130','AT5G20150','AT2G26660','AT1G68740','AT2G32830','AT2G38940','AT1G20860',
          'AT3G58810','AT5G56080','AT5G04950','AT1G23020','AT4G19690','AT3G46900','AT4G19680',
          'AT1G09790','AT3G58060','AT5G03570','AT2G28160','AT2G03260','AT1G74770','AT4G09110',
          'AT5G06490','AT1G49390','AT1G72200','AT2G35000','AT3G12900','AT4G30120') 


#' ### Preparation of the expression table for the shortlisted genes
AvgTOR <- Avg[match(sel1,rownames(Avg)),]
SDTOR <- SD[match(sel1,rownames(SD)),]

#' ### Graphical representations
#' #### Combined barplots
barcenters <- barplot(AvgTOR,beside=T, cex.names=0.7,legend.text = sel1)
arrows(x0=barcenters, x1=barcenters, y0=AvgTOR-SDTOR, y1=AvgTOR+SDTOR,
       lwd=1.5, angle=90, length=0.05,code=3)

#' #### Separated barplots
for (i in 1:length(AvgTOR[,1])){
    barcenters <- barplot(AvgTOR[i,],beside=T, cex.names=0.7,legend.text = sel1[i], las=2)
    arrows(x0=barcenters, x1=barcenters, y0=AvgTOR[i,]-SDTOR[i,], y1=AvgTOR[i,]+SDTOR[i,],
           lwd=1.5, angle=90, length=0.2,code=3)
}


#' #### Individual graphs without colors
for (i in 1:length(AvgTOR[,1]))
{
    plot(AvgTOR[i,],ylim=c(min(AvgTOR[i,]-SDTOR[i,]),max(AvgTOR[i,]+SDTOR[i,])),
         xaxt="n",xlab="",ylab=sprintf("gene = %s",rownames(AvgTOR)[i]))
    axis(side = 1,at=1:13, colnames(AvgTOR),cex.axis=0.7,las=2)
    arrows(x0=1:length(AvgTOR[i,]), x1=1:length(AvgTOR[i,]),
           y0=AvgTOR[i,]-SDTOR[i,], y1=AvgTOR[i,]+SDTOR[i,],
           lwd=1.5, angle=90, length=0.05,code=3)
}

#' #### Individual graphs with colors
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

#' #### Heatmap
AvgTOR[AvgTOR == 0] <- 0.000001
AvgTOR_norm <- log2(AvgTOR[,] / AvgTOR[,1])
AvgTOR_norm[AvgTOR_norm < -3] <- -3; AvgTOR_norm[AvgTOR_norm > 3] <- 3


heatmap.2(AvgTOR_norm,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,12))


#' ## Genes deregulated in Nicolai et al., 2006
#' ### List of the genes
sel1 <- c('AT2G03090','AT3G10950','AT3G25250','AT4G30280','AT5G13210','AT5G65110','AT1G12780','AT4G02520','AT3G12970','AT3G47540',
          'AT4G11650','AT4G30490','AT2G29490','AT2G38870','AT3G45970','AT5G22920','AT5G39320','AT1G75380','AT5G39580','AT5G57655',
          'AT1G54100','AT1G68440','AT1G76870','AT2G33150','AT3G04720','AT3G44300','AT3G59270','AT1G51400','AT2G30860','AT3G28180',
          'AT4G37610','AT5G07440','AT5G49360') 



#' ### Preparation of the expression table for the shortlisted genes
AvgTOR <- Avg[match(sel1,rownames(Avg)),]
SDTOR <- SD[match(sel1,rownames(SD)),]

#' ### Graphical representations
#' #### Combined barplots
barcenters <- barplot(AvgTOR,beside=T, cex.names=0.7,legend.text = sel1)
arrows(x0=barcenters, x1=barcenters, y0=AvgTOR-SDTOR, y1=AvgTOR+SDTOR,
       lwd=1.5, angle=90, length=0.05,code=3)

#' #### Separated barplots
for (i in 1:length(AvgTOR[,1])){
    barcenters <- barplot(AvgTOR[i,],beside=T, cex.names=0.7,legend.text = sel1[i], las=2)
    arrows(x0=barcenters, x1=barcenters, y0=AvgTOR[i,]-SDTOR[i,], y1=AvgTOR[i,]+SDTOR[i,],
           lwd=1.5, angle=90, length=0.2,code=3)
}


#' #### Individual graphs without colors
for (i in 1:length(AvgTOR[,1]))
{
    plot(AvgTOR[i,],ylim=c(min(AvgTOR[i,]-SDTOR[i,]),max(AvgTOR[i,]+SDTOR[i,])),
         xaxt="n",xlab="",ylab=sprintf("gene = %s",rownames(AvgTOR)[i]))
    axis(side = 1,at=1:13, colnames(AvgTOR),cex.axis=0.7,las=2)
    arrows(x0=1:length(AvgTOR[i,]), x1=1:length(AvgTOR[i,]),
           y0=AvgTOR[i,]-SDTOR[i,], y1=AvgTOR[i,]+SDTOR[i,],
           lwd=1.5, angle=90, length=0.05,code=3)
}

#' #### Individual graphs with colors
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

#' #### Heatmap
AvgTOR[AvgTOR == 0] <- 0.000001
AvgTOR_norm <- log2(AvgTOR[,] / AvgTOR[,1])
AvgTOR_norm[AvgTOR_norm < -3] <- -3; AvgTOR_norm[AvgTOR_norm > 3] <- 3


heatmap.2(AvgTOR_norm,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,12))









#' ## Preparation of other lists
#' * Preparation of the hexokinase list
sel1 <- c("AT1G05205","AT1G47840","AT1G47845","AT2G19860","AT4G29130")
#' * Preparation of the AtHXK1
sel1 <- c("AT1G47845","AT4G29130")



#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```

