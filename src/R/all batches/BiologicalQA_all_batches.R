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
#' Set the working dir
setwd("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/all-batches")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/all-batches")
#' ```

#' Load libraries
suppressPackageStartupMessages({
    library(data.table)
    library(DESeq2)
    library(gplots)
    library(here)
    library(hyperSpec)
    library(limma)
    library(parallel)
    library(pander)
    library(plotly)
    library(RColorBrewer)
    library(scatterplot3d)
    library(tidyverse)
    library(tximport)
    library(vsn)
})


#' Source some helper functions
source("~/Git/UPSCb/UPSCb-common/src/R/plot.multidensity.R")
source("~/Git/UPSCb/UPSCb-common/src/R/featureSelection.R")

#' Create palettes
pal <- c(brewer.pal(8,"Dark2"),1)
cols <- rainbow(17)
hpal <- colorRampPalette(c("blue","white","red"))(100)

#' Register the default plot margin
mar <- par("mar")

#' # Raw data
#' ## Loading
#' Read the sample information
samples <- read.csv("~/Git/UPSCb/projects/arabidopsis-nutrition-TOR/doc/samples3.csv")
samples <- read.csv("~/Git/UPSCb/projects/arabidopsis-nutrition-TOR/doc/samples4.csv")

# The data

filenames <- list.files("../Salmon", 
                    recursive = TRUE, 
                    pattern = "quant.sf",
                    full.names = TRUE)

#' name them
#' 
names(filenames) <- sub("_S.*","",sapply(strsplit(filenames, "/"), .subset, 3))

#' Reorder the sample data.frame according to the way
#' the results were read in and accounting for tech reps
samples <- samples[match(names(filenames),samples$SciLifeID),]
samples$Conditions <- factor(paste(samples$Timepoint,samples$Nutrition,samples$AZD,sep="_"))
samples$Batch <- factor(substr(samples$SciLifeID,1,8))
    
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

#' ### Display the samples mean raw counts distribution
#' 
#' i.e. the mean raw count of every gene across samples is calculated and displayed on a log10 scale.
#'
plot(density(log10(rowMeans(kg))),col=pal[1],
     main="gene mean raw counts distribution",
     xlab="mean raw counts (log10)")

#' ### Display all samples raw counts distribution
plot.multidensity(log10(kg+1),
                  col=cols[as.integer(samples$Conditions)],
                  legend.x="topright",
                  legend=levels(samples$SampleName),
                  legend.col=cols,
                  legend.lwd=1,legend.cex=.7,
                  main="samples raw counts distribution",
                  xlab="per gene raw counts (log10)")

#' # Data normalisation 
#'  For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample type 
dds <- DESeqDataSetFromMatrix(
    countData = kg,
    colData = samples,
    design = ~ Conditions)

#' Check the size factors (i.e. the sequencing library size effect)
#' There is no big variation, a Variance Stabilizing Transformation can
#' be used (over a Relative Log Transformation)
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
names(sizes) <- colnames(kg)
sizes
boxplot(sizes, main="Sequencing libraries size factor")

#' split by batch effect
boxplot(split(sizes,samples$Batch))

#' ## VST
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' Validation
meanSdPlot(vst)

#' Write out
#write.csv(vst,"/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/analysis_Tom/library-size-normalized_variance-stabilized_data_nutrition.csv")

#' ## Biological QA
#' 
#' ### MDS
#' Verification of the sample distribution by a multivariate approach on the raw data
#' 
#' * Samples
plotMDS(kg,cex=.6,col=cols[as.integer(samples$Conditions)])

#' * Suspect samples
suspect <- c("P11554_155","P11554_237","P13406_105")
sel <- samples$Conditions %in% samples[match(suspect,samples$SciLifeID),"Conditions"]
pander(samples[sel,])

plotMDS(kg[,sel],cex=.6,col=cols[as.integer(samples$Conditions)][sel])

#' * Suspect Batch
plotMDS(kg[,sel],cex=.6,labels=samples$Batch[sel],col=cols[as.integer(samples$Conditions)][sel])

#' ### PCA
#' 
#' Principal Component Analysis on the normalized data
#' 
#' * Colour coded by Nutrition and AZD treatment - all data
conditions1 <- factor(paste(samples$AZD,samples$Nutrition,sep="_"))
pch<-as.integer(factor(samples$Timepoint))+14
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100);percent

#' * 0 hours
t.sel <- samples$Timepoint == 0
conds <- droplevels(conditions1[t.sel])
plot(pc$x[t.sel,1],
     pc$x[t.sel,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(conds)],
     pch=c(8,15,16,17)[as.integer(samples$Batch)[t.sel]],
     main="T0")

legend("top",pch=19,
       col=pal[1:nlevels(conds)],
       legend=levels(conds))

legend("bottomleft",pch=c(8,15,16,17),
       col="black",
       legend=c("Batch1","Batch2","Batch3","Batch4"))

#' * 6 hours
t.sel <- samples$Timepoint == 6
conds <- droplevels(conditions1[t.sel])
plot(pc$x[t.sel,1],
     pc$x[t.sel,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(conds)],
     pch=c(8,15,16,17)[as.integer(samples$Batch)[t.sel]],
     main="T6")

legend("bottom",pch=19,
       col=pal[1:nlevels(conds)],
       legend=levels(conds))

legend("bottomright",pch=c(8,15,16,17),
       col="black",
       legend=c("Batch1","Batch2","Batch3","Batch4"))

#' * 24 hours
t.sel <- samples$Timepoint == 24
conds <- droplevels(conditions1[t.sel])
plot(pc$x[t.sel,1],
     pc$x[t.sel,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(conds)],
     pch=c(8,15,16,17)[as.integer(samples$Batch)[t.sel]],
     main="T24")

legend("top",pch=19,
       col=pal[1:nlevels(conds)],
       legend=levels(conds))

legend("bottomright",pch=c(8,15,16,17),
       col="black",
       legend=c("Batch1","Batch2","Batch3","Batch4"))

#' * 24 hours or 0 hrs
t.sel <- samples$Timepoint %in% c(0,24)
conds <- droplevels(conditions1[t.sel])
plot(pc$x[t.sel,1],
     pc$x[t.sel,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(conds)],
     pch=c(8,15,16,17)[as.integer(samples$Batch)[t.sel]],
     main="T24 and 0hrs")

legend("top",pch=19,
       col=pal[1:nlevels(conds)],
       legend=levels(conds))

legend("bottomright",pch=c(8,15,16,17),
       col="black",
       legend=c("Batch1","Batch2","Batch3","Batch4"))

#' * All timepoints
t.sel <- samples$Timepoint %in% c(0,6,24)
conds <- droplevels(conditions1[t.sel])
plot(pc$x[t.sel,1],
     pc$x[t.sel,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(conds)],
     pch=c(8,15,16,17)[as.integer(samples$Batch)[t.sel]],
     main="All timepoints")

legend("bottomright",pch=19,
       col=pal[1:nlevels(conds)],
       legend=levels(conds))

legend("bottomleft",pch=c(8,15,16,17),
       col="black",
       legend=c("Batch1","Batch2","Batch3","Batch4"))

#' * Batch 2, 3 and 4 only
t.sel <- samples$Batch %in% c("P11554_2","P13406_1", "P14066_1")
conds <- droplevels(conditions1[t.sel])
plot(pc$x[t.sel,1],
     pc$x[t.sel,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(conds)],
     pch=c(8,15,16,17)[as.integer(samples$Batch)[t.sel]],
     cex=c(0.5,1,2)[as.factor(samples$Timepoint)[t.sel]],
     main="All timepoints for only batch2, 3 and 4")


legend("bottomright",pch=19,
       col=pal[1:nlevels(conds)],
       legend=levels(conds))

legend("bottomleft",pch=c(15,16,17),
       col="black",
       legend=c("Batch2","Batch3","Batch4"))

legend("topleft",pch=19,
       cex=c(0.5,1,2),
       col="black",
       legend=c("T0","T6","T24"))





#' # Heatmap
ssel <- samples$SciLifeID[samples$Timepoint %in% c(0,6) &
                              samples$Nutrition %in% c("NPS","NS","T0") &
                              substr(colnames(tx$counts),1,8) %in% c("P11554_2","P13406_1","P14066_1")]
counts <- tx$counts[,substr(colnames(tx$counts),1,8) %in% c("P11554_2","P13406_1")]
counts <- tx$counts[,colnames(tx$counts) %in% factor(ssel)]
sel <- order(apply(log2(counts),1,sd),decreasing=TRUE)[1:1000]
heatmap.2(log2(counts[sel,]+1),labRow = NA,trace = "none")

ssel <- samples$SciLifeID[samples$Timepoint %in% c(0,6) &
                              samples$Nutrition %in% c("NPS","NS") &
                              substr(colnames(tx$counts),1,8) %in% c("P11554_2","P13406_1","P14066_1")]
counts <- tx$counts[,substr(colnames(tx$counts),1,8) %in% c("P11554_2","P13406_1")]
counts <- tx$counts[,colnames(tx$counts) %in% factor(ssel)]
sel <- order(apply(log2(counts),1,sd),decreasing=TRUE)[1:1000]
heatmap.2(log2(counts[sel,]+1),labRow = NA,trace = "none")

#' # Selection of technical replicates
#' ## T0
t.sel <- samples$Timepoint == 0 & samples$Batch %in% c("P11554_2","P13406_1","P14066_1")
conds <- droplevels(conditions1[t.sel])
plot(pc$x[t.sel,1],
     pc$x[t.sel,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     type="n",
     main="All timepoints")

text(pc$x[t.sel,1],
     pc$x[t.sel,2],
     samples$SampleName[t.sel],
     col=pal[as.integer(samples$Batch)[t.sel]])

legend("bottomleft",col=pal[1:nlevels(samples$Batch)],
       pch=15,
       legend=c("Batch1","Batch2","Batch3","Batch4"))

pander(samples[t.sel,][order(samples$SampleName[t.sel]),])

#' ## T6
t.sel <- samples$Timepoint == 6 & samples$Batch %in% c("P11554_2","P13406_1","P14066_1")
conds <- droplevels(conditions1[t.sel])
plot(pc$x[t.sel,1],
     pc$x[t.sel,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     type="n",
     main="All timepoints")

text(pc$x[t.sel,1],
     pc$x[t.sel,2],
     samples$SampleName[t.sel],
     col=pal[as.integer(samples$Batch)[t.sel]])

legend("topleft",col=pal[1:nlevels(samples$Batch)],
       pch=15,
       legend=c("Batch1","Batch2","Batch3","Batch4"))

pander(samples[t.sel,][order(samples$SampleName[t.sel]),])

#' ## T24
t.sel <- samples$Timepoint == 24 & samples$Batch %in% c("P11554_2","P13406_1","P14066_1")
conds <- droplevels(conditions1[t.sel])
plot(pc$x[t.sel,1],
     pc$x[t.sel,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     type="n",
     main="All timepoints")

text(pc$x[t.sel,1],
     pc$x[t.sel,2],
     samples$SampleName[t.sel],
     col=pal[as.integer(samples$Batch)[t.sel]])

legend("bottomleft",col=pal[1:nlevels(samples$Batch)],
       pch=15,
       legend=c("Batch1","Batch2","Batch3","Batch4"))

pander(samples[t.sel,][order(samples$SampleName[t.sel]),])

#' # All samples after removing the technical replicates
samples <- samples[- samples$SciLifeID %in% "P13406_101","P14066_128"]
#' * Selection of the samples
samples <- samples[-c(samples$SciLifeID=="P13406_101",
           samples$SciLifeID=="P13406_102",
           samples$SciLifeID=="P14066_128",
           samples$SciLifeID=="P14066_133",
           samples$SciLifeID=="P14066_131"),]







#' Heatmap normalized to the T0
#norm_data <- read.csv("Tom-normalizedToT0-data-vst-blind.csv",row.names = 1)
#ssel <- samples$SciLifeID[samples$Timepoint %in% c(0,6) &
#                              samples$Nutrition %in% c("NPS","NS","T0") &
#                              substr(colnames(tx$counts),1,8) %in% c("P11554_2","P13406_1","P14066_1")]
#counts <- norm_data[,colnames(norm_data) %in% factor(ssel)]
#sel <- order(apply(log2(counts),1,sd),decreasing=TRUE)[1:1000]
#heatmap.2(log2(counts[sel,]+1),labRow = NA,trace = "none")
#
#ssel <- samples$SciLifeID[samples$Timepoint %in% c(0,6) &
#                              samples$Nutrition %in% c("NPS","NS") &
#                              substr(colnames(tx$counts),1,8) %in% c("P11554_2","P13406_1")]
#counts <- tx$counts[,substr(colnames(tx$counts),1,8) %in% c("P11554_2","P13406_1")]
#counts <- tx$counts[,colnames(tx$counts) %in% factor(ssel)]
#sel <- order(apply(log2(counts),1,sd),decreasing=TRUE)[1:1000]
#heatmap.2(log2(counts[sel,]+1),labRow = NA,trace = "none")









#zero <- samples$SciLifeID[samples$Nutrition == "T0"]
#count_zero <- tx$counts[,colnames(tx$counts) %in% factor(zero) &
#                                      substr(colnames(tx$counts),1,8) %in% c("P11554_2","P13406_1")]
#c <- apply(count_zero,1,mean)
#count_zero <- cbind(count_zero,average=apply(count_zero,1,mean))
#counts2 <- counts / c
#counts3 <- log2(counts2)
#count_zero <- count_zero / c
#sel <- order(apply(log2(counts2),1,sd),decreasing=TRUE)[1:1000]
#hmcol = colorRampPalette(brewer.pal(9, "RdBu"))(100)
#heatmap.2(log2(counts2[sel,]+1),
#          col = hmcol,
#          scale="row",
#          labRow = NA,
#          trace = "none")
#write.table(counts3,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/analysis_Tom/log2_normalized_to_T0.csv")

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
