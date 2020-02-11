#' ---
#' title: "Biological QA of the 6hr samples"
#' author: "Iryna Shutava, Thomas Dobrenel and Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR")
#' ```

#' Load libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(limma)) 
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(vsn))

#' Source some helper functions
source("~/Git/UPSCb/src/R/plot.multidensity.R")
source("~/Git/UPSCb/src/R/featureSelection.R")

#' Create palettes
pal <- c(brewer.pal(8,"Dark2"),1)
cols <- rainbow(17)
hpal <- colorRampPalette(c("blue","white","red"))(100)

#' Register the default plot margin
mar <- par("mar")

#' # Raw data
#' ## Loading
#' Read the sample information
samples <- read.csv("~/Git/UPSCb/projects/arabidopsis-nutrition-TOR/doc/samples2.csv")
#' Remove the samples of the 1st batch
samples <- samples[! substr(samples$SciLifeID,1,8) == "P11554_1",]
#' Remove the samples of the 24hr timepoint
samples <- samples[! samples$Timepoint == 24,]

# The data
filenames <- list.files("Salmon", 
                    recursive = TRUE, 
                    pattern = "quant.sf",
                    full.names = TRUE)

#' name them
#' 
names(filenames) <- sub("_S.*","",sapply(strsplit(filenames, "/"), .subset, 2))

#' Reorder the sample data.frame according to the way
#' the results were read in and accounting for tech reps
samples <- samples[match(names(filenames),samples$SciLifeID),]
samples
#' Remove the NAs
samples <- na.omit(samples)
#' Rename rows
rownames(samples) <- samples$SciLifeID

samples$Conditions <- factor(paste(samples$Timepoint,samples$Nutrition,samples$AZD,sep="_"))
samples$Batch <- factor(substr(samples$SciLifeID,1,8))

#' Clean-up the sample list in the filenames
filenames <- filenames[match(samples$SciLifeID,names(filenames))]
    
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

#' ### PCA
#' 
#' Principal Component Analysis on the normalized data
#' 
#' * Colour coded by Nutrition and AZD treatment - all data
conditions1 <- factor(paste(samples$AZD,samples$Nutrition,sep="_"))
pch<-as.integer(factor(samples$Timepoint))+14
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100);percent



#' * Drawing the PCA
t.sel <- samples$Timepoint %in% c(0,24)
conds <- droplevels(conditions1[t.sel])
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(samples$Nutrition)],
     pch=c(8,17,19)[as.integer(samples$AZD)],
     main="T6 and 0hrs")

legend("bottomright",pch=19,
       col=pal[1:nlevels(samples$Nutrition)],
       legend=levels(samples$Nutrition))

legend("top",pch=c(8,17,19),
       col="black",
       legend=levels(samples$AZD))

#' # Heatmaps
#' ## All nutritions
counts <- tx$counts
sel <- order(apply(log2(counts),1,sd),decreasing=TRUE)[1:1000]
heatmap.2(log2(counts[sel,]+1),labRow = NA,trace = "none")

#' ## Only +/- phosphorus with T0
ssel <- samples$SciLifeID[samples$Timepoint %in% c(0,6) &
                              samples$Nutrition %in% c("NPS","NS","T0") &
                              substr(colnames(tx$counts),1,8) %in% c("P11554_2","P13406_1")]
counts <- tx$counts[,substr(colnames(tx$counts),1,8) %in% c("P11554_2","P13406_1")]
counts <- tx$counts[,colnames(tx$counts) %in% factor(ssel)]
sel <- order(apply(log2(counts),1,sd),decreasing=TRUE)[1:1000]
heatmap.2(log2(counts[sel,]+1),labRow = NA,trace = "none")
hm <- heatmap.2(log2(counts[sel,]+1),labRow = NA,trace = "none")

#' ## Only +/- phosphorus without T0
ssel <- samples$SciLifeID[samples$Timepoint %in% c(0,6) &
                              samples$Nutrition %in% c("NPS","NS") &
                              substr(colnames(tx$counts),1,8) %in% c("P11554_2","P13406_1")]
counts <- tx$counts[,substr(colnames(tx$counts),1,8) %in% c("P11554_2","P13406_1")]
counts <- tx$counts[,colnames(tx$counts) %in% factor(ssel)]
sel <- order(apply(log2(counts),1,sd),decreasing=TRUE)[1:1000]
heatmap.2(log2(counts[sel,]+1),labRow = NA,trace = "none")















##' Heatmap normalized to the T0
#norm_data <- read.csv("Tom-normalizedToT0-data-vst-blind.csv",row.names = 1)
#ssel <- samples$SciLifeID[samples$Timepoint %in% c(0,6) &
#                              samples$Nutrition %in% c("NPS","NS","T0") &
#                              substr(colnames(tx$counts),1,8) %in% c("P11554_2","P13406_1")]
#counts <- norm_data[,colnames(norm_data) %in% factor(ssel)]
#sel <- order(apply(log2(counts),1,sd),decreasing=TRUE)[1:1000]
#heatmap.2(log2(counts[sel,]+1),labRow = NA,trace = "none")

#ssel <- samples$SciLifeID[samples$Timepoint %in% c(0,6) &
#                              samples$Nutrition %in% c("NPS","NS") &
#                              substr(colnames(tx$counts),1,8) %in% c("P11554_2","P13406_1")]
#counts <- tx$counts[,substr(colnames(tx$counts),1,8) %in% c("P11554_2","P13406_1")]
#counts <- tx$counts[,colnames(tx$counts) %in% factor(ssel)]
#sel <- order(apply(log2(counts),1,sd),decreasing=TRUE)[1:1000]
#heatmap.2(log2(counts[sel,]+1),labRow = NA,trace = "none")






#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
