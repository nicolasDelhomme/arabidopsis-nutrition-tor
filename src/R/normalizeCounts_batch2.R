#' ---
#' title: "Analysis of RNAseq"
#' author: "Thomas Dobrenel (adapted from Nicolas Delhomme & Iryna Shutava)"
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
#' Libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(LSD))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(VennDiagram))
suppressMessages(source("~/Git/UPSCb/src/R/volcanoPlot.R"))

#' Source some helper functions
source("~/Git/UPSCb/src/R/plot.multidensity.R")
source("~/Git/UPSCb/src/R/featureSelection.R")

#' Create a palette
pal <- brewer.pal(9,"Paired")

#' Register the default plot margin
mar <- par("mar")

#' # Get the count data

#' ## Raw data
#' ### Loading
#' Read the sample information

kg <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/analysis_Tom/batch2-unormalised-gene-expression_data.csv", header = TRUE)
load("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/analysis_Tom/batch2_counts.rda")

#' summarise to genes
#tx2gene <- data.frame(TXID=rownames(tx$counts),
#                      GENEID=sub("\\.[0-9]+","",rownames(tx$counts)))

#gx <- summarizeToGene(tx,tx2gene=tx2gene)

#kg <- round(gx$counts) 

#' Sanity check
stopifnot(all(colnames(kg) == samples$SciLifeID))

#' ### Preliminary validations
#' #### Check for the genes that are never expressed
sel <- rowSums(kg) == 0 
# str(sel)= logical vector
# sums the counts for one gene in each samples
# show TRUE or FALSE for each gene if the gene is never expressed (has 0 counts in every time samples)

#' #### Check for the number of genes never expressed
sprintf("%s%% (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(kg),digits=1),
        sum(sel),
        nrow(kg))
#in our experiment the number of genes and transcript are the same (one transcript for one gene)

#' #### Display the samples mean raw counts distribution
#' 
#' i.e. the mean raw count of every gene across samples is calculated and displayed on a log10 scale.
#'
plot(density(log10(rowMeans(kg))),col=pal[1],
     main="gene mean raw counts distribution",
     xlab="mean raw counts (log10)")
#density=likelihood, chance
# plot the likelihood that the mean of genes has this mean of raw counts 
# give an idea of the read depth, here the majority of gene maps ~30 reads (1.5 log10) 

#' #### Display all samples raw counts distribution
cols <- sample(pal,nlevels(samples$SampleName),replace = TRUE)
plot.multidensity(mclapply(1:ncol(kg),function(k){log10(kg)[,k]},mc.cores=16L),
                  col=cols[as.integer(samples$User.ID)],
                  legend.x="topright",
                  legend=levels(samples$User.ID),
                  legend.col=cols,
                  legend.lwd=2,
                  main="samples raw counts distribution",
                  xlab="per gene raw counts (log10)")



#'# Verification of the sample distribution by a multivariate approach
#'## Use of a Multidimensional scaling plot of distances
plotMDS(kg, cex=0.7)

#'## PCA analysis
#'### Establishment of the PCA
conditions1 <- factor(paste(samples$AZD,samples$Nutrition))
pch<-as.integer(factor(samples$Timepoint))+14
pc <- prcomp(t(kg))

percent <- round(summary(pc)$importance[2,]*100);percent

#'### Verification of working in Batch2 only

factor(substr(sub(".*_","",samples$SciLifeID),1,1))

plot(0,0,type="n",xlim=range(pc$x[,1]),ylim=range(pc$x[,2]))
text(pc$x[,1],pc$x[,2],cex=.7,substr(sub(".*_","",samples$SciLifeID),1,1),adj = c(1,-1))

#' Looking for conditions coming from different batches
samples[samples$AZD=="DMSO" & samples$Nutrition == "NP",]
samples[samples$AZD=="DMSO" & samples$Nutrition == "NPS",]
samples[samples$AZD=="DMSO" & samples$Nutrition == "PKS",]
samples[samples$AZD=="DMSO" & samples$Nutrition == "NS",]
samples[samples$AZD=="AZD" & samples$Nutrition == "NP",]
samples[samples$AZD=="AZD" & samples$Nutrition == "NPS",]
samples[samples$AZD=="AZD" & samples$Nutrition == "PKS",]
samples[samples$AZD=="AZD" & samples$Nutrition == "NS",]
samples[samples$AZD=="0",]

#'### PCA for biological meaning
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(factor(conditions1))],
     pch=pch)
legend("bottomleft",pch=19,
       col=pal[as.integer(levels(factor(as.integer(conditions1))))],
       legend=levels(factor(conditions1)))
legend("topleft",pch=as.integer(levels(factor(pch))),
       col="black",
       legend=levels(factor(samples$Timepoint)))


mar=c(5.1,4.1,4.1,2.1)
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=pal[as.integer(factor(conditions1))],
              pch=pch)
legend("topleft",pch=19,
       col=pal[as.integer(levels(factor(as.integer(conditions1))))],
       legend=levels(factor(conditions1)))
legend("topright",pch=as.integer(levels(factor(pch))),
       col="black",
       legend=levels(factor(samples$Timepoint)))

par(mar=c(5.1, 4.1, 4.1, 2.1))

#'# Biological analysis


#' For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate
#'
#' Create the dds object, without giving any prior on the design
#
#'## Data normalisation to T0
#'### Normalization procedure
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#---------------------------------- Analysis with T0 as a reference ----------------------------------------
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
conditions <- factor(paste(samples$Timepoint,".",samples$Nutrition,".",samples$AZD))

conditions <- relevel(conditions,ref = "0 . T0 . 0")

dds.kg <- DESeqDataSetFromMatrix(
    countData = kg,
    colData = data.frame(condition=conditions),
    design = ~ condition)

#' Check the size factors (i.e. the sequencing library size effect)
dds.kg <- estimateSizeFactors(dds.kg)
sizes.kg <- sizeFactors(dds.kg)
# names(sizes.kt) <- colnames(kt)
# not useful sizes.kt already have the right names
pander(sizes.kg)
boxplot(sizes.kg, main="Sequencing libraries size factor")
abline(h=1,lty=2,col="gray")
# There is a small variation between the mediane and the value 1.0

#' #### Variance Stabilising Transformation
#' ##### Blind
system.time(vsd.kg <- varianceStabilizingTransformation(dds.kg, blind=TRUE))
vst_blind <- assay(vsd.kg)

#colnames(vst.kt) <- colnames(kt)
#already have the right names
vst_blind <- vst_blind - min(vst_blind) 
# substracts from all values the smaller value, so the smaller value is now 0

#' Validate the VST 
meanSdPlot(vst_blind[rowSums(vst_blind)>0,]) #mean variance stabilized between 0.5 and 1
meanSdPlot(log2(kg+1)) #without Variance Stabilizing Transformation: variation around 1
meanSdPlot(log2(counts(dds.kg,normalized=TRUE)+1)) #normalized read depth: variance is not stabilized ("curvy")
#VST validated: mean variance stabilized around 0.5


write.csv(vst_blind,"/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Tom-normalizedToT0-data-vst-blind.csv")

#save(kg, vst_blind, file="vst_blind.rda")
save(kg, vst_blind, file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/vst_blind.rda")

#' ##### vst aware (calculated accordingly to samplingDates)

#' Calculate vst aware with DESeq2 taking the date as a parameter
#
#' ##### Model-aware
system.time(vsd2 <- varianceStabilizingTransformation(dds.kg, blind=FALSE))
vst_aware <- assay(vsd2)

#colnames(vst.kt) <- colnames(kt)
#already have the right names
vst_aware <- vst_aware - min(vst_aware) 

# substracts from all values the smaller value, so the smaller value is now 0

#' Validate the VST 
meanSdPlot(vst_aware[rowSums(vst_aware)>0,]) #mean variance stabilized between 0 and 0.5
meanSdPlot(log2(kg+1)) #without Variance Stabilizing Transformation: variation around 1
meanSdPlot(log2(counts(dds.kg,normalized=TRUE)+1)) #normalized read depth: variance is not stabilized ("curvy")

# vst_aware validated: mean variance is way lower when vst is calculated "aware of the date factor" compare to before when vst was calculated "blind"
# stabilized around 0.2 instead of 0.5 when vst was blind


write.csv(vst_aware,"/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Tom-normalized-data-vst-aware.csv")
#save(kg, vst_aware, file="vst_aware.rda")
#load("vst_aware.rda")

#' Differential expression
dds.kg <- DESeq(dds.kg)

#' which results
resultsNames(dds.kg)

#' ### Analysis of 6hrs timepoint
#--------------------------------------------------------------------------------------------------------------------
#-------------------------------------------- Analysis of 6hrs ------------------------------------------------------
#' #### Effect of phosphorus starvation
# ------------------------------ Extract results from 6hrs of Phosphorus starvation vs T0 ---------------------------
res <- results(dds.kg,name="condition_6...NS...DMSO_vs_0...T0...0")

alpha=0.01

#' ##### The log2 fold-change range is also rather tight.
plotMA(res)

#' ##### Plot the log odds vs. log2 fold change
#' 
#' The volcano plot is very flat, cross-validating the results of the MA plot
volcanoPlot(res,alpha=alpha)

#' ##### Plot the adjusted p-value histogram
#' 
#' Which is as expected, over-enriched for adjusted p-values of 1
hist(res$padj,breaks=seq(0,1,.01))

# select genes that are significant
alpha=0.01
l2fc = 0.5
sel <- ! is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= l2fc

# list of DE genes that are significant
six_NS_DMSO_vs_T0 <- rownames(res[sel,])

#' ##### DE genes data export
#' 
write.csv(res,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/6_NS_DMSO_vs_T0-DEgenes_all.csv")
write.csv(res[sel,],file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/6_NS_DMSO_vs_T0-DEgenes__significant.csv")

#' #### Effect of full medium resupply
# ------------------------------ Extract results from 6hrs of No starvation vs T0 ---------------------------
res <- results(dds.kg,name="condition_6...NPS...DMSO_vs_0...T0...0")

alpha=0.01

#' ##### The log2 fold-change range is also rather tight.
plotMA(res)

#' ##### Plot the log odds vs. log2 fold change
#' 
#' The volcano plot is very flat, cross-validating the results of the MA plot
volcanoPlot(res,alpha=alpha)

#' ##### Plot the adjusted p-value histogram
#' 
#' Which is as expected, over-enriched for adjusted p-values of 1
hist(res$padj,breaks=seq(0,1,.01))

# select genes that are significant
alpha=0.01
l2fc = 0.5
sel <- ! is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= l2fc

# list of DE genes that are significant
six_NPS_DMSO_vs_T0 <- rownames(res[sel,])

#' ##### DE genes data export
#' 
write.csv(res,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/6_NPS_DMSO_vs_T0-DEgenes_all.csv")
write.csv(res[sel,],file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/6_NPS_DMSO_vs_T0-DEgenes__significant.csv")

#' #### Effect of phosphorus starvation and AZD treatment
# ------------------------------ Extract results from 6hrs of AZD and Phosphorus starvation vs T0 ---------------------------
res <- results(dds.kg,name="condition_6...NS...AZD_vs_0...T0...0")

alpha=0.01

#' ##### The log2 fold-change range is also rather tight.
plotMA(res)

#' ##### Plot the log odds vs. log2 fold change
#' 
#' The volcano plot is very flat, cross-validating the results of the MA plot
volcanoPlot(res,alpha=alpha)

#' ##### Plot the adjusted p-value histogram
#' 
#' Which is as expected, over-enriched for adjusted p-values of 1
hist(res$padj,breaks=seq(0,1,.01))

# select genes that are significant
alpha=0.01
l2fc = 0.5
sel <- ! is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= l2fc

# list of DE genes that are significant
six_NS_AZD_vs_T0 <- rownames(res[sel,])

#' ##### DE genes data export
#' 
write.csv(res,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/6_NS_AZD_vs_T0-DEgenes_all.csv")
write.csv(res[sel,],file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/6_NS_AZD_vs_T0-DEgenes__significant.csv")

#' #### Effect of AZD treatment
# ------------------------------ Extract results from 6hrs of AZD vs T0 ---------------------------
res <- results(dds.kg,name="condition_6...NPS...AZD_vs_0...T0...0")

alpha=0.01

#' ##### The log2 fold-change range is also rather tight.
plotMA(res)

#' ##### Plot the log odds vs. log2 fold change
#' 
#' The volcano plot is very flat, cross-validating the results of the MA plot
volcanoPlot(res,alpha=alpha)

#' ##### Plot the adjusted p-value histogram
#' 
#' Which is as expected, over-enriched for adjusted p-values of 1
hist(res$padj,breaks=seq(0,1,.01))

# select genes that are significant
alpha=0.01
l2fc = 0.5
sel <- ! is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= l2fc

# list of DE genes that are significant
six_NPS_AZD_vs_T0 <- rownames(res[sel,])

#' ##### DE genes data export
#' 
write.csv(res,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/6_NPS_AZD_vs_T0-DEgenes_all.csv")
write.csv(res[sel,],file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/6_NPS_AZD_vs_T0-DEgenes__significant.csv")

#' #### Effect of nitrogen starvation
# ------------------------------ Extract results from 6hrs of nitrogen starvation vs T0 ---------------------------
res <- results(dds.kg,name="condition_6...PKS...DMSO_vs_0...T0...0")

alpha=0.01

#' ##### The log2 fold-change range is also rather tight.
plotMA(res)

#' ##### Plot the log odds vs. log2 fold change
#' 
#' The volcano plot is very flat, cross-validating the results of the MA plot
volcanoPlot(res,alpha=alpha)

#' ##### Plot the adjusted p-value histogram
#' 
#' Which is as expected, over-enriched for adjusted p-values of 1
hist(res$padj,breaks=seq(0,1,.01))

# select genes that are significant
alpha=0.01
l2fc = 0.5
sel <- ! is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= l2fc

# list of DE genes that are significant
six_PKS_DMSO_vs_T0 <- rownames(res[sel,])

#' ##### DE genes data export
#' 
write.csv(res,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/6_PKS_DMSO_vs_T0-DEgenes_all.csv")
write.csv(res[sel,],file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/6_PKS_DMSO_vs_T0-DEgenes__significant.csv")

#' #### Effect of sucrose starvation
# ------------------------------ Extract results from 6hrs of sucrose starvation vs T0 ---------------------------
res <- results(dds.kg,name="condition_6...NP...DMSO_vs_0...T0...0")

alpha=0.01

#' ##### The log2 fold-change range is also rather tight.
plotMA(res)

#' ##### Plot the log odds vs. log2 fold change
#' 
#' The volcano plot is very flat, cross-validating the results of the MA plot
volcanoPlot(res,alpha=alpha)

#' ##### Plot the adjusted p-value histogram
#' 
#' Which is as expected, over-enriched for adjusted p-values of 1
hist(res$padj,breaks=seq(0,1,.01))

# select genes that are significant
alpha=0.01
l2fc = 0.5
sel <- ! is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= l2fc

# list of DE genes that are significant
six_NP_DMSO_vs_T0 <- rownames(res[sel,])

#' ##### DE genes data export
#' 
write.csv(res,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/6_NP_DMSO_vs_T0-DEgenes_all.csv")
write.csv(res[sel,],file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/6_NP_DMSO_vs_T0-DEgenes__significant.csv")

#' #### Effect of nitrogen starvation and AZD treatment
# ------------------------------ Extract results from 6hrs of PKS & AZD vs T0 ---------------------------
res <- results(dds.kg,name="condition_6...PKS...AZD_vs_0...T0...0")

alpha=0.01

#' ##### The log2 fold-change range is also rather tight.
plotMA(res)

#' ##### Plot the log odds vs. log2 fold change
#' 
#' The volcano plot is very flat, cross-validating the results of the MA plot
volcanoPlot(res,alpha=alpha)

#' ##### Plot the adjusted p-value histogram
#' 
#' Which is as expected, over-enriched for adjusted p-values of 1
hist(res$padj,breaks=seq(0,1,.01))

# select genes that are significant
alpha=0.01
l2fc = 0.5
sel <- ! is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= l2fc

# list of DE genes that are significant
six_PKS_AZD_vs_T0 <- rownames(res[sel,])

#' ##### DE genes data export
#' 
write.csv(res,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/6_PKS_AZD_vs_T0-DEgenes_all.csv")
write.csv(res[sel,],file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/6_PKS_AZD_vs_T0-DEgenes__significant.csv")

#' #### Effect of sucrose starvation and AZD treatment
# ------------------------------ Extract results from 6hrs of NP & AZD vs T0 ---------------------------
res <- results(dds.kg,name="condition_6...NP...AZD_vs_0...T0...0")

alpha=0.01

#' ##### The log2 fold-change range is also rather tight.
plotMA(res)

#' ##### Plot the log odds vs. log2 fold change
#' 
#' The volcano plot is very flat, cross-validating the results of the MA plot
volcanoPlot(res,alpha=alpha)

#' ##### Plot the adjusted p-value histogram
#' 
#' Which is as expected, over-enriched for adjusted p-values of 1
hist(res$padj,breaks=seq(0,1,.01))

# select genes that are significant
alpha=0.01
l2fc = 0.5
sel <- ! is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= l2fc

# list of DE genes that are significant
six_NP_AZD_vs_T0 <- rownames(res[sel,])

#' ##### DE genes data export
#' 
write.csv(res,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/6_NP_AZD_vs_T0-DEgenes_all.csv")
write.csv(res[sel,],file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/6_NP_AZD_vs_T0-DEgenes__significant.csv")



#' #### Venn diagrams
#' ##### VennDiagram of +/- Phoshorus +/- AZD

SixH <- VennDiagram::calculate.overlap(list(six_NS_DMSO_vs_T0, six_NS_AZD_vs_T0,six_NPS_AZD_vs_T0,six_NPS_DMSO_vs_T0))
grid.newpage()
svg("Venn1.svg", width = 7, height = 7)
grid.draw(venn.diagram(list(six_NS_DMSO_vs_T0, six_NS_AZD_vs_T0,six_NPS_AZD_vs_T0,six_NPS_DMSO_vs_T0),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("PStarv","PStarv_AZD","AZD","NoStarv")))
dev.off()
grid.draw(venn.diagram(list(six_NS_DMSO_vs_T0, six_NS_AZD_vs_T0,six_NPS_AZD_vs_T0,six_NPS_DMSO_vs_T0),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("PStarv","PStarv_AZD","AZD","NoStarv")))

#' ##### VennDiagram of all nutrition in DMSO

SixH <- VennDiagram::calculate.overlap(list(six_NS_DMSO_vs_T0,
                                            six_NPS_DMSO_vs_T0,
                                            six_PKS_DMSO_vs_T0,
                                            six_NP_DMSO_vs_T0))
grid.newpage()
svg("Venn2.svg", width = 7, height = 7)
grid.draw(venn.diagram(list(six_NS_DMSO_vs_T0,
                            six_NPS_DMSO_vs_T0,
                            six_PKS_DMSO_vs_T0,
                            six_NP_DMSO_vs_T0),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("Pstarv","full","Nstarv","Sstarv")))
dev.off()
grid.draw(venn.diagram(list(six_NS_DMSO_vs_T0,
                            six_NPS_DMSO_vs_T0,
                            six_PKS_DMSO_vs_T0,
                            six_NP_DMSO_vs_T0),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("Pstarv","full","Nstarv","Sstarv")))

#' ##### VennDiagram of all nutrition in AZD


SixH <- VennDiagram::calculate.overlap(list(six_NS_AZD_vs_T0,
                                            six_NPS_AZD_vs_T0,
                                            six_PKS_AZD_vs_T0,
                                            six_NP_AZD_vs_T0))
grid.newpage()
svg("Venn3.svg", width = 7, height = 7)
grid.draw(venn.diagram(list(six_NS_AZD_vs_T0,
                            six_NPS_AZD_vs_T0,
                            six_PKS_AZD_vs_T0,
                            six_NP_AZD_vs_T0),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("Pstarv & AZD","full & AZD","Nstarv & AZD","Sstarv & AZD")))
dev.off()
grid.draw(venn.diagram(list(six_NS_AZD_vs_T0,
                            six_NPS_AZD_vs_T0,
                            six_PKS_AZD_vs_T0,
                            six_NP_AZD_vs_T0),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("Pstarv & AZD","full & AZD","Nstarv & AZD","Sstarv & AZD")))

#'## Data normalisation to NPS_DMSO
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#---------------------------------- Analysis with NPS_DMSO as a reference ----------------------------------------
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
conditions <- factor(paste(samples$Timepoint,".",samples$Nutrition,".",samples$AZD))

conditions <- relevel(conditions,ref = "6 . NPS . DMSO")

dds.kg <- DESeqDataSetFromMatrix(
    countData = kg,
    colData = data.frame(condition=conditions),
    design = ~ condition)

#' Check the size factors (i.e. the sequencing library size effect)
dds.kg <- estimateSizeFactors(dds.kg)
sizes.kg <- sizeFactors(dds.kg)
# names(sizes.kt) <- colnames(kt)
# not useful sizes.kt already have the right names
pander(sizes.kg)
boxplot(sizes.kg, main="Sequencing libraries size factor")
abline(h=1,lty=2,col="gray")
# There is a small variation between the mediane and the value 1.0

#' ### Variance Stabilising Transformation
#' #### Blind
system.time(vsd.kg <- varianceStabilizingTransformation(dds.kg, blind=TRUE))
vst_blind <- assay(vsd.kg)

#colnames(vst.kt) <- colnames(kt)
#already have the right names
vst_blind <- vst_blind - min(vst_blind) 
# substracts from all values the smaller value, so the smaller value is now 0

#' Validate the VST 
meanSdPlot(vst_blind[rowSums(vst_blind)>0,]) #mean variance stabilized between 0.5 and 1
meanSdPlot(log2(kg+1)) #without Variance Stabilizing Transformation: variation around 1
meanSdPlot(log2(counts(dds.kg,normalized=TRUE)+1)) #normalized read depth: variance is not stabilized ("curvy")
#VST validated: mean variance stabilized around 0.5


write.csv(vst_blind,"/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Tom-normalizedTo6_NPS_DMSO-data-vst-blind.csv")

#save(kg, vst_blind, file="vst_blind.rda")
save(kg, vst_blind, file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/vst_blind.rda")

#' #### vst aware (calculated accordingly to samplingDates)

#' Calculate vst aware with DESeq2 taking the date as a parameter
#
#' #### Model-aware
system.time(vsd2 <- varianceStabilizingTransformation(dds.kg, blind=FALSE))
vst_aware <- assay(vsd2)

#colnames(vst.kt) <- colnames(kt)
#already have the right names
vst_aware <- vst_aware - min(vst_aware) 

# substracts from all values the smaller value, so the smaller value is now 0

#' Validate the VST 
meanSdPlot(vst_aware[rowSums(vst_aware)>0,]) #mean variance stabilized between 0 and 0.5
meanSdPlot(log2(kg+1)) #without Variance Stabilizing Transformation: variation around 1
meanSdPlot(log2(counts(dds.kg,normalized=TRUE)+1)) #normalized read depth: variance is not stabilized ("curvy")

# vst_aware validated: mean variance is way lower when vst is calculated "aware of the date factor" compare to before when vst was calculated "blind"
# stabilized around 0.2 instead of 0.5 when vst was blind


write.csv(vst_aware,"/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Tom-normalizedTo_6_NPS_DMSO-data-vst-aware.csv")
#save(kg, vst_aware, file="vst_aware.rda")
#load("vst_aware.rda")

#' Differential expression
dds.kg <- DESeq(dds.kg)

#' which results
resultsNames(dds.kg)

#' ### Effect of AZD treatment in NPS
# ------------------------------ Extract results from 6hrs of AZD vs DMSO ---------------------------
res <- results(dds.kg,name="condition_6...NPS...AZD_vs_6...NPS...DMSO")

alpha=0.01

#' #### The log2 fold-change range is also rather tight.
plotMA(res)

#' #### Plot the log odds vs. log2 fold change
#
#' The volcano plot is very flat, cross-validating the results of the MA plot
volcanoPlot(res,alpha=alpha)

#' #### Plot the adjusted p-value histogram
#' 
#' Which is as expected, over-enriched for adjusted p-values of 1
hist(res$padj,breaks=seq(0,1,.01))

# select genes that are significant
alpha=0.01
l2fc = 0.5
sel <- ! is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= l2fc

# list of DE genes that are significant
AZDvsDMSO_in_NPS <- rownames(res[sel,])

#' #### DE genes data export
#' 
write.csv(res,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/AZDvsDMSO_in_NPS-DEgenes_all.csv")
write.csv(res[sel,],file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/AZDvsDMSO_in_NPS-DEgenes__significant.csv")

#' ### Effect of Phosphorus starvation
# ------------------------------ Extract results from 6hrs of NS vs NPS ---------------------------
res <- results(dds.kg,name="condition_6...NS...DMSO_vs_6...NPS...DMSO")

alpha=0.01

#' #### The log2 fold-change range is also rather tight.
plotMA(res)

#' #### Plot the log odds vs. log2 fold change
#
#' The volcano plot is very flat, cross-validating the results of the MA plot
volcanoPlot(res,alpha=alpha)

#' #### Plot the adjusted p-value histogram
#' 
#' Which is as expected, over-enriched for adjusted p-values of 1
hist(res$padj,breaks=seq(0,1,.01))

# select genes that are significant
alpha=0.01
l2fc = 0.5
sel <- ! is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= l2fc

# list of DE genes that are significant
NSvsNPS_in_DMSO <- rownames(res[sel,])

#' #### DE genes data export
#' 
write.csv(res,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/NSvsNPS_in_DMSO-DEgenes_all.csv")
write.csv(res[sel,],file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/NSvsNPS_in_DMSO-DEgenes__significant.csv")


#' ### Effect of Nitrogen starvation
# ------------------------------ Extract results from 6hrs of PKS vs NPS ---------------------------
res <- results(dds.kg,name="condition_6...PKS...DMSO_vs_6...NPS...DMSO")

alpha=0.01

#' #### The log2 fold-change range is also rather tight.
plotMA(res)

#' #### Plot the log odds vs. log2 fold change
#
#' The volcano plot is very flat, cross-validating the results of the MA plot
volcanoPlot(res,alpha=alpha)

#' #### Plot the adjusted p-value histogram
#' 
#' Which is as expected, over-enriched for adjusted p-values of 1
hist(res$padj,breaks=seq(0,1,.01))

# select genes that are significant
alpha=0.01
l2fc = 0.5
sel <- ! is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= l2fc

# list of DE genes that are significant
PKSvsNPS_in_DMSO <- rownames(res[sel,])

#' #### DE genes data export
#' 
write.csv(res,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/PKSvsNPS_in_DMSO-DEgenes_all.csv")
write.csv(res[sel,],file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/APKSvsNPS_in_DMSO-DEgenes__significant.csv")


#' ### Effect of sugar starvation
# ------------------------------ Extract results from 6hrs of NP vs NPS ---------------------------
res <- results(dds.kg,name="condition_6...NP...DMSO_vs_6...NPS...DMSO")

alpha=0.01

#' #### The log2 fold-change range is also rather tight.
plotMA(res)

#' #### Plot the log odds vs. log2 fold change
#
#' The volcano plot is very flat, cross-validating the results of the MA plot
volcanoPlot(res,alpha=alpha)

#' #### Plot the adjusted p-value histogram
#' 
#' Which is as expected, over-enriched for adjusted p-values of 1
hist(res$padj,breaks=seq(0,1,.01))

# select genes that are significant
alpha=0.01
l2fc = 0.5
sel <- ! is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= l2fc

# list of DE genes that are significant
NPvsNPS_in_DMSO <- rownames(res[sel,])

#' #### DE genes data export
#' 
write.csv(res,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/NPvsNPS_in_DMSO-DEgenes_all.csv")
write.csv(res[sel,],file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/NPvsNPS_in_DMSO-DEgenes__significant.csv")


#' ### Effect of AZD treatment and Phosphorus starvation compared to NPS_DMSO
# ------------------------------ Extract results from 6hrs of NS_AZD vs NPS_DMSO ---------------------------
res <- results(dds.kg,name="condition_6...NS...AZD_vs_6...NPS...DMSO")

alpha=0.01

#' #### The log2 fold-change range is also rather tight.
plotMA(res)

#' #### Plot the log odds vs. log2 fold change
#
#' The volcano plot is very flat, cross-validating the results of the MA plot
volcanoPlot(res,alpha=alpha)

#' #### Plot the adjusted p-value histogram
#' 
#' Which is as expected, over-enriched for adjusted p-values of 1
hist(res$padj,breaks=seq(0,1,.01))

# select genes that are significant
alpha=0.01
l2fc = 0.5
sel <- ! is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= l2fc

# list of DE genes that are significant
NS.AZD_vs_NPS.DMSO <- rownames(res[sel,])

#' #### DE genes data export
#' 
write.csv(res,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/NS.AZD_vs_NPS.DMSO-DEgenes_all.csv")
write.csv(res[sel,],file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/NS.AZD_vs_NPS.DMSO-DEgenes__significant.csv")


#' ### Venn Diagram of Phosphorus starvation and AZD treatment compared to NPS_DMSO


SixH <- VennDiagram::calculate.overlap(list(NSvsNPS_in_DMSO,
                                            AZDvsDMSO_in_NPS,
                                            NS.AZD_vs_NPS.DMSO))
grid.newpage()
svg("Venn5.svg", width = 7, height = 7)
grid.draw(venn.diagram(list(NSvsNPS_in_DMSO,
                            AZDvsDMSO_in_NPS,
                            NS.AZD_vs_NPS.DMSO),
                       filename=NULL,
                       col=pal[1:3],
                       category.names=c("Pstarv","AZD","Pstarv & AZD")))
dev.off()
grid.draw(venn.diagram(list(NSvsNPS_in_DMSO,
                            AZDvsDMSO_in_NPS,
                            NS.AZD_vs_NPS.DMSO),
                       filename=NULL,
                       col=pal[1:3],
                       category.names=c("Pstarv","AZD","Pstarv & AZD")))



#' ## Data normalisation to PKS_DMSO
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#---------------------------------- Analysis with PKS_DMSO as a reference ----------------------------------------
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
conditions <- factor(paste(samples$Timepoint,".",samples$Nutrition,".",samples$AZD))

conditions <- relevel(conditions,ref = "6 . PKS . DMSO")

dds.kg <- DESeqDataSetFromMatrix(
    countData = kg,
    colData = data.frame(condition=conditions),
    design = ~ condition)

#' Check the size factors (i.e. the sequencing library size effect)
dds.kg <- estimateSizeFactors(dds.kg)
sizes.kg <- sizeFactors(dds.kg)
# names(sizes.kt) <- colnames(kt)
# not useful sizes.kt already have the right names
pander(sizes.kg)
boxplot(sizes.kg, main="Sequencing libraries size factor")
abline(h=1,lty=2,col="gray")
# There is a small variation between the mediane and the value 1.0

#' ### Variance Stabilising Transformation
#' #### Blind
system.time(vsd.kg <- varianceStabilizingTransformation(dds.kg, blind=TRUE))
vst_blind <- assay(vsd.kg)

#colnames(vst.kt) <- colnames(kt)
#already have the right names
vst_blind <- vst_blind - min(vst_blind) 
# substracts from all values the smaller value, so the smaller value is now 0

#' Validate the VST 
meanSdPlot(vst_blind[rowSums(vst_blind)>0,]) #mean variance stabilized between 0.5 and 1
meanSdPlot(log2(kg+1)) #without Variance Stabilizing Transformation: variation around 1
meanSdPlot(log2(counts(dds.kg,normalized=TRUE)+1)) #normalized read depth: variance is not stabilized ("curvy")
#VST validated: mean variance stabilized around 0.5


write.csv(vst_blind,"/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Tom-normalizedTo6_PKS_DMSO-data-vst-blind.csv")

#save(kg, vst_blind, file="vst_blind.rda")
save(kg, vst_blind, file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/vst_blind.rda")

#' #### vst aware (calculated accordingly to samplingDates)

#' Calculate vst aware with DESeq2 taking the date as a parameter
#
#' #### Model-aware
system.time(vsd2 <- varianceStabilizingTransformation(dds.kg, blind=FALSE))
vst_aware <- assay(vsd2)

#colnames(vst.kt) <- colnames(kt)
#already have the right names
vst_aware <- vst_aware - min(vst_aware) 

# substracts from all values the smaller value, so the smaller value is now 0

#' Validate the VST 
meanSdPlot(vst_aware[rowSums(vst_aware)>0,]) #mean variance stabilized between 0 and 0.5
meanSdPlot(log2(kg+1)) #without Variance Stabilizing Transformation: variation around 1
meanSdPlot(log2(counts(dds.kg,normalized=TRUE)+1)) #normalized read depth: variance is not stabilized ("curvy")

# vst_aware validated: mean variance is way lower when vst is calculated "aware of the date factor" compare to before when vst was calculated "blind"
# stabilized around 0.2 instead of 0.5 when vst was blind


write.csv(vst_aware,"/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Tom-normalizedTo_6_PKS_DMSO-data-vst-aware.csv")
#save(kg, vst_aware, file="vst_aware.rda")
#load("vst_aware.rda")

#' Differential expression
dds.kg <- DESeq(dds.kg)

#' which results
resultsNames(dds.kg)
#' ### Effect of AZD treatment in PKS
# ------------------------------ Extract results from 6hrs of AZD vs DMSO in PKS ---------------------------
res <- results(dds.kg,name="condition_6...PKS...AZD_vs_6...PKS...DMSO")

alpha=0.01

#' #### The log2 fold-change range is also rather tight.
plotMA(res)

#' #### Plot the log odds vs. log2 fold change
#' 
#' The volcano plot is very flat, cross-validating the results of the MA plot
volcanoPlot(res,alpha=alpha)

#' #### Plot the adjusted p-value histogram
#' 
#' Which is as expected, over-enriched for adjusted p-values of 1
hist(res$padj,breaks=seq(0,1,.01))

# select genes that are significant
alpha=0.01
l2fc = 0.5
sel <- ! is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= l2fc

# list of DE genes that are significant
AZDvsDMSO_in_PKS <- rownames(res[sel,])

#' #### DE genes data export
#' 
write.csv(res,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/AZDvsDMSO_in_PKS-DEgenes_all.csv")
write.csv(res[sel,],file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/AZDvsDMSO_in_PKS-DEgenes__significant.csv")

#'## Data normalisation to NS_DMSO
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#---------------------------------- Analysis with NS_DMSO as a reference ----------------------------------------
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
conditions <- factor(paste(samples$Timepoint,".",samples$Nutrition,".",samples$AZD))

conditions <- relevel(conditions,ref = "6 . NS . DMSO")

dds.kg <- DESeqDataSetFromMatrix(
    countData = kg,
    colData = data.frame(condition=conditions),
    design = ~ condition)

#' Check the size factors (i.e. the sequencing library size effect)
dds.kg <- estimateSizeFactors(dds.kg)
sizes.kg <- sizeFactors(dds.kg)
# names(sizes.kt) <- colnames(kt)
# not useful sizes.kt already have the right names
pander(sizes.kg)
boxplot(sizes.kg, main="Sequencing libraries size factor")
abline(h=1,lty=2,col="gray")
# There is a small variation between the mediane and the value 1.0

#' ### Variance Stabilising Transformation
#' #### Blind
system.time(vsd.kg <- varianceStabilizingTransformation(dds.kg, blind=TRUE))
vst_blind <- assay(vsd.kg)

#colnames(vst.kt) <- colnames(kt)
#already have the right names
vst_blind <- vst_blind - min(vst_blind) 
# substracts from all values the smaller value, so the smaller value is now 0

#' Validate the VST 
meanSdPlot(vst_blind[rowSums(vst_blind)>0,]) #mean variance stabilized between 0.5 and 1
meanSdPlot(log2(kg+1)) #without Variance Stabilizing Transformation: variation around 1
meanSdPlot(log2(counts(dds.kg,normalized=TRUE)+1)) #normalized read depth: variance is not stabilized ("curvy")
#VST validated: mean variance stabilized around 0.5


write.csv(vst_blind,"/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Tom-normalizedTo6_NS_DMSO-data-vst-blind.csv")

#save(kg, vst_blind, file="vst_blind.rda")
save(kg, vst_blind, file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/vst_blind.rda")

#' #### vst aware (calculated accordingly to samplingDates)

#' Calculate vst aware with DESeq2 taking the date as a parameter
#' #### Model-aware
system.time(vsd2 <- varianceStabilizingTransformation(dds.kg, blind=FALSE))
vst_aware <- assay(vsd2)

#colnames(vst.kt) <- colnames(kt)
#already have the right names
vst_aware <- vst_aware - min(vst_aware) 

# substracts from all values the smaller value, so the smaller value is now 0

#' Validate the VST 
meanSdPlot(vst_aware[rowSums(vst_aware)>0,]) #mean variance stabilized between 0 and 0.5
meanSdPlot(log2(kg+1)) #without Variance Stabilizing Transformation: variation around 1
meanSdPlot(log2(counts(dds.kg,normalized=TRUE)+1)) #normalized read depth: variance is not stabilized ("curvy")

# vst_aware validated: mean variance is way lower when vst is calculated "aware of the date factor" compare to before when vst was calculated "blind"
# stabilized around 0.2 instead of 0.5 when vst was blind


write.csv(vst_aware,"/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Tom-normalizedTo_6_NS_DMSO-data-vst-aware.csv")
#save(kg, vst_aware, file="vst_aware.rda")
#load("vst_aware.rda")

#' Differential expression
dds.kg <- DESeq(dds.kg)

#' which results
resultsNames(dds.kg)

#' ### Effect of AZD treatment in NS medium
# ------------------------------ Extract results from 6hrs of AZD vs DMSO in NS ---------------------------
res <- results(dds.kg,name="condition_6...NS...AZD_vs_6...NS...DMSO")

alpha=0.01

#' #### The log2 fold-change range is also rather tight.
plotMA(res)

#' #### Plot the log odds vs. log2 fold change
#' 
#' The volcano plot is very flat, cross-validating the results of the MA plot
volcanoPlot(res,alpha=alpha)

#' #### Plot the adjusted p-value histogram
#' 
#' Which is as expected, over-enriched for adjusted p-values of 1
hist(res$padj,breaks=seq(0,1,.01))

# select genes that are significant
alpha=0.01
l2fc = 0.5
sel <- ! is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= l2fc

# list of DE genes that are significant
AZDvsDMSO_in_NS <- rownames(res[sel,])

#' #### DE genes data export
#' 
write.csv(res,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/AZDvsDMSO_in_NS-DEgenes_all.csv")
write.csv(res[sel,],file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/AZDvsDMSO_in_NS-DEgenes__significant.csv")

#'## Data normalisation to NP_DMSO
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#---------------------------------- Analysis with NP_DMSO as a reference ----------------------------------------
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
conditions <- factor(paste(samples$Timepoint,".",samples$Nutrition,".",samples$AZD))

conditions <- relevel(conditions,ref = "6 . NP . DMSO")

dds.kg <- DESeqDataSetFromMatrix(
    countData = kg,
    colData = data.frame(condition=conditions),
    design = ~ condition)

#' Check the size factors (i.e. the sequencing library size effect)
dds.kg <- estimateSizeFactors(dds.kg)
sizes.kg <- sizeFactors(dds.kg)
# names(sizes.kt) <- colnames(kt)
# not useful sizes.kt already have the right names
pander(sizes.kg)
boxplot(sizes.kg, main="Sequencing libraries size factor")
abline(h=1,lty=2,col="gray")
# There is a small variation between the mediane and the value 1.0

#' ### Variance Stabilising Transformation
#' #### Blind
system.time(vsd.kg <- varianceStabilizingTransformation(dds.kg, blind=TRUE))
vst_blind <- assay(vsd.kg)

#colnames(vst.kt) <- colnames(kt)
#already have the right names
vst_blind <- vst_blind - min(vst_blind) 
# substracts from all values the smaller value, so the smaller value is now 0

#' Validate the VST 
meanSdPlot(vst_blind[rowSums(vst_blind)>0,]) #mean variance stabilized between 0.5 and 1
meanSdPlot(log2(kg+1)) #without Variance Stabilizing Transformation: variation around 1
meanSdPlot(log2(counts(dds.kg,normalized=TRUE)+1)) #normalized read depth: variance is not stabilized ("curvy")
#VST validated: mean variance stabilized around 0.5


write.csv(vst_blind,"/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Tom-normalizedTo6_NP_DMSO-data-vst-blind.csv")

#save(kg, vst_blind, file="vst_blind.rda")
save(kg, vst_blind, file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/vst_blind.rda")

#' #### vst aware (calculated accordingly to samplingDates)

#' Calculate vst aware with DESeq2 taking the date as a parameter
#
#' #### Model-aware
system.time(vsd2 <- varianceStabilizingTransformation(dds.kg, blind=FALSE))
vst_aware <- assay(vsd2)

#colnames(vst.kt) <- colnames(kt)
#already have the right names
vst_aware <- vst_aware - min(vst_aware) 

# substracts from all values the smaller value, so the smaller value is now 0

#' Validate the VST 
meanSdPlot(vst_aware[rowSums(vst_aware)>0,]) #mean variance stabilized between 0 and 0.5
meanSdPlot(log2(kg+1)) #without Variance Stabilizing Transformation: variation around 1
meanSdPlot(log2(counts(dds.kg,normalized=TRUE)+1)) #normalized read depth: variance is not stabilized ("curvy")

# vst_aware validated: mean variance is way lower when vst is calculated "aware of the date factor" compare to before when vst was calculated "blind"
# stabilized around 0.2 instead of 0.5 when vst was blind


write.csv(vst_aware,"/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Tom-normalizedTo_6_NP_DMSO-data-vst-aware.csv")
#save(kg, vst_aware, file="vst_aware.rda")
#load("vst_aware.rda")

#' Differential expression
dds.kg <- DESeq(dds.kg)

#' which results
resultsNames(dds.kg)

#' ### Effect of AZD treatment in NP medium
# ------------------------------ Extract results from 6hrs of AZD vs DMSO in NP ---------------------------
res <- results(dds.kg,name="condition_6...NP...AZD_vs_6...NP...DMSO")

alpha=0.01

#' #### The log2 fold-change range is also rather tight.
plotMA(res)

#' #### Plot the log odds vs. log2 fold change
#' 
#' The volcano plot is very flat, cross-validating the results of the MA plot
volcanoPlot(res,alpha=alpha)

#' #### Plot the adjusted p-value histogram
#' 
#' Which is as expected, over-enriched for adjusted p-values of 1
hist(res$padj,breaks=seq(0,1,.01))

# select genes that are significant
alpha=0.01
l2fc = 0.5
sel <- ! is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= l2fc

# list of DE genes that are significant
AZDvsDMSO_in_NP <- rownames(res[sel,])

#' #### DE genes data export
#' 
write.csv(res,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/AZDvsDMSO_in_NP-DEgenes_all.csv")
write.csv(res[sel,],file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/AZDvsDMSO_in_NP-DEgenes__significant.csv")


#' ## VennDiagram of +/- AZD in different nutritional conditions
#'

SixH <- VennDiagram::calculate.overlap(list(AZDvsDMSO_in_NPS,
                                            AZDvsDMSO_in_PKS,
                                            AZDvsDMSO_in_NS,
                                            AZDvsDMSO_in_NP))
grid.newpage()
svg("Venn4.svg", width = 7, height = 7)
grid.draw(venn.diagram(list(AZDvsDMSO_in_NPS,
                            AZDvsDMSO_in_PKS,
                            AZDvsDMSO_in_NS,
                            AZDvsDMSO_in_NP),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("AZD in NPS","AZD in PKS","AZD in NS","AZD in NP")))
dev.off()
grid.draw(venn.diagram(list(AZDvsDMSO_in_NPS,
                            AZDvsDMSO_in_PKS,
                            AZDvsDMSO_in_NS,
                            AZDvsDMSO_in_NP),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("AZD in NPS","AZD in PKS","AZD in NS","AZD in NP")))







#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```


