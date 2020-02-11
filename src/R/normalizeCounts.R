#' ---
#' title: "Biological QA"
#' author: "Nicolas Delhomme & Iryna Shutava"
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
pal <- brewer.pal(8,"Dark2")

#' Register the default plot margin
mar <- par("mar")

#'#### Get the count data

#' # Raw data
#' ## Loading
#' Read the sample information
samples <- read.csv("~/Git/UPSCb/projects/arabidopsis-nutrition-TOR/doc/samples.csv")
#goi <- read.delim("~/Git/UPSCb/projects/aspen-FTL1-growth-cessation/doc/list_of_genes.txt", header = FALSE)
kg <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/analysis_Tom/nutrition-unormalised-gene-expression_data.csv", header = TRUE)

#kg <- read.csv("analysis/salmon/Domenique-unormalised-gene-expression_data.csv", header = TRUE)
#load("analysis/salmon/counts.rda")

#' Check for the genes that are never expressed
sel <- rowSums(kg) == 0 
# str(sel)= logical vector
# sums the counts for one gene in each samples
# show TRUE or FALSE for each gene if the gene is never expressed (has 0 counts in every time samples)

#' Check for the number of genes never expressed
sprintf("%s%% (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(kg),digits=1),
        sum(sel),
        nrow(kg))
#in our experiment the number of genes and transcript are the same (one transcript for one gene)

#' Display the samples mean raw counts distribution
#' 
#' i.e. the mean raw count of every 
#' gene across samples is calculated
#' and displayed on a log10 scale.
#'
plot(density(log10(rowMeans(kg))),col=pal[1],
     main="gene mean raw counts distribution",
     xlab="mean raw counts (log10)")
#density=likelihood, chance
# plot the likelihood that the mean of genes has this mean of raw counts 
# give an idea of the read depth, here the majority of gene maps ~30 reads (1.5 log10) 

#' Display all samples raw counts distribution
cols <- sample(pal,nlevels(samples$User.ID),replace = TRUE)
plot.multidensity(mclapply(1:ncol(kg),function(k){log10(kg)[,k]},mc.cores=16L),
                  col=cols[as.integer(samples$User.ID)],
                  legend.x="topright",
                  legend=levels(samples$User.ID),
                  legend.col=cols,
                  legend.lwd=2,
                  main="samples raw counts distribution",
                  xlab="per gene raw counts (log10)")

load("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/analysis_Tom/counts.rda")

#'##### Data normalisation 

#' ### vst blind

#' For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate
#'
#' Create the dds object, without giving any prior on the design


#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#---------------------------------- Analysis with T0 as a reference ----------------------------------------
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#conditions <- factor(rep(c("MIMIC172","OX156","OX172","T89"),each=6)) # to add the 
conditions <- factor(paste(samples$Timepoint,".",samples$Nutrition,".",samples$AZD))
 
conditions <- relevel(conditions,ref = "0 . 0 . 0")

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

#' ## Variance Stabilising Transformation
#' ### Blind
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

#write.csv(vst_blind,"analysis/kallisto/Domenique-normalized-data-vst-blind.csv")
write.csv(vst_blind,"/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Tom-normalized-data-vst-blind.csv")

#save(kg, vst_blind, file="vst_blind.rda")
save(kg, vst_blind, file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/vst_blind.rda")

#' ### vst aware (calculated accordingly to samplingDates)

#' Calculate vst aware with DESeq2 taking the date as a parameter
#' ### Model-aware
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

#write.csv(vst_aware,"analysis/kallisto/Domenique-normalized-data-vst-aware.csv")
write.csv(vst_aware,"/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Tom-normalized-data-vst-aware.csv")
#save(kg, vst_aware, file="vst_aware.rda")
#load("vst_aware.rda")

#' Differential expression
dds.kg <- DESeq(dds.kg)

#' which results
resultsNames(dds.kg)


#--------------------------------------------------------------------------------------------------------------------
#-------------------------------------------- Analysis of 24hrs ------------------------------------------------------

# ------------------------------ Extract results from 24hrs of Phosphorus starvation vs T0 ---------------------------
#' # Extract results from 6hrs of Phosphorus starvation vs T0
res <- results(dds.kg,name="condition_24...NS...DMSO_vs_0...0...0")

alpha=0.01

#' ## The log2 fold-change range is also rather tight.
plotMA(res)

#' ## Plot the log odds vs. log2 fold change
#' 
#' The volcano plot is very flat, cross-validating the results of the MA plot
volcanoPlot(res,alpha=alpha)

#' ## Plot the adjusted p-value histogram
#' 
#' Which is as expected, over-enriched for adjusted p-values of 1
hist(res$padj,breaks=seq(0,1,.01))

# select genes that are significant
alpha=0.01
l2fc = 0.5
sel <- ! is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= l2fc

# list of DE genes that are significant
twentyfour_NS_DMSO_vs_T0 <- rownames(res[sel,])

#' ## DE genes data export
#' 
#write.csv(res,file="analysis/salmon/MIMIC172_vs_T89-DEgenes_all.csv")
write.csv(res,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/24_NS_DMSO_vs_T0-DEgenes_all.csv")
#write.csv(res[sel,],file="analysis/salmon/MIMIC172_vs_T89-DEgenes_significant.csv")
write.csv(res[sel,],file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/24_NS_DMSO_vs_T0-DEgenes__significant.csv")

# ------------------------------ Extract results from 24hrs of No starvation vs T0 ---------------------------
#' # Extract results from 6hrs of Phosphorus starvation vs T0
res <- results(dds.kg,name="condition_24...NPS...DMSO_vs_0...0...0")

alpha=0.01

#' ## The log2 fold-change range is also rather tight.
plotMA(res)

#' ## Plot the log odds vs. log2 fold change
#' 
#' The volcano plot is very flat, cross-validating the results of the MA plot
volcanoPlot(res,alpha=alpha)

#' ## Plot the adjusted p-value histogram
#' 
#' Which is as expected, over-enriched for adjusted p-values of 1
hist(res$padj,breaks=seq(0,1,.01))

# select genes that are significant
alpha=0.01
l2fc = 0.5
sel <- ! is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= l2fc

# list of DE genes that are significant
twentyfour_NPS_DMSO_vs_T0 <- rownames(res[sel,])

#' ## DE genes data export
#' 
#write.csv(res,file="analysis/salmon/MIMIC172_vs_T89-DEgenes_all.csv")
write.csv(res,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/24_NPS_DMSO_vs_T0-DEgenes_all.csv")
#write.csv(res[sel,],file="analysis/salmon/MIMIC172_vs_T89-DEgenes_significant.csv")
write.csv(res[sel,],file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/24_NPS_DMSO_vs_T0-DEgenes__significant.csv")



# ------------------------------ Extract results from 24hrs of AZD and Phosphorus starvation vs T0 ---------------------------
#' # Extract results from 6hrs of Phosphorus starvation vs T0
res <- results(dds.kg,name="condition_24...NS...AZD_vs_0...0...0")

alpha=0.01

#' ## The log2 fold-change range is also rather tight.
plotMA(res)

#' ## Plot the log odds vs. log2 fold change
#' 
#' The volcano plot is very flat, cross-validating the results of the MA plot
volcanoPlot(res,alpha=alpha)

#' ## Plot the adjusted p-value histogram
#' 
#' Which is as expected, over-enriched for adjusted p-values of 1
hist(res$padj,breaks=seq(0,1,.01))

# select genes that are significant
alpha=0.01
l2fc = 0.5
sel <- ! is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= l2fc

# list of DE genes that are significant
twentyfour_NS_AZD_vs_T0 <- rownames(res[sel,])

#' ## DE genes data export
#' 
#write.csv(res,file="analysis/salmon/MIMIC172_vs_T89-DEgenes_all.csv")
write.csv(res,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/24_NS_AZD_vs_T0-DEgenes_all.csv")
#write.csv(res[sel,],file="analysis/salmon/MIMIC172_vs_T89-DEgenes_significant.csv")
write.csv(res[sel,],file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/24_NS_AZD_vs_T0-DEgenes__significant.csv")

# ------------------------------ Extract results from 24hrs of No starvation vs T0 ---------------------------
#' # Extract results from 6hrs of Phosphorus starvation vs T0
res <- results(dds.kg,name="condition_24...NPS...AZD_vs_0...0...0")

alpha=0.01

#' ## The log2 fold-change range is also rather tight.
plotMA(res)

#' ## Plot the log odds vs. log2 fold change
#' 
#' The volcano plot is very flat, cross-validating the results of the MA plot
volcanoPlot(res,alpha=alpha)

#' ## Plot the adjusted p-value histogram
#' 
#' Which is as expected, over-enriched for adjusted p-values of 1
hist(res$padj,breaks=seq(0,1,.01))

# select genes that are significant
alpha=0.01
l2fc = 0.5
sel <- ! is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= l2fc

# list of DE genes that are significant
twentyfour_NPS_AZD_vs_T0 <- rownames(res[sel,])

#' ## DE genes data export
#' 
#write.csv(res,file="analysis/salmon/MIMIC172_vs_T89-DEgenes_all.csv")
write.csv(res,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/24_NPS_AZD_vs_T0-DEgenes_all.csv")
#write.csv(res[sel,],file="analysis/salmon/MIMIC172_vs_T89-DEgenes_significant.csv")
write.csv(res[sel,],file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/24_NPS_AZD_vs_T0-DEgenes__significant.csv")



#' ## goi- list of genes of interests

#list_goi_6hPstarv <- res[rownames(res) %in% goi$V1,]

#res1 <- res[sel,]
#sel1 <- rownames(res[sel,]) %in% goi$V1
#list_goi_MIMIC172_significant <- res1[sel1,]

#rownames(list_goi_MIMIC172)
#rownames(list_goi_MIMIC172_significant)

#write.csv(res,file="analysis/salmon/MIMIC172_vs_T89-DEgenes_list_goi.csv")
#write.csv(list_goi_MIMIC172,file="/mnt/picea/home/ishutava/MIMIC172_vs_T89-DEgenes_list_goi.csv")
#write.csv(list_goi_MIMIC172_significant,file="/mnt/picea/home/ishutava/MIMIC172_vs_T89-DEgenes_list_goi_significant.csv")




#' ## DE genes data export
#' 
#' We have no DE genes for the OX172 mutants!
#' 
#write.csv(res,file="analysis/salmon/OX172_vs_T89-DEgenes_all.csv")
#write.csv(res,file="/mnt/picea/home/ishutava/OX172_vs_T89-DEgenes_all.csv")
#write.csv(res[sel,],file="analysis/salmon/OX172_vs_T89-DEgenes_significant.csv")
#write.csv(res[sel,],file="/mnt/picea/home/ishutava/OX172_vs_T89-DEgenes_significant.csv")



#' # VennDiagram of all comparison
#' 
#' 
twentyfourH <- VennDiagram::calculate.overlap(list(twentyfour_NS_DMSO_vs_T0, twentyfour_NS_AZD_vs_T0,twentyfour_NPS_AZD_vs_T0,twentyfour_NPS_DMSO_vs_T0))
grid.newpage()
grid.draw(venn.diagram(list(twentyfour_NS_DMSO_vs_T0, twentyfour_NS_AZD_vs_T0,twentyfour_NPS_AZD_vs_T0,twentyfour_NPS_DMSO_vs_T0),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("PStarv","PStarv_AZD","AZD","NoStarv")))




#--------------------------------------------------------------------------------------------------------------------
#-------------------------------------------- Analysis of 6hrs ------------------------------------------------------

# ------------------------------ Extract results from 6hrs of Phosphorus starvation vs T0 ---------------------------
#' # Extract results from 6hrs of Phosphorus starvation vs T0
res <- results(dds.kg,name="condition_6...NS...DMSO_vs_0...0...0")

alpha=0.01

#' ## The log2 fold-change range is also rather tight.
plotMA(res)

#' ## Plot the log odds vs. log2 fold change
#' 
#' The volcano plot is very flat, cross-validating the results of the MA plot
volcanoPlot(res,alpha=alpha)

#' ## Plot the adjusted p-value histogram
#' 
#' Which is as expected, over-enriched for adjusted p-values of 1
hist(res$padj,breaks=seq(0,1,.01))

# select genes that are significant
alpha=0.01
l2fc = 0.5
sel <- ! is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= l2fc

# list of DE genes that are significant
six_NS_DMSO_vs_T0 <- rownames(res[sel,])

#' ## DE genes data export
#' 
#write.csv(res,file="analysis/salmon/MIMIC172_vs_T89-DEgenes_all.csv")
write.csv(res,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/6_NS_DMSO_vs_T0-DEgenes_all.csv")
#write.csv(res[sel,],file="analysis/salmon/MIMIC172_vs_T89-DEgenes_significant.csv")
write.csv(res[sel,],file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/6_NS_DMSO_vs_T0-DEgenes__significant.csv")

# ------------------------------ Extract results from 6hrs of No starvation vs T0 ---------------------------
#' # Extract results from 6hrs of Phosphorus starvation vs T0
res <- results(dds.kg,name="condition_6...NPS...DMSO_vs_0...0...0")

alpha=0.01

#' ## The log2 fold-change range is also rather tight.
plotMA(res)

#' ## Plot the log odds vs. log2 fold change
#' 
#' The volcano plot is very flat, cross-validating the results of the MA plot
volcanoPlot(res,alpha=alpha)

#' ## Plot the adjusted p-value histogram
#' 
#' Which is as expected, over-enriched for adjusted p-values of 1
hist(res$padj,breaks=seq(0,1,.01))

# select genes that are significant
alpha=0.01
l2fc = 0.5
sel <- ! is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= l2fc

# list of DE genes that are significant
Six_NPS_DMSO_vs_T0 <- rownames(res[sel,])

#' ## DE genes data export
#' 
#write.csv(res,file="analysis/salmon/MIMIC172_vs_T89-DEgenes_all.csv")
write.csv(res,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/6_NPS_DMSO_vs_T0-DEgenes_all.csv")
#write.csv(res[sel,],file="analysis/salmon/MIMIC172_vs_T89-DEgenes_significant.csv")
write.csv(res[sel,],file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/6_NPS_DMSO_vs_T0-DEgenes__significant.csv")



# ------------------------------ Extract results from 6hrs of AZD and Phosphorus starvation vs T0 ---------------------------
#' # Extract results from 6hrs of Phosphorus starvation vs T0
res <- results(dds.kg,name="condition_6...NS...AZD_vs_0...0...0")

alpha=0.01

#' ## The log2 fold-change range is also rather tight.
plotMA(res)

#' ## Plot the log odds vs. log2 fold change
#' 
#' The volcano plot is very flat, cross-validating the results of the MA plot
volcanoPlot(res,alpha=alpha)

#' ## Plot the adjusted p-value histogram
#' 
#' Which is as expected, over-enriched for adjusted p-values of 1
hist(res$padj,breaks=seq(0,1,.01))

# select genes that are significant
alpha=0.01
l2fc = 0.5
sel <- ! is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= l2fc

# list of DE genes that are significant
six_NS_AZD_vs_T0 <- rownames(res[sel,])

#' ## DE genes data export
#' 
#write.csv(res,file="analysis/salmon/MIMIC172_vs_T89-DEgenes_all.csv")
write.csv(res,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/6_NS_AZD_vs_T0-DEgenes_all.csv")
#write.csv(res[sel,],file="analysis/salmon/MIMIC172_vs_T89-DEgenes_significant.csv")
write.csv(res[sel,],file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/6_NS_AZD_vs_T0-DEgenes__significant.csv")

# ------------------------------ Extract results from 6hrs of No starvation vs T0 ---------------------------
#' # Extract results from 6hrs of Phosphorus starvation vs T0
res <- results(dds.kg,name="condition_6...NPS...AZD_vs_0...0...0")

alpha=0.01

#' ## The log2 fold-change range is also rather tight.
plotMA(res)

#' ## Plot the log odds vs. log2 fold change
#' 
#' The volcano plot is very flat, cross-validating the results of the MA plot
volcanoPlot(res,alpha=alpha)

#' ## Plot the adjusted p-value histogram
#' 
#' Which is as expected, over-enriched for adjusted p-values of 1
hist(res$padj,breaks=seq(0,1,.01))

# select genes that are significant
alpha=0.01
l2fc = 0.5
sel <- ! is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= l2fc

# list of DE genes that are significant
Six_NPS_AZD_vs_T0 <- rownames(res[sel,])

#' ## DE genes data export
#' 
#write.csv(res,file="analysis/salmon/MIMIC172_vs_T89-DEgenes_all.csv")
write.csv(res,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/6_NPS_AZD_vs_T0-DEgenes_all.csv")
#write.csv(res[sel,],file="analysis/salmon/MIMIC172_vs_T89-DEgenes_significant.csv")
write.csv(res[sel,],file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/6_NPS_AZD_vs_T0-DEgenes__significant.csv")



#' ## goi- list of genes of interests

#list_goi_6hPstarv <- res[rownames(res) %in% goi$V1,]

#res1 <- res[sel,]
#sel1 <- rownames(res[sel,]) %in% goi$V1
#list_goi_MIMIC172_significant <- res1[sel1,]

#rownames(list_goi_MIMIC172)
#rownames(list_goi_MIMIC172_significant)

#write.csv(res,file="analysis/salmon/MIMIC172_vs_T89-DEgenes_list_goi.csv")
#write.csv(list_goi_MIMIC172,file="/mnt/picea/home/ishutava/MIMIC172_vs_T89-DEgenes_list_goi.csv")
#write.csv(list_goi_MIMIC172_significant,file="/mnt/picea/home/ishutava/MIMIC172_vs_T89-DEgenes_list_goi_significant.csv")




#' ## DE genes data export
#' 
#' We have no DE genes for the OX172 mutants!
#' 
#write.csv(res,file="analysis/salmon/OX172_vs_T89-DEgenes_all.csv")
#write.csv(res,file="/mnt/picea/home/ishutava/OX172_vs_T89-DEgenes_all.csv")
#write.csv(res[sel,],file="analysis/salmon/OX172_vs_T89-DEgenes_significant.csv")
#write.csv(res[sel,],file="/mnt/picea/home/ishutava/OX172_vs_T89-DEgenes_significant.csv")

#' # VennDiagram of all comparison
#'
#'
mimic172_vs_ox156 <- VennDiagram::calculate.overlap(list(mimic172,ox156))
grid.newpage()
draw.pairwise.venn(length(mimic172_vs_ox156$a1), length(mimic172_vs_ox156$a2), length(mimic172_vs_ox156$a3),
                   c("MIMIC172", "OX156"), fill = rainbow(2), alpha = 0.6)

#' # VennDiagram of all comparison include gene list of interest 
#' 
#' 
SixH <- VennDiagram::calculate.overlap(list(six_NS_DMSO_vs_T0, six_NS_AZD_vs_T0,Six_NPS_AZD_vs_T0,Six_NPS_DMSO_vs_T0))
grid.newpage()
grid.draw(venn.diagram(list(six_NS_DMSO_vs_T0, six_NS_AZD_vs_T0,Six_NPS_AZD_vs_T0,Six_NPS_DMSO_vs_T0),
                       filename=NULL,
                       col=pal[1:4],
                       category.names=c("PStarv","PStarv_AZD","AZD","NoStarv")))












#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#----------------------------- Analysis with NPS_DMSO-6h as a reference ------------------------------------
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#conditions <- factor(rep(c("MIMIC172","OX156","OX172","T89"),each=6)) # to add the 
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

#' ## Variance Stabilising Transformation
#' ### Blind
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

#write.csv(vst_blind,"analysis/kallisto/Domenique-normalized-data-vst-blind.csv")
write.csv(vst_blind,"/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Tom-normalized_to_6h_NPS_DMSO-data-vst-blind.csv")

#save(kg, vst_blind, file="vst_blind.rda")
save(kg, vst_blind, file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/vst_blind.rda")

#' ### vst aware (calculated accordingly to samplingDates)

#' Calculate vst aware with DESeq2 taking the date as a parameter
#' ### Model-aware
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

#write.csv(vst_aware,"analysis/kallisto/Domenique-normalized-data-vst-aware.csv")
write.csv(vst_aware,"/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Tom-normalized-data-6_NPS_DMSO-vst-aware.csv")
#save(kg, vst_aware, file="vst_aware.rda")
#load("vst_aware.rda")

#' Differential expression
dds.kg <- DESeq(dds.kg)

#' which results
resultsNames(dds.kg)

# ------------------------------ Extract results from 6hrs of Phosphorus starvation vs six_NPS_DMSO ---------------------------
#' # Extract results from 6hrs of Phosphorus starvation vs T0
res <- results(dds.kg,name="condition_6...NS...DMSO_vs_6...NPS...DMSO")

alpha=0.01

#' ## The log2 fold-change range is also rather tight.
plotMA(res)

#' ## Plot the log odds vs. log2 fold change
#' 
#' The volcano plot is very flat, cross-validating the results of the MA plot
volcanoPlot(res,alpha=alpha)

#' ## Plot the adjusted p-value histogram
#' 
#' Which is as expected, over-enriched for adjusted p-values of 1
hist(res$padj,breaks=seq(0,1,.01))

# select genes that are significant
alpha=0.01
l2fc = 0.5
sel <- ! is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= l2fc

# list of DE genes that are significant
six_NS_DMSO_vs_six_NPS_DMSO <- rownames(res[sel,])

#' ## DE genes data export
#' 
#write.csv(res,file="analysis/salmon/MIMIC172_vs_T89-DEgenes_all.csv")
write.csv(res,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/6_NS_DMSO_vs_six_NPS_DMSO-DEgenes_all.csv")
#write.csv(res[sel,],file="analysis/salmon/MIMIC172_vs_T89-DEgenes_significant.csv")
write.csv(res[sel,],file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/6_NS_DMSO_vs_six_NPS_DMSO-DEgenes__significant.csv")



# ------------------------------ Extract results from 6hrs of AZD and Phosphorus starvation vs 6_NPS_DMSO ---------------------------
#' # Extract results from 6hrs of Phosphorus starvation vs T0
res <- results(dds.kg,name="condition_6...NS...AZD_vs_6...NPS...DMSO")

alpha=0.01

#' ## The log2 fold-change range is also rather tight.
plotMA(res)

#' ## Plot the log odds vs. log2 fold change
#' 
#' The volcano plot is very flat, cross-validating the results of the MA plot
volcanoPlot(res,alpha=alpha)

#' ## Plot the adjusted p-value histogram
#' 
#' Which is as expected, over-enriched for adjusted p-values of 1
hist(res$padj,breaks=seq(0,1,.01))

# select genes that are significant
alpha=0.01
l2fc = 0.5
sel <- ! is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= l2fc

# list of DE genes that are significant
six_NS_AZD_vs_six_NPS_DMSO <- rownames(res[sel,])

#' ## DE genes data export
#' 
#write.csv(res,file="analysis/salmon/MIMIC172_vs_T89-DEgenes_all.csv")
write.csv(res,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/6_NS_AZD_vs_six_NPS_DMSO-DEgenes_all.csv")
#write.csv(res[sel,],file="analysis/salmon/MIMIC172_vs_T89-DEgenes_significant.csv")
write.csv(res[sel,],file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/6_NS_AZD_vs_six_NPS_DMSO-DEgenes__significant.csv")

# ------------------------------ Extract results from 6hrs of No starvation vs 6_NPS_DMSO ---------------------------
#' # Extract results from 6hrs of Phosphorus starvation vs T0
res <- results(dds.kg,name="condition_6...NPS...AZD_vs_6...NPS...DMSO")

alpha=0.01

#' ## The log2 fold-change range is also rather tight.
plotMA(res)

#' ## Plot the log odds vs. log2 fold change
#' 
#' The volcano plot is very flat, cross-validating the results of the MA plot
volcanoPlot(res,alpha=alpha)

#' ## Plot the adjusted p-value histogram
#' 
#' Which is as expected, over-enriched for adjusted p-values of 1
hist(res$padj,breaks=seq(0,1,.01))

# select genes that are significant
alpha=0.01
l2fc = 0.5
sel <- ! is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= l2fc

# list of DE genes that are significant
Six_NPS_AZD_vs_six_NPS_DMSO <- rownames(res[sel,])

#' ## DE genes data export
#' 
#write.csv(res,file="analysis/salmon/MIMIC172_vs_T89-DEgenes_all.csv")
write.csv(res,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/6_NPS_AZD_vs_six_NPS_DMSO-DEgenes_all.csv")
#write.csv(res[sel,],file="analysis/salmon/MIMIC172_vs_T89-DEgenes_significant.csv")
write.csv(res[sel,],file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/6_NPS_AZD_vs_six_NPS_DMSO-DEgenes__significant.csv")



#' ## goi- list of genes of interests

#list_goi_6hPstarv <- res[rownames(res) %in% goi$V1,]

#res1 <- res[sel,]
#sel1 <- rownames(res[sel,]) %in% goi$V1
#list_goi_MIMIC172_significant <- res1[sel1,]

#rownames(list_goi_MIMIC172)
#rownames(list_goi_MIMIC172_significant)

#write.csv(res,file="analysis/salmon/MIMIC172_vs_T89-DEgenes_list_goi.csv")
#write.csv(list_goi_MIMIC172,file="/mnt/picea/home/ishutava/MIMIC172_vs_T89-DEgenes_list_goi.csv")
#write.csv(list_goi_MIMIC172_significant,file="/mnt/picea/home/ishutava/MIMIC172_vs_T89-DEgenes_list_goi_significant.csv")




#' ## DE genes data export
#' 
#' We have no DE genes for the OX172 mutants!
#' 
#write.csv(res,file="analysis/salmon/OX172_vs_T89-DEgenes_all.csv")
#write.csv(res,file="/mnt/picea/home/ishutava/OX172_vs_T89-DEgenes_all.csv")
#write.csv(res[sel,],file="analysis/salmon/OX172_vs_T89-DEgenes_significant.csv")
#write.csv(res[sel,],file="/mnt/picea/home/ishutava/OX172_vs_T89-DEgenes_significant.csv")

#' # VennDiagram of all comparison
#'
#'
mimic172_vs_ox156 <- VennDiagram::calculate.overlap(list(mimic172,ox156))
grid.newpage()
draw.pairwise.venn(length(mimic172_vs_ox156$a1), length(mimic172_vs_ox156$a2), length(mimic172_vs_ox156$a3),
                   c("MIMIC172", "OX156"), fill = rainbow(2), alpha = 0.6)

#' # VennDiagram of all comparison include gene list of interest 
#' 
#' 
SixH <- VennDiagram::calculate.overlap(list(six_NS_DMSO_vs_six_NPS_DMSO, six_NS_AZD_vs_six_NPS_DMSO,Six_NPS_AZD_vs_six_NPS_DMSO))
grid.newpage()
grid.draw(venn.diagram(list(six_NS_DMSO_vs_six_NPS_DMSO, six_NS_AZD_vs_six_NPS_DMSO,Six_NPS_AZD_vs_six_NPS_DMSO),
                       filename=NULL,
                       col=pal[1:3],
                       category.names=c("PStarv","PStarv_AZD","AZD")))





#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#----------------------------- Analysis with NPS_DMSO-24h as a reference ------------------------------------
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#conditions <- factor(rep(c("MIMIC172","OX156","OX172","T89"),each=6)) # to add the 
conditions <- factor(paste(samples$Timepoint,".",samples$Nutrition,".",samples$AZD))

conditions <- relevel(conditions,ref = "24 . NPS . DMSO")

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

#' ## Variance Stabilising Transformation
#' ### Blind
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

#write.csv(vst_blind,"analysis/kallisto/Domenique-normalized-data-vst-blind.csv")
write.csv(vst_blind,"/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Tom-normalized_to_24h_NPS_DMSO-data-vst-blind.csv")

#save(kg, vst_blind, file="vst_blind.rda")
save(kg, vst_blind, file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/vst_blind.rda")

#' ### vst aware (calculated accordingly to samplingDates)

#' Calculate vst aware with DESeq2 taking the date as a parameter
#' ### Model-aware
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

#write.csv(vst_aware,"analysis/kallisto/Domenique-normalized-data-vst-aware.csv")
write.csv(vst_aware,"/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Tom-normalized-data-24_NPS_DMSO-vst-aware.csv")
#save(kg, vst_aware, file="vst_aware.rda")
#load("vst_aware.rda")

#' Differential expression
dds.kg <- DESeq(dds.kg)

#' which results
resultsNames(dds.kg)

# ------------------------------ Extract results from 24_NS_DMSO vs 24_NPS_DMSO ---------------------------
#' # Extract results from 24hrs of Phosphorus starvation vs 24_NPS_DMSO
res <- results(dds.kg,name="condition_24...NS...DMSO_vs_24...NPS...DMSO")

alpha=0.01

#' ## The log2 fold-change range is also rather tight.
plotMA(res)

#' ## Plot the log odds vs. log2 fold change
#' 
#' The volcano plot is very flat, cross-validating the results of the MA plot
volcanoPlot(res,alpha=alpha)

#' ## Plot the adjusted p-value histogram
#' 
#' Which is as expected, over-enriched for adjusted p-values of 1
hist(res$padj,breaks=seq(0,1,.01))

# select genes that are significant
alpha=0.01
l2fc = 0.5
sel <- ! is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= l2fc

# list of DE genes that are significant
twentyfour_NS_DMSO_vs_twentyfour_NPS_DMSO <- rownames(res[sel,])

#' ## DE genes data export
#' 
#write.csv(res,file="analysis/salmon/MIMIC172_vs_T89-DEgenes_all.csv")
write.csv(res,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/24_NS_DMSO_vs_24_NPS_DMSO-DEgenes_all.csv")
#write.csv(res[sel,],file="analysis/salmon/MIMIC172_vs_T89-DEgenes_significant.csv")
write.csv(res[sel,],file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/24_NS_DMSO_vs_24_NPS_DMSO-DEgenes__significant.csv")



# ------------------------------ Extract results from 24hrs of AZD and Phosphorus starvation vs 24_NPS_DMSO ---------------------------
#' # Extract results from 24_NS_AZD vs 24_NPS_DMSO
res <- results(dds.kg,name="condition_24...NS...AZD_vs_24...NPS...DMSO")

alpha=0.01

#' ## The log2 fold-change range is also rather tight.
plotMA(res)

#' ## Plot the log odds vs. log2 fold change
#' 
#' The volcano plot is very flat, cross-validating the results of the MA plot
volcanoPlot(res,alpha=alpha)

#' ## Plot the adjusted p-value histogram
#' 
#' Which is as expected, over-enriched for adjusted p-values of 1
hist(res$padj,breaks=seq(0,1,.01))

# select genes that are significant
alpha=0.01
l2fc = 0.5
sel <- ! is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= l2fc

# list of DE genes that are significant
twentyfour_NS_AZD_vs_twentyfour_NPS_DMSO <- rownames(res[sel,])

#' ## DE genes data export
#' 
#write.csv(res,file="analysis/salmon/MIMIC172_vs_T89-DEgenes_all.csv")
write.csv(res,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/24_NS_AZD_vs_24_NPS_DMSO-DEgenes_all.csv")
#write.csv(res[sel,],file="analysis/salmon/MIMIC172_vs_T89-DEgenes_significant.csv")
write.csv(res[sel,],file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/24_NS_AZD_vs_24_NPS_DMSO-DEgenes__significant.csv")

# ------------------------------ Extract results from 24hrs of AZD vs 24_NPS_DMSO ---------------------------
#' # Extract results from 24hrs of Phosphorus starvation vs 24_NPS_DMSO
res <- results(dds.kg,name="condition_24...NPS...AZD_vs_24...NPS...DMSO")

alpha=0.01

#' ## The log2 fold-change range is also rather tight.
plotMA(res)

#' ## Plot the log odds vs. log2 fold change
#' 
#' The volcano plot is very flat, cross-validating the results of the MA plot
volcanoPlot(res,alpha=alpha)

#' ## Plot the adjusted p-value histogram
#' 
#' Which is as expected, over-enriched for adjusted p-values of 1
hist(res$padj,breaks=seq(0,1,.01))

# select genes that are significant
alpha=0.01
l2fc = 0.5
sel <- ! is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= l2fc

# list of DE genes that are significant
twentyfour_NPS_AZD_vs_twentyfour_NPS_DMSO <- rownames(res[sel,])

#' ## DE genes data export
#' 
#write.csv(res,file="analysis/salmon/MIMIC172_vs_T89-DEgenes_all.csv")
write.csv(res,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/24_NPS_AZD_vs_24_NPS_DMSO-DEgenes_all.csv")
#write.csv(res[sel,],file="analysis/salmon/MIMIC172_vs_T89-DEgenes_significant.csv")
write.csv(res[sel,],file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/24_NPS_AZD_vs_24_NPS_DMSO-DEgenes__significant.csv")



#' ## goi- list of genes of interests

#list_goi_6hPstarv <- res[rownames(res) %in% goi$V1,]

#res1 <- res[sel,]
#sel1 <- rownames(res[sel,]) %in% goi$V1
#list_goi_MIMIC172_significant <- res1[sel1,]

#rownames(list_goi_MIMIC172)
#rownames(list_goi_MIMIC172_significant)

#write.csv(res,file="analysis/salmon/MIMIC172_vs_T89-DEgenes_list_goi.csv")
#write.csv(list_goi_MIMIC172,file="/mnt/picea/home/ishutava/MIMIC172_vs_T89-DEgenes_list_goi.csv")
#write.csv(list_goi_MIMIC172_significant,file="/mnt/picea/home/ishutava/MIMIC172_vs_T89-DEgenes_list_goi_significant.csv")




#' ## DE genes data export
#' 
#' We have no DE genes for the OX172 mutants!
#' 
#write.csv(res,file="analysis/salmon/OX172_vs_T89-DEgenes_all.csv")
#write.csv(res,file="/mnt/picea/home/ishutava/OX172_vs_T89-DEgenes_all.csv")
#write.csv(res[sel,],file="analysis/salmon/OX172_vs_T89-DEgenes_significant.csv")
#write.csv(res[sel,],file="/mnt/picea/home/ishutava/OX172_vs_T89-DEgenes_significant.csv")

#' # VennDiagram of all comparison

#' 
#' 
twentyfourH <- VennDiagram::calculate.overlap(list(twentyfour_NS_DMSO_vs_twentyfour_NPS_DMSO, twentyfour_NS_AZD_vs_twentyfour_NPS_DMSO,twentyfour_NPS_AZD_vs_twentyfour_NPS_DMSO))
grid.newpage()
grid.draw(venn.diagram(list(twentyfour_NS_DMSO_vs_twentyfour_NPS_DMSO, twentyfour_NS_AZD_vs_twentyfour_NPS_DMSO,twentyfour_NPS_AZD_vs_twentyfour_NPS_DMSO),
                       filename=NULL,
                       col=pal[1:3],
                       category.names=c("PStarv","PStarv_AZD","AZD")))

















#---------
##' ### Filter data for expressed enough genes 

##' ## Select the genes that are expressed 

#sels <- sapply(1:13,function(i){
#    featureSelect(vst_blind,conditions = factor(paste0(samples$SciLifeID,samples$User.ID)),
#                  exp=i)})


#plot(colSums(sels),type="l",xlab="vst cutoff",
#     main="number of genes selected at cutoff",ylab="number of genes")

#sel <- sels[,2]



##' Hierarchical clustering of the data
#plot(hclust(dist(t(vst_aware[]))),labels=samples$User.ID)


##' Create a heatmap
#hpal <- colorRampPalette(c("blue","white","red"))(100)

#heatmap.2(vst_aware[],trace="none",col=hpal)

#s.vst <- t(scale(t(vst_aware)))

#library(hyperSpec)

#heatmap.2(s.vst[],distfun = pearson.dist,
#          hclustfun = function(X){hclust(X,method="ward.D")},
#          trace="none",col=hpal,labRow = FALSE,
#          labCol = paste(samples$Genotype,samples$Time,sep="-"))


#hc <- hclust(pearson.dist(s.vst[]),method = "ward.D")

#tc <- cutree(hc,k=4)

#tc

#nams <- split(names(tc),tc) 
#library(IRanges)
#barplot(elementNROWS(nams))

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
