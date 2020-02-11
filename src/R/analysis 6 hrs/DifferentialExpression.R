#' ---
#' title: "Differential Expression analyses"
#' author: "Thomas Dobrenel and Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' * Working dir
setwd("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR")
#' ```

#' * Libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(hyperSpec))
suppressPackageStartupMessages(library(LSD))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyverse))

#' * Helpers
suppressPackageStartupMessages(source("~/Git/UPSCb/src/R/gopher.R"))

#' * Graphics
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")
suppressPackageStartupMessages(source("~/Git/UPSCb/src/R/volcanoPlot.R"))

#' * Cutoff
#' As of Schurch et al., RNA, 2016
lfc <- 0.5
FDR <- 0.01

#' * Raw data
load("analysis_Tom/counts.rda")

#' * Selection of the validated samples
samples %<>% filter(!grepl("P11554_1",SciLifeID)) %>% 
    filter(Timepoint!=24) %>% mutate(Timepoint,Timepoint=factor(paste0("T",Timepoint))) %>% 
    mutate(Nutrition,Nutrition=relevel(Nutrition,"NPS")) %>% 
    mutate(AZD,AZD=relevel(AZD,"DMSO"))

kg <- kg[,match(samples$SciLifeID,colnames(kg))]

#' * Validation
stopifnot(all(samples$SciLifeID == colnames(kg)))
stopifnot(all(samples$Batch == "B2"))

#' # Differential expression
#' ## Nutrition and Treatment at T6
sel_T6 <- samples$Timepoint == "T6"
suppressMessages(dds_T6 <- DESeqDataSetFromMatrix(
    countData = kg[,sel_T6],
    colData = samples[sel_T6,],
    design = ~ Nutrition * AZD ))

dds_T6 <- DESeq(dds_T6)

#' ## Variance Stabilising Transformation
#' Perform a VST for plotting
vst_T6 <- varianceStabilizingTransformation(dds_T6,blind=FALSE)
vsd_T6 <- assay(vst_T6)
vsd_T6 <- vsd_T6 - min(vsd_T6)

#' Calculate fold change
colAvg_T6 <- sapply(split.data.frame(t(vsd_T6),f = droplevels(samples$Conditions[sel_T6])),colMeans)
fc_T6 <- colAvg_T6[,- match("6_NPS_DMSO",colnames(colAvg_T6))] - colAvg_T6[,"6_NPS_DMSO"]
fc2_T6 <- vsd_T6 - colAvg_T6[,"6_NPS_DMSO"]
range(fc_T6)
fc_T6[fc_T6>2] <- 2 ; fc_T6[fc_T6 < -2] <- -2
fc2_T6[fc2_T6>2] <- 2 ; fc2_T6[fc2_T6 < -2] <- -2
#range(fc_T6)

#' ### Contrasts on NS
#' The contrast by default is the first one (not Intercept)
resultsNames(dds_T6)

#' #### Nutrition effect
res_T6 <- results(dds_T6,name = "Nutrition_NS_vs_NPS")

# res <- results(dds,contrast = c(Nutrition","NP","NS"))
# res <- results(dds,contrast = list("Nutrition_NP_vs_NPS","Nutrition_NS_vs_NPS"))
# res <- results(dds,contrast = c(0,1,-1,0,0,0,0,0,0) # NP vs. NS

cutoffs <- abs(res_T6$log2FoldChange) >= lfc & ! is.na(res_T6$padj) & res_T6$padj <= FDR

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs)))

#' #### MA plot
DESeq2::plotMA(res_T6)

#' #### Volcano plot 
volcanoPlot(res_T6)

#' #### Heatmap
#' ##### Z-score
hm <- heatmap.2(t(scale(t(vsd_T6[cutoffs,]))),col=hpal,
          Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
          distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel_T6],
          hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))

#' ##### Fold-change
hmfc <- heatmap.2(fc_T6[cutoffs,],
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                col=hpal, margins = c(6,5),cexCol = 0.8)

hmfc <- heatmap.2(fc_T6[cutoffs,],
                  trace = "none",labRow = FALSE,
                  col=hpal, margins = c(6,5),cexCol = 0.8)

hmfc <- heatmap.2(fc_T6[cutoffs,],
                  trace = "none",labRow = FALSE,
                  col=hpal, margins = c(6,5),cexCol = 0.8,
                  distfun = pearson.dist)

hmfc <- heatmap.2(fc2_T6[cutoffs,],
                  trace = "none",labRow = FALSE,
                  col=hpal, margins = c(6,5),
                  labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel_T6],cexCol = 0.8)

#' #### hclust
tree <- cutree(as.hclust(hm$rowDendrogram),6)
c1 <- names(tree)[tree==1]
message(sprintf("There are %s genes in the cluster1",length(c1)))
heatmap.2(t(scale(t(vsd_T6[c1,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel_T6],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))
c2 <- names(tree)[tree==2]
message(sprintf("There are %s genes in the cluster2",length(c2)))
heatmap.2(t(scale(t(vsd_T6[c2,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel_T6],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))
c3 <- names(tree)[tree==3]
message(sprintf("There are %s genes in the cluster3",length(c3)))
heatmap.2(t(scale(t(vsd_T6[c3,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel_T6],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))
c4 <- names(tree)[tree==4]
message(sprintf("There are %s genes in the cluster4",length(c4)))
heatmap.2(t(scale(t(vsd_T6[c4,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel_T6],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))
c5 <- names(tree)[tree==5]
message(sprintf("There are %s genes in the cluster5",length(c5)))
heatmap.2(t(scale(t(vsd_T6[c5,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel_T6],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))
c6 <- names(tree)[tree==6]
message(sprintf("There are %s genes in the cluster6",length(c6)))
heatmap.2(t(scale(t(vsd_T6[c6,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel_T6],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))

goi <- names(tree)[tree==3]

ggplot(melt(vsd_T6[goi,]),aes(x=Var2,y=value,group=Var2)) + geom_boxplot() +
    theme(axis.text.x = element_text(angle=90)) +
    scale_x_discrete(label=paste(samples$Nutrition,samples$AZD,sep="-")[sel_T6],element_blank()) +
    scale_y_continuous("expression (VST)")

goi2 <- names(tree)[tree==2]

ggplot(melt(vsd_T6[goi,]),aes(x=Var2,y=value,group=Var2)) + geom_boxplot() +
    theme(axis.text.x = element_text(angle=90)) +
    scale_x_discrete(label=paste(samples$Nutrition,samples$AZD,sep="-")[sel_T6],element_blank()) +
    scale_y_continuous("expression (VST)")



#' ##### Analysis of GOI
#' ###### Enrichment
#enr <- gopher(goi,background=rownames(vsd_T6)[rowSums(vsd_T6)>0],url="athaliana")

#enr$go$name

#enr$kegg

#enr$pfam$name

#' ###### Expression profile
#' * Differential expression in T0 and T6
sel_T0 <- samples$Timepoint %in% c("T6","T0")
dds_T0 <- DESeqDataSetFromMatrix(
    countData = kg[,sel_T0],
    colData = samples[sel_T0,],
    design = ~ Conditions )

dds_T0 <- DESeq(dds_T0)

#' * Variance Stabilising Transformation and Perform a VST for plotting
vst_T0 <- varianceStabilizingTransformation(dds_T0,blind=FALSE)
vsd_T0 <- assay(vst_T0)
vsd_T0 <- vsd_T0 - min(vsd_T0)

#' * Calculate fold change
colAvg_T0 <- sapply(split.data.frame(t(vsd_T0),f = droplevels(samples$Conditions[sel_T0])),colMeans)
fc2_T0 <- vsd_T0- colAvg_T0[,"0_T0_0"]
range(fc2_T0)
fc2_T0[fc2_T0>2] <- 2 ; fc2_T0[fc2_T0 < -2] <- -2

#' * Contrasts on T0
#' The contrast by default is the first one (not Intercept)
resultsNames(dds_T0)

message(sprintf("There are %s genes of interest",length(goi)))

#' #### Heatmap
#' ##### Z-score
hm <- heatmap.2(t(scale(t(vsd_T0[goi,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel_T0],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))

#' ##### Fold-change
hmfc <- heatmap.2(fc2_T0[goi,],
                  Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                  col=hpal, margins = c(6,5),cexCol = 0.8)

hmfc <- heatmap.2(fc2_T0[goi,],
                  trace = "none",labRow = FALSE,
                  col=hpal, margins = c(6,5),cexCol = 0.8)

hmfc <- heatmap.2(fc2_T0[goi,],
                  trace = "none",labRow = FALSE,
                  col=hpal, margins = c(6,5),cexCol = 0.8,
                  distfun = pearson.dist)

hmfc <- heatmap.2(fc2_T0[goi,],
                  trace = "none",labRow = FALSE,
                  col=hpal, margins = c(6,5),
                  labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel_T0],cexCol = 0.8)






#' # Expression level of the TOR complex
#' * Preparation of the complex member list
sel1 <- c("AT1G50030","AT3G18140","AT2G22040","AT3G08850","AT5G01770")
expr <- vsd_T6[match(rownames(vsd_T6),sel1),]
expr <- na.omit(expr)
plot(expr[1,])
plot(expr[3,])
AvgTOR <- sapply(split.data.frame(t(expr),f = droplevels(samples$Conditions[sel_T6])),colMeans)
barplot(AvgTOR,beside=T, cex.names=0.7,legend.text = sel1)
plot(AvgTOR[1,])

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```

