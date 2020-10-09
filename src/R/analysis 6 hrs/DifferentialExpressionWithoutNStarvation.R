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
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(VennDiagram))


#' * Helpers
suppressPackageStartupMessages(source("~/Git/UPSCb/UPSCb-common/src/R/gopher.R"))

#' * Graphics
hpal <- colorRampPalette(c("blue","white","red"))(100)
pal <- brewer.pal(9,"Paired") #require package RColorBrewer
mar <- par("mar")
suppressPackageStartupMessages(source("~/Git/UPSCb/UPSCb-common/src/R/volcanoPlot.R"))

#' * Cutoff
#' As of Schurch et al., RNA, 2016
lfc <- 0.5
FDR <- 0.01

#' * Raw data
load("analysis_Tom/counts.rda")

#' * Selection of the validated samples
samples %<>% filter(!grepl("P11554_1",SciLifeID)) %>% 
    filter(Timepoint!=24) %>% 
    filter(Nutrition!="PKS") %>% 
    mutate(Timepoint,Timepoint=factor(paste0("T",Timepoint))) %>% 
    mutate(Nutrition,Nutrition=relevel(Nutrition,"NPS")) %>% 
    mutate(AZD,AZD=relevel(AZD,"DMSO"))

kg <- kg[,match(samples$SciLifeID,colnames(kg))]

#' * Validation
stopifnot(all(samples$SciLifeID == colnames(kg)))
stopifnot(all(samples$Batch == "B2"))

#' # Differential expression based on the nutrition and treatment at T6
#' ## Filtration of samples based on timepoint
sel_T6 <- samples$Timepoint == "T6"
suppressMessages(dds_T6 <- DESeqDataSetFromMatrix(
    countData = kg[,sel_T6],
    colData = samples[sel_T6,],
    design = ~ Nutrition * AZD ))

#' ## Differential expression analysis
dds_T6 <- DESeq(dds_T6)

#' ## Variance Stabilising Transformation
#' ### Perform a Variance Stabilizing Transformation for plotting
vst_T6 <- varianceStabilizingTransformation(dds_T6,blind=FALSE)
vsd_T6 <- assay(vst_T6)
vsd_T6 <- vsd_T6 - min(vsd_T6)

#' ### Calculation of fold change compared to NPS_DMSO at 6 hrs
colAvg_T6 <- sapply(split.data.frame(t(vsd_T6),f = droplevels(samples$Conditions[sel_T6])),colMeans)
fc_T6 <- colAvg_T6[,- match("6_NPS_DMSO",colnames(colAvg_T6))] - colAvg_T6[,"6_NPS_DMSO"]
fc2_T6 <- vsd_T6 - colAvg_T6[,"6_NPS_DMSO"]
range(fc_T6)
fc_T6[fc_T6>2] <- 2 ; fc_T6[fc_T6 < -2] <- -2
fc2_T6[fc2_T6>2] <- 2 ; fc2_T6[fc2_T6 < -2] <- -2
#range(fc_T6)

#' ## Contrasts to NPS_DMSO
#' The contrast by default is the first one (not Intercept)
resultsNames(dds_T6)

#' ### Nutrition effect of carbon starvation
#' #### Extraction of the results from the DESeq analysis
res_Suc <- results(dds_T6,name = "Nutrition_NP_vs_NPS")

# res <- results(dds,contrast = c(Nutrition","NP","NS"))
# res <- results(dds,contrast = list("Nutrition_NP_vs_NPS","Nutrition_NS_vs_NPS"))
# res <- results(dds,contrast = c(0,1,-1,0,0,0,0,0,0) # NP vs. NS

cutoffs_Suc <- abs(res_Suc$log2FoldChange) >= lfc & ! is.na(res_Suc$padj) & res_Suc$padj <= FDR
cutoff1_Suc <- res_Suc$log2FoldChange >= lfc & ! is.na(res_Suc$padj) & res_Suc$padj <= FDR
cutoff2_Suc <- res_Suc$log2FoldChange <= -lfc & ! is.na(res_Suc$padj) & res_Suc$padj <= FDR
SucStarvEffect <- rownames(res_Suc[cutoffs_Suc,])
SucStarvLow <- rownames(res_Suc[cutoff2_Suc,])
SucStarvHigh <- rownames(res_Suc[cutoff1_Suc,])

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs_Suc)))
message(sprintf("There are %s genes that are induced",sum(cutoff1_Suc)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2_Suc)))

#' #### MA plot
DESeq2::plotMA(res_Suc)

#' #### Volcano plot 
volcanoPlot(res_Suc)
svg(filename="VolcanoSucStarv.svg", 
    width=10, 
    height=8, 
    pointsize=12)
volcanoPlot(res_Suc)
dev.off()

#' #### Heatmap
hm <- heatmap.2(t(scale(t(vsd_T6[cutoffs_Suc,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel_T6],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))
svg(filename="HeatMapSucStarv.svg", 
    width=10, 
    height=8, 
    pointsize=12)
hm <- heatmap.2(t(scale(t(vsd_T6[cutoffs_Suc,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel_T6],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))
dev.off()

#' #### Enrichment analysis
#' ##### Downregulated genes
EnrSucLow <- gopher(SucStarvLow,background=rownames(vsd_T6)[rowSums(vsd_T6)>0],url="athaliana")
EnrSucLow$go$name
EnrSucLow$kegg$id
EnrSucLow$pfam$name
write.table(data.frame(EnrSucLow$go$id, EnrSucLow$go$padj), file = "GeneOntologySucLow.txt", sep = "\t",
            row.names = FALSE)

#' ##### Upregulated genes
EnrSucHigh <- gopher(SucStarvHigh,background=rownames(vsd_T6)[rowSums(vsd_T6)>0],url="athaliana")
EnrSucHigh$go$name
EnrSucHigh$kegg$id
EnrSucHigh$pfam$name
write.table(data.frame(EnrSucHigh$go$id, EnrSucHigh$go$padj), file = "GeneOntologySucHigh.txt", sep = "\t",
            row.names = FALSE)







#' ### Nutrition effect of phosphorus starvation
#' #' #### Extraction of the results from the DESeq analysis
res_Pi <- results(dds_T6,name = "Nutrition_NS_vs_NPS")

# res <- results(dds,contrast = c(Nutrition","NP","NS"))
# res <- results(dds,contrast = list("Nutrition_NP_vs_NPS","Nutrition_NS_vs_NPS"))
# res <- results(dds,contrast = c(0,1,-1,0,0,0,0,0,0) # NP vs. NS

cutoffs_Pi <- abs(res_Pi$log2FoldChange) >= lfc & ! is.na(res_Pi$padj) & res_Pi$padj <= FDR
cutoff1_Pi <- res_Pi$log2FoldChange >= lfc & ! is.na(res_Pi$padj) & res_Pi$padj <= FDR
cutoff2_Pi <- res_Pi$log2FoldChange <= -lfc & ! is.na(res_Pi$padj) & res_Pi$padj <= FDR
PiStarvEffect <- rownames(res_Pi[cutoffs_Pi,])
PiStarvLow <- rownames(res_Pi[cutoff2_Pi,])
PiStarvHigh <- rownames(res_Pi[cutoff1_Pi,])

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs_Pi)))
message(sprintf("There are %s genes that are induced",sum(cutoff1_Pi)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2_Pi)))

#' #### MA plot
DESeq2::plotMA(res_Pi)

#' #### Volcano plot 
volcanoPlot(res_Pi)
svg(filename="VolcanoPiStarv.svg", 
    width=10, 
    height=8, 
    pointsize=12)
volcanoPlot(res_Pi)
dev.off()

#' #### Heatmap
#' ##### Z-score
hm <- heatmap.2(t(scale(t(vsd_T6[cutoffs_Pi,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel_T6],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))
svg(filename="HeatMapPiStarv.svg", 
    width=10, 
    height=8, 
    pointsize=12)
hm <- heatmap.2(t(scale(t(vsd_T6[cutoffs_Pi,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel_T6],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))
dev.off()

#' #### Enrichment analysis
#' ##### Downregulated genes
EnrPiLow <- gopher(PiStarvLow,background=rownames(vsd_T6)[rowSums(vsd_T6)>0],url="athaliana")
EnrPiLow$go$name
EnrPiLow$kegg$id
EnrPiLow$pfam$name
write.table(data.frame(EnrPiLow$go$id, EnrPiLow$go$padj), file = "GeneOntologyPiLow.txt", sep = "\t",
            row.names = FALSE)

#' ##### Upregulated genes
EnrPiHigh <- gopher(PiStarvHigh,background=rownames(vsd_T6)[rowSums(vsd_T6)>0],url="athaliana")
EnrPiHigh$go$name
EnrPiHigh$kegg$id
EnrPiHigh$pfam$name
write.table(data.frame(EnrPiHigh$go$id, EnrPiHigh$go$padj), file = "GeneOntologyPiHigh.txt", sep = "\t",
            row.names = FALSE)




#' ### AZD effect
#' #' #### Extraction of the results from the DESeq analysis
res_AZD <- results(dds_T6,name = "AZD_AZD_vs_DMSO")

cutoffs_AZD <- abs(res_AZD$log2FoldChange) >= lfc & ! is.na(res_AZD$padj) & res_AZD$padj <= FDR
cutoff1_AZD <- res_AZD$log2FoldChange >= lfc & ! is.na(res_AZD$padj) & res_AZD$padj <= FDR
cutoff2_AZD <- res_AZD$log2FoldChange <= -lfc & ! is.na(res_AZD$padj) & res_AZD$padj <= FDR
AZDEffect <- rownames(res_AZD[cutoffs_AZD,])
AZDLow <- rownames(res_AZD[cutoff2_AZD,])
AZDHigh <- rownames(res_AZD[cutoff1_AZD,])

message(sprintf("There are %s genes that are differentially expressed",sum(cutoffs_AZD)))
message(sprintf("There are %s genes that are induced",sum(cutoff1_AZD)))
message(sprintf("There are %s genes that are repressed",sum(cutoff2_AZD)))

#' ##### MA plot
DESeq2::plotMA(res_AZD)

#' ##### Volcano plot 
volcanoPlot(res_AZD)
svg(filename="VolcanoAZD.svg", 
    width=10, 
    height=8, 
    pointsize=12)
volcanoPlot(res_AZD)
dev.off()

#' #### Enrichment analysis
#' ##### Downregulated genes
EnrAZDLow <- gopher(AZDLow,background=rownames(vsd_T6)[rowSums(vsd_T6)>0],url="athaliana")
EnrAZDLow$go$name
EnrAZDLow$kegg$id
EnrAZDLow$pfam$name
write.table(data.frame(EnrAZDLow$go$id, EnrAZDLow$go$padj), file = "GeneOntologyAZDLow.txt", sep = "\t",
            row.names = FALSE)

#' ##### Upregulated genes
EnrAZDHigh <- gopher(AZDHigh,background=rownames(vsd_T6)[rowSums(vsd_T6)>0],url="athaliana")
EnrAZDHigh$go$name
EnrAZDHigh$kegg$id
EnrAZDHigh$pfam$name
write.table(data.frame(EnrAZDHigh$go$id, EnrAZDHigh$go$padj), file = "GeneOntologyAZDHigh.txt", sep = "\t",
            row.names = FALSE)

#' #### Heatmap
#' ##### Z-score
hm <- heatmap.2(t(scale(t(vsd_T6[cutoffs_AZD,]))),col=hpal,
          Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
          distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel_T6],
          hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))
svg(filename="HeatMapAZD.svg", 
    width=10, 
    height=8, 
    pointsize=12)
hm <- heatmap.2(t(scale(t(vsd_T6[cutoffs_AZD,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel_T6],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))
dev.off()


#' ##### Heatmaps with Fold-change
#' ###### Without NPS_DMSO, average of replicates, conditions ordered alphabetically
hmfc <- heatmap.2(fc_T6[cutoffs_AZD,],
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                col=hpal, margins = c(6,5),cexCol = 0.8)

#' ###### Without NPS_DMSO, average of replicates, conditions ordered by dendrogram
hmfc <- heatmap.2(fc_T6[cutoffs_AZD,],
                  trace = "none",labRow = FALSE,
                  col=hpal, margins = c(6,5),cexCol = 0.8)

#' ###### Like above but genes organized by Pearson correlation
hmfc <- heatmap.2(fc_T6[cutoffs_AZD,],
                  trace = "none",labRow = FALSE,
                  col=hpal, margins = c(6,5),cexCol = 0.8,
                  distfun = pearson.dist)
#' ###### With NPS_DMSO, all replicates, conditions ordered by dendrogram
hmfc <- heatmap.2(fc2_T6[cutoffs_AZD,],
                  trace = "none",labRow = FALSE,
                  col=hpal, margins = c(6,5),
                  labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel_T6],cexCol = 0.8)

#' ###### With NPS_DMSO, all replicates, conditions ordered according to the experimental design
hmfc <- heatmap.2(fc2_T6[cutoffs_AZD,],
                  Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                  col=hpal, margins = c(6,5),
                  labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel_T6],cexCol = 0.8)

#' #### Cluster analysis
tree <- cutree(as.hclust(hm$rowDendrogram),6)
c1 <- names(tree)[tree==1]
message(sprintf("There are %s genes in the cluster1",length(c1)))
heatmap.2(t(scale(t(vsd_T6[c1,]))),col=hpal,
                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel_T6],
                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))
svg(filename="HeatMapCl1.svg", 
    width=10, 
    height=8, 
    pointsize=12)
heatmap.2(t(scale(t(vsd_T6[c1,]))),col=hpal,
          Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
          distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel_T6],
          hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))
dev.off()
cluster1 <- names(tree)[tree==1]
write.csv(cluster1,"cluster1.csv")

ggplot(melt(vsd_T6[cluster1,]),aes(x=Var2,y=value,group=Var2)) + geom_boxplot() +
    theme(axis.text.x = element_text(angle=90)) +
    scale_x_discrete(label=paste(samples$Nutrition,samples$AZD,sep="-")[sel_T6],element_blank()) +
    scale_y_continuous("expression (VST)")

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

svg(filename="HeatMapCl4.svg", 
    width=10, 
    height=8, 
    pointsize=12)
heatmap.2(t(scale(t(vsd_T6[c4,]))),col=hpal,
          Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
          distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel_T6],
          hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))
dev.off()

cluster4 <- names(tree)[tree==4]
write.csv(cluster4,"cluster4.csv")

ggplot(melt(vsd_T6[cluster4,]),aes(x=Var2,y=value,group=Var2)) + geom_boxplot() +
    theme(axis.text.x = element_text(angle=90)) +
    scale_x_discrete(label=paste(samples$Nutrition,samples$AZD,sep="-")[sel_T6],element_blank()) +
    scale_y_continuous("expression (VST)")

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


#' ### Analysis of cluster composition
#' #### Enrichment
#' ##### Enrichment for cluster1
cl1 <- gopher(cluster1,background=rownames(vsd_T6)[rowSums(vsd_T6)>0],url="athaliana")
cl1$go$name
cl1$kegg
cl1$pfam$name
write.csv(data.frame(cl1$go$id, cl1$go$padj),"GOofC1.csv")

#' ##### Enrichment for cluster4
cl4 <- gopher(cluster4,task=c("go","pfam","kegg"),background=rownames(vsd_T6)[rowSums(vsd_T6)>0],url="athaliana")
cl4$go$name
cl4$kegg
cl4$pfam$name
write.csv(data.frame(cl4$go$id, cl4$go$padj),"GOofC4.csv")

#' ### Comparisons of GOI lists
#' #### Comparison of the DEGs in the different conditions
#' ##### Venn diagram
# Require package VennDiagram
grid.draw(venn.diagram(list(AZDEffect, PiStarvEffect, SucStarvEffect),
                       filename=NULL,
                       col=pal[1:3],
                       category.names=c("AZD","Phos. Starv.","Suc. Starv.")))

#' #### Comparison of genes induced in the different conditions
grid.draw(venn.diagram(list(AZDHigh, PiStarvHigh, SucStarvHigh),
                       filename=NULL,
                       col=pal[1:3],
                       category.names=c("AZD+","Phos. Starv.+","Suc. Starv.+")))
#' ##### GOI overlapping lists
UpOverlap <- intersect(AZDHigh,PiStarvHigh)
UpAzdSpecific <- setdiff(AZDHigh,PiStarvHigh)
UpPiSpecific <- setdiff(PiStarvHigh,AZDHigh)
#' ##### GO enrichment analysis
message(sprintf("There are %s genes commonly induced",length(UpOverlap)))
EnrUpOverlap <- gopher(UpOverlap,background=rownames(vsd_T6)[rowSums(vsd_T6)>0],url="athaliana")
EnrUpOverlap$go$name
EnrUpOverlap$kegg$id
EnrUpOverlap$pfam$name
write.table(data.frame(EnrUpOverlap$go$id, EnrUpOverlap$go$padj), file = "GeneOntologyUpOverlap.txt", sep = "\t",
            row.names = FALSE)

message(sprintf("There are %s genes specifically induced by AZD",length(UpAzdSpecific)))
EnrUpAZD <- gopher(UpAzdSpecific,background=rownames(vsd_T6)[rowSums(vsd_T6)>0],url="athaliana")
EnrUpAZD$go$name
EnrUpAZD$kegg$id
EnrUpAZD$pfam$name
write.table(data.frame(EnrUpAZD$go$id, EnrUpAZD$go$padj), file = "GeneOntologyUpAZD.txt", sep = "\t",
            row.names = FALSE)

message(sprintf("There are %s genes specifically induced by phosphate starvation",length(UpPiSpecific)))
EnrUpPi <- gopher(UpPiSpecific,background=rownames(vsd_T6)[rowSums(vsd_T6)>0],url="athaliana")
EnrUpPi$go$name
EnrUpPi$kegg$id
EnrUpPi$pfam$name
write.table(data.frame(EnrUpPi$go$id, EnrUpPi$go$padj), file = "GeneOntologyUpPi.txt", sep = "\t",
            row.names = FALSE)

#' #### Comparison of genes repressed in the different conditions
grid.draw(venn.diagram(list(AZDLow, PiStarvLow, SucStarvLow),
                       filename=NULL,
                       col=pal[1:3],
                       category.names=c("AZDlow","Phos. Starv.low","Suc. Starv.low")))
#' ##### GOI overlapping lists
DownOverlap <- intersect(AZDLow,PiStarvLow)
DownAzdSpecific <- setdiff(AZDLow,PiStarvLow)
DownPiSpecific <- setdiff(PiStarvLow,AZDLow)

#' ##### GO enrichment analysis
message(sprintf("There are %s genes commonly repressed",length(DownOverlap)))
EnrDownOverlap <- gopher(DownOverlap,background=rownames(vsd_T6)[rowSums(vsd_T6)>0],url="athaliana")
EnrDownOverlap$go$name
EnrDownOverlap$kegg$id
EnrDownOverlap$pfam$name
write.table(data.frame(EnrDownOverlap$go$id, EnrDownOverlap$go$padj), file = "GeneOntologyDownOverlap.txt", sep = "\t",
            row.names = FALSE)

message(sprintf("There are %s genes specifically repressed by AZD",length(DownAzdSpecific)))
EnrDownAZD <- gopher(DownAzdSpecific,background=rownames(vsd_T6)[rowSums(vsd_T6)>0],url="athaliana")
EnrDownAZD$go$name
EnrDownAZD$kegg$id
EnrDownAZD$pfam$name
write.table(data.frame(EnrDownAZD$go$id, EnrDownAZD$go$padj), file = "GeneOntologyDownAZD.txt", sep = "\t",
            row.names = FALSE)

message(sprintf("There are %s genes specifically repressed by phosphate starvation",length(DownPiSpecific)))
EnrDownPi <- gopher(DownPiSpecific,background=rownames(vsd_T6)[rowSums(vsd_T6)>0],url="athaliana")
EnrDownPi$go$name
EnrDownPi$kegg$id
EnrDownPi$pfam$name
write.table(data.frame(EnrDownPi$go$id, EnrDownPi$go$padj), file = "GeneOntologyDownPi.txt", sep = "\t",
            row.names = FALSE)

#' ### Comparisons of treatments regarding the Cell cycle GO term
# Create list of important GO terms based on the GO terms enriched in the genes downregulated after Pi starvation
myList <- c("GO:0007049", "GO:0044770", "GO:0016458", "GO:0000278", "GO:0006720", "GO:0007017", "GO:0040029", 
            "GO:0022900", "GO:0008283", "GO:0051301")

myList <- c("GO:0007049","GO:0090066","GO:0055086","GO:0048522","GO:0048519","GO:0009259","GO:0051726","GO:0006553",
            "GO:0044770","GO:0009825","GO:0006090","GO:0022900","GO:0007059","GO:0060249","GO:0051301","GO:0019637",
            "GO:0007267","GO:0050789","GO:0007017","GO:0010608","GO:0040029","GO:0046471","GO:0046490","GO:0040008",
            "GO:0008652","GO:0008610","GO:0016458","GO:0006720","GO:0051171")


AZDLow <- AZDLow
PiStarvLow <- PiStarvLow

# Create list of genes in each GO term
geneList <- gopher(myList, task="go", url="athaliana", endpoint="term-to-gene")

#AZD_CellCycle <- lapply(geneList, function(x){intersect(x, AZDLow)})
#Pi_CellCycle <- lapply(geneList, function(x){intersect(x, PiStarvLow)})

# Remove the GO term having one or less gene
geneList <- geneList[lapply(geneList,length)>2]
geneList
length(geneList)
lengths(geneList)

lapply(geneList, function(x){length(x)})
lapply(geneList, function(x){length(intersect(x,AZDLow))})
lapply(geneList, function(x){length(intersect(x,AZDHigh))})

lapply(geneList, function(x){
    goi <- res_AZD[rownames(res_AZD) %in% x,]
    volcanoPlot(goi)
    goi2 <- t(scale(t(vsd_T6[rownames(res_AZD) %in% x,])))
    goi2 <- na.omit(goi2)
    hm <- heatmap.2(goi2,col=hpal,
                    Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                    distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel_T6],
                    hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))
})


lapply(geneList, function(x){length(x)})
lapply(geneList, function(x){length(intersect(x,PiStarvLow))})
lapply(geneList, function(x){length(intersect(x,PiStarvHigh))})

lapply(geneList, function(x){
    goi <- res_Pi[rownames(res_Pi) %in% x,]
    volcanoPlot(goi)
    goi2 <- t(scale(t(vsd_T6[rownames(res_Pi) %in% x,])))
    goi2 <- na.omit(goi2)
    hm <- heatmap.2(goi2,col=hpal,
                    Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
                    distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel_T6],
                    hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))
})
























#' ###### Expression profile
#' * Differential expression in T0 and T6
sel_T0 <- samples$Timepoint %in% c("T6","T0")
suppressMessages(dds_T0 <- DESeqDataSetFromMatrix(
    countData = kg[,sel_T0],
    colData = samples[sel_T0,],
    design = ~ Conditions ))

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
#resultsNames(dds_T0)

#message(sprintf("There are %s genes of interest",length(goi)))

#' #### Heatmap for goi1
#' ##### Z-score
#hm <- heatmap.2(t(scale(t(vsd_T0[goi,]))),col=hpal,
#                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
#                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel_T0],
#                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))

#' ##### Fold-change
#hmfc <- heatmap.2(fc2_T0[goi,],
#                  Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
#                  labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel_T0],
#                  col=hpal, margins = c(6,5),cexCol = 0.8)

#' #### Heatmap for goi2
#' ##### Z-score
#hm <- heatmap.2(t(scale(t(vsd_T0[goi2,]))),col=hpal,
#                Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
#                distfun = pearson.dist, labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel_T0],
#                hclustfun = function(X){hclust(X,method="ward.D2")},margins = c(6,5))

#' ##### Fold-change
#hmfc <- heatmap.2(fc2_T0[goi2,],
#                  Colv=FALSE,dendrogram = "row",trace = "none",labRow = FALSE,
#                  labCol = paste(samples$Nutrition,samples$AZD,sep="-")[sel_T0],
#                  col=hpal, margins = c(6,5),cexCol = 0.8)











#' # Expression level of the TOR complex
#' * Preparation of the complex member list
sel1 <- c("AT1G50030","AT3G18140","AT2G22040","AT3G08850","AT5G01770")
#' * Preparation of the hexokinase list
sel2 <- c("AT1G05205","AT1G47840","AT2G19860","AT4G29130")
#' * Preparation of the AtHXK1
sel3 <- c("AT1G47840","AT4G29130")


j <- sel1
    expr <- vsd_T0[match(j,rownames(vsd_T0)),]
    #expr <- vsd_T0[match(rownames(vsd_T0),sel1),]
    #expr <- na.omit(expr)
    AvgTOR <- sapply(split.data.frame(t(expr),f = droplevels(samples$Conditions[sel_T0])),colMeans)
    SDTOR <- sapply(split.data.frame(t(expr),f = droplevels(samples$Conditions[sel_T0])),colSds)
    barplot(AvgTOR,beside=T, cex.names=0.7,legend.text = j,ylim=c(0,3.5))
    barcenters <- barplot(AvgTOR,beside=T, cex.names=0.7,legend.text = j,ylim=c(0,3.5))
    arrows(x0=barcenters, x1=barcenters, y0=AvgTOR-SDTOR, y1=AvgTOR+SDTOR,
           lwd=1.5, angle=90, length=0.05,code=3)
    for (i in 1:length(AvgTOR[,1]))
        {
        plot(AvgTOR[i,],ylim=c(min(AvgTOR[i,]-SDTOR[i,]),max(AvgTOR[i,]+SDTOR[i,])),
             xaxt="n",xlab="",ylab=sprintf("gene = %s",rownames(AvgTOR)[i]))
        axis(side = 1,at=1:7, colnames(AvgTOR),cex.axis=0.7,las=2)
        arrows(x0=1:length(AvgTOR[i,]), x1=1:length(AvgTOR[i,]),
           y0=AvgTOR[i,]-SDTOR[i,], y1=AvgTOR[i,]+SDTOR[i,],
               lwd=1.5, angle=90, length=0.05,code=3)
        }

    j <- sel2
    expr <- vsd_T0[match(j,rownames(vsd_T0)),]
    #expr <- vsd_T0[match(rownames(vsd_T0),sel1),]
    #expr <- na.omit(expr)
    AvgTOR <- sapply(split.data.frame(t(expr),f = droplevels(samples$Conditions[sel_T0])),colMeans)
    SDTOR <- sapply(split.data.frame(t(expr),f = droplevels(samples$Conditions[sel_T0])),colSds)
    barplot(AvgTOR,beside=T, cex.names=0.7,legend.text = j,ylim=c(0,3.5))
    barcenters <- barplot(AvgTOR,beside=T, cex.names=0.7,legend.text = j,ylim=c(0,3.5))
    arrows(x0=barcenters, x1=barcenters, y0=AvgTOR-SDTOR, y1=AvgTOR+SDTOR,
           lwd=1.5, angle=90, length=0.05,code=3)
    for (i in 1:length(AvgTOR[,1]))
    {
        plot(AvgTOR[i,],ylim=c(min(AvgTOR[i,]-SDTOR[i,]),max(AvgTOR[i,]+SDTOR[i,])),
             xaxt="n",xlab="",ylab=sprintf("gene = %s",rownames(AvgTOR)[i]))
        axis(side = 1,at=1:7, colnames(AvgTOR),cex.axis=0.7,las=2)
        arrows(x0=1:length(AvgTOR[i,]), x1=1:length(AvgTOR[i,]),
               y0=AvgTOR[i,]-SDTOR[i,], y1=AvgTOR[i,]+SDTOR[i,],
               lwd=1.5, angle=90, length=0.05,code=3)
    }

    j <- sel3
    expr <- vsd_T0[match(j,rownames(vsd_T0)),]
    #expr <- vsd_T0[match(rownames(vsd_T0),sel1),]
    #expr <- na.omit(expr)
    AvgTOR <- sapply(split.data.frame(t(expr),f = droplevels(samples$Conditions[sel_T0])),colMeans)
    SDTOR <- sapply(split.data.frame(t(expr),f = droplevels(samples$Conditions[sel_T0])),colSds)
    barplot(AvgTOR,beside=T, cex.names=0.7,legend.text = j,ylim=c(0,3.5))
    barcenters <- barplot(AvgTOR,beside=T, cex.names=0.7,legend.text = j,ylim=c(0,3.5))
    arrows(x0=barcenters, x1=barcenters, y0=AvgTOR-SDTOR, y1=AvgTOR+SDTOR,
           lwd=1.5, angle=90, length=0.05,code=3)
    for (i in 1:length(AvgTOR[,1]))
    {
        plot(AvgTOR[i,],ylim=c(min(AvgTOR[i,]-SDTOR[i,]),max(AvgTOR[i,]+SDTOR[i,])),
             xaxt="n",xlab="",ylab=sprintf("gene = %s",rownames(AvgTOR)[i]))
        axis(side = 1,at=1:7, colnames(AvgTOR),cex.axis=0.7,las=2)
        arrows(x0=1:length(AvgTOR[i,]), x1=1:length(AvgTOR[i,]),
               y0=AvgTOR[i,]-SDTOR[i,], y1=AvgTOR[i,]+SDTOR[i,],
               lwd=1.5, angle=90, length=0.05,code=3)
    }

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```

