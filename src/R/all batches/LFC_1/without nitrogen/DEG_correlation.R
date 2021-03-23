#' ---
#' title: "DEG analysis, GO analysis and bubble plot preparation"
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
    library(plyr)
    library(RColorBrewer)
    library(reshape2)
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
#' Different from Schurch et al., RNA, 2016
lfc <- 0.5
FDR <- 0.01

#' * Create palettes
pal <- c(brewer.pal(8,"Dark2"),1)
pal2 <- brewer.pal(9,"Paired") #require package RColorBrewer
cols <- rainbow(17)
hpal <- colorRampPalette(c("blue","white","red"))(100)

#' * Register the default plot margin
mar <- par("mar")

#' # Expression data
#' ## Loading of Fold changes and adjusted p-values
#' * Read the expression data for 6 hours
T6_NP_DMSO <- read.csv("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/Conditions_T6_NP_DMSO_vs_T6_NPS_DMSO.csv",sep=";",dec=",")
T6_NS_DMSO <- read.csv("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/Conditions_T6_NS_DMSO_vs_T6_NPS_DMSO.csv",sep=";",dec=",")
T6_NPS_AZD <- read.csv("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/Conditions_T6_NPS_AZD_vs_T6_NPS_DMSO.csv",sep=";",dec=",")
T6_NP_AZD <- read.csv("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/Conditions_T6_NP_AZD_vs_T6_NPS_DMSO.csv",sep=";",dec=",")
T6_NS_AZD <- read.csv("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/Conditions_T6_NS_AZD_vs_T6_NPS_DMSO.csv",sep=";",dec=",")

#' * Read the expression data for 24 hours
T24_NP_DMSO <- read.csv("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/Conditions_T24_NP_DMSO_vs_T24_NPS_DMSO.csv",sep=";",dec=",")
T24_NS_DMSO <- read.csv("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/Conditions_T24_NS_DMSO_vs_T24_NPS_DMSO.csv",sep=";",dec=",")
T24_NPS_AZD <- read.csv("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/Conditions_T24_NPS_AZD_vs_T24_NPS_DMSO.csv",sep=";",dec=",")
T24_NP_AZD <- read.csv("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/Conditions_T24_NP_AZD_vs_T24_NPS_DMSO.csv",sep=";",dec=",")
T24_NS_AZD <- read.csv("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/Conditions_T24_NS_AZD_vs_T24_NPS_DMSO.csv",sep=";",dec=",")

#' * Read Wouter's data
Seedling <- read.csv("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Wouter_short.csv",sep=";",dec=",")
Seedling <- cbind(Seedling["x"],Seedling["Coef.acol...col"],Seedling["p.value.adj.acol...col"])
colnames(Seedling) <- c("rownames.res.","res.log2FoldChange","res.padj")

#' ## Combine the files while keeping the non significant genes
#' * Add the experimental annotation to the files
T6_NP_DMSO <- cbind(Description='T6_NP_DMSO', T6_NP_DMSO)
T6_NS_DMSO <- cbind(Description='T6_NS_DMSO', T6_NS_DMSO)
T6_NPS_AZD <- cbind(Description='T6_NPS_AZD', T6_NPS_AZD)
T6_NP_AZD <- cbind(Description='T6_NP_AZD', T6_NP_AZD)
T6_NS_AZD <- cbind(Description='T6_NS_AZD', T6_NS_AZD)
T24_NP_DMSO <- cbind(Description='T24_NP_DMSO', T24_NP_DMSO)
T24_NS_DMSO <- cbind(Description='T24_NS_DMSO', T24_NS_DMSO)
T24_NPS_AZD <- cbind(Description='T24_NPS_AZD', T24_NPS_AZD)
T24_NP_AZD <- cbind(Description='T24_NP_AZD', T24_NP_AZD)
T24_NS_AZD <- cbind(Description='T24_NS_AZD', T24_NS_AZD)
Seedling <- cbind(Description='Seedling', Seedling)

#' * Rbind the files
combined1 <- rbind.fill(T6_NP_DMSO,T6_NS_DMSO,T6_NPS_AZD,T6_NP_AZD,T6_NS_AZD,T24_NP_DMSO,T24_NS_DMSO,T24_NPS_AZD,
                       T24_NP_AZD,T24_NS_AZD,Seedling)

data1 <- dcast(combined, rownames.res. ~ Description, mean, value.var = "res.log2FoldChange") #Require reshape2 package
rownames(data1) <- data1[,1]
data1 <- data1[,2:length(data1)]

#' ## Remove non significant genes
T6_NP_DMSO$res.log2FoldChange[T6_NP_DMSO$res.padj >= FDR | is.na(T6_NP_DMSO$res.padj)] <- NA 
T6_NS_DMSO$res.log2FoldChange[T6_NS_DMSO$res.padj >= FDR | is.na(T6_NS_DMSO$res.padj)] <- NA 
T6_NPS_AZD$res.log2FoldChange[T6_NPS_AZD$res.padj >= FDR | is.na(T6_NPS_AZD$res.padj)] <- NA 
T6_NP_AZD$res.log2FoldChange[T6_NP_AZD$res.padj >= FDR | is.na(T6_NP_AZD$res.padj)] <- NA 
T6_NS_AZD$res.log2FoldChange[T6_NS_AZD$res.padj >= FDR | is.na(T6_NS_AZD$res.padj)] <- NA 

T24_NP_DMSO$res.log2FoldChange[T24_NP_DMSO$res.padj >= FDR | is.na(T24_NP_DMSO$res.padj)] <- NA 
T24_NS_DMSO$res.log2FoldChange[T24_NS_DMSO$res.padj >= FDR | is.na(T24_NS_DMSO$res.padj)] <- NA 
T24_NPS_AZD$res.log2FoldChange[T24_NPS_AZD$res.padj >= FDR | is.na(T24_NPS_AZD$res.padj)] <- NA 
T24_NP_AZD$res.log2FoldChange[T24_NP_AZD$res.padj >= FDR | is.na(T24_NP_AZD$res.padj)] <- NA 
T24_NS_AZD$res.log2FoldChange[T24_NS_AZD$res.padj >= FDR | is.na(T24_NS_AZD$res.padj)] <- NA 

Seedling$res.log2FoldChange[Seedling$res.padj >= FDR | is.na(Seedling$res.padj)] <- NA 

#' ## Combine the files together after removal of the non-significant genes

#' * Rbind the files
combined <- rbind.fill(T6_NP_DMSO,T6_NS_DMSO,T6_NPS_AZD,T6_NP_AZD,T6_NS_AZD,T24_NP_DMSO,T24_NS_DMSO,T24_NPS_AZD,
                          T24_NP_AZD,T24_NS_AZD,Seedling)

#' * Create the matrix
data <- dcast(combined, rownames.res. ~ Description, mean, value.var = "res.log2FoldChange") #Require reshape2 package
rownames(data) <- data[,1]
data <- data[,2:length(data)]

#' * Create the plots
plot(data)
plot(data$T24_NS_DMSO ~ data$T24_NS_AZD,xlim=c(-4,4),ylim=c(-4,4),)

plot(data$T6_NP_DMSO ~ data$T6_NP_AZD,
     xlim=c(min(data, na.rm=T),max(data,na.rm=T)),
     ylim=c(min(data, na.rm=T),max(data,na.rm=T)),
     pch=20,cex=1)

plot(data,
     xlim=c(min(data, na.rm=T),max(data,na.rm=T)),
     ylim=c(min(data, na.rm=T),max(data,na.rm=T)),
     pch=20,cex=1)


#' ## Calculate the correlations
cor(data$T6_NS_DMSO,data$T6_NS_AZD,method = "pearson", use= "pairwise.complete.obs")

#' # Identify proportions of deregulated genes in specific GOs
GO0006468 <- read.csv("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/GO0006468.csv",sep=";",dec=",")
GO0006468 <- GO0006468$X..gene_id
GO <- combined1[combined1$rownames.res. %in% GO0006468,]
ggplot(GO, aes(x=Description, y=res.log2FoldChange)) + 
    geom_violin(trim=FALSE, fill="gray")+
    labs(title="GO:0006468",x="Comparisons", y = "Log2FC")+
    geom_boxplot(width=0.1)+
    theme_classic()

#' # Identify proportions of deregulated genes in specific gene lists
#' ## Load gene lists
PiSpecLow6hrs <- read.csv("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/PiSpecLow6hrs.csv", sep="\t")
PiSpecLow24hrs <- read.csv("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/PiSpecLow24hrs.csv", sep="\t")
PiSpecHigh6hrs <- read.csv("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/PiSpecHigh6hrs.csv", sep="\t")
PiSpecHigh24hrs <- read.csv("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/PiSpecHigh24hrs.csv", sep="\t")

AzdSpecLow6hrs <- read.csv("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/AzdSpecLow6hrs.csv", sep="\t")
AzdSpecLow24hrs <- read.csv("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/AzdSpecLow24hrs.csv", sep="\t")
AzdSpecHigh6hrs <- read.csv("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/AzdSpecHigh6hrs.csv", sep="\t")
AzdSpecHigh24hrs <- read.csv("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/AzdSpecHigh24hrs.csv", sep="\t")

#' ## Select the list in the all gene list
GO <- combined1[combined1$rownames.res. %in% union(AzdSpecHigh6hrs$x,AzdSpecLow6hrs$x),]
GO <- combined[combined$rownames.res. %in% AzdSpecHigh24hrs$x,]
list <- data[rownames(data) %in% union(AzdSpecHigh24hrs$x,AzdSpecLow24hrs$x),]
plot(list$T24_NS_AZD,list$T24_NS_DMSO)

#' ## Make the violin plot
ggplot(GO, aes(x=Description, y=res.log2FoldChange)) + 
    geom_violin(trim=T, fill="gray")+
    labs(title="GO:0006468",x="Comparisons", y = "Log2FC")+
    geom_boxplot(width=0.1)+
    theme_classic()

#' # Identification of the genes that are differentially deregulated between +/-AZD in T24-NS samples
#' ## Removing non-significant genes based on the LFC and adjusted p-values
T24_NS_DMSO$res.log2FoldChange[T24_NS_DMSO$res.padj >= FDR | abs(T24_NS_DMSO$res.log2FoldChange) < lfc  | is.na(T24_NS_DMSO$res.padj)] <- NA 
T24_NS_AZD$res.log2FoldChange[T24_NS_AZD$res.padj >= FDR | abs(T24_NS_AZD$res.log2FoldChange) < lfc  | is.na(T24_NS_AZD$res.padj)] <- NA 
plot(T24_NS_DMSO$res.log2FoldChange,T24_NS_AZD$res.log2FoldChange, xlim=c(-4,4),ylim=c(-4,4))

#' ## Select the genes presenting opposite patterns
HighInAzd <- intersect(T24_NS_DMSO$rownames.res.[T24_NS_DMSO$res.log2FoldChange<0],T24_NS_AZD$rownames.res.[T24_NS_AZD$res.log2FoldChange>0])
LowInAzd <- intersect(T24_NS_DMSO$rownames.res.[T24_NS_DMSO$res.log2FoldChange>0],T24_NS_AZD$rownames.res.[T24_NS_AZD$res.log2FoldChange<0])
HighInAzd <- HighInAzd[2:length(HighInAzd)]
LowInAzd <- LowInAzd[2:length(LowInAzd)]

#' ## Make a GO enrichment
GO <- gopher(HighInAzd, url="athaliana", alpha=2)
GO_High <- GO$go
GO <- gopher(LowInAzd, url="athaliana", alpha=2)
GO_Low <- GO$go



#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```

