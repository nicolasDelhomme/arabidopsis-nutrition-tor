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

#' # Preparing the data
#' ## Loading the data

counts <- read.csv("Averaged_Normalized_Counts.csv", sep = ";", dec=",")
genes <- read.csv("Averaged_Normalized_Counts_genes.csv", sep = ";", dec=",")
stopifnot(nrow(counts)==nrow(genes))

rownames(counts) <- as.vector(genes[,1])


#' ## Normalization by the average
avg <- rowMeans(counts)
norm <- log2(counts/avg)

#' # Heatmap of the cell cycle related genes
#' ## Call the list of genes (Menges et al., 2003) > table 2
CC <- read.csv("Menges2.csv", sep = ";", dec=",")

#' ## Call the list of genes (Menges et al., 2003) > supplemental tables 5 and 6
CC <- read.csv("Menges.csv", sep = ";", dec=",")

#' ## Separate the lists according to the phase
Mphase <- CC$AGI[CC$Phase=="M"]
G1phase <- CC$AGI[CC$Phase=="G1"]
G2phase <- CC$AGI[CC$Phase=="G2"]
Sphase <- CC$AGI[CC$Phase=="S"]

#' ## Split the normalized count list per phase
M <- as.matrix(norm[rownames(norm) %in% Mphase,])
G1 <- as.matrix(norm[rownames(norm) %in% G1phase,])
G2 <- as.matrix(norm[rownames(norm) %in% G2phase,])
S <- as.matrix(norm[rownames(norm) %in% Sphase,])

#' ## Correction of the tables
#' ### Removal of the infinte value
M[!is.finite(M)] <- 0 ; M[is.na(M)] <- 0
G1[!is.finite(G1)] <- 0 ; G1[is.na(G1)] <- 0
G2[!is.finite(G2)] <- 0 ; G2[is.na(G2)] <- 0
S[!is.finite(S)] <- 0 ; S[is.na(S)] <- 0

#' ### Correction to fit a (-3;3) range
M[M < -3] <- -3; M[M > 3] <- 3
G1[G1 < -3] <- -3; G1[G1 > 3] <- 3
G2[G2 < -3] <- -3; G2[G2 > 3] <- 3
S[S < -3] <- -3; S[S > 3] <- 3




#' ### Heatmap
heatmap.2(M,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

heatmap.2(M,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})

heatmap.2(G1,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

heatmap.2(G1,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})

heatmap.2(G2,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

heatmap.2(G2,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})

heatmap.2(S,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

heatmap.2(S,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})


#' ## Normalization by the T0
T0 <- counts[,1]
normT0 <- log2(counts/T0)

#' # Heatmap of the cell cycle related genes
#' ## Call the list of genes
CC <- read.csv("Menges2.csv", sep = ";", dec=",")

#' ## Separate the lists according to the phase
Mphase <- CC$AGI[CC$Phase=="M"]
G1phase <- CC$AGI[CC$Phase=="G1"]
G2phase <- CC$AGI[CC$Phase=="G2"]
Sphase <- CC$AGI[CC$Phase=="S"]

#' ## Split the normalized count list per phase
M <- as.matrix(normT0[rownames(normT0) %in% Mphase,])
G1 <- as.matrix(normT0[rownames(normT0) %in% G1phase,])
G2 <- as.matrix(normT0[rownames(normT0) %in% G2phase,])
S <- as.matrix(normT0[rownames(normT0) %in% Sphase,])

#' ## Correction of the tables
#' ### Removal of the infinte value
M[!is.finite(M)] <- 0 ; M[is.na(M)] <- 0
G1[!is.finite(G1)] <- 0 ; G1[is.na(G1)] <- 0
G2[!is.finite(G2)] <- 0 ; G2[is.na(G2)] <- 0
S[!is.finite(S)] <- 0 ; S[is.na(S)] <- 0

#' ### Correction to fit a (-3;3) range
M[M < -3] <- -3; M[M > 3] <- 3
G1[G1 < -3] <- -3; G1[G1 > 3] <- 3
G2[G2 < -3] <- -3; G2[G2 > 3] <- 3
S[S < -3] <- -3; S[S > 3] <- 3




#' ### Heatmap
heatmap.2(M,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

heatmap.2(M,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})

heatmap.2(G1,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

heatmap.2(G1,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})

heatmap.2(G2,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

heatmap.2(G2,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})

heatmap.2(S,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

heatmap.2(S,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})


#' ## Normalization by the T6_NPS_DMSO
T6 <- counts[,2]
normT6 <- log2(counts/T6)

#' # Heatmap of the cell cycle related genes
#' ## Call the list of genes
CC <- read.csv("Menges2.csv", sep = ";", dec=",")

#' ## Separate the lists according to the phase
Mphase <- CC$AGI[CC$Phase=="M"]
G1phase <- CC$AGI[CC$Phase=="G1"]
G2phase <- CC$AGI[CC$Phase=="G2"]
Sphase <- CC$AGI[CC$Phase=="S"]

#' ## Split the normalized count list per phase
M <- as.matrix(normT6[rownames(normT6) %in% Mphase,2:7])
G1 <- as.matrix(normT6[rownames(normT6) %in% G1phase,2:7])
G2 <- as.matrix(normT6[rownames(normT6) %in% G2phase,2:7])
S <- as.matrix(normT6[rownames(normT6) %in% Sphase,2:7])

#' ## Correction of the tables
#' ### Removal of the infinte value
M[!is.finite(M)] <- 0 ; M[is.na(M)] <- 0
G1[!is.finite(G1)] <- 0 ; G1[is.na(G1)] <- 0
G2[!is.finite(G2)] <- 0 ; G2[is.na(G2)] <- 0
S[!is.finite(S)] <- 0 ; S[is.na(S)] <- 0

#' ### Correction to fit a (-3;3) range
M[M < -3] <- -3; M[M > 3] <- 3
G1[G1 < -3] <- -3; G1[G1 > 3] <- 3
G2[G2 < -3] <- -3; G2[G2 > 3] <- 3
S[S < -3] <- -3; S[S > 3] <- 3




#' ### Heatmap
heatmap.2(M,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

heatmap.2(M,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})

heatmap.2(G1,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

heatmap.2(G1,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})

heatmap.2(G2,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

heatmap.2(G2,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})

heatmap.2(S,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

heatmap.2(S,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})


#' ## Normalization by the T24
T24 <- counts[,8]
normT24 <- log2(counts/T24)

#' # Heatmap of the cell cycle related genes
#' ## Call the list of genes
CC <- read.csv("Menges2.csv", sep = ";", dec=",")

#' ## Separate the lists according to the phase
Mphase <- CC$AGI[CC$Phase=="M"]
G1phase <- CC$AGI[CC$Phase=="G1"]
G2phase <- CC$AGI[CC$Phase=="G2"]
Sphase <- CC$AGI[CC$Phase=="S"]

#' ## Split the normalized count list per phase
M <- as.matrix(normT24[rownames(normT24) %in% Mphase,8:13])
G1 <- as.matrix(normT24[rownames(normT24) %in% G1phase,8:13])
G2 <- as.matrix(normT24[rownames(normT24) %in% G2phase,8:13])
S <- as.matrix(normT24[rownames(normT24) %in% Sphase,8:13])

#' ## Correction of the tables
#' ### Removal of the infinte value
M[!is.finite(M)] <- 0 ; M[is.na(M)] <- 0
G1[!is.finite(G1)] <- 0 ; G1[is.na(G1)] <- 0
G2[!is.finite(G2)] <- 0 ; G2[is.na(G2)] <- 0
S[!is.finite(S)] <- 0 ; S[is.na(S)] <- 0

#' ### Correction to fit a (-3;3) range
M[M < -3] <- -3; M[M > 3] <- 3
G1[G1 < -3] <- -3; G1[G1 > 3] <- 3
G2[G2 < -3] <- -3; G2[G2 > 3] <- 3
S[S < -3] <- -3; S[S > 3] <- 3




#' ### Heatmap
heatmap.2(M,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

heatmap.2(M,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})

heatmap.2(G1,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

heatmap.2(G1,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})

heatmap.2(G2,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

heatmap.2(G2,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})

heatmap.2(S,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

heatmap.2(S,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})












#' # Heatmap of the genes related to the cell cycle and TOR
#' ## Call the list of genes (Xiong et al, 2013) > supplemental table 6 full
CC <- read.csv("Xiong_Nature_complete.csv", sep = ";", dec=",")

#' ## Call the list of genes (Xiong et al, 2013) > supplemental table 6 only genes with |log (FC)| > 0,5
CC <- read.csv("Xiong_Nature_lfc05.csv", sep = ";", dec=",")
norm <- cbind(norm, AGI=rownames(norm))
norm <- as.matrix(merge(norm,CC,by="AGI"))

#' ## Separate the lists according to the phase
Mphase <- norm[norm[,15]=="M",] ; rownames(Mphase) <- Mphase[,1] ; Mphase <- Mphase[,-15] ; Mphase <- Mphase[,-1]
G1phase <- norm[norm[,15]=="G1",] ; rownames(G1phase) <- G1phase[,1] ; G1phase <- G1phase[,-15] ; G1phase <- G1phase[,-1]
G2phase <- norm[norm[,15]=="G2",] ; rownames(G2phase) <- G2phase[,1] ; G2phase <- G2phase[,-15] ; G2phase <- G2phase[,-1]
Sphase <- norm[norm[,15]=="S",] ; rownames(Sphase) <- Sphase[,1] ; Sphase <- Sphase[,-15] ; Sphase <- Sphase[,-1]

#' ## Prepare the tables as matrix of numeric values
M <- matrix(as.numeric(Mphase), ncol = ncol(Mphase)) ; colnames(M) <- colnames(Mphase); rownames(M) <- rownames(Mphase)
G1 <- matrix(as.numeric(G1phase), ncol = ncol(G1phase)) ; colnames(G1) <- colnames(G1phase); rownames(G1) <- rownames(G1phase)
G2 <- matrix(as.numeric(G2phase), ncol = ncol(G2phase)) ; colnames(G2) <- colnames(G2phase); rownames(G2) <- rownames(G2phase)
S <- matrix(as.numeric(Sphase), ncol = ncol(Sphase)) ; colnames(S) <- colnames(Sphase); rownames(S) <- rownames(Sphase)

#' ## Correction of the tables
#' ### Removal of the infinte value
M[!is.finite(M)] <- 0 ; M[is.na(M)] <- 0
G1[!is.finite(G1)] <- 0 ; G1[is.na(G1)] <- 0
G2[!is.finite(G2)] <- 0 ; G2[is.na(G2)] <- 0
S[!is.finite(S)] <- 0 ; S[is.na(S)] <- 0

#' ### Correction to fit a (-3;3) range
M[M < -3] <- -3; M[M > 3] <- 3
G1[G1 < -3] <- -3; G1[G1 > 3] <- 3
G2[G2 < -3] <- -3; G2[G2 > 3] <- 3
S[S < -3] <- -3; S[S > 3] <- 3


#' ### Heatmap
heatmap.2(M,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

heatmap.2(M,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})

heatmap.2(G1,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

heatmap.2(G1,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})

heatmap.2(G2,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

heatmap.2(G2,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})

heatmap.2(S,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

heatmap.2(S,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})


#' ## Normalization by the T0
T0 <- counts[,1]
normT0 <- log2(counts/T0)

#' # Heatmap of the cell cycle related genes
#' ## Call the list of genes
CC <- read.csv("Menges2.csv", sep = ";", dec=",")

#' ## Separate the lists according to the phase
Mphase <- CC$AGI[CC$Phase=="M"]
G1phase <- CC$AGI[CC$Phase=="G1"]
G2phase <- CC$AGI[CC$Phase=="G2"]
Sphase <- CC$AGI[CC$Phase=="S"]

#' ## Split the normalized count list per phase
M <- as.matrix(normT0[rownames(normT0) %in% Mphase,])
G1 <- as.matrix(normT0[rownames(normT0) %in% G1phase,])
G2 <- as.matrix(normT0[rownames(normT0) %in% G2phase,])
S <- as.matrix(normT0[rownames(normT0) %in% Sphase,])

#' ## Correction of the tables
#' ### Removal of the infinte value
M[!is.finite(M)] <- 0 ; M[is.na(M)] <- 0
G1[!is.finite(G1)] <- 0 ; G1[is.na(G1)] <- 0
G2[!is.finite(G2)] <- 0 ; G2[is.na(G2)] <- 0
S[!is.finite(S)] <- 0 ; S[is.na(S)] <- 0

#' ### Correction to fit a (-3;3) range
M[M < -3] <- -3; M[M > 3] <- 3
G1[G1 < -3] <- -3; G1[G1 > 3] <- 3
G2[G2 < -3] <- -3; G2[G2 > 3] <- 3
S[S < -3] <- -3; S[S > 3] <- 3




#' ### Heatmap
heatmap.2(M,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

heatmap.2(M,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})

heatmap.2(G1,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

heatmap.2(G1,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})

heatmap.2(G2,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

heatmap.2(G2,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})

heatmap.2(S,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

heatmap.2(S,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})


#' ## Normalization by the T6_NPS_DMSO
T6 <- counts[,2]
normT6 <- log2(counts/T6)

#' # Heatmap of the cell cycle related genes
#' ## Call the list of genes
CC <- read.csv("Menges2.csv", sep = ";", dec=",")

#' ## Separate the lists according to the phase
Mphase <- CC$AGI[CC$Phase=="M"]
G1phase <- CC$AGI[CC$Phase=="G1"]
G2phase <- CC$AGI[CC$Phase=="G2"]
Sphase <- CC$AGI[CC$Phase=="S"]

#' ## Split the normalized count list per phase
M <- as.matrix(normT6[rownames(normT6) %in% Mphase,2:7])
G1 <- as.matrix(normT6[rownames(normT6) %in% G1phase,2:7])
G2 <- as.matrix(normT6[rownames(normT6) %in% G2phase,2:7])
S <- as.matrix(normT6[rownames(normT6) %in% Sphase,2:7])

#' ## Correction of the tables
#' ### Removal of the infinte value
M[!is.finite(M)] <- 0 ; M[is.na(M)] <- 0
G1[!is.finite(G1)] <- 0 ; G1[is.na(G1)] <- 0
G2[!is.finite(G2)] <- 0 ; G2[is.na(G2)] <- 0
S[!is.finite(S)] <- 0 ; S[is.na(S)] <- 0

#' ### Correction to fit a (-3;3) range
M[M < -3] <- -3; M[M > 3] <- 3
G1[G1 < -3] <- -3; G1[G1 > 3] <- 3
G2[G2 < -3] <- -3; G2[G2 > 3] <- 3
S[S < -3] <- -3; S[S > 3] <- 3




#' ### Heatmap
heatmap.2(M,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

heatmap.2(M,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})

heatmap.2(G1,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

heatmap.2(G1,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})

heatmap.2(G2,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

heatmap.2(G2,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})

heatmap.2(S,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

heatmap.2(S,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})


#' ## Normalization by the T24
T24 <- counts[,8]
normT24 <- log2(counts/T24)

#' # Heatmap of the cell cycle related genes
#' ## Call the list of genes
CC <- read.csv("Menges2.csv", sep = ";", dec=",")

#' ## Separate the lists according to the phase
Mphase <- CC$AGI[CC$Phase=="M"]
G1phase <- CC$AGI[CC$Phase=="G1"]
G2phase <- CC$AGI[CC$Phase=="G2"]
Sphase <- CC$AGI[CC$Phase=="S"]

#' ## Split the normalized count list per phase
M <- as.matrix(normT24[rownames(normT24) %in% Mphase,8:13])
G1 <- as.matrix(normT24[rownames(normT24) %in% G1phase,8:13])
G2 <- as.matrix(normT24[rownames(normT24) %in% G2phase,8:13])
S <- as.matrix(normT24[rownames(normT24) %in% Sphase,8:13])

#' ## Correction of the tables
#' ### Removal of the infinte value
M[!is.finite(M)] <- 0 ; M[is.na(M)] <- 0
G1[!is.finite(G1)] <- 0 ; G1[is.na(G1)] <- 0
G2[!is.finite(G2)] <- 0 ; G2[is.na(G2)] <- 0
S[!is.finite(S)] <- 0 ; S[is.na(S)] <- 0

#' ### Correction to fit a (-3;3) range
M[M < -3] <- -3; M[M > 3] <- 3
G1[G1 < -3] <- -3; G1[G1 > 3] <- 3
G2[G2 < -3] <- -3; G2[G2 > 3] <- 3
S[S < -3] <- -3; S[S > 3] <- 3




#' ### Heatmap
heatmap.2(M,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

heatmap.2(M,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})

heatmap.2(G1,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

heatmap.2(G1,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})

heatmap.2(G2,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

heatmap.2(G2,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})

heatmap.2(S,
          col=hpal,
          dendrogram="none",
          trace="none",
          Colv=FALSE,
          Rowv=FALSE,
          margins=c(8,8))

heatmap.2(S,
          trace="none",
          col=hpal,
          Colv=FALSE,
          dendrogram = "row",
          margins=c(8,8),
          hclustfun = function(X){hclust(X,method="ward.D2")})




#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```

