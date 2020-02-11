#' ---
#' title: "Biological QA"
#' author: "Thomas Dobrenel (adapted from Iryna Shutava and Nicolas Delhomme to only consider a subset of the dataset)"
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
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(vsn))

#' Source some helper functions
source("~/Git/UPSCb/src/R/plot.multidensity.R")
source("~/Git/UPSCb/src/R/featureSelection.R")

#' Create palettes
pal <- brewer.pal(8,"Dark2")
#pal12 <- brewer.pal(12,"Paired")
#tol15rainbow=c("#114477", "#4477AA", "#77AADD", "#117755", "#44AA88", "#99CCBB", "#777711", "#AAAA44", "#DDDD77", "#771111", "#AA4444", "#DD7777", "#771144", "#AA4477", "#DD77AA")
hpal <- colorRampPalette(c("blue","white","red"))(100)

#' Register the default plot margin
mar <- par("mar")

#' # Raw data
#' ## Loading
#' Read the sample information
samples <- read.csv("~/Git/UPSCb/projects/arabidopsis-nutrition-TOR/doc/samples.csv")

#' ### Original data
orig <- list.files("Salmon", 
                    recursive = TRUE, 
                    pattern = "quant.sf",
                    full.names = TRUE)
orig <- orig[!grepl("unpaired",orig)]

#' name them
names(orig) <- sub("_S.*.","",sapply(strsplit(orig, "/"), .subset, 2))

#' Reorder the sample data.frame according to the way
#' the results were read in and accounting for tech reps
samples <- samples[match(names(orig),samples$SciLifeID),]

#' Remove samples from batch 1
samples <- samples[(factor(substr(sub(".*_","",samples$SciLifeID),1,1)))==2,]
orig <- orig[match(samples$SciLifeID,names(orig))]


#' Read the expression at the transcript level
tx <- suppressMessages(tximport(files = orig, 
                                type = "salmon", 
                                txOut = TRUE))
                                #countsFromAbundance = "scaledTPM"))



#' # Data normalisation 
#'  For visualization, the data is
#' submitted to a variance stabilization
#' transformation using DESeq2. The 
#' dispersion is estimated independently
#' of the sample type and sample experiment

sel <- colnames(tx$counts) %in% samples$SciLifeID
kt <- round(tx$counts[,sel])

#' summarise to genes
tx2gene <- data.frame(TXID=rownames(kt),
                      GENEID=sub("\\.[0-9]+","",rownames(kt)))

# Variant for the lefting of the splicing forms of genes
#tx2gene <- data.frame(TXID=rownames(kt),
#                      GENEID=rownames(kt))

gx <- summarizeToGene(tx,tx2gene=tx2gene)

kg <- round(gx$counts[,sel]) 

#' ## Raw data export

dir.create(file.path("analysis_Tom","Salmon"),showWarnings=FALSE,recursive=TRUE)
write.table(kg,file="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/analysis_Tom/batch2-unormalised-gene-expression_data.csv")
save(kg, samples, file = "/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/analysis_Tom/batch2_counts.rda")


#' # Data normalisation 
#'  For visualization, the data is
#' submitted to a variance stabilization
#' transformation using DESeq2. The 
#' dispersion is estimated independently
#' of the sample type 
#' TRANSCRIPTS


conditions <- as.integer(samples$Timepoint) + as.integer(samples$Nutrition) + as.integer(samples$AZD)

dds <- DESeqDataSetFromMatrix(
  countData = kt,
  colData = data.frame(condition=conditions),
  design = ~ condition)


#' Check the size factors (i.e. the sequencing library size effect)
#' There is no big variation, a Variance Stabilizing Transformation can
#' be used (over a Relative Log Transformation)
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
names(sizes) <- colnames(kt)
sizes
boxplot(sizes, main="Sequencing libraries size factor")

#' Perform the VST
colData(dds)$condition <- factor(colData(dds)$condition,
                                 levels=unique(conditions))
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
colnames(vst) <- colnames(kt)
vst <- vst - min(vst)

#' # Data normalisation 
#'  For visualization, the data is
#' submitted to a variance stabilization
#' transformation using DESeq2. The 
#' dispersion is estimated independently
#' of the sample type 
#' GENES


conditions <- as.integer(samples$Timepoint) + as.integer(samples$Nutrition) + as.integer(samples$AZD)

dds <- DESeqDataSetFromMatrix(
    countData = kg,
    colData = data.frame(condition=conditions),
    design = ~ condition)


#' Check the size factors (i.e. the sequencing library size effect)
#' There is no big variation, a Variance Stabilizing Transformation can
#' be used (over a Relative Log Transformation)
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
names(sizes) <- colnames(kg)
sizes
boxplot(sizes, main="Sequencing libraries size factor")

#' Perform the VST
colData(dds)$condition <- factor(colData(dds)$condition,
                                 levels=unique(conditions))
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
colnames(vst) <- colnames(kg)
vst <- vst - min(vst)


write.csv(vst,"/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/analysis_Tom/library-size-normalized_variance-stabilized_data_batch2.csv")

#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```

