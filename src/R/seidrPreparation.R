#' ---
#' title: "Somatic Embryogenesis Network data preparation"
#' author: "Iryna Shutava and Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
})

#' Helper
source(here("UPSCb-common/src/R/featureSelection.R"))

#' # Data
vst <- read.csv(here("data/seidr/ForSeidr.csv"))
rownames(vst) <- vst$X
vst <- vst[,2:length(vst)]

samples <- read.csv("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/analysis_Tom/samplelist.csv")


stopifnot(all(samples$SampleID == colnames(vst)))

#' # Filter
sels <- rangeFeatureSelect(counts=as.matrix(vst),
                           conditions=samples$Conditions,
                           nrep=2)

vst.cutoff <- 1

#' # Export
dir.create(here("data/analysis/seidr"),showWarnings=FALSE)

#' * gene by column, without names matrix
write.table(t(vst[sels[[vst.cutoff+1]],]),
            file=here("data/analysis/seidr/headless.tsv"),
            col.names=FALSE,
            row.names=FALSE,
            sep="\t",quote=FALSE)

#' * gene names, one row
write.table(t(sub("\\.1$","",rownames(vst)[sels[[vst.cutoff+1]]])),
            file=here("data/analysis/seidr/genes.tsv"),
            col.names=FALSE,
            row.names=FALSE,
            sep="\t",quote=FALSE)

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
