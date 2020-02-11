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
vst <- read_tsv(here("data/analysis/salmon/model-aware_variance-stabilised-data.tsv"),
                col_types=cols(.default=col_double(),rowname=col_character())) %>% 
  column_to_rownames("rowname")

samples <- read_csv(here("doc/droughtneedles.csv"),
                    col_types=cols(
                      col_character(),
                      col_character(),
                      col_factor()
                    )) %>% 
  mutate(SampleID=sub("_.*","",SciLifeID))

stopifnot(all(samples$SampleID == colnames(vst)))

#' # Filter
sels <- rangeFeatureSelect(counts=as.matrix(vst),
                           conditions=samples$Level,
                           nrep=3)

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
