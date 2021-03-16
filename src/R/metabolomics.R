library(tidyverse)
library(here)
library(matrixStats)

#' # Preparation of the metabolomics data
#' ## Load the metabolomics data results
metabolomics <- read_csv2(here("data/metabolomics/metabolo.csv"),
                          col_types=cols(Sample=col_character(),
                                         .default=col_double())) %>% 
  select(-starts_with("N",ignore.case=FALSE)) %>% column_to_rownames("Sample")

#' ## Verification if a variance stabilization is needed
plot(density(as.matrix(log2(metabolomics))))

plot(density(rowSds(as.matrix(log2(metabolomics)))))

plot(density(colSds(as.matrix(log2(metabolomics)))))

ord <- order(colMeans(as.matrix(log2(metabolomics))),decreasing=FALSE)

plot(colSds(as.matrix(log2(metabolomics)))[ord],type="l")

#' ## Associate the data.matrix with the sample names
#' * Load the sample list
samples_metabo <- read.csv(here("doc/samples_metabolo.csv"), sep=";")
samples_metabo <- cbind(samples_metabo,
                        Tp_Id = factor(paste(samples_metabo$Tp,
                                      samples_metabo$Id,
                                      sep=".")))
#' * Verification step
stopifnot(all(rownames(metabolomics) == samples_metabo$Sample))
#' * Change the names of the samples
rownames(metabolomics) <- samples_metabo$Tp_Id



#' # Preparation of the lipidomics data
#' ## Load the lipidomics data results

lipids <- read.csv2(here("data/lipidomics/lipid.csv"),dec=".",row.names = 1) 
lipids <- t(lipids)
lipids[is.na(lipids)] <- 0
lipids <- lipids +1

#' ## Verification if a variance stabilization is needed                  
plot(density(as.matrix(log2(lipids))))

plot(density(rowSds(as.matrix(log2(lipids)))))

plot(density(colSds(as.matrix(log2(lipids)))))

ord <- order(colMeans(as.matrix(log2(lipids))),decreasing=FALSE)

plot(colSds(as.matrix(log2(lipids)))[ord],type="l")

#' ## Associate the data.matrix with the sample names
#' * Load the sample list
samples_lipido <- read.csv(here("doc/samples_lipido.csv"), sep=";")
samples_lipido <- cbind(samples_lipido,
                        Tp_Id = factor(paste(samples_lipido$Timepoint,
                                             samples_lipido$Id,
                                             sep=".")))
#' * Remove the samples from the SC experiment
samples_lipido <- samples_lipido[(!grepl("SC",samples_lipido$Sample_Id)),]
lipids <- lipids[match(samples_lipido$Tube,rownames(lipids)),]


#' * Verification step
stopifnot(all(rownames(lipids) == samples_lipido$Tube))
#' * Change the names of the samples
rownames(lipids) <- samples_lipido$Tp_Id

#' # Preparation of the transcriptomic data
transcripto <- read.csv("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/analysis_Tom/library-size-normalized_variance-stabilized_data_nutrition.csv")

#' # Selection of the samples found in the transcriptomic analysis
#' * For the metabolomics data
colnames(transcripto)
rownames(metabolomics)
metabolomics <- metabolomics[rownames(metabolomics) %in% colnames(transcripto),]

#' * For the lipidomics data
colnames(transcripto)
rownames(lipids)
lipids <- lipids[rownames(lipids) %in% colnames(transcripto),]






#' # Creation of the table for SEIDR
#' ## Prepare the transcriptomics data table
rownames(transcripto) <- transcripto$X
transcripto <- transcripto[2:length(transcripto)]
transcripto <- transcripto[,sort(colnames(transcripto))]

#' ## Prepare the metabolomics data table
metabolomics <- log2(metabolomics)
metabolomics <- t(metabolomics)
metabolomics <- metabolomics[,sort(colnames(metabolomics))]

#' ## Prepare the lipidomics data table
lipids <- log2(lipids)
lipids <- t(lipids)
lipids <- lipids[,sort(colnames(lipids))]
#lipids[is.infinite(lipids)] <- 0

#' ## Verification
colnames(transcripto)
colnames(metabolomics)
colnames(lipids)
stopifnot(all(colnames(transcripto) == colnames(metabolomics)))
stopifnot(all(colnames(transcripto) == colnames(lipids)))

#' ## Preparation of the table
table <- rbind(transcripto,metabolomics,lipids)
dir.create(here("data/seidr"),showWarnings=FALSE,recursive=TRUE)
write.csv(table,here("data/seidr/ForSeidr.csv"),row.names=TRUE)

