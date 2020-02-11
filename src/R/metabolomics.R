library(tidyverse)
library(here)
library(matrixStats)

metabolomics <- read_csv2(here("data/metabolomics/metabolo.csv"),
                          col_types=cols(Sample=col_character(),
                                         .default=col_double())) %>% 
  select(-starts_with("N",ignore.case=FALSE)) %>% column_to_rownames("Sample")

                          
plot(density(as.matrix(log2(metabolomics))))

plot(density(rowSds(as.matrix(log2(metabolomics)))))

plot(density(colSds(as.matrix(log2(metabolomics)))))

ord <- order(colMeans(as.matrix(log2(metabolomics))),decreasing=FALSE)

plot(colSds(as.matrix(log2(metabolomics)))[ord],type="l")

lipids <- log2(lipids)

lipids[is.infinite(lipids)] <- 0

