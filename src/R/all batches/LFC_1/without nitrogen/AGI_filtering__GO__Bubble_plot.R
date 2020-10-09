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


#' ## Apply filters for the DEGs
#' ### For the 6hrs timepoint
res <- T6_NP_DMSO
cutoff1 <- res$res.log2FoldChange >= lfc & ! is.na(res$res.padj) & res$res.padj <= FDR
cutoff2 <- res$res.log2FoldChange <= -lfc & ! is.na(res$res.padj) & res$res.padj <= FDR
SucLow6hrs <- res$rownames.res.[cutoff2]
SucHigh6hrs <- res$rownames.res.[cutoff1]

res <- T6_NS_DMSO
cutoff1 <- res$res.log2FoldChange >= lfc & ! is.na(res$res.padj) & res$res.padj <= FDR
cutoff2 <- res$res.log2FoldChange <= -lfc & ! is.na(res$res.padj) & res$res.padj <= FDR
PiLow6hrs <- res$rownames.res.[cutoff2]
PiHigh6hrs <- res$rownames.res.[cutoff1]

res <- T6_NPS_AZD
cutoff1 <- res$res.log2FoldChange >= lfc & ! is.na(res$res.padj) & res$res.padj <= FDR
cutoff2 <- res$res.log2FoldChange <= -lfc & ! is.na(res$res.padj) & res$res.padj <= FDR
AzdLow6hrs <- res$rownames.res.[cutoff2]
AzdHigh6hrs <- res$rownames.res.[cutoff1]

res <- T6_NP_AZD
cutoff1 <- res$res.log2FoldChange >= lfc & ! is.na(res$res.padj) & res$res.padj <= FDR
cutoff2 <- res$res.log2FoldChange <= -lfc & ! is.na(res$res.padj) & res$res.padj <= FDR
AzdNPLow6hrs <- res$rownames.res.[cutoff2]
AzdNPHigh6hrs <- res$rownames.res.[cutoff1]

res <- T6_NS_AZD
cutoff1 <- res$res.log2FoldChange >= lfc & ! is.na(res$res.padj) & res$res.padj <= FDR
cutoff2 <- res$res.log2FoldChange <= -lfc & ! is.na(res$res.padj) & res$res.padj <= FDR
AzdNSLow6hrs <- res$rownames.res.[cutoff2]
AzdNSHigh6hrs <- res$rownames.res.[cutoff1]

#' ### For the 24hrs timepoint
res <- T24_NP_DMSO
cutoff1 <- res$res.log2FoldChange >= lfc & ! is.na(res$res.padj) & res$res.padj <= FDR
cutoff2 <- res$res.log2FoldChange <= -lfc & ! is.na(res$res.padj) & res$res.padj <= FDR
SucLow24hrs <- res$rownames.res.[cutoff2]
SucHigh24hrs <- res$rownames.res.[cutoff1]

res <- T24_NS_DMSO
cutoff1 <- res$res.log2FoldChange >= lfc & ! is.na(res$res.padj) & res$res.padj <= FDR
cutoff2 <- res$res.log2FoldChange <= -lfc & ! is.na(res$res.padj) & res$res.padj <= FDR
PiLow24hrs <- res$rownames.res.[cutoff2]
PiHigh24hrs <- res$rownames.res.[cutoff1]

res <- T24_NPS_AZD
cutoff1 <- res$res.log2FoldChange >= lfc & ! is.na(res$res.padj) & res$res.padj <= FDR
cutoff2 <- res$res.log2FoldChange <= -lfc & ! is.na(res$res.padj) & res$res.padj <= FDR
AzdLow24hrs <- res$rownames.res.[cutoff2]
AzdHigh24hrs <- res$rownames.res.[cutoff1]

res <- T24_NP_AZD
cutoff1 <- res$res.log2FoldChange >= lfc & ! is.na(res$res.padj) & res$res.padj <= FDR
cutoff2 <- res$res.log2FoldChange <= -lfc & ! is.na(res$res.padj) & res$res.padj <= FDR
AzdNPLow24hrs <- res$rownames.res.[cutoff2]
AzdNPHigh24hrs <- res$rownames.res.[cutoff1]

res <- T24_NS_AZD
cutoff1 <- res$res.log2FoldChange >= lfc & ! is.na(res$res.padj) & res$res.padj <= FDR
cutoff2 <- res$res.log2FoldChange <= -lfc & ! is.na(res$res.padj) & res$res.padj <= FDR
AzdNSLow24hrs <- res$rownames.res.[cutoff2]
AzdNSHigh24hrs <- res$rownames.res.[cutoff1]

#' ### For Wouter's data
res <- as.data.frame(Seedling)
cutoff1 <- res$res.log2FoldChange >= lfc & ! is.na(res$res.padj) & res$res.padj <= FDR
cutoff2 <- res$res.log2FoldChange <= -lfc & ! is.na(res$res.padj) & res$res.padj <= FDR
Seedling_Low <- res$rownames.res.[cutoff2]
Seedling_High <- res$rownames.res.[cutoff1]

#' # GO analysis
#' ## For the 6hr timepoint
GO <- gopher(SucLow6hrs,url="athaliana", alpha=2)
GO_SucLow6hrs <- GO$go

GO <- gopher(SucHigh6hrs,url="athaliana", alpha=2)
GO_SucHigh6hrs <- GO$go

GO <- gopher(PiLow6hrs, url="athaliana", alpha=2)
GO_PiLow6hrs <- GO$go

GO <- gopher(PiHigh6hrs, url="athaliana", alpha=2)
GO_PiHigh6hrs <- GO$go

GO <- gopher(AzdLow6hrs, url="athaliana", alpha=2)
GO_AzdLow6hrs <- GO$go

GO <- gopher(AzdHigh6hrs, url="athaliana", alpha=2)
GO_AzdHigh6hrs <- GO$go

GO <- gopher(AzdNPLow6hrs, url="athaliana", alpha=2)
GO_AzdNPLow6hrs <- GO$go

GO <- gopher(AzdNPHigh6hrs, url="athaliana", alpha=2)
GO_AzdNPHigh6hrs <- GO$go

GO <- gopher(AzdNSLow6hrs, url="athaliana", alpha=2)
GO_AzdNSLow6hrs <- GO$go

GO <- gopher(AzdNSHigh6hrs, url="athaliana", alpha=2)
GO_AzdNSHigh6hrs <- GO$go

#' ## For the 24hr timepoint
GO <- gopher(SucLow24hrs, url="athaliana", alpha=2)
GO_SucLow24hrs <- GO$go

GO <- gopher(SucHigh24hrs, url="athaliana", alpha=2)
GO_SucHigh24hrs <- GO$go

GO <- gopher(PiLow24hrs, url="athaliana", alpha=2)
GO_PiLow24hrs <- GO$go

GO <- gopher(PiHigh24hrs, url="athaliana", alpha=2)
GO_PiHigh24hrs <- GO$go

GO <- gopher(AzdLow24hrs, url="athaliana", alpha=2)
GO_AzdLow24hrs <- GO$go

GO <- gopher(AzdHigh24hrs, url="athaliana", alpha=2)
GO_AzdHigh24hrs <- GO$go

GO <- gopher(AzdNPLow24hrs, url="athaliana", alpha=2)
GO_AzdNPLow24hrs <- GO$go

GO <- gopher(AzdNPHigh24hrs, url="athaliana", alpha=2)
GO_AzdNPHigh24hrs <- GO$go

GO <- gopher(AzdNSLow24hrs, url="athaliana", alpha=2)
GO_AzdNSLow24hrs <- GO$go

GO <- gopher(AzdNSHigh24hrs, url="athaliana", alpha=2)
GO_AzdNSHigh24hrs <- GO$go

#' ## For Wouter's data
GO <- gopher(Seedling_Low, url="athaliana", alpha=2)
GO_Seedling_Low <- GO$go

GO <- gopher(Seedling_High, url="athaliana", alpha=2)
GO_Seedling_High <- GO$go

#' ## Add the experiment name in the GO tables
GO_AzdLow6hrs <- cbind(Description='GO_Azd6hrs', GO_AzdLow6hrs)
GO_AzdHigh6hrs <- cbind(Description='GO_Azd6hrs', GO_AzdHigh6hrs)
GO_SucHigh6hrs <- cbind(Description='GO_Suc6hrs', GO_SucHigh6hrs)
GO_SucLow6hrs <- cbind(Description='GO_Suc6hrs', GO_SucLow6hrs)
GO_PiHigh6hrs <- cbind(Description='GO_Pi6hrs', GO_PiHigh6hrs)
GO_PiLow6hrs <- cbind(Description='GO_Pi6hrs', GO_PiLow6hrs)
GO_AzdNPHigh6hrs <- cbind(Description='GO_AzdNP6hrs', GO_AzdNPHigh6hrs)
GO_AzdNPLow6hrs <- cbind(Description='GO_AzdNP6hrs', GO_AzdNPLow6hrs)
GO_AzdNSHigh6hrs <- cbind(Description='GO_AzdNS6hrs', GO_AzdNSHigh6hrs)
GO_AzdNSLow6hrs <- cbind(Description='GO_AzdNS6hrs', GO_AzdNSLow6hrs)
GO_AzdHigh24hrs <- cbind(Description='GO_Azd24hrs', GO_AzdHigh24hrs)
GO_AzdLow24hrs <- cbind(Description='GO_Azd24hrs', GO_AzdLow24hrs)
GO_SucHigh24hrs <- cbind(Description='GO_Suc24hrs', GO_SucHigh24hrs)
GO_SucLow24hrs <- cbind(Description='GO_Suc24hrs', GO_SucLow24hrs)
GO_PiHigh24hrs <- cbind(Description='GO_Pi24hrs', GO_PiHigh24hrs)
GO_PiLow24hrs <- cbind(Description='GO_Pi24hrs', GO_PiLow24hrs)
GO_AzdNPHigh24hrs <- cbind(Description='GO_AzdNP24hrs', GO_AzdNPHigh24hrs)
GO_AzdNPLow24hrs <- cbind(Description='GO_AzdNP24hrs', GO_AzdNPLow24hrs)
GO_AzdNSHigh24hrs <- cbind(Description='GO_AzdNS24hrs', GO_AzdNSHigh24hrs)
GO_AzdNSLow24hrs <- cbind(Description='GO_AzdNS24hrs', GO_AzdNSLow24hrs)
GO_Seedling_High <- cbind(Description='GO_Seedling', GO_Seedling_High)
GO_Seedling_Low <- cbind(Description='GO_Seedling', GO_Seedling_Low)

#' # Combine the GO tables
#' ## Prepare the tables
combined_up <- rbind.fill(GO_AzdHigh6hrs, GO_AzdNPHigh6hrs, GO_AzdNSHigh6hrs, GO_PiHigh6hrs, GO_SucHigh6hrs,GO_AzdHigh24hrs, GO_AzdNPHigh24hrs, GO_AzdNSHigh24hrs, GO_PiHigh24hrs, GO_SucHigh24hrs,GO_Seedling_High)
combined_down <- rbind.fill(GO_AzdLow6hrs,GO_SucLow6hrs,GO_PiLow6hrs,GO_AzdNPLow6hrs,GO_AzdNSLow6hrs,GO_AzdLow24hrs,GO_SucLow24hrs,GO_PiLow24hrs,GO_AzdNPLow24hrs,GO_AzdNSLow24hrs,GO_Seedling_Low)
#' ## Keep only the biological processes
combined_up <- combined_up[combined_up$namespace == "BP",]
combined_down <- combined_down[combined_down$namespace == "BP",]

#' ## Combining Up and Down-regulated genes together 
combined_up[,"Direction"]  <- "Induced"
combined_down[,"Direction"]  <- "Repressed"

combined_all_up_and_down <- rbind(combined_up, combined_down)

#' * Export the obtained data.frame containing also the non-significant GOs
write.table(combined_all_up_and_down,file="GOs_also_NON_sgnificant_GOs.txt", quote = F ,sep = "\t", row.names=FALSE)


#' ## Removing the root terms
exclude <- c("GO:0008150", "GO:0008283", "GO:0032502", "GO:0044848", "GO:0032501", "GO:0098743", "GO:0008152", "GO:0051704", "GO:0007610", "GO:0000003", "GO:0009987",
             'GO:0006791', 'GO:0001906', 'GO:0002376', 'GO:0006794', 'GO:0009758', 'GO:0015976', 'GO:0019740', 'GO:0022414', 
             'GO:0022610', 'GO:0023052', 'GO:0040007', 'GO:0040011', 'GO:0043473', 'GO:0048511', 'GO:0048518', 'GO:0048519', 'GO:0050789', 'GO:0050896', 'GO:0051179', 'GO:0065007', 'GO:0071840', 'GO:0098754', 'GO:0110148')
combined_all_up_and_down <- combined_all_up_and_down[! combined_all_up_and_down$id %in% exclude,]

# ' ## Change the p-values to change all the ones that are non-significant
combined_all_up_and_down$padj[combined_all_up_and_down$padj > 0.05] <- NA #Change this value as wished (normally 0.05)

#' # Create a list of interesting GO
#' * Create a list of cell cycle related GO enriched in Repressed genes after 24hrs of Pi starvation
C_Cycle <- c("GO:0007049","GO:0051726","GO:0044770","GO:0007051","GO:0007059","GO:0042549","GO:0051301","GO:0051304",
             "GO:0009061","GO:0007017","GO:0046490","GO:0006720","GO:0032196")

#' * Pyrimidine-containing compound metabolism
Pyrimidine <- c("GO:0072527",'GO:0006304','GO:0055086','GO:0018193','GO:0016072','GO:0090501','GO:0006354',
                'GO:1901564','GO:1901566','GO:0006139','GO:0000375','GO:0043603','GO:0032446','GO:0019318','GO:0006259',
                'GO:0006260','GO:0090305','GO:0006518','GO:0006553','GO:0006090','GO:0006089',
                'GO:0008213','GO:0006396','GO:0034641','GO:0005991','GO:0009451','GO:0034660')

#' * DNA replication
Replication <- c('GO:0006260','GO:0006304','GO:0090066','GO:0055086','GO:0018193','GO:0048519','GO:0006793',
                 'GO:0031323','GO:0045814','GO:0006270','GO:0006275','GO:0018205','GO:0016072','GO:0044260',
                 'GO:0090501','GO:0009069','GO:0050789','GO:1901576','GO:0046483','GO:1901564','GO:0051338',
                 'GO:0006139','GO:0043603','GO:0032446','GO:0006725','GO:0019318','GO:0006259',
                 'GO:0010498','GO:0006518','GO:0010467','GO:0033559','GO:0006090','GO:0043412','GO:0008213',
                 'GO:0060249','GO:0040029','GO:0009892','GO:0040008','GO:0072527','GO:0034641',
                 'GO:0005991','GO:0009451','GO:0034660')

#' * Lipopolysaccharide metabolism (repressed by sugar starvation after 24 hrs)
LPS <- c("GO:1903509","GO:1901137","GO:0006643","GO:0006629")

#' * Response to heat (repressed by sugar starvation after 24 hrs)
heat <- c("GO:0009408","GO:0071496","GO:0009642","GO:0016036","GO:0006950","GO:0006979","GO:1901700","GO:0002215")

#' * ncRNA (induced by 24hrs of sugar starvation)
ncRNA <- c('GO:0034660','GO:0042537','GO:0009768','GO:0055070','GO:0043900','GO:0051246','GO:0018198','GO:0017014',
           'GO:0016054','GO:0016070','GO:0016072','GO:0042592','GO:0022900','GO:0055114','GO:0006081','GO:0006082',
           'GO:0008219','GO:0009069','GO:0009308','GO:0001505','GO:0006396','GO:1901605','GO:0030258','GO:0044283',
           'GO:0009850','GO:0044281','GO:0005991','GO:0005997','GO:0042430','GO:0051174','GO:0065008','GO:1901657','GO:0006766')


myList <- list(C_Cycle,Pyrimidine,Replication,LPS, heat, ncRNA)

#' ## List of GO groups overrepresented among genes induced after 24hrs of Suc starvation
aging <- c("GO:0007568","GO:0090558","GO:0048366","GO:0010618","GO:0007275")
response_to_red_light <- c("GO:0010114","GO:1901700","GO:0009743","GO:0006290","GO:0080021","GO:0017085","GO:0043617","GO:0034285","GO:0009637","GO:0009646","GO:0009611","GO:0010033")
photosynthesis <- "GO:0015979"
response_to_stimulus <- "GO:0050896"
localization <- "GO:0051179"
multi_organism_process <- "GO:0051704"
biological_regulation <- "GO:0065007"
catabolism <- "GO:0009056"
plastid_organization <- c("GO:0009657","GO:0043933","GO:0072657","GO:0022613","GO:0008643","GO:0044085","GO:0010256")
isoprenoid_transport <- c("GO:0046864","GO:0055085","GO:0006811","GO:0015850","GO:0032879","GO:0015853","GO:0072511","GO:0010232","GO:0015672","GO:0015718")
drug_metabolism <- "GO:0017144"
multidimensional_cell_growth <- "GO:0009825"
generation_of_precursor_metabolites_and_energy <- "GO:0006091"
secondary_metabolism <- c("GO:0019748","GO:0015977")
ROS_metabolism <- "GO:0072593"
sulfur_compound_biosynthesis <- "GO:0044272"
dephosphorylation <- c("GO:0016311","GO:0006733","GO:0006739")
sulfur_compound_metabolism <- "GO:0006790"
ncRNA_metabolism <- c("GO:0034660","GO:0042537","GO:0009768","GO:0055070","GO:0043900","GO:0006520","GO:0051246","GO:0018198","GO:0017014","GO:0048583","GO:0016070","GO:0016072","GO:0042592","GO:0033559","GO:0022900","GO:0055114","GO:0006081","GO:0006082","GO:0008219","GO:0009311","GO:0009069","GO:0009308","GO:0001505","GO:0006396","GO:0009850","GO:0044283","GO:0044281","GO:0009607","GO:0005991","GO:0044093","GO:0005997","GO:0042430","GO:0051174","GO:0065008","GO:0042221","GO:0006766")
cofactor_metabolism <- c("GO:0051186","GO:0006793")
carbohydrate_metabolism <- "GO:0005975"
ureide_metabolism <- "GO:0010135"
oxidative_photosynthetic_carbon_pathway <- "GO:0009854"


myList <- list(
    aging = aging,
    response_to_red_light= response_to_red_light,
    photosynthesis=photosynthesis,
    catabolism=catabolism,
    plastid_organization=plastid_organization,
    isoprenoid_transport=isoprenoid_transport,
    drug_metabolism=drug_metabolism,
    multidimensional_cell_growth=multidimensional_cell_growth,
    generation_of_precursor_metabolites_and_energy=generation_of_precursor_metabolites_and_energy,
    secondary_metabolism=secondary_metabolism,
    ROS_metabolism=ROS_metabolism,
    sulfur_compound_biosynthesis=sulfur_compound_biosynthesis,
    dephosphorylation=dephosphorylation,
    sulfur_compound_metabolism=sulfur_compound_metabolism,
    ncRNA_metabolism=ncRNA_metabolism,
    cofactor_metabolism=cofactor_metabolism,
    carbohydrate_metabolism=carbohydrate_metabolism,
    ureide_metabolism=ureide_metabolism,
    oxidative_photosynthetic_carbon_pathway=oxidative_photosynthetic_carbon_pathway)

#' ## List of GO groups overrepresented under the different conditions
LowSuc6hr <- GO_SucLow6hrs$id[GO_SucLow6hrs$padj <= 0.05]
LowPi6hr <- GO_PiLow6hrs$id[GO_PiLow6hrs$padj <= 0.05]
LowAzd6hr <- GO_AzdLow6hrs$id[GO_AzdLow6hrs$padj <= 0.05]
LowAzdSuc6hr <- GO_AzdNPLow6hrs$id[GO_AzdNPLow6hrs$padj <= 0.05]
LowAzdPi6hr <- GO_AzdNSLow6hrs$id[GO_AzdNSLow6hrs$padj <= 0.05]

HighSuc6hr <- GO_SucHigh6hrs$id[GO_SucHigh6hrs$padj <= 0.05]
HighPi6hr <- GO_PiHigh6hrs$id[GO_PiHigh6hrs$padj <= 0.05]
HighAzd6hr <- GO_AzdHigh6hrs$id[GO_AzdHigh6hrs$padj <= 0.05]
HighAzdSuc6hr <- GO_AzdNPHigh6hrs$id[GO_AzdNPHigh6hrs$padj <= 0.05]
HighAzdPi6hr <- GO_AzdNSHigh6hrs$id[GO_AzdNSHigh6hrs$padj <= 0.05]

LowSuc24hr <- GO_SucLow24hrs$id[GO_SucLow24hrs$padj <= 0.05]
LowPi24hr <- GO_PiLow24hrs$id[GO_PiLow24hrs$padj <= 0.05]
LowAzd24hr <- GO_AzdLow24hrs$id[GO_AzdLow24hrs$padj <= 0.05]
LowAzdSuc24hr <- GO_AzdNPLow24hrs$id[GO_AzdNPLow24hrs$padj <= 0.05]
LowAzdPi24hr <- GO_AzdNSLow24hrs$id[GO_AzdNSLow24hrs$padj <= 0.05]

HighSuc24hr <- GO_SucHigh24hrs$id[GO_SucHigh24hrs$padj <= 0.05]
HighPi24hr <- GO_PiHigh24hrs$id[GO_PiHigh24hrs$padj <= 0.05]
HighAzd24hr <- GO_AzdHigh24hrs$id[GO_AzdHigh24hrs$padj <= 0.05]
HighAzdSuc24hr <- GO_AzdNPHigh24hrs$id[GO_AzdNPHigh24hrs$padj <= 0.05]
HighAzdPi24hr <- GO_AzdNSHigh24hrs$id[GO_AzdNSHigh24hrs$padj <= 0.05]

LowSeedling <- GO_Seedling_Low$id[GO_Seedling_Low$padj <= 0.05]
HighSeedling <- GO_Seedling_High$id[GO_Seedling_High$padj <= 0.05]

myList <- list(LowSuc6hr = LowSuc6hr,
               LowPi6hr = LowPi6hr,
               LowAzd6hr= LowAzd6hr,
               LowAzdSuc6hr= LowAzdSuc6hr,
               LowAzdPi6hr= LowAzdPi6hr,
               HighSuc6hr= HighSuc6hr,
               HighPi6hr= HighPi6hr,
               HighAzd6hr= HighAzd6hr,
               HighAzdSuc6hr= HighAzdSuc6hr,
               HighAzdPi6hr= HighAzdPi6hr,
               LowSuc24hr= LowSuc24hr,
               LowPi24hr= LowPi24hr,
               LowAzd24hr= LowAzd24hr,
               LowAzdSuc24hr= LowAzdSuc24hr,
               LowAzdPi24hr= LowAzdPi24hr,
               HighSuc24hr= HighSuc24hr,
               HighPi24hr= HighPi24hr,
               HighAzd24hr= HighAzd24hr,
               HighAzdSuc24hr= HighAzdSuc24hr,
               HighAzdPi24hr= HighAzdPi24hr,
               LowSeedling = LowSeedling,
               HighSeedling=HighSeedling
)

#' ## List of interesting GOs found in the genes specifically induced after Pi+AZD but not found deregulated after Pi only treatment
Phosphate_N_lipid <- c("GO:0006629","GO:0006793","GO:0006796","GO:0016310","GO:0006468","GO:0044255")

#' * Make the GO plot
go <- combined_all_up_and_down[combined_all_up_and_down$id %in% Phosphate_N_lipid,]
ggplot(go,aes(x=Description,y=id))+ geom_point(aes(size=npat, color=-log10(padj)))+ 
    facet_wrap(~Direction, nrow=1) +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("phosphate and lipid")+scale_shape_discrete(solid=F)+
    scale_colour_gradientn(colours=c("red2", "purple2", "blue2")) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

#' # Create a sublist of GO with only the GOs of interest
go <- combined_all_up_and_down[combined_all_up_and_down$id %in% myList,]
go <- lapply(myList,function(x)
{
    combined_all_up_and_down[combined_all_up_and_down$id %in% x,]
})

#' # Make the bubble plot
lapply(names(go),function(nam,l)
{
    x <- l[[nam]]
    ggplot(x,aes(x=Description,y=id))+ geom_point(aes(size=npat, color=-log10(padj)))+ 
        facet_wrap(~Direction, nrow=1) +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        ggtitle(nam)+scale_shape_discrete(solid=F)+
        scale_colour_gradientn(colours=c("red2", "purple2", "blue2")) +
        theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
}, go)


#' # Identification of genes differentially deregulated specifically after Pi starvation or specifically after the combination of Pi starvation and AZD treatment
PiSpecLow6hrs <-setdiff(PiLow6hrs,AzdNSLow6hrs)
AzdSpecLow6hrs <- setdiff(AzdNSLow6hrs,PiLow6hrs)
PiSpecLow24hrs <-setdiff(PiLow24hrs,AzdNSLow24hrs)
AzdSpecLow24hrs <- setdiff(AzdNSLow24hrs,PiLow24hrs)

PiSpecHigh6hrs <-setdiff(PiHigh6hrs,AzdNSHigh6hrs)
AzdSpecHigh6hrs <- setdiff(AzdNSHigh6hrs,PiHigh6hrs)
PiSpecHigh24hrs <-setdiff(PiHigh24hrs,AzdNSHigh24hrs)
AzdSpecHigh24hrs <- setdiff(AzdNSHigh24hrs,PiHigh24hrs)


View(PiSpecLow6hrs); View(PiSpecLow24hrs); View(PiSpecHigh6hrs); View(PiSpecHigh24hrs)
View(AzdSpecLow6hrs); View(AzdSpecLow24hrs); View(AzdSpecHigh6hrs); View(AzdSpecHigh24hrs)

write.table(PiSpecLow6hrs, "/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/PiSpecLow6hrs.csv", sep="\t")
write.table(PiSpecLow24hrs, "/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/PiSpecLow24hrs.csv", sep="\t")
write.table(PiSpecHigh6hrs, "/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/PiSpecHigh6hrs.csv", sep="\t")
write.table(PiSpecHigh24hrs, "/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/PiSpecHigh24hrs.csv", sep="\t")

write.table(AzdSpecLow6hrs, "/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/AzdSpecLow6hrs.csv", sep="\t")
write.table(AzdSpecLow24hrs, "/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/AzdSpecLow24hrs.csv", sep="\t")
write.table(AzdSpecHigh6hrs, "/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/AzdSpecHigh6hrs.csv", sep="\t")
write.table(AzdSpecHigh24hrs, "/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/AzdSpecHigh24hrs.csv", sep="\t")

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```

