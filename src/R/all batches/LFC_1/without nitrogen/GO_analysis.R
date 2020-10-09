#' ---
#' title: "Creation of GO bubble plots"
#' author: "Thomas Dobrenel and Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup

#' * Set the working dir
setwd("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/GO")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/GO")
#' ```

#' * Load libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(plyr); library(dplyr) # Library 'plyr' should be loaded before library 'dplyr'
})

#' * Source some helper functions
suppressPackageStartupMessages({
    source("~/Git/UPSCb/UPSCb-common/src/R/gopher.R")
})


#' * Create palettes


#' * Register the default plot margin
mar <- par("mar")

#' # Load the Gopher outputs
#' ## Analyzed by conditions in contrast to the T0
#' ### 6hrs of treatment

#' * DMSO samples
GO_T6NPSDMSO_High <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/T0/GO_T6NPSDMSO_High.csv", sep= ";", head=T)
GO_T6NPSDMSO_Low <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/T0/GO_T6NPSDMSO_Low.csv", sep= ";", head=T)
GO_T6NSDMSO_High <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/T0/GO_T6NSDMSO_High.csv", sep= ";", head=T)
GO_T6NSDMSO_Low <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/T0/GO_T6NSDMSO_Low.csv", sep= ";", head=T)
GO_T6NPDMSO_High <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/T0/GO_T6NPDMSO_High.csv", sep= ";", head=T)
GO_T6NPDMSO_Low <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/T0/GO_T6NPDMSO_Low.csv", sep= ";", head=T)

#' * AZD samples
GO_T6NPSAZD_High <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/T0/GO_T6NPSAZD_High.csv", sep= ";", head=T)
GO_T6NPSAZD_Low <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/T0/GO_T6NPSAZD_Low.csv", sep= ";", head=T)
GO_T6NSAZD_High <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/T0/GO_T6NSAZD_High.csv", sep= ";", head=T)
GO_T6NSAZD_Low <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/T0/GO_T6NSAZD_Low.csv", sep= ";", head=T)
GO_T6NPAZD_High <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/T0/GO_T6NPAZD_High.csv", sep= ";", head=T)
GO_T6NPAZD_Low <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/T0/GO_T6NPAZD_Low.csv", sep= ";", head=T)

#' ### 24hrs of treatment

#' * DMSO samples
GO_T24NPSDMSO_High <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/T0/GO_T24NPSDMSO_High.csv", sep= ";", head=T)
GO_T24NPSDMSO_Low <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/T0/GO_T24NPSDMSO_Low.csv", sep= ";", head=T)
GO_T24NSDMSO_High <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/T0/GO_T24NSDMSO_High.csv", sep= ";", head=T)
GO_T24NSDMSO_Low <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/T0/GO_T24NSDMSO_Low.csv", sep= ";", head=T)
GO_T24NPDMSO_High <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/T0/GO_T24NPDMSO_High.csv", sep= ";", head=T)
GO_T24NPDMSO_Low <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/T0/GO_T24NPDMSO_Low.csv", sep= ";", head=T)

#' * AZD samples
GO_T24NPSAZD_High <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/T0/GO_T24NPSAZD_High.csv", sep= ";", head=T)
GO_T24NPSAZD_Low <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/T0/GO_T24NPSAZD_Low.csv", sep= ";", head=T)
GO_T24NSAZD_High <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/T0/GO_T24NSAZD_High.csv", sep= ";", head=T)
GO_T24NSAZD_Low <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/T0/GO_T24NSAZD_Low.csv", sep= ";", head=T)
GO_T24NPAZD_High <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/T0/GO_T24NPAZD_High.csv", sep= ";", head=T)
GO_T24NPAZD_Low <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/T0/GO_T24NPAZD_Low.csv", sep= ";", head=T)

#' ## Analyzed by conditions in contrast to the respective timepoints
#' ### 6hrs of treatment

#' * AZD treatment
GO_AzdHigh6hrs <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/GO_AzdHigh6hrs.csv", sep= ";", head=T)
GO_AzdLow6hrs <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/GO_AzdLow6hrs.csv", sep= ";", head=T)

#' * Carbon starvation
GO_SucHigh6hrs <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/GO_SucHigh6hrs.csv", sep= ";", head=T)
GO_SucLow6hrs <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/GO_SucLow6hrs.csv", sep= ";", head=T)

#' * Phosphorus starvation
GO_PiHigh6hrs <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/GO_PiHigh6hrs.csv", sep= ";", head=T)
GO_PiLow6hrs <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/GO_PiLow6hrs.csv", sep= ";", head=T)

#' * Carbon starvation + AZD treatment
GO_AzdNPHigh6hrs <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/GO_AzdNPHigh6hrs.csv", sep= ";", head=T)
GO_AzdNPLow6hrs <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/GO_AzdNPLow6hrs.csv", sep= ";", head=T)

#' * Phosphorus starvation + AZD treatment
GO_AzdNSHigh6hrs <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/GO_AzdNSHigh6hrs.csv", sep= ";", head=T)
GO_AzdNSLow6hrs <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/GO_AzdNSLow6hrs.csv", sep= ";", head=T)

#' ### 24hrs of treatment

#' * AZD treatment
GO_AzdHigh24hrs <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/GO_AzdHigh24hrs.csv", sep= ";", head=T)
GO_AzdLow24hrs <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/GO_AzdLow24hrs.csv", sep= ";", head=T)

#' * Carbon starvation
GO_SucHigh24hrs <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/GO_SucHigh24hrs.csv", sep= ";", head=T)
GO_SucLow24hrs <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/GO_SucLow24hrs.csv", sep= ";", head=T)

#' * Phosphorus starvation
GO_PiHigh24hrs <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/GO_PiHigh24hrs.csv", sep= ";", head=T)
GO_PiLow24hrs <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/GO_PiLow24hrs.csv", sep= ";", head=T)

#' * Carbon starvation + AZD treatment
GO_AzdNPHigh24hrs <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/GO_AzdNPHigh24hrs.csv", sep= ";", head=T)
GO_AzdNPLow24hrs <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/GO_AzdNPLow24hrs.csv", sep= ";", head=T)

#' * Phosphorus starvation + AZD treatment
GO_AzdNSHigh24hrs <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/GO_AzdNSHigh24hrs.csv", sep= ";", head=T)
GO_AzdNSLow24hrs <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Conditions/GO_AzdNSLow24hrs.csv", sep= ";", head=T)

#' ## Analyzed by nutrition*AZD in contrast to the respective timepoints
#' ### 6hrs of treatment

#' * Pi.starv_AZD vs. Pi.starv + AZD
GO_Azd_NSHigh6hrs <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Nutrition_AZD/GO_Azd_NSHigh6hrs.csv", sep= ";", head=T)
GO_Azd_NSLow6hrs <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Nutrition_AZD/GO_Azd_NSHigh6hrs.csv", sep= ";", head=T)

#' * Suc.starv_AZD vs. Suc.starv + AZD
GO_Azd_NPHigh6hrs <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Nutrition_AZD/GO_Azd_NPHigh6hrs.csv", sep= ";", head=T)
GO_Azd_NPLow6hrs <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Nutrition_AZD/GO_Azd_NPHigh6hrs.csv", sep= ";", head=T)

#' ### 24hrs of treatment

#' * Pi.starv_AZD vs. Pi.starv + AZD
GO_Azd_NSHigh24hrs <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Nutrition_AZD/GO_Azd_NSHigh24hrs.csv", sep= ";", head=T)
GO_Azd_NSLow24hrs <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Nutrition_AZD/GO_Azd_NSHigh24hrs.csv", sep= ";", head=T)

#' * Suc.starv_AZD vs. Suc.starv + AZD
GO_Azd_NPHigh24hrs <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Nutrition_AZD/GO_Azd_NPHigh24hrs.csv", sep= ";", head=T)
GO_Azd_NPLow24hrs <- read.table("/mnt/picea/projects/arabidopsis/jhanson/arabidopsis-nutrition-TOR/Nutrition_AZD/GO_Azd_NPHigh24hrs.csv", sep= ";", head=T)





#' # Add the headers as a description in the tables
ls()[substr(ls(),1,3) == "GO_"]


GO_T6NPSDMSO_High <- cbind(Description='GO_T6NPSDMSO', GO_T6NPSDMSO_High)
GO_T6NPSDMSO_Low <- cbind(Description='GO_T6NPSDMSO', GO_T6NPSDMSO_Low)
GO_T6NSDMSO_High <- cbind(Description='GO_T6NSDMSO', GO_T6NSDMSO_High)
GO_T6NSDMSO_Low <- cbind(Description='GO_T6NSDMSO', GO_T6NSDMSO_Low)
GO_T6NPDMSO_High <- cbind(Description='GO_T6NPDMSO', GO_T6NPDMSO_High)
GO_T6NPDMSO_Low <- cbind(Description='GO_T6NPDMSO', GO_T6NPDMSO_Low)
GO_T6NPSAZD_High <- cbind(Description='GO_T6NPSAZD', GO_T6NPSAZD_High)
GO_T6NPSAZD_Low <- cbind(Description='GO_T6NPSAZD', GO_T6NPSAZD_Low)
GO_T6NSAZD_High <- cbind(Description='GO_T6NSAZD', GO_T6NSAZD_High)
GO_T6NSAZD_Low <- cbind(Description='GO_T6NSAZD', GO_T6NSAZD_Low)
GO_T6NPAZD_High <- cbind(Description='GO_T6NPAZD', GO_T6NPAZD_High)
GO_T6NPAZD_Low <- cbind(Description='GO_T6NPAZD', GO_T6NPAZD_Low)
GO_T24NPSDMSO_High <- cbind(Description='GO_T24NPSDMSO', GO_T24NPSDMSO_High)
GO_T24NPSDMSO_Low <- cbind(Description='GO_T24NPSDMSO', GO_T24NPSDMSO_Low)
GO_T24NSDMSO_High <- cbind(Description='GO_T24NSDMSO', GO_T24NSDMSO_High)
GO_T24NSDMSO_Low <- cbind(Description='GO_T24NSDMSO', GO_T24NSDMSO_Low)
GO_T24NPDMSO_High <- cbind(Description='GO_T24NPDMSO', GO_T24NPDMSO_High)
GO_T24NPDMSO_Low <- cbind(Description='GO_T24NPDMSO', GO_T24NPDMSO_Low)
GO_T24NPSAZD_High <- cbind(Description='GO_T24NPSAZD', GO_T24NPSAZD_High)
GO_T24NPSAZD_Low <- cbind(Description='GO_T24NPSAZD', GO_T24NPSAZD_Low)
GO_T24NSAZD_High <- cbind(Description='GO_T24NSAZD', GO_T24NSAZD_High)
GO_T24NSAZD_Low <- cbind(Description='GO_T24NSAZD', GO_T24NSAZD_Low)
GO_T24NPAZD_High <- cbind(Description='GO_T24NPAZD', GO_T24NPAZD_High)
GO_T24NPAZD_Low <- cbind(Description='GO_T24NPAZD', GO_T24NPAZD_Low)
GO_AzdHigh6hrs <- cbind(Description='GO_Azd6hrs', GO_AzdHigh6hrs)
GO_AzdLow6hrs <- cbind(Description='GO_Azd6hrs', GO_AzdLow6hrs)
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
GO_Azd_NSHigh6hrs <- cbind(Description='GO_Azd_NS6hrs', GO_Azd_NSHigh6hrs)
GO_Azd_NSLow6hrs <- cbind(Description='GO_Azd_NS6hrs', GO_Azd_NSLow6hrs)
GO_Azd_NPHigh6hrs <- cbind(Description='GO_Azd_NP6hrs', GO_Azd_NPHigh6hrs)
GO_Azd_NPLow6hrs <- cbind(Description='GO_Azd_NP6hrs', GO_Azd_NPLow6hrs)
GO_Azd_NSHigh24hrs <- cbind(Description='GO_Azd_NS24hrs', GO_Azd_NSHigh24hrs)
GO_Azd_NSLow24hrs <- cbind(Description='GO_Azd_NS24hrs', GO_Azd_NSLow24hrs)
GO_Azd_NPHigh24hrs <- cbind(Description='GO_Azd_NP24hrs', GO_Azd_NPHigh24hrs)
GO_Azd_NPLow24hrs <- cbind(Description='GO_Azd_NP24hrs', GO_Azd_NPLow24hrs)


#' # Combine all the tables
#' ## Analyzed in contrast to T0
#' ### Prepare the tables
combined_up_T0 <- rbind.fill(GO_T6NPSDMSO_High, GO_T6NSDMSO_High, GO_T6NPDMSO_High, GO_T6NPSAZD_High, GO_T6NSAZD_High, GO_T6NPAZD_High, GO_T24NPSDMSO_High, GO_T24NSDMSO_High, GO_T24NPDMSO_High, GO_T24NPSAZD_High, GO_T24NSAZD_High, GO_T24NPAZD_High)
combined_down_T0 <- rbind.fill(GO_T6NPSDMSO_Low, GO_T6NSDMSO_Low, GO_T6NPDMSO_Low, GO_T6NPSAZD_Low, GO_T6NSAZD_Low, GO_T6NPAZD_Low, GO_T24NPSDMSO_Low, GO_T24NSDMSO_Low, GO_T24NPDMSO_Low, GO_T24NPSAZD_Low, GO_T24NSAZD_Low, GO_T24NPAZD_Low)
#' ### Keep only the biological processes
combined_up_T0 <- combined_up_T0[combined_up_T0$namespace == "BP",]
combined_down_T0 <- combined_down_T0[combined_down_T0$namespace == "BP",]

#' ## Analyzed in contrast to NPS-DMSO-6hrs
#' ### Prepare the tables
combined_up_T0 <- rbind.fill(GO_AzdHigh6hrs, GO_AzdNPHigh6hrs, GO_AzdNSHigh6hrs, GO_PiHigh6hrs, GO_SucHigh6hrs)
combined_down_T0 <- rbind.fill(GO_AzdLow6hrs,GO_SucLow6hrs,GO_PiLow6hrs,GO_AzdNPLow6hrs,GO_AzdNSLow6hrs)
#' ### Keep only the biological processes
combined_up_T0 <- combined_up_T0[combined_up_T0$namespace == "BP",]
combined_down_T0 <- combined_down_T0[combined_down_T0$namespace == "BP",]

#' ## Analyzed in contrast to NPS-DMSO-24hrs
#' ### Prepare the tables
combined_up_T0 <- rbind.fill(GO_AzdHigh24hrs, GO_AzdNPHigh24hrs, GO_AzdNSHigh24hrs, GO_PiHigh24hrs, GO_SucHigh24hrs)
combined_down_T0 <- rbind.fill(GO_AzdLow24hrs,GO_SucLow24hrs,GO_PiLow24hrs,GO_AzdNPLow24hrs,GO_AzdNSLow24hrs)
#' ### Keep only the biological processes
combined_up_T0 <- combined_up_T0[combined_up_T0$namespace == "BP",]
combined_down_T0 <- combined_down_T0[combined_down_T0$namespace == "BP",]


#' ## Combining 6hrs and 24hrs
#' ### Prepare the tables
combined_up_T0 <- rbind.fill(GO_SucHigh6hrs,
                             GO_SucHigh24hrs,
                             GO_PiHigh6hrs,
                             GO_PiHigh24hrs,
                             GO_AzdHigh6hrs, 
                             GO_AzdHigh24hrs,
                             GO_AzdNPHigh6hrs, 
                             GO_AzdNPHigh24hrs, 
                             GO_AzdNSHigh6hrs, 
                             GO_AzdNSHigh24hrs)
combined_down_T0 <- rbind.fill(GO_SucLow6hrs,
                               GO_SucLow24hrs,
                               GO_PiLow6hrs,
                               GO_PiLow24hrs,
                               GO_AzdLow6hrs, 
                               GO_AzdLow24hrs,
                               GO_AzdNPLow6hrs, 
                               GO_AzdNPLow24hrs, 
                               GO_AzdNSLow6hrs, 
                               GO_AzdNSLow24hrs)
#' ### Keep only the biological processes
combined_up_T0 <- combined_up_T0[combined_up_T0$namespace == "BP",]
combined_down_T0 <- combined_down_T0[combined_down_T0$namespace == "BP",]

#' # Make a first plot with all the GO terms to see if it works
ggplot(combined_up_T0,aes(x=Description,y=id))+ geom_point(aes(size=nt, color=padj))+ 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("GO enrichments Up-regulated genes \n")+scale_shape_discrete(solid=F)+
    scale_colour_gradientn(colours=c("royalblue4", "lightskyblue1")) 

#' # Combining Up and Down-regulated genes together 
combined_up_T0[,"Direction"]  <- "Induced"
combined_down_T0[,"Direction"]  <- "Repressed"

combined_all_up_and_down_T0 <- rbind(combined_up_T0, combined_down_T0)


#' * Export the obtained data.frame containing also the non-significant GOs
write.table(combined_all_up_and_down_T0,file="combined_all_up_and_down_with_also_NON_sgnificant_GOs.txt", quote = F ,sep = "\t", row.names=FALSE)


#' # Removing the root terms
exclude <- c("GO:0008150", "GO:0008283", "GO:0032502", "GO:0044848", "GO:0032501", "GO:0098743", "GO:0008152", "GO:0051704", "GO:0007610", "GO:0000003", "GO:0009987",
             'GO:0006791', 'GO:0001906', 'GO:0002376', 'GO:0006794', 'GO:0009758', 'GO:0015976', 'GO:0019740', 'GO:0022414', 
             'GO:0022610', 'GO:0023052', 'GO:0040007', 'GO:0040011', 'GO:0043473', 'GO:0048511', 'GO:0048518', 'GO:0048519', 'GO:0050789', 'GO:0050896', 'GO:0051179', 'GO:0065007', 'GO:0071840', 'GO:0098754', 'GO:0110148')
combined_all_up_and_down_T0 <- combined_all_up_and_down_T0[! combined_all_up_and_down_T0$id %in% exclude,]

# ' # Change the p-values to change all the ones that are non-significant
combined_all_up_and_down_T0$padj[combined_all_up_and_down_T0$padj > 0.05] <- NA #Change this value as wished (normally 0.05)

#' # Make a list of GO of interest
myList <- GO_T6NPSDMSO_High[GO_T6NPSDMSO_High$padj <= 0.00001,3]
combined_all_up_and_down_T0 <- combined_all_up_and_down_T0[order(combined_all_up_and_down_T0$padj),]
myList <- combined_all_up_and_down_T0[1:80,3]
myList <- c(as.character(GO_T6NPSAZD_High[1:8,3]),as.character(GO_T6NPSAZD_Low[1:8,3]),as.character(GO_T6NPAZD_High[1:8,3]),as.character(GO_T6NPAZD_Low[1:8,3]),as.character(GO_T6NSAZD_High[1:8,3]),as.character(GO_T6NSAZD_Low[1:8,3]),
            as.character(GO_T6NPSDMSO_High[1:8,3]),as.character(GO_T6NPSDMSO_Low[1:8,3]),as.character(GO_T6NPDMSO_High[1:8,3]),as.character(GO_T6NPDMSO_Low[1:8,3]),as.character(GO_T6NSDMSO_High[1:8,3]),as.character(GO_T6NSDMSO_Low[1:8,3]),
            as.character(GO_T24NPSAZD_High[1:8,3]),as.character(GO_T24NPSAZD_Low[1:8,3]),as.character(GO_T24NPAZD_High[1:8,3]),as.character(GO_T24NPAZD_Low[1:8,3]),as.character(GO_T24NSAZD_High[1:8,3]),as.character(GO_T24NSAZD_Low[1:8,3]),
            as.character(GO_T24NPSDMSO_High[1:8,3]),as.character(GO_T24NPSDMSO_Low[1:8,3]),as.character(GO_T24NPDMSO_High[1:8,3]),as.character(GO_T24NPDMSO_Low[1:8,3]),as.character(GO_T24NSDMSO_High[1:8,3]),as.character(GO_T24NSDMSO_Low[1:8,3]))
myList <- GO_AzdLow6hrs[1:30,3]
myList <- c("GO:0022402","GO:0022403","GO:0016070","GO:0006260","GO:0051726","GO:0042254","GO:0006412")
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

myList <- list(C_Cycle,Pyrimidine,Replication)

#' # Create a sublist of GO with only the GOs of interest
go <- combined_all_up_and_down_T0[combined_all_up_and_down_T0$id %in% myList,]
go <- lapply(myList,function(x)
{
    combined_all_up_and_down_T0[combined_all_up_and_down_T0$id %in% x,]
})



#' # Make the bubble plot
ggplot(go,aes(x=Description,y=id))+ geom_point(aes(size=npat, color=padj))+ 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("GO enrichments Up-regulated genes \n")+scale_shape_discrete(solid=F)+
    scale_colour_gradientn(colours=c("royalblue4", "lightskyblue1")) 



#' # Make the bubble plot with only the selected GOs and in gray the ones with a non-significant p-value
ggplot(go,aes(x=Description,y=id))+ geom_point(aes(size=npat, color=padj))+ 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    ggtitle("GO enrichments Up and Down-regulated genes\n")+scale_shape_discrete(solid=F)+
    scale_colour_gradientn(colours=c("red2", "orange2"))

#' # Separate induced and repressed genes
ggplot(go,aes(x=Description,y=id))+ geom_point(aes(size=npat, color=padj))+ 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + facet_wrap(~Direction, nrow=1) +
    ggtitle("GO enrichments Up and Down-regulated genes\n")+scale_shape_discrete(solid=F)+
    scale_colour_gradientn(colours=c("red2", "orange2"))


#' # With a grid in the background
ggplot(go,aes(x=Description,y=id))+ geom_point(aes(size=npat, color=-log10(padj)))+ 
    facet_wrap(~Direction, nrow=1) +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("DNA replication (npat)\n")+scale_shape_discrete(solid=F)+
    scale_colour_gradientn(colours=c("red2", "purple2", "blue2")) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 



lapply(go,function(x)
{
    ggplot(x,aes(x=Description,y=id))+ geom_point(aes(size=npat, color=-log10(padj)))+ 
        facet_wrap(~Direction, nrow=1) +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        ggtitle("DNA replication (npat)\n")+scale_shape_discrete(solid=F)+
        scale_colour_gradientn(colours=c("red2", "purple2", "blue2")) +
        theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
})





#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```

