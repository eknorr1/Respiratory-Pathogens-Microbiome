# Canine nosomes: Merge panel and 16S data
#
# April 22, 2019
#Edited June 11 by Liz Knorr

rm(list = ls())
graphics.off()


library(metacoder)
library(phyloseq)
library(RColorBrewer)



# Paths -------------------------------------------------------------------

pathToAssembledData <- "~/Users/elizabethknorr/Respiratory-Pathogens_Microbiome/Microbiome Code"




# Load data ---------------------------------------------------------------

# sxs (for '16S') is the microbiome data
# sam (for 'sample') is the sample data 

setwd(pathToAssembledData)
load("16S.Rdata")
load("sample.Rdata")





# Augment sample data in sxs ----------------------------------------------

samsxs <- data.frame(sample_data(sxs))
samsxs_aug <- base::merge(samsxs, sam, all.x = T, by = 'SampleIndex') 
rownames(samsxs_aug) <- samsxs_aug$SampleID                                   #preserve rownames as SampleID because phyloseq wants that
samsxs_aug$SampleType[is.na(samsxs_aug$SampleType)] <- "control"              #default covers samples introduced during library prep
samsxs_aug$Symptomatic[is.na(samsxs_aug$Symptomatic)] <- FALSE                #default covers samples introduced during library prep
sample_data(sxs) <- samsxs_aug





# Save --------------------------------------------------------------------

save(sxs, file = "16S_combined.Rdata")

 


