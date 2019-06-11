# Canine nosomes: Assemble respiratory panel and host data
#
# Feb 8, 2019

rm(list=ls())
graphics.off()

# Script parameters -------------------------------------------------------

InputDirectory <- "~/Dropbox/Research/Active/Canine nosomes/Data/Pilot data from Liz Jan2019/"
OutputDirectory <- "~/Dropbox/Research/Active/Canine nosomes/Data/Assembled/"




# Import data -------------------------------------------------------------

# Note: 'Well' is a sequential sample index used by us and Cornell, refering to the order in which we processed the samples
# Below we will rename 'Well' as 'SampleIndex', which is more informative
# Note this is not the same as sampleid, which also includes 'samples' introduced at cgrb ie blanks

setwd(InputDirectory)

host <- read.csv(file = "SMB_Pilot_Data_1.csv", stringsAsFactors = F)
keep <- c('Well','PetID','PetName','Origin','CollectionDate','Symptomatic','TotalDNA','TotalRNA','VialID','SampleType')
host <- host[,keep]

panel <- read.csv(file = "paneldata_1.csv", stringsAsFactors = F)
panel <- panel[,-which(names(panel)%in%c("X","MS2","SZ","Type"))]
panel <- subset(panel, !is.na(panel$Well)) #NA well numbers are internal controls used by Cornell



# Clean up host data ------------------------------------------------------

host$Well[host$Well=='145b'] <- 145    #fix a typo
host$Well <- as.numeric(host$Well)

host <- subset(host, !is.na(host$Well))

suppressWarnings(host$TotalDNA <- as.numeric(host$TotalDNA))
suppressWarnings(host$TotalRNA <- as.numeric(host$TotalRNA))

host$SampleType[host$SampleType == ""] <- "control"





# Merge host and panel data -----------------------------------------------
# Here we rename 'Well' as 'SampleIndex', which is more informative

sam <- merge(panel, host, by = "Well")
n <- names(sam)
n[n=='Well'] <- 'SampleIndex'
names(sam) <- n


# Adjust the classes of some variables ------------------------------------

sam$Symptomatic <- sam$Symptomatic=="Y"
sam$Origin <- as.factor(sam$Origin)
sam$CollectionDate <- strptime(sam$CollectionDate, "%m/%d/%y")






# Manually fix some errors ------------------------------------------------

# Date coded wrong for November samples
i <- which(as.character(sam$CollectionDate)=="2016-01-11")
sam$CollectionDate[i] <- strptime("2016-11-01", "%Y-%m-%d")


# In origin, shall we say OHS? means OHS 
sam$Origin[sam$Origin=='OHS?'] <- "OHS"
sam$Origin <- factor(sam$Origin)





# Update some variables by looking at individuals -------------------------

sam$SamplingSequence <- NA
sam$EverSymptomatic <- FALSE
sam$BecameSymptomatic <- FALSE

uniqueID <- unique(sam$PetID)
uniqueID <- uniqueID[!is.na(uniqueID)]
nID <- length(uniqueID)

for(i in 1:nID){
  
  j <- which(sam$PetID==uniqueID[i])
  
  ss <- sam[j,]
  ss <- ss[order(ss$CollectionDate),]       #ordering by time...
  
  ss$Origin <- ss$Origin[1]                 #set Origin to what it was on their first sample
  
  ss$SamplingSequence <- 1:length(j)        
  ss$EverSymptomatic <- any(ss$Symptomatic)
  ss$BecameSymptomatic <- any(ss$Symptomatic) & !ss$Symptomatic[1]
  
  sam[j,] <- ss                           #...they go back in a different order than they went in
  
}

rm(ss)




# Save --------------------------------------------------------------------

setwd(OutputDirectory)
save(sam, file = 'sample.Rdata')
