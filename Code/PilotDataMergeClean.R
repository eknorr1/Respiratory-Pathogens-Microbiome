
# Respiratory Pathogen Community Analysis
# Initial data merging and cleaning of unnecessary columns 
# Liz Knorr
# created 4/24/2019

######## clear workspace #######

rm(list = ls())
graphics.off()

# Load Libraries -------------------------------------------------------

# None

# Script Parameters ----------------------------------------------------
#set to where data files are located

DataDirectory <- "~/Desktop/Pilot Data/Raw Data"
CTtarget <- c("CADEN", "BCOR", "BORD", "MCYN", "PINF", "PNVPCR")


# Load and merge host data and pathogen panel data ----------------------

setwd(DataDirectory)

host <- read.csv(file = "SMB_Pilot_Data_1.csv", stringsAsFactors = F)
keep <- c('Well','PetID','SampleType','PetName','Origin','CollectionDate','Symptomatic', 'VialID')
host <- host[,keep]


panel <- read.csv(file = "paneldata_1.csv", stringsAsFactors = F)
discard <- which(names(panel)%in%c("X","MS2","SZ","Type"))  #discard useless columns
panel <- panel[,-discard]

panel <- merge(panel, host, by = "Well")
panel <- subset(panel, SampleType == 'nasal')   #subset only nasal samples

panel$Symptomatic <- panel$Symptomatic=="Y"
panel$Origin <- as.factor(panel$Origin)
panel$CollectionDate <- strptime(panel$CollectionDate, "%m/%d/%y")

rm(host, keep, discard)



# Manually fix some errors ------------------------------------------------

# Date coded wrong for November samples
i <- which(as.character(panel$CollectionDate)=="2016-01-11")
panel$CollectionDate[i] <- strptime("2016-11-01", "%Y-%m-%d")


# In origin, shall we say OHS? means OHS 
panel$Origin[panel$Origin=='OHS?'] <- "OHS"
panel$Origin <- factor(panel$Origin)



# Add some variables ------------------------------------------------------

CTcolumn <- which(names(panel) %in% CTtarget)
x <- panel[,CTcolumn]
panel$Richness <- rowSums(!is.na(x))


# Update some variables by looking at individuals -------------------------

panel$SamplingSequence <- NA
panel$EverSymptomatic <- FALSE
panel$BecameSymptomatic <- FALSE

uniqueID <- unique(panel$PetID)
nID <- length(uniqueID)

for(i in 1:nID){
  
  j <- which(panel$PetID==uniqueID[i])
  
  ss <- panel[j,]
  ss <- ss[order(ss$CollectionDate),]       #ordering by time...
  
  ss$Origin <- ss$Origin[1]                 #set Origin to what it was on their first sample
  
  ss$SamplingSequence <- 1:length(j)        
  ss$EverSymptomatic <- any(ss$Symptomatic)
  ss$BecameSymptomatic <- any(ss$Symptomatic) & !ss$Symptomatic[1]
  
  panel[j,] <- ss                           #...they go back in a different order than they went in
  
}

rm(ss)

# Subest to get just nasal samples 

panel <- subset(panel, !is.na(panel$nasal))  #doesn't work 
panel <- subset(panel, !is.na(panel$PetID))   #remove controls

panel <- panel[,-which(names(panel)=='DIS')]  #remove distemper
