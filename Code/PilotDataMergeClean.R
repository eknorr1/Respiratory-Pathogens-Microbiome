
# Respiratory Pathogen Community Analysis
# Initial data merging and cleaning of unnecessary columns 
# Liz Knorr
# created 4/24/2019

######## clear workspace #######

rm(list = ls())
graphics.off()

# Load Libraries -------------------------------------------------------

library(vegan)

# Script Parameters ----------------------------------------------------

#set to where data files are located

DataDirectory <- "~/Desktop/Pilot Data/Raw Data"
CTtarget <- c("CADEN", "BCOR", "BORD", "MCYN", "PINF", "PNVPCR")



# Load and merge datasets -------------------------------------------------

setwd(DataDirectory)

host <- read.csv(file = "SMB_Pilot_Data_1.csv", stringsAsFactors = F)
keep <- c('Well','PetID','PetName','Origin','CollectionDate','Symptomatic','TotalDNA','TotalRNA')
host <- host[,keep]

suppressWarnings(host$TotalDNA <- as.numeric(host$TotalDNA))
suppressWarnings(host$TotalRNA <- as.numeric(host$TotalRNA))

host$Symptomatic <- host$Symptomatic=="Y"

panel <- read.csv(file = "paneldata_1.csv", stringsAsFactors = F)
discard <- which(names(panel)%in%c("X","MS2","SZ","Type"))
panel <- panel[,-discard]

panel <- merge(panel, host, by = "Well")

rm(host, keep, discard)


# Subset ------------------------------------------------------------------

panel <- subset(panel, !is.na(panel$PetID))   #remove controls

panel <- panel[,-which(names(panel)=='DIS')]  #remove distemper



# Invert CT score ---------------------------------------------------------

# !!! Check colnames that 2 - 7 are "CADEN", "BCOR", "BORD", "MCYN", "PINF", "PNVPCR"
print(colnames(panel))


CTcolumns <- 2:7

panel[,CTcolumns] <- 1/panel[,CTcolumns]
panel[,CTcolumns][is.na(panel[,CTcolumns])] <- 0

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

# Subset ------------------------------------------------------------------

panel <- subset(panel, !is.na(panel$PetID))   #remove controls

panel <- panel[,-which(names(panel)=='DIS')]  #remove distemper



# Invert CT score ---------------------------------------------------------

CTcolumns <- 2:7
panel[,CTcolumns] <- 1/panel[,CTcolumns]
panel[,CTcolumns][is.na(panel[,CTcolumns])] <- 0


# Manually fix some errors ------------------------------------------------

# Date coded wrong for November samples
i <- which(as.character(panel$CollectionDate)=="2016-01-11")
panel$CollectionDate[i] <- strptime("2016-11-01", "%Y-%m-%d")


# In origin, shall we say OHS? means OHS 
panel$Origin[panel$Origin=='OHS?'] <- "OHS"
panel$Origin <- factor(panel$Origin)



# Add some variables ------------------------------------------------------

# Total load 

panel$Total_Load <- panel$CADEN + panel$BCOR + panel$BORD + panel$MCYN 
  + panel$PINF + panel$PNVPCR 

#Number of infections

panel$Coinfection<- apply(panel[,2:7]>0,1,sum)


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

# State: Cali vs Oregon --------------------------------------------------

panel$State <- ifelse(panel$Origin=="Madera", "California",
                      ifelse(panel$Origin == "Hayward", "California",
                             ifelse(panel$Origin == "Merced", "California",
                                    ifelse(panel$Origin == "Sacramento", "California",
                                           ifelse(panel$Origin == "Fresno", "California",
                                                  ifelse(panel$Origin=="Oakland", "California",
                                                         ifelse(panel$Origin=="OHS", "Oregon",
                                                                NA )))))))



# Shannon diversity ------------------------------------------------------


panel$Shannon <- diversity(panel[,2:7], index = "shannon")

# Evenness ---------------------------------------------------------------

panel$Evenness <- (panel$Shannon/ log(panel$Coinfection))


even <- panel$Evenness=="0"
panel$Evenness[even] <- "NaN"


