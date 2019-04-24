# Canine nosomes: Assemble and plot panel data
# creates figures for temporal analysis of pathogens present in symptmatic and asymptomatic dogs
# Feb 8, 2019

rm(list=ls())
graphics.off()

# Script parameters -------------------------------------------------------

DataDirectory <- "~/Dropbox/Research/Active/Canine nosomes/Data/Pilot data from Liz Jan2019/"
CTtarget <- c("CADEN", "BCOR", "BORD", "MCYN", "PINF", "PNVPCR")



# Load and merge datasets -------------------------------------------------

setwd(DataDirectory)

host <- read.csv(file = "SMB_Pilot_Data_1.csv", stringsAsFactors = F)
keep <- c('Well','PetID','PetName','Origin','CollectionDate','Symptomatic','TotalDNA','TotalRNA')
host <- host[,keep]

suppressWarnings(host$TotalDNA <- as.numeric(host$TotalDNA))
suppressWarnings(host$TotalRNA <- as.numeric(host$TotalRNA))

panel <- read.csv(file = "paneldata_1.csv", stringsAsFactors = F)
discard <- which(names(panel)%in%c("X","MS2","SZ","Type"))  #discard useless columns
panel <- panel[,-discard]

panel <- merge(panel, host, by = "Well")
panel <- subset(panel, !is.na(panel$PetID))                  #remove controls

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







# Exploratory plots -------------------------------------------------------



# CT score by symptomatic and sampling seqeunce
quartz(h=6,w=6, title = "Shedding by sampling sequence and disease status")
par(mfrow=c(3,2))
par(mar=c(3,4,3,2))

CTcolumn <- which(names(panel) %in% CTtarget)
symp <- rep("N", nrow(panel))
symp[panel$Symptomatic] <- "Y"
sampxsymp <- as.factor(paste(panel$SamplingSequence,symp, sep = ""))
nl <- nlevels(sampxsymp)

cols <- c('blue','red')
at <- c(1,1.75,3,3.75,5,5.75,7,7.75,9:11)


for(i in 1:length(CTcolumn)){
  x <- panel[,CTcolumn[i]]
  y <- 1/x
  y[!is.finite(y)] <- 0 
  boxplot(y~sampxsymp, at = at,
          varwidth = T, col = cols, ylim = 1.5*range(y), lty=1, boxcol = cols,
          cex.axis = 0.8,
          xlim = c(0.5,8), boxwex = 1, main = names(panel)[CTcolumn[i]],ylab = "1/CT")
  
  n <- as.numeric(table(sampxsymp))
  text(at,rep(1.35*max(y),nl),paste("n=",n,sep=""),cex=0.7)
  
  
  if(i == 2){
    
    legend('topright', pch = rep(22,2), col = cols, pt.bg = cols,  pt.cex = 2, 
           bty='n', inset = c(0,.15),
           legend = c('No signs of CIRD', 'Signs of CIRD'))
    
  }
  
  
}

rm(CTcolumn,i,j,n,nID,nl,sampxsymp,uniqueID,x,y, symp)




# CT score by shelter and sampling sequence, for selected locations and times
quartz(h=6,w=6, title = "Shedding by origin and sampling sequence for selected locations")
par(mfrow=c(3,2))
par(mar=c(3,4,3,2))

cols <- c('darkblue', 'lightblue')

ss <- subset(panel, panel$Origin %in% c("Fresno", "Madera", "Merced", "Oakland"))
ss <- subset(ss, SamplingSequence<3)

ss$Origin <- factor(ss$Origin)

locxsamp <- paste(ss$Origin, ss$SamplingSequence, sep="")

CTcolumn <- which(names(ss) %in% CTtarget)
for(i in 1:length(CTcolumn)){
  
  x <- ss[,CTcolumn[i]]
  y <- 1/x
  y[!is.finite(y)] <- 0 
  
  at <- 1:8 - rep(c(0,0.2),4)
  
  boxplot(y~locxsamp, at = at,
          ylim = 1.5*range(y), lty=1, varwidth=T,
          main = names(ss)[CTcolumn[i]],ylab = "1/CT",
          col = cols,
          boxcol = cols,
          xaxt='n')
  
  axis(1, at, labels = NA)
  axis(1, seq(1,8,2)+0.4, levels(ss$Origin), tick = F)
  
  n <- as.numeric(table(locxsamp))
  text(at,rep(1.35*max(y),8),paste("n=",n,sep=""),cex=0.7)
  
  if(i == 2){
    
    legend('topright', pch = rep(22,2), col = cols, pt.bg = cols,  pt.cex = 2, 
           bty='n', inset = c(0,.15),
           legend = c('Sample 1', 'Sample 2'))
    
  }
  
  
}







# Questions

# STUDY DESIGN
# [ ] how many dogs?
# [ ] breakdown of where they came from
# [ ] when we sampled
# [ ] how many repeat samples?
# [ ] (infographic of what we did to the samples)


# DISEASE (less emphasis on disease per se because they can monitor this themselves)
# [ ] how many sampling events had associated signs of CIRD?
# [ ] what was a dog's chance of developing CIRD? (overall, and looking at transition from not having it)


# RESPIRATORY PANEL DATA
# [X] heatmap pathogen prevalence across samples
# [X] initial and secondary sample shedding by shelter of origin


# PANEL DATA AND DISEASE
# [X] shedding by sampling sequence and disease status 
# [ ] future symptoms based on past panels? eg myco now, versus future symptoms






