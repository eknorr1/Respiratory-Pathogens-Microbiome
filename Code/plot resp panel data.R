# Canine nosomes: Visualize respiratory panel data by host factors
#
# Feb 8, 2019

rm(list=ls())
#graphics.off()


# Script parameters -------------------------------------------------------

DoRandomization <- T   #shuffle the inverse CT scores for each pathogen?
DataDirectory <- "~/Desktop/Pilot Data/Raw Data"




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

CTcolumns <- 2:7
panel[,CTcolumns] <- 1/panel[,CTcolumns]
panel[,CTcolumns][is.na(panel[,CTcolumns])] <- 0



# Possibly shuffle  -------------------------------------------------------

randomize <- function(x){
  
  x <- as.matrix(x)
  for(i in 1:dim(x)[2]){
    ord <- order(runif(dim(x)[1]))
    x[,i] <- x[ord,i]
  }
  
  return(x)
  
}

if(DoRandomization){
  x <- panel[,CTcolumns]
  x <- randomize(x)
  panel[,CTcolumns] <- x
}


# Calculate proportions and load ------------------------------------------

x <- panel[,CTcolumns]
panel$Load <- rowSums(x)
panel$Richness <- rowSums(x>0)

p <- x
for(i in 1:nrow(panel)){
  p[i,] <- p[i,]/panel$Load[i]
}

n <- names(p)
n <- paste('p', n, sep='')
names(p) <- n

panel$Shannon <- -rowSums(p*log(p),na.rm=T)

panel <- cbind(panel,p)

nrep <- 1000
Richness0 <- matrix(NA,nrow(x),nrep)
Load0 <- matrix(NA,nrow(x),nrep)
Shannon0 <- matrix(NA,nrow(x),nrep)

for(i in 1:nrep){
  
  y <- randomize(x)
  Load0[,i] <- rowSums(y)
  Richness0[,i] <- rowSums(y>0)
  
  p <- y
  for(j in 1:nrow(p)){
    p[j,] <- p[j,]/Load0[j,i]
  }
  
  Shannon0[,i] <- -rowSums(p*log(p),na.rm=T)
  
}



rm(n,x,p,i,y)




# Shannon versus total load -----------------------------------------------

quartz()
plot(panel$Shannon, panel$Load)









# Works in progress... ----------------------------------------------------


if(0){


# Barplot -----------------------------------------------------------------

#quartz(h=4,w=10)
graphics.off()
symp <- subset(panel, panel$Symptomatic)
asymp <- subset(panel, !panel$Symptomatic)

x <- panel[,CTcolumns]
x <- x[order(panel$Load),]
x <- x[,c('MCYN', 'PINF', 'BCOR','BORD',"CADEN",'PNVPCR')]

pal <- 2:(ncol(x)+1)

barplot(t(as.matrix(x)),col = pal, border = pal, 
        names.arg = rep(NA,nrow(x)),
        ylab = 'Load')

legend('topleft',
       col = rev(pal),
       pt.bg = rev(pal),
       pch = rep(22,length(pal)),
       bty = 'n',
       legend = rev(names(x)))

#-----------------------------------
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



DogID<- c('195349','195356','195353','195319','195352', '195635')

# Barplot -----------------------------------------------------------------

panelTiny <- panel[which(panel$PetID==195349),]

#quartz(h=4,w=10)
graphics.off()
symp <- subset(panel, panel$Symptomatic)
asymp <- subset(panel, !panel$Symptomatic)

x <- panelTiny[,CTcolumns]
x <- x[order(panelTiny$CollectionDate),]
x <- x[,c('MCYN','BCOR','BORD',"CADEN",'PNVPCR')]

pal <- 2:(ncol(x)+1)

barplot(t(as.matrix(x)),col = pal, border = pal, 
        names.arg = rep(NA,nrow(x)),
        ylab = 'Load')

legend('topleft',
       col = rev(pal),
       pt.bg = rev(pal),
       pch = rep(22,length(pal)),
       bty = 'n',
       legend = rev(names(x)))

#######################################################
panelMichelle <- panel[which(panel$PetID==195356),]

x <- panelMichelle[,CTcolumns]
x <- x[order(panelMichelle$CollectionDate),]
x <- x[,c('MCYN','BCOR','BORD',"CADEN",'PNVPCR')]

pal <- 2:(ncol(x)+1)

barplot(t(as.matrix(x)),col = pal, border = pal, 
        names.arg = rep(NA,nrow(x)),
        ylab = 'Load')

legend('topleft',
       col = rev(pal),
       pt.bg = rev(pal),
       pch = rep(22,length(pal)),
       bty = 'n',
       legend = rev(names(x)))

#######################################################
panelOsborne <- panel[which(panel$PetID==195353),]

x <- panelOsborne[,CTcolumns]
x <- x[order(panelOsborne$CollectionDate),]
x <- x[,c('MCYN','BCOR','BORD',"CADEN",'PNVPCR')]

pal <- 2:(ncol(x)+1)

barplot(t(as.matrix(x)),col = pal, border = pal, 
        names.arg = rep(NA,nrow(x)),
        ylab = 'Load')

legend('topleft',
       col = rev(pal),
       pt.bg = rev(pal),
       pch = rep(22,length(pal)),
       bty = 'n',
       legend = rev(names(x)))

#######################################################
panelSqueak <- panel[which(panel$PetID==195319),]

x <- panelSqueak[,CTcolumns]
x <- x[order(panelSqueak$CollectionDate),]
x <- x[,c('MCYN','BCOR','BORD',"CADEN",'PNVPCR')]

pal <- 2:(ncol(x)+1)

barplot(t(as.matrix(x)),col = pal, border = pal, 
        names.arg = rep(NA,nrow(x)),
        ylab = 'Load')

legend('topleft',
       col = rev(pal),
       pt.bg = rev(pal),
       pch = rep(22,length(pal)),
       bty = 'n',
       legend = rev(names(x)))

#######################################################
panelBismark <- panel[which(panel$PetID==195352),]

x <- panelBismark[,CTcolumns]
x <- x[order(panelBismark$CollectionDate),]
x <- x[,c('MCYN','BCOR','BORD',"CADEN",'PNVPCR')]

pal <- 2:(ncol(x)+1)

barplot(t(as.matrix(x)),col = pal, border = pal, 
        names.arg = rep(NA,nrow(x)),
        ylab = 'Load')

legend('topleft',
       col = rev(pal),
       pt.bg = rev(pal),
       pch = rep(22,length(pal)),
       bty = 'n',
       legend = rev(names(x)))


#######################################################
panelSally <- panel[which(panel$PetID==195635),]

x <- panelSally[,CTcolumns]
x <- x[order(panelSally$CollectionDate),]
x <- x[,c('MCYN','BCOR','BORD',"CADEN",'PNVPCR')]

pal <- 2:(ncol(x)+1)

barplot(t(as.matrix(x)),col = pal, border = pal, 
        names.arg = rep(NA,nrow(x)),
        ylab = 'Load')

legend('topleft',
       col = rev(pal),
       pt.bg = rev(pal),
       pch = rep(22,length(pal)),
       bty = 'n',
       legend = rev(names(x)))







