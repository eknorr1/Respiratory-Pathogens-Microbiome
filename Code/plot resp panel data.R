# Canine nosomes: Visualize respiratory panel data by host factors
#
# Feb 8, 2019

rm(list=ls())
#graphics.off()


# Script parameters -------------------------------------------------------

DoRandomization <- T   #shuffle the inverse CT scores for each pathogen?
DataDirectory <- "~/Dropbox/Research/Active/Canine nosomes/Data/Pilot data from Liz Jan2019/"




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

quartz(h=4,w=10)

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


# Load versus Shannon -----------------------------------------------------

quartz(h=6,w=8)
par(mfrow=c(2,3))
plot(jitter(panel$Richness), panel$MCYN)
plot(jitter(panel$Richness), panel$PNVPCR)
plot(jitter(panel$Richness), panel$PINF)

plot(panel$Shannon, panel$MCYN)
plot(panel$Shannon, panel$PNVPCR)
plot(panel$Shannon, panel$PINF)



}