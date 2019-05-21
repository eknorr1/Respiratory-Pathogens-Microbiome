
# Respiratory Pathogen Microbiome Pilot Analysis 
# Respiratory Pathogen Community Analysis
# Liz Knorr
# created 4/24/2019

######## clear workspace #######

rm(list = ls())
graphics.off()

# Load Libraries -------------------------------------------------------

library(vegan)
library(corrplot)
library(factoextra)
library(FactoMineR)

# Script Parameters ----------------------------------------------------

#set to where data files are located
DataDirectory <- "~/Desktop/Pilot Data/Raw Data"


# Load and merge datasets -------------------------------------------------

setwd(DataDirectory)

panel <- read.csv(file = "panel.csv", stringsAsFactors = F)


# Analysis Starts -------------------------------------------------


#levels(panel$Symptomatic)[1]<- "0"
#levels(panel$Symptomatic)[2]<- "1"


# Species richness plot ------------------------------------------------------

Species_Richness <- sort(panel$Coinfection,decreasing = FALSE)

plot(Species_Richness)

# Shannon diversity ---------------------------------------------------------
Shannon <- diversity(panel[,3:8], index = "shannon")  # note hard coding for the columns of inv Ct scores
hist(Shannon) 

Shan_increasing <- sort(panel$Shannon,decreasing = FALSE)
plot(Shan_increasing)


plot(panel$Shannon,panel$Total_Load, col=panel$Coinfection)

##########


#EVENNESS CALCULATION HERE
panel$Evenness <- (panel$Shannon/ log(panel$Coinfection))


even <- panel$Evenness=="0"
panel$Evenness[even] <- "NaN"
plot(panel$Evenness, panel$Total_Load, col = panel$Coinfection)

plot(panel$Evenness, panel$Total_Load, col = panel$Coinfection, 
     pch = 16, xlab = '', ylab = '',cex=3, cex.lab=3, cex.axis=3, cex.main=3, cex.sub=3)

legend(x="bottomleft", 
       legend = paste(c('2','3','4','5','6')),
       col = c('red','green','blue','cyan','magenta'),
       title = 'Infections',
       pch = 16,
       inset = 0.02,
       cex = 1.5,
       pt.cex=2)


plot(panel$Evenness, panel$Total_Load, col = panel$Symptomatic, 
     pch = 16, xlab = '', ylab = '',cex=3, cex.lab=3, cex.axis=3, cex.main=3, cex.sub=3)

load_myco <- lm(Total_Load ~MCYN, data = panel)
summary(load_myco)


Nasal_Pilot_Data$Index<- seq.int(nrow(Nasal_Pilot_Data))

Nasal_Pilot_Data<-Nasal_Pilot_Data[order(Nasal_Pilot_Data$MCYN_CT_invCT),]
Nasal_Pilot_Data$Index <- as.numeric(rownames(Nasal_Pilot_Data))
Nasal_Pilot_Data$Index <- seq.int(nrow(Nasal_Pilot_Data))


ggplot(data=Nasal_Pilot_Data) +
  geom_point(aes(x=Nasal_Pilot_Data$Index,y=Nasal_Pilot_Data$MCYN_CT_invCT, colour = IsSymptomatic), shape=19, size=5) +
  labs(x="Nasal Sample") +
  labs(y="Mycoplasma Shedding") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black", size = 0.5, linetype = "solid")) +
  theme(axis.text.x = element_text(face="plain", color="black", size=14, angle=0), axis.text.y = element_text(face="plain", color="black", size=14, angle=0)) +
  theme(axis.title = element_text(color="black", face="bold", size=18)) +
  theme(legend.position="bottom") +
  theme(legend.title = element_text(colour="black", size=14, face="bold")) +
  theme(legend.text = element_text(colour="black", size=14, face="plain")) +
  scale_x_continuous(limits=c(0, 105)) +
  scale_y_continuous(limits=c(0,0.5)) +
  scale_color_manual(values = c("darkgrey","red3")
  )


### Needs Work 

par(mfrow=c(1,2))
plot(Nasal_Pilot_Data$Evenness[Nasal_Pilot_Data$MCYN_CT_invCT>0],Nasal_Pilot_Data$MCYN_CT_invCT[Nasal_Pilot_Data$MCYN_CT_invCT>0], col = Nasal_Pilot_Data$IsSymptomatic, pch = 16, xlab = "Evenness",
     ylab = "Mycoplasma 1/CT",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(Nasal_Pilot_Data$Evenness[Nasal_Pilot_Data$PNVPCR_CT_invCT>0],Nasal_Pilot_Data$PNVPCR_CT_invCT[Nasal_Pilot_Data$PNVPCR_CT_invCT>0], col = Nasal_Pilot_Data$IsSymptomatic, pch = 16, xlab = "Evenness",
     ylab = "Pneumovirus 1/CT", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)

graphics.off()
plot(Nasal_Pilot_Data$Evenness, Nasal_Pilot_Data$Total_Load, pch = 16, xlab = "Evenness", ylab = "Total Load (1/CT)",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, col = Nasal_Pilot_Data$Coinfection)



graphics.off()


plot(Nasal_Pilot_Data$Evenness,Nasal_Pilot_Data$MCYN_CT_invCT, col = c("black","red")[Nasal_Pilot_Data$IsSymptomatic], 
     pch = 19, xlim = c(0.4,1), ylim = c(0,0.43), xlab = "", ylab = "", cex.lab=3, cex.axis=3, cex.main=3, cex.sub=3, cex = 4)


plot(Nasal_Pilot_Data$Evenness,Nasal_Pilot_Data$PNVPCR_CT_invCT, col = c("black","red")[Nasal_Pilot_Data$IsSymptomatic], 
     pch = 19, xlim = c(0.4,1), ylim = c(0,0.43), xlab = "", 
     ylab = "",cex.lab=3, cex.axis=3, cex.main=3, cex.sub=3, cex = 4, lwd = 3)
legend(x="topright", 
       legend = paste(c('Asymptomatic','Symptomatic')),
       col = c('black','red'),
       title = 'Health Status',
       pch = 16,
       inset = 0.02,
       cex = 1.75,
       pt.cex=2.25)

dog_uniq<- unique(Nasal_Pilot_Data$PetName)
length(dog_uniq)

infected<- length(which(Nasal_Pilot_Data$Total_Load>0))
print(infected)







#####    GLMs   ########

Nasal_Pilot_Data$Origin <- factor(as.character(Nasal_Pilot_Data$Origin))

Nasal_Pilot_Data$IsSymptomatic <- factor(Nasal_Pilot_Data$Symptomatic==1)

fit <- glm(IsSymptomatic ~ Total_Load, data = Nasal_Pilot_Data, family = binomial(link = "logit"))
summary(fit)

myco_total_load <- glm(IsSymptomatic ~ Total_Load + MCYN_CT_invCT, data = Nasal_Pilot_Data, family = binomial(link = "logit"))
summary(myco_total_load)


Sympt_fit2 <- glm(IsSymptomatic ~ CADEN_CT_invCT + MS2_CT_invCT + BCOR_CT_invCT + BORD_CT_invCT + 
                    DIS_CT_invCT + MCYN_CT_invCT + PINF_CT_invCT + PNVPCR_CT_invCT + SZ_CT_invCT, 
                  data = Nasal_Pilot_Data, family = binomial(link = "logit"))
summary(Sympt_fit2)

fit3 <- glm(IsSymptomatic ~ MCYN_CT_invCT, data = Nasal_Pilot_Data, family = binomial(link = "logit"))
summary(fit3)








#####     Boxplots     #####
graphics.off()
boxplot(Coinfection~Origin, data = Nasal_Pilot_Data, ylab = "Number of Coinfections", xlab = "Shelter of Origin")

boxplot(Coinfection~Symptomatic, data = Nasal_Pilot_Data, ylab = "", names=c("Asymptomatic","Symptomatic"),cex.lab=2.5,
        cex.axis=2.5, cex.main=2.5, cex.sub=2.5, col=c("blue","red"), lwd = 3)

boxplot(Coinfection~State, data = Nasal_Pilot_Data, ylab = "Number of Coinfections", xlab = "Shelter of Origin")

boxplot(Nasal_Pilot_Data$MCYN_CT_invCT[Nasal_Pilot_Data$MCYN_CT_invCT>0]~IsSymptomatic[Nasal_Pilot_Data$MCYN_CT_invCT>0],
        data = Nasal_Pilot_Data, ylab = "Mycoplasma Load", xlab = "Symptomatic Status", names=c("",""),cex.lab=2.5, cex.axis=2.5,
        cex.main=2.5, cex.sub=2.5, col=c("blue","red"), lwd = 2.7)

boxplot(Nasal_Pilot_Data$Coinfection[Nasal_Pilot_Data$Total_Load>0]~IsSymptomatic[Nasal_Pilot_Data$Total_Load>0],
        data = Nasal_Pilot_Data, ylab = "Number of Infections", xlab = "Symptomatic Status",names=c("",""),cex.lab=2.5,
        cex.axis=2.5, cex.main=2.5, cex.sub=2.5, col=c("blue","red"), lwd = 2.7)

boxplot(MCYN_CT_invCT~IsSymptomatic, data = Nasal_Pilot_Data, ylab = "Mycoplasma 1/CT", xlab = "Symptomatic Status")

boxplot(Nasal_Pilot_Data$MCYN_CT_invCT[Nasal_Pilot_Data$MCYN_CT_invCT>0]~IsSymptomatic[Nasal_Pilot_Data$MCYN_CT_invCT>0],
        data = Nasal_Pilot_Data, ylab = "Mycoplasma Load", xlab = "Symptomatic Status", names=c("",""),cex.lab=2.5, cex.axis=2.5,
        cex.main=2.5, cex.sub=2.5, col=c("blue","red"), lwd = 2.7)

boxplot(Nasal_Pilot_Data$Coinfection[Nasal_Pilot_Data$Total_Load>0]~IsSymptomatic[Nasal_Pilot_Data$Total_Load>0],
        data = Nasal_Pilot_Data, ylab = "Number of Infections", xlab = "Symptomatic Status",names=c("",""),cex.lab=2.5,
        cex.axis=2.5, cex.main=2.5, cex.sub=2.5, col=c("blue","red"), lwd = 2.7)

boxplot(MCYN_CT_invCT ~Symptomatic, data = Nasal_Pilot_Data, ylab = "", names=c("",""),cex.lab=2.5,
        cex.axis=2.5, cex.main=2.5, cex.sub=2.5, col=c("blue","red"), lwd = 2.7)

boxplot(Nasal_Pilot_Data$MCYN_CT_invCT[Nasal_Pilot_Data$MCYN_CT_invCT>0]~IsSymptomatic[Nasal_Pilot_Data$MCYN_CT_invCT>0],
        data = Nasal_Pilot_Data, ylab = "", xlab = "", names=c("",""),cex.lab=2.5, cex.axis=2.5,
        cex.main=2.5, cex.sub=2.5, col=c("royalblue3","red3"), lwd = 2.7)

boxplot(Nasal_Pilot_Data$Coinfection[Nasal_Pilot_Data$Total_Load>0]~IsSymptomatic[Nasal_Pilot_Data$Total_Load>0],
        data = Nasal_Pilot_Data, ylab = "", xlab = "",names=c("",""),cex.lab=2.5,
        cex.axis=2.5, cex.main=2.5, cex.sub=2.5, col=c("royalblue3","red3"), lwd = 2.7)

boxplot(Nasal_Pilot_Data$MCYN_CT_invCT[Nasal_Pilot_Data$Total_Load>0]~IsSymptomatic[Nasal_Pilot_Data$Total_Load>0],
        data = Nasal_Pilot_Data, ylab = "", xlab = "Symptomatic Status", names=c("",""),cex.lab=2.5, cex.axis=2.5,
        cex.main=2.5, cex.sub=2.5, col=c("royalblue3","red3"), lwd = 2.7)

par(mfrow=c(1,2))
par(mar=c(3,4,3,2), mai = c(1.25, 1.25, 1.25, 1.25))

boxplot(Nasal_Pilot_Data$Coinfection[Nasal_Pilot_Data$Total_Load>0]~IsSymptomatic[Nasal_Pilot_Data$Total_Load>0],
        data = Nasal_Pilot_Data, ylab = "", xlab = "",names=c("",""),cex.lab=2.5,
        cex.axis=2.5, cex.main=2.5, cex.sub=2.5, col=c("royalblue3","red3"), lwd = 2.7)
boxplot(Nasal_Pilot_Data$Total_Load[Nasal_Pilot_Data$Total_Load>0]~IsSymptomatic[Nasal_Pilot_Data$Total_Load>0],
        data = Nasal_Pilot_Data, ylab = "", xlab = "", names=c("",""), cex.lab=2.5, 
        cex.axis=2.5, cex.main=2.5, cex.sub=2.5, col=c("royalblue3","red3"), lwd = 2.7)


######         Histograms        ######


graphics.off()
par(mfrow=c(2,3), mar=c(7,6,6,5))
Mycoplasma<-Nasal_Pilot_Data$MCYN_CT_invCT
hist(Mycoplasma, xlab = "", xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, col='darkgrey')
Pneumovirus<-Nasal_Pilot_Data$PNVPCR_CT_invCT
hist(Pneumovirus, xlab = "1/CT (Shedding)", ylab = "", xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, col='darkgrey')
Betacoronavirus<-Nasal_Pilot_Data$BCOR_CT_invCT
hist(Betacoronavirus, xlab = "", ylab = "", xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, col='darkgrey')
Bordetella<-Nasal_Pilot_Data$BORD_CT_invCT
hist(Bordetella, xlab = "", xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, col='darkgrey')
Parainfluenza<-Nasal_Pilot_Data$PINF_CT_invCT
hist(Parainfluenza, xlab = "1/CT (Shedding)",ylab = "", xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, col='darkgrey')
Adenovirus<-Nasal_Pilot_Data$CADEN_CT_invCT
hist(Adenovirus, xlab = "",ylab="", xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, col='darkgrey')






#........Histograms without zeros and colored by low, moderate, and high positives

graphics.off()
par(mfrow=c(2,3))

Mycoplasma<-Nasal_Pilot_Data$MCYN_CT_invCT
myco_hist<- hist(Mycoplasma[!Mycoplasma<0.060887], xlab = "", xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=3, cex.sub=2, main = "Mycoplasma", col='darkgrey')
myco_col<- cut(myco_hist$breaks, c(-Inf, 0.12177302, Inf))

Pneumovirus<-Nasal_Pilot_Data$PNVPCR_CT_invCT
pneumo_hist<- hist(Pneumovirus[!Pneumovirus<0.05445], xlab = "1/CT (Shedding)", ylab = "", xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=3, cex.sub=2, col='darkgrey', main = "Pneumovirus")
pneumo_col<- cut(pneumo_hist$breaks, c(-Inf, 0.10890043, Inf))

Betacoronavirus<-Nasal_Pilot_Data$BCOR_CT_invCT
beta_hist<-hist(Betacoronavirus[!Betacoronavirus<0.04918], xlab = "", ylab = "", xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=3, cex.sub=2, col='darkgrey', main = "Betacoronavirus")
beta_col<- cut(beta_hist$breaks, c(-Inf, 0.09834776, Inf))

Bordetella<-Nasal_Pilot_Data$BORD_CT_invCT
bord_hist<-hist(Bordetella[!Bordetella<0.051806], xlab = "", xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=3, cex.sub=2, col='darkgrey', main = "Bordetella")
bord_col<- cut(bord_hist$breaks, c(-Inf, 0.10361191, Inf))

Parainfluenza<-Nasal_Pilot_Data$PINF_CT_invCT
para_hist<-hist(Parainfluenza[!Parainfluenza<0.063173], xlab = "1/CT (Shedding)",ylab = "", xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=3, cex.sub=2, col='darkgrey', main = "Parainfluenza")
para_col<- cut(para_hist$breaks, c(-Inf, 0.12634558, Inf))

Adenovirus<-Nasal_Pilot_Data$CADEN_CT_invCT
adeno_hist<-hist(Adenovirus[!Adenovirus<0.047086], xlab = "",ylab="", xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=3, cex.sub=2, col='darkgrey', main = "Adenovirus")
adeno_col<- cut(adeno_hist$breaks, c(-Inf, 0.09417083, Inf))


graphics.off()
par(mfrow=c(2,3))
plot(myco_hist,col = c("gold","darkorange1","red2")[myco_col],xlab = "",xlim = c(0,0.45), cex.lab=2, cex.axis=2, cex.main=3, cex.sub=2, main = "Mycoplasma")
plot(pneumo_hist,col = c("gold","darkorange1","red2")[pneumo_col],xlab = "",xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=3, cex.sub=2, main = "Pneumovirus")
plot(beta_hist,col = c("gold","darkorange1","red2")[beta_col],xlab = "",xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=3, cex.sub=2, main = "Betacoronavirus")
plot(bord_hist,col = c("gold","darkorange1","red2")[bord_col],xlab = "",xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=3, cex.sub=2, main = "Bordetella")
plot(para_hist,col = c("gold","darkorange1","red2")[para_col],xlab = "",xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=3, cex.sub=2, main = "Parainfluenza")
plot(adeno_hist,col = c("gold","darkorange1","red2")[adeno_col],xlab = "",xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=3, cex.sub=2, main = "Adenovirus")


graphics.off()
par(mfrow=c(2,3))
plot(myco_hist,col = c("darkgrey","red")[myco_col],xlab = "",xlim = c(0,0.45), cex.lab=2, cex.axis=2, cex.main=3, cex.sub=2, main = "Mycoplasma")
plot(pneumo_hist,col = c("darkgrey","red")[pneumo_col],xlab = "",xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=3, cex.sub=2, main = "Pneumovirus")
plot(beta_hist,col = c("darkgrey","red")[beta_col],xlab = "",xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=3, cex.sub=2, main = "Betacoronavirus")
plot(bord_hist,col = c("darkgrey","red")[bord_col],xlab = "",xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=3, cex.sub=2, main = "Bordetella")
plot(para_hist,col = c("darkgrey","red")[para_col],xlab = "",xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=3, cex.sub=2, main = "Parainfluenza")
plot(adeno_hist,col = c("darkgrey","red")[adeno_col],xlab = "",xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=3, cex.sub=2, main = "Adenovirus")

graphics.off()











##   Extras   ##

Nasal_Pilot_Data$Shannon <- diversity(Nasal_Pilot_Data[,33:39], index = "shannon")

Nasal_Pilot_Data$Index<- seq.int(nrow(Nasal_Pilot_Data))

Nasal_Pilot_Data<-Nasal_Pilot_Data[order(Nasal_Pilot_Data$Shannon),]
Nasal_Pilot_Data$Index <- as.numeric(rownames(Nasal_Pilot_Data))
Nasal_Pilot_Data$Index <- seq.int(nrow(Nasal_Pilot_Data))

ggplot(data=Nasal_Pilot_Data) +
  geom_point(aes(x=Nasal_Pilot_Data$Index,y=Nasal_Pilot_Data$Shannon, colour=Origin), shape=19, size=2) +
  labs(x="Nasal Sample") +
  labs(y="Shannon") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black", size = 0.5, linetype = "solid")) +
  theme(axis.text.x = element_text(face="plain", color="black", size=14, angle=0), axis.text.y = element_text(face="plain", color="black", size=14, angle=0)) +
  theme(axis.title = element_text(color="black", face="bold", size=18)) +
  theme(legend.position="bottom") +
  theme(legend.title = element_text(colour="black", size=14, face="bold")) +
  theme(legend.text = element_text(colour="black", size=14, face="plain")) +
  scale_x_continuous(limits=c(0, 105)) +
  scale_y_continuous(limits=c(0,2))


par(mfrow=c(1,2))
plot(Nasal_Pilot_Data$Shannon,Nasal_Pilot_Data$MCYN_CT_invCT, col = Nasal_Pilot_Data$IsSymptomatic, pch = 16, xlab = "Diversity",
     ylab = "Mycoplasma 1/CT", ylim=c(0,0.44), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(Nasal_Pilot_Data$Shannon,Nasal_Pilot_Data$PNVPCR_CT_invCT, col = Nasal_Pilot_Data$IsSymptomatic, pch = 16, xlab = "Diversity",
     ylab = "Pneumovirus 1/CT", ylim=c(0,0.44), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)

graphics.off()
plot(Nasal_Pilot_Data$Shannon,Nasal_Pilot_Data$Total_Load, col=Nasal_Pilot_Data$Coinfection, xlab = "Diversity", ylab="Total Load", pch = 16,
     cex.lab=1.5, cex.axis=1.5, cex.main=2, cex.sub=2)

plot(Nasal_Pilot_Data$Shannon,Nasal_Pilot_Data$Total_Load, col=Nasal_Pilot_Data$IsSymptomatic, xlab = "Shannon Diversity", ylab="Total Load", pch = 16,
     cex.lab=1.5, cex.axis=1.5, cex.main=2, cex.sub=2)


plot(Nasal_Pilot_Data$Shannon,Nasal_Pilot_Data$MCYN_CT_invCT, col=Nasal_Pilot_Data$IsSymptomatic, xlab = "Shannon Diversity", ylab="Mycoplasma 1/CT", pch = 16,
     cex.lab=1.5, cex.axis=1.5, cex.main=2, cex.sub=2)

quartz("Shannon Diversity by Pneumovirus Load", 6,6)
plot(Nasal_Pilot_Data$Shannon,Nasal_Pilot_Data$PNVPCR_CT_invCT, col=Nasal_Pilot_Data$IsSymptomatic, xlab = "Shannon Diversity", ylab="Pneumovirus 1/CT", pch = 16,
     cex.lab=1.5, cex.axis=1.5, cex.main=2, cex.sub=2)
graphics.off()

plot(Nasal_Pilot_Data$Shannon,Nasal_Pilot_Data$MCYN_CT_invCT, col = Nasal_Pilot_Data$IsSymptomatic, pch = 15, xlim = c(0,1.5), ylim = c(0,0.43), xlab = "", ylab = "",cex.lab=1.5, cex.axis=1.5, cex.main=2, cex.sub=2)
par(new=TRUE)
plot(Nasal_Pilot_Data$Shannon,Nasal_Pilot_Data$PNVPCR_CT_invCT, col = Nasal_Pilot_Data$IsSymptomatic, pch = 4, xlim = c(0,1.5), ylim = c(0,0.43), xlab = "Diversity", 
     ylab = "1/CT",cex.lab=1.5, cex.axis=1.5, cex.main=2, cex.sub=2)
legend(x="topright", 
       legend = paste("Symptomatic", levels(Nasal_Pilot_Data$IsSymptomatic)),
       col = c('black','red'),
       pch = 16,
       cex = 0.7)


legend("topright",legend=paste(rep(c("Symptomatic","Medium","Large","Big"),times=2),rep(c("Site 1","Site 2"),each=4),sep=", "),col=rep(cols,times=2),pch=rep(c(16,18),each=4),bty="n",
       ncol=1,cex=0.7,pt.cex=0.7)

plot(Nasal_Pilot_Data$Evenness,Nasal_Pilot_Data$PNVPCR_CT_invCT, col = c("black","red")[Nasal_Pilot_Data$IsSymptomatic], 
     pch = 19, xlim = c(0.4,1), ylim = c(0,0.43), xlab = "Evenness", 
     ylab = "",cex.lab=3, cex.axis=3, cex.main=3, cex.sub=3, cex = 4, lwd = 3)

par(new=TRUE)


plot(Nasal_Pilot_Data$Evenness,Nasal_Pilot_Data$PNVPCR_CT_invCT, col = c("blue","red")[Nasal_Pilot_Data$IsSymptomatic], 
     pch = 7, xlim = c(0.4,1), ylim = c(0,0.43), xlab = "Evenness", 
     ylab = "1/CT",cex.lab=3, cex.axis=3, cex.main=3, cex.sub=3, cex = 4, lwd = 3)



#### Correlation matrices #####

#without MS2 and SZ because they were NA because of singularity
cor_matrix <- cor(Nasal_Pilot_Data[,33:39], y = NULL)
cor_matrix
res1 <- cor.mtest(cor_matrix, conf.level = .95)

corrplot(cor_matrix,na.label = "NA",  method = "color", p.mat = res1$p, insig = "p-value", sig.level = -1)

res_all <- cor.mtest(cor_matrix, conf.level = .95)
corrplot(cor_matrix, method = "color",  p.mat = res_all$p, insig = "p-value", sig.level = -1)


#### Subset each shelter for individual correlation matrices ######

#Madera

Madera_subset <- subset(Nasal_Pilot_Data, Origin == 'Madera', select = Well:IsSymptomatic)
madera_matrix <- cor(Madera_subset[,c(33:35,37:39)], y = NULL)
madera_matrix
corrplot(madera_matrix, method = "color")

res_madera <- cor.mtest(madera_matrix, conf.level = .95)
corrplot(madera_matrix, method = "color",  p.mat = res_madera$p, insig = "p-value", sig.level = -1)


#Merced

Merced_subset <- subset(Nasal_Pilot_Data, Origin == 'Merced', select = Well:IsSymptomatic)
merced_matrix <- cor(Merced_subset[,c(33:35,37:39)], y = NULL)
merced_matrix

res_merced <- cor.mtest(merced_matrix, conf.level = .95)
corrplot(merced_matrix, method = "color",  p.mat = res_merced$p, insig = "p-value", sig.level = -1)

# Oakland 

Oakland_subset <- subset(Nasal_Pilot_Data, Origin == 'Oakland', select = Well:IsSymptomatic)
oakland_matrix <- cor(Oakland_subset[,c(34:37,39)], y = NULL)
oakland_matrix

colnames(Oakland_subset)
res_oak <- cor.mtest(oakland_matrix, conf.level = .95)
corrplot(oakland_matrix, method = "color",  p.mat = res_oak$p, insig = "p-value", sig.level = -1)

#Fresno 

Fresno_subset <- subset(Nasal_Pilot_Data, Origin == 'Fresno', select = Well:IsSymptomatic)
fresno_matrix <- cor(Fresno_subset[,33:39], y = NULL)
fresno_matrix

res_fresno <- cor.mtest(fresno_matrix, conf.level = .95)
corrplot(fresno_matrix, method = "color", na.label = "NA",  p.mat = res_fresno$p, insig = "p-value", sig.level = -1)

#Sacramento

Sacramento_subset <- subset(Nasal_Pilot_Data, Origin == 'Sacramento', select = Well:IsSymptomatic)
sacramento_matrix <- cor(Sacramento_subset[,c(34:37,39)], y = NULL)
sacramento_matrix

res_sac <- cor.mtest(sacramento_matrix, conf.level = .95)
corrplot(sacramento_matrix, method = "color", na.label = "NA",  p.mat = res_sac$p, insig = "p-value", sig.level = -1)


#Hayward

Hayward_subset <- subset(Nasal_Pilot_Data, Origin == 'Hayward', select = Well:IsSymptomatic)
hayward_matrix <- cor(Hayward_subset[,33:39], y = NULL)
hayward_matrix
corrplot(hayward_matrix, method = "color")

#Theres only 3 dogs from Hayward, 1 with infection
#View(Hayward_subset)

# OHS

OHS_subset <- subset(Nasal_Pilot_Data, Origin == 'OHS', select = Well:IsSymptomatic)
ohs_matrix <- cor(OHS_subset[,33:39], y = NULL)
ohs_matrix

res_ohs <- cor.mtest(ohs_matrix, conf.level = .95)
corrplot(ohs_matrix, method = "color", na.label = "NA",  p.mat = res_ohs$p, insig = "p-value", sig.level = -1)

# Oregon

oregon_subset <- subset(Nasal_Pilot_Data, State == 'Oregon', select = Well:IsSymptomatic)
oregon_matrix <- cor(oregon_subset[,33:39], y = NULL)
oregon_matrix

res_ore <- cor.mtest(oregon_matrix, conf.level = .95)
corrplot(oregon_matrix,na.label = "NA",  method = "color",  p.mat = res_ore$p, insig = "p-value", sig.level = -1)


# California

cali_subset <- subset(Nasal_Pilot_Data, State == 'California', select = Well:IsSymptomatic)
cali_matrix <- cor(cali_subset[,33:39], y = NULL)
cali_matrix

corrplot(cali_matrix, method = "color")

res_cali <- cor.mtest(cali_matrix, conf.level = .95)
corrplot(cali_matrix,na.label = "NA",  method = "color",  p.mat = res_cali$p, insig = "p-value", sig.level = -1)



graphics.off()
par(mfrow=c(2,3))

Mycoplasma<-Nasal_Pilot_Data$MCYN_CT_invCT
myco_hist<- hist(Mycoplasma[!Mycoplasma==0], xlab = "", xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=3, cex.sub=2, main = "Mycoplasma", col='darkgrey')
myco_col<- cut(myco_hist$breaks, c(-Inf,0.060887, 0.12177302, Inf))

Pneumovirus<-Nasal_Pilot_Data$PNVPCR_CT_invCT
pneumo_hist<- hist(Pneumovirus[!Pneumovirus==0], xlab = "1/CT (Shedding)", ylab = "", xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=3, cex.sub=2, col='darkgrey', main = "Pneumovirus")
pneumo_col<- cut(pneumo_hist$breaks, c(-Inf,0.05445, 0.10890043, Inf))

Betacoronavirus<-Nasal_Pilot_Data$BCOR_CT_invCT
beta_hist<-hist(Betacoronavirus[!Betacoronavirus==0], xlab = "", ylab = "", xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=3, cex.sub=2, col='darkgrey', main = "Betacoronavirus")
beta_col<- cut(beta_hist$breaks, c(-Inf,0.04918, 0.09834776, Inf))

Bordetella<-Nasal_Pilot_Data$BORD_CT_invCT
bord_hist<-hist(Bordetella[!Bordetella==0], xlab = "", xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=3, cex.sub=2, col='darkgrey', main = "Bordetella")
bord_col<- cut(bord_hist$breaks, c(-Inf,0.051806, 0.10361191, Inf))

Parainfluenza<-Nasal_Pilot_Data$PINF_CT_invCT
para_hist<-hist(Parainfluenza[!Parainfluenza==0], xlab = "1/CT (Shedding)",ylab = "", xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=3, cex.sub=2, col='darkgrey', main = "Parainfluenza")
para_col<- cut(para_hist$breaks, c(-Inf,0.063173, 0.12634558, Inf))

Adenovirus<-Nasal_Pilot_Data$CADEN_CT_invCT
adeno_hist<-hist(Adenovirus[!Adenovirus==0], xlab = "",ylab="", xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=3, cex.sub=2, col='darkgrey', main = "Adenovirus")
adeno_col<- cut(adeno_hist$breaks, c(-Inf,0.047086, 0.09417083, Inf))


graphics.off()
par(mfrow=c(2,3))
plot(myco_hist,col = c("gold","darkorange1","red2")[myco_col],xlab = "",ylab="",xlim = c(0,0.45), cex.lab=2, cex.axis=2, cex.main=3, cex.sub=2, main = "Mycoplasma")
plot(pneumo_hist,col = c("gold","darkorange1","red2")[pneumo_col],xlab = "",ylab="",xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=3, cex.sub=2, main = "Pneumovirus")
plot(beta_hist,col = c("gold","darkorange1","red2")[beta_col],xlab = "",ylab="",xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=3, cex.sub=2, main = "Betacoronavirus")
# #legend(x="topright", 
#        legend = paste(c('Low Positive','Moderate Positive','High Positive')),
#        col = c("gold","darkorange1","red2"),
#        pch = 15,
#        inset = 0.02,
#        cex = 1.75,
#        pt.cex=2.5)
plot(bord_hist,col = c("gold","darkorange1","red2")[bord_col],xlab = "",ylab="",xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=3, cex.sub=2, main = "Bordetella")
plot(para_hist,col = c("gold","darkorange1","red2")[para_col],xlab = "",ylab="",xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=3, cex.sub=2, main = "Parainfluenza")
plot(adeno_hist,col = c("gold","darkorange1","red2")[adeno_col],xlab = "",ylab="",xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=3, cex.sub=2, main = "Adenovirus")




infected<- length(which(Nasal_Pilot_Data$IsSymptomatic==TRUE))
print(infected)




panel <- Nasal_Pilot_Data


################################




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