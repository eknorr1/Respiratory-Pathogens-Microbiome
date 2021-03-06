# Respiratory Pathogen Microbiome Pilot Analysis 
# Respiratory Pathogen Community Analysis
# Liz Knorr
# created 4/24/2019

# Clear workspace ------------------------------------------------------

rm(list = ls())
graphics.off()
palette("default")


# Load Libraries -------------------------------------------------------

library(vegan)
library(corrplot)
library(factoextra)
library(FactoMineR)

# Script Parameters ----------------------------------------------------

#set to where data files are located
DataDirectory <- "~/Desktop/Pilot Data/Raw Data"


# Load dataset ------------------------------------------------------------

setwd(DataDirectory)

panel <- read.csv(file = "panel.csv", stringsAsFactors = F)

panel$IsSymptomatic <- as.integer(panel$IsSymptomatic)



# Species richness plot ------------------------------------------------------

Species_Richness <- sort(panel$Coinfection,decreasing = FALSE)

plot(Species_Richness)



# Shannon diversity plots ---------------------------------------------------------

Shannon <- diversity(panel[,3:8], index = "shannon")  # note hard coding for the columns of inv Ct scores
hist(Shannon) 

Shan_increasing <- sort(panel$Shannon,decreasing = FALSE)
plot(Shan_increasing)

par(fin = c(4,4))
par(mai = c(0.8,0.8,0.8,0.8))
plot(panel$Shannon,panel$Total_Load, col=panel$Coinfection)





# Evenness Plots ----------------------------------------------------------------

# Evenness - Total Load colored by coinfections

#quartz("",5,5)
par(fin = c(4,4))
par(mai = c(0.8,0.8,0.8,0.8))
plot(panel$Evenness, panel$Total_Load, col = panel$Coinfection, 
     pch = 16, xlab = 'Evenness', ylab = 'Total Load')

legend(x="bottomleft", 
       legend = paste(c('2','3','4','5','6')),
       col = c('red','green','blue','cyan','magenta'),
       title = 'Infections',
       pch = 16,
       inset = 0.02,
       cex = 0.9
       )

# Evenness - Total Load colored by symptomatic status
graphics.off()
quartz("",7,6)
par(fin = c(4,4))
par(mai = c(0.8,0.8,0.8,0.8))
ggplot(data=panel) +
  geom_point(aes(x=panel$Evenness,y=panel$Total_Load, colour = Symptomatic), shape=19, size=5) +
  labs(x="Evenness") +
  labs(y="Total Load") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black", size = 0.5, linetype = "solid")) +
  theme(axis.text.x = element_text(face="plain", color="black", size=14, angle=0), 
        axis.text.y = element_text(face="plain", color="black", size=14, angle=0)) +
  theme(axis.title = element_text(color="black", face="bold", size=18)) +
  theme(legend.position="bottom") +
  theme(panel.grid = element_blank()) +
  theme(legend.title = element_text(colour="black", size=14, face="bold")) +
  theme(legend.text = element_text(colour="black", size=14, face="plain")) +
  scale_x_continuous(limits=c(0.5, 1)) +
  scale_y_continuous(limits=c(0,0.7)) +
  scale_color_manual(values = c("black","red3")
  )




# Mycoplasma Plots ---------------------------------------------------------------------

panel$Index<- seq.int(nrow(panel))

panel<-panel[order(panel$MCYN),]
panel$Index <- as.numeric(rownames(panel))
panel$Index <- seq.int(nrow(panel))

# Mycoplasma abundance in each sample colored by symptomatic
quartz("",7,6)
ggplot(data=panel) +
  geom_point(aes(x=panel$Index,y=panel$MCYN, colour = Symptomatic), shape=19, size=5) +
  labs(x="Nasal Sample") +
  labs(y="Mycoplasma Shedding") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black", size = 0.5, linetype = "solid")) +
  theme(axis.text.x = element_text(face="plain", color="black", size=14, angle=0), 
        axis.text.y = element_text(face="plain", color="black", size=14, angle=0)) +
  theme(axis.title = element_text(color="black", face="bold", size=18)) +
  theme(legend.position="bottom") +
  theme(panel.grid = element_blank()) +
  theme(legend.title = element_text(colour="black", size=14, face="bold")) +
  theme(legend.text = element_text(colour="black", size=14, face="plain")) +
  scale_x_continuous(limits=c(0, 105)) +
  scale_y_continuous(limits=c(0,0.5)) +
  scale_color_manual(values = c("darkgrey","red3")
  )

# Mycoplasma - Evenness colored by symptomatic 
quartz("",7,6)
ggplot(data=panel) +
  geom_point(aes(x=panel$Evenness,y=panel$MCYN, colour = Symptomatic), shape=19, size=5) +
  labs(x="Evenness") +
  labs(y="Mycoplasma") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black", size = 0.5, linetype = "solid")) +
  theme(axis.text.x = element_text(face="plain", color="black", size=14, angle=0), 
        axis.text.y = element_text(face="plain", color="black", size=14, angle=0)) +
  theme(axis.title = element_text(color="black", face="bold", size=18)) +
  theme(legend.position="bottom") +
  theme(panel.grid = element_blank()) +
  theme(legend.title = element_text(colour="black", size=14, face="bold")) +
  theme(legend.text = element_text(colour="black", size=14, face="plain")) +
  scale_x_continuous(limits=c(0.5, 1)) +
  scale_y_continuous(limits=c(0,0.44)) +
  scale_color_manual(values = c("black","red3")
  )




# Pneumovirus Evenness Plot - colored by symptomatic -------------------------------------------

ggplot(data=panel) +
  geom_point(aes(x=panel$Evenness,y=panel$PNVPCR, colour = Symptomatic), shape=19, size=5) +
  labs(x="Evenness") +
  labs(y="Pneumovirus") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black", size = 0.5, linetype = "solid")) +
  theme(axis.text.x = element_text(face="plain", color="black", size=14, angle=0), 
        axis.text.y = element_text(face="plain", color="black", size=14, angle=0)) +
  theme(axis.title = element_text(color="black", face="bold", size=18)) +
  theme(legend.position="bottom") +
  theme(panel.grid = element_blank()) +
  theme(legend.title = element_text(colour="black", size=14, face="bold")) +
  theme(legend.text = element_text(colour="black", size=14, face="plain")) +
  scale_x_continuous(limits=c(0.5, 1)) +
  scale_y_continuous(limits=c(0,0.44)) +
  scale_color_manual(values = c("black","red3")
  )




# Boxplots --------------------------------------------------------------------------------

graphics.off()

# Number of coinfections in symptomatic and asymptomatic dogs all samples
par(fin = c(4,4))
par(mai = c(0.8,0.8,0.8,0.8))
boxplot(Coinfection~Symptomatic, data = panel, ylab = "", names=c("Asymptomatic","Symptomatic"), col=c("blue","red"))

# Difference in coinfections from symptomatic/asymptomatic dogs shedding any pathogen
par(fin = c(4,4))
par(mai = c(0.8,0.8,0.8,0.8))
boxplot(panel$Coinfection[panel$Total_Load>0]~IsSymptomatic[panel$Total_Load>0],
        data = panel, ylab = "", xlab = "",names=c("Symptomatic","Asymptomatic"), col=c("royalblue3","red3"))

graphics.off()
boxplot(panel$Coinfection[panel$Total_Load>0]~IsSymptomatic[panel$Total_Load>0],
        data = panel, ylab = "", xlab = "",names=c("Symptomatic","Asymptomatic"), col=c("royalblue3","red3"))

#Difference in mycoplasma loads from symptomatic and asymptomatic dogs shedding myco 
par(fin = c(4,4))
par(mai = c(0.7,0.7,0.7,0.7))
boxplot(panel$MCYN[panel$MCYN>0]~IsSymptomatic[panel$MCYN>0],
        data = panel, ylab = "", xlab = "", names=c("a","b"), col=c("royalblue3","red3"))

#Difference in mycoplasma loads from symptomatic and asymptomatic dogs (ALL)
par(mgp=c(5,1,0))
par(mar=c(2,5,2,2))
boxplot(panel$MCYN[panel$Total_Load>0]~IsSymptomatic[panel$Total_Load>0],
        data = panel, ylab = "Mycoplasma Load", xlab = "", 
        names=c("Asymptomatic","Symptomatic"), col=c("royalblue3","red3"))




# Histograms without zeros and colored by low, moderate, and high positives ------------------------

graphics.off()
par(mfrow=c(2,3))

Mycoplasma<-panel$MCYN
myco_hist<- hist(Mycoplasma[!Mycoplasma==0], xlab = "", xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, main = "Mycoplasma", col='darkgrey')
myco_col<- cut(myco_hist$breaks, c(-Inf,0.060887, 0.12177302, Inf))

Pneumovirus<-panel$PNVPCR
pneumo_hist<- hist(Pneumovirus[!Pneumovirus==0], xlab = "1/CT (Shedding)", ylab = "", xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, col='darkgrey', main = "Pneumovirus")
pneumo_col<- cut(pneumo_hist$breaks, c(-Inf,0.05445, 0.10890043, Inf))

Betacoronavirus<-panel$BCOR
beta_hist<-hist(Betacoronavirus[!Betacoronavirus==0], xlab = "", ylab = "", xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, col='darkgrey', main = "Betacoronavirus")
beta_col<- cut(beta_hist$breaks, c(-Inf,0.04918, 0.09834776, Inf))

Bordetella<-panel$BORD
bord_hist<-hist(Bordetella[!Bordetella==0], xlab = "", xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, col='darkgrey', main = "Bordetella")
bord_col<- cut(bord_hist$breaks, c(-Inf,0.051806, 0.10361191, Inf))

Parainfluenza<-panel$PINF
para_hist<-hist(Parainfluenza[!Parainfluenza==0], xlab = "1/CT (Shedding)",ylab = "", xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, col='darkgrey', main = "Parainfluenza")
para_col<- cut(para_hist$breaks, c(-Inf,0.063173, 0.12634558, Inf))

Adenovirus<-panel$CADEN
adeno_hist<-hist(Adenovirus[!Adenovirus==0], xlab = "",ylab="", xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, col='darkgrey', main = "Adenovirus")
adeno_col<- cut(adeno_hist$breaks, c(-Inf,0.047086, 0.09417083, Inf))


graphics.off()
par(mfrow=c(2,3))
plot(myco_hist,col = c("gold","darkorange1","red2")[myco_col],xlab = "",ylab="",xlim = c(0,0.45), cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, main = "Mycoplasma")
plot(pneumo_hist,col = c("gold","darkorange1","red2")[pneumo_col],xlab = "",ylab="",xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, main = "Pneumovirus")
plot(beta_hist,col = c("gold","darkorange1","red2")[beta_col],xlab = "",ylab="",xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, main = "Betacoronavirus")
# #legend(x="topright",
#        legend = paste(c('Low Positive','Moderate Positive','High Positive')),
#        col = c("gold","darkorange1","red2"),
#        pch = 15,
#        inset = 0.02,
#        cex = 1.75,
#        pt.cex=2.5)
plot(bord_hist,col = c("gold","darkorange1","red2")[bord_col],xlab = "",ylab="",xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, main = "Bordetella")
plot(para_hist,col = c("gold","darkorange1","red2")[para_col],xlab = "",ylab="",xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, main = "Parainfluenza")
plot(adeno_hist,col = c("gold","darkorange1","red2")[adeno_col],xlab = "",ylab="",xlim = c(0,0.43), cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, main = "Adenovirus")




#GLMs -------------------------------------------------------------------------------------------

panel$Origin <- factor(as.character(panel$Origin))

fit <- glm(Symptomatic ~ Total_Load, data = panel, family = binomial(link = "logit"))
summary(fit)

myco_total_load <- glm(Symptomatic ~ Total_Load + MCYN, data = panel,
                       family = binomial(link = "logit"))
summary(myco_total_load)

Sympt_fit2 <- glm(IsSymptomatic ~ CADEN + BCOR + BORD + MCYN + PINF + PNVPCR, 
                  data = panel, family = binomial(link = "logit"))
summary(Sympt_fit2)

fit3 <- glm(IsSymptomatic ~ MCYN, data = panel, family = binomial(link = "logit"))
summary(fit3)





# Some functions for information on samples ---------------------------------------------------

dog_uniq<- unique(panel$PetName)
length(dog_uniq)

infected<- length(which(panel$Total_Load>0))
print(infected)



# Correlation matrices -------------------------------------------------------------------------
 
# #without MS2 and SZ: NA because of singularity

graphics.off()
cor_matrix <- cor(panel[,3:8], y = NULL)
cor_matrix
res_all <- cor.mtest(cor_matrix, conf.level = .95)
corrplot(cor_matrix, method = "color",  p.mat = res_all$p, insig = "p-value", sig.level = -1)


# Subset each shelter for individual correlation matrices ----------------------------------------

#Madera

Madera_subset <- subset(panel, Origin == 'Madera', select = Well:PNVPCR)
madera_matrix <- cor(Madera_subset[,c(2:7)], y = NULL)
madera_matrix
res_madera <- cor.mtest(madera_matrix, conf.level = .95)
corrplot(madera_matrix, method = "color",  p.mat = res_madera$p, insig = "p-value", sig.level = -1)


#Merced

Merced_subset <- subset(panel, Origin == 'Merced', select = Well:PNVPCR)
merced_matrix <- cor(Merced_subset[,c(2:7)], y = NULL)
merced_matrix
res_merced <- cor.mtest(merced_matrix, conf.level = .95)
corrplot(merced_matrix, method = "color",  p.mat = res_merced$p, insig = "p-value", sig.level = -1)

# Oakland

Oakland_subset <- subset(panel, Origin == 'Oakland', select = Well:PNVPCR)
oakland_matrix <- cor(Oakland_subset[,c(2:7)], y = NULL)
oakland_matrix
res_oak <- cor.mtest(oakland_matrix, conf.level = .95)
corrplot(oakland_matrix, method = "color",  p.mat = res_oak$p, insig = "p-value", sig.level = -1)

#Fresno

Fresno_subset <- subset(panel, Origin == 'Fresno', select = Well:PNVPCR)
fresno_matrix <- cor(Fresno_subset[,2:7], y = NULL)
fresno_matrix
res_fresno <- cor.mtest(fresno_matrix, conf.level = .95)
corrplot(fresno_matrix, method = "color", na.label = "NA",  p.mat = res_fresno$p, insig = "p-value", sig.level = -1)

#Sacramento

Sacramento_subset <- subset(panel, Origin == 'Sacramento', select = Well:PNVPCR)
sacramento_matrix <- cor(Sacramento_subset[,c(2:7)], y = NULL)
sacramento_matrix
res_sac <- cor.mtest(sacramento_matrix, conf.level = .95)
corrplot(sacramento_matrix, method = "color", na.label = "NA",  p.mat = res_sac$p, insig = "p-value", sig.level = -1)


#Hayward only 5 samples

Hayward_subset <- subset(panel, Origin == 'Hayward', select = Well:PNVPCR)
hayward_matrix <- cor(Hayward_subset[,2:7], y = NULL)
hayward_matrix
corrplot(hayward_matrix, method = "color")













# ##   Extras   ## -----------

# Extra for when there are samples from California and Oregon -----------------------------------

# OHS
# 
# OHS_subset <- subset(panel, Origin == 'OHS', select = Well:PNVPCR)
# ohs_matrix <- cor(OHS_subset[,2:7], y = NULL)
# ohs_matrix
# 
# res_ohs <- cor.mtest(ohs_matrix, conf.level = .95)
# corrplot(ohs_matrix, method = "color", na.label = "NA",  p.mat = res_ohs$p, insig = "p-value", sig.level = -1)
# 

# # Oregon
# 

# oregon_subset <- subset(panel, State == 'Oregon', select = Well:PNVPCR)
# oregon_matrix <- cor(oregon_subset[,2:7], y = NULL)
# oregon_matrix
# 
# res_ore <- cor.mtest(oregon_matrix, conf.level = .95)
# corrplot(oregon_matrix,na.label = "NA",  method = "color",  p.mat = res_ore$p, insig = "p-value", sig.level = -1)
# 
# 
# # California
# 
# cali_subset <- subset(panel, State == 'California', select = Well:PNVPCR)
# cali_matrix <- cor(cali_subset[,2:7], y = NULL)
# cali_matrix
# 
# corrplot(cali_matrix, method = "color")
# 
# res_cali <- cor.mtest(cali_matrix, conf.level = .95)
# corrplot(cali_matrix,na.label = "NA",  method = "color",  p.mat = res_cali$p, insig = "p-value", sig.level = -1)



