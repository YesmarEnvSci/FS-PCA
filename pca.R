#USFS habitat and fish pop analysis

#----
#CT & DV
fs_fish <- read.csv("fs_fish.csv")
fish<-fs_fish[c(1:3,9, 17:18)]
rm(fs_fish)
fish <-na.omit(fish)

library(dplyr)
pop <- fish%>% 
  group_by(Site, Fish.Species)%>%
  summarise(pop.ave=mean(Pop_Est))

library(tidyr)
pop <- spread(pop, Fish.Species, pop.ave)
par(mfrow=c(1,2))
hist(pop$CT)
hist(pop$DV)
pop <-pop[-c(1)]
pop.l <-log(pop+1)
par(mfrow=c(1,2))
hist(pop.l$CT)
hist(pop.l$DV)


#Fish length
fs_lgh <- read.csv("fs_fish_length.csv")
length <- fs_lgh[-c(2:3,6:7)]#2009 fish
fk.l <-length%>%
  group_by(site, species)%>%
  summarise(fk.ave=mean(fk.length..mm.))
rm(fs_lgh)
rm(length)

#Summary
fs_smry <- read.csv("fs_fish_smry.csv")
smry <-fs_smry[c(1,4,6:7,11:15)]
rm(fs_smry)

#Habitat
library(dplyr)
fs_hab <- read.csv("fs_habitat.csv")
hab <-fs_hab[c(1,6:9,12:14,16:18,20:21,26:28,30:31)]
hab <-na.omit(hab)
rm(fs_hab)
ave.hab <-hab%>%
  group_by(Site)%>%
  summarise_each(funs(mean))
  
hab.1 <-ave.hab[-c(1,2)]
hab.l <-log(hab.1+1)
boxplot(scale(hab.1))

#-------
boxplot(scale(hbfish))
boxplot(scale(log(hbfish+1)))
cor.matrix(scale(hb))

cor(log(hbfish+1))

lshap <- lapply(hbfish.lg, shapiro.test) #shapiro test on log transformed data
lres <- sapply(lshap, `[`, c("statistic","p.value"))
t(lres)

#----------
#PCA 
require(MASS) #loads the PCA package
pca <- princomp(hab.l, cor=TRUE) #creates a PC matrix using the correlation matrix
biplot(pca, expand = 1.05,main = "Biplot", xlab = "Comp.1 (38.0%)", ylab = "Comp.2 (19.0%)")
#Scale for sites(PC matrix-pca$scores) on top, scale for variables (vectors-loadings) along bottom
summary(pca) #proportion of variance is eigenvalues for each PC

library(vegan)
screeplot(pca, bstick = TRUE, main="PCA") #inertia= variance n PCA

round(loadings(pca),2) #Check eigenvectors: length of vector is relative variance and how much it contributes to the PC
#Principal component loading (pg 50).  The further from zero, the greater the contribution.
round(loadings(pca)[,c(1:2)],2) #Loading for PC1 & 2 only

round((pca$scores),2) #PC matrix showing site scores for all PCs. How far each is(SD) from the the grand centroid
#This is the distribution of PC1 and PC2 site scores (top scale).  Each variable for each component. 
#In this case due to broken stick, PC1 and PC2

#-----
#channel type plot vs fish plot
par(mfrow=c(2,2), bty="u",fg="black", bg="white")
plot(pca$scores[,1], pca$scores[,2],xlab="PC 1", ylab="PC 2", pch=5) # plot the first two PCs
plot(pca$scores[,1], pca$scores[,2],xlab="PC 1", ylab="PC 2", type='n') # plot the first two PCs.
text(pca$scores[,1], pca$scores[,2],labels=smry$Chan.Type, lwd=2)
symbols(pca$scores[,1], pca$scores[,2],circles=pop[,1],inches=0.2, main="Cutthroat Trout") #show relative concentrations of CT among sites
symbols(pca$scores[,1], pca$scores[,2],squares=pop[,2], inches=0.2, main="Dolly varden Charr") #show relative concentrations of DV among sites

#---------
#create Shepard diagram
euc<-dist(scale(hab.l)) #Calculate Euclidian distance among sites scale=centered to Z-score (multidimentional spcae). Check transformation and matrix used.
euc.1<-dist(pca$scores[,c(1,2)]) #calculate Euclidian distance among sites in PCA space using only first 2 PCs (reduced space).
plot(euc,euc.1,main="PC=2", xlab="Distance in Multidimensional space", ylab="Distance in Reduced space") #x=euc, y=euc.1  

#------
#RDA
hab <- hab.l
fsh <- pop.l
library(vegan)
library(MASS)
mod <-rda(fsh~., data=hab, scale=T) #Full model RDA

mod #Inertia is variance (Lec 6)
summary(mod) #Species scores is eigenvectors
plot(mod, type="n") #Biplot
text(mod, dis="cn", col="red")
points(mod, pch=21, col="darkgreen", bg="lightgrey", cex=1.2)
text(mod, "species", col="black", cex=0.8)
#Fish species is response variable(black), habitat explanatory variable is (red) and sites are green/grey.

#Is the RDA model relationship between the two matrices significant?  (anova.cca: global permutation test)
anova.cca(mod) #No, significant with p-value > 0.05

#Can the RDA model be reduced?(step selection with AIC plus VIF)

vif.cca(mod)#Redundancy among species (Veriance inflation factor)
# VIF > 4 or 5 suggests multi-collinearity; VIF > 10 is strong evidence 
#that collinearity is affecting the regression coefficients.

#Selection procedure (AIC approach)- Hybrid approach, search method that compares models sequentially

#Full model (with all Xs):
rda.fs<-rda(fsh ~.,data=hab,scale=T)
#Null model (with no Xs):
rda.0<-rda(fsh ~1, data=hab, scale=T)
#Hybrid selection:
rda.1<-step(rda.0, scope=formula(rda.fs))

#Run VIF again on new selection.  See if there are any multi-col > 5
#If yes, regress and drop higher ones that are correlated
hab.1 <-hab[-c(3)]
rda.1<-rda(fsh ~.,data=hab.1,scale=T)
vif.cca(rda.1) #Still high

#Full model (with all Xs):
rda.1<-rda(fsh ~.,data=hab,scale=T)
#Null model (with no Xs):
rda.0<-rda(fsh ~1, data=hab, scale=T)
#Hybrid selection:
rda.1<-step(rda.0, scope=formula(rda.fs))

#Run VIF again to see if all are < 5
vif.cca(rda.ca.2)

#Test if RDA is significant
anova.cca(rad.ca, step=1000)

#Test if RDA axis is significant
anova.cca(rad.ca,by=‘axis’, step=1000)





