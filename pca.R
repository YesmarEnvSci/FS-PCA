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
fs_hab <- read.csv("fs_habitat.csv")
hab <-fs_hab[c(1,6:9,12:14,16:18,20:21,26:28,30:31)]
hab <-na.omit(hab)
rm(fs_hab)
ave.hab <-hab%>%
  group_by(Site)%>%
  summarise_each(funs(mean))
  
hab.1 <-ave.hab[-c(1,2)]
hab.l <-scale(log(hab.1+1))
boxplot(scale(hab.1))
boxplot(hab.l)

#-------
#Data for analysis that includes relative fish abundance
hb_fish <-left_join(hab,fish,by="Site")
hbfish <-hb_fish[-c(1,3,18:25)]
rm(hb_fish)
hbfish <-na.omit(hbfish)
hbfish.lg <- log(hbfish+1) #log transformed data
  
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

plot(pca, main="Scree Plot") #Scree plot
broken.stick(16) #After comparing, keep comp 1 & 2

round(loadings(pca),2) #Check eigenvectors: length of vector is relative variance and how much it contributes to the PC
#Principal component loading (pg 50).  The further from zero, the greater the contribution.
round(loadings(pca)[,c(1:2)],2) #Loading for PC1 & 2 only

round((pca$scores),2) #PC matrix showing site scores for all PCs. How far each is(SD) from the the grand centroid
#This is the distribution of PC1 and PC2 site scores (top scale).  Each variable for each component. 
#In this case due to broken stick, PC1 and PC2

#-----
#channel type plot
plot(pca$scores[,1], pca$scores[,2],xlab="PC 1", ylab="PC 2", type='n') # plot the first two PCs.  Can you interpret the plot?
text(pca$scores[,1], pca$scores[,2],labels=smry$Chan.Type, lwd=2)
plot(pca$scores[,1], pca$scores[,2],xlab="PC 1", ylab="PC 2", cex=3) # plot the first two PCs

#---------
#create Shepard diagram
head(round(pca$scores,2))
euc<-dist(hab.l)#Computes the distance matrix in the log transformed data (multidimensional space)
euc.1<-dist(pca$scores[,c(1,2)]) #Computes the distance matrix in the PC matrix (reduced space)
plot(euc,euc.1,main="Shepards Diagram (PC=2)", xlab="Distance in Multidimensional space", ylab="Distance in Reduced space")   #Shepard diagram

