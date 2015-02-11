#USFS habitat and fish pop analysis

fs_fish <- read.csv("fs_fish.csv")
fs_lgh <- read.csv("fs_fish_length.csv")
fs_smry <- read.csv("fs_fish_smry.csv")
fs_hab <- read.csv("fs_habitat.csv")

#CT & DV
fish<-fs_fish[c(1:3,17:18)]
rm(fs_fish)
fish <-na.omit(fish)
library(dplyr)
ave.rel.den <- fish%>% # combines species for mean relative fish density
  group_by(Site)%>%
  summarise(ave.den=mean(Species_Density..fish.100sq.m.))

#Fish length
length <- fs_lgh[-c(2:3,6:7)]#2009 fish
rm(fs_lgh)

#Habitat
hab <-fs_hab[c(1,7:9,12:14,16:18,20:21,26:28,30:31)]
hab <-na.omit(hab)
rm(fs_hab)
str(hab)
ave.hab <-hab%>%
  group_by(Site)%>%
  summarise_each(funs(mean))
rm(hab)  
hb <-ave.hab[-c(1)]

boxplot(scale(hb))
boxplot(scale(log(hb+1)))
cor.matrix(scale(hb))

cor(log(hb+1))
hb.lg <- log(hb+1)

#Summary
smry <-fs_smry[c(1,4,6:7,11:15)]
rm(fs_smry)
all <-left_join(smry, ave.rel.den) #includes ave density of fish
rm(smry)

lshap <- lapply(hb.lg, shapiro.test) #shapiro test on log transformed data
lres <- sapply(lshap, `[`, c("statistic","p.value"))
t(lres)



#----------

#PCA 
require(MASS) #loads the PCA package
pca <- princomp(scale(hb.lg)) #creates a PC matrix using the correlation matrix
biplot(pca, expand = 1.05,main = "Biplot", xlab = "Comp.1 (37.9%)", ylab = "Comp.2 (20.4%)")
#Scale for sites(PC matrix-pca$scores) on top, scale for variables (vectors-loadings) along bottom
summary(pca) #proportion of variance is eigenvalues for each PC

plot(pca, main="Scree Plot") #Scree plot
broken.stick(18) #After comparing, keep comp 1 & 2

round(loadings(pca),2) #Check eigenvectors: length of vector is relative variance and how much it contributes to the PC
#Principal component loading (pg 50).  The further from zero, the greater the contribution.
round(loadings(pca)[,c(1:2)],2) #Loading for PC1 & 2 only

round((pca$scores),2) #PC matrix showing site scores for all PCs. How far each is(SD) from the the grand centroid
#This is the distribution of PC1 and PC2 site scores (top scale).  Each variable for each component. 
#In this case due to broken stick, PC1 and PC2

#---------
#create shepard diagram
wtr.d<-round(var(scale(hb.lg,scale=F)),0)  #calculate variance-covariance matrix and save it to 'wtr.d'
e.1<-eigen(wtr.d) #eigen-analysis
pc.matrix.1<-(as.matrix(scale(hb.lg, scale=F)))%*%e.1$vectors #calculate pc.matrix
euc<-dist(scale(hb.lg, scale=F))  #calculate Euclidian distance among site for centered original data
round(euc,2)
euc.1<-dist(pc.matrix.1[,c(1,2)])  #calculate Euclidian distance among sites in PCA space using only first 2 PCs
round(euc.1,2)
plot(euc,euc.1,main="Shepard diagram", xlab="Distance in Multidimensional space", ylab="Distance in Reduced space")
#Farther from y=x the more distortion
