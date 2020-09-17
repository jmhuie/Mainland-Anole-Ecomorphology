rm(list = ls(all=TRUE))
source("MainlandAnole_Functions.R")
library(tidyverse)
library(readxl)
library(stats)
library(phytools)
library(factoextra)



# Data Cleaning -----------------------------------------------------------

MorphoData <- read_excel("AnoleDataClean.xlsx", sheet = 1) %>% #Read in full data
#MorphoData <- read.csv("AnoleDataClean.csv", header = TRUE, sep = ",", na.strings = "") %>%
  filter(Sex == "M") %>% #Keeps only Males
  filter(is.na(Omit)) %>% #removes species marked to omit
  as.data.frame
MorphoData <- MorphoData[,c(1:2,4,7:28)]

#log the data
MorphoData[MorphoData$Species == "onca",c("Finger.II.Lam","Finger.III.Lam","Toe.II.Lam","Toe.III.Lam")] <- 1 # change onca toepad count from 0 to 1
MorphoData[,4:25] <- apply(MorphoData[,4:25], 2, log) # logs all trait values


#size-correct the data
residuals <- MorphoData
for (y in 5:21) {
  fit <-lm(MorphoData[which(!is.na(MorphoData[,y])),y]~MorphoData[which(!is.na(MorphoData[,y])),4])
  residuals[which(!is.na(MorphoData[,y])),y] <- as.array(fit$residuals)
  remove(fit)
}

#for (y in 5:25) {
#  residuals[which(!is.na(MorphoData[,y])),y] <- MorphoData[which(!is.na(MorphoData[,y])),y]/MorphoData[which(!is.na(MorphoData[,y])),"SVL"]
#}

#average species trait data
NewData <- residuals %>% 
  dplyr::group_by(Species) %>% 
  summarise_at(vars(SVL:Toe.III.Lam), list(mean), na.rm = TRUE) %>% # finds the group mean for each trait
  filter(!is.na(TL)) %>%
  as.data.frame
rownames(NewData) <- NewData$Species


NewData <- NewData[,c(1:3,5:15,17,19)]


# Read Tree ---------------------------------------------------------------

tree <- read.nexus("Trees/Poe_2017_MCC_names.txt")
tree <- ladderize(tree, right = FALSE)
tree <- drop.tip(tree, setdiff(tree$tip.label, unique(NewData$Species))) #keeps tip labels that matches vector
NewData[which(is.na((match(NewData$Species,tree$tip.label)))),1] #checks that tree tip labels match datasets species list
NewData <- NewData[tree$tip.label,] # subsets data to tree species and reorders (I think)
all(as.factor(unique(NewData[,1])) == tree$tip.label) #checks that tree tip labels match datasets species list




# PCA ---------------------------------------------------------------------

#pca <- prcomp(residuals[,8:20], scale = TRUE)
pca <- prcomp(NewData[,2:16])
pca$L
plot(pca)


fviz_pca_biplot(pca, 
                axes = c(2,3),
                geom.ind = "point",
                pointshape = 21,
                pointsize = 2.5,
                #fill.ind = NewData$Genus,
                col.ind = "black",
                #legend.title = list(fill = "Genus", color = "Clusters"),
                repel = TRUE,
                addEllipses = FALSE,
                mean.point = FALSE,
                title = "", palette = "Paired",
                xlab= paste0("PC",1," (", round(summary(pca)$importance[2,1]*100, digits = 1), "% explained var.)"), 
                ylab=paste0("PC",2," (", round(summary(pca)$importance[2,2]*100, digits = 1), "% explained var.)"))


