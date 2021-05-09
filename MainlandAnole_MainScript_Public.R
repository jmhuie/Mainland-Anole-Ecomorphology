rm(list = ls(all=TRUE))
source("MainlandAnole_Functions_Public.R")
library(tidyverse)
library(readxl)
library(phytools)
library(plyr)
library(ggpubr)
library(MASS)
library(geiger)
library(l1ou)
library(pvclust)
library(NbClust)


# Read Tree ---------------------------------------------------------------

#read in species list with ecomorph and region assignments
EcoData <- read_excel("Data/MainlandAnole_SpeciesList.xlsx")%>%
  as.data.frame
rownames(EcoData) <- EcoData$Species 
EcoData <- EcoData[sort(EcoData$Species),]

tree <- read.nexus("Trees/Poe_2017_MCC_names.txt")
#tree <- read.nexus("Trees/Poe_2017_Time_names.tre") # remove tag to use time tree

tree$tip.label[which(tree$tip.label=="cybotes")] <- "hispaniolae"
tree$tip.label[which(tree$tip.label=="tropidonotus")] <- "mccraniei"
tree$tip.label[which(tree$tip.label=="polylepis")] <- "osa"
#tree$tip.label[which(tree$tip.label=="argillaceus")] <- "ruibali" # remove tag when using time tree


#tree <- ladderize(tree, right = T)
tree <- drop.tip(tree, setdiff(tree$tip.label, unique(EcoData$Species))) #keeps tip labels that matches vector
if (is.ultrametric(tree) == FALSE) {
  tree <- force.ultrametric(tree) # forces tree to be ultrametric. Tree is ultrametic but R did not read that way due to rounding
}


# Morphology Data Cleaning ------------------------------------------------

MorphoData <- read_excel("Data/MainlandAnole_MorphoData_Public.xlsx", sheet = 1, na = "NA") # read in Morphology Data

#log the data
MorphoData[MorphoData$Species == "onca",c("Finger.II.Lam","Finger.III.Lam","Toe.II.Lam","Toe.III.Lam")] <- 1 # change onca toepad count from 0 to 1
MorphoData[,4:25] <- apply(MorphoData[,4:25], 2, log) # logs all trait values

# read in substitute tail data from Poe and Anderson 2019 and replace the NAs
TailSub <- read_excel("Data/MainlandAnole_TailSub.xlsx", sheet = 1) %>% as.data.frame
for (i in 1:nrow(TailSub)) {
  MorphoData[which(MorphoData$Species == TailSub$Species[i]),"TL"] <- log(TailSub$TL[i])
}
remove(TailSub)


# Size-Correction ---------------------------------------------------------

#size-correct the tail data
tail <- MorphoData[!is.na(MorphoData$TL),] %>%
  dplyr::group_by(Species) %>% 
  summarise_at(vars(SVL:TL), list(mean), na.rm = TRUE)  %>% # finds the species mean for each trait
  as.data.frame
rownames(tail) <- tail$Species
tail <- tail[tree$tip.label,]
tail.resid <- phyl.resid(tree, x = setNames(tail$SVL,tail$Species), Y = setNames(tail$TL,tail$Species)) # phylo size-correction

#size-correct the other morpho data
MorphoData <- MorphoData %>%
  dplyr::group_by(Species) %>% 
  summarise_at(vars(SVL:Toe.III.Lam), list(mean), na.rm = TRUE) %>%# finds the species mean for each trait
  as.data.frame
rownames(MorphoData) <- MorphoData$Species
#MorphoData<-MorphoData[tree$tip.label,]

# create new data frame with the size-corrected data
NewData <- MorphoData[tree$tip.label,]
all(match(NewData$Species, tree$tip.label))
phyl.resid <- phyl.resid(tree, x = setNames(NewData$SVL,NewData$Species), Y = NewData[,4:19]) # phylo size-correction

NewData[,3] <-tail.resid$resid
NewData[,4:19]<-phyl.resid$resid
NewData <- NewData[,c(1:3,5:9,11:13,15,17,19,21,23)]
NewData <- cbind(EcoData[tree$tip.label,-c(4,7:8)],NewData[tree$tip.label,-1])
remove(tail,tail.resid,phyl.resid,MorphoData)

# PCA ---------------------------------------------------------------------

phylopca <- phyl.pca(tree, NewData[tree$tip.label,6:18], method = "BM", mode = "cov")
plot.phyl.pca(phylopca, type = "lines")

col <- setNames(c("Blue","cyan", "Gold", "grey73", "ForestGreen", "tan4", "Red", "Purple4", "Darkorange1"),sort(unique(NewData$Ground)))
col2 <- setNames(c("Blue", "Gold", "grey73", "ForestGreen", "tan4", "Red", "Purple4", "Darkorange1"),sort(unique(NewData$Ecomorph)))

species <- rownames(phylopca$S)
#species <- rownames(phylopca$S[which(!NewData$Ground == "M" & !NewData$Ground == "U"),]) # use to plot only the ecomorph species
eco <- setNames(NewData[species,"Ground"], species)
phylopca$S <- phylopca$S[species,]
phylopca$S[,c(3)] <-  phylopca$S[,c(3)] *-1
phylopca$S[,c(1)] <-  phylopca$S[,c(1)] *-1
phylopca$S[,c(5)] <-  phylopca$S[,c(5)] *1
ggarrange(ggplot.pca(pca = phylopca, axis1 = 1, axis2 =2, species = species, 
                groups = eco, labels = F, 
                region = NewData[species,"Region"]),
          ggplot.pca(pca = phylopca, axis1 = 1, axis2 =3, species = species, 
                groups = eco, labels = F, 
                region = NewData[species,"Region"]),
          ggplot.pca(pca = phylopca, axis1 = 4, axis2 =3, species = species, 
                groups = eco, labels = F, 
                region = NewData[species,"Region"]),
          ggplot.pca(pca = phylopca, axis1 = 4, axis2 =5, species = species, 
                groups = eco, labels = F, 
                region = NewData[species,"Region"]),
          ncol =2,nrow =2)
dev.off()

#### . ####
# Randomization Tests -----------------------------------------------------

# perform randomization tests with Caribbean ecomorphs
# uses pPC scores but size-corretected traits yield same results
Ecomorph.Scores <- cbind("Ecomorph" = NewData$Ecomorph,as.data.frame(phylopca$S)) %>%
  filter(., !Ecomorph == "M" & !Ecomorph == "U")

bon <- posthoc.cross(Ecomorph.Scores, axes = 5, fun = "random.dist",
                     p.adj = "bonferroni", nsim = 1000)
bon

# perform randomization tests with Caribbean ecomorphs + ground ecomorph
# uses pPC scores but size-corretected traits yield same results
Ecomorph.Scores <- cbind("Ecomorph" = NewData$Ground,as.data.frame(phylopca$S)) %>%
  filter(., !Ecomorph == "M" & !Ecomorph == "U")
bon2 <- posthoc.cross(Ecomorph.Scores, axes = 5, fun = "random.dist", 
                      p.adj = "bonferroni", nsim = 1000)
bon2


# Cluster Analysis --------------------------------------------------------

# perform cluster analysis on Caribbean ecomorphs
EcomorphData <- NewData[which(!NewData$Ecomorph == "M" & !NewData$Ecomorph == "U"),]
clust1 <- hclust(dist((EcomorphData[,6:18]),method = "euclidean"), method = "ward.D2", members = NULL)
clust1 <- as.phylo(clust1)
plotTree(clust1,cex =0.5, offset = .3, fsize = 0.5, ftype = "i")
tiplabels(bg = col[setNames(EcomorphData$Ecomorph,EcomorphData$Species)], pch = 21)

# perform cluster analysis on Caribbean ecomorphs + ground ecomorph
GroundData <- NewData[which(!NewData$Ground == "M" & !NewData$Ground == "U"),]
clust2 <- hclust(dist((GroundData[,6:18]),method = "euclidean"), method = "ward.D2", members = NULL)
clust2 <- as.phylo(clust2)
plotTree(clust2,cex =0.5, offset = .3, fsize = 0.5, ftype = "i")
tiplabels(bg = col[setNames(GroundData$Ground,GroundData$Species)], pch = 21)



# DFA w/ Caribbean Ecomorphs ----------------------------------------------

EcomorphData <- NewData[which(!NewData$Ecomorph == "M" & !NewData$Ecomorph == "U"),]
EcomorphData <- droplevels(EcomorphData)

# MANOVA for signifiance testing
summary(manova(as.matrix(EcomorphData[,c(6:18)])~EcomorphData$Ecomorph), test = "Wilks")

# perform DFA on Caribbean ecomorph species assigned a priori
DFA <- lda(EcomorphData$Ecomorph ~ ., 
           data = EcomorphData[,c(6:18)],
           prior = c(rep(1/6, 6)), CV = FALSE) # initial DFA with subsitution
assessment <- predict(DFA)
acc <- table(EcomorphData$Ecomorph, assessment$class) # actual vs predicted without CV
acc
sum(acc[row(acc) == col(acc)]) / sum(acc) # Accuracy of the DFA with subsitution
EcomorphData$Species[which(!assessment$class == EcomorphData$Ecomorph)] # species that got missclassified

DFA_cv <- lda(EcomorphData$Ecomorph ~ ., 
           data = EcomorphData[,c(6:18)],
           prior = c(rep(1/6, 6)), CV = TRUE) # initial DFA with cross validation
acc <- table(EcomorphData$Ecomorph, DFA_cv$class) # actual vs predicted with CV
acc
sum(acc[row(acc) == col(acc)]) / sum(acc) # Accuracy of the DFA with cross validation
EcomorphData$Species[which(!DFA_cv$class == EcomorphData$Ecomorph)] # species that got missclassified


# DFA w/ Ground Ecomorph ----------------------------------------------

GroundData <- NewData[which(!NewData$Ground == "M" & !NewData$Ground == "U"),]
GroundData <- droplevels(GroundData)

# MANOVA for signifiance testing
summary(manova(as.matrix(GroundData[,c(6:18)])~GroundData$Ground), test = "Wilks")

# perform DFA on Caribbean + Ground ecomorph species assigned a priori
DFA2 <- lda(GroundData$Ground ~ ., 
           data = GroundData[,c(6:18)],
           prior = c(rep(1/7, 7)), CV = FALSE) # initial DFA with subsitution
assessment <- predict(DFA2)
acc <- table(GroundData$Ground, assessment$class) # actual vs predicted without CV
acc
sum(acc[row(acc) == col(acc)]) / sum(acc) # Accuracy of the DFA with subsitution
GroundData$Species[which(!assessment$class == GroundData$Ground)] # species that got missclassified

DFA_cv2 <- lda(GroundData$Ground ~ ., 
              data = GroundData[,c(6:18)],
              prior = c(rep(1/7, 7)), CV = TRUE) # initial DFA with cross validation
acc <- table(GroundData$Ground, DFA_cv2$class) # actual vs predicted with CV
acc
sum(acc[row(acc) == col(acc)]) / sum(acc) # Accuracy of the DFA with cross validation
GroundData$Species[which(!DFA_cv2$class == GroundData$Ground)] # species that got missclassified




#### . ####
# DFA Criteria w/ Caribbean Ecomorphs ---------------------------------

# re-do DFA
EcomorphData <- NewData[which(!NewData$Ecomorph == "M" & !NewData$Ecomorph == "U"),]
EcomorphData <- droplevels(EcomorphData)
DFA <- lda(EcomorphData$Ecomorph ~ ., 
data = EcomorphData[,c(6:18)],
prior = c(rep(1/6, 6)), CV = FALSE) # initial DFA with subsitution

# predict Draconura and unclassified Caribbean ecomorphs with DFA
predictions <- predict(DFA, as.data.frame(NewData[which(NewData$Ecomorph == "M" | NewData$Ecomorph == "U"),]))

# create table with all of DFA posterior probability results
DFA_predictions<- cbind(Species = rownames(predictions$x), 
                      Ecomorph = as.character(NewData$Ecomorph[which(NewData$Ecomorph == "M" | NewData$Ecomorph == "U")]),
                      round(predictions$posterior,3), 
                      Pred.Eco = as.character(predictions$class)) %>% data.frame()

# removes species with posterior probability >=90
DFA_criteria <- as.matrix(DFA_predictions)
for (i in 1:nrow(DFA_criteria)) {
  tmp <- max(DFA_criteria[i,4:ncol(DFA_criteria)-1])
  if (tmp < 0.90) {
    DFA_criteria[i,"Pred.Eco"] <- NA
  }
}
DFA_criteria <- data.frame(DFA_criteria[which(!is.na(DFA_criteria[,"Pred.Eco"])),])

# how many species got classified
tapply(DFA_criteria$Pred.Eco[which(DFA_criteria$Ecomorph == "M" | DFA_criteria$Ecomorph == "U")],
       DFA_criteria$Ecomorph[which(DFA_criteria$Ecomorph == "M" | DFA_criteria$Ecomorph == "U")],length)

# break down of mainland and island Draconura species that got classified with DFA
tapply(DFA_criteria$Pred.Eco[which(DFA_criteria$Ecomorph == "M")],
       DFA_criteria$Pred.Eco[which(DFA_criteria$Ecomorph == "M")],length)

# break down of previously unclassified Caribbean species that got classified with DFA
tapply(DFA_criteria$Pred.Eco[which(DFA_criteria$Ecomorph == "U")],
       DFA_criteria$Pred.Eco[which(DFA_criteria$Ecomorph == "U")],length)


# ED Criteria w/ Caribbean Ecomorphs ----------------------------------

# predict ecomorphs based on Euclidean Distance criteria
# take note of the inputs
ED_predicted <- ED.predict(scores = phylopca$S[tree$tip.label,1:5], 
                           species = tree$tip.label, 
                           groups = setNames(NewData[tree$tip.label,"Ecomorph"],tree$tip.label),
                           hard.mode = T, all.species = F)

# print result breakdown of ED criterion 1
ED_criteria1 <- ED_predicted$criteria1
tapply(ED_criteria1$Pred.Eco[which(ED_criteria1$Ecomorph=="M")],
       ED_criteria1$Pred.Eco[which(ED_criteria1$Ecomorph=="M")],length) # Draconura
tapply(ED_criteria1$Pred.Eco[which(ED_criteria1$Ecomorph=="U")],
       ED_criteria1$Pred.Eco[which(ED_criteria1$Ecomorph=="U")],length) # Unclassified Caribbean

# print result breakdown of ED criterion 2
ED_criteria2 <- ED_predicted$criteria2
tapply(ED_criteria2$Pred.Eco[which(ED_criteria2$Ecomorph=="M")],
       ED_criteria2$Pred.Eco[which(ED_criteria2$Ecomorph=="M")],length) # Draconura
tapply(ED_criteria2$Pred.Eco[which(ED_criteria2$Ecomorph=="U")],
       ED_criteria2$Pred.Eco[which(ED_criteria2$Ecomorph=="U")],length) # Unclassified Caribbean

# print result breakdown of ED criterion 3
ED_criteria3 <- ED_predicted$criteria3
tapply(ED_criteria3$Pred.Eco[which(ED_criteria3$Ecomorph=="U")],
       ED_criteria3$Pred.Eco[which(ED_criteria3$Ecomorph=="U")],length) # Draconura
tapply(ED_criteria3$Pred.Eco[which(ED_criteria3$Ecomorph=="M")],
       ED_criteria3$Pred.Eco[which(ED_criteria3$Ecomorph=="M")],length) # Unclassified Caribbean



# Final Classification w/Caribbean Ecomorphs ------------------------------

# compile the DFA and ED criteria results into final classification
# take note of the input
compile <- synth.compile(dfa = DFA_predictions, ed = ED_predicted, 
                         lower.cut = 0.9, upper.cut = 0.95, hard.mode =  T)

# how many species got classified
tapply(compile$Predicted,compile$Ecomorph,length)

# how many island and mainland Draconura species got classified
tapply(compile$Predicted[which(compile$Ecomorph == "M")],compile$Predicted[which(compile$Ecomorph == "M")],length)

# how many previosuly unclassified Caribbean species got classified
tapply(compile$Predicted[which(compile$Ecomorph == "U")],compile$Predicted[which(compile$Ecomorph == "U")],length)

# save this for later
# table with species and their update ecomorph classes 
Pred.Eco <- data.frame(Species = NewData$Species, Pred.Eco = NewData$Ecomorph)
rownames(Pred.Eco) <- Pred.Eco$Species
Pred.Eco[!is.na(match(Pred.Eco$Species,compile$Species)),2] <- compile$Predicted


# Intermediate w/Caribbean ------------------------------------------------

# determine intermediate ecomorph classifications
inter <- inter(DFA_predictions, compile)

# breakdown of intermediate ecomorph classifications for Draconura species
tapply(inter$inter.value[which(inter$inter$Ecomorph == "M"), "Match"],inter$inter.value[which(inter$inter$Ecomorph == "M"), "Match"],length)

# breakdown of intermediate ecomorph classifications for previously unclassified Caribbean 
tapply(inter$inter.value[which(inter$inter$Ecomorph == "U"), "Match"],inter$inter.value[which(inter$inter$Ecomorph == "U"), "Match"],length)


# create the prior probabilities for future simmap analysis
# uses final and intermediate classification
# intermediate species priors use the dfa posteriors
simmap.eco <- setNames(NewData$Ecomorph,NewData$Species)
simmap.eco[which(simmap.eco == "M")] <- "U"
simmap.eco[compile$Species] <- compile$Predicted
simmap.col <- setNames(c("Blue", "Gold",  "ForestGreen", "tan4","Red", "Purple4","grey73"),sort(unique(simmap.eco)))
simmap.eco<-to.matrix(simmap.eco,levels(as.factor(simmap.eco)))

for (i in 1:nrow(inter$inter)) {
  simmap.eco[inter$inter$Species[i],-ncol(simmap.eco)] <- as.numeric(DFA_predictions[inter$inter$Species[i],3:(ncol(DFA_predictions)-1)])
  simmap.eco[inter$inter$Species[i],"U"] <- 0
}
simmap.eco


#### . ####
# DFA Criteria w/ Ground Ecomorph -----------------------------------------

# re-do DFA
GroundData <- NewData[which(!NewData$Ground == "M" & !NewData$Ground == "U"),]
DFA2 <- lda(GroundData$Ground ~ ., 
            data = GroundData[,c(6:18)],
            prior = c(rep(1/7, 7)), CV = FALSE) # initial DFA with subsitution

# predict Draconura and unclassified Caribbean ecomorphs with DFA
predictions <- predict(DFA2, as.data.frame(NewData[which(NewData$Ground == "M" | NewData$Ground == "U"),]))

# create table with all of DFA posterior probability results
DFA_predictions2 <- cbind(Species = rownames(predictions$x), 
                      Ecomorph = as.character(NewData$Ground[which(NewData$Ground == "M" | NewData$Ground == "U")]),
                      round(predictions$posterior,3), 
                      Pred.Eco = as.character(predictions$class))

#removes species with post prob <90
DFA_criteria2 <- DFA_predictions2
for (i in 1:nrow(DFA_criteria2)) {
  tmp <- max(DFA_criteria2[i,4:ncol(DFA_criteria2)-1])
  if (tmp < .90) {
    DFA_criteria2[i,"Pred.Eco"] <- NA
  }
}
DFA_criteria2 <- as.data.frame(DFA_criteria2[which(!is.na(DFA_criteria2[,"Pred.Eco"])),])

# how many species got classified
tapply(DFA_criteria2$Pred.Eco[which(DFA_criteria2$Ecomorph == "M" | DFA_criteria2$Ecomorph == "U")],
       DFA_criteria2$Ecomorph[which(DFA_criteria2$Ecomorph == "M" | DFA_criteria2$Ecomorph == "U")],length)

# break down of mainland and island Draconura species that got classified with DFA
tapply(DFA_criteria2$Pred.Eco[which(DFA_criteria2$Ecomorph == "M")],
       DFA_criteria2$Pred.Eco[which(DFA_criteria2$Ecomorph == "M")],length)

# break down of previously unclassified Caribbean species that got classified with DFA
tapply(DFA_criteria2$Pred.Eco[which(DFA_criteria2$Ecomorph == "U")],
       DFA_criteria2$Pred.Eco[which(DFA_criteria2$Ecomorph == "U")],length)


# ED Criteria w/ Ground ---------------------------------------------------

# predict ecomorphs based on certain criteria

# predict ecomorphs based on Euclidean Distance criteria
# take note of the inputs
ED_predicted2 <- ED.predict(scores = phylopca$S[tree$tip.label,1:5], 
                           species = tree$tip.label, 
                           groups = setNames(NewData[tree$tip.label,"Ground"],tree$tip.label),
                           hard.mode = T, all.species = F)

# print result breakdown of ED criterion 1
ED_criteria1 <- ED_predicted2$criteria1
tapply(ED_criteria1$Pred.Eco[which(ED_criteria1$Ecomorph=="M")],
       ED_criteria1$Pred.Eco[which(ED_criteria1$Ecomorph=="M")],length) # Draconura
tapply(ED_criteria1$Pred.Eco[which(ED_criteria1$Ecomorph=="U")],
       ED_criteria1$Pred.Eco[which(ED_criteria1$Ecomorph=="U")],length) # Unclassified Caribbean

# print result breakdown of ED criterion 2
ED_criteria2 <- ED_predicted2$criteria2
tapply(ED_criteria2$Pred.Eco[which(ED_criteria2$Ecomorph=="M")],
       ED_criteria2$Pred.Eco[which(ED_criteria2$Ecomorph=="M")],length) # Draconura
tapply(ED_criteria2$Pred.Eco[which(ED_criteria2$Ecomorph=="U")],
       ED_criteria2$Pred.Eco[which(ED_criteria2$Ecomorph=="U")],length) # Unclassified Caribbean

# print result breakdown of ED criterion 3
ED_criteria3 <- ED_predicted2$criteria3
tapply(ED_criteria3$Pred.Eco[which(ED_criteria3$Ecomorph=="U")],
       ED_criteria3$Pred.Eco[which(ED_criteria3$Ecomorph=="U")],length) # Draconura
tapply(ED_criteria3$Pred.Eco[which(ED_criteria3$Ecomorph=="M")],
       ED_criteria3$Pred.Eco[which(ED_criteria3$Ecomorph=="M")],length) # Unclassified Caribbean



# Final Classifications w/ Ground ---------------------------------------

# compile the DFA and ED criteria results into final classification
# take note of the input
compile2 <- synth.compile(dfa = DFA_predictions2, ed = ED_predicted2, 
                         lower.cut = 0.9, upper.cut = 0.95, hard.mode =  T)

# how many species got classified
tapply(compile2$Predicted,compile2$Ecomorph,length)

# how many island and mainland Draconura species got classified
tapply(compile2$Predicted[which(compile2$Ecomorph == "M")],compile2$Predicted[which(compile2$Ecomorph == "M")],length)

# how many previosuly unclassified Caribbean species got classified
tapply(compile2$Predicted[which(compile2$Ecomorph == "U")],compile2$Predicted[which(compile2$Ecomorph == "U")],length)

# figure out which species that used to be classified w/out the ground ecomorph
# but became unclassified with it 
compile$Species[which(is.na(match(compile$Species, compile2$Species)))]

# reconcile species that were not classified with ground ecomorph 
# but were previously classified without it
# not really a function just using it to not create too many variables
comp.diff<-function(){
  diff <- compile$Species[which(is.na(match(compile$Species, compile2$Species)))]
  diff <- diff[which(is.na(match(diff,EcomorphData$Species)))]
  diff # vector with the species that became unclassified
  
  old <- compile[!is.na(match(compile$Species,diff)),]; old # old classifications
  rownames(old) <- old$Species
  new <- as.data.frame(DFA_predictions[diff,]);new # new classifications
  just <- cbind(old[,c("Species","Predicted")],New.LDA = new[,"Pred.Eco"], C1.Predicted = ED_criteria1[diff,"Pred.Eco"])
  need.c1 <- ED.predict(scores = phylopca$S[,1:5], species = tree$tip.label, 
                groups = NewData$Ground,
                hard.mode = F, all.species = F)
  just[which(is.na(just$C1.Predicted)),"C1.Predicted"] <- need.c1$criteria1[just[which(is.na(just$C1.Predicted)),"Species"],"Pred.Eco"]
  just
  synth <- old[which(just$C1.Predicted == just$New.LDA & just$C1.Predicted == just$Predicted),]
  for (i in 1:nrow(synth)) {
    synth[i,"LDA"] <- new[synth$Species[i],synth$Predicted[i]]
  }
  synth <- synth[which(synth$LDA >= 0.75),]
  trans <- just$Species[which(just$New.LDA == "G" & just$C1.Predicted == just$New.LDA)]
  total <- synth.compile(dfa = DFA_predictions, ed = ED_predicted, hard.mode = F)
  trans <- total[!is.na(match(total$Species,trans)),]
  trans <- trans[which(trans$LDA >= 0.75),]
  rownames(trans) <- trans$Species
  
  new.compile2 <- rbind(compile2,synth,trans)
  new.compile2 <- new.compile2[sort(new.compile2$Species),]
  
  return(new.compile2)
}
new.compile2 <- comp.diff()

compile$Species[which(is.na(match(compile$Species, new.compile2$Species)))] # Species that still did not get classified

# how many species got classified
tapply(new.compile2$Predicted,new.compile2$Ecomorph,length)

# how many island and mainland Draconura species got classified
tapply(new.compile2$Predicted[which(new.compile2$Ecomorph == "M")],new.compile2$Predicted[which(new.compile2$Ecomorph == "M")],length)

# how many previosuly unclassified Caribbean species got classified
tapply(new.compile2$Predicted[which(new.compile2$Ecomorph == "U")],new.compile2$Predicted[which(new.compile2$Ecomorph == "U")],length)


# save this for later
# table with species and their update ecomorph classes 
Pred.Eco2 <- data.frame(Species = NewData$Species, Pred.Eco = NewData$Ground)
rownames(Pred.Eco2) <- Pred.Eco2$Species
Pred.Eco2[!is.na(match(Pred.Eco2$Species,new.compile2$Species)),2] <- new.compile2$Predicted

# Intermediate w/Ground --------------------------------------------------

# determine intermediate ecomorph classifications
inter2 <- inter(DFA_predictions2, new.compile2)

# breakdown of intermediate ecomorph classifications for Draconura species
tapply(inter2$inter.value[which(inter2$inter$Ecomorph == "M"), "Match"],inter2$inter.value[which(inter2$inter$Ecomorph == "M"), "Match"],length)

# breakdown of intermediate ecomorph classifications for previously unclassified Caribbean 
tapply(inter2$inter.value[which(inter2$inter$Ecomorph == "U"), "Match"],inter2$inter.value[which(inter2$inter$Ecomorph == "U"), "Match"],length)


# create the prior probabilities for future simmap analysis
# uses final and intermediate classification
# intermediate species priors use the dfa posteriors
simmap.eco2 <- setNames(NewData$Ground,NewData$Species)
simmap.eco2[which(simmap.eco2 == "M")] <- "U"
simmap.eco2[new.compile2$Species] <- new.compile2$Predicted
simmap.col2 <- setNames(c("Blue", "cyan","Gold",  "ForestGreen", "tan4","Red", "Purple4","grey73"),sort(unique(simmap.eco2)))
simmap.eco2<-to.matrix(simmap.eco2,levels(as.factor(simmap.eco2)))

for (i in 1:nrow(inter2$inter)) {
  #split <- str_split(intermediate[i,"Match"], "/")
  #simmap.eco[intermediate$Species[i],c(split)[[1]]] <- 0.5
  simmap.eco2[inter2$inter$Species[i],-ncol(simmap.eco2)] <- as.numeric(DFA_predictions2[inter2$inter$Species[i],3:(ncol(DFA_predictions2)-1)])
  simmap.eco2[inter2$inter$Species[i],"U"] <- 0
}
simmap.eco2



# . -----------------------------------------------------------------------
# SIMMAP Analyses w/Caribbean  --------------------------------------------

# ER -273.6212 df 1
#> SYM -246.9038 df 21
#> ER vs SYM LR 64.1824 df 20 
#-2*(-273.6212 - (-246.9038))
#quantile(rchisq(10000, 20),.95) #SYM IS BETTER!!!! YAY
# ARD -236.2189 df 42
# SYM vs ARD LR 17.272 df 21
#-2*(-246.9038 - (-236.2189))
#quantile(rchisq(10000, 21),.95) #ARD IS NOT BETTER!!!! GO FOR SYM

#simmap<-make.simmap(tree,simmap.eco,model = "SYM",nsim = 1000)

# Post Burn In Sim
hundotrees <- read.tree("Trees/Poe_2017_MCC_postburnin_names_100.txt")
for ( i in 1:100) {
  hundotrees[[i]]$tip.label[which(hundotrees[[i]]$tip.label=="cybotes")] <- "hispaniolae"
  hundotrees[[i]]$tip.label[which(hundotrees[[i]]$tip.label=="tropidonotus")] <- "mccraniei"
  hundotrees[[i]]$tip.label[which(hundotrees[[i]]$tip.label=="polylepis")] <- "osa"
}
hundotrees<-lapply(hundotrees,keep.tip,tree$tip.label)
class(hundotrees) <- "multiPhylo"
Mainland.tree <- drop.tip(tree,setdiff(tree$tip.label,NewData$Species[which(NewData$Region == "Mainland")]))
hundotrees <- lapply(hundotrees, ladderize, right = FALSE)
class(hundotrees) <- "multiPhylo"

#define which set of simmap species and col to use
simeco <- simmap.eco2
simeco.col <- simmap.col2

post.sim<-make.simmap(hundotrees,simeco,model = "SYM",nsim = 100)
write.simmap(post.sim, "Outputs/TooBig/SIMMAP_SYM_Ground_postburnin_NEW.txt", append = TRUE, version = 1.0)
#write.simmap(post.sim1, "Outputs/SIMMAP_SYM_Ground_MCC.txt", append = TRUE, version = 1.0)

post.sim <- read.simmap("Outputs/TooBig/SIMMAP_SYM_Caribbean_postburnin.txt", format="phylip",version =1.0)
#class(post.sim) <- "multiPhylo"

simmap.tree<-lapply(post.sim,drop.tip.simmap,NewData$Species[which(!NewData$Region == "Mainland")])
class(simmap.tree) <- "multiPhylo"
pd<-describe.simmap(simmap.tree, plot = FALSE, ref.tree = Mainland.tree)
#write_rds(pd, "Outputs/TooBig/SIMMAP_SYM_Ground_post_PD_NEW.rds")
pd <- read_rds("Outputs/TooBig/SIMMAP_SYM_Ground_post_PD_NEW.rds")


transitions <- pd$count
transition.type <- (",G$")
tapply(apply(transitions[,grep(transition.type, colnames(transitions))],1,sum),
       apply(transitions[,grep(transition.type, colnames(transitions))],1,sum),length)

#plot(simmap.tree[[1]], col = simeco.col, fsize = 0.5, offset = 0.3)
plotTree(Mainland.tree, fsize = 0.5, lwd = 2, offset = 0.3)
nodelabels(pie=pd$ace,piecol=simeco.col[colnames(pd$ace)],cex=.35)
tiplabels(pie=simeco[Mainland.tree$tip.label,],piecol=simeco.col,cex=0.2)



# Wheatsheaf's Index ------------------------------------------------------

dat <- cbind(NewData[,c("Species", "Ground")],as.data.frame(phylopca$S)[,1:5])
dat[which(dat[,"Ground"] == "G"),2] <- 1
dat[new.compile2[which(new.compile2$Predicted == "G"),"Species"],2] <- 1
dat[which(!dat[,"Ground"] == 1),2] <- 0
colnames(dat) <- c("species","focal","PC1","PC2","PC3","PC4","PC5")

conv <- windex(dat[tree$tip.label,], tree, c(3:7), SE = FALSE);conv
cov.text <- test.windex(dat[tree$tip.label,], tree, c(3:7),reps = 100, SE = FALSE);cov.text


dat <- cbind(NewData[,c("Species", "Ground")],as.data.frame(phylopca$S)[,1:5])
dat[compile[which(compile$Predicted == "Tw"),"Species"],2] <- 1
dat[which(!dat[,"Ground"] == 1),2] <- 0
dat <- dat[Mainland.tree$tip.label,]
colnames(dat) <- c("species","focal","PC1","PC2","PC3","PC4","PC5")

conv <- windex(dat[Mainland.tree$tip.label,], Mainland.tree, c(3:7), SE = FALSE);conv
cov.text <- test.windex(dat[Mainland.tree$tip.label,], Mainland.tree, c(3:7),reps = 1000, SE = FALSE);cov.text


# C Tests -----------------------------------------------------------------

traits <- as.matrix(phylopca$S[which(Pred.Eco$Region == "Mainland"),1:5])
convtips <- Pred.Eco2[which(Pred.Eco$Region == "Mainland" & Pred.Eco2$Pred.Eco == "GB"),"Species"]
convSig(Mainland.tree,traits,convtips,nsim =1000)

# Species for Species Matching --------------------------------------------

FiveIsland <- NewData %>%
  filter(Country == "DR" | Country == "Haiti" | Country == "Puerto Rico" | 
           Country == "Cuba" | Country == "Jamaica"| Region ==  "Mainland")
FiveIsland <- cbind(FiveIsland[,1:5], phylopca$S[rownames(FiveIsland),])
FiveIsland[which(FiveIsland$Country == "DR" | FiveIsland$Country == "Haiti"),"Country"] <- "Hispaniola"
FiveIsland[which(FiveIsland$Region == "Caribbean"),"Region"] <- FiveIsland[which(FiveIsland$Region == "Caribbean"),"Country"]

euc <- as.matrix(dist(FiveIsland[,6:10], method = "euclidean"))
euc[which(euc == 0)] <- NA
MainNND <- data.frame(FiveIsland[,c(1,4)], 
                      Cuba.NND = NA,
                      Jam.NND = NA,
                      Hisp.NND = NA,
                      PR.NND = NA,
                      Main.NND = NA)
for (i in 1:nrow(MainNND)) {
  MainNND[i,"Cuba.NND"] <- min(euc[i,which(FiveIsland$Region == "Cuba")],na.rm = T)
  MainNND[i,"Jam.NND"] <- min(euc[i,which(FiveIsland$Region == "Jamaica")],na.rm = T)
  MainNND[i,"Hisp.NND"] <- min(euc[i,which(FiveIsland$Region == "Hispaniola")],na.rm = T)
  MainNND[i,"PR.NND"] <- min(euc[i,which(FiveIsland$Region == "Puerto Rico")],na.rm = T)
  MainNND[i,"Main.NND"] <- min(euc[i,which(FiveIsland$Region == "Mainland")],na.rm = T)
}
Cuba <- mean(unlist(MainNND[which(MainNND$Region == "Cuba"),c(4,5,6,7)]))
#Cuba <- mean(unlist(MainNND[which(MainNND$Region == "Cuba"),c(7)]))
Jam <- mean(unlist(MainNND[which(MainNND$Region == "Jamaica"),c(3,5,6,7)]))
#Jam <- mean(unlist(MainNND[which(MainNND$Region == "Jamaica"),c(7)]))
Hisp <- mean(unlist(MainNND[which(MainNND$Region == "Hispaniola"),c(3,4,6,7)]))
#Hisp <- mean(unlist(MainNND[which(MainNND$Region == "Hispaniola"),c(7)]))
PR <- mean(unlist(MainNND[which(MainNND$Region == "Puerto Rico"),c(3,4,5,7)]))
#PR <- mean(unlist(MainNND[which(MainNND$Region == "Puerto Rico"),c(7)]))
Main <- mean(unlist(MainNND[which(MainNND$Region == "Mainland"),c(3,4,5,6)]))
#Main <- apply((MainNND[which(MainNND$Region == "Mainland"),c(3,4,5,6)]),2,mean)

totalmean <- mean(c(Cuba,Jam,Hisp,PR));totalmean
Carmean <- mean(c(Cuba,Jam,Hisp,PR))
Mainmean <- mean(Main)
#totalmean <- mean(c(Cuba,Jam,Hisp,PR,apply(MainNND[,3:6],2,mean)))
#totalmean <- Main



fit1 <- fitContinuous(tree,setNames(phylopca$S[,1],rownames(phylopca$S))[tree$tip.label],model = "BM")
fit2 <- fitContinuous(tree,setNames(phylopca$S[,2],rownames(phylopca$S))[tree$tip.label],model = "BM")
fit3 <- fitContinuous(tree,setNames(phylopca$S[,3],rownames(phylopca$S))[tree$tip.label],model = "BM")
fit4 <- fitContinuous(tree,setNames(phylopca$S[,4],rownames(phylopca$S))[tree$tip.label],model = "BM")
fit5 <- fitContinuous(tree,setNames(phylopca$S[,5],rownames(phylopca$S))[tree$tip.label],model = "BM")
totalsim <-matrix(NA,1000,1)
Carsim <-matrix(NA,1000,1)
Mainsim <-matrix(NA,1000,1)
for(y in 1:nrow(totalsim)) {
  simNND <- data.frame(FiveIsland[,c(1,4)], 
                       Cuba.NND = NA,
                       Jam.NND = NA,
                       Hisp.NND = NA,
                       PR.NND = NA,
                       Main.NND = NA)
  simdata <- cbind(fastBM(tree, sig2 = fit1$opt$sigsq, nsim = 1),
                   fastBM(tree, sig2 = fit2$opt$sigsq, nsim = 1),
                   fastBM(tree, sig2 = fit3$opt$sigsq, nsim = 1),
                   fastBM(tree, sig2 = fit4$opt$sigsq, nsim = 1),
                   fastBM(tree, sig2 = fit5$opt$sigsq, nsim = 1)) # simulate data under BM
  sim.euc <- as.matrix(dist(simdata[FiveIsland$Species,], method = "euclidean"))
  sim.euc[which(sim.euc == 0)] <- NA
  for (i in 1:nrow(sim.euc)) {
    simNND[i,"Cuba.NND"] <- min(sim.euc[i,which(FiveIsland$Region == "Cuba")], na.rm = T)
    simNND[i,"Jam.NND"] <- min(sim.euc[i,which(FiveIsland$Region == "Jamaica")], na.rm = T)
    simNND[i,"Hisp.NND"] <- min(sim.euc[i,which(FiveIsland$Region == "Hispaniola")], na.rm = T)
    simNND[i,"PR.NND"] <- min(sim.euc[i,which(FiveIsland$Region == "Puerto Rico")], na.rm = T)
    simNND[i,"Main.NND"] <- min(sim.euc[i,which(FiveIsland$Region ==  "Mainland")], na.rm = T)
  }
  simCuba <- mean(unlist(simNND[which(simNND$Region == "Cuba"), c(4,5,6)]))
  #simCuba <- mean(unlist(simNND[which(simNND$Region == "Cuba"), c(7)]))
  simJam <- mean(unlist(simNND[which(simNND$Region == "Jamaica"), c(3,5,6)]))
  #simJam <- mean(unlist(simNND[which(simNND$Region == "Jamaica"), c(7)]))
  simHisp <- mean(unlist(simNND[which(simNND$Region == "Hispaniola"), c(3,4,6)]))
  #simHisp <- mean(unlist(simNND[which(simNND$Region == "Hispaniola"), c(7)]))
  simPR <- mean(unlist(simNND[which(simNND$Region == "Puerto Rico"), c(3,4,5)]))
  #simPR <- mean(unlist(simNND[which(simNND$Region == "Puerto Rico"), c(7)]))
  simMain <- mean(unlist(simNND[which(simNND$Region == "Mainland"), c(3,4,5,6)]))
  #simMain <- apply((simNND[which(simNND$Region == "Mainland"), c(3,4,5,6)]),2,mean)
  
  
  totalsim[y] <- mean(c(simCuba,simJam,simHisp,simPR))
  Carsim[y] <- mean(c(simCuba,simJam,simHisp,simPR))
  Mainsim[y] <- mean(simMain)
}
quantile(c(totalsim),.05)
totalmean
which(!is.na(match(sort(c(totalmean,totalsim)),totalmean)))/1000

hist(c(totalsim,totalmean))
abline(v = totalmean, col = "red", lwd = 2)

quantile(c(Carsim),.05)
Carmean
which(!is.na(match(sort(c(Carmean,Carsim)),Carmean)))/1000

quantile(c(Mainsim),.05)
Mainmean
which(!is.na(match(sort(c(Mainmean,Mainsim)),Mainmean)))/1000

#RESULTS - it really matters whether you compare a species to the same species or a random one.
# Former is super sig but the later is really not sig
# I think former is correct 



# L1OU --------------------------------------------------------------------

lizard <- adjust_data(tree,phylopca$S[tree$tip.label,1:5])

#fit_ind <- estimate_shift_configuration(lizard$tree,lizard$Y, nCores = 6, criterion = "AIC")
#saveRDS(fit_ind, "Outputs/shift_config_AIC_MCC.rds")
fit_ind <- readRDS("Outputs/shift_config_AIC_MCC.rds")

#fit_conv <- estimate_convergent_regimes(fit_ind, criterion = "AIC", nCores = 4)
#saveRDS(fit_conv, "Outputs/convergent_regime_AIC_MCC.rds")
fit_conv <- readRDS("Outputs/convergent_regime_AIC_MCC.rds")

#fit_ind_AIC_bootstrap <- l1ou_bootstrap_support(fit_ind,nItrs = 100, multicore = T, nCores=4)
#saveRDS(fit_ind_AIC_bootstrap, "Outputs/shift_bootstrap_AIC_MCC.rds")
bootstrap <- readRDS("Outputs/shift_bootstrap_AIC_MCC.rds")

BS_constrained_pBIC_all_data <- round(bootstrap[[1]] * 100,digits=1)
BS_constrained_pBIC_all_data <- ifelse(BS_constrained_pBIC_all_data>50,  paste0(BS_constrained_pBIC_all_data,""), NA)
BS_constrained_pBIC_all_data[setdiff(1:nrow(fit_conv$tree$edge),as.numeric(fit_conv$shift.configuration))] <- NA

ann <- setNames(bootstrap$detection.rate[fit_ind$shift.configuration], fit_ind$shift.configuration)


pal <-c(fit_ind$shift.configuration,'darkgrey')
pal[!is.na(match(pal,396))] <- "#009EFF"
pal[!is.na(match(pal,268))] <- "Red"
pal[!is.na(match(pal,c(182, 355)))] <- "#FFBB05"
pal[!is.na(match(pal,c(402, 334, 399)))] <- "Blue"
pal[!is.na(match(pal,c(161, 117)))] <- "#3D6CD7" #col_vector[1]
pal[!is.na(match(pal,c(123, 176)))] <- "cyan"
pal[!is.na(match(pal,c(94, 345)))] <- "#31FEAD" #col_vector[2]
pal[!is.na(match(pal,c(192)))] <- "#FF96B4" #col_vector[3]
pal[!is.na(match(pal,c(341, 383)))] <- "#C47DF9"
pal[!is.na(match(pal,c(405, 305, 346, 203)))] <-  "#560097"
pal[!is.na(match(pal,c(319, 260, 213)))] <- "Gold"
pal[!is.na(match(pal,c(206)))] <-"#CFFD56" #col_vector[4]
pal[!is.na(match(pal,c(384)))] <- "#E9DF0E"
pal[!is.na(match(pal,c(320, 153)))] <- "#9D19FF"
pal[!is.na(match(pal,c(286, 389)))] <- "forestgreen"
pal[!is.na(match(pal,c(179, 14)))] <-"#FF7A24" #col_vector[5]



plot(fit_conv, edge.ann.cex = .5, edge.label.ann=F,
     cex=0.25,label.offset=0.01, 
     show.tip.label=T, plot.bar=F, edge.shift.ann=F,bar.axis=F, 
     edge.label.adj=1.5, asterisk=F, edge.width=1.75, palette = pal)
group <- setNames(Pred.Eco2$Pred.Eco,Pred.Eco2$Species)
group[which(group == "U")] <- "M"
tiplabels(pch=21,frame="none",bg=col[group[fit_conv$tree$tip.label]],cex=.4)

Z = l1ou:::generate_design_matrix(fit_conv$tree, type = "apprX")
for (idx in as.numeric(names(which(ann >=0.95)))) {
  pos = max(Z[, idx])
  edgelabels(edge = idx, pch = 21, cex = 1.2, col = "black", bg = "black",
             adj = c(0.5, 0.8), frame = "none", date = pos)
}
for (idx in as.numeric(names(which(ann < 0.95 & ann >=0.75)))) {
  pos = max(Z[, idx])
  edgelabels(edge = idx, pch = 21, cex = 1.2, col = "black", bg = "darkgrey",
             adj = c(0.5, 0.8), frame = "none", date = pos)
}
for (idx in as.numeric(names(which(ann < 0.75)))) {
  pos = max(Z[, idx])
  edgelabels(edge = idx, pch = 21, cex = 1.2, col = "black", bg = "white",
             adj = c(0.5, 0.8), frame = "none", date = pos)
}


#regime 1 - luteogularis
#regime 10 - trunk anoles
#regime 11 - olssoni & auratus
#regime 12 - barahone & nobeli & vermiculatus
#regime 13 - petersii & biporcatus
#regime 14 - chrysolepis & gracilipes
#regime 15 - second major Draconura shift and etheridgei
#regime 16 - first major Draconura shift
#regime 2 - darlingtoni & big twig
#regime 3 - 4 twig regimes
#regime 4 - 3 GB regimes
#regime 5 - pinchoti
#regime 6 - other GB with koopmani
#regime 7 - mainland twig & angusticeps
#regime 8 - two TC regimes
#regime 9 - onca & lionotus/oxylophus


fit_ind <- readRDS("Outputs/shift_config_AIC_MRCT.rds")

#fit_conv <- estimate_convergent_regimes(fit_ind, criterion = "AIC", nCores = 4)
#saveRDS(fit_conv, "Outputs/convergent_regime_AIC_MRCT.rds")
fit_conv <- readRDS("Outputs/convergent_regime_AIC_MRCT.rds")

#fit_ind_AIC_bootstrap <- l1ou_bootstrap_support(fit_ind,nItrs = 100, multicore = T, nCores=4)
#saveRDS(fit_ind_AIC_bootstrap, "Outputs/shift_bootstrap_AIC_MRCT.rds")
bootstrap <- readRDS("Outputs/shift_bootstrap_AIC_MRCT.rds")

BS_constrained_pBIC_all_data <- round(bootstrap[[1]] * 100,digits=1)
BS_constrained_pBIC_all_data <- ifelse(BS_constrained_pBIC_all_data>50,  paste0(BS_constrained_pBIC_all_data,""), NA)
BS_constrained_pBIC_all_data[setdiff(1:nrow(fit_conv$tree$edge),as.numeric(fit_conv$shift.configuration))] <- NA

ann <- setNames(bootstrap$detection.rate[fit_ind$shift.configuration], fit_ind$shift.configuration)


pal <-c(fit_ind$shift.configuration,'darkgrey')
pal[!is.na(match(pal,c(353,75,116)))] <- "ForestGreen"
pal[!is.na(match(pal,c(299,189,223,332)))] <- "Gold"
pal[!is.na(match(pal,c(274)))] <- 3
pal[!is.na(match(pal,c(11)))] <- "#FF7A24"
pal[!is.na(match(pal,c(76)))] <- "#31FEAD" #col_vector[1]
pal[!is.na(match(pal,c(183,369,279,330)))] <- "Purple4"
pal[!is.na(match(pal,c(354)))] <- "#E9DF0E" #col_vector[2]
pal[!is.na(match(pal,c(356)))] <- "Blue" #col_vector[3]
pal[!is.na(match(pal,c(139)))] <- "pink"
pal[!is.na(match(pal,c(171,103)))] <-  "cyan"
pal[!is.na(match(pal,c(321, 341)))] <- "#C47DF9"
pal[!is.na(match(pal,c(311,241,283)))] <-"Red" #col_vector[4]
pal[!is.na(match(pal,c(300,135)))] <- "#9D19FF"


plot(fit_conv, edge.ann.cex = .5, edge.label.ann=F,
     cex=0.25,label.offset=0.01, 
     show.tip.label=T, plot.bar=F, edge.shift.ann=F,bar.axis=F, 
     edge.label.adj=1.5, asterisk=F, edge.width=1.75, palette = pal)
group <- setNames(Pred.Eco2$Pred.Eco,Pred.Eco2$Species)
group[which(group == "U")] <- "M"
tiplabels(pch=21,frame="none",bg=col[group[fit_conv$tree$tip.label]],cex=.4)

Z = l1ou:::generate_design_matrix(fit_conv$tree, type = "apprX")
for (idx in as.numeric(names(which(ann >=0.95)))) {
  pos = max(Z[, idx])
  edgelabels(edge = idx, pch = 21, cex = 1.2, col = "black", bg = "black",
             adj = c(0.5, 0.8), frame = "none", date = pos)
}
for (idx in as.numeric(names(which(ann < 0.95 & ann >=0.75)))) {
  pos = max(Z[, idx])
  edgelabels(edge = idx, pch = 21, cex = 1.2, col = "black", bg = "darkgrey",
             adj = c(0.5, 0.8), frame = "none", date = pos)
}
for (idx in as.numeric(names(which(ann < 0.75)))) {
  pos = max(Z[, idx])
  edgelabels(edge = idx, pch = 21, cex = 1.2, col = "black", bg = "white",
             adj = c(0.5, 0.8), frame = "none", date = pos)
}



#regime 1 - omiltemanus & ortonii & coelestinus
#regime 10 - 4 GB regimes
#regime 11 - 1 TC regime
#regime 12 - lionotus
#regime 13 - major Draconura radiation
#regime 2 - 4 twig regimes
#regime 3 - 1 GB regime
#regime 4 - nobeli
#regime 5 - salvini
#regime 6 - chrysolepis and gracilipes
#regime 7 - big twig and darlingtoni
#regime 8 - two trunk + eugenegrahami
#regime 9 - angusticeps and mainland twig
#
optima <- data.frame(fit_conv$optima)
colnames(optima) <- c("PC1","PC2","PC3","PC4","PC5")
optima<-distinct(optima, PC1,.keep_all = T)
optima$col <- NewData[rownames(optima),"Ground"]
optima["vermiculatus","col"] <- "CG"
optima["altavelensis","col"] <- "Tr"
optima["barbatus","col"] <- "Tw"
optima2 <- optima[which(!rownames(optima) == "luteogularis" & !rownames(optima) == "pinchoti") ,]
#optima2 <- optima[which(!rownames(optima) == "salvini" & !rownames(optima) == "noblei") ,]
#optima2["auratus","col"] <- "G"
#optima2["omiltemanus","col"] <- "TC"




ggplot(as.data.frame(phylopca$S), aes(x = PC1,y=PC2)) +
  geom_point(aes(fill = NewData$Ground), col = "black",shape =21, size = 2) +
  scale_fill_manual(values = col) +
  geom_point(data = as.data.frame(optima2), aes(x = PC1, y = PC2, fill = col), size = 4, shape = 21)+
  theme_classic()



####PLOT SIMMAP####
simmap.eco.max <- setNames(c(1:nrow(simmap.eco)), rownames(simmap.eco))
for (i in 1:length(simmap.eco.max)) {
  simmap.eco.max[i] <- colnames(simmap.eco)[which(simmap.eco[i,] ==  max(simmap.eco[i,]))]
}

CG<- names(which(simmap.eco.max[Mainland.tree$tip.label] == "CG"))
GB<- names(which(simmap.eco.max[Mainland.tree$tip.label] == "GB"))
TR<- names(which(simmap.eco.max[Mainland.tree$tip.label] == "Tr"))
TC<- names(which(simmap.eco.max[Mainland.tree$tip.label] == "TC"))
TG<- names(which(simmap.eco.max[Mainland.tree$tip.label] == "TG"))
TW<- names(which(simmap.eco.max[Mainland.tree$tip.label] == "Tw"))

coltree <- Mainland.tree2
coltree$maps <- NULL
tt<-paintBranches(coltree,edge=sapply(CG,match,coltree$tip.label),
                  state="CG",anc.state="U")
tt<-paintBranches(tt,edge=sapply(GB,match,coltree$tip.label),
                  state="GB")
tt<-paintBranches(tt,edge=sapply(TR,match,coltree$tip.label),
                  state="Tr")
tt<-paintBranches(tt,edge=sapply(TC,match,coltree$tip.label),
                  state="TC")
tt<-paintBranches(tt,edge=sapply(TG,match,coltree$tip.label),
                  state="TG")
tt<-paintBranches(tt,edge=sapply(TW,match,coltree$tip.label),
                  state="Tw")


simmap.eco2.max <- setNames(c(1:nrow(simmap.eco2)), rownames(simmap.eco2))
for (i in 1:length(simmap.eco2.max)) {
  simmap.eco2.max[i] <- colnames(simmap.eco2)[which(simmap.eco2[i,] ==  max(simmap.eco2[i,]))]
}

CG2<- names(which(simmap.eco2.max[Mainland.tree$tip.label] == "CG"))
GB2<- names(which(simmap.eco2.max[Mainland.tree$tip.label] == "GB"))
G2<- names(which(simmap.eco2.max[Mainland.tree$tip.label] == "G"))
TR2<- names(which(simmap.eco2.max[Mainland.tree$tip.label] == "Tr"))
TC2<- names(which(simmap.eco2.max[Mainland.tree$tip.label] == "TC"))
TG2<- names(which(simmap.eco2.max[Mainland.tree$tip.label] == "TG"))
TW2<- names(which(simmap.eco2.max[Mainland.tree$tip.label] == "Tw"))

coltree2 <- Mainland.tree2
coltree2$maps <- NULL
tt2<-paintBranches(coltree,edge=sapply(CG2,match,coltree2$tip.label),
                   state="CG",anc.state="U")
tt2<-paintBranches(tt2,edge=sapply(GB2,match,coltree2$tip.label),
                   state="GB")
tt2<-paintBranches(tt2,edge=sapply(TR2,match,coltree2$tip.label),
                   state="Tr")
tt2<-paintBranches(tt2,edge=sapply(TC2,match,coltree2$tip.label),
                   state="TC")
tt2<-paintBranches(tt2,edge=sapply(TG2,match,coltree2$tip.label),
                   state="TG")
tt2<-paintBranches(tt2,edge=sapply(TW2,match,coltree2$tip.label),
                   state="Tw")
tt2<-paintBranches(tt2,edge=sapply(G2,match,coltree2$tip.label),
                   state="G")
plot(tt2,colors=simmap.col2,lwd=3,split.vertical=TRUE,ftype="i",fsize = 0.6, offset = 0.3)
nodelabels(pie=pd3$ace,piecol=simeco.col2,cex=.6)
tiplabels(pie=simmap.eco2[Mainland.tree$tip.label,],piecol=simmap.col2,cex=0.2)


is_tip <- Mainland.tree2$edge[,2] <= length(Mainland.tree2$tip.label)
ordered_tips <- Mainland.tree2$edge[is_tip, 2]
Mainland.tree2$tip.label[ordered_tips]

pdf(width = 7.5, height = 7.8)
layout(matrix(1:3,1,3),widths=c(0.47,0.1,0.47))
plot(tt,colors=simmap.col,lwd=3,split.vertical=TRUE,ftype = "off")
nodelabels(pie=pd$ace,piecol=simmap.col,cex=.6)
tiplabels(pie=simmap.eco[Mainland.tree2$tip.label,],piecol=simmap.col,cex=0.4)
#legend("topleft", c("Crown-Giant","Grass-Bush","Ground","Trunk","Trunk-Crown","Trunk-Ground","Twig","Unclassified"), pch = 21, pt.bg = simmap.col2[c(1,3,2,6,4,5,7,8)], 
#       pt.cex = 2.2, col = "black", bty ="n", cex = 1.1, y.intersp = 1, x.intersp = .8)
plot.new()
plot.window(xlim=c(0.1,-0.1),ylim=c(1, length(Mainland.tree2$tip.label)))
par(cex=1)
text(rep(0,length(Mainland.tree2$tip.label)), 1:length(Mainland.tree2$tip.label),Mainland.tree$tip.label[ordered_tips], cex = 0.5, font = 3)
plot(tt2,colors=simmap.col2,lwd=3,split.vertical=TRUE,ftype="off",direction = "leftwards")
nodelabels(pie=pd2$ace,piecol=simmap.col2,cex=.6)
tiplabels(pie=simmap.eco2[Mainland.tree2$tip.label,],piecol=simmap.col2,cex=0.4)
dev.off()


# PLOT BOXPLOT ------------------------------------------------------------

EcomorphData <- NewData[which(!NewData$Ground == "M" & !NewData$Ground == "U"),]
ggarrange(ggbox(trait = "SVL", ylab = "Snout-vent Length"),
          ggbox(trait = "TL", ylab = "Tail Length"),
          ggbox(trait = "HL", ylab = "Head Length"),
          ggbox(trait = "SL", ylab = "Snout Length"),
          ggbox(trait = "HW", ylab = "Head Width"),
          ggbox(trait = "Humerus", ylab = "Humerus Length"),
          ggbox(trait = "Radius", ylab = "Radius Length"),
          ggbox(trait = "Hand", ylab = "Hand Length"),
          ggbox(trait = "Femur", ylab = "Femur Length"),
          ggbox(trait = "Tibia", ylab = "Tibia Length"),
          ggbox(trait = "Foot", ylab = "Foot Length"),
          ggbox(trait = "Finger.III.W", ylab = "Fingerpad Width"),
          ggbox(trait = "Toe.III.W", ylab = "Toepad Width"),
          ggbox(trait = "Finger.III.Lam", ylab = "Fingerpad Lamellae"),
          ggbox(trait = "Toe.III.Lam", ylab = "Toepad Lamellae"),
          ncol = 3, nrow = 5)


# PLOT PCA ROCK -----------------------------------------------------------


col <- setNames(c("Blue","cyan", "Gold", "grey73", "ForestGreen", "tan4", "Red", "Purple4", "Darkorange1","Black","Pink"),c(sort(unique(NewData$Ground)),"RW","Rock"))

rock.TG <- c("argenteolus","bartschi","lucius","longitibialis","taylori","gadovii","strahmi","armouri","shrevei")
rock.G <- c("rupinae","monticola","barkeri","aquaticus","rivalis")

phylopca <- phyl.pca(tree, NewData[tree$tip.label,6:18], method = "BM", mode = "cov")

species <- rownames(phylopca$S[c(which(!NewData$Ground == "M" & !NewData$Ground == "U"),
                                 which(!is.na(match(NewData$Species,rock.TG))),
                                 which(!is.na(match(NewData$Species,rock.G)))),]) # use to plot only the ecomorph species
eco <- setNames(NewData[species,"Ground"], species)
eco[rock.TG] <- "RW"
eco[rock.G] <- "Rock"
phylopca$S <- phylopca$S[species,]
phylopca$S[,c(3)] <-  phylopca$S[,c(3)] *-1
phylopca$S[,c(1)] <-  phylopca$S[,c(1)] *-1
phylopca$S[,c(5)] <-  phylopca$S[,c(5)] *1
ggarrange(ggplot.pca(pca = phylopca, axis1 = 1, axis2 =2, species = species, 
                     groups = eco, labels = F, 
                     region = NewData[species,"Region"]),
          ggplot.pca(pca = phylopca, axis1 = 1, axis2 =3, species = species, 
                     groups = eco, labels = F, 
                     region = NewData[species,"Region"]),
          ggplot.pca(pca = phylopca, axis1 = 4, axis2 =3, species = species, 
                     groups = eco, labels = F, 
                     region = NewData[species,"Region"]),
          ggplot.pca(pca = phylopca, axis1 = 4, axis2 =5, species = species, 
                     groups = eco, labels = F, 
                     region = NewData[species,"Region"]),
          ncol =2,nrow =2)

bon <- posthoc.cross(rockscores, axes = 5, fun = "random.dist", p.adj = "bonferroni", nsim = 1000)
rockscores <- cbind("Ecomorph"=eco,as.data.frame(phylopca$S[,1:5]))
rockscores<-rockscores %>%
  filter(Ecomorph == "G" | Ecomorph == "TG" | Ecomorph == "Rock" | Ecomorph == "RW")



# Tanglegram --------------------------------------------------------------



hm <- hclust(dist((NewData[,6:18]),method = "euclidean"), method = "ward.D2", members = NULL)
hm2 <- as.phylo(hm)
plotTree(hm2,cex =0.5, offset = .3, fsize = 0.3, ftype = "i")
var <- rect.hclust(hm, k=(10), border="red")

tangle <- cophylo(ladderize(tree),ladderize(hm2))
connect <- setNames(NewData$Region,NewData$Species)
connect.col <- setNames(c(4,2,2),unique(connect))
thick <- setNames(rep(1.5,205),NewData$Species); thick[which(NewData$Region == "Mainland")] <- 1.5
line <- setNames(rep("solid",205),NewData$Species); line[which(NewData$Region == "Mainland")] <- "solid"


left<-rep(4,nrow(tangle$trees[[1]]$edge))
nodes<-getDescendants(tangle$trees[[1]],212)
left[sapply(nodes,function(x,y) which(y==x),y=tangle$trees[[1]]$edge[,2])]<-
  2
nodes<-getDescendants(tangle$trees[[1]],309)
left[sapply(nodes,function(x,y) which(y==x),y=tangle$trees[[1]]$edge[,2])]<-
  2
left[200] <- 2

right<-rep("black",nrow(tangle$trees[[2]]$edge))
nodes<-getDescendants(tangle$trees[[2]],392)
right[sapply(nodes,function(x,y) which(y==x),y=tangle$trees[[2]]$edge[,2])]<-
  "Blue"
nodes<-getDescendants(tangle$trees[[2]],376)
right[sapply(nodes,function(x,y) which(y==x),y=tangle$trees[[2]]$edge[,2])]<-
  "Purple4"
nodes<-getDescendants(tangle$trees[[2]],368)
right[sapply(nodes,function(x,y) which(y==x),y=tangle$trees[[2]]$edge[,2])]<-
  "ForestGreen"
nodes<-getDescendants(tangle$trees[[2]],334)
right[sapply(nodes,function(x,y) which(y==x),y=tangle$trees[[2]]$edge[,2])]<-
  "Red"
nodes<-getDescendants(tangle$trees[[2]],342)
right[sapply(nodes,function(x,y) which(y==x),y=tangle$trees[[2]]$edge[,2])]<-
  "ForestGreen"
nodes<-getDescendants(tangle$trees[[2]],300)
right[sapply(nodes,function(x,y) which(y==x),y=tangle$trees[[2]]$edge[,2])]<-
  "tan4"
nodes<-getDescendants(tangle$trees[[2]],269)
right[sapply(nodes,function(x,y) which(y==x),y=tangle$trees[[2]]$edge[,2])]<-
  "Gold"
nodes<-getDescendants(tangle$trees[[2]],211)
right[sapply(nodes,function(x,y) which(y==x),y=tangle$trees[[2]]$edge[,2])]<-
  "cyan"
nodes<-getDescendants(tangle$trees[[2]],233)
right[sapply(nodes,function(x,y) which(y==x),y=tangle$trees[[2]]$edge[,2])]<-
  "cyan"


edge.col<-list(left=left,right=right)
group <- setNames(Pred.Eco2$Pred.Eco,Pred.Eco2$Species)
group[which(group == "U")] <- "M"
#group <- setNames(NewData$Ground,NewData$Species)



plot.cophylo(tangle, fsize = 0.3, 
             link.lty = line[tangle$assoc[,1]], 
             link.col = connect.col[connect[tangle$assoc[,1]]], 
             link.lwd = thick[tangle$assoc[,1]],
             tip.lty = "blank",
             edge.col = edge.col, lwd = 1)
tiplabels.cophylo(pch=19,frame="none",col=connect.col[connect[tangle$trees[[1]]$tip.label]],cex=.3)
#tiplabels.cophylo(pie=simmap.eco2[tangle$trees[[2]]$tip.label,],piecol=simmap.col2,cex=0.07, which = "right")
tiplabels.cophylo(pch=19,frame="none",col=col[group[tangle$trees[[2]]$tip.label]],cex=.3, which = "right")
#nodelabels.cophylo(which="right", cex = 0.3)
