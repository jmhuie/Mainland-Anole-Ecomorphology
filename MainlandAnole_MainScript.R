rm(list = ls(all=TRUE))
source("MainlandAnole_Functions.R")
library(tidyverse)
library(readxl)
library(phytools)
library(factoextra)
library(geometry)
library(plyr)
library(ggpubr)
library(MASS)
library(geiger)


# Perch Data --------------------------------------------------------

EcoData <- read_excel("Data/MainlandAnole_SpeciesList.xlsx")%>%
  as.data.frame
rownames(EcoData) <- EcoData$Species 
EcoData <- EcoData[sort(EcoData$Species),]

RawPerch <- read_excel("Data/MainlandAnole_EcoData.xlsx") %>%
  filter(., Quality == "S" | Quality == "AS" ) %>%
  dplyr::select(.,Species,'Height N','Perch Height (m)','Diameter N', 'Perch Diameter (cm)') %>%
  as.data.frame
colnames(RawPerch) <- c("Species", "PH.N", "PH", "PD.N", "PD")
RawPerch[is.na(RawPerch)] <- 1
PerchData <- matrix(nrow =length(unique(RawPerch$Species)), ncol = 3) %>% as.data.frame
colnames(PerchData) <- c("Species", "PH", "PD")
PerchData[,1] <- unique(RawPerch$Species)

for (i in 1:nrow(PerchData)) {
  species <- as.character(PerchData[i,1])
  PerchData[i,2] <- round(weighted.mean(x = as.numeric(unlist(filter(RawPerch, Species == species)[,"PH"])),
                                  w = as.numeric(unlist(as.numeric(filter(RawPerch, Species == species)[,"PH.N"])/
                                  sum(as.numeric(filter(RawPerch, Species == species)[,"PH.N"]))))),2) 
                
  PerchData[i,3] <- round(weighted.mean(x = as.numeric(unlist(filter(RawPerch, Species == species)[,"PD"])),
                                  w = as.numeric(unlist(as.numeric(filter(RawPerch, Species == species)[,"PD.N"])/
                                  sum(as.numeric(filter(RawPerch, Species == species)[,"PD.N"]))))),2)
  remove(species)
}

EcoData <- left_join(EcoData,PerchData)
rownames(EcoData) <- EcoData$Species
remove(PerchData,RawPerch)

#plot(EcoData$PH~EcoData$PD, col = col[EcoData$Ground], pch = 19, cex = 1)
#text(EcoData$PH~EcoData$PD, col = col[EcoData$Ground], label = EcoData$Species)

# Read Tree ---------------------------------------------------------------

tree <- read.nexus("Trees/Poe_2017_MCC_names.txt")
#tree <- read.nexus("Trees/Poe_2017_Time_names.tre")
tree <- ladderize(tree, right = FALSE)
tree <- drop.tip(tree, setdiff(tree$tip.label, unique(EcoData$Species))) #keeps tip labels that matches vector
plotTree(tree, size = .5) # plot tree. looks pretty ultrametric but R doesn't think so
tree <- force.ultrametric(tree) # forces tree to be ultrametric
dev.off()

# Morphology Data Cleaning ------------------------------------------------

MorphoData <- read_excel("Data/MainlandAnole_MorphoData.xlsx", sheet = 1) %>% #Read in full data
#MorphoData <- read.csv("AnoleDataClean.csv", header = TRUE, sep = ",", na.strings = "") %>%
  filter(Sex == "M") %>% #Keeps only Males
  filter(is.na(Omit)) %>% #removes species marked to omit
  dplyr::select(., Collection:Toe.III.Lam, -Trip,-Omit,-Sex) %>% # select relevant columns
  as.data.frame
MorphoData[,4:20] <- round(MorphoData[4:20],2)
MorphoData <- MorphoData[!is.na(match(MorphoData$Species,EcoData$Species)),]


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

#size-correct the data
MorphoData <- MorphoData %>%
  dplyr::group_by(Species) %>% 
  summarise_at(vars(SVL:Toe.III.Lam), list(mean), na.rm = TRUE) %>%# finds the species mean for each trait
  as.data.frame
rownames(MorphoData) <- MorphoData$Species
MorphoData<-MorphoData[tree$tip.label,]


NewData <- MorphoData[tree$tip.label,]
all(match(NewData$Species, tree$tip.label))
phyl.resid <- phyl.resid(tree, x = setNames(NewData$SVL,NewData$Species), Y = NewData[,4:19]) # phylo size-correction

NewData[,3] <-tail.resid$resid
NewData[,4:19]<-phyl.resid$resid
NewData <- NewData[,c(1:3,5:9,11:13,15,17,19,21,23)]
NewData <- cbind(EcoData[tree$tip.label,-c(4,6,7)],NewData[tree$tip.label,-1])
remove(tail,tail.resid,phyl.resid,MorphoData)

# PCA ---------------------------------------------------------------------

phylopca <- phyl.pca(tree, NewData[tree$tip.label,5:17], method = "BM", mode = "cov")
phylopca$L
plot.phyl.pca(phylopca, type = "lines")

col <- setNames(c("Blue","cyan", "Gold", "grey73", "ForestGreen", "tan4", "Red", "Purple4", "Darkorange1"),sort(unique(NewData$Ground)))
col2 <- setNames(c("Blue", "Gold", "grey73", "ForestGreen", "tan4", "Red", "Purple4", "Darkorange1"),sort(unique(NewData$Ecomorph)))

species <- rownames(phylopca$S)
eco <- setNames(NewData[species,"Ground"], species)
phylopca$S[,c(2)] <-  phylopca$S[,c(2)] *-1
ggarrange(ggplot.pca(pca = phylopca, axis1 = 1, axis2 =2, species = species, 
                groups = eco, labels = F, 
                region = NewData[species,"Region"]),
          ggplot.pca(pca = phylopca, axis1 = 1, axis2 =3, species = species, 
                groups = eco, labels = F, 
                region = NewData[species,"Region"]),
          ggplot.pca(pca = phylopca, axis1 = 4, axis2 =3, species = species, 
                groups = eco, labels = F, 
                region = NewData[species,"Region"]) + 
                scale_y_continuous(labels = scales::number_format(accuracy = 0.1)),
          ggplot.pca(pca = phylopca, axis1 = 4, axis2 =5, species = species, 
                groups = eco, labels = F, 
                region = NewData[species,"Region"]),
          ncol =2,nrow =2)
#ggsave(tail.phy.pc, filename = "PCA2.pdf",  bg = "transparent", height = 15, width = 15)
#dev.off()



# DFA w/ Caribbean --------------------------------------------------------

EcomorphData <- NewData[which(!NewData$Ecomorph == "M" & !NewData$Ecomorph == "U"),]
EcomorphData <- droplevels(EcomorphData)

summary(manova(as.matrix(EcomorphData[,c(5:17)])~EcomorphData$Ecomorph), test = "Wilks")

LDA <- lda(EcomorphData$Ecomorph ~ ., 
           data = EcomorphData[,c(5:17)],
           prior = c(rep(1/6, 6)), CV = FALSE) # initial LDA

assessment <- predict(LDA)
#acc <- table(EcomorphData$Ecomorph, LDA$class) # actual vs predicted with CV
acc <- table(EcomorphData$Ecomorph, assessment$class) # actual vs predicted without CV
acc
sum(acc[row(acc) == col(acc)]) / sum(acc)

EcomorphData$Species[which(!assessment$class == EcomorphData$Ecomorph)]
#EcomorphData$Species[which(!LDA$class == EcomorphData$Ecomorph)]
#LDA$class[which(!LDA$class == EcomorphData$Ecomorph)]

predictions <- predict(LDA, as.data.frame(NewData[which(NewData$Ecomorph == "M" | NewData$Ecomorph == "U"),]))

criteria.lda <- cbind(Species = rownames(predictions$x), 
                      Ecomorph = as.character(NewData$Ecomorph[which(NewData$Ecomorph == "M" | NewData$Ecomorph == "U")]),
                      round(predictions$posterior,3), 
                      Pred.Eco = as.character(predictions$class))
#removes species with post prob >90
criteria.lda.trim <- criteria.lda
for (i in 1:nrow(criteria.lda.trim)) {
  tmp <- max(criteria.lda.trim[i,4:ncol(criteria.lda.trim)-1])
  if (tmp < 0.90) {
    criteria.lda.trim[i,"Pred.Eco"] <- NA
  }
}
criteria.lda.trim <- as.data.frame(criteria.lda.trim[which(!is.na(criteria.lda.trim[,"Pred.Eco"])),])

tapply(criteria.lda.trim$Pred.Eco[which(criteria.lda.trim$Ecomorph == "M" | criteria.lda.trim$Ecomorph == "U")],
       criteria.lda.trim$Ecomorph[which(criteria.lda.trim$Ecomorph == "M" | criteria.lda.trim$Ecomorph == "U")],length)
tapply(criteria.lda.trim$Pred.Eco[which(criteria.lda.trim$Ecomorph == "M")],
       criteria.lda.trim$Pred.Eco[which(criteria.lda.trim$Ecomorph == "M")],length)



# ED Criteria w/ Caribbean ------------------------------------------------

# predict ecomorphs based on certain criteria

predicted <- ED.predict(scores = phylopca$S[tree$tip.label,1:5], 
                           species = tree$tip.label, 
                           groups = setNames(NewData[tree$tip.label,"Ecomorph"],tree$tip.label),
                           hard.mode = T, all.species = F)
criteria1 <- predicted$criteria1
tapply(criteria1$Pred.Eco[which(criteria1$Ecomorph=="M")],criteria1$Pred.Eco[which(criteria1$Ecomorph=="M")],length)
criteria2 <- predicted$criteria2
tapply(criteria2$Pred.Eco[which(criteria2$Ecomorph=="M")],criteria2$Pred.Eco[which(criteria2$Ecomorph=="M")],length)
criteria3 <- predicted$criteria3
tapply(criteria3$Pred.Eco[which(criteria3$Ecomorph=="M")],criteria3$Pred.Eco[which(criteria3$Ecomorph=="M")],length)


# Compile w/Caribbean -----------------------------------------------------

compile <- synth.compile(lda = criteria.lda, predicted = predicted, hard.mode =  T,)
compile
tapply(compile$Predicted,compile$Ecomorph,length)
tapply(compile$Predicted[which(compile$Ecomorph == "M")],compile$Predicted[which(compile$Ecomorph == "M")],length)


Pred.Eco <- data.frame(Species = NewData$Species, Region = NewData$Region, Pred.Eco = NewData$Ecomorph)
rownames(Pred.Eco) <- Pred.Eco$Species
Pred.Eco[!is.na(match(Pred.Eco$Species,compile$Species)),"Pred.Eco"] <- compile$Predicted

ggplot.pca(pca = phylopca, axis1 = 1, axis2 =2, species = tree$tip.label, 
           groups = Pred.Eco$Pred.Eco, labels = FALSE)


# Intermediate w/Caribbean ------------------------------------------------
noclass <- as.data.frame(criteria.lda[is.na(match(criteria.lda[,1],compile[,1])),])
intermediate <- cbind(noclass[,c("Species","Ecomorph")],"DFA"=NA,"C1"=NA,"C2"=NA,"C3"=NA,"Match" =NA)

predicted <- ED.predict(scores = phylopca$S[tree$tip.label,1:5], 
                        species = tree$tip.label, 
                        groups = setNames(NewData[tree$tip.label,"Ecomorph"],tree$tip.label),
                        hard.mode = T, all.species = F)
criteria1 <- predicted$criteria1
criteria2 <- predicted$criteria2
criteria3 <- predicted$criteria3

for (i in 1:nrow(intermediate)) {
  species <- intermediate[i,"Species"]
  tmp <- sort((noclass[i,4:ncol(noclass)-1]),decreasing = TRUE)
  tmp2 <- sort((criteria1[species,4:ncol(noclass)-1]),decreasing = F)
  tmp3 <- sort((criteria2[species,4:ncol(noclass)-1]),decreasing = F)
  tmp4 <- sort((criteria3[species,4:ncol(noclass)-1]),decreasing = F)
  if (as.numeric(tmp[1]) >= 0.9) {
    intermediate$DFA[i] <- paste0(names(tmp[1]))
  #  intermediate$C1[i] <- paste0(names(tmp2[1]))
  #  intermediate$C2[i] <- paste0(names(tmp3[1]))
  #  intermediate$C3[i] <- paste0(names(tmp4[1]))
  }
  if ((as.numeric(tmp[1]) + as.numeric(tmp[2])) >= 0.900 & as.numeric(tmp[1]) < 0.9) {
    intermediate$DFA[i] <- paste0(names(tmp[1]),"/",names(tmp[2]))
    if (length(tmp2) >=2){
      intermediate$C1[i] <- paste0(names(tmp2[1]),"/",names(tmp2[2]))
    }
    if (length(tmp2) ==1){
      intermediate$C1[i] <- paste0(names(tmp2[1]))
    }
    if (length(tmp3) >=2){
      intermediate$C2[i] <- paste0(names(tmp3[1]),"/",names(tmp3[2]))
    }
    if (length(tmp3) ==1){
      intermediate$C2[i] <- paste0(names(tmp3[1]))
    }
    if (length(tmp4) >=2){
      intermediate$C3[i] <- paste0(names(tmp4[1]),"/",names(tmp4[2]))
    }
    if (length(tmp4) ==1){
      intermediate$C3[i] <- paste0(names(tmp4[1]))
    }
  }
  #if ((as.numeric(tmp[1]) + as.numeric(tmp[2]) + as.numeric(tmp[3])) >= 0.900 & as.numeric(tmp[1]) + as.numeric(tmp[2])< 0.9) {
    #intermediate$DFA[i] <- paste0(names(tmp[1]),"/",names(tmp[2]),"/",names(tmp[3]))
    #intermediate$C1[i] <- paste0(names(tmp2[1]),"/",names(tmp2[2]),"/",names(tmp2[3]))
    #intermediate$C2[i] <- paste0(names(tmp3[1]),"/",names(tmp3[2]),"/",names(tmp3[3]))
    #intermediate$C3[i] <- paste0(names(tmp4[1]),"/",names(tmp4[2]),"/",names(tmp4[3]))
  #}
  match <- any(match(str_split(intermediate[i,], "/")[[3]],str_split(intermediate[i,], "/")[[4]]))
  if(isTRUE(match) == TRUE) {
    intermediate$Match[i] <- intermediate$DFA[i]
  } else {
    intermediate$Match[i] <- NA
  }
}
intermediate <- filter(intermediate, str_detect(intermediate$DFA, "/"))
intermediate <- intermediate[!is.na(intermediate$Match),c("Species","Ecomorph","Match")]

# DFA w/ Ground -----------------------------------------------------------

EcomorphData <- NewData[which(!NewData$Ground == "M" & !NewData$Ground == "U"),]
EcomorphData <- droplevels(EcomorphData)

summary(manova(as.matrix(EcomorphData[,c(5:17)])~EcomorphData$Ground), test = "Wilks")

LDA <- lda(EcomorphData$Ground ~ ., 
           data = EcomorphData[,c(5:17)],
           prior = c(rep(1/7, 7)), CV = FALSE) # initial LDA

assessment <- predict(LDA)
#acc <- table(EcomorphData$Ground, LDA$class) # actual vs predicted with CV
acc <- table(EcomorphData$Ground, assessment$class) # actual vs predicted without CV
acc
sum(acc[row(acc) == col(acc)]) / sum(acc)

EcomorphData$Species[which(!assessment$class == EcomorphData$Ground)]
#EcomorphData$Species[which(!LDA$class == EcomorphData$Ground)]
#LDA$class[which(!LDA$class == EcomorphData$Ground)]

predictions <- predict(LDA, as.data.frame(NewData[which(NewData$Ground == "M" | NewData$Ground == "U"),]))

criteria.lda <- cbind(Species = rownames(predictions$x), 
                      Ecomorph = as.character(NewData$Ground[which(NewData$Ground == "M" | NewData$Ground == "U")]),
                      round(predictions$posterior,2), 
                      Pred.Eco = as.character(predictions$class))

#removes species with post prob <90
criteria.lda.trim <- criteria.lda
for (i in 1:nrow(criteria.lda.trim)) {
  tmp <- max(criteria.lda.trim[i,4:ncol(criteria.lda.trim)-1])
  if (tmp > .90) {
    criteria.lda.trim[i,"Pred.Eco"] <- NA
  }
}
criteria.lda.trim <- as.data.frame(criteria.lda.trim[which(!is.na(criteria.lda.trim[,"Pred.Eco"])),])

tapply(criteria.lda.trim$Pred.Eco[which(criteria.lda.trim$Ecomorph == "M" | criteria.lda.trim$Ecomorph == "U")],
       criteria.lda.trim$Ecomorph[which(criteria.lda.trim$Ecomorph == "M" | criteria.lda.trim$Ecomorph == "U")],length)
tapply(criteria.lda.trim$Pred.Eco[which(criteria.lda.trim$Ecomorph == "M")],
       criteria.lda.trim$Pred.Eco[which(criteria.lda.trim$Ecomorph == "M")],length)



# ED Criteria w/ Ground ---------------------------------------------------

# predict ecomorphs based on certain criteria

predicted <- ED.predict(scores = phylopca$S[,1:5], species = tree$tip.label, 
                           groups = NewData[tree$tip.label,"Ground"],
                           hard.mode = T, all.species = F)
criteria1 <- predicted$criteria1
tapply(criteria1$Pred.Eco[which(criteria1$Ecomorph=="M")],criteria1$Pred.Eco[which(criteria1$Ecomorph=="M")],length)
criteria2 <- predicted$criteria2
tapply(criteria2$Pred.Eco[which(criteria2$Ecomorph=="M")],criteria2$Pred.Eco[which(criteria2$Ecomorph=="M")],length)
criteria3 <- predicted$criteria3
tapply(criteria3$Pred.Eco[which(criteria3$Ecomorph=="M")],criteria3$Pred.Eco[which(criteria3$Ecomorph=="M")],length)


# Compile w/ Ground -------------------------------------------------------

compile2 <- synth.compile(lda = criteria.lda, predicted = predicted, upper.cut = 0.95, lower.cut = 0.9,
                        hard.mode = T)

tapply(compile2$Predicted,compile2$Ecomorph,length)
tapply(compile2$Predicted[which(compile2$Ecomorph == "M")],
       compile2$Predicted[which(compile2$Ecomorph == "M")],length)
compile$Species[which(is.na(match(compile$Species, compile2$Species)))]

#reconcile many of the species not classified with ground ecomorph but were previously classified
comp.diff<-function(){
  diff <- compile$Species[which(is.na(match(compile$Species, compile2$Species)))]
  diff <- diff[which(is.na(match(diff,EcomorphData$Species)))]
  diff
  
  old <- compile[!is.na(match(compile$Species,diff)),]; old
  rownames(old) <- old$Species
  new <- as.data.frame(criteria.lda[diff,]);new
  criteria1[diff,]
  just <- cbind(old[,c("Species","Predicted")],New.LDA = new[,"Pred.Eco"], C1.Predicted = criteria1[diff,"Pred.Eco"])
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
  trans <-just$Species[which(just$New.LDA == "G" & just$C1.Predicted == just$New.LDA)]
  total <- synth.compile(lda = criteria.lda, predicted = predicted, hard.mode = F)
  trans <- total[!is.na(match(total$Species,trans)),]
  trans <- trans[which(trans$LDA >= 0.75),]
  rownames(trans) <- trans$Species
  
  new.compile2 <- rbind(compile2,synth,trans)
  new.compile2 <- new.compile2[sort(new.compile2$Species),]
  
  return(new.compile2)
}
new.compile2 <- comp.diff()

compile$Species[which(is.na(match(compile$Species, new.compile2$Species)))]

tapply(new.compile2$Predicted,new.compile2$Ecomorph,length)
tapply(new.compile2$Predicted[which(new.compile2$Ecomorph == "M")],
       new.compile2$Predicted[which(new.compile2$Ecomorph == "M")],length)

Pred.Eco2 <- data.frame(Species = NewData$Species, Pred.Eco = NewData$Ground)
rownames(Pred.Eco2) <- Pred.Eco2$Species
Pred.Eco2[!is.na(match(Pred.Eco2$Species,compile2$Species)),2] <- compile2$Predicted

ggplot.pca(pca = phylopca, axis1 = 1, axis2 =3, species = tree$tip.label, 
           groups = Pred.Eco2$Pred.Eco, labels = FALSE)

# Intermediate w/Ground --------------------------------------------------
noclass <- as.data.frame(criteria.lda[is.na(match(criteria.lda[,1],new.compile2[,1])),])
intermediate <- cbind(noclass[,c("Species","Ecomorph")],"DFA"=NA,"C1"=NA,"C2"=NA,"C3"=NA,"Match" =NA)


predicted <- ED.predict(scores = phylopca$S[tree$tip.label,1:5], 
                        species = tree$tip.label, 
                        groups = setNames(NewData[tree$tip.label,"Ground"],tree$tip.label),
                        hard.mode = F, all.species = F)
criteria1 <- predicted$criteria1
criteria2 <- predicted$criteria2
criteria3 <- predicted$criteria3

for (i in 1:nrow(intermediate)) {
  species <- intermediate[i,"Species"]
  tmp <- sort((noclass[i,4:ncol(noclass)-1]),decreasing = TRUE)
  tmp2 <- sort((criteria1[species,4:ncol(noclass)-1]),decreasing = F)
  tmp3 <- sort((criteria2[species,4:ncol(noclass)-1]),decreasing = F)
  tmp4 <- sort((criteria3[species,4:ncol(noclass)-1]),decreasing = F)
  if (as.numeric(tmp[1]) >= 0.9) {
    intermediate$DFA[i] <- paste0(names(tmp[1]))
    #  intermediate$C1[i] <- paste0(names(tmp2[1]))
    #  intermediate$C2[i] <- paste0(names(tmp3[1]))
    #  intermediate$C3[i] <- paste0(names(tmp4[1]))
  }
  if ((as.numeric(tmp[1]) + as.numeric(tmp[2])) >= 0.900 & as.numeric(tmp[1]) < 0.9) {
    intermediate$DFA[i] <- paste0(names(tmp[1]),"/",names(tmp[2]))
    if (length(tmp2) >=2){
      intermediate$C1[i] <- paste0(names(tmp2[1]),"/",names(tmp2[2]))
    }
    if (length(tmp2) ==1){
      intermediate$C1[i] <- paste0(names(tmp2[1]))
    }
    if (length(tmp3) >=2){
      intermediate$C2[i] <- paste0(names(tmp3[1]),"/",names(tmp3[2]))
    }
    if (length(tmp3) ==1){
      intermediate$C2[i] <- paste0(names(tmp3[1]))
    }
    if (length(tmp4) >=2){
      intermediate$C3[i] <- paste0(names(tmp4[1]),"/",names(tmp4[2]))
    }
    if (length(tmp4) ==1){
      intermediate$C3[i] <- paste0(names(tmp4[1]))
    }
  }
  #if ((as.numeric(tmp[1]) + as.numeric(tmp[2]) + as.numeric(tmp[3])) >= 0.900 & as.numeric(tmp[1]) + as.numeric(tmp[2])< 0.9) {
  #intermediate$DFA[i] <- paste0(names(tmp[1]),"/",names(tmp[2]),"/",names(tmp[3]))
  #intermediate$C1[i] <- paste0(names(tmp2[1]),"/",names(tmp2[2]),"/",names(tmp2[3]))
  #intermediate$C2[i] <- paste0(names(tmp3[1]),"/",names(tmp3[2]),"/",names(tmp3[3]))
  #intermediate$C3[i] <- paste0(names(tmp4[1]),"/",names(tmp4[2]),"/",names(tmp4[3]))
  #}
  match <- all(match(str_split(intermediate[i,], "/")[[3]],str_split(intermediate[i,], "/")[[4]]))
  if(isTRUE(match) == TRUE) {
    intermediate$Match[i] <- intermediate$DFA[i]
  } else {
    intermediate$Match[i] <- NA
  }
}
intermediate <- filter(intermediate, str_detect(intermediate$DFA, "/"))
intermediate <- intermediate[!is.na(intermediate$Match),c("Species","Ecomorph","Match")]
intermediate[is.na(match(intermediate$Species,new.compile2$Species)),]


# Randomization Tests -----------------------------------------------------

Ecomorph.Scores <- cbind("Ecomorph" = NewData$Ground,as.data.frame(phylopca$S)) %>%
  filter(., !Ecomorph == "M" & !Ecomorph == "U")
Ecomorph.Scores <- NewData[,-c(1:2,4,18:19)] %>%
  filter(., !Ground == "M" & !Ground == "U")

bon <- posthoc.cross(Ecomorph.Scores, axes = 13, fun = "random.dist", p.adj = "bonferroni", nsim = 1000)

# Simulated Trait Data ----------------------------------------------------

# all ecomorph comparisons with multi comparison correction
bon <- posthoc.cross(Ecomorph.Scores, axes = 13, fun = "sim.dist", p.adj = "holm", nsim = 1000)

# simulate data and test compare distance to centroid and NND against null
var <- sim.ED(tree, scores = Ecomorph.Scores, axes = 13, nsim = 1000)
p.adjust(var[,4], method = "bonferroni")


# SIMMAP Analyses w/Caribbean  --------------------------------------------

simmap.eco <- setNames(NewData$Ecomorph,NewData$Species)
simmap.eco[which(simmap.eco == "M")] <- "U"
simmap.eco[compile$Species] <- compile$Predicted
simeco.col <- setNames(c("Blue", "Gold",  "ForestGreen", "tan4","Red", "Purple4","grey73"),sort(unique(simmap.eco)))
simmap.eco<-to.matrix(simmap.eco,levels(as.factor(simmap.eco)))

for (i in 1:nrow(intermediate)) {
  #split <- str_split(intermediate[i,"Match"], "/")
  #simmap.eco[intermediate$Species[i],c(split)[[1]]] <- 0.5
  simmap.eco[intermediate$Species[i],-ncol(simmap.eco)] <- as.numeric(criteria.lda[intermediate$Species[i],3:(ncol(criteria.lda)-1)])
  simmap.eco[intermediate$Species[i],"U"] <- 0
}
# ER -276.4185 df 1
#> SYM -244.3273 df 21
#> ER vs SYM LR 64.1824 df 20 
#-2*(-276.4185 - (-244.3273))
#quantile(rchisq(10000, 20),.95) #SYM IS BETTER!!!! YAY
# ARD -235.6913 df 42
# SYM vs ARD LR 17.272 df 21
#-2*(-244.3273 - (-235.6913))
#quantile(rchisq(10000, 21),.95) #ARD IS NOT BETTER!!!! GO FOR SYM

simmap<-make.simmap(tree,simmap.eco,model = "SYM",nsim = 1)
simmap.trim <- lapply(simmap,drop.tip.simmap,NewData$Species[which(!NewData$Region == "Mainland")])
class(simmap.trim) <- "multiPhylo"
mainland.tree <- drop.tip(tree, setdiff(tree$tip.label,NewData$Species[which(NewData$Region == "Mainland")]))
pd<-describe.simmap(simmap.trim[[1]], plot = FALSE, ref.tree = mainland.tree)
plot(simmap.trim, fsize = 0.5, col = simeco.col)
nodelabels(pie=pd$ace,piecol = simeco.col,cex=.4)
tiplabels(pie=simmap.eco[simmap.trim[[1]]$tip.label,],piecol = simeco.col,cex=0.2)

transitions <- pd$count
transition.type <- (",TG$")
tapply(apply(transitions[,grep(transition.type, colnames(transitions))],1,sum),
       apply(transitions[,grep(transition.type, colnames(transitions))],1,sum),length)


# L1ou attempt ------------------------------------------------------------
tree <- force.ultrametric(tree)
lizard <- adjust_data(tree,as.matrix(phylopca$S[,1:5]))
fit_ind_AIC <- estimate_shift_configuration(lizard$tree,lizard$Y, nCores = 7, criterion = "AICc")
fit_ind_AIC_bootstrap <- l1ou_bootstrap_support(fit_ind_AIC,nItrs = 100, multicore = T, nCores=4)

fit_conv_AIC <- estimate_convergent_regimes(fit_ind_AIC, criterion = "AICc", nCores = 7)
plot(fit_conv_AIC, show.data	= FALSE)

fit_conv_AIC$shift.configuration
plot(fit_ind_AIC, edge.ann.cex=1, cex=0.5, label.offset=0.02, edge.label.ann = FALSE)


nEdges <- Nedge(lizard$tree) # total number of edges
ew <- rep(1,nEdges)  # to set default edge width of 1
ew[fit_ind_AIC$shift.configuration] <- 3   # to widen edges with a shift 
plot(fit_ind_AIC, cex=0.5, label.offset=0.02, edge.width=ew)

