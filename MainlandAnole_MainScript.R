rm(list = ls(all=TRUE))
source("MainlandAnole_Functions.R")
library(tidyverse)
library(readxl)
library(phytools)
library(factoextra)
library(geometry)
library(plyr)

# Read Tree ---------------------------------------------------------------

EcoData <- read_excel("MainlandAnole_SpeciesList.xlsx")%>%
  as.data.frame
rownames(EcoData) <- EcoData$Species 
EcoData <- EcoData[sort(EcoData$Species),]

tree <- read.nexus("Trees/Poe_2017_MCC_names.txt")
tree <- ladderize(tree, right = FALSE)
tree <- drop.tip(tree, setdiff(tree$tip.label, unique(EcoData$Species))) #keeps tip labels that matches vector

#plotTree(tree, size = .5)

# Morphology Data Cleaning -----------------------------------------------------------

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

NewData <- MorphoData[tree$tip.label,]
all(match(NewData$Species, tree$tip.label))
phyl.resid <- phyl.resid(tree, x = setNames(NewData$SVL,NewData$Species), Y = NewData[,4:19]) # phylo size-correction

NewData[,3] <-tail.resid$resid
NewData[,4:19]<-phyl.resid$resid
NewData <- NewData[,c(1:3,5:9,11:13,15,17,19,21,23)]
NewData <- cbind(EcoData,NewData[EcoData$Species,-1])
remove(tail,tail.resid,phyl.resid)


# PCA ---------------------------------------------------------------------

phylopca <- phyl.pca(tree, NewData[,5:17], method = "BM", mode = "cov")
plot.phyl.pca(phylopca)
summary(phylopca)

col <- setNames(c("Blue", "Gold","cyan","grey73", "Red", "ForestGreen", "tan4", "Purple4", "Darkorange1"),sort(unique(NewData$Ground)))
col2 <- setNames(c("Blue","Gold","grey73", "Red", "ForestGreen", "tan4", "Purple4", "Darkorange1"),sort(unique(NewData$Ecomorph)))

species = rownames(phylopca$S)
eco <- setNames(NewData[species,"Ground"], species)
phylopca$S[,c(2)] <-  phylopca$S[,c(2)] *-1
ggarrange(ggplot.pca(pca = phylopca, axis1 = 1, axis2 =2, species = species, 
                groups = eco, labels = FALSE, 
                region = NewData[species,"Region"]),
          ggplot.pca(pca = phylopca, axis1 = 1, axis2 =3, species = species, 
                groups = eco, labels = FALSE, 
                region = NewData[species,"Region"]),
          ggplot.pca(pca = phylopca, axis1 = 4, axis2 =3, species = species, 
                groups = eco, labels = FALSE, 
                region = NewData[species,"Region"]) + 
                scale_y_continuous(labels = scales::number_format(accuracy = 0.1)),
          ggplot.pca(pca = phylopca, axis1 = 4, axis2 =5, species = species, 
                groups = eco, labels = FALSE, 
                region = NewData[species,"Region"]),
          ncol =2,nrow =2)
#ggsave(tail.phy.pc, filename = "PCA2.pdf",  bg = "transparent", height = 15, width = 15)
#dev.off()



# DFA w/ Caribbean --------------------------------------------------------

EcomorphData <- NewData[which(!NewData$Ecomorph == "Mainland" & !NewData$Ecomorph == "Unknown"),]
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

predictions <- predict(LDA, as.data.frame(NewData[which(NewData$Ecomorph == "Mainland" | NewData$Ecomorph == "Unknown"),]))

criteria.lda <- cbind(Species = rownames(predictions$x), 
                      Ecomorph = as.character(NewData$Ecomorph[which(NewData$Ecomorph == "Mainland" | NewData$Ecomorph == "Unknown")]),
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

tapply(criteria.lda.trim$Pred.Eco[which(criteria.lda.trim$Ecomorph == "Mainland" | criteria.lda.trim$Ecomorph == "Unknown")],
       criteria.lda.trim$Ecomorph[which(criteria.lda.trim$Ecomorph == "Mainland" | criteria.lda.trim$Ecomorph == "Unknown")],length)
tapply(criteria.lda.trim$Pred.Eco[which(criteria.lda.trim$Ecomorph == "Mainland")],
       criteria.lda.trim$Pred.Eco[which(criteria.lda.trim$Ecomorph == "Mainland")],length)



# ED Criteria w/ Caribbean ------------------------------------------------

# predict ecomorphs based on certain criteria

predicted <- predict.class(scores = phylopca$S[tree$tip.label,1:5], species = tree$tip.label, 
                           groups = setNames(NewData[tree$tip.label,"Ecomorph"],tree$tip.label))
criteria1 <- predicted$criteria1
criteria2 <- predicted$criteria2
tapply(criteria2$Pred.Eco,criteria2$Pred.Eco,length)
criteria3 <- predicted$criteria3
tapply(criteria3$Pred.Eco,criteria3$Pred.Eco,length)
criteria4 <- predicted$criteria4
tapply(criteria4$Pred.Eco,criteria4$Pred.Eco,length)



# Compile w/Caribbean -----------------------------------------------------

compile <- compilation(lda = criteria.lda, c1=criteria1, c2 = criteria2, 
                       c3 = criteria3, c4 = criteria4)
tapply(compile[[1]]$Predicted,compile[[1]]$Ecomorph,length)
tapply(compile[[1]]$Predicted[which(compile[[1]]$Ecomorph == "Mainland")],compile[[1]]$Predicted[which(compile[[1]]$Ecomorph == "Mainland")],length)


Pred.Eco <- data.frame(Species = NewData$Species, Pred.Eco = NewData$Ecomorph)
#all(match(Pred.Eco[!is.na(match(Pred.Eco$Species,compile[[1]]$Species,compile[[1]]$Species))
Pred.Eco[!is.na(match(Pred.Eco$Species,compile[[1]]$Species)),2] <- compile[[1]]$Predicted

ggplot.pca(pca = phylopca, axis1 = 1, axis2 =2, species = tree$tip.label, 
           groups = Pred.Eco$Pred.Eco, labels = FALSE)




# DFA w/ Ground -----------------------------------------------------------

EcomorphData <- NewData[which(!NewData$Ground == "Mainland" & !NewData$Ground == "Unknown"),]
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

predictions <- predict(LDA, as.data.frame(NewData[which(NewData$Ground == "Mainland" | NewData$Ground == "Unknown"),]))

criteria.lda <- cbind(Species = rownames(predictions$x), 
                      Ecomorph = as.character(NewData$Ground[which(NewData$Ground == "Mainland" | NewData$Ground == "Unknown")]),
                      round(predictions$posterior,3), 
                      Pred.Eco = as.character(predictions$class))

#removes species with post prob <90
criteria.lda.trim <- criteria.lda
for (i in 1:nrow(criteria.lda.trim)) {
  tmp <- max(criteria.lda.trim[i,4:ncol(criteria.lda.trim)-1])
  if (tmp < 0.90) {
    criteria.lda.trim[i,"Pred.Eco"] <- NA
  }
}
criteria.lda.trim <- as.data.frame(criteria.lda.trim[which(!is.na(criteria.lda.trim[,"Pred.Eco"])),])

tapply(criteria.lda.trim$Pred.Eco[which(criteria.lda.trim$Ecomorph == "Mainland" | criteria.lda.trim$Ecomorph == "Unknown")],
       criteria.lda.trim$Ecomorph[which(criteria.lda.trim$Ecomorph == "Mainland" | criteria.lda.trim$Ecomorph == "Unknown")],length)
tapply(criteria.lda.trim$Pred.Eco[which(criteria.lda.trim$Ecomorph == "Mainland")],
       criteria.lda.trim$Pred.Eco[which(criteria.lda.trim$Ecomorph == "Mainland")],length)



# ED Criteria w/ Ground ---------------------------------------------------

# predict ecomorphs based on certain criteria

predicted <- predict.class(scores = phylopca$S[,1:6], species = tree$tip.label, 
                           groups = NewData$Ground)
criteria1 <- predicted$criteria1
criteria2 <- predicted$criteria2
tapply(criteria2$Pred.Eco,criteria2$Pred.Eco,length)
criteria3 <- predicted$criteria3
tapply(criteria3$Pred.Eco,criteria3$Pred.Eco,length)
criteria4 <- predicted$criteria4
tapply(criteria4$Pred.Eco,criteria4$Pred.Eco,length)



# Compile w/ Ground -------------------------------------------------------

compile2 <- compilation(lda = criteria.lda, c1=criteria1, c2 = criteria2, 
                        c3 = criteria3, c4 = criteria4)
compile[[1]]$Species[which(is.na(match(compile[[1]]$Species, compile2[[1]]$Species)))]


tapply(compile2[[1]]$Predicted,compile2[[1]]$Ecomorph,length)
tapply(compile2[[1]]$Predicted[which(compile2[[1]]$Ecomorph == "Mainland")],
       compile2[[1]]$Predicted[which(compile2[[1]]$Ecomorph == "Mainland")],length)
#Pred.ground <- Pred.ground[which(Pred.ground$Ecomorph == "Unknown"),]



