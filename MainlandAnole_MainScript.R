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

predicted <- predict.class(scores = phylopca$S[tree$tip.label,1:5], species = tree$tip.label, 
                           groups = setNames(NewData[tree$tip.label,"Ecomorph"],tree$tip.label),
                           hard.mode = T, all.species = F)
criteria1 <- predicted$criteria1
criteria2 <- predicted$criteria2
tapply(criteria2$Pred.Eco[which(criteria2$Ecomorph=="M")],criteria2$Pred.Eco[which(criteria2$Ecomorph=="M")],length)
criteria3 <- predicted$criteria3
tapply(criteria3$Pred.Eco[which(criteria3$Ecomorph=="M")],criteria3$Pred.Eco[which(criteria3$Ecomorph=="M")],length)
criteria4 <- predicted$criteria4
tapply(criteria4$Pred.Eco[which(criteria4$Ecomorph=="M")],criteria4$Pred.Eco[which(criteria4$Ecomorph=="M")],length)



# Compile w/Caribbean -----------------------------------------------------

compile <- compilation(lda = criteria.lda, predicted = predicted, upper.cut = 0.95, lower.cut = 0.9)
compile
tapply(compile$Predicted,compile$Ecomorph,length)
tapply(compile$Predicted[which(compile$Ecomorph == "M")],compile$Predicted[which(compile$Ecomorph == "M")],length)


Pred.Eco <- data.frame(Species = NewData$Species, Region = NewData$Region, Pred.Eco = NewData$Ecomorph)
rownames(Pred.Eco) <- Pred.Eco$Species
#all(match(Pred.Eco[!is.na(match(Pred.Eco$Species,compile[[1]]$Species,compile[[1]]$Species))
Pred.Eco[!is.na(match(Pred.Eco$Species,compile$Species)),"Pred.Eco"] <- compile$Predicted
#Pred.Eco[criteria3[which(!is.na(criteria3$Pred.Eco)),1],"Pred.Eco"] <- criteria3[criteria3[which(!is.na(criteria3$Pred.Eco)),1],"Pred.Eco"] #Irshick and Losos 97
#Pred.Eco[criteria2[which(!is.na(criteria2$Pred.Eco)),1],"Pred.Eco"] <- criteria2[criteria2[which(!is.na(criteria2$Pred.Eco)),1],"Pred.Eco"] #Irshick and Losos 97
#tapply(filter(Pred.Eco,Region == "Mainland")$Pred.Eco,filter(Pred.Eco,Region == "Mainland")$Pred.Eco,length) #Irshick and Losos 97


ggplot.pca(pca = phylopca, axis1 = 1, axis2 =2, species = tree$tip.label, 
           groups = Pred.Eco$Pred.Eco, labels = FALSE)

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
  if (tmp < .90) {
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

predicted <- predict.class(scores = phylopca$S[,1:5], species = tree$tip.label, 
                           groups = NewData$Ground,
                           hard.mode = T, all.species = F)
criteria1 <- predicted$criteria1
criteria2 <- predicted$criteria2
tapply(criteria2$Pred.Eco[which(criteria2$Ecomorph=="M")],criteria2$Pred.Eco[which(criteria2$Ecomorph=="M")],length)
criteria3 <- predicted$criteria3
tapply(criteria3$Pred.Eco[which(criteria3$Ecomorph=="M")],criteria3$Pred.Eco[which(criteria3$Ecomorph=="M")],length)
criteria4 <- predicted$criteria4
tapply(criteria4$Pred.Eco[which(criteria4$Ecomorph=="M")],criteria4$Pred.Eco[which(criteria4$Ecomorph=="M")],length)



# Compile w/ Ground -------------------------------------------------------

compile2 <- compilation(lda = criteria.lda, predicted = predicted, upper.cut = 0.95, lower.cut = 0.9,
                        hard.mode = T)

tapply(compile2$Predicted,compile2$Ecomorph,length)
tapply(compile2$Predicted[which(compile2$Ecomorph == "M")],
       compile2$Predicted[which(compile2$Ecomorph == "M")],length)

#reconcile many of the species not classified with ground ecomorph but were previously classified
comp.diff<-function(){
  rownames(compile2) <- compile2$Species
  
  diff <- compile$Species[which(is.na(match(compile$Species, compile2$Species)))]
  diff <- diff[which(is.na(match(diff,EcomorphData$Species)))]
  diff
  
  old <- compile[!is.na(match(compile$Species,diff)),]; old
  rownames(old) <- old$Species
  new <- as.data.frame(criteria.lda[diff,]);new
  criteria2[diff,]
  just <- cbind(old[,c("Species","Predicted")],New.LDA = new[,"Pred.Eco"], C2.Predicted = criteria2[diff,"Pred.Eco"])
  #just["laeviventris","C2.Predicted"] <- "TC"
  #just["cusuco","C2.Predicted"] <- "TC"
  just
  synth <- old[which(just$C2.Predicted == just$New.LDA & just$C2.Predicted == just$Predicted),]
  for (i in 1:nrow(synth)) {
    synth[i,"LDA"] <- new[synth$Species[i],synth$Predicted[i]]
  }
  synth <- synth[which(synth$LDA >= 0.8),]
  trans <-just$Species[which(just$New.LDA == "G" & just$C2.Predicted == just$New.LDA)]
  total <- compilation(lda = criteria.lda, predicted = predicted, hard.mode = F)
  trans <- total[!is.na(match(total$Species,trans)),]
  trans <- trans[which(trans$LDA >= 0.8),]
  rownames(trans) <- trans$Species
  
  new.compile2 <- rbind(compile2,synth,trans)
  new.compile2 <- new.compile2[sort(new.compile2$Species),]
  
  return(new.compile2)
}
new.compile2 <- comp.diff()

compile$Species[which(is.na(match(compile$Species, compile2$Species)))]

tapply(new.compile2$Predicted,new.compile2$Ecomorph,length)
tapply(new.compile2$Predicted[which(new.compile2$Ecomorph == "M")],
       new.compile2$Predicted[which(new.compile2$Ecomorph == "M")],length)
#Pred.ground <- Pred.ground[which(Pred.ground$Ecomorph == "Unknown"),]

Pred.Eco2 <- data.frame(Species = NewData$Species, Pred.Eco = NewData$Ground)
rownames(Pred.Eco2) <- Pred.Eco2$Species
#all(match(Pred.Eco[,1]!is.na(match(Pred.Eco$Species,compile[[1]]$Species,compile[[1]]$Species))
Pred.Eco2[!is.na(match(Pred.Eco2$Species,compile2$Species)),2] <- compile2$Predicted
#Pred.Eco2[criteria3[which(!is.na(criteria3$Pred.Eco)),1],2] <- criteria3$Pred.Eco[which(!is.na(criteria3$Pred.Eco))] #Irshick and Losos 97
#Pred.Eco2[criteria2[which(!is.na(criteria2$Pred.Eco)),1],2] <- criteria2$Pred.Eco[which(!is.na(criteria2$Pred.Eco))] #Irshick and Losos 97


ggplot.pca(pca = phylopca, axis1 = 1, axis2 =3, species = tree$tip.label, 
           groups = Pred.Eco2$Pred.Eco, labels = FALSE)


# Randomization Tests -----------------------------------------------------

Ecomorph.Scores <- cbind("Ecomorph" = NewData$Ground,as.data.frame(phylopca$S)) %>%
  filter(., !Ecomorph == "M" & !Ecomorph == "U")
Ecomorph.Scores <- NewData[,-c(1:2,4:7,21:22)] %>%
  filter(., !Ground == "M" & !Ground == "U")

# perform randomization
random.dist <- function(scores, axes, group1, group2, nsim = 1000) {
  scores <- scores[which(scores[,1] == group1 | scores[,1] == group2),1:(axes+1)] # trims scores matrix to include only relevant species and axes
  group <- scores[,1]

  #calculate actual centroid distance
  group.centroids <- (aggregate(scores[,-1], list(group), mean)) # calculates centroid for each  group
  #rownames(group.centroids) <- group.centroids[,1]
  group.centroids <- group.centroids %>%
    dplyr::select(., 2:ncol(group.centroids)) %>%
    as.data.frame
  dfa <- as.matrix(dist(group.centroids, method = "euclidean")) # calculates euclidean distances between centroids
  dist <- dfa[2,1] #isolate centroid dist
  
  #perform randomization tests and calculate centroid distances
    dummy.dist <- c(1:nsim)
    for (i in 1:nsim) {
      group1.count <- tapply(scores[,1],scores[,1],length)[1] # how many species are in group 1
      dummy.scores <- scores
      sample <- sample(1:nrow(dummy.scores),group1.count)
      dummy.scores[sample,1] <- group1
      dummy.scores[is.na(match(1:nrow(scores),sample)),1] <- group2
      
      dummy.group <- dummy.scores[,1]
      
      dummy.group.centroids <- (aggregate(dummy.scores[,-1], list(dummy.group), mean)) # calculates centroid for each  group
      dummy.group.centroids <- dummy.group.centroids %>%
        dplyr::select(., 2:ncol(dummy.group.centroids)) %>%
        as.data.frame
      dummy.dfa.table <- as.matrix(dist(dummy.group.centroids, method = "euclidean")) # calculates euclidean distances between centroids
      dummy.dist[i] <- dummy.dfa.table[2,1] #isolate centroid dist
    }
  
  pval <- length(which(dist < (dummy.dist)))/nsim

  return(list(dist = dist, sim = dummy.dist, pval = pval))
  
}


# do multiple corrections for either randomization tests or simulated test
posthoc.cross <- function(scores, axes, nsim = 1000, fun, p.adj = "holm") {
  bon<- matrix(NA,length(unique(scores[,1])),length(unique(scores[,1])))
  colnames(bon) <- sort(unique(scores[,1]));rownames(bon) <- sort(unique(scores[,1]))
  dat <- bon
  cross <- crossing(colnames(bon),rownames(bon))
  cross <- cross[!duplicated(t(apply(cross, 1, sort))),]
  cross <- cross[which(!cross[,1] == cross[,2]),]
  cross <- as.data.frame(cross)
  
  if (fun == "random.dist") {
    for (i in 1:nrow(cross)){
      dist <- random.dist(scores = scores, axes = axes, group1 = cross[i,2], group2 = cross[i,1], nsim = nsim)
      dat[cross[i,2],cross[i,1]] <- dist$dist
      bon[cross[i,2],cross[i,1]] <- dist$pval
    }
  }
  
  if (fun == "sim.dist") {
    for (i in 1:nrow(cross)){
      dist <- sim.dist(tree, scores = scores, axes = axes, group1 = cross[i,2], group2 = cross[i,1], nsim = nsim)
      dat[cross[i,2],cross[i,1]] <- dist$dist
      bon[cross[i,2],cross[i,1]] <- dist$pval
    }
  }
  posthoc <- matrix((p.adjust(bon, method = p.adj)),ncol(bon),nrow(bon))
  dimnames(posthoc) <- dimnames(bon)
  return(list("dist" = dat, "Pf" = bon, "Pt" = posthoc))
}
  bon <- posthoc.cross(Ecomorph.Scores, axes = 13, fun = "random.dist", p.adj = "holm", nsim = 1000)

# Simulated Trait Data ----------------------------------------------------

#simulate data for 5 axes under BM and test differences between ecomorph centroids
sim.dist <- function(tree = tree, scores, axes, group1, group2, nsim = 1000) {
  scores <- scores[which(scores[,1] == group1 | scores[,1] == group2),1:(axes+1)] # trims scores matrix to include only relevant species and axes
  group <- scores[,1]
  trials <- nsim
  #calculate actual centroid distance
  group.centroids <- (aggregate(scores[,-1], list(group), mean)) # calculates centroid for each  group
  #rownames(group.centroids) <- group.centroids[,1]
  group.centroids <- group.centroids %>%
    dplyr::select(., 2:ncol(group.centroids)) %>%
    as.data.frame
  dfa <- as.matrix(dist(group.centroids, method = "euclidean")) # calculates euclidean distances between centroids
  dist <- dfa[2,1] #isolate centroid dist
  
  #perform centroid comparison with simulated data
  dummy.dist <- c(1:trials)
  for (i in 1:trials) {
    simdata <- fastBM(tree = tree,nsim = axes) # simulate data under BM
    simdata <- cbind(group,as.data.frame(simdata[rownames(scores),])) # prep simulated data frame
    
    #calculate actual centroid distance
    dummy.group.centroids <- (aggregate(simdata[,-1], list(simdata[,1]), mean)) # calculates centroid for each group
    #rownames(group.centroids) <- group.centroids[,1]
    dummy.group.centroids <- dummy.group.centroids %>%
      dplyr::select(., 2:ncol(dummy.group.centroids)) %>%
      as.data.frame
    dummy.dfa <- as.matrix(dist(dummy.group.centroids, method = "euclidean")) # calculates euclidean distances between centroids
    dummy.dist[i] <- dummy.dfa[2,1] #isolate centroid dist
  }
  
  pval <- length(which(dist < (dummy.dist)))/trials
  
  return(list(dist = dist, sim = dummy.dist, pval = pval))
  
}

# all ecomorph comparisons with multi comparison correction
bon <- posthoc.cross(Ecomorph.Scores, axes = 13, fun = "sim.dist", p.adj = "holm", nsim = 1000)

# simulate data and test compare distance to centroid and NND against null
sim.ED <- function(tree = tree, scores, axes, nsim = 1000) {
  #scores <- scores[which(scores[,1] == group1 | scores[,1] == group2),1:(axes+1)] # trims scores matrix to include only relevant species and axes
  scores <- scores[,1:(axes + 1)]
  groups <- scores[,1]
  trials <- nsim
  #calculate actual centroid distance
  group.centroids <- (aggregate(scores[,-1], list(groups), mean)) # calculates centroid for each  group
  rownames(group.centroids) <- group.centroids[,1]
  colnames(group.centroids) <- colnames(scores)
  df <- rbind(scores,group.centroids)
  df <- df %>%
    dplyr::select(., 2:ncol(df)) %>%
    as.data.frame
  dfa <- as.matrix(dist(df, method = "euclidean")) # calculates euclidean distances between centroids
  matrix <- cbind("Species" = rownames(dfa), 
                  Ecomorph = c(groups,sort(unique(groups))), 
                  as.data.frame(dfa))
  matrix[matrix == 0] <- NA
  rownames(matrix) <- matrix[,1]
  colnames(matrix) <- c("Species", "Ecomorph", as.character(matrix[,1]))
  Centroid.matrix <- dfa[1:nrow(scores),(ncol(dfa)-length(unique(groups))+1):ncol(dfa)]
  
  Ecomorph.dist <- matrix(NA, length(sort(unique(groups))), 4)
  rownames(Ecomorph.dist) <- (sort(unique(groups)))
  colnames(Ecomorph.dist) <- c("Mean Centroid Dist", "Centroid Pval","Mean NDD","NND Pval")
  for (i in 1:nrow(Ecomorph.dist)) {
    class <- sort(unique(groups))[i]
    Ecomorph.dist[class,1] <- mean(Centroid.matrix[which(groups == class),class])
    members <- rownames(Centroid.matrix)[which(groups == class)]
    help <- as.matrix(matrix[members,members])
    Ecomorph.dist[class,3] <- mean(apply(help, 1, FUN = mean, na.rm = TRUE)) #calculates average NND for each ecomorph
    }

  sim.centroid.dist <- matrix(NA, trials, length(sort(unique(groups)))); colnames(sim.centroid.dist) <- sort(unique(groups))
  sim.NND.dist <- matrix(NA, trials, length(sort(unique(groups)))); colnames(sim.NND.dist) <- sort(unique(groups))
  #perform centroid comparison with simulated data
    for (i in 1:trials) {
      simdata <- fastBM(tree = tree,nsim = axes) # simulate data under BM
      simdata <- cbind("Ecomorph" = groups,as.data.frame(simdata[rownames(scores),])) # prep simulated data frame
      
      sim.group.centroids <- (aggregate(simdata[,-1], list(groups), mean)) # calculates centroid for each  group
      rownames(sim.group.centroids) <- sim.group.centroids[,1]
      colnames(sim.group.centroids) <- colnames(simdata)
      sim.df <- rbind(simdata,sim.group.centroids)
      sim.df <- sim.df %>%
        dplyr::select(., 2:ncol(sim.df)) %>%
        as.data.frame
      sim.dfa <- as.matrix(dist(sim.df, method = "euclidean")) # calculates euclidean distances between centroids
      sim.matrix <- cbind("Species" = rownames(sim.dfa), 
                      "Ecomorph" = c(groups,sort(unique(groups))), 
                      as.data.frame(sim.dfa))
      sim.matrix[sim.matrix == 0] <- NA
      rownames(sim.matrix) <- sim.matrix[,1]
      #colnames(sim.matrix) <- c("Species", "Ecomorph", as.character(sim.matrix[,1]))
      sim.Centroid.matrix <- sim.dfa[1:nrow(simdata),(ncol(sim.dfa)-length(unique(groups))+1):ncol(sim.dfa)]
      
      for (y in 1:ncol(sim.centroid.dist)) {
        sim.class <- sort(unique(groups))[y]
        sim.centroid.dist[i,sim.class] <- mean(sim.Centroid.matrix[which(groups == sim.class),sim.class])
        members <- rownames(sim.Centroid.matrix)[which(groups == sim.class)]
        help <- as.matrix(sim.matrix[members,members])
        sim.NND.dist[i,sim.class] <- mean(apply(help, 1, FUN = mean, na.rm = TRUE)) #calculates average NND for each ecomorph
      }
    }
    for (i in 1:nrow(Ecomorph.dist)){
      Ecomorph.dist[i,2] <- 1 - (length(which(Ecomorph.dist[i,1] < (sim.centroid.dist[,i])))/trials)
      Ecomorph.dist[i,4] <- 1 - (length(which(Ecomorph.dist[i,3] < (sim.NND.dist[,i])))/trials)
    }
  
  
  
  
  return(Ecomorph.dist)
  
}
var <- sim.ED(tree, scores = Ecomorph.Scores, axes = 13, nsim = 1000)
p.adjust(var[,4], method = "bonferroni")



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


# PhylogeneticEM ----------------------------------------------------------

library(PhylogeneticEM)
