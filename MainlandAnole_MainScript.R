rm(list = ls(all=TRUE))
source("MainlandAnole_Functions.R")
library(tidyverse)
library(readxl)
library(phytools)
library(plyr)
library(ggpubr)
library(MASS)
library(geiger)
library(l1ou)
library(windex)


# Read Tree ---------------------------------------------------------------

#read in species list with ecomorph and region assignments
EcoData <- read_excel("Data/MainlandAnole_SpeciesList.xlsx")%>%
  as.data.frame
rownames(EcoData) <- EcoData$Species 
EcoData <- EcoData[sort(EcoData$Species),]

tree <- read.nexus("Trees/Poe_2017_MCC_names.txt")
#tree<-chronopl(tree,lambda = 0 ,age.min = c(70,17), age.max = c(72,23), node = c(384,407))
#tree <- read.nexus("Trees/Poe_2017_Time_names.tre")
tree <- ladderize(tree, right = FALSE)
tree <- drop.tip(tree, setdiff(tree$tip.label, unique(EcoData$Species))) #keeps tip labels that matches vector
#if (is.ultrametric(tree) == FALSE) {
#  tree <- force.ultrametric(tree) # forces tree to be ultrametric
#}
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
NewData <- cbind(EcoData[tree$tip.label,-c(4,7:8)],NewData[tree$tip.label,-1])
remove(tail,tail.resid,phyl.resid,MorphoData)

# PCA ---------------------------------------------------------------------

phylopca <- phyl.pca(tree, NewData[tree$tip.label,6:18], method = "BM", mode = "cov")
phylopca$L
plot.phyl.pca(phylopca, type = "lines")

col <- setNames(c("Blue","cyan", "Gold", "grey73", "ForestGreen", "tan4", "Red", "Purple4", "Darkorange1"),sort(unique(NewData$Ground)))
col2 <- setNames(c("Blue", "Gold", "grey73", "ForestGreen", "tan4", "Red", "Purple4", "Darkorange1"),sort(unique(NewData$Ecomorph)))

species <- rownames(phylopca$S)
eco <- setNames(NewData[species,"Ground"], species)
phylopca$S[,c(3)] <-  phylopca$S[,c(3)] *-1
phylopca$S[,c(5)] <-  phylopca$S[,c(5)] *-1
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
Cuba <- mean(unlist(MainNND[which(MainNND$Region == "Cuba"),c(7)]))
Jam <- mean(unlist(MainNND[which(MainNND$Region == "Jamaica"),c(7)]))
Hisp <- mean(unlist(MainNND[which(MainNND$Region == "Hispaniola"),c(7)]))
PR <- mean(unlist(MainNND[which(MainNND$Region == "Puerto Rico"),c(7)]))
Main <- (apply(MainNND[which(MainNND$Region == "Mainland"),c(3,4,5,6)],2,mean))
totalmean <- mean(c(Cuba,Jam,Hisp,PR,Main));totalmean
#totalmean <- mean(c(Cuba,Jam,Hisp,PR,apply(MainNND[,3:6],2,mean)))
#totalmean <- Main



fit1 <- fitContinuous(tree,setNames(phylopca$S[,1],rownames(phylopca$S))[tree$tip.label],model = "BM")
fit2 <- fitContinuous(tree,setNames(phylopca$S[,2],rownames(phylopca$S))[tree$tip.label],model = "BM")
fit3 <- fitContinuous(tree,setNames(phylopca$S[,3],rownames(phylopca$S))[tree$tip.label],model = "BM")
fit4 <- fitContinuous(tree,setNames(phylopca$S[,4],rownames(phylopca$S))[tree$tip.label],model = "BM")
fit5 <- fitContinuous(tree,setNames(phylopca$S[,5],rownames(phylopca$S))[tree$tip.label],model = "BM")
totalsim <-matrix(NA,1000,1)
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
  simCuba <- mean(unlist(simNND[which(simNND$Region == "Cuba"), c(7)]))
  simJam <- mean(unlist(simNND[which(simNND$Region == "Jamaica"), c(7)]))
  simHisp <- mean(unlist(simNND[which(simNND$Region == "Hispaniola"), c(7)]))
  simPR <- mean(unlist(simNND[which(simNND$Region == "Puerto Rico"), c(7)]))
  simMain <- (apply(simNND[which(simNND$Region == "Mainland"), c(3,4,5,6)],2,mean))
  #totalsim[y] <- mean(c(simCuba,simJam,simHisp,simPR))
  totalsim[y] <- mean(c(simCuba,simJam,simHisp,simMain))
}
quantile(c(totalsim),.05)
totalmean
which(!is.na(match(sort(c(totalmean,totalsim)),totalmean)))/1000

hist(c(totalsim,totalmean))
abline(v = totalmean, col = "red", lwd = 2)


#RESULTS - it really matters whether you compare a species to the same species or a random one.
# Former is super sig but the later is really not sig
# I think former is correct 


# L1OU --------------------------------------------------------------------

lizard <- adjust_data(tree,phylopca$S[tree$tip.label,1:5])

fit_ind <- estimate_shift_configuration(lizard$tree,lizard$Y, nCores = 6, criterion = "pBIC")
#saveRDS(fit_ind, "Outputs/shift_config_pBIC_MRCT.rds")
fit_ind <- readRDS("Outputs/shift_config_AIC_MRCT.rds")

fit_conv <- estimate_convergent_regimes(fit_ind, criterion = "AIC", nCores = 4)
#saveRDS(fit_conv, "Outputs/convergent_regime_AIC_MRCT.rds")
fit_conv <- readRDS("Outputs/convergent_regime_AIC_MCC.rds")

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

dev.off()
plot(fit_conv,edge.shift.ann = F,plot.bar = F,  asterisk = T, cex = 0.5)
edgelabels(edge = c(300, 135),cex =0.5,frame = "none")

fit_ind_AIC_bootstrap <- l1ou_bootstrap_support(fit_ind,nItrs = 100, multicore = T, nCores=6)
saveRDS(fit_ind_AIC_bootstrap, "Outputs/shift_bootstrap_AIC_MCC.rds")
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

#regime 1 - omiltemanus & ortonii & coelestinus
#regime 10 - 4 GB regimes
#regime 11 - 1 TC regime
#regime 12 - lionotus
#regime 13 - major Draconura radiation
#regime 2 - 4 twig regimes
#regime 3 - 4 twig regimes
#regime 4 - 1 GB regime
#regime 5 - salvini
#regime 6 - chrysolepis and gracilipes
#regime 7 - big twig and darlingtoni
#regime 8 - two trunk + eugenegrahami
#regime 9 - angusticeps and mainland twig


# DFA w/ Caribbean --------------------------------------------------------

EcomorphData <- NewData[which(!NewData$Ecomorph == "M" & !NewData$Ecomorph == "U"),]
EcomorphData <- droplevels(EcomorphData)

summary(manova(as.matrix(EcomorphData[,c(6:18)])~EcomorphData$Ecomorph), test = "Wilks")

LDA <- lda(EcomorphData$Ecomorph ~ ., 
           data = EcomorphData[,c(6:18)],
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
  match2 <- any(match(str_split(intermediate[i,], "/")[[3]],str_split(intermediate[i,], "/")[[6]][1]))
  if(isTRUE(match) == TRUE & isTRUE(match2) == TRUE) {
    intermediate$Match[i] <- intermediate$DFA[i]
  } else {
    intermediate$Match[i] <- NA
  }
}
intermediate <- filter(intermediate, str_detect(intermediate$DFA, "/"))
intermediate <- intermediate[!is.na(intermediate$Match),c("Species","Ecomorph","DFA","C1","C3","Match")]
intermediate

tapply(intermediate[which(intermediate$Ecomorph == "M"), "Match"],intermediate[which(intermediate$Ecomorph == "M"), "Match"],length)

# this is for later to do the simmap
simmap.eco <- setNames(NewData$Ecomorph,NewData$Species)
simmap.eco[which(simmap.eco == "M")] <- "U"
simmap.eco[compile$Species] <- compile$Predicted
simmap.col <- setNames(c("Blue", "Gold",  "ForestGreen", "tan4","Red", "Purple4","grey73"),sort(unique(simmap.eco)))
simmap.eco<-to.matrix(simmap.eco,levels(as.factor(simmap.eco)))

for (i in 1:nrow(intermediate)) {
  #split <- str_split(intermediate[i,"Match"], "/")
  #simmap.eco[intermediate$Species[i],c(split)[[1]]] <- 0.5
  simmap.eco[intermediate$Species[i],-ncol(simmap.eco)] <- as.numeric(criteria.lda[intermediate$Species[i],3:(ncol(criteria.lda)-1)])
  simmap.eco[intermediate$Species[i],"U"] <- 0
}

# DFA w/ Ground -----------------------------------------------------------

EcomorphData <- NewData[which(!NewData$Ground == "M" & !NewData$Ground == "U"),]
EcomorphData <- droplevels(EcomorphData)

summary(manova(as.matrix(EcomorphData[,c(6:18)])~EcomorphData$Ground), test = "Wilks")

LDA <- lda(EcomorphData$Ground ~ ., 
           data = EcomorphData[,c(6:18)],
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

ggplot.pca(pca = phylopca, axis1 = 1, axis2 =2, species = tree$tip.label, 
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
  match2 <- any(match(str_split(intermediate[i,], "/")[[3]],str_split(intermediate[i,], "/")[[6]][1]))
  if(isTRUE(match) == TRUE & isTRUE(match2) == TRUE) {
    intermediate$Match[i] <- intermediate$DFA[i]
  } else {
    intermediate$Match[i] <- NA
  }
}
intermediate <- filter(intermediate, str_detect(intermediate$DFA, "/"))
intermediate <- intermediate[!is.na(intermediate$Match),c("Species","Ecomorph","DFA","C1","C3","Match")]
intermediate <- intermediate[is.na(match(intermediate$Species,new.compile2$Species)),]
intermediate

tapply(intermediate[which(intermediate$Ecomorph == "M"), "Match"],intermediate[which(intermediate$Ecomorph == "M"), "Match"],length)


# this is for later to do the simmap
simmap.eco2 <- setNames(NewData$Ground,NewData$Species)
simmap.eco2[which(simmap.eco2 == "M")] <- "U"
simmap.eco2[new.compile2$Species] <- new.compile2$Predicted
simmap.col2 <- setNames(c("Blue", "cyan","Gold",  "ForestGreen", "tan4","Red", "Purple4","grey73"),sort(unique(simmap.eco2)))
simmap.eco2<-to.matrix(simmap.eco2,levels(as.factor(simmap.eco2)))

for (i in 1:nrow(intermediate)) {
  #split <- str_split(intermediate[i,"Match"], "/")
  #simmap.eco[intermediate$Species[i],c(split)[[1]]] <- 0.5
  simmap.eco2[intermediate$Species[i],-ncol(simmap.eco2)] <- as.numeric(criteria.lda[intermediate$Species[i],3:(ncol(criteria.lda)-1)])
  simmap.eco2[intermediate$Species[i],"U"] <- 0
}


# Randomization Tests -----------------------------------------------------

Ecomorph.Scores <- cbind("Ecomorph" = NewData$Ground,as.data.frame(phylopca$S)) %>%
  filter(., !Ecomorph == "M" & !Ecomorph == "U")
#Ecomorph.Scores <- NewData[,-c(1:2,4,18:19)] %>%
#  filter(., !Ground == "M" & !Ground == "U")

bon <- posthoc.cross(Ecomorph.Scores, axes = 5, fun = "random.dist", p.adj = "bonferroni", nsim = 1000)

# Simulated Trait Data ----------------------------------------------------


Ecomorph.Scores <- cbind("Ecomorph" = NewData$Ground,as.data.frame(phylopca$S))
# all ecomorph comparisons with multi comparison correction
sim.dist(tree = tree, scores = Ecomorph.Scores, axes =5, group1 = "CG", group2 = "TG", nsim = 1000)
bon2 <- posthoc.cross(Ecomorph.Scores, axes = 5, fun = "sim.dist", p.adj = "bonferroni", nsim = 1000)

# simulate data and test compare distance to centroid and NND against null
var <- sim.ED(tree, scores = Ecomorph.Scores, axes = 5, nsim = 100)
p.adjust(var[,2], method = "bonferroni")


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
hundotrees<-lapply(hundotrees,keep.tip,tree$tip.label)
class(hundotrees) <- "multiPhylo"
Mainland.tree <- drop.tip(tree,setdiff(tree$tip.label,NewData$Species[which(NewData$Region == "Mainland")]))
hundotrees <- lapply(hundotrees, ladderize, right = FALSE)
class(hundotrees) <- "multiPhylo"

#define which set of simmap species and col to use
simeco <- simmap.eco2
simeco.col <- simmap.col2

post.sim<-make.simmap(tree,simmap.eco,model = "SYM",nsim = 1000)
write.simmap(post.sim, "Outputs/SIMMAP_SYM_Caribbean_MRCT.txt", append = TRUE, version = 1.0)
#write.simmap(post.sim1, "Outputs/SIMMAP_SYM_Ground_MCC.txt", append = TRUE, version = 1.0)

#post.sim <- read.simmap("Outputs/SIMMAP_SYM_Ground_postburnin.txt", format="phylip",version =1.0)
#class(post.sim) <- "multiPhylo"

simmap.tree<-lapply(post.sim,drop.tip.simmap,NewData$Species[which(!NewData$Region == "Mainland")])
class(simmap.tree) <- "multiPhylo"
pd<-describe.simmap(simmap.tree, plot = FALSE, ref.tree = Mainland.tree)

transitions <- pd$count
transition.type <- (",CG$")
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

