library(pvclust)

Substrate <- read_excel("Data/MainlandAnole_EcoData.xlsx") 
Substrate$Percentage <- round(apply(Substrate[,4:27],1, FUN = sum,  na.rm = TRUE),1)
Substrate[,4:27] <- Substrate[,4:27]*Substrate$"Substrate N"/100
Substrate <- Substrate %>%
  filter(!is.na(Substrate$'Substrate N')) %>%
  filter(Percentage >= 99) %>%
  dplyr::select(Species:Fern) %>%
  group_by(Species) %>%
  summarise_at(vars("Substrate N":Fern), sum, na.rm = TRUE)
Substrate <- Substrate %>%
  dplyr::select(., !Crevice)
Substrate[,3:24] <- round(Substrate[,3:24]/Substrate$"Substrate N",2)


main <- Substrate[tree$tip.label[which(NewData[tree$tip.label,"Region"]=="Mainland")],] %>%
  filter(.,!is.na(Species))

sub.traits <- read_excel("Data/MainlandAnole_EcoData.xlsx", sheet = 2) %>%
  as.data.frame
rownames(sub.traits) <- sub.traits[,1]
scaled.data <- scale(sub.traits[,-1])
wss <- (nrow(scaled.data)-1)*sum(apply(scaled.data,2,var)) #calculate within groups sum of squares
for (i in 2:15) wss[i] <- sum(kmeans(scaled.data, iter.max = 1000,
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

sub.data <- t(sub.traits[,-1])
prey.fit.p <- pvclust(sub.data, method.hclust="ward.D2",
                      method.dist="euclidean", nboot = 200)
plot(prey.fit.p)
var <- rect.hclust(prey.fit.p$hclust, k=(10), border="red")

comb.sub <- Substrate %>%
  mutate(.,"Stick" =Stick+Root+Deadfall,.keep = "unused") %>%
  mutate(.,"Ground" = Ground+`Leaf Litter`,.keep = "unused") %>%
  mutate(.,"Bushy" = Shrub+Bush+Fern, .keep = "unused") %>%
  mutate(.,"Rock" = Rock +Log+Pipe, .keep = "unused") %>%
  mutate(.,"RockWall" = Post+`Rock face`+`Artificial Wall`, .keep = "unused") %>%
  mutate(.,"Trunk" = Stump+Trunk+Stem, .keep = "unused")

comb.sub$Max <- NA
comb.sub <- relocate(comb.sub, Max, .after = Species)
for (i in 1:nrow(comb.sub)) {
  max <- max(comb.sub[i,4:ncol(comb.sub)])
  if (max > 0.75) {
    comb.sub$Max[i] <- colnames(sort(comb.sub[i,4:ncol(comb.sub)],decreasing = TRUE)[1])
  } else {
    comb.sub$Max[i] <- paste0(colnames(sort(comb.sub[i,4:ncol(comb.sub)],decreasing = TRUE)[1]),"_",
                              colnames(sort(comb.sub[i,4:ncol(comb.sub)],decreasing = TRUE)[2]))
    
    
  }
  
}  

