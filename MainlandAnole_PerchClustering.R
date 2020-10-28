library(pvclust)

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

