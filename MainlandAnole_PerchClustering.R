

Substrate <- read_excel("Data/MainlandAnole_EcoData.xlsx") 
Substrate$Percentage <- round(apply(Substrate[,4:27],1, FUN = sum,  na.rm = TRUE),1)
Substrate[,4:27] <- Substrate[,4:27]*Substrate$"Substrate N"/100
Substrate <- Substrate %>%
  filter(!is.na(Substrate$'Substrate N')) %>%
  filter(Percentage >= 99) %>%
  dplyr::select(Species:Fruit) %>%
  group_by(Species) %>%
  summarise_at(vars("Substrate N":Fruit), sum, na.rm = TRUE)
Substrate[,3:26] <- round(Substrate[,3:26]/Substrate$"Substrate N",2)


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
                      method.dist="euclidean", nboot = 100)
plot(prey.fit.p)
var <- rect.hclust(prey.fit.p$hclust, k=(8), border="red")

comb.sub <- Substrate %>%
  mutate(.,"Ground" =Ground+`Leaf Litter`+Stick+Water,.keep = "unused") %>%
  mutate(.,"Leafy_Vegetation" = Leaf+Shrub+Bush+Fern,.keep = "unused") %>%
  mutate(.,"Root" = Root+Deadfall, .keep = "unused") %>%
  mutate(.,"Stem" = Post+Stump+Log+Stem, .keep = "unused") %>%
  mutate(.,"Trunk_Rock" = Trunk+Rock+`Rock face`+Artifical, .keep = "unused") %>%
  mutate(.,"Grass_Vine" = Grass+Vine, .keep = "unused") %>%
  mutate(.,"Branch_Twig" = Branch+Twig, .keep = "unused") %>%
  mutate(.,"Other" = Cervice+Fruit, .keep = "unused")

for (i in 1:nrow(comb.sub)) {
  max <- which(comb.sub[i,] == max(comb.sub[i,3:10]))
  if (length(max) > 1) {
    comb.sub$Max[i] <- paste0(colnames(comb.sub)[max][1],"_",colnames(comb.sub)[max][2])
  } else {
    comb.sub$Max[i] <- colnames(comb.sub)[max]
  }
  
}  
