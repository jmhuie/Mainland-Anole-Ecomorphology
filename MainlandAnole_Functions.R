##### phylo scree plot #####
plot.phyl.pca<- function(x,...){
  if(hasArg(main)) main<-list(...)$main
  else main="screeplot"
  x$sdev<-sqrt(diag(x$Eval))
  screeplot(x,main=main)
}

##### ggplot pca #####
ggplot.pca <- function(pca, species, groups, region = NewData$Region, axis1, axis2, labels) {
  scores <- as.data.frame(pca$S)
  scores$Species <- species
  scores$Group <- groups
  scores$Region <- region
  scaleFUN <- function(x) sprintf("%.1f", x)
  hulls <- ddply(scores, .(Group), function(scores) scores[chull(scores[,axis1], scores[,axis2]), ])
  hulls <- hulls[!hulls$Group=="M" & !hulls$Group=="U",] 
  plot <- ggplot(scores, aes(x = scores[,axis1], y = scores[,axis2], col = Group))+
    geom_polygon(data=hulls, aes(x = hulls[,axis1], y = hulls[,axis2], group=Group, fill = Group), alpha = 0.4, show.legend = FALSE) +# size = 0.75)+
    geom_point(aes(fill=Group, shape=(Region)), colour = "black", size = 1.5) +
    labs(fill = "Ecomorph") +
    scale_shape_manual(values = setNames(c(21,24,22),unique(NewData$Region))) +
    scale_colour_manual(values = col) +
    scale_fill_manual(values = col,
                      guide = guide_legend(override.aes = aes(shape = 21, color = "black"))) +  
    xlab(paste0("pPC",axis1," (", round(diag(pca$Eval)/sum(pca$Eval)*100, digits = 1)[axis1], "% explained var.)")) +
    ylab(paste0("pPC",axis2," (", round(diag(pca$Eval)/sum(pca$Eval)*100, digits = 1)[axis2], "% explained var.)")) +
    theme_classic2()+ 
    theme(legend.position = "None",
          axis.title.x=element_text(size = 8, face = "bold"),
          axis.title.y=element_text(size = 8, face = "bold"),
          axis.text.x=element_text(size = 7,face = "bold"),
          axis.text.y=element_text(size = 7,face = "bold"),
          plot.margin=unit(c(.5,.5,.5,.5),"cm"))+
    scale_y_continuous(labels=scaleFUN) +
    scale_x_continuous(labels=scaleFUN) 
  if (labels == TRUE) {
    plot <- plot + geom_text(aes(label=c(rownames(scores)), col = groups),hjust=0, vjust=0, size = 3)
  } else {
    plot
  }
  return(plot)
}

###### Euclidean Distance Predict Class #####
ED.predict <- function(scores, species, groups, hard.mode = T, all.species = F) {
  # scores <- pca scores
  # creates a matrix of Euclidean distancs
  df <- rbind(scores, 
              (aggregate(scores, list(groups), mean))[2:(ncol(scores)+1)]) 
  dfa <- as.matrix(dist(df, method = "euclidean")) # calculates euclidean distances
  matrix <- cbind(Species = c(as.character(species),sort(unique(groups))), 
                  Ecomorph = c(as.character(groups),rep("Centroid",length(unique(groups)))), 
                  as.data.frame(dfa))
  matrix[matrix == 0] <- NA
  rownames(matrix) <- matrix[,1]
  colnames(matrix) <- c("Species", "Ecomorph", as.character(matrix[,1]))
  Centroid.matrix <- matrix[1:nrow(scores),(nrow(scores)+3):ncol(matrix)] # matrix of each species and dist to each centroid
  Centroid.matrix <- Centroid.matrix[,which(!colnames(Centroid.matrix) == "M" & !colnames(Centroid.matrix)=="U")]
  
  Ecomorphs <- unique(groups)[-which(unique(groups) == "M" | unique(groups) == "U")]
  Ecomorphs <- sort(Ecomorphs)
  
  # CRITERIA 1 - if closer to a centroid than the furthest member of that centroid
  # start by creating a matrix that resembles the output
  criteria1 <- data.frame(Species = rownames(scores))
  criteria1$Ecomorph <- groups
  criteria1 <- cbind(criteria1,round(Centroid.matrix,3))
  EcomorphSpecies <- criteria1[-which(groups == "M" | groups == "U"),1:2]
  EcomorphSpecies[,3] <- as.vector(1:nrow(EcomorphSpecies))
  
  # determine the max centroid distance for each ecomorph among ecomorph members
  for (i in 1:nrow(EcomorphSpecies)) {
    taxon <- EcomorphSpecies[i,1]
    Ecomorph <- EcomorphSpecies[i,2]
    EcomorphSpecies[i,3] <- criteria1[as.character(taxon),as.character(Ecomorph)]
  }
  EcomorphSpecies <- droplevels(EcomorphSpecies)
  compare <- as.matrix(tapply(as.numeric(EcomorphSpecies[,3]), INDEX = EcomorphSpecies[,2], FUN = max))
  
  # determine which is ecomorph centroid each species is closest too
  # must be closer that the furthest member of that ecomorph to count
  criteria1$Pred.Eco <- NA
  for (i in 1:nrow(criteria1)) {
    test <- compare
    for (y in 1:nrow(compare)) {
      test[y] <- isTRUE(criteria1[i,2+y] <= compare[y])
    }
    
    if (hard.mode == TRUE) {
      criteria1[i,rownames(compare)][which(test < 1)] <- NA 
      min <- names(which.min(criteria1[i,rownames(compare)][which(test == 1)]))
      if (!is.null(min)) {
        criteria1$Pred.Eco[i] <- min
      }
    } else {
      criteria1$Pred.Eco[i] <- names(which.min(criteria1[i,rownames(compare)]))
    }
  }
  
  # naturally subsets out the a priori ecomorph species 
  if (all.species == T) {
  } else {
    criteria1 <- criteria1[which(groups == "M" | groups == "U"),]
  }
  
  # CRITERIA 2 - UNKNOWN MPD </= MAX MPD
  criteria2 <- data.frame(Species = rownames(scores))
  criteria2$Ecomorph <- groups
  criteria2 <- cbind(criteria2,round(Centroid.matrix,3))
  EcomorphSpecies <- criteria2[-which(groups == "M" | groups == "U"),1:2]
  EcomorphSpecies[,3] <- as.vector(1:nrow(EcomorphSpecies))
  criteria3 <- criteria2
  
  # determine the mean pairwise distance for all ecomorphs
  compareMPD <-compare
  compareNND <-compare
  for (i in 1:length(compareMPD)) {
    members <- rownames(matrix)[which(rownames(compare)[i] == matrix[,2])] 
    help <- as.matrix(matrix[members,members])
    compareMPD[i] <- max(apply(help, 1, FUN = mean, na.rm = TRUE)) #calculates average NND for each ecomorph
    numb <- c()
    for (y in 1:nrow(help)) {
      numb <- c(numb,min(help[y,], na.rm = T))
    }
    compareNND[i] <- max(numb)
  }
  
  # determined the MPD distances for each species being classified
  for (i in 1:nrow(criteria2)) {
    all.MPD <- matrix[which(!matrix[,2] == "M" & !matrix[,2] == "U" & !matrix[,2] == "Centroid"),c(1:2,i+2)] 
    for (y in 3:length(colnames(criteria2))) {
      MPD <- mean(all.MPD[which(all.MPD[,2] == colnames(criteria2)[y]),][,3],na.rm = TRUE)
      criteria2[i,y] <- round(MPD,3)
    }
  }
  
  # determine whether the MPD criteria is satisfied and the nearest ecomorph
  criteria2$Pred.Eco <- NA
  for (i in 1:nrow(criteria2)) {
    test <- compareMPD
    for (y in 1:nrow(compareMPD)) {
      test[y] <- isTRUE(criteria2[i,2+y] <= compareMPD[y])
    }
    
    if (hard.mode == T){
      criteria2[i,rownames(compareMPD)][which(test < 1)] <- NA
      min <- names(which.min(criteria2[i,rownames(compareMPD)][which(test == 1)]))
      if (!is.null(min)) {
        criteria2$Pred.Eco[i] <- min
      }
    } else {
      criteria2$Pred.Eco[i] <- names(which.min(criteria2[i,rownames(compareMPD)]))
    }
  }
  
  # naturally subsets out the a priori ecomorph species 
  if (all.species == TRUE) {
  } else {
    criteria2 <- criteria2[which(criteria2[,2] == "M" | criteria2[,2] == "U"),]
  }
  
  
  # CRITERIA 3 - UNKNOWN NND </= MAX NND
  # matrix already made
  for (i in 1:nrow(criteria3)) {
    all.NND <- matrix[which(!matrix[,2] == "M" & !matrix[,2] == "U" & !matrix[,2] == "Centroid"),c(1:2,i+2)] 
    for (y in 3:length(colnames(criteria3))) {
      NND <- min(all.NND[which(all.NND[,2] == colnames(criteria3)[y]),][,3],na.rm = TRUE)
      criteria3[i,y] <- round(NND,3)
    }
  }
  
  # determine whether the NND criteria is satisfied and the nearest ecomorph
  criteria3$Pred.Eco <- NA
  for (i in 1:nrow(criteria3)) {
    test <- compareNND
    for (y in 1:nrow(compareNND)) {
      test[y] <- isTRUE(criteria3[i,2+y] <= compareNND[y])
    }
    
    if (hard.mode == T){
      criteria3[i,rownames(compareNND)][which(test < 1)] <- NA
      min <- names(which.min(criteria3[i,rownames(compareNND)][which(test == 1)]))
      if (!is.null(min)) {
        criteria3$Pred.Eco[i] <- min
      }
    } else {
      criteria3$Pred.Eco[i] <- names(which.min(criteria3[i,rownames(compareNND)]))
    }
  }
  
  # naturally subsets out the a priori ecomorph species 
  if (all.species == TRUE) {
  } else {
    criteria3 <- criteria3[which(criteria3[,2] == "M" | criteria3[,2] == "U"),]
  }
  
  
  comb.compare <- round(cbind(compare,compareMPD, compareNND),3)
  colnames(comb.compare) <- c("Max CD", "Max MPD", "Max NND")
  print(comb.compare)
  
  return(list(criteria1 = criteria1, 
              criteria2 = criteria2,
              criteria3 = criteria3,
              compare = comb.compare))
}

###### Synthetic Compile #####
synth.compile <- function(lda = criteria.lda, predicted, hard.mode = T, upper.cut = 0.95, lower.cut = 0.9) {
  compilation <- matrix(0,0,7)
  other <- matrix(0,0,6)
  criteria.lda <- as.matrix(criteria.lda)
  
  compare <- predicted$compare 
  criteria1 <- predicted$criteria1
  criteria2 <- predicted$criteria2
  criteria3 <- predicted$criteria3
  
  for (i in 1:nrow(criteria.lda)) {
    temp.species <- criteria.lda[i,1] #reads species name
    temp.region <- criteria.lda[i,2] #reads region (mainland or unknown)
    temp.predicted <- criteria.lda[i,ncol(criteria.lda)] # reads LDA predicted class
    
    # looks at DFA probabilities
    if (!is.na(temp.predicted)){
      lda <- criteria.lda[i,temp.predicted] # reads LDA posterior probability
      # determines if the LDA post prob is > 0.95, > 0.9, or < 0.90
      if (lda >= upper.cut) { 
        temp.lda <- 1
      } 
      if (lda >= lower.cut & lda < upper.cut) {
        temp.lda <- 2
      } 
      if (lda < lower.cut){
        temp.lda <- 3
      }
    }
    
    # determines if it satisfies criteria1 for the predicted ecomorph
    c1.val <- criteria1[which(criteria2$Species == temp.species),temp.predicted]
    if (!is.na(c1.val) == TRUE){
      temp.centroid <- 1
    } else {
      temp.centroid <- 0
    }
    
    # determines if it satisfies criteria2 for the predicted ecomorph
    c2.val <- criteria2[which(criteria2$Species == temp.species),temp.predicted]
    if (!is.na(c2.val) == TRUE){
      temp.pairwise <- 1
    } else {
      temp.pairwise <- 0
    }
    
    c3.val <- criteria3[which(criteria3$Species == temp.species),temp.predicted]
    if (!is.na(c3.val) == TRUE){
      temp.neighbor <- 1
    } else {
      temp.neighbor <- 0
    }
    
    temp.total <- temp.centroid + temp.pairwise + temp.neighbor
    
    temp.final <- data.frame(temp.species,temp.region,temp.predicted,lda,c1.val,c2.val,c3.val)
    # if LDA post prob > 0.95 and it satisfies any of the criteria then it gets added to matrix
    if (temp.lda == 1 & temp.total >= 1){
      compilation <- rbind(compilation,temp.final)
    } 
    if (temp.lda == 1 & temp.total < 1){
      other <- rbind(other,temp.final)
    }
    # if LDA post prob <0.95 and >0.9 and it satisfies at least one centroid and one neighbor criteria then it gets added to matrix
    if (temp.lda == 2 & temp.total >= 2){
      compilation <- rbind(compilation,temp.final)
    } 
    if (temp.lda == 2 & temp.total < 2){
      other <- rbind(other,temp.final)
    }
    #if LDA post prob <0.9 and it satisfies at least one centroid and one neighbor criteria then it gets added to matrix
    #this is just every species
    if (temp.lda == 3){
      other <- rbind(other,temp.final)
    }
  }
  colnames(compilation) <- c("Species","Ecomorph","Predicted","LDA","C1","C2","C3")
  rownames(compilation) <- compilation$Species
  colnames(other) <- c("Species","Ecomorph","Predicted","LDA","C1","C2","C3")
  rownames(other) <- other$Species
  
  if (hard.mode == F) {
    final <- other
  } else {
    final <- compilation
  }
  return(final)
}

###### BOXPLOT #######
ggbox <- function(data = EcomorphData, group = "Ground", trait = "SVL", ylab = "Snout-vent Length") {
  ggplot(data,aes(x = data[,group], y = data[,trait], fill = data[,group])) +
    geom_boxplot(color = "black", size = .5) +
    scale_fill_manual(values=col) +
    xlab("") +
    ylab(ylab) +
    scale_x_discrete(labels = c("CG", "G","GB","TC", "TG","Tr", "Tw")) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    theme(axis.line.x = element_line(size = .65, colour = "black"),
          axis.line.y = element_line(size = .65, colour = "black"),
          panel.background = element_rect(fill = "white"), # bg of the panel 
          plot.background = element_rect(fill = "white"), # bg of the plot
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          legend.background = element_rect(fill = "white"), # get rid of legend bg
          legend.box.background = element_rect(fill = "white"),
          text=element_text(colour="black", size = 10, face = "bold"),
          axis.text.x=element_text(colour="black", size = 10),
          axis.text.y=element_text(colour="black", size = 8),
          # axis.title.x=element_text(colour="black", size = 12),
          axis.ticks = element_line(colour = "black", size = 0.8),
          axis.ticks.x = element_blank(),
          legend.position = "none")
}

###### RANDOMIZATION TESTS #######
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

###### SIMULATE DIST BETWEEN CENTROIDS TESTS #######
sim.dist <- function(tree = tree, scores, axes, group1, group2, nsim = 1000) {
  allscores <- scores
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
  fit1 <- fitContinuous(tree,setNames(allscores[,"PC1"],rownames(allscores)),model = "BM")
  fit2 <- fitContinuous(tree,setNames(allscores[,"PC2"],rownames(allscores)),model = "BM")
  fit3 <- fitContinuous(tree,setNames(allscores[,"PC3"],rownames(allscores)),model = "BM")
  fit4 <- fitContinuous(tree,setNames(allscores[,"PC4"],rownames(allscores)),model = "BM")
  fit5 <- fitContinuous(tree,setNames(allscores[,"PC5"],rownames(allscores)),model = "BM")
  for (i in 1:trials) {
    simdata <- cbind(fastBM(tree, sig2 = fit1$opt$sigsq, nsim = 1),
                     fastBM(tree, sig2 = fit2$opt$sigsq, nsim = 1),
                     fastBM(tree, sig2 = fit3$opt$sigsq, nsim = 1),
                     fastBM(tree, sig2 = fit4$opt$sigsq, nsim = 1),
                     fastBM(tree, sig2 = fit5$opt$sigsq, nsim = 1)) # simulate data under BM
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

###### POST HOC FOR RANDOMIZATION OR SIMULATION TESTS #######
posthoc.cross <- function(scores, axes, nsim = 1000, fun, p.adj = "holm") {
  allscores <- scores
  eco.scores <- allscores %>%
    filter(., !Ecomorph == "M" & !Ecomorph == "U")
  bon<- matrix(NA,length(unique(eco.scores[,1])),length(unique(eco.scores[,1])))
  colnames(bon) <- sort(unique(eco.scores[,1]));rownames(bon) <- sort(unique(eco.scores[,1]))
  dat <- bon
  cross <- tidyr::crossing(colnames(bon),rownames(bon))
  cross <- cross[!duplicated(t(apply(cross, 1, sort))),]
  cross <- cross[which(!cross[,1] == cross[,2]),]
  cross <- as.data.frame(cross)
  
  if (fun == "random.dist") {
    for (i in 1:nrow(cross)){
      dist <- random.dist(scores = eco.scores, axes = axes, group1 = cross[i,2], group2 = cross[i,1], nsim = nsim)
      dat[cross[i,2],cross[i,1]] <- dist$dist
      bon[cross[i,2],cross[i,1]] <- dist$pval
    }
  }
  
  if (fun == "sim.dist") {
    for (i in 1:nrow(cross)){
      dist <- sim.dist(tree, scores = allscores, axes = axes, group1 = cross[i,2], group2 = cross[i,1], nsim = nsim)
      dat[cross[i,2],cross[i,1]] <- dist$dist
      bon[cross[i,2],cross[i,1]] <- dist$pval
    }
  }
  posthoc <- matrix((p.adjust(bon, method = p.adj)),ncol(bon),nrow(bon))
  dimnames(posthoc) <- dimnames(bon)
  return(list("dist" = dat, "Pf" = bon, "Pt" = posthoc))
}

###### SIMULATE DIST TO CENTROID AND NND #######
sim.ED <- function(tree = tree, scores, axes, nsim = 1000) {
  #scores <- scores[which(scores[,1] == group1 | scores[,1] == group2),1:(axes+1)] # trims scores matrix to include only relevant species and axes
  allscores <- scores
  scores <- scores %>% filter(!Ecomorph == "M" & !Ecomorph =="U")
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
  fit1 <- fitContinuous(tree,setNames(allscores[,"PC1"],rownames(allscores)),model = "BM")
  fit2 <- fitContinuous(tree,setNames(allscores[,"PC2"],rownames(allscores)),model = "BM")
  fit3 <- fitContinuous(tree,setNames(allscores[,"PC3"],rownames(allscores)),model = "BM")
  fit4 <- fitContinuous(tree,setNames(allscores[,"PC4"],rownames(allscores)),model = "BM")
  fit5 <- fitContinuous(tree,setNames(allscores[,"PC5"],rownames(allscores)),model = "BM")
  for (i in 1:trials) {
    simdata <- cbind(fastBM(tree, sig2 = fit1$opt$sigsq, nsim = 1),
                     fastBM(tree, sig2 = fit2$opt$sigsq, nsim = 1),
                     fastBM(tree, sig2 = fit3$opt$sigsq, nsim = 1),
                     fastBM(tree, sig2 = fit4$opt$sigsq, nsim = 1),
                     fastBM(tree, sig2 = fit5$opt$sigsq, nsim = 1)) # simulate data under BM
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

