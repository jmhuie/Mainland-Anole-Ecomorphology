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
  hulls <- ddply(scores, .(Group), function(scores) scores[chull(scores[,axis1], scores[,axis2]), ])
  hulls <- hulls[!hulls$Group=="M" & !hulls$Group=="U",] 
  plot <- ggplot(scores, aes(x = setNames(scores[,axis1],species), y = setNames(scores[,axis2],species), col = Group))+
    geom_polygon(data=hulls, aes(x = hulls[,axis1], y = hulls[,axis2], group=Group, fill = Group), alpha = 0.4, show.legend = FALSE) +# size = 0.75)+
    geom_point(aes(fill=Group, shape=(Region)), colour = "black", size = 1.5) +
    labs(fill = "Ecomorph") +
    scale_shape_manual(values = c(21,23,24)) +
    scale_colour_manual(values = col) +
    scale_fill_manual(values = col,
                      guide = guide_legend(override.aes = aes(shape = 21, color = "black"))) +  
    xlab(paste0("pPC",axis1," (", round(diag(pca$Eval)/sum(pca$Eval)*100, digits = 1)[axis1], "% explained var.)")) +
    ylab(paste0("pPC",axis2," (", round(diag(pca$Eval)/sum(pca$Eval)*100, digits = 1)[axis2], "% explained var.)")) +
    theme(axis.line.x = element_line(size = .5, colour = "black"),
          axis.line.y = element_line(size = .5, colour = "black"),
          panel.background = element_rect(fill = "white"), # bg of the panel 
          plot.background = element_rect(fill = "white"), # bg of the plot
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          legend.background = element_rect(fill = "white"), # get rid of legend bg
          legend.box.background = element_rect(fill = "white"),
          legend.key = element_blank(),
          legend.title = element_blank(),
          text=element_text(colour="black", size = 8, face = "bold"),
          axis.text.x=element_text(colour="black", size = 7),
          axis.text.y=element_text(colour="black", size = 7),
          axis.title.x=element_text(colour="black", size = 8, face = "bold"),
          axis.ticks = element_line(colour = "black", size = .5),
          legend.position = "none",
          plot.margin=unit(c(.5,.5,.5,.5),"cm"))
  if (labels == TRUE) {
    plot + geom_text(aes(label=c(rownames(scores)), col = groups),hjust=0, vjust=0, size = 3)
  } else {
    plot
  }
}

###### predict class Euclidean #####
predict.class <- function(scores, species, groups, hard.mode = T, all.species = F) {
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
  
  # criteria 1 - uses convex hulls to assign classes
  {criteria1 <- data.frame(Species = rownames(scores))
    criteria1[,2] <- groups
    for (i in 1:length(Ecomorphs)) {
      data <- as.matrix(scores[which(groups == Ecomorphs[i]),])
      hv <- convhulln(data, return.non.triangulated.facets = TRUE)
      set.seed(20)
      list <- inhulln(hv, as.matrix(scores))
      criteria1[which(list),i+2] <- 1
    } 
    
    criteria1[,ncol(criteria1)+1] <- groups
    colnames(criteria1) <- c("Species","Ecomorph",Ecomorphs, "Pred.Eco")
    rownames(criteria1) <- criteria1[,1]
    criteria1 <- droplevels(criteria1)
    
    for (i in 1:nrow(criteria1)) {
      class <- which(!is.na(criteria1[i,c(3:(ncol(criteria1)-1))]))
      if (length(class) == 0) {
        criteria1[i,ncol(criteria1)] <- NA
      }
      if (length(class) == 1) {
        criteria1[i,ncol(criteria1)] <- colnames(criteria1[,c(3:(ncol(criteria1)-1))])[class]
      }
      if (length(class) > 1) {
        cent.list <- setNames(Centroid.matrix[i,class],colnames(Centroid.matrix[,class]))
        cent.min <- min(cent.list)
        criteria1[i,ncol(criteria1)] <- names(cent.list[which(cent.list == cent.min)])
      }
    }
    criteria1 = criteria1[which(groups == "M" | groups == "U"),]
  } 
  
  # criteria 2 - if closer to a centroid than the furthest member of that centroid
  # start by creating a matrix that resembles the output
  criteria2 <- data.frame(Species = rownames(scores))
    criteria2$Ecomorph <- groups
    criteria2 <- cbind(criteria2,round(Centroid.matrix,3))
    EcomorphSpecies <- criteria2[-which(groups == "M" | groups == "U"),1:2]
    EcomorphSpecies[,3] <- as.vector(1:nrow(EcomorphSpecies))
    
    # determine the max centroid distance for each ecomorph among ecomorph members
    for (i in 1:nrow(EcomorphSpecies)) {
      taxon <- EcomorphSpecies[i,1]
      Ecomorph <- EcomorphSpecies[i,2]
      EcomorphSpecies[i,3] <- criteria2[as.character(taxon),as.character(Ecomorph)]
    }
    EcomorphSpecies <- droplevels(EcomorphSpecies)
    compare <- as.matrix(tapply(as.numeric(EcomorphSpecies[,3]), INDEX = EcomorphSpecies[,2], FUN = max))
    
    # determine which is ecomorph centroid each species is closest too
    # must be closer that the furthest member of that ecomorph to count
    criteria2$Pred.Eco <- NA
    for (i in 1:nrow(criteria2)) {
      test <- compare
      for (y in 1:nrow(compare)) {
        test[y] <- isTRUE(criteria2[i,2+y] <= compare[y])
      }
      
      if (hard.mode == TRUE) {
        criteria2[i,rownames(compare)][which(test < 1)] <- NA 
        min <- names(which.min(criteria2[i,rownames(compare)][which(test == 1)]))
        if (!is.null(min)) {
          criteria2$Pred.Eco[i] <- min
        }
      } else {
        criteria2$Pred.Eco[i] <- names(which.min(criteria2[i,rownames(compare)]))
      }
    }
    
    # naturally subsets out the a priori ecomorph species 
    if (all.species == T) {
    } else {
      criteria2 <- criteria2[which(groups == "M" | groups == "U"),]
    }
    
    # criteria 3 - looks at NND distances. Must be closer to an ecomorph species than the smallest ecomorph average
    # start with creating a matrix that resembles the output
    criteria3 <- data.frame(Species = rownames(scores))
    criteria3$Ecomorph <- groups
    criteria3 <- cbind(criteria3,round(Centroid.matrix,3))
    EcomorphSpecies <- criteria3[-which(groups == "M" | groups == "U"),1:2]
    EcomorphSpecies[,3] <- as.vector(1:nrow(EcomorphSpecies))
    
    # determine the mean NND for all ecomorphs
    compareNND <-compare
    for (i in 1:length(compareNND)) {
      members <- rownames(matrix)[which(rownames(compare)[i] == matrix[,2])] 
      help <- as.matrix(matrix[members,members])
      compareNND[i] <- min(apply(help, 1, FUN = mean, na.rm = TRUE)) #calculates average NND for each ecomorph
    }
    
    # determined the NND distances for each species being classified
    for (i in 1:nrow(criteria3)) {
      all.NND <- matrix[which(!matrix[,2] == "M" & !matrix[,2] == "U" & !matrix[,2] == "Centroid"),c(1:2,i+2)] 
      for (y in 3:length(colnames(criteria3))) {
        NND <- min(all.NND[which(all.NND[,2] == colnames(criteria3)[y]),][,3],na.rm = TRUE)
        criteria3[i,y] <- round(NND,3)
      }
    }
    criteria4 <- criteria3
    
    # determine whether the NND criteria is satisfied and the nearest ecomorph
    criteria3$Pred.Eco <- NA
    for (i in 1:nrow(criteria3)) {
      test <- compareNND
      for (y in 1:nrow(compareNND)) {
        test[y] <- isTRUE(criteria3[i,2+y] <= min(compareNND))
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
    
    
    # criteria 4 - looks at NND and classifies species if they closer to an eco species than the mean NND for that eco
    # already created a matrix based on criteria 3
    criteria4$Pred.Eco <- NA
    for (i in 1:nrow(criteria4)) {
      test <- compareNND
      for (y in 1:nrow(compareNND)) {
        test[y] <- isTRUE(criteria4[i,2+y] <= min(compareNND[y]))
      }
      
      if (hard.mode == TRUE){
        criteria4[i,rownames(compareNND)][which(test < 1)] <- NA
        min <- names(which.min(criteria4[i,rownames(compareNND)][which(test == 1)]))
        if (!is.null(min)) {
          criteria4$Pred.Eco[i] <- min
        }
      } else {
        criteria4$Pred.Eco[i] <- names(which.min(criteria4[i,rownames(compareNND)]))
      }
    }
    
    # naturally subsets out the a priori ecomorph species 
    if (all.species == TRUE) {
    } else {
      criteria4 <- criteria4[which(criteria4[,2] == "M" | criteria4[,2] == "U"),]
    }
   
  
  comb.compare <- round(cbind(compare,compareNND),3)
  colnames(comb.compare) <- c("Max Centroid Dist", "Mean NDD")
  print(comb.compare)
  
  return(list(criteria1 = criteria1, 
              criteria2 = criteria2,
              criteria3 = criteria3, 
              criteria4 = criteria4,
              compare = comb.compare))
}

##### COMPILE CLASSIFICATION RESULTS #####
compilation <- function(lda = criteria.lda, predicted, upper.cut = 0.95, lower.cut = 0.9,
                        hard.mode = T) {
  compilation <- matrix(0,0,8)
  other <- matrix(0,0,8)
  criteria.lda <- as.matrix(criteria.lda)
  
  compare <- predicted$compare 
  criteria1 <- predicted$criteria1
  criteria2 <- predicted$criteria2
  criteria3 <- predicted$criteria3
  criteria4 <- predicted$criteria4
  
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
    c1.val <- criteria1[which(criteria1$Species == temp.species),temp.predicted]
    if (!is.na(c1.val) == TRUE){
      temp.c1 <- 1 # yes
    } else {
      temp.c1 <- 0 # no
    }
    
    # determines if it satisfies criteria2 for the predicted ecomorph
    c2.val <- criteria2[which(criteria2$Species == temp.species),temp.predicted]
    if (!is.na(c2.val) == TRUE){
      temp.c2 <- 1
    } else {
      temp.c2 <- 0
    }
    
    # determines if it satisfies criteria3 for the predicted ecomorph
    c3.val <- criteria3[which(criteria3$Species == temp.species),temp.predicted]
    if (!is.na(c3.val) == TRUE){
      temp.c3 <- 1
    } else {
      temp.c3 <- 0
    }
    
    # determines if it satisfies criteria3 for the predicted ecomorph
    c4.val <- criteria4[which(criteria4$Species == temp.species),temp.predicted]
    if (!is.na(c4.val) == TRUE){
      temp.c4 <- 1
    } else {
      temp.c4 <- 0
    }
    
    # determines if it satisfies any of the centroid criteria (1 or 2)
    if (temp.c1 | temp.c2 >= 1) {
      temp.centroid <- 1
    } else {
      temp.centroid <- 0
    }
    
    # determines if it satisfies any of the neighbor criteria (3 or 4)
    if (temp.c3 | temp.c4 >= 1) {
      temp.neighbor <- 1
    } else {
      temp.neighbor <- 0
    }
    
    temp.final <- data.frame(temp.species,temp.region,temp.predicted,lda,temp.c1,c2.val,c3.val,c4.val)
    # if LDA post prob > 0.95 and it satisfies any of the criteria then it gets added to matrix
    if (temp.lda == 1 & (temp.centroid + temp.neighbor) >= 1){
      compilation <- rbind(compilation,temp.final)
    }
    # if LDA post prob <0.95 and >0.9 and it satisfies at least one centroid and one neighbor criteria then it gets added to matrix
    if (temp.lda == 2 & (temp.centroid + temp.neighbor) >= 2){
      compilation <- rbind(compilation,temp.final)
    }
    #if LDA post prob <0.9 and it satisfies at least one centroid and one neighbor criteria then it gets added to matrix
    #this is just every species
    if (temp.lda == 1 | temp.lda == 2 | temp.lda == 3){
      other <- rbind(other,temp.final)
    }
  }
  colnames(compilation) <- c("Species","Ecomorph","Predicted","LDA","C1","C2","C3","C4")
  colnames(other) <- c("Species","Ecomorph","Predicted","LDA","C1","C2","C3","C4")
  
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
    scale_x_discrete(labels = c("CG", "GB","G", "TR", "TC", "TG", "Tw")) +
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




