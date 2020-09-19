##### ggplot pca #####
ggplot.pca <- function(pca, species, groups, region = NewData$Region, axis1, axis2, labels) {
  scores <- as.data.frame(pca$S)
  scores$Species <- species
  scores$Group <- groups
  scores$Region <- region
  hulls <- ddply(scores, .(Group), function(scores) scores[chull(scores[,axis1], scores[,axis2]), ])
  hulls <- hulls[!hulls$Group=="Mainland" & !hulls$Group=="Unknown",] 
  plot <- ggplot(scores, aes(x = setNames(scores[,axis1],species), y = setNames(scores[,axis2],species), col = Group))+
    geom_polygon(data=hulls, aes(x = hulls[,axis1], y = hulls[,axis2], group=Group, fill = Group), alpha = 0.4, show.legend = FALSE) +# size = 0.75)+
    geom_point(aes(fill=Group, shape=(Region)), colour = "black", size = 1.5) +
    labs(fill = "Ecomorph") +
    scale_shape_manual(values = c(21,24)) +
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

##### plot phylomorphspace #####
plot.phylopca <- function(tree, pca, axis1, axis2, group, group.col = rep("black",length(group)), 
                          pch = 20, cex =1, node.size = c(0,0), edge.col = "black" ,label = "off") {
  scores <- pca$S
  rep("black", length(tree$tip.label))
  treeTraits = matrix(c(scores[,axis1], scores[,axis2]), nrow = length(tree$tip.label), ncol =2) 
  row.names(treeTraits) = tree$tip.label
  tips.col <- rep(group.col[group], length(tree$tip.label)) # change depending on classification
  names(tips.col) <- tree$tip.label
  phylo.col<-c(tips.col[tree$tip.label],rep("black",tree$Nnode))
  col.edge<-setNames(rep(edge.col,nrow(tree$edge)),as.character(tree$edge[,2]))
  names(phylo.col)<-1:(length(tree$tip)+tree$Nnode)
  phylomorphospace(tree, treeTraits, label = label, control=list(col.node=phylo.col, col.edge = col.edge), node.size = node.size,
                   xlab="", ylab="", fsize = .75, ftype = "i", lwd = 1)
  if ((pch[1] >= 21) == TRUE) {
    points(scores[,axis1],scores[,axis2],col="black",bg = phylo.col, pch=pch,cex=cex)
  } else {
    points(scores[,axis1],scores[,axis2],col=phylo.col, pch=pch,cex=cex)
  }
  title(xlab=paste0("PC",axis1, " (",round(summary(pca)$importance[2,axis1]*100,1), "% explained var.)"),
        ylab=paste0("PC",axis2, " (",round(summary(pca)$importance[2,axis2]*100,1), "% explained var.)"), cex.lab = 1)
}

###### predict class Euclidean #####
predict.class <- function(scores, species, groups) {
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
  Centroid.matrix <- Centroid.matrix[,which(!colnames(Centroid.matrix) == "Mainland" & !colnames(Centroid.matrix)=="Unknown")]
  
  Ecomorphs <- unique(groups)[-which(unique(groups) == "Mainland" | unique(groups) == "Unknown")]
  Ecomorphs <- sort(Ecomorphs)
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
    criteria1 = criteria1[which(groups == "Mainland" | groups == "Unknown"),]
  } # criteria 1 - uses convex hulls to assign classes
  
  {criteria2 <- data.frame(Species = rownames(scores))
    criteria2$Ecomorph <- groups
    criteria2 <- cbind(criteria2,round(Centroid.matrix,5))
    EcomorphSpecies <- criteria2[-which(groups == "Mainland" | groups == "Unknown"),1:2]
    EcomorphSpecies[,3] <- as.vector(1:nrow(EcomorphSpecies))
    
    for (i in 1:nrow(EcomorphSpecies)) {
      taxon <- EcomorphSpecies[i,1]
      Ecomorph <- EcomorphSpecies[i,2]
      EcomorphSpecies[i,3] <- criteria2[as.character(taxon),as.character(Ecomorph)]
    }
    EcomorphSpecies <- droplevels(EcomorphSpecies)
    compare <- as.matrix(tapply(as.numeric(EcomorphSpecies[,3]), INDEX = EcomorphSpecies[,2], FUN = max))
    
    for(i in 3:ncol(criteria2)) {
      test <- which(!criteria2[,i] <= compare[(i-2),])
      criteria2[test,i] <- NA
    }
    criteria2$Pred.Eco <- NA
    for (i in 1:nrow(criteria2)) {
      min <- names(which.min(criteria2[i,3:(ncol(criteria2)-1)]))
      if (!is.null(min) == TRUE) {
        criteria2$Pred.Eco[i] <- min
      }
    }
    
    criteria2 <- criteria2[which(groups == "Mainland" | groups == "Unknown"),]
  } # criteria 2 - if closer to a centroid than the furthest member of that centroid
  
  {criteria3 <- data.frame(Species = rownames(scores))
    criteria3$Ecomorph <- groups
    criteria3 <- cbind(criteria3,round(Centroid.matrix,5))
    EcomorphSpecies <- criteria3[-which(groups == "Mainland" | groups == "Unknown"),1:2]
    EcomorphSpecies[,3] <- as.vector(1:nrow(EcomorphSpecies))
    
    compareNND <-compare
    for (i in 1:length(compareNND)) {
      members <- rownames(matrix)[which(rownames(compare)[i] == matrix[,2])] 
      help <- as.matrix(matrix[members,members])
      compareNND[i] <- min(apply(help, 1, FUN = mean, na.rm = TRUE)) #calculates average NND for each ecomorph
    }
    
    for (i in 1:nrow(criteria3)) {
      all.NND <- matrix[which(!matrix[,2] == "Mainland" & !matrix[,2] == "Unknown" & !matrix[,2] == "Centroid"),c(1:2,i+2)] 
      for (y in 3:length(colnames(criteria3))) {
        NND <- min(all.NND[which(all.NND[,2] == colnames(criteria3)[y]),][,3],na.rm = TRUE)
        if (NND <= min(compareNND)) {
          criteria3[i,y] <- NND
        } else {
          criteria3[i,y] <- NA
        }
      }
    }
    
    criteria3$Pred.Eco <- NA
    for (i in 1:nrow(criteria3)) {
      min <- names(which.min(criteria3[i,3:(ncol(criteria3)-1)]))
      if (!is.null(min) == TRUE) {
        criteria3$Pred.Eco[i] <- min
      }
    }
    criteria3 <- criteria3[which(criteria3[,2] == "Mainland" | criteria3[,2] == "Unknown"),]
  } # criteria 3 - strict
  
  {criteria4 <- data.frame(Species = rownames(scores))
    criteria4$Ecomorph <- groups
    criteria4 <- cbind(criteria4,round(Centroid.matrix,5))
    EcomorphSpecies <- criteria4[-which(groups == "Mainland" | groups == "Unknown"),1:2]
    EcomorphSpecies[,3] <- as.vector(1:nrow(EcomorphSpecies))
    
    for (i in 1:nrow(criteria4)) {
      all.NND <- matrix[which(!matrix[,2] == "Mainland" & !matrix[,2] == "Unknown" & !matrix[,2] == "Centroid"),c(1:2,i+2)] 
      for (y in 3:length(colnames(criteria4))) {
        NND <- min(all.NND[which(all.NND[,2] == colnames(criteria4)[y]),][,3],na.rm = TRUE)
        if (NND <= compareNND[y-2]) {
          criteria4[i,y] <- NND
        } else {
          criteria4[i,y] <- NA
        }
      }
    }
    
    criteria4$Pred.Eco <- NA
    for (i in 1:nrow(criteria4)) {
      min <- names(which.min(criteria4[i,3:(ncol(criteria4)-1)]))
      if (!is.null(min) == TRUE) {
        criteria4$Pred.Eco[i] <- min
      }
    }
    criteria4 <- criteria4[which(criteria4[,2] == "Mainland" | criteria4[,2] == "Unknown"),]
  } # criteria 4 - liberal 
  
  comb.compare <- (cbind(compare,compareNND))
  colnames(comb.compare) <- c("Max Centroid Dist", "Mean NDD")
  print(comb.compare)
  
  return(list(criteria1 = criteria1, 
              criteria2 = criteria2,
              criteria3 = criteria3, 
              criteria4 = criteria4,
              Pred.Eco = cbind(criteria1[,1:2],criteria1$Pred.Eco,criteria2$Pred.Eco,criteria3$Pred.Eco,criteria4$Pred.Eco)))
}

##### COMPILE CLASSIFICATION RESULTS #####
compilation <- function(lda = criteria.lda, c1=criteria1, c2 = criteria2, c3 = criteria3, c4 = criteria4) {
  compilation <- matrix(0,0,8)
  other <- matrix(0,0,8)
  criteria.lda <- as.matrix(criteria.lda)
  for (i in 1:nrow(criteria.lda)) {
    temp.species <- criteria.lda[i,1] #reads species name
    temp.region <- criteria.lda[i,2] #reads region (mainland or unknown)
    temp.predicted <- criteria.lda[i,ncol(criteria.lda)] # reads LDA predicted class
    if (!is.na(temp.predicted)){
      lda <- criteria.lda[i,temp.predicted] # reads LDA posterior probability
      # determines if the LDA post prob is > 0.95, > 0.9, or < 0.90
      if (lda >= 0.95) { 
        temp.lda <- 1
      } 
      if (lda >= 0.90 & lda < 0.95) {
        temp.lda <- 2
      } 
      if (lda < 0.9){
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
    if (temp.lda == 1 | temp.lda == 2 | temp.lda == 3){
      other <- rbind(other,temp.final)
    }
  }
  colnames(compilation) <- c("Species","Ecomorph","Predicted","LDA","C1","C2","C3","C4")
  colnames(other) <- c("Species","Ecomorph","Predicted","LDA","C1","C2","C3","C4")
  return(list(compilation, other))
}


###### ground compile #####
winnow <- function(compile.eco = compile, compile.eco.ground = compile2, lda = criteria.lda, c1=criteria1, c2 = criteria2, 
                   c3 = criteria3, c4 = criteria4) {
  # species that went from class with 6 to unclass with 7
  main.species <- compile[[1]]$Species[which(compile[[1]]$Ecomorph == "Mainland")][which(is.na(match(compile[[1]]$Species[which(compile[[1]]$Ecomorph == "Mainland")],compile2[[1]]$Species[which(compile2[[1]]$Ecomorph == "Mainland")])))]
  
    a <- setNames(c(1,1,1,1,1),c("C1","C2","C3","C4","Rank"))
    b <- setNames(c(1,1,1,0,2),c("C1","C2","C3","C4","Rank"))
    c <- setNames(c(1,1,0,1,4),c("C1","C2","C3","C4","Rank"))
    d <- setNames(c(1,0,1,1,2),c("C1","C2","C3","C4","Rank"))
    e <- setNames(c(1,0,1,0,3),c("C1","C2","C3","C4","Rank"))
    f <- setNames(c(1,0,0,1,5),c("C1","C2","C3","C4","Rank"))
    g <- setNames(c(1,1,0,0,5),c("C1","C2","C3","C4","Rank"))
    h <- setNames(c(1,0,0,0,7),c("C1","C2","C3","C4","Rank"))
    i <- setNames(c(0,1,1,1,4),c("C1","C2","C3","C4","Rank"))
    j <- setNames(c(0,1,1,0,5),c("C1","C2","C3","C4","Rank"))
    k <- setNames(c(0,1,0,1,6),c("C1","C2","C3","C4","Rank"))
    l <- setNames(c(0,1,0,0,8),c("C1","C2","C3","C4","Rank"))
    m <- setNames(c(0,0,1,1,5),c("C1","C2","C3","C4","Rank"))
    n <- setNames(c(0,0,1,0,7),c("C1","C2","C3","C4","Rank"))
    o <- setNames(c(0,0,0,1,8),c("C1","C2","C3","C4","Rank"))
    p <- setNames(c(0,0,0,0,9),c("C1","C2","C3","C4","Rank"))
    priority <- as.data.frame(rbind(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p))

  
  Pred.Eco.Ground <- compile.eco.ground[[1]]
  for (i in 1:length(main.species)){
      compile.eco1.ED <- compile.eco[[1]][which(compile[[1]]$Species == as.character(main.species[i])),5:8]
      compile.eco1.ED[which(is.na(compile.eco1.ED))] <- 0
      compile.eco1.ED[which(compile.eco1.ED > 0)] <- 1
       if (is.na(criteria1[as.character(main.species[i]),"Ground"])){
        Cr1 <- 0
      } else {
        Cr1 <- 1
      }
      if (is.na(criteria2[as.character(main.species[i]),"Ground"])){
        Cr2 <- 0
      } else {
        Cr2 <- 1
      }
      if (is.na(criteria3[as.character(main.species[i]),"Ground"])){
        Cr3 <- 0
      } else {
        Cr3 <- 1
      }
      if (is.na(criteria4[as.character(main.species[i]),"Ground"])){
        Cr4 <- 0
      } else {
        Cr4 <- 1
      }
      compile.eco.ground.ED <- data.frame(Cr1,Cr2,Cr3,Cr4)
      colnames(compile.eco.ground.ED) <- colnames(compile.eco1.ED)
      #determines support rank of classification with 6 ecomorphs
      for (p in 1:nrow(priority)) {
        if (all(compile.eco1.ED == priority[p,1:4])) {
          rank1 <- priority$Rank[p]
        }
      }
      #determines support rank of classification with 7 ecomorphs
      for (p in 1:nrow(priority)) {
        if (all(compile.eco.ground.ED == priority[p,1:4])){
          rank2 <- priority$Rank[p]
        }
      }
      
      if (rank1 < rank2){
        starter <- as.data.frame(compile.eco[[1]][which(compile.eco[[1]]$Species == as.character(main.species[i])),1:3])
        finisher <- compile.eco[[1]][which(compile[[1]]$Species == as.character(main.species[i])),5:8]
        Pred.Eco.Ground <- rbind(Pred.Eco.Ground,as.data.frame(c(starter, LDA = criteria.lda[as.character(starter$Species),as.character(starter$Predicted)],finisher)))
      }
      if (rank1 > rank2){
        starter <- as.data.frame(compile.eco[[1]][which(compile.eco[[1]]$Species == as.character(main.species[i])),1:2])
        finisher <- as.matrix(c(criteria1[as.character(main.species[i]),"Ground"],criteria2[as.character(main.species[i]),"Ground"],criteria3[as.character(main.species[i]),"Ground"],criteria4[as.character(main.species[i]),"Ground"]))
        finisher <- as.data.frame(t(finisher))
        colnames(finisher) <- c("C1", "C2", "C3","C4")
        Pred.Eco.Ground <- rbind(Pred.Eco.Ground,as.data.frame(c(starter, Predicted = "Ground", LDA = criteria.lda[as.character(starter$Species),"Ground"],finisher)))
      }
  }

  Pred.Eco.Ground <- Pred.Eco.Ground[order(as.character(Pred.Eco.Ground$Species)),]
  rownames(Pred.Eco.Ground) <- 1:nrow(Pred.Eco.Ground)
  return(Pred.Eco.Ground)
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


### consolidate ####
consolidate <- function(compile.eco = compile, compile.eco.ground = compile2, lda = criteria.lda, c1=criteria1, c2 = criteria2, 
                        c3 = criteria3, c4 = criteria4, Ecomorphspp = EcomorphData$Species) {
  conflict <- compile.eco[[1]]$Species[which(is.na(match(compile.eco[[1]]$Species,compile.eco.ground[[1]]$Species)))]
  conflict <- conflict[which(is.na(match(conflict,Ecomorphspp)))]
  Pred.Eco.Ground <- compile.eco.ground[[1]]
  for (i in 1:length(conflict)) {
    name <- conflict[i]
    lda.max <- max(lda[as.character(name),3:9])
    lda.col <-which(lda[as.character(name),3:9] == max(lda[as.character(name),3:9]))
    lda.eco <- colnames(lda)[lda.col+2]
    if (all(is.na(c2[as.character(name),3:9]))) {
      c2.eco <- NA
    } else {
      c2.min <- min(c2[as.character(name),3:9],na.rm = TRUE)
      c2.col <- which(c2[as.character(name),3:9] == min(c2[as.character(name),3:9], na.rm = TRUE))
      c2.eco <- colnames(c2)[c2.col+2]
    }
    if(is.na(c1[(name),lda.col])){
      c1.val <- 0
    } else {
      c1.val <- 1
    }
    c3.val <- c3[as.character(name),lda.col+2]
    c4.val <- c4[as.character(name),lda.col+2]
    if (is.na(match(lda.eco,c2.eco)) == TRUE) {
      row <- cbind(as.character(name), as.character(lda[(name),2]), as.character(lda.eco), as.numeric(lda.max),
                   as.numeric(c1.val), as.numeric(c2.min), (c3.val), (c4.val))
      colnames(row) <- colnames(compile.eco.ground[[1]])
      Pred.Eco.Ground <- rbind(Pred.Eco.Ground,row)
    }
  }
  match(sort(levels(Pred.Eco.Ground$Species)), Pred.Eco.Ground$Species)
  Pred.Eco.Ground <- Pred.Eco.Ground[match(sort(levels(Pred.Eco.Ground$Species)), Pred.Eco.Ground$Species),]
  return(Pred.Eco.Ground)
  
}
