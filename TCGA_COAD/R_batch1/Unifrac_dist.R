# UniFrac Distance
library(ape)
library(geiger)

UniFrac <- function (datP, tree) {
  # Calculate Weighted/Unweighted UniFrac distances
  #	
  # Args:
  #		datP: immune cell composition data after merge 
  #         the propotion to the node
  #		tree: rooted phylogenetic tree of R class "phylo"
  #
  # Returns:
  # 	unifracs: three dimensional array containing the 
  #             weight/unweighted UniFrac distance 
  #
  
  if (!is.rooted(tree)) stop("Rooted phylogenetic tree required!")
  
  # Convert into proportions
  datP <- as.matrix(datP)
  n <- nrow(datP)
  
  # Construct the returning array
  if (is.null(rownames(datP))) {
    rownames(datP) <- paste("comm", 1:n, sep="_")
  }
  # d_UW: unweighted UniFrac, d_VAW: weighted UniFrac
  dimname3 <- c("d_UW", "d_W")
  unifracs <- array(0, c(n, n, 2),
                    dimnames=list(rownames(datP), rownames(datP), dimname3))
  
  datP <- datP[, c(tree$tip.label, tree$node.label)]
  br.len <- tree$edge.length
  edge <- tree$edge
  
  # caculate propotion of immune cell descending from brunch l
  branch_P <- matrix(0, ncol =  nrow(edge),nrow = n) 
  # for(l in 1:nrow(edge)){
  #   branch_P[, l] <- datP[, edge[l, 1]] - datP[, edge[l,2]] # wrong version
  # }
  for(l in 1:nrow(edge)){
    branch_P[, l] <- datP[, edge[l,2]]
  }
  
  # Calculate various UniFrac distances
  for (i in 2:n) {
    for (j in 1:(i-1)) {
      tmp1 <- branch_P[i, ] # sample i
      tmp2 <- branch_P[j, ] # sample j
      if(any(tmp1 != tmp2)){
      ind <- (tmp1 + tmp2) != 0
      tmp1 <- tmp1[ind]
      tmp2 <- tmp2[ind]		
      br.len2 <- br.len[ind]		
      
      # Weighted UniFrac Distance
      unifracs[i, j, "d_W"] <- unifracs[j, i, "d_W"] <- 
        sum(br.len2 * abs(tmp1 - tmp2)) / sum(br.len2 * (tmp1 + tmp2)) 
      
      #	Unweighted UniFrac Distance
      tmp1 <- (tmp1 != 0)
      tmp2 <- (tmp2 != 0)			
      unifracs[i, j, "d_UW"] <- unifracs[j, i, "d_UW"] <- 
        sum(br.len2 * abs(tmp1 - tmp2)) / sum(br.len2)
      }
      else{
        unifracs[i, j, "d_UW"] <- unifracs[j, i, "d_UW"] <- 
          unifracs[i, j, "d_W"] <- unifracs[j, i, "d_W"] <- 0
      }
    }
  }
  return(list(unifracs=unifracs))
}


GUniFrac <- function(datP, tree, alpha) {
  # Calculate Generilized UniFrac distances
  #	
  # Args:
  #		datP: immune cell composition data after merge 
  #         the propotion to the node
  #		tree: rooted phylogenetic tree of R class "phylo"
  #   alpha: attenuate the contribution from the proportions, 
  #          varies from 0 to 1, 1 weighted
  #
  # Returns:
  # 	unifracs: three dimensional array containing the 
  #             weight/unweighted UniFrac distance 
  #
  
  if (!is.rooted(tree)) stop("Rooted phylogenetic tree required!")
  
  # Convert into proportions
  datP <- as.matrix(datP)
  n <- nrow(datP)
  
  # Construct the returning array
  if (is.null(rownames(datP))) {
    rownames(datP) <- paste("comm", 1:n, sep="_")
  }
  
  unifracs <- matrix(0, nrow = n, ncol = n,
                    dimnames=list(rownames(datP), rownames(datP)))
  
  datP <- datP[, c(tree$tip.label, tree$node.label)]
  br.len <- tree$edge.length
  edge <- tree$edge
  
  # caculate propotion of immune cell descending from brunch l
  branch_P <- matrix(0, ncol =  nrow(edge),nrow = n) 
  for(l in 1:nrow(edge)){
    branch_P[, l] <- datP[, edge[l,2]]
  }
  
  
  # Calculate various UniFrac distances
  for (i in 2:n) {
    for (j in 1:(i-1)) {
      tmp1 <- branch_P[i, ] # sample i
      tmp2 <- branch_P[j, ] # sample j
      if(any(tmp1 != tmp2)){
        ind <- (tmp1 + tmp2) != 0
        tmp1 <- tmp1[ind]
        tmp2 <- tmp2[ind]		
        br.len2 <- br.len[ind]		
        
        # Weighted UniFrac Distance
        unifracs[i, j]  <- unifracs[j, i] <-
          sum(br.len2 *(tmp1 + tmp2)^alpha * abs((tmp1 - tmp2)/(tmp1 + tmp2) )) / sum(br.len2 * (tmp1 + tmp2)^alpha) 
        
      }
      else{
        unifracs[i, j] <- unifracs[j, i] <- 0
      }
    }
  }
  return(unifracs)
}

NodeProp22 <- function(datC, tree){
  nms = c('B_naive', 'B_memory', 'Plasma', 'CD8+_T', 'CD4+T_naive', 
          'CD4+T_resting', 'CD4+T_activated', 'Tfh', 'Tregs', 'Tgd',
          'NK_resting', 'NK_activated', 'Monocytes', 'Macrophages_M0',
          'Macrophages_M1','Macrophages_M2','Dendritic_resting', 'Dendritic_activated', 'Mast_resting',
          'Mast_activated','Eosinophils','Neutrophils')
  colnms = c('B cells naive' , 'B cells memory' , 'Plasma cells' , 'T cells CD8' , 'T cells CD4 naive' ,
             'T cells CD4 memory resting' , 'T cells CD4 memory activated' , 'T cells follicular helper' ,
             'T cells regulatory (Tregs)' , 'T cells gamma delta' , 'NK cells resting' , 'NK cells activated' ,
             'Monocytes' , 'Macrophages M0' , 'Macrophages M1' , 'Macrophages M2' , 'Dendritic cells resting' ,
             'Dendritic cells activated' , 'Mast cells resting' , 'Mast cells activated' , 'Eosinophils' , 'Neutrophils')
  colnames(datC) = nms[match(colnames(datC), colnms)]
  
  n = nrow(datC)
  edge = tree$edge
  edge2 = edge[,2]
  nbr = nrow(edge)
  ntip = length(tree$tip.label)
  labs = c(tree$tip.label, tree$node.label)
  
  Mono = datC[, 'Monocytes']
  datC = datC[, tree$tip.label]
  cum <- matrix(0, n, nbr + 1)							# Branch abundance matrix
  colnames(cum) <- c(labs[edge2], "Multipoential_Hematopoietic_Stem_Cell")
  rownames(cum) <- rownames(datC)
  for (i in 1:ntip) {
    tip.loc <- which(edge2 == i)
    cum[,tip.loc] <- cum[,tip.loc] + datC[, i]	
    node <- edge[tip.loc, 1]						# Assume the direction of edge 
    node.loc <- which(edge2 == node)
    while (length(node.loc)) {
      cum[,node.loc] <- cum[,node.loc] + datC[, i]
      node <- edge[node.loc, 1]
      node.loc <- which(edge2 == node)
    }
  }
  cum[, 'Monocytes'] <- cum[, 'Monocytes'] + Mono
  cum[, 'Myeloid'] <- cum[, 'Myeloid'] + Mono
  cum[, "Multipoential_Hematopoietic_Stem_Cell"] <- 1
  return(cum)
}

CellType12 <- function(datC){
  # collapse composition results to large group and rowsum to 1
  datC = data.matrix(datC)
  ctL = list()
  ctL[["T cells follicular helper"]] = "T cells follicular helper"
  ctL[["T cells regulatory (Tregs)"]] =  "T cells regulatory (Tregs)"
  ctL[["T CD4"]] =c("T cells CD4 memory activated", "T cells CD4 memory resting", "T cells CD4 naive")
  ctL[["CD8+_T"]] = c("T cells CD8")
  ctL[["B lineage"]]      = c("B cells naive", "B cells memory", "Plasma cells")
  ctL[["NK"]]     = c("NK cells resting", "NK cells activated")
  ctL[["Dendritic"]]   = c("Dendritic cells activated", "Dendritic cells resting")
  ctL[["Macrophages"]] = c("Macrophages M0", "Macrophages M1", "Macrophages M2")
  ctL[["Mast"]]        = c("Mast cells activated", "Mast cells resting")
  ctL[["Monocytes"]]   = "Monocytes"
  ctL[["Eosinophils"]] = "Eosinophils"
  ctL[["Neutrophils"]] = "Neutrophils"
  
  
  datC1 = matrix(NA, nrow = nrow(datC), ncol = length(ctL))
  colnames(datC1) <- names(ctL)
  rownames(datC1) <- rownames(datC)
  for(i in names(ctL)){
    datC1[,i] = apply(as.matrix(datC[, ctL[[i]]]), 1, sum)
  }
  return(datC1/rowSums(datC1))
}

NodeProp <- function(datC, tree){
  # Calculate Proportion for all the nodes
  #	
  # Args:
  #		datC: immune cell composition data 
  #		tree: rooted phylogenetic tree of R class "phylo"
  #
  # Returns:
  # 	datP: immune cell composition data after merge 
  #         the propotion to the node
  if (!is.rooted(tree)) stop("Rooted phylogenetic tree required!")
  
  # Convert into proportions
  
  n <- nrow(datC)
  datP <- matrix(0, nrow = n, ncol = length(c(tree$tip.label, tree$node.label)))
  colnames(datP) = c(tree$tip.label, tree$node.label)
  rownames(datP) = rownames(datC)
  
  datC = data.matrix(datC)
  if( ! all(colnames(datC) %in% c(tree$tip.label, tree$node.label))){  
    index = match(c("T cells follicular helper",  "T cells regulatory (Tregs)", "T CD4", "B lineage" ), colnames(datC) )
    colnames(datC)[index[1:3]] = c('Tfh','Tregs', 'other_CD4+T')
    colnames(datC)[index[4]] = "B"
  }
  datP[, colnames(datC)] = datC
  
  edge <- as.matrix(tree$edge[order(tree$edge[,2]),])
  
  datP[, 'Monocytes'] = apply(datC[, c("Macrophages", "Dendritic", "Monocytes")], 1, sum)
  datP[, 'Myeloid'] = apply(datP[, c("Mast",'Neutrophils','Eosinophils', 'Monocytes')], 1, sum)
  datP[, "CD4+_T"] = apply(datP[ , c('Tfh','Tregs', 'other_CD4+T')], 1, sum)
  datP[, 'T'] = apply(datP[, c("CD4+_T",'CD8+_T') ],1 , sum)
  datP[, 'Lymphoid'] = apply(datP[, c("T",'B', "NK")], 1, sum )
  datP[, 'Multipoential_Hematopoietic_Stem_Cell'] = apply(datP[, c("Myeloid",'Lymphoid')], 1, sum )
  return(datP)
  
}



