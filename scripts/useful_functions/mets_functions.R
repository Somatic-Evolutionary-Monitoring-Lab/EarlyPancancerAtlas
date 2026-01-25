# Functions for metastasis analyses - clonality, seeding etc.

##################################################################################################
#########################                    Functions                    #########################
####################################################################################################

get.treeStructure <- function(tree, trunk) {
  tree.structure.final <- c()
  tree.trunk.final     <- as.numeric(trunk)
  treeInfo.tmp         <- tree
  # if only one node (trunk), return trunk as sole component of tree
  if (nrow(treeInfo.tmp) == 1 & tree.trunk.final == treeInfo.tmp[1,1] & tree.trunk.final == treeInfo.tmp[1,2]) {
    return(list(trunk = tree.trunk.final, tips = tree.trunk.final, structure = tree.trunk.final, allTreeClones = tree.trunk.final))
  }
  # all nodes
  allNodes             <- unique(as.numeric(tree))
  # all nodes only in second col, leaf nodes
  tree.tips.final      <- allNodes[which(!allNodes %in% treeInfo.tmp[, 1])]
  # start from root, find all branches leading to leaf nodes
  # collapse each path from trunk to leaf node into a string (e.g., 1:2:10:5)
  # once branch is fully traversed, remove edges from treeInfo.tmp
  while(length(which(treeInfo.tmp[, 1] == tree.trunk.final)) != 0) {
    idx <- idx.rm <- min(which(treeInfo.tmp[, 1] == tree.trunk.final))
    node <- treeInfo.tmp[idx, 2]
    while (!node %in% tree.tips.final) {
      idx <- c(idx, min(which(treeInfo.tmp[, 1] == node)))
      if (length(which(treeInfo.tmp[, 1] == node)) > 1) {
        idx.rm <- idx[length(idx)]
      } else {
        idx.rm <- c(idx.rm, min(which(treeInfo.tmp[, 1] == node)))
      }
      node <- treeInfo.tmp[min(which(treeInfo.tmp[, 1] == node)), 2]
    }
    ### collapsing by row so it's easier to remove duplicate node entries
    tree.structure.final <- c(tree.structure.final, paste0(apply(treeInfo.tmp[idx,, drop = FALSE], 1, paste0, collapse = ":"), collapse = ":"))    
    ### removing branches already covered
    treeInfo.tmp <- treeInfo.tmp[-idx.rm, , drop = FALSE]
  }
  ### removing duplicate nodes
  tree.structure.final <- unlist(lapply(lapply(strsplit(tree.structure.final, split = ":"), unique), paste0, collapse = ":"))
  return(list(trunk = tree.trunk.final, tips = tree.tips.final, structure = tree.structure.final, allTreeClones = allNodes))
}


get.tumLevel.clonality <- function(clonalities, regions.to.use = NULL) {
  ### only use regions of interest
  ### if regions.to.use == NULL then assume using all regions in the clonality table
  if (!is.null(regions.to.use)) {
    clonalities <- clonalities[regions.to.use]
  }
  ### sanity checks
  if (class(clonalities) != "data.frame") stop("Clonalities have wrong format.")
  if (nrow(clonalities) == 0) return(NULL)
  
  tmp <- sapply(seq(1, nrow(clonalities)), function(i) {
    if (all(clonalities[i,] == "clonal")) return("clonal")
    if (any(clonalities[i,] == "subclonal")) return("subclonal")
    if (all(clonalities[i,] == "absent")) return("absent")
    if (any(clonalities[i,] == "absent")) return("subclonal") 
    return(NA)
  })
  names(tmp) <- rownames(clonalities)
  return(tmp)
} 

get.seedingPattern.region <- function(tree, clonality, met.region, primary.regions) {
  ### subset clonality table to only include tree clones
  clonality <- clonality[rownames(clonality) %in% unique(c(tree)), ]
  primary.clonality <- get.tumLevel.clonality(clonality, primary.regions)
  
  clusters.met <- setNames(as.character(clonality[, met.region]), rownames(clonality))
  # clusters.in.met <- clusters.met[which(clusters.met != "absent")]
  clusters.met.subclonal <- names(which(clusters.met == "subclonal"))
  
  # if no subclonal clusters in met, met must be monoclonal origin
  if (length(clusters.met.subclonal) == 0) {
    return("monoclonal")
  } else {
    # check primary origin of subclonal clusters
    # if all subclonal clusters in met are absent from primary, still monoclonal
    if (all(primary.clonality[clusters.met.subclonal] == "absent")) {
      return("monoclonal")
    } else {
      # if any subclonal clusters from met are present in primary, its polyclonal
      return("polyclonal")
    }
  }
}

get.overallClonality <- function(tree, clonality, met.regions, primary.regions) {
  clonality <- clonality[rownames(clonality) %in% unique(c(tree)), ]
  primary.clonality <- get.tumLevel.clonality(clonality, primary.regions)
  
  # for each region, find shared clusters with primary
  allSharedClusters <- lapply(met.regions, function(met.region) {
    clusters.met <- setNames(as.character(clonality[, met.region]), rownames(clonality))
    intersect(names(which(primary.clonality != "absent")), names(which(clusters.met != "absent")))
  })
  
  names(allSharedClusters) <- met.regions
  
  # get list of clusters that all mets share with primary
  sharedClusters <- Reduce(intersect, allSharedClusters)
  
  # get any clusters that are not universally shared with the primary by all mets
  # e.g., if met 1 shares clusters 1,2,3 and met 2 shares clusters 1,2,3,4, 
  # this would return cluster 4
  additionalClusters <- sapply(allSharedClusters, function(x) {
    length(setdiff(x, sharedClusters)) > 0
  })
  
  # If some metastases have additional mutations (not universally shared with prim),
  # the spread is polyclonal
  if (any(additionalClusters == TRUE)) {
    return("polyclonal")
  } else {
    return("monoclonal")
  }
}

get.overallPhyletic <- function(tree, trunk, clonality, met.regions, primary.regions) {
  clonality <- clonality[rownames(clonality) %in% unique(c(tree)), ]
  primary.clonality <- get.tumLevel.clonality(clonality, primary.regions)
  
  ### all clusters shared between primary and met regions
  allSharedClusters.full <- lapply(met.regions, function(met.region) {
    clusters.met <- setNames(as.character(clonality[, met.region]), rownames(clonality))
    intersect(names(which(primary.clonality != "absent")), names(which(clusters.met != "absent")))
  })
  names(allSharedClusters.full) <- met.regions
  treeStructure <- get.treeStructure(tree, trunk)
  branches <- treeStructure$structure
  
  if (length(branches) == 1) {
    if (branches == trunk) {
      return("monophyletic")
    }
  }
  
  # for each met region, check which branches in the tree contain the shared clusters
  tmp <- lapply(allSharedClusters.full, function(y) {
    tmp <- lapply(y, function(x) {
      tmp.branches <- grep(paste0("^", x, ":|:", x, ":|:", x, "$"), branches, value = TRUE)
    })
    overlapBranches <- Reduce(intersect, tmp)
    return(overlapBranches)
  })
  
  # get common branches where shared clusters exist
  overlapBranches <- Reduce(intersect, tmp)
  
  # no overlapping branches = "polyphyletic"
  if (length(overlapBranches) == 0) {
    return("polyphyletic")
  }
  
  allSharedClusters <- unique(unlist(allSharedClusters.full))
  overlapedClusters <- unique(unlist(strsplit(overlapBranches, split = ":")))
  # all shared clusters are found within the overlapping branches = "monophyletic"
  if (all(allSharedClusters %in% overlapedClusters)) {
    return("monophyletic")
  } else {
    return("check")
  }
}


get.seedingClones <- function(clonality, sharedClones, treeBranches, prim = "PrimaryClonality", met = "MetClonality") {
  seedingClones <- c()
  while(any(sapply(treeBranches, length) > 0)) {
    parentNodes <- unique(sapply(treeBranches, `[`, 1))
    parentNodes <- parentNodes[!is.na(parentNodes)]
    for (x in parentNodes) {
      if (x %in% sharedClones) {
        tmpStructure <- treeBranches[which(sapply(treeBranches, `[`, 1) == x)]
        childNodes <- unique(sapply(tmpStructure, `[`, 2))    
        child.of.interest <- c()
        for (y in childNodes) {
          if (is.na(y)) {
            child.of.interest <- c(child.of.interest, FALSE)
          } else {
            if (y %in% sharedClones) {
              child.of.interest <- c(child.of.interest, TRUE)
            } else {
              child.of.interest <- c(child.of.interest, FALSE)
              if (clonality[which(rownames(clonality) == y), prim] == "absent") {
                if (clonality[which(rownames(clonality) == y), met] %in% c("subclonal", "clonal")) {
                  seedingClones <- c(seedingClones, x)
                  child.of.interest <- TRUE
                }
              }
            }
          }
        }
        if (any(child.of.interest)) {
          next
        } else {
          seedingClones <- c(seedingClones, x)
        }
      }
    }
    treeBranches <- lapply(treeBranches, `[`, -1)
  }
  return(seedingClones)
}

f.sameBranch <- function(x, y, treeBranches) {
  if (class(treeBranches) == "list") {
    branches <- intersect(grep(paste0("^", x, "$"), treeBranches), grep(paste0("^", y, "$"), treeBranches))
    if (length(branches) == 0) return(FALSE)
    return(TRUE)
  } else if (class(treeBranches) == "character") {
    branches <- intersect(grep(paste0("^", x, ":|:", x, ":|:", x, "$"), treeBranches), grep(paste0("^", y, ":|:", y, ":|:", y, "$"), treeBranches))
    if (length(branches) == 0) return(FALSE)
    return(TRUE)
  }
}

f.findParent <- function(x, y, treeBranches) {
  if (class(treeBranches) == "list") {
    branches <- intersect(grep(paste0("^", x, "$"), treeBranches), grep(paste0("^", y, "$"), treeBranches))
    tmp <- sapply(branches, function(z) {
      if (grep(paste0("^", x, "$"), treeBranches[[z]]) < grep(paste0("^", y, "$"), treeBranches[[z]])) {
        return(x)
      } else {
        return(y)
      }
    })
  } else if (class(treeBranches) == "character") {
    branches <- intersect(grep(paste0("^", x, ":|:", x, ":|:", x, "$"), treeBranches), grep(paste0("^", y, ":|:", y, ":|:", y, "$"), treeBranches))
    tmp.branches <- strsplit(treeBranches, split = ":")
    tmp <- sapply(branches, function(z) {
      if (grep(paste0("^", x, "$"), tmp.branches[[z]]) < grep(paste0("^", y, "$"), tmp.branches[[z]])) {
        return(x)
      } else {
        return(y)
      }
    })  
  }
  if (length(unique(tmp)) == 1) {
    return(unique(tmp))
  } else {
    return(NA)
  }
}

######################################################################################
## Functions to plot trees ##
######################################################################################

gettingClusterType.mets <- function(patient, clonalities, sampleOverview) {
  regions         <- sampleOverview %>% filter(tumour_id %in% patient)
  primary.regions <- as.character(regions %>% filter(sample_type == "primary") %>% pull(sample_name_hash))
  met.regions     <- as.character(regions %>% filter(sample_type == "metastasis") %>% pull(sample_name_hash))
  
  df <- data.frame(Primary = get.tumLevel.clonality(clonalities, primary.regions), Mets = get.tumLevel.clonality(clonalities, met.regions), stringsAsFactors = FALSE)
  
  df$type <- ifelse(df$Primary %in% c("clonal", "subclonal"), 
                    ifelse(df$Mets %in% c("clonal", "subclonal"), "shared", "primary"),
                    ifelse(df$Mets %in% c("clonal", "subclonal"), "metastasis", "NA"))
  return(df)
}

plottingTxTree <- function(tree, trunk, edgelength, tumorID = NULL, seedingCluster, add.title = FALSE, add.subtitle = FALSE, subtitle.label = NULL, color.mapping = "single", col.palette = c("#1f78b4"), node.shape = 21, node.size = 13, node.label.custom = NULL, node.label.type = "none", segment.size = 1, add.edge.labels = FALSE, edge.labels = NULL, edge.highlights = NULL, includeEdgeLengths = FALSE, addGermlineNode = TRUE) {
  if (addGermlineNode) {
    tree <- rbind(c("0", trunk), tree)
    edgelength <- c(setNames(0, "0"), edgelength)
  }
  g.tree <- igraph::as.undirected(igraph::graph_from_data_frame(tree))
  node.indx <- igraph::V(g.tree)$name
  l.tree <- igraph::layout_as_tree(g.tree, root = '0', flip.y = FALSE)
  if (includeEdgeLengths) {
    newPositions <- l.tree
    for (node in node.indx) {
      if (node == 0) next
      if (node == trunk) {
        newPositions[which(node.indx == node),] <- c(0, as.numeric(edgelength[node]) / sum(edgelength) * 10)
        next
      }
      
      ### get parent of current node
      parent <- as.character(tree[which(tree[,2] == node), 1])
      
      ### calculate new value for y based on mutations
      y.new <- sum(edgelength[find.pathTree(tree, trunk, node)]) / sum(edgelength) * 10
      
      ### get current positions of node and parend
      x1 <- l.tree[which(node.indx == node), 1]
      y1 <- l.tree[which(node.indx == node), 2]
      x2 <- l.tree[which(node.indx == parent), 1]
      y2 <- l.tree[which(node.indx == parent), 2]
      
      ### if x positions are the same then keep current x position
      if (x1 == x2) {
        newPositions[which(node.indx == node),] <- c(newPositions[which(node.indx == parent), 1], y.new)
        next
      }
      
      ### otherwise, calculate slope and intersect
      m <- (y2 - y1) / (x2 - x1)
      b <- y1 - m * x1
      
      ### shift the x value so that it lies on the previous edge just at the new y value
      x.new <- (y.new - b) / m
      
      newPositions[which(node.indx == node),] <- c(x.new, y.new)
    }
    l.tree <- newPositions
  }
  
  ### setting colors of nodes
  tmp <- col.palette[1]
  if (color.mapping == "single") {
    if (length(col.palette) > 1) {
      warning("Single color specified, but palette includes more colors. Using only first color.")
    }
    col.palette <- setNames(tmp, "standard")
    color <- rep("standard", length(node.indx) - 1)
  } else if (color.mapping == "perNode") {
    if (length(col.palette) != length(unique(as.numeric(tree))) - 1) {
      warning("Length of color palette and number of nodes do not match. Using only first color.")
      col.palette <- setNames(tmp, "standard")
      color <- rep("standard", length(node.indx) - 1)
    } else {
      if (is.null(names(col.palette))) {
        warning("Color palette not named. Using only first color.")
        col.palette <- setNames(tmp, "standard")
        color <- rep("standard", length(node.indx) - 1)
      } else {
        col.palette <- col.palette[node.indx[-1]]
        if (any(is.na(col.palette))) {
          warning("Names of color palette did not match nodes. Using only first color.")
          col.palette <- setNames(tmp, "standard")
          color <- rep("standard", length(node.indx) - 1)
        } else {
          color <- names(col.palette)
        }
      }
    }
  } else {
    warning(paste0("The functionality color mapping ", color.mapping, " does not exists. Reverting to single color."))
  }
  
  if (node.label.type == "custom") {
    if (!is.null(node.label.custom)) {
      node.label.info <- setNames(rep("", length(node.indx) - 1), node.indx[-1])
      node.label.info[names(node.label.custom)] <- node.label.custom
    } else {
      warning("No custom labels submitted.")
      node.label.info <- setNames(rep("", length(node.indx) - 1), node.indx[-1])
    }
  } else {
    node.label.info <- setNames(rep("", length(node.indx) - 1), node.indx[-1])
  }
  
  ### setting size of nodes
  if (length(node.size) == 1) {
    node.size <- rep(node.size, length(node.indx) - 1)
  } else if (length(node.size) == length(node.indx) - 1) {
    if (is.null(names(node.size))) {
      warning("Size vector not named. Reverting to single size.")
      node.size <- rep(13, length(node.indx) - 1) 
    } else {
      node.size <- node.size[node.indx[-1]]
      if (any(is.na(node.size))) {
        warning("Names of sizes did not match nodes. Reverting to single size.")
        node.size <- rep(13, length(node.indx) - 1)
      }
    }
  }
  
  ### setting shape of nodes
  if (length(node.shape) == 1) {
    node.shape <- rep(node.shape, length(node.indx) - 1)
  } else if (length(node.shape) == length(node.indx) - 1) {
    if (is.null(names(node.shape))) {
      warning("Shape vector not named. Reverting to single shape.")
      node.shape <- rep(21, length(node.indx) - 1) 
    } else {
      node.shape <- node.shape[node.indx[-1]]
      if (any(is.na(node.shape))) {
        warning("Names of shapes did not match nodes. Reverting to single shape.")
        node.shape <- rep(21, length(node.indx) - 1)
      }
    }
  }
  
  
  tree.plot.nodes <- data.frame(Node = node.indx[-1], 
                                x = l.tree[-1, 1], 
                                y = l.tree[-1, 2],
                                size = node.size,
                                shape = node.shape,
                                color = color,
                                node.label = node.label.info)
  tree.plot.edges <- data.frame(Edge = apply(tree, 1, function(x) paste0(x, collapse = "-")),
                                xstart = l.tree[match(tree[,1], node.indx), 1],
                                xend = l.tree[match(tree[,2], node.indx), 1],
                                ystart = l.tree[match(tree[,1], node.indx), 2],
                                yend = l.tree[match(tree[,2], node.indx), 2])
  
  if (add.edge.labels) {
    if (is.null(edge.labels)) {
      warning("No edge labels supplied. Ignoring edge labelling.")
      add.edge.labels <- FALSE
    } else if (length(grep(c("data.frame|matrix|data.table"), class(edge.labels))) > 0) {
      colnames(edge.labels) <- c("Edge", "label")
      tree.plot.edges <- plyr::join(tree.plot.edges, edge.labels, type = "left", by = "Edge")
      tree.plot.edges$label[is.na(tree.plot.edges$label)] <- ""
      tree.plot.edges$label.x <- (tree.plot.edges$xstart + tree.plot.edges$xend) / 2
      tree.plot.edges$label.y <- (tree.plot.edges$ystart + tree.plot.edges$yend) / 2
      ### commented out when using labels on edges
      tree.plot.edges$label.y[which(tree.plot.edges$label.x < 0)] <- tree.plot.edges$label.y[which(tree.plot.edges$label.x < 0)] + 0.1
      # tree.plot.edges$label.y[which(tree.plot.edges$label.x > 0)] <- tree.plot.edges$label.y[which(tree.plot.edges$label.x > 0)] - 0.05
    }
  }
  
  if (!is.null(edge.highlights)) {
    if (length(grep(c("data.frame|matrix|data.table"), class(edge.highlights))) > 0) {
      colnames(edge.highlights) <- c("Edge", "label", "color")
      tree.plot.edges <- plyr::join(tree.plot.edges, edge.highlights, type = "left", by = "Edge")
      tree.plot.edges$color[is.na(tree.plot.edges$label)] <- "transparent"
      tree.plot.edges$label[is.na(tree.plot.edges$label)] <- ""
    }
  }
  
  p <- ggplot()
  
  if (!is.null(edge.highlights)) {
    p <- p + geom_segment(data = tree.plot.edges, aes(x = xstart, y = ystart, xend = xend, yend = yend), size = segment.size * 5, color = tree.plot.edges$color)
  }
  
  p <- p + geom_segment(data = tree.plot.edges, aes(x = xstart, y = ystart, xend = xend, yend = yend), size = segment.size) +
    geom_point(data = tree.plot.nodes, aes(x = x, y = y), fill = "white", size = tree.plot.nodes$size, shape = tree.plot.nodes$shape, stroke = segment.size) +
    geom_point(data = tree.plot.nodes, aes(x = x, y = y, fill = color), size = tree.plot.nodes$size, shape = tree.plot.nodes$shape, stroke = segment.size) +
    scale_fill_manual(values = col.palette, guide = "none") +
    # Highlight seeding cluster with a dark circle
    geom_point(data = tree.plot.nodes %>% filter(Node %in% seedingCluster), aes(x = x, y = y), size = tree.plot.nodes %>% filter(Node %in% seedingCluster) %>% pull(size), shape = 1, color = 'black', stroke = 2) +
    theme_void() +
    scale_y_reverse()
  
  if (add.title & !is.null(tumorID)) {
    p <- p + ggtitle(tumorID) + theme(plot.title = element_text(hjust = 0.5))
  }
  
  if (add.subtitle) {
    p <- p + labs(subtitle = subtitle.label) + theme(plot.subtitle = element_text(hjust = 0.5))
  }
  
  if (add.edge.labels) {
    p <- p + geom_text(data = tree.plot.edges, aes(x = label.x, y = label.y, label = stringr::str_wrap(label, 20)), hjust = 0, nudge_x = 0.05, vjust = 1, size = 3)
    # p <- p + geom_text(data = tree.plot.edges, aes(x = label.x, y = label.y, label = stringr::str_wrap(label, 20)), hjust = 0.5, vjust = 0.5, size = 20)
  }
  
  if (node.label.type == "none") {
    return(p)
  } else if (node.label.type == "nodeNumber") {
    p <- p + geom_text(data = tree.plot.nodes, aes(x = x, y = y, label = Node), size = 7)
    return(p)
  } else if (node.label.type == "custom") {
    p <- p + geom_point(data = tree.plot.nodes %>% filter(node.label != ""), aes(x = x, y = y), fill = "white", color = "white", shape = 24, size = tree.plot.nodes %>% filter(node.label != "") %>% pull(size) * 0.5)
    return(p)
  } else {
    warning(paste0("Node label type ", node.label.type, " unknown. Returning basic plot."))
    return(p)
  }
}


