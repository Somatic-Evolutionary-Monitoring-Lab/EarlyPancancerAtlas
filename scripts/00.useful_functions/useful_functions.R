# Useful functions for personalis ctdna analysis


# Function to find matching tumour id (genomic) given a sample hash and, optionally,
# change the sample name in the mets table to match the tree ouput clonality table
find_patient_tumour_id <- function(cruk, sample_name_hash, full_tree_output_list) {
  pattern <- paste0('^', cruk, '_')
  matching_lists <- grep(pattern, names(full_tree_output_list), value = TRUE)
  split_name <- strsplit(sample_name_hash, "--")[[1]]  
  split_name[1] <- gsub("-", ".", split_name[1]) 
  split_name[1] <- sub("^LTX[0-9]+", cruk, split_name[1])
  sample <- split_name[1]
  hash <- split_name[2]
  
  for (list_name in matching_lists) {
    clonality_table <- full_tree_output_list[[list_name]]$clonality_out$clonality_table_corrected
    
    if (!is.null(clonality_table)) {
      #matching_col <- grep(sample, colnames(clonality_table), value = TRUE)
      matching_col <- grep(paste(sample, hash, sep = "|"), colnames(clonality_table), value = TRUE)
  
      if (length(matching_col) > 0) {
        return(list_name)
      }
    }
  }
  return(NA)
}


# change the sample ids in the met table to match those in the tree output clonality tables
harmonize_sample_names <- function(met_table, full_tree_output_list) {
  for (i in seq_len(nrow(met_table))) {
    cruk_id <- met_table$cruk_id[i]
    cruk <- met_table$cruk[i]
    sample_name_hash <- met_table$sample_name_hash[i]
    split_name <- strsplit(sample_name_hash, "--")[[1]]  
    split_name[1] <- gsub("-", ".", split_name[1]) 
    split_name[1] <- sub("^LTX[0-9]+", cruk, split_name[1])
    sample <- split_name[1]
    hash <- split_name[2]
    
    if (!is.na(cruk_id) && cruk_id %in% names(full_tree_output_list)) {
      clonality_table <- full_tree_output_list[[cruk_id]]$clonality_out$clonality_table_corrected
      
      if (!is.null(clonality_table)) {
        #matching_col <- grep(sample_hash, colnames(clonality_table), value = TRUE)
        #matching_col <- grep(sample, colnames(clonality_table), value = TRUE)
        matching_col <- grep(paste(sample, hash, sep = "|"), colnames(clonality_table), value = TRUE)
        
        if (length(matching_col) > 0) {
          met_table$new_sample_name_hash[i] <- matching_col
        }
      }
    }
  }
  return(met_table)
}


convert_zero_pad <- function(patient_id, zero_pad) {
  # Extract prefix and numeric id
  prefix <- substr(patient_id, 1, 3)
  numeric_id <- as.numeric(substr(patient_id, 4, nchar(patient_id)))
  
  # Format the numeric id to have the desired zero padding
  padded_numeric_id <- sprintf(paste0("%0", zero_pad, "d"), numeric_id)
  
  # Combine prefix and padded numeric id
  new_id <- paste0(prefix, padded_numeric_id)
  
  return(new_id)
}


# Plot a tree
plot_tree <- function( tree, clone_colours, edgelength = 0.9, node_size = 20, node.shape = 21,
                       boarder_type = NA, line.thickness = 2, dot_type = NA, tracked_muts = NA,
                       plot_names = FALSE, patient = NA){
  
  tree <- rbind(c("0", find_root(tree)), tree)
  edgelength <- c(setNames(0, "0"), edgelength)
  tree_graph <- igraph::as.undirected(igraph::graph.data.frame(tree))
  node.indx <- igraph::V(tree_graph)$name
  l.tree <- igraph::layout_as_tree(tree_graph, root = '0', flip.y = FALSE)
  col.palette <- clone_colours[node.indx[-1]]
  color <- names(col.palette)
  node.size <- rep(node_size, length(node.indx) - 1)
  
  tree.plot.nodes <- data.table(Node = node.indx[-1], 
                                x = l.tree[-1, 1], 
                                y = l.tree[-1, 2],
                                size = node.size,
                                shape = node.shape,
                                color = color)
  tree.plot.edges <- data.table(Edge = apply(tree, 1, function(x) paste0(x, collapse = "-")),
                                xstart = l.tree[match(tree[,1], node.indx), 1],
                                xend = l.tree[match(tree[,2], node.indx), 1],
                                ystart = l.tree[match(tree[,1], node.indx), 2],
                                yend = l.tree[match(tree[,2], node.indx), 2])
  
  if( !all(is.na(boarder_type)) ){
    thickness <- c(2,10)[as.numeric(!boarder_type == 'absent')+1]
    colour <- c('absent','present')[as.numeric(boarder_type == 'present_relapse')+1]
    tree.plot.nodes[, bord.col := colour[ match(Node, names(boarder_type)) ] ]
    tree.plot.nodes[, bord.thick := thickness[ match(Node, names(boarder_type)) ] ]
  } else {
    tree.plot.nodes[, bord.col := 'absent' ]
    tree.plot.nodes[, bord.thick := 1 ]
  }
  colours <- c('#000000','#CC0000')
  names(colours) <- c('absent', 'present')
  
  # Assign border color and thickness based on boarder_type
  #tree.plot.nodes[, bord.col := ifelse(boarder_type[Node] == "met_unique", 'red', 'black')]
  #tree.plot.nodes[, bord.thick := ifelse(boarder_type[Node] == "met_unique", 5, 2)]
  
  #colours <- c('black' = '#000000', 'red' = '#CC0000')
  
  x_pad <- 0.1
  y_pad <- 0.3
  
  minx <- mfloor( tree.plot.edges[, min(c(xstart, xend))], x_pad )
  maxx <- mceiling( tree.plot.edges[, max(c(xstart, xend))], x_pad )
  miny <- mfloor( tree.plot.edges[, min(c(ystart, yend))], y_pad )
  maxy <- mceiling( tree.plot.edges[, max(c(ystart, yend))], y_pad )
  
  plot <- ggplot() + 
    geom_segment(data = tree.plot.edges, size = line.thickness, 
                 aes(x = xstart, y = ystart, xend = xend, yend = yend)) +
    geom_point(data = tree.plot.nodes, aes(x = x, y = y), fill = "white", 
               size = tree.plot.nodes$size, shape = tree.plot.nodes$shape) +
    geom_point(data = tree.plot.nodes, size = tree.plot.nodes$size, shape = tree.plot.nodes$shape,
               aes(x = x, y = y, fill = color, stroke = bord.thick, colour = bord.col)) +
    scale_fill_manual(values = col.palette, guide = "none") +
    scale_colour_manual(values = colours, guide = "none") +
    scale_x_continuous( breaks = seq(minx - x_pad, maxx + x_pad, x_pad ), limits = c(minx - x_pad, maxx + x_pad) ) +
    scale_y_reverse( breaks = rev(seq(miny, maxy + y_pad, y_pad )), limits = rev(c(miny, maxy + y_pad)) ) +
    theme_void()
  
  if(!all(is.na(tracked_muts))){
    tree.plot.nodes[, num_tracked := tracked_muts[ match(Node, names(tracked_muts)) ] ]
    plot <- plot + 
      geom_text(data = tree.plot.nodes, aes(x = x, y = y, label = num_tracked ), size = 10)
  }
  
  if(plot_names){
    plot <- plot + 
      geom_text(data = tree.plot.nodes, aes(x = x, y = y, label = Node ), size = 10)
  }
  
  print( plot )
  
}