#'---
#'title: Assign Lineages Function
#'author: Kathryn Campbell
#'date: 04/09/2020
#'---
#'
#'
#'**REMOVE NODES**
#'=================
#'
#'File: ancestor_difference.r
#'=====================


assign_lineages <- function(nodes, tree) {
  nodes$overlaps <- NA 
  for (i in 1:length(nodes$Node)) {
    nodes$overlaps[i] <- length(which((allDescendants(tree)[[(nodes[i,1])]]) %in% nodes[,1]))
  } 
  
  # Create a data frame for lineage assignments. Add the tip labels, and a column ready to add the lineage they're assigned to
  lineage_assignments <- data.frame(tip = tree$tip.label, cluster = NA) 
  
  # Order the nodes of interest by the number of times they overlap the other nodes of interest (descending)
  nodes_diff <- nodes[order(-nodes$overlaps),]
  
  # Add a column called cluster and label the clusters
  nodes_diff$cluster <- c(1:(length(nodes_diff$Node)))
  
  
  for (i in 1:(length(nodes_diff$Node))) {
    lineage_assignments[which(lineage_assignments[,1] %in% clade.members(nodes_diff[i,1], tree, include.nodes = F, tip.labels = T)), 2] <- nodes_diff[i,5]
  }
  # For each sequence, see if it's a member of a lineage. If yes, put the number of the cluster in it's lineage assignment
  # Do this in order of the node with the most overlaps to the least, to ensure the assignment is at the lowest possible level
  # E.g. if a sequence is in clusters 1-7, it will appear as 7 
  
  summary <- data.frame(cluster = nodes_diff$cluster, count = NA)
  
  for (i in 1:(length(summary$cluster))) {
    summary$count[i] <- length(which(lineage_assignments$cluster == summary$cluster[i]))
  }
  # Count the number of sequences assigned to each lineage
  
  nodes_diff <- nodes_diff[-c(which(nodes_diff$cluster %in% summary$cluster[(which(summary$count < 2))])),]
  # If any lineages have no sequences in them, remove them as an option from the nodes_diff table
  
  min <- min(summary$count)
  
  while (min < 2){
    nodes_diff <- nodes_diff[order(-nodes_diff$overlaps),]
    nodes_diff$cluster <-c(1:(length(nodes_diff$Node)))
    lineage_assignments$cluster <- NA
    for (i in c(1:(length(nodes_diff$Node)))) {
      lineage_assignments[which(lineage_assignments[,1] %in% clade.members((nodes_diff[i,1]), tree, include.nodes = F, tip.labels = T)),2]<-nodes_diff[i,5]
    }
    summary <- data.frame(cluster = nodes_diff$cluster, count = NA)
    
    for (i in 1:(length(summary$cluster))) {
      summary$count[i] <- length(which(lineage_assignments$cluster == summary$cluster[i]))
    }
    
    min <- min(summary$count)
    
    if (min == 2) {
      print("done")
    } else {
      nodes_diff<-nodes_diff[-c(which(nodes_diff$cluster %in% summary$cluster[(which(summary$count < 2))])), ]
    }
    
    things<-list(lineage_assignments, nodes_diff)
    
    return(things)
  }
}

# Repeat the above steps until there are no clusters with 0 sequences left 
