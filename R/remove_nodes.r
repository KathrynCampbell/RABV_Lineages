#'---
#'title: Remove <95% Coverage Sequence Nodes Function
#'author: Kathryn Campbell
#'date: 04/09/2020
#'---
#'
#'
#'**REMOVE NODES**
#'=================
#'
#'File: remove_nodes.r
#'=====================

remove_nodes <- function(tree, sequence.data, nodes) {
  nodes_remove <- Ancestors(tree, 
                            (which(tree$tip.label 
                                   %in% (sequence.data$ID[which(sequence.data$Length_after < (sequence.data$Length_before * 0.95))])
                            )), 
                            'all')
  # Identify seqs with less than 95% coverage and corresponding to tip numbers
  # List the ancestor nodes for each of these tip numbers
  
  removes <- nodes_remove[[1]]
  for (i in 2:(length(nodes_remove))) {
    removes <- c(removes, nodes_remove[[i]])
  }
  remove_counts <- data.frame(table(removes))
  # Make a table to count the number the removed sequences descended from each node (e.g. for the deeper nodes, all 10 are descended)
  
  names(remove_counts) <-c('Node', 'freq')
  # Change the names
  
  remove_counts$Node <- as.integer(levels(remove_counts$Node))
  # Need to change this, or it creates many levels and causes errors
  
  new_remove <- remove_counts[which(remove_counts[,1] %in% nodes[,1]),]; new_remove
  # Not all nodes are included in the nodes_5 data (some are already excluded) 
  # Get rid of the nodes not in the nodes_5 data 
  
  nodes_new<-nodes
  
  for (i in new_remove$Node) {
    nodes_new[which(nodes_new == i), 2] <-(nodes[which(nodes == i), 2] - (new_remove[which(new_remove == i), 2]))
  }
  # Take away the number of removed tips from the previous total number of tips calculated for each node
  
  nodes_5 <- nodes_new[(which(nodes_new[,2] >= 5)),] # Redo this to remove any that now have less than 5, and write over the old nodes_5 so this is updated with the new tip numbers
  
  nodes_5
}
