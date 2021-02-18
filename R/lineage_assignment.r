#'---
#'title: Lineage Assignment Function
#'author: Kathryn Campbell
#'date: 16/09/2020
#'---
#'
#'
#'**TOTAL LINEAGE ASSIGNMENT**
#'=================
#'
#'File: lineage_assignment.r
#'=====================

lineage_assignment <- function(tree, min.support, max.support, alignment, metadata, ancestral) {
  alignment_matrix <- as.matrix.alignment(alignment)
  ancestral_matrix <- as.matrix.alignment(ancestral)
  sequences <- 10
  # Need it as a matrix for later analyses
  
  #############################################
  #            BOOTSTRAP SUPPORT              #
  #############################################
  # Identify nodes with a bootstrap of over 70 (why would the first ~570 nodes be NA?)
  nodes_70 <- which(tree$node.comment > min.support | tree$node.comment == max.support)
  nodes_70 <- nodes_70 + length(tree$tip.label)
  
  node_data <- data.frame(Node = nodes_70, n_tips = NA)
  # Make a dataframe ready for values to be put in
  # Fill the first column with the numbers of the nodes identified in the previous steps
  
  for(i in 1:length(nodes_70)) {
    node_data[i,2] <- length(Descendants(tree, nodes_70[i], type = "tips")[[1]])
  }
  # For each node identified in the previous step, count the number of tips descended from that node
  
  nodes_5 <- node_data[(which(node_data[,2]>= sequences)),]
  # Only carry forwards nodes which have more than 5 tips descended from it
  # This has been identified as the definition for a cluster in previous studies
  
  #############################################
  #            95% COVERAGE WGS               #
  #############################################
  # Make a dataframe ready to fill with info about number of gaps and N bases, and length of the alignment and sequence
  seq_data <- data.frame(ID = alignment$nam, N = NA, "gap" = NA,
                         Length_before = nchar(alignment$seq[[1]]), Length_after = NA)
  
  
  for (i in 1:length(alignment$seq)) {
    seq_data$N[i] <- str_count(alignment$seq[[i]], pattern = 'n')
    seq_data$gap[i] <- str_count(alignment$seq[[i]], pattern = '-')
    seq_data$Length_after[i] <- (seq_data$Length_before[i] - seq_data$N[i] - seq_data$gap[i])
  }
  # For each sequence, count the number of n bases and gaps
  # Calculate the length after removing these 
  
  nodes_remove <- Ancestors(tree, 
                            (which(tree$tip.label 
                                   %in% (seq_data$ID[which(seq_data$Length_after < (seq_data$Length_before * 0.95))])
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
  
  new_remove <- remove_counts[which(remove_counts[,1] %in% nodes_5[,1]),]; new_remove
  # Not all nodes are included in the nodes_5 data (some are already excluded) 
  # Get rid of the nodes not in the nodes_5 data 
  
  nodes_new<-nodes_5
  
  for (i in new_remove$Node) {
    nodes_new[which(nodes_new == i), 2] <-(nodes_5[which(nodes_5 == i), 2] - (new_remove[which(new_remove == i), 2]))
  }
  # Take away the number of removed tips from the previous total number of tips calculated for each node
  
  nodes_5 <- nodes_new[(which(nodes_new[,2] >= sequences)),] # Redo this to remove any that now have less than 5, and write over the old nodes_5 so this is updated with the new tip numbers
  
  
  #############################################
  #         DIFFERENCE FROM ANCESTOR          #
  #############################################
  seq_data$Year <- NA # Add another column to the seq data ready to fill in dates
  # Add collection year of each sequence to the table (Use latest, as exact collection not always filled in)
  for (i in 1:length(alignment$seq)) {
    seq_data$Year[i] <- metadata$year[which(metadata$ID == seq_data$ID[i])]
  } 
  
  nodes_5$diff <- NA # Add a column in nodes_5 to count the number of nucleotide differences each cluster has from the old seq
  
  # For each node of interest, find all the tips
  # Make a note of the differences between the oldest seq in the each cluster/lineage and one of the seqs in the lineage
  # Which differences between the old seq and each seq are shared between all the seqs in the lineage
  # E.g. which lineages show one or more shared nucleotides differences from the ancestor
  # Count these differences and add them to the table to be analysed further (may just be n's)
  
  nodes_reduced <- data.frame(Nodes = (nodes_5$Node - (1+length(tree$tip.label))))
  
  for (i in 1:length(nodes_5$Node)) {
    cm <- clade.members(nodes_5$Node[i], tree, include.nodes = F, tip.labels = T)
    seq_cm <- which(seq_data$ID %in% cm)
    old <- which(row.names(ancestral_matrix) == paste("NODE_", (sprintf("%07d", nodes_reduced$Nodes[i])), sep=""))
    
    tips <- which(row.names(ancestral_matrix) %in% cm)
    x <- which(ancestral_matrix[old,] != ancestral_matrix[(tips[1]),])
    
    for (j in tips[-c(1)]) {
      x <- x[which(x %in% (which(ancestral_matrix[old,] != ancestral_matrix[j,])))]
      print(x)
      nodes_5$diff[i] <- length(x)
    }
  }
  
  nodes_diff <- nodes_5[(which(nodes_5[,3]!=0)),] # Get rid of the ones with no differences straight away 
  

  #############################################
  #         OVERLAPPING TIPS REMOVAL          #
  #############################################
  # Add a column to nodes_diff and for each node, count how many of the other nodes of interest are descended from it
  nodes_diff$overlaps <- NA 
  for (i in 1:length(nodes_diff$Node)) {
    nodes_diff$overlaps[i] <- length(which((allDescendants(tree)[[(nodes_diff[i,1])]]) %in% nodes_diff[,1]))
  } 
  
  # Create a data frame for lineage assignments. Add the tip labels, and a column ready to add the lineage they're assigned to
  lineage_assignments <- data.frame(tip = tree$tip.label, cluster = NA) 
  
  # Order the nodes of interest by the number of times they overlap the other nodes of interest (descending)
  nodes_diff <- nodes_diff[order(-nodes_diff$overlaps),]
  
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
  }
  # Repeat the above steps until there are no clusters with 0 sequences left 
  
  for(i in 1:length(seq_data$ID)){
    seq_data$cluster[i]<-lineage_assignments$cluster[which(lineage_assignments$tip == seq_data$ID[i])]
  }
  
  
  data <- list(nodes_diff, seq_data)
  return(data)
}