#'---
#'title: Difference From Ancestor Function
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


ancestor_difference <- function(nodes, alignment.matrix, 
                                sequence.data, tree) {
  nodes$diff <- NA # Add a column in nodes to count the number of nucleotide differences each cluster has from the old seq
  
  # For each node of interest, find all the tips
  # Make a note of the differences between the oldest seq in the each cluster/lineage and one of the seqs in the lineage
  # Which differences between the old seq and each seq are shared between all the seqs in the lineage
  # E.g. which lineages show one or more shared nucleotides differences from the ancestor
  # Count these differences and add them to the table to be analysed further (may just be n's)
  
  for (i in 1:length(nodes$Node)) {
    cm <- clade.members(nodes[i,1], tree, include.nodes = F, tip.labels = T)
    seq_cm <- which(sequence.data$ID %in% cm)
    
    old <- which(row.names(alignment.matrix) %in% (
      sequence.data$ID[seq_cm[which(sequence.data$Year[seq_cm] == min(sequence.data$Year[seq_cm]))]] # This row is still a little confusing!
    ))
    old <- old[1]
    
    tips <- which(row.names(alignment.matrix) %in% cm)
    tips <- tips[-c(which(tips == old))]
    x <- which(alignment.matrix[old,] != alignment.matrix[(tips[1]),])
    
    for (j in tips[-c(1)]) {
      x <- x[which(x %in% (which(alignment.matrix[old,] != alignment.matrix[j,])))]
      print(x)
      nodes$diff[i] <- length(x)
    }
  }
  
  nodes_diff <- nodes[(which(nodes[,3]!=0)),] # Get rid of the ones with no differences straight away 
  
  nodes_diff
}
