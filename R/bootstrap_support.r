#'---
#'title: Bootstrap Support Function
#'author: Kathryn Campbell
#'date: 04/09/2020
#'---
#'
#'
#'**BOOTSTRAP SUPPORT**
#'=================
#'
#'File: bootstrap_support.r
#'=====================

bootstrap_support <- function(tree) {
  # Identify nodes with a bootstrap of over 70
  nodes_70 <- which(tree$node.comment > 70 | tree$node.comment == 100)
  
  node_data <- data.frame(Node = nodes_70, n_tips = NA)
  # Make a dataframe ready for values to be put in
  # Fill the first column with the numbers of the nodes identified in the previous steps
  
  for(i in 1:length(nodes_70)) {
    node_data[i,2] <- length(Descendants(tree, nodes_70[i], type = "tips")[[1]])
  }
  # For each node identified in the previous step, count the number of tips descended from that node
  
  nodes_5 <- node_data[ (which(node_data[,2]>=5)), ]
  # Only carry forwards nodes which have more than 5 tips descended from it
  # This has been identified as the definition for a cluster in previous studies
}