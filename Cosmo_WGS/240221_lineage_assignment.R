rm(list = ls())

library(phylotate)
tree <- read_annotated("Cosmo_WGS/Trees/test.nex", format = "nexus")
data <- read.csv("Cosmo_WGS//Outputs/Cosmo_WGS_sequence_data.csv")

sequence_table<-data.frame(ID=tree$tip.label, cluster=NA)

for (i in 1:length(sequence_table$ID)) {
  sequence_table$cluster[i]<-data$cluster[which(data$ID == sequence_table$ID[i])]
}

tree$tip.cluster<-sequence_table$cluster

calculate_mode <- function(x) {
  uniqx<-unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}

assignment_table<-data.frame(ID=tree$tip.label, assignment=NA)

for (i in 1:length(assignment_table$assignment)) {
  assignment_table$assignment[i]<-calculate_mode(
    tree$tip.cluster[c(as.numeric(unlist(
    Descendants(tree, (Ancestors(tree, i, type = "parent")), type = "tips"))))]
  )
}




