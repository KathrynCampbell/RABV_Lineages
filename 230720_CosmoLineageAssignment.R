#'---
#'title: Cosmopolitan Lineage Assignment
#'author: Kathryn Campbell
#'date: 22/07/2020
#'---

rm(list = ls())

# load packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("treeio")
# BiocManager::install("ggtree")

library(seqinr)
library(ape)
library(dplyr)
library(TreeTools)
library(adephylo)
library(phangorn)
library(phylotate)
library(caper)
library(stringr)
library(pracma)
library(ggrepel)
library(phytools)
library(treeio)
library(ggtree)

#############################################
#            IMPORT THE DATA                #
#############################################

# Import the tree
Cosmotree <- read_annotated(file="Trees/230720_Cosmo_copy.nex.txt") # GLUE sequences

# Sequence names got messed up in MAFFT, need to fix these
Cosmotree$tip.label <- sub("(?<=\\.).*$", "", Cosmotree$tip.label, perl = T)
Cosmotree$tip.label <- gsub("\\.", "", Cosmotree$tip.label, perl = T)

# Import the metadata
Cosmometa <- read.csv("Sequences/220720_GLUE_CosmoMeta.csv")

# Import the alignment
Cosmoalign <- read.alignment("Sequences/220720_GLUE_CosmoSeqs_align.fasta", format = "fasta")

# Sequence names got messed up in MAFFT, need to fix these
Cosmoalign$nam <- sub("(?<=\\.).*$", "", Cosmotree$tip.label, perl = T) 
Cosmoalign$nam <- gsub("\\.", "", Cosmotree$tip.label, perl = T)

alignment_matrix <- as.matrix.alignment(Cosmoalign)
# Need it as a matrix for later analyses

#############################################
#            PLOT THE TREE                  #
#############################################
# NEED TO COME BACK AND MAKE THIS NICER
ggtree(Cosmotree) + geom_nodelab()

#############################################
#            BOOTSTRAP SUPPORT              #
#############################################
# Identify nodes with a bootstrap of over 70 (why would the first ~570 nodes be NA?)
nodes_70 <- which(Cosmotree$node.comment > 70 | Cosmotree$node.comment == 100); nodes_70

node_data <- data.frame(Node = nodes_70, n_tips = NA)
# Make a dataframe ready for values to be put in
# Fill the first column with the numbers of the nodes identified in the previous steps

for(i in 1:length(nodes_70)) {
  node_data[i,2] <- length(Descendants(Cosmotree, nodes_70[i], type = "tips")[[1]])
}
# For each node identified in the previous step, count the number of tips descended from that node

nodes_5 <- node_data[(which(node_data[,2]>=5)),]
# Only carry forwards nodes which have more than 5 tips descended from it
# This has been identified as the definition for a cluster in previous studies

#############################################
#            95% COVERAGE WGS               #
#############################################
# Make a dataframe ready to fill with info about number of gaps and N bases, and length of the alignment and sequence
seq_data <- data.frame(ID = Cosmoalign$nam, N = NA, "-" = NA,
                       Length_before = nchar(Cosmoalign$seq[[1]]), Length_after = NA)


for (i in 1:length(Cosmoalign$seq)) {
  seq_data$N[i] <- str_count(Cosmoalign$seq[[i]], pattern = 'n')
  seq_data$`-`[i] <- str_count(Cosmoalign$seq[[i]], pattern = '-')
  seq_data$Length_after[i] <- (seq_data$Length_before[i] - seq_data$N[i] - seq_data$`-`[i])
}
# For each sequence, count the number of n bases and gaps
# Calculate the length after removing these 

nodes_remove <- Ancestors(Cosmotree, 
                        (which(Cosmotree$tip.label 
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

nodes_new # Works!!! Checked a few manually 

nodes_5 <- nodes_new[(which(nodes_new[,2] >= 5)),] # Redo this to remove any that now have less than 5, and write over the old nodes_5 so this is updated with the new tip numbers


#############################################
#         DIFFERENCE FROM ANCESTOR          #
#############################################
seq_data$Year <- NA # Add another column to the seq data ready to fill in dates
# Add collection year of each sequence to the table (Use latest, as exact collection not always filled in)
for (i in 1:length(Cosmoalign$seq)) {
  seq_data$Year[i] <- Cosmometa$sequence.latest_collection_year[which(Cosmometa$sequence.sequenceID == seq_data$ID[i])]
} 

nodes_5$diff <- NA # Add a column in nodes_5 to count the number of nucleotide differences each cluster has from the old seq

# For each node of interest, find all the tips
# Make a note of the differences between the oldest seq in the each cluster/lineage and one of the seqs in the lineage
# Which differences between the old seq and each seq are shared between all the seqs in the lineage
# E.g. which lineages show one or more shared nucleotides differences from the ancestor
# Count these differences and add them to the table to be analysed further (may just be n's)

for (i in 1:length(nodes_5$Node)) {
  cm <- clade.members(nodes_5[i,1], Cosmotree, include.nodes = F, tip.labels = T)
  seq_cm <- which(seq_data$ID %in% cm)
  
  old <- which(row.names(alignment_matrix) %in% (
    seq_data$ID[seq_cm[which(seq_data$Year[seq_cm] == min(seq_data$Year[seq_cm]))]] # This row is still a little confusing!
  ))
  old <- old[1]
  
  tips <- which(row.names(alignment_matrix) %in% cm)
  tips <- tips[-c(which(tips == old))]
  x <- which(alignment_matrix[old,] != alignment_matrix[(tips[1]),])
  
  for (j in tips[-c(1)]) {
    x <- x[which(x %in% (which(alignment_matrix[old,] != alignment_matrix[j,])))]
    print(x)
    nodes_5$diff[i] <- length(x)
  }
}

nodes_diff <- nodes_5[(which(nodes_5[,3]!=0)),] # Get rid of the ones with no differences straight away 

# Test the ones with only a few differences to check these aren't just n's
  old<-which(row.names(alignment_matrix) %in% (
    seq_data$ID[
      which(seq_data$ID %in% clade.members(573, Cosmotree, include.nodes = F, tip.labels = T))[
        which((seq_data$Year[
          which(seq_data$ID %in% clade.members(573, Cosmotree, include.nodes = F, tip.labels = T))]) == min(
            seq_data$Year[which(seq_data$ID %in% clade.members(573, Cosmotree, include.nodes = F, tip.labels = T))]))
      ]
    ]
  ))
  old<-old[1]
  
  tips<-which(row.names(alignment_matrix) %in% clade.members(573, Cosmotree, include.nodes = F, tip.labels = T))
  tips<-tips[-c(which(tips == old))]
  x<-which(alignment_matrix[old,] != alignment_matrix[(tips[1]),])
  
  for (j in tips[-c(1)]) {
    x<-x[which(x %in% (which(alignment_matrix[old,] != alignment_matrix[j,])))]
    print(x)
    nodes_5$diff[i]<-length(x)
  }

tips
old
alignment_matrix[171, 76]
alignment_matrix[175, 76]
# Checked 2; looked fine. Worked fine in previous scripts.

#############################################
#         OVERLAPPING TIPS REMOVAL          #
#############################################
# Add a column to nodes_diff and for each node, count how many of the other nodes of interest are descended from it
nodes_diff$overlaps <- NA 
for (i in 1:length(nodes_diff$Node)) {
  nodes_diff$overlaps[i] <- length(which((allDescendants(Cosmotree)[[(nodes_diff[i,1])]]) %in% nodes_diff[,1]))
} 

# Create a data frame for lineage assignments. Add the tip labels, and a column ready to add the lineage they're assigned to
lineage_assignments <- data.frame(tip = Cosmotree$tip.label, cluster = NA) 

# Order the nodes of interest by the number of times they overlap the other nodes of interest (descending)
nodes_diff <- nodes_diff[order(-nodes_diff$overlaps),]

# Add a column called cluster and label the clusters
nodes_diff$cluster <- c(1:(length(nodes_diff$Node)))


for (i in 1:(length(nodes_diff$Node))) {
  lineage_assignments[which(lineage_assignments[,1] %in% clade.members(nodes_diff[i,1], Cosmotree, include.nodes = F, tip.labels = T)), 2] <- nodes_diff[i,5]
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
    lineage_assignments[which(lineage_assignments[,1] %in% clade.members((nodes_diff[i,1]), Cosmotree, include.nodes = F, tip.labels = T)),2]<-nodes_diff[i,5]
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

#############################################
#              PLOT THE TREE                #
#############################################
lineage_assignments$cluster <- as.factor(lineage_assignments$cluster)

lineage_assignments$previous <- NA
for (i in 1:567) {
  lineage_assignments$previous[i]<-
    Cosmometa$alignment.displayName[which(Cosmometa$sequence.sequenceID == lineage_assignments$tip[i])]
}

tree<-ggtree(Cosmotree)
# initial plot of the tree; need to see this to understand the lineage names

nodes_diff$cluster[54]<-"AF1a_A1"
nodes_diff$cluster[84]<-"AF1a_B1"
nodes_diff$cluster[85]<-"AF1a_C1"
nodes_diff$cluster[4]<-"AF1b_A1"
nodes_diff$cluster[65]<-"AF4_A1"
nodes_diff$cluster[115]<-"AM2a_A1"
nodes_diff$cluster[112]<-"AM3a_A1"
nodes_diff$cluster[30]<-"CA1_A1"
nodes_diff$cluster[105]<-"CA1_B1"
nodes_diff$cluster[106]<-"CA1_C1"
nodes_diff$cluster[22]<-"CA2_A1"
nodes_diff$cluster[58]<-"CE_A1"
nodes_diff$cluster[99]<-"CE_B1"
nodes_diff$cluster[44]<-"EE_A1"
nodes_diff$cluster[7]<-"ME1a_A1"
nodes_diff$cluster[20]<-"ME2_A1"
nodes_diff$cluster[62]<-"ME2_A1"
nodes_diff$cluster[25]<-"NEE_A1"
nodes_diff$cluster[100]<-"WE_A1"
nodes_diff$cluster[101]<-"WE_A1"
# Look at the previous lineage assigned to each sequence
# Gone through by hand and identified which previous lineages are present in each new cluster - probably a way to do this in R
# The first time a previous lineage appears in a cluster on it's own, write it down and check that everything descended from it is in that previous lineage
# If so, assign that cluster as the previous lineage (seen below), and then assign further from here
# Many previous lineages are missing from this (not in a cluster on their own)

for (i in c(1:8)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.01*i, offset.text = 0)
}
for (i in c(67:115)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.15, offset.text = 0)
}
for (i in c(48:66)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.14, offset.text = 0)
}
for (i in c(38:47)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.12, offset.text = 0)
}
for (i in c(18, 19, 24)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = -0.1+0.01*i, offset.text = 0)
}
for (i in c(9:17)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = -0.04+0.01*i, offset.text = 0)
}
for (i in c(20:23)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.06, offset.text = 0)
}
for (i in c(25, 26)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.07, offset.text = 0)
}
for (i in c(27:30, 32, 33, 36, 37)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.1, offset.text = 0)
}
for (i in c(31, 34, 35)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.11, offset.text = 0)
}
tree
# Plot with everything on it!



#############################################
#         RENAME THE LINEAGES               #
#############################################

# Just doing this by eye currently, must be a better way
nodes_diff$cluster[1]<-c("A1")
nodes_diff$cluster[2]<-"B1"
nodes_diff$cluster[3]<-"A1.1"
nodes_diff$cluster[4]<-"A1.1.1"
nodes_diff$cluster[5]<-"AF1b_A1"
nodes_diff$cluster[6]<-"ME1a_A1"
nodes_diff$cluster[7]<-"AF1b_B1"
nodes_diff$cluster[8]<-"C1"
nodes_diff$cluster[9]<-"AF1b_B1.1"
nodes_diff$cluster[10]<-"C1.1"
nodes_diff$cluster[11]<-"B1.1"
nodes_diff$cluster[12]<-"D1"
nodes_diff$cluster[13]<-"ME2_A1"
nodes_diff$cluster[14]<-"AF1b_B1.1.1"
nodes_diff$cluster[15]<-"AF1a_A1"
nodes_diff$cluster[16]<-"CA1_A1"
nodes_diff$cluster[17]<-"AF1b_A1.1"
nodes_diff$cluster[18]<-"AF1b_A1.2"
nodes_diff$cluster[19]<-"AF1b_A1.2.1"
nodes_diff$cluster[20]<-"AF1b_C1"
nodes_diff$cluster[21]<-"AF1b_A1.1.1"
nodes_diff$cluster[22]<-"B1.2"
nodes_diff$cluster[23]<-"AF1b_D1"
nodes_diff$cluster[24]<-"AF1b_E1"
nodes_diff$cluster[25]<-"AF1b_A1.3"
nodes_diff$cluster[26]<-"E1"
nodes_diff$cluster[27]<-"AF1b_D1.1"
nodes_diff$cluster[28]<-"AF1b_A1.1.2"
nodes_diff$cluster[29]<-"AF1b_E1.1"
nodes_diff$cluster[30]<-"AF1b_F1"
nodes_diff$cluster[31]<-"AF1b_A1.2.2"
nodes_diff$cluster[32]<-"AF1b_A1.3.1"
nodes_diff$cluster[33]<-"AF1b_A1.4"
nodes_diff$cluster[34]<-"AF1b_A1.5"
nodes_diff$cluster[35]<-"B1.3"
nodes_diff$cluster[36]<-"AF1b_C1.1"
nodes_diff$cluster[37]<-"AF1b_G1"
nodes_diff$cluster[38]<-"AF1b_B1.1.2"
nodes_diff$cluster[39]<-"AF1b_B1.1.3"
nodes_diff$cluster[40]<-"AF1b_B1.2"
nodes_diff$cluster[41]<-"AF1b_H1"
nodes_diff$cluster[42]<-"AF1b_I1"
nodes_diff$cluster[43]<-"AF1b_J1"
nodes_diff$cluster[44]<-"AF1a_A1.1"
nodes_diff$cluster[45]<-"AF1a_A1.2"
nodes_diff$cluster[46]<-"AF1a_A1.3"
nodes_diff$cluster[47]<-"B1.1.1"
nodes_diff$cluster[48]<-"ME1a_A1.1"
nodes_diff$cluster[49]<-"ME1a_A1.2"
nodes_diff$cluster[50]<-"ME1a_A1.3"
nodes_diff$cluster[51]<-"ME1a_A1.4"
nodes_diff$cluster[52]<-"ME1a_A1.5"
nodes_diff$cluster[53]<-"ME1a_A1.6"
nodes_diff$cluster[54]<-"ME1a_A1.7"
nodes_diff$cluster[55]<-"ME1a_A1.8"
nodes_diff$cluster[56]<-"ME1a_A1.9"
nodes_diff$cluster[57]<-"ME1a_A1.10"
nodes_diff$cluster[58]<-"CE_A1"
nodes_diff$cluster[59]<-"C1.1.1"
nodes_diff$cluster[60]<-"WE_A1"
nodes_diff$cluster[61]<-"C1.1.2"
nodes_diff$cluster[62]<-"EE_A1"
nodes_diff$cluster[63]<-"NEE_A1"
nodes_diff$cluster[64]<-"CA1_A1.1"
nodes_diff$cluster[65]<-"CA1_A1.2"
nodes_diff$cluster[66]<-"D1.1"
nodes_diff$cluster[67]<-"ME2_A1.1"
nodes_diff$cluster[68]<-"ME2_A1.2"
nodes_diff$cluster[69]<-"A1.1.2"
nodes_diff$cluster[70]<-"A1.1.3"
nodes_diff$cluster[71]<-"CA2_A1"
nodes_diff$cluster[72]<-"AM3a_A1"
nodes_diff$cluster[73]<-"AF4_A1"
nodes_diff$cluster[74]<-"E1.1"
nodes_diff$cluster[75]<-"AM2a_A1"

# Renamed according to Rambaut et al (2020) with A1.1.1.1 becoming a new letter (e.g. C1)

lineage_assignments$new_cluster <- NA
# For some reason it won't just overwrite the 'cluster' column like in the other scripts, so do this instead
for (i in 1:75) {
  lineage_assignments$new_cluster[which(lineage_assignments$cluster == i)] <- nodes_diff$cluster[i]
}
# Rename the lineages in the sequence assignment table

names(lineage_assignments) <- c('tip', 'number', 'previous', 'cluster')
names(lineage_assignments)

#############################################
#             REPLOT THE TREE               #
#############################################
attach(Cosmotree)

# Make a nice figure to save! Similar to before but with text sizes, legend etc
tree<-ggtree(Cosmotree) %<+% lineage_assignments +
  geom_tippoint(na.rm = T, aes(colour = (cluster))) +
  theme(legend.position = c(0.17, 0.83),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) +
  guides(colour=guide_legend(override.aes=list(alpha=1, size=5)))

for (i in 1:5) {
  tree<-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.006*i, offset.text = 0)
}
for (i in 6) {
  tree<-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.018+0.003*i, offset.text = 0)
}
for (i in 7:10) {
  tree<-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.008+0.003*i, offset.text = 0)
}
for (i in 11) {
  tree<-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.014, offset.text = 0)
}
for (i in 14) {
  tree<-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.008+0.003*i, offset.text = 0)
}
for (i in c(12, 13, 15, 16, 17, 18, 21, 23,24)) {
  tree<-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.0025*i, offset.text = 0)
}
for (i in c(19,20)) {
  tree<-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.0027*i, offset.text = 0)
}
for (i in c(22, 25)) {
  tree<-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.0017*i, offset.text = 0)
}
for (i in 26) {
  tree<-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.015, offset.text = 0)
}
for (i in 27:46) {
  tree<-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.07, offset.text = 0)
}
for (i in 47:75) {
  tree<-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.05, offset.text = 0)
}
# The clade bars all need to be offset by a different amount so can't see a way to automate this

tree
# Plot with everything on it!


ggsave("figures/LineageTree.png", 
       plot = last_plot(),
       height = 15, width = 25)
# Save it