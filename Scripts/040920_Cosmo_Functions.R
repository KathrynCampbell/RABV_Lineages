#'---
#'title: Cosmopolitan Lineage Assignment
#'author: Kathryn Campbell
#'date: 22/07/2020
#'---

rm(list = ls())

# Packages installed without problem
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

# Difficult packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
library(phytools)
library(treeio)
library(ggtree)

source("R/bootstrap_support.r")
source("R/sequence_data.r")
source("R/remove_nodes.r")
source("R/ancestor_difference.r")


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
ggtree(Cosmotree) + geom_nodelab()

#############################################
#            BOOTSTRAP SUPPORT              #
#############################################
# Function in bootstrap_support.r
# Identifies nodes with at least 70% bootstrap support
# Counts the number of tips descended from each of these nodes
# Carries forwards only the nodes with 5 or more tips descended 

nodes_5 <- bootstrap_support(tree = Cosmotree)

#############################################
#            95% COVERAGE WGS               #
#############################################
# Function in sequence_data.r
# Makes a dataframe with the number of -'s and n's in the alignment for each sequence
# Calculates the length of each sequence if these were removed 

seq_data <- sequence_data(alignment = Cosmoalign)

# Function in remove_nodes.r
# Identify seqs with less than 95% coverage and corresponding to tip numbers
# List the ancestor nodes for each of these tip numbers
# Take away the number of removed tips from the previous total number of tips calculated for each node
# Again, remove any nodes with less than 5 tips descended

nodes_5 <- remove_nodes(tree = Cosmotree, sequence.data = seq_data, nodes = nodes_5)

#############################################
#         DIFFERENCE FROM ANCESTOR          #
#############################################
seq_data$Year <- NA
# Add another column to the seq data ready to fill in dates

for (i in 1:length(Cosmoalign$seq)) {
  seq_data$Year[i] <- Cosmometa$sequence.latest_collection_year[which(Cosmometa$sequence.sequenceID == seq_data$ID[i])]
}
# Add the collection year of each sequence to the table
# Use latest, as exact collection not always filled in

# Function in ancestor_difference.r
# For each node of interest, find all the tips
# Make a note of the differences between the oldest seq in the each cluster/lineage and one of the seqs in the lineage
# Which differences between the old seq and each seq are shared between all the seqs in the lineage
# E.g. which lineages show one or more shared nucleotides differences from the ancestor
# Count these differences and add them to the table to be analysed further (may just be n's)
# Carry forward those with differences

nodes_diff <- ancestor_difference(nodes = nodes_5,
                                  alignment.matrix = alignment_matrix,
                                  sequence.data = seq_data,
                                  tree = Cosmotree)


#############################################
#         OVERLAPPING TIPS REMOVAL          #
#############################################
nodes_diff$overlaps <- NA
for (i in c(1:(length(nodes_diff$Node)))) {
  nodes_diff$overlaps[i]<-length(which((allDescendants(Cosmotree)[[(nodes_diff[i,1])]]) %in% nodes_diff[,1]))
}
# Add a column to nodes_diff and for each node, count how many of the other nodes of interest are descended from it

m <- matrix(nrow = length(Cosmoalign$seq), ncol = 2)
lineage_assignments<-data.frame(m)
names(lineage_assignments)<-c("tip", "cluster")
lineage_assignments$tip<-Cosmotree$tip.label
# Create a data frame for lineage assignments. Add the tip labels, and a column ready to add the lineage they're assigned to

nodes_diff <- nodes_diff[order(-nodes_diff$overlaps),]
# Order the nodes of interest by the number of times they overlap the other nodes of interest (descending)

nodes_diff$cluster<-c(1:(length(nodes_diff$Node)))
# Add a column called cluster and label the clusters


for (i in c(1:(length(nodes_diff$Node)))) {
  lineage_assignments[which(lineage_assignments[,1] %in% clade.members((nodes_diff[i,1]), Cosmotree, include.nodes = F, tip.labels = T)),2]<-nodes_diff[i,5]
}
# For each sequence, see if it's a member of a lineage. If it is, put the number of the cluster in it's lineage assignment
# Do this in order of the node with the most overlaps to the least, to ensure the assignment is at the lowest possible level
# E.g. if a sequence is in clusters 1-7, it will appear as 7 

m <- matrix(nrow = length(nodes_diff$Node), ncol = 2)
summary <- data.frame(m)
names(summary) <- c("cluster", "count")
summary$cluster <- nodes_diff$cluster
            
for (i in 1:(length(summary$cluster))) {
  summary$count[i] <- length(which(lineage_assignments$cluster == summary$cluster[i]))
}
# Count the number of sequences assigned to each lineage

nodes_diff<-nodes_diff[-c(which(nodes_diff$cluster %in% summary$cluster[(which(summary$count < 5))])),]
# If any lineages have less than 5 sequences in them, remove them as an option from the nodes_diff table

min <- min(summary$count)

while (min < 5){
  nodes_diff <- nodes_diff[order(-nodes_diff$overlaps),]
  nodes_diff$cluster <-c(1:(length(nodes_diff$Node)))
  lineage_assignments$cluster <- NA
  for (i in c(1:(length(nodes_diff$Node)))) {
    lineage_assignments[which(lineage_assignments[,1] %in% clade.members((nodes_diff[i,1]), Cosmotree, include.nodes = F, tip.labels = T)),2]<-nodes_diff[i,5]
  }
  m <- matrix(nrow = length(nodes_diff$Node), ncol = 2)
  summary <- data.frame(m)
  names(summary) <- c("cluster", "count")
  summary$cluster <- nodes_diff$cluster

  for (i in 1:(length(summary$cluster))) {
    summary$count[i] <- length(which(lineage_assignments$cluster == summary$cluster[i]))
  }
  
  min <- min(summary$count)

  if (min == 5) {
    print("done")
  } else {
    nodes_diff<-nodes_diff[-c(which(nodes_diff$cluster %in% summary$cluster[(which(summary$count < 5))])), ]
  }
}



# Repeat the above steps until there are no clusters with less than 5 sequences left 

#############################################
#              PLOT THE TREE                #
#############################################
lineage_assignments$cluster <- as.factor(lineage_assignments$cluster)

tree<-ggtree(Cosmotree) %<+% lineage_assignments +
  geom_tippoint(aes(colour = (cluster)))
# initial plot of the tree; need to see this to understand the lineage names

lineage_assignments$previous<-NA
for (i in 1:567) {
  lineage_assignments$previous[i]<-
    Cosmometa$alignment.displayName[which(Cosmometa$sequence.sequenceID == lineage_assignments$tip[i])]
}
# Look at the previous lineage assigned to each sequence
# Gone through by hand and identified which previous lineages are present in each new cluster - probably a way to do this in R
# The first time a previous lineage appears in a cluster on it's own, write it down and check that everything descended from it is in that previous lineage
# If so, assign that cluster as the previous lineage (seen below), and then assign further from here
# Many previous lineages are missing from this (not in a cluster on their own)

nodes_diff$cluster[15]<-"AF1a_A1"
nodes_diff$cluster[5]<-"AF1b_A1"
nodes_diff$cluster[73]<-"AF4_A1"
nodes_diff$cluster[75]<-"AM2a_A1"
nodes_diff$cluster[72]<-"AM3a_A1"
nodes_diff$cluster[16]<-"CA1_A1"
nodes_diff$cluster[71]<-"CA2_A1"
nodes_diff$cluster[58]<-"CE_A1"
nodes_diff$cluster[62]<-"EE_A1"
nodes_diff$cluster[6]<-"ME1a_A1"
nodes_diff$cluster[13]<-"ME2_A1"
nodes_diff$cluster[63]<-"NEE_A1"
nodes_diff$cluster[60]<-"WE_A1"

for (i in 1:57) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.006, offset.text = 0)
}
# The clade bars all need to be offset by a different amount so can't see a way to automate this

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

lineage_assignments$new_cluster<-NA
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