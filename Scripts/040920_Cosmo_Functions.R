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

source("R/bootstrap_support.r")
source("R/sequence_data.r")
source("R/remove_nodes.r")
source("R/ancestor_difference.r")
source("R/assign_lineages.r")


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
#              ASSIGN LINEAGES              #
#############################################
# Function in assign_lineages.r
# For each sequence, see if it's a member of a lineage. If yes, put the number of the cluster in it's lineage assignment
# Do this in order of the node with the most overlaps to the least, to ensure the assignment is at the lowest possible level
# E.g. if a sequence is in clusters 1-7, it will appear as 7 
# Count the number of sequences assigned to each lineage, remove any with only 1 sequence as an option
# Reassign the lineages and repeat until no lineages with only 1 seq 
# Returns a dataframe listing each tip and it's lineage assignment, and a data frame with information about each node of interest

lineage_assignments <- assign_lineages(tree = Cosmotree, nodes = nodes_diff)[[1]]
nodes_diff <- assign_lineages(tree = Cosmotree, nodes = nodes_diff)[[2]]


#############################################
#              PLOT THE TREE                #
#############################################
lineage_assignments$cluster <- as.factor(lineage_assignments$cluster)

lineage_assignments$previous <- NA
for (i in 1:567) {
  lineage_assignments$previous[i]<-
    Cosmometa$alignment.displayName[which(Cosmometa$sequence.sequenceID == lineage_assignments$tip[i])]
}

#---------------------------------
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
nodes_diff$cluster[62]<-"ME2_B1"
nodes_diff$cluster[25]<-"NEE_A1"
nodes_diff$cluster[100]<-"WE_A1"
nodes_diff$cluster[101]<-"WE_A1"
#---------------------------------
# Look at the previous lineage assigned to each sequence
# Gone through by hand and identified which previous lineages are present in each new cluster - probably a way to do this in R
# The first time a previous lineage appears in a cluster on it's own, write it down and check that everything descended from it is in that previous lineage
# If so, assign that cluster as the previous lineage (seen below), and then assign further from here
# Many previous lineages are missing from this (not in a cluster on their own)

tree<-ggtree(Cosmotree)
# initial plot of the tree; need to see this to understand the lineage names

# Plot each clade bar
# ---------------------------------------------------------------------------------------------
# GROUP 1
for (i in c(1)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.01, offset.text = 0)
}
for (i in c(46, 47)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.012, offset.text = 0)
}
for (i in c(65)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.023, offset.text = 0)
}
for (i in c(115)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.023, offset.text = 0)
}

# GROUP 2
for (i in c(2)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.03, offset.text = 0)
}
for (i in c(84)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.028, offset.text = 0)
}
for (i in c(3)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.026, offset.text = 0)
}
for (i in c(22)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.039, offset.text = 0)
}
for (i in c(85)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.037, offset.text = 0)
}
for (i in c(54)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.042, offset.text = 0)
}
for (i in c(41, 64)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.015, offset.text = 0)
}
for (i in c(113, 66)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.028, offset.text = 0)
}

# GROUP 3
for (i in c(4)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.032, offset.text = 0)
}
for (i in c(109,26)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.046, offset.text = 0)
}
for (i in c(86,62)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.052, offset.text = 0)
}
for (i in c(6, 20)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.041, offset.text = 0)
}
for (i in c(5, 112, 114)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.036, offset.text = 0)
}
for (i in c(55)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.021, offset.text = 0)
}

#  GROUP 4
for (i in c(99, 109, 37)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.06, offset.text = 0)
}
for (i in c(9, 10, 82)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.045, offset.text = 0)
}
for (i in c(23, 81, 83, 25)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.05, offset.text = 0)
}
for (i in c(7, 58, 18, 36, 108)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.053, offset.text = 0)
}
for (i in c(87)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.032, offset.text = 0)
}
for (i in c(100, 101, 44)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.056, offset.text = 0)
}

#  GROUP 5
for (i in c(53, 59, 63)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.075, offset.text = 0)
}
for (i in c(27, 68, 50, 51, 8, 19, 107)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.07, offset.text = 0)
}
for (i in c(21, 11, 29)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.065, offset.text = 0)
}
for (i in c(73, 75, 98)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.07, offset.text = 0)
}

#  GROUP 6
for (i in c(31,70,72,73,74,75,13,42,43,95,45,105,106)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.095, offset.text = 0)
}
for (i in c(49, 71, 12)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.091, offset.text = 0)
}
for (i in c(53, 102, 30)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.1, offset.text = 0)
}
for (i in c(87,64)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.078, offset.text = 0)
}
for (i in c(86, 110)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.105, offset.text = 0)
}
for (i in c(57,111)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.1, offset.text = 0)
}

#  GROUP 7
for (i in c(38,69,32,33,14,94,60)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.12, offset.text = 0)
}
for (i in c(39, 28)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.115, offset.text = 0)
}
for (i in c(80,97,35)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.124, offset.text = 0)
}

#  GROUP 8
for (i in c(48,76,52,78,15,103)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.145, offset.text = 0)
}
for (i in c(34)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.135, offset.text = 0)
}
for (i in c(61)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.15, offset.text = 0)
}

#  GROUP 9
for (i in c(48,16,104)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.17, offset.text = 0)
}
for (i in c(77)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.175, offset.text = 0)
}
for (i in c(40)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.165, offset.text = 0)
}

#  GROUP 10
for (i in c(67,79,17)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.19, offset.text = 0)
}

#  GROUP 11
for (i in c(24,91,92)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.2, offset.text = 0)
}

#  GROUP 12
for (i in c(88,89,90,93)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.215, offset.text = 0)
}
#-------------------------------------------------------------------------------
tree
# Plot with everything on it!



#############################################
#         RENAME THE LINEAGES               #
#############################################

# Just doing this by eye currently, must be a better way
# ------------------------------
nodes_diff$cluster[1]<-"A1"
nodes_diff$cluster[2]<-"A1.1"
nodes_diff$cluster[3]<-"A1.2"
nodes_diff$cluster[4]<-"AF1b_A1"
nodes_diff$cluster[5]<-"A1.1.1"
nodes_diff$cluster[6]<-"A1.1.2"
nodes_diff$cluster[7]<-"ME1a_A1"
nodes_diff$cluster[8]<-"ME1a_A1.1"
nodes_diff$cluster[9]<-"AF1b_A1.1"
nodes_diff$cluster[10]<-"AF1b_A1.2"
nodes_diff$cluster[11]<-"AF1b_A1.2.1"
nodes_diff$cluster[12]<-"AF1b_B1"
nodes_diff$cluster[13]<-"ME1a_A1.1.1"
nodes_diff$cluster[14]<-"ME1a_B1"
nodes_diff$cluster[15]<-"ME1a_B1.1"
nodes_diff$cluster[16]<-"ME1a_B1.1.1"
nodes_diff$cluster[17]<-"ME1a_C1"
nodes_diff$cluster[18]<-"B1"
nodes_diff$cluster[19]<-"B1.1"
nodes_diff$cluster[20]<-"ME2_A1"
nodes_diff$cluster[21]<-"AF1b_A1.1.1"
nodes_diff$cluster[22]<-"CA2_A1"
nodes_diff$cluster[23]<-"AF1b_A1.3"
nodes_diff$cluster[24]<-"ME1a_C1.1"
nodes_diff$cluster[25]<-"NEE_A1"
nodes_diff$cluster[26]<-"CA2_A1.1"
nodes_diff$cluster[27]<-"AF1b_A1.1.2"
nodes_diff$cluster[28]<-"AF1b_B1.1"
nodes_diff$cluster[29]<-"NEE_A1.1"
nodes_diff$cluster[30]<-"CA1_A1"
nodes_diff$cluster[31]<-"AF1b_C1"
nodes_diff$cluster[32]<-"AF1b_B1.2"
nodes_diff$cluster[33]<-"AF1b_B1.3"
nodes_diff$cluster[34]<-"AF1b_B1.1.1"
nodes_diff$cluster[35]<-"CA1_A1.1"
nodes_diff$cluster[36]<-"ME2_A1.1"
nodes_diff$cluster[37]<-"CA2_A1.1.1"
nodes_diff$cluster[38]<-"AF1b_C1.1"
nodes_diff$cluster[39]<-"AF1b_B1.4"
nodes_diff$cluster[40]<-"AF1b_D1"
nodes_diff$cluster[41]<-"C1"
nodes_diff$cluster[42]<-"ME1a_A1.1.2"
nodes_diff$cluster[43]<-"ME1a_A1.1.3"
nodes_diff$cluster[44]<-"EE_A1"
nodes_diff$cluster[45]<-"NEE_A1.1.1"
nodes_diff$cluster[46]<-"D1"
nodes_diff$cluster[47]<-"E1"
nodes_diff$cluster[48]<-"AF1b_C1.1.1"
nodes_diff$cluster[49]<-"AF1b_E1"
nodes_diff$cluster[50]<-"AF1b_A1.1.3"
nodes_diff$cluster[51]<-"AF1b_A1.3.1"
nodes_diff$cluster[52]<-"AF1b_B1.3.1"
nodes_diff$cluster[53]<-"AF1b_A1.2.2"
nodes_diff$cluster[54]<-"AF1a_A1"
nodes_diff$cluster[55]<-"C1.1"
nodes_diff$cluster[56]<-"ME1a_C1.2"
nodes_diff$cluster[57]<-"ME1a_A1.1.4"
nodes_diff$cluster[58]<-"CE_A1"
nodes_diff$cluster[59]<-"EE_A1.1"
nodes_diff$cluster[60]<-"NEE_B1"
nodes_diff$cluster[61]<-"CA1_A1.1.1"
nodes_diff$cluster[62]<-"ME2_B1"
nodes_diff$cluster[63]<-"CA2_B1"
nodes_diff$cluster[64]<-"D1.1"
nodes_diff$cluster[65]<-"AF4_A1"
nodes_diff$cluster[66]<-"E1.1"
nodes_diff$cluster[67]<-"AF1b_F1"
nodes_diff$cluster[68]<-"AF1b_A1.1.4"
nodes_diff$cluster[69]<-"AF1b_E1.1"
nodes_diff$cluster[70]<-"AF1b_G1"
nodes_diff$cluster[71]<-"AF1b_H1"
nodes_diff$cluster[72]<-"AF1b_I1"
nodes_diff$cluster[73]<-"AF1b_A1.1.5"
nodes_diff$cluster[74]<-"AF1b_J1"
nodes_diff$cluster[75]<-"AF1b_A1.4"
nodes_diff$cluster[76]<-"AF1b_B1.2.1"
nodes_diff$cluster[77]<-"AF1b_K1"
nodes_diff$cluster[78]<-"AF1b_B1.4.1"
nodes_diff$cluster[79]<-"AF1b_D1.1"
nodes_diff$cluster[80]<-"AF1b_L1"
nodes_diff$cluster[81]<-"AF1b_A1.5"
nodes_diff$cluster[82]<-"AF1b_A1.6"
nodes_diff$cluster[83]<-"AF1b_A1.7"
nodes_diff$cluster[84]<-"AF1a_B1"
nodes_diff$cluster[85]<-"AF1a_C1"
nodes_diff$cluster[86]<-"AF1a_A1.1"
nodes_diff$cluster[87]<-"C1.1.1"
nodes_diff$cluster[88]<-"ME1a_C1.1.1"
nodes_diff$cluster[89]<-"ME1a_C1.1.2"
nodes_diff$cluster[90]<-"ME1a_C1.1.3"
nodes_diff$cluster[91]<-"ME1a_C1.3"
nodes_diff$cluster[92]<-"ME1a_C1.4"
nodes_diff$cluster[93]<-"ME1a_C1.2.1"
nodes_diff$cluster[94]<-"ME1a_D1"
nodes_diff$cluster[95]<-"ME1a_A1.1.5"
nodes_diff$cluster[96]<-"ME1a_E1"
nodes_diff$cluster[97]<-"ME1a_F1"
nodes_diff$cluster[98]<-"CE_A1.1"
nodes_diff$cluster[99]<-"CE_B1"
nodes_diff$cluster[100]<-"WE_A1"
nodes_diff$cluster[101]<-"WE_A1"
nodes_diff$cluster[102]<-"EE_A1.1.1"
nodes_diff$cluster[103]<-"NEE_B1.1"
nodes_diff$cluster[104]<-"CA1_D1"
nodes_diff$cluster[105]<-"CA1_B1"
nodes_diff$cluster[106]<-"CA1_C1"
nodes_diff$cluster[107]<-"ME2_A1.1.1"
nodes_diff$cluster[108]<-"ME2_A1.2"
nodes_diff$cluster[109]<-"ME2_A1.3"
nodes_diff$cluster[110]<-"ME2_B1.1"
nodes_diff$cluster[111]<-"CA2_B1.1"
nodes_diff$cluster[112]<-"AM3a_A1"
nodes_diff$cluster[113]<-"AF4_A1.1"
nodes_diff$cluster[114]<-"E1.1.1"
nodes_diff$cluster[115]<-"AM2a_A1"
#-------------------------------------------------------------------------------
# Renamed according to Rambaut et al (2020) with A1.1.1.1 becoming a new letter (e.g. C1)

lineage_assignments$new_cluster <- NA
# For some reason it won't just overwrite the 'cluster' column like in the other scripts, so do this instead
for (i in 1:length(nodes_diff$cluster)) {
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
  theme(legend.position = c(0.1, 0.83),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8)) +
  guides(colour=guide_legend(override.aes=list(alpha=1, size=5)))

# Plot each clade bar
# ---------------------------------------------------------------------------------------------
# GROUP 1
for (i in c(1)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.01, offset.text = 0)
}
for (i in c(46, 47)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.012, offset.text = 0)
}
for (i in c(65)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.023, offset.text = 0)
}
for (i in c(115)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.023, offset.text = 0)
}

# GROUP 2
for (i in c(2)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.03, offset.text = 0)
}
for (i in c(84)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.028, offset.text = 0)
}
for (i in c(3)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.026, offset.text = 0)
}
for (i in c(22)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.039, offset.text = 0)
}
for (i in c(85)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.037, offset.text = 0)
}
for (i in c(54)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.042, offset.text = 0)
}
for (i in c(41, 64)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.015, offset.text = 0)
}
for (i in c(113, 66)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.028, offset.text = 0)
}

# GROUP 3
for (i in c(4)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.032, offset.text = 0)
}
for (i in c(109,26)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.046, offset.text = 0)
}
for (i in c(86,62)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.052, offset.text = 0)
}
for (i in c(6, 20)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.041, offset.text = 0)
}
for (i in c(5, 112, 114)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.036, offset.text = 0)
}
for (i in c(55)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.021, offset.text = 0)
}

#  GROUP 4
for (i in c(99, 109, 37)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.06, offset.text = 0)
}
for (i in c(9, 10, 82)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.045, offset.text = 0)
}
for (i in c(23, 81, 83, 25)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.05, offset.text = 0)
}
for (i in c(7, 58, 18, 36, 108)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.053, offset.text = 0)
}
for (i in c(87)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.032, offset.text = 0)
}
for (i in c(100, 101, 44)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.056, offset.text = 0)
}

#  GROUP 5
for (i in c(53, 59, 63)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.075, offset.text = 0)
}
for (i in c(27, 68, 50, 51, 8, 19, 107)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.07, offset.text = 0)
}
for (i in c(21, 11, 29)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.065, offset.text = 0)
}
for (i in c(73, 75, 98)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.07, offset.text = 0)
}

#  GROUP 6
for (i in c(31,70,72,73,74,75,13,42,43,95,45,105,106)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.095, offset.text = 0)
}
for (i in c(49, 71, 12)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.091, offset.text = 0)
}
for (i in c(53, 102, 30)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.1, offset.text = 0)
}
for (i in c(87,64)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.078, offset.text = 0)
}
for (i in c(86, 110)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.105, offset.text = 0)
}
for (i in c(57,111)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.1, offset.text = 0)
}

#  GROUP 7
for (i in c(38,69,32,33,14,94,60)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.12, offset.text = 0)
}
for (i in c(39, 28)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.115, offset.text = 0)
}
for (i in c(80,97,35)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.124, offset.text = 0)
}

#  GROUP 8
for (i in c(48,76,52,78,15,103)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.143, offset.text = 0)
}
for (i in c(34)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.133, offset.text = 0)
}
for (i in c(61)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.147, offset.text = 0)
}

#  GROUP 9
for (i in c(48,16,104)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.165, offset.text = 0)
}
for (i in c(77)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.171, offset.text = 0)
}
for (i in c(40)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.162, offset.text = 0)
}

#  GROUP 10
for (i in c(67,79,17)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.184, offset.text = 0)
}

#  GROUP 11
for (i in c(24,91,92)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.194, offset.text = 0)
}

#  GROUP 12
for (i in c(88,89,90,93)) {
  tree <-tree +
    geom_cladelabel(nodes_diff$Node[i], nodes_diff$cluster[i], offset = 0.21, offset.text = 0)
}
#-------------------------------------------------------------------------------

tree
# Plot with everything on it!


ggsave("figures/LineageTree.png", 
       plot = last_plot(),
       height = 15, width = 30)
# Save it

write.csv(sequence_data, "Outputs/sequence_data_cosmo.csv")
write.csv(node_data, "Outputs/node_data_cosmo.csv")