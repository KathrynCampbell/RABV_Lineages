#'---------------------------------------------------------
#'title: Cosmopolitan Lineage Assignment - Total Automation
#'author: Kathryn Campbell
#'date: 16/09/2020
#'---------------------------------------------------------

rm(list=ls())

#############################################
#            INSTALL PACKAGES               #
#############################################
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
#          SOURCE THE FUNCTION              #
#############################################

source("R/lineage_assignment.r")

#############################################
#            IMPORT THE DATA                #
#############################################
#'
#'**TREE**
#'========================================================================================================
#' The tree must contain the element 'node.comment' which contains the bootstrap support/posterior support
#' And the element 'tip.label' which lists all the sequence ID's
#' These sequence ID's must match the sequence ID's in the metadata and alignment
#'=========================================================================================================
tree <- read_annotated(file="Trees/230720_Cosmo_copy.nex.txt") # GLUE sequences

# Sequence names got messed up in MAFFT, need to fix these so they match metadata and alignment
tree$tip.label <- sub("(?<=\\.).*$", "", tree$tip.label, perl = T)
tree$tip.label <- gsub("\\.", "", tree$tip.label, perl = T)


#'**METADATA**
#'========================================================================================================
#' The metadata must contain the element 'year' which lists the collection year for each sequence 
#' And the element 'ID' which lists all the sequence ID's
#' These sequence ID's must match the sequence ID's in the tree and alignment
#'=========================================================================================================
metadata <- read.csv("Sequences/220720_GLUE_CosmoMeta.csv")

# Need to edit column names so they match what is required
metadata<- metadata %>%
  rename(ID = sequence.sequenceID,
         year = sequence.latest_collection_year)


#'**ALIGNMENT**
#'========================================================================================================
#' The alignment must contain the element 'seq' which contains the sequences
#' And the element 'nam' which lists all the sequence ID's
#' These sequence ID's must match the sequence ID's in the metadata and tree
#'=========================================================================================================
alignment <- read.alignment("Sequences/220720_GLUE_CosmoSeqs_align.fasta", format = "fasta")

# Sequence names got messed up in MAFFT, need to fix these so they match metadata and alignment
alignment$nam <- sub("(?<=\\.).*$", "", tree$tip.label, perl = T) 
alignment$nam <- gsub("\\.", "", tree$tip.label, perl = T)


#############################################
#           RUN THE ASSIGNMENT              #
#############################################
#'========================================================================================================
#' Function in 'R/lineage_assignment.r'
#' This function will split the tree into lineages, and assign a lineage number to each sequence
#' It may take a little while!
#' 
#' The lineages are defined according to Rambaut et al. (2020) in which there must be:
#'    - 70% support at the defining node
#'    - At least 5 genomes with 95% coverage
#'    - 1 or more shared nucleotide differences from the ancestral lineage
#'    - At least 1 shared nucleotide change
#' 
#' The function will return 2 elements:
#'    - A data frame with information about each sequence, including the number of 'n's, the number of '-'s,
#'      the length of the alignment, the length after removing n's and -'s, the collection year and 
#'      the assigned lineage number (called sequence_data)
#'    - A data frame with information about each node the lineages are defined from including the node number,
#'      the number of tips descended from that node, the number of shared differences between sequences in 
#'      that lineage from the ancestor, the number of other lineages that lineage overlaps with and the 
#'      lineage number assigned to that node (called node_data)
#'      
#' The function uses the arguments tree, alignment and metadata which were specified earlier
#' It also uses the min.support and max.support arguments:
#'    - If using bootstrap support, use min.support = 70 and max.support = 100
#'    - If using posterior support, use min.support = 0.7 and max.support = 1.0
#'=========================================================================================================

sequence_data <- lineage_assignment(tree, min.support = 70, max.support = 100, alignment, metadata)[[2]]
node_data <- lineage_assignment(tree, min.support = 70, max.support = 100, alignment, metadata)[[1]]

#############################################
#              PLOT THE TREE                #
#############################################
#'========================================================================================================
#' Plot an initial tree to colour each tip according to it's lineage assignment
#' Add clade bars to show where each lineage falls on the tree by editing the commented out section
#' This will be very messy and need some manual editing (especially of 'offset') - but will give and initial idea!
#' It may take a little while
#'=========================================================================================================

# sequence_data$cluster <- as.factor(sequence_data$cluster)
# 
# attach(tree)
# 
# plot_tree<-ggtree(tree) %<+% sequence_data +
#   geom_tippoint(na.rm = T, aes(colour = (cluster))) +
#   theme(legend.position = c(0.12, 0.83),
#         legend.text = element_text(size = 8),
#         legend.title = element_text(size = 8)) +
#   guides(colour=guide_legend(override.aes=list(alpha=1, size=5)))
# 
# 
# # for (i in c(1:length(sequence_data$ID))) {
# #   plot_tree <-plot_tree +
# #     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.01*i)
# # }
# 
# 
# plot_tree
#
#---------------------------------------------------------------------------------------
#
# Everything above this is part of the lineage assignment script
# Everything below is added extras specific to the cosmopolitan lineages - not essential
# but helps to make it neat, and with informative naming
#
#---------------------------------------------------------------------------------------

#############################################
#         RENAME THE LINEAGES               #
#############################################

# Just doing this by eye currently, must be a better way
node_data$cluster[1]<-"A1"
node_data$cluster[2]<-"A1.1"
node_data$cluster[3]<-"A1.2"
node_data$cluster[4]<-"AF1b_A1"
node_data$cluster[5]<-"A1.1.1"
node_data$cluster[6]<-"A1.1.2"
node_data$cluster[7]<-"ME1a_A1"
node_data$cluster[8]<-"ME1a_A1.1"
node_data$cluster[9]<-"AF1b_A1.1"
node_data$cluster[10]<-"AF1b_A1.2"
node_data$cluster[11]<-"AF1b_A1.2.1"
node_data$cluster[12]<-"AF1b_B1"
node_data$cluster[13]<-"ME1a_A1.1.1"
node_data$cluster[14]<-"ME1a_B1"
node_data$cluster[15]<-"ME1a_B1.1"
node_data$cluster[16]<-"ME1a_B1.1.1"
node_data$cluster[17]<-"ME1a_C1"
node_data$cluster[18]<-"B1"
node_data$cluster[19]<-"B1.1"
node_data$cluster[20]<-"ME2_A1"
node_data$cluster[21]<-"AF1b_A1.1.1"
node_data$cluster[22]<-"CA2_A1"
node_data$cluster[23]<-"AF1b_A1.3"
node_data$cluster[24]<-"ME1a_C1.1"
node_data$cluster[25]<-"NEE_A1"
node_data$cluster[26]<-"CA2_A1.1"
node_data$cluster[27]<-"AF1b_A1.1.2"
node_data$cluster[28]<-"AF1b_B1.1"
node_data$cluster[29]<-"NEE_A1.1"
node_data$cluster[30]<-"CA1_A1"
node_data$cluster[31]<-"AF1b_C1"
node_data$cluster[32]<-"AF1b_B1.2"
node_data$cluster[33]<-"AF1b_B1.3"
node_data$cluster[34]<-"AF1b_B1.1.1"
node_data$cluster[35]<-"CA1_A1.1"
node_data$cluster[36]<-"ME2_A1.1"
node_data$cluster[37]<-"CA2_A1.1.1"
node_data$cluster[38]<-"AF1b_C1.1"
node_data$cluster[39]<-"AF1b_B1.4"
node_data$cluster[40]<-"AF1b_D1"
node_data$cluster[41]<-"C1"
node_data$cluster[42]<-"ME1a_A1.1.2"
node_data$cluster[43]<-"ME1a_A1.1.3"
node_data$cluster[44]<-"EE_A1"
node_data$cluster[45]<-"NEE_A1.1.1"
node_data$cluster[46]<-"D1"
node_data$cluster[47]<-"E1"
node_data$cluster[48]<-"AF1b_C1.1.1"
node_data$cluster[49]<-"AF1b_E1"
node_data$cluster[50]<-"AF1b_A1.1.3"
node_data$cluster[51]<-"AF1b_A1.3.1"
node_data$cluster[52]<-"AF1b_B1.3.1"
node_data$cluster[53]<-"AF1b_A1.2.2"
node_data$cluster[54]<-"AF1a_A1"
node_data$cluster[55]<-"C1.1"
node_data$cluster[56]<-"ME1a_C1.2"
node_data$cluster[57]<-"ME1a_A1.1.4"
node_data$cluster[58]<-"CE_A1"
node_data$cluster[59]<-"EE_A1.1"
node_data$cluster[60]<-"NEE_B1"
node_data$cluster[61]<-"CA1_A1.1.1"
node_data$cluster[62]<-"ME2_B1"
node_data$cluster[63]<-"CA2_B1"
node_data$cluster[64]<-"D1.1"
node_data$cluster[65]<-"AF4_A1"
node_data$cluster[66]<-"E1.1"
node_data$cluster[67]<-"AF1b_F1"
node_data$cluster[68]<-"AF1b_A1.1.4"
node_data$cluster[69]<-"AF1b_E1.1"
node_data$cluster[70]<-"AF1b_G1"
node_data$cluster[71]<-"AF1b_H1"
node_data$cluster[72]<-"AF1b_I1"
node_data$cluster[73]<-"AF1b_A1.1.5"
node_data$cluster[74]<-"AF1b_J1"
node_data$cluster[75]<-"AF1b_A1.4"
node_data$cluster[76]<-"AF1b_B1.2.1"
node_data$cluster[77]<-"AF1b_K1"
node_data$cluster[78]<-"AF1b_B1.4.1"
node_data$cluster[79]<-"AF1b_D1.1"
node_data$cluster[80]<-"AF1b_L1"
node_data$cluster[81]<-"AF1b_A1.5"
node_data$cluster[82]<-"AF1b_A1.6"
node_data$cluster[83]<-"AF1b_A1.7"
node_data$cluster[84]<-"AF1a_B1"
node_data$cluster[85]<-"AF1a_C1"
node_data$cluster[86]<-"AF1a_A1.1"
node_data$cluster[87]<-"C1.1.1"
node_data$cluster[88]<-"ME1a_C1.1.1"
node_data$cluster[89]<-"ME1a_C1.1.2"
node_data$cluster[90]<-"ME1a_C1.1.3"
node_data$cluster[91]<-"ME1a_C1.3"
node_data$cluster[92]<-"ME1a_C1.4"
node_data$cluster[93]<-"ME1a_C1.2.1"
node_data$cluster[94]<-"ME1a_D1"
node_data$cluster[95]<-"ME1a_A1.1.5"
node_data$cluster[96]<-"ME1a_E1"
node_data$cluster[97]<-"ME1a_F1"
node_data$cluster[98]<-"CE_A1.1"
node_data$cluster[99]<-"CE_B1"
node_data$cluster[100]<-"WE_A1"
node_data$cluster[101]<-"WE_A1"
node_data$cluster[102]<-"EE_A1.1.1"
node_data$cluster[103]<-"NEE_B1.1"
node_data$cluster[104]<-"CA1_D1"
node_data$cluster[105]<-"CA1_B1"
node_data$cluster[106]<-"CA1_C1"
node_data$cluster[107]<-"ME2_A1.1.1"
node_data$cluster[108]<-"ME2_A1.2"
node_data$cluster[109]<-"ME2_A1.3"
node_data$cluster[110]<-"ME2_B1.1"
node_data$cluster[111]<-"CA2_B1.1"
node_data$cluster[112]<-"AM3a_A1"
node_data$cluster[113]<-"AF4_A1.1"
node_data$cluster[114]<-"E1.1.1"
node_data$cluster[115]<-"AM2a_A1"



# Renamed according to Rambaut et al (2020) with A1.1.1.1 becoming a new letter (e.g. C1)

for (i in 1:length(node_data$cluster)) {
  sequence_data$cluster[which(sequence_data$cluster == i)] <- node_data$cluster[i]
}
# Rename the lineages in the sequence assignment table


#############################################
#             REPLOT THE TREE               #
#############################################
attach(tree)

# Make a nice figure to save! Similar to before but with text sizes, legend etc
plot_tree<-ggtree(tree) %<+% sequence_data +
  geom_tippoint(na.rm = T, aes(colour = (cluster))) +
  theme(legend.position = c(0.12, 0.83),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8)) +
  guides(colour=guide_legend(override.aes=list(alpha=1, size=5)))


for (i in c(1:8)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.01*i, offset.text = 0)
}
for (i in c(67:115)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.15, offset.text = 0)
}
for (i in c(48:66)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.14, offset.text = 0)
}
for (i in c(38:47)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.12, offset.text = 0)
}
for (i in c(18, 19, 24)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = -0.1+0.01*i, offset.text = 0)
}
for (i in c(9:17)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = -0.04+0.01*i, offset.text = 0)
}
for (i in c(20:23)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.06, offset.text = 0)
}
for (i in c(25, 26)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.07, offset.text = 0)
}
for (i in c(27:30, 32, 33, 36, 37)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.1, offset.text = 0)
}
for (i in c(31, 34, 35)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.11, offset.text = 0)
}
# The clade bars all need to be offset by a different amount so can't see a way to automate this

plot_tree
# Plot with everything on it!


