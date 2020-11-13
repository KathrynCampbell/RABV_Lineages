#'---------------------------------------------------------
#'title: Cosmopolitan Lineage Assignment - Total Automation
#'author: Kathryn Campbell
#'date: 16/09/2020
#'---------------------------------------------------------

#KB added suggestions 08/10/20

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
library(gridExtra)

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
# #KB- can replace above 2 lines with this:
tree$tip.label <- gsub("\\..*", "", tree$tip.label, perl = T)

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
# #KB- can replace above 2 lines with this:
alignment$nam <- gsub("\\..*", "", alignment$nam, perl = T)


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

sequence_data <- lineage_assignment(tree, min.support = 95, max.support = 100, alignment, metadata, sequences = 10)[[2]]
node_data <- lineage_assignment(tree, min.support = 95, max.support = 100, alignment, metadata, sequences = 10)[[1]]

#---------------------------------------------------------------------------------------
#
# Everything above this is part of the lineage assignment script
# Everything below is added extras specific to the cosmopolitan lineages - 
# Makes a plot with informative naming
#
#-------------------------------------------------------------------------------


# #############################################
# #         RENAME THE LINEAGES               #
# #############################################
# 
# # Just doing this by eye currently, must be a better way
# # ------------------------------
node_data$cluster[1]<-"A1"
node_data$cluster[2]<-"A1.1"
node_data$cluster[3]<-"A1.1.1"
node_data$cluster[4]<-"A1.1.2"
node_data$cluster[5]<-"AF1b_A1"
node_data$cluster[6]<-"AF1b_A1.1"
node_data$cluster[7]<-"AF1b_A1.1.1"
node_data$cluster[8]<-"B1"
node_data$cluster[9]<-"ME1a_A1"
node_data$cluster[10]<-"ME1a_A1.1"
node_data$cluster[11]<-"C1"
node_data$cluster[12]<-"B1.1"
node_data$cluster[13]<-"ME1a_A1.1.1"
node_data$cluster[14]<-"AF1b_B1"
node_data$cluster[15]<-"ME1a_B1"
node_data$cluster[16]<-"ME2_A1"
node_data$cluster[17]<-"ME1a_B1.1"
node_data$cluster[18]<-"ME2_A1.1"
node_data$cluster[19]<-"AF1b_B1.1"
node_data$cluster[20]<-"AF1b_A1.2"
node_data$cluster[21]<-"ME1a_B1.1.1"
node_data$cluster[22]<-"B1.1.1"
node_data$cluster[23]<-"ME2_A1.1.1"
node_data$cluster[24]<-"AF1b_B1.1.1"
node_data$cluster[25]<-"AF1b_C1"
node_data$cluster[26]<-"AF1b_A1.2.1"
node_data$cluster[27]<-"ME1a_C1"
node_data$cluster[28]<-"C1"
node_data$cluster[29]<-"ME2_B1"
node_data$cluster[30]<-"AF1b_D1"
node_data$cluster[31]<-"AF1b_E1"
node_data$cluster[32]<-"ME1a_C1.1"
node_data$cluster[33]<-"NEE_A1"
node_data$cluster[34]<-"CA1_A1"
node_data$cluster[35]<-"A1.1.3"
node_data$cluster[36]<-"AF1b_D1.1"
node_data$cluster[37]<-"AF1b_A1.3"
node_data$cluster[38]<-"AF1a_A1"
node_data$cluster[39]<-"NEE_A1.1"
node_data$cluster[40]<-"CA1_A1.1"
node_data$cluster[41]<-"D1"
node_data$cluster[42]<-"AF1b_D1.1.1"
node_data$cluster[43]<-"AF1b_E1.1"
node_data$cluster[44]<-"ME1a_C1.1.1"
node_data$cluster[45]<-"CA1_A1.1.1"
node_data$cluster[46]<-"ME2_B1.1"
node_data$cluster[47]<-"AF1b_F1"
node_data$cluster[48]<-"AF1b_C1.1"
node_data$cluster[49]<-"AF1b_E1.1.1"
node_data$cluster[50]<-"AF1b_A1.3.1"
node_data$cluster[51]<-"AF1a_A1.1"
node_data$cluster[52]<-"A1.1.4"
node_data$cluster[53]<-"ME1a_D1"
node_data$cluster[54]<-"ME1a_C1.1.1"
node_data$cluster[55]<-"ME1a_E1"
node_data$cluster[56]<-"ME1a_A1.1.1"
node_data$cluster[57]<-"CE_A1"
node_data$cluster[58]<-"EE_A1"
node_data$cluster[59]<-"NEE_A1.1.1"
node_data$cluster[60]<-"CA1_B1"
node_data$cluster[61]<-"ME2_B1.1.1"
node_data$cluster[62]<-"D1.1"
node_data$cluster[63]<-"A1.2"
node_data$cluster[64]<-"AM2a_A1"

#-------------------------------------------------------------------------------
# # Renamed according to Rambaut et al (2020) with A1.1.1.1 becoming a new letter (e.g. C1)
# 
for (i in 1:length(node_data$cluster)) {
  sequence_data$cluster[which(sequence_data$cluster == i)] <- node_data$cluster[i]
}
# # Rename the lineages in the sequence assignment table
# 
# # Add the country of origin for each sequence for comparison later
# for (i in 1:length(sequence_data$ID)) {
#   sequence_data$country[i] <- metadata$sequence.m49_country.display_name[(which(metadata$ID == sequence_data$ID[i]))]
# }

#############################################
#             REPLOT THE TREE               #
#############################################
sequence_data$previous <- NA
for (i in 1:length(sequence_data$ID)) {
  sequence_data$previous[i]<-
    metadata$alignment.displayName[which(metadata$ID == sequence_data$ID[i])]
}

attach(tree)

sequence_data$cluster <- as.factor(sequence_data$cluster)
# Make a nice figure to save! Similar to before but with text sizes, legend etc
plot_tree<-ggtree(tree) %<+% sequence_data +
  geom_tippoint(na.rm = T, aes(colour = (cluster))) +
  theme(legend.position = c(0.1, 0.83),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8)) +
  guides(colour=guide_legend(override.aes=list(alpha=1, size=2))) +
  ggtitle("WGS")

# Plot each clade bar
# ---------------------------------------------------------------------------------------------
# GROUP 1
for (i in c(1)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.02*i, offset.text = 0)
}
# GROUP 2
for (i in c(2)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.02*i, offset.text = 0)
}
for (i in c(35)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.072, offset.text = 0)
}
for (i in c(64)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.051, offset.text = 0)
}
for (i in c(63)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.04, offset.text = 0)
}
# GROUP 3
for (i in c(3, 38)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.063, offset.text = 0)
}
for (i in c(4)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.06, offset.text = 0)
}
# GROUP 4
for (i in c(5)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.07, offset.text = 0)
}
for (i in c(52)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.06, offset.text = 0)
}
for (i in c(11, 51)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.075, offset.text = 0)
}
for (i in c(8, 16)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.08, offset.text = 0)
}
for (i in c(41)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.083, offset.text = 0)
}
# GROUP 5
for (i in c(6, 20, 37)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.08, offset.text = 0)
}
for (i in c(9, 18)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.09, offset.text = 0)
}
for (i in c(62)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.096, offset.text = 0)
}
for (i in c(12)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.085, offset.text = 0)
}
# GROUP 6
for (i in c(10, 22, 23, 58)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.1, offset.text = 0)
}
for (i in c(7, 26, 50)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.093, offset.text = 0)
}
for (i in c(33, 57)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.096, offset.text = 0)
}
# GROUP 7
for (i in c(14, 31)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.105, offset.text = 0)
}
for (i in c(13, 28, 29)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.115, offset.text = 0)
}
for (i in c(56)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.119, offset.text = 0)
}
for (i in c(25, 39)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.11, offset.text = 0)
}
# GROUP 8
for (i in c(19, 43, 48)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.12, offset.text = 0)
}
for (i in c(15, 34, 46, 55, 59)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.13, offset.text = 0)
}
# GROUP 9
for (i in c(24, 49)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.14, offset.text = 0)
}
for (i in c(17, 40, 61)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.145, offset.text = 0)
}
# GROUP 10
for (i in c(30)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.155, offset.text = 0)
}
for (i in c(21, 45)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.16, offset.text = 0)
}
# GROUP 11
for (i in c(36, 27, 60)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.17, offset.text = 0)
}
# GROUP 12
for (i in c(32, 42)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.19, offset.text = 0)
}
# GROUP 13
for (i in c(44, 47, 54)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.22, offset.text = 0)
}
# GROUP 14
for (i in c(53)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.24, offset.text = 0)
}
#-------------------------------------------------------------------------------

plot_tree
# Plot with everything on it!


ggsave("figures/Lineageplot_tree_10.png", 
       plot = last_plot(),
       height = 15, width = 30)
# Save it

#KB added row.names=F to avoid a column of row numbers
write.csv(sequence_data, "Outputs/sequence_data_cosmo10.csv", row.names=F)
write.csv(node_data, "Outputs/node_data_cosmo10.csv", row.names=F)

#############################################
#             COMBINE FIGURES               #
#############################################
tree_N <- read_annotated(filename = "Trees/051020_GLUE_CosmoSeqs_N_align.fasta.nex")
sequence_data_N <- read.csv(file = "Outputs/sequence_data_cosmoN10.csv")
node_data_N <- read.csv(file = "Outputs/node_data_cosmoN10.csv")

sequence_data_N$previous <- NA
for (i in 1:length(sequence_data_N$ID)) {
  sequence_data_N$previous[i]<-
    metadata$alignment.displayName[which(metadata$ID == sequence_data_N$ID[i])]
}

node_data_N$cluster[1]<-"A1"
node_data_N$cluster[2]<-"A1.1"
node_data_N$cluster[3]<-"A1.1.1"
node_data_N$cluster[4]<-"B1"
node_data_N$cluster[5]<-"AF1b_A1"
node_data_N$cluster[6]<-"ME1a_A1"
node_data_N$cluster[7]<-"ME2_A1"
node_data_N$cluster[8]<-"A1.2"
node_data_N$cluster[9]<-"B1.1"
node_data_N$cluster[10]<-"ME2_A1.1"
node_data_N$cluster[11]<-"A1.2.1"
node_data_N$cluster[12]<-"CA1_A1"
node_data_N$cluster[13]<-"ME2_A1.1.1"
node_data_N$cluster[14]<-"ME1a_A1.1"
node_data_N$cluster[15]<-"AF1b_A1.1"
node_data_N$cluster[16]<-"NEE_A1"
node_data_N$cluster[17]<-"CA1_A1.1"
node_data_N$cluster[18]<-"ME2_B1"
node_data_N$cluster[19]<-"ME1a_A1.1.1"
node_data_N$cluster[20]<-"C1"
node_data_N$cluster[21]<-"CA1_A1.1.1"
node_data_N$cluster[22]<-"ME2_B1.1"
node_data_N$cluster[23]<-"ME1a_B1"
node_data_N$cluster[24]<-"A1.1.3"
node_data_N$cluster[25]<-"CA2_A1"
node_data_N$cluster[26]<-"AF1b_A1.1.1"
node_data_N$cluster[27]<-"AF1b_A1.2"
node_data_N$cluster[28]<-"CE_A1"
node_data_N$cluster[29]<-"WE_A1"
node_data_N$cluster[30]<-"EE_A1"
node_data_N$cluster[31]<-"NEE_A1.1"
node_data_N$cluster[32]<-"CA1_B1"
node_data_N$cluster[33]<-"ME2_B1.1.1"
node_data_N$cluster[34]<-"ME1a_B1.1"
node_data_N$cluster[35]<-"ME1a_A1.2"
node_data_N$cluster[36]<-"ME1a_A1.3"
node_data_N$cluster[37]<-"CA2_A1.1"
node_data_N$cluster[38]<-"C1.1"
node_data_N$cluster[39]<-"AM2a_A1"
node_data_N$cluster[40]<-"D1"
node_data_N$cluster[41]<-"AF1b_B1"
node_data_N$cluster[42]<-"AF1b_A1.3"
node_data_N$cluster[43]<-"AF1b_A1.2.1"
node_data_N$cluster[44]<-"E1"

for (i in 1:length(node_data_N$cluster)) {
  sequence_data_N$cluster[which(sequence_data_N$cluster == i)] <- node_data_N$cluster[i]
}
# Rename the lineages in the sequence assignment table


sequence_data_N$cluster <- as.factor(sequence_data_N$cluster)
# Make a nice figure to save! Similar to before but with text sizes, legend etc
plot_tree_N<-ggtree(tree_N) %<+% sequence_data_N +
  geom_tippoint(na.rm = T, aes(colour = (cluster))) +
  theme(legend.position = c(0.05, 0.83),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8)) +
  guides(colour=guide_legend(override.aes=list(alpha=1, size=2))) +
  ggtitle("N gene")
# Plot each clade bar
# ---------------------------------------------------------------------------------------------
# GROUP 1
for (i in c(1, 44)) {
  plot_tree_N <-plot_tree_N +
    geom_cladelabel(node_data_N$Node[i], node_data_N$cluster[i], offset = 0, offset.text = 0)
}
# GROUP 2
for (i in c(2, 5)) {
  plot_tree_N <-plot_tree_N +
    geom_cladelabel(node_data_N$Node[i], node_data_N$cluster[i], offset = 0.055, offset.text = 0)
}
for (i in c(8)) {
  plot_tree_N <-plot_tree_N +
    geom_cladelabel(node_data_N$Node[i], node_data_N$cluster[i], offset = 0.012, offset.text = 0)
}
for (i in c(25)) {
  plot_tree_N <-plot_tree_N +
    geom_cladelabel(node_data_N$Node[i], node_data_N$cluster[i], offset = 0.074, offset.text = 0)
}
# GROUP 3
for (i in c(3, 6, 24)) {
  plot_tree_N <-plot_tree_N +
    geom_cladelabel(node_data_N$Node[i], node_data_N$cluster[i], offset = 0.065, offset.text = 0)
}
for (i in c(27)) {
  plot_tree_N <-plot_tree_N +
    geom_cladelabel(node_data_N$Node[i], node_data_N$cluster[i], offset = 0.078, offset.text = 0)
}
for (i in c(7, 15)) {
  plot_tree_N <-plot_tree_N +
    geom_cladelabel(node_data_N$Node[i], node_data_N$cluster[i], offset = 0.069, offset.text = 0)
}
for (i in c(42)) {
  plot_tree_N <-plot_tree_N +
    geom_cladelabel(node_data_N$Node[i], node_data_N$cluster[i], offset = 0.066, offset.text = 0)
}
for (i in c(11)) {
  plot_tree_N <-plot_tree_N +
    geom_cladelabel(node_data_N$Node[i], node_data_N$cluster[i], offset = 0.022, offset.text = 0)
}
for (i in c(37)) {
  plot_tree_N <-plot_tree_N +
    geom_cladelabel(node_data_N$Node[i], node_data_N$cluster[i], offset = 0.09, offset.text = 0)
}
# GROUP 4
for (i in c(4)) {
  plot_tree_N <-plot_tree_N +
    geom_cladelabel(node_data_N$Node[i], node_data_N$cluster[i], offset = 0.08, offset.text = 0)
}
for (i in c(10)) {
  plot_tree_N <-plot_tree_N +
    geom_cladelabel(node_data_N$Node[i], node_data_N$cluster[i], offset = 0.085, offset.text = 0)
}
for (i in c(14, 26, 35, 36)) {
  plot_tree_N <-plot_tree_N +
    geom_cladelabel(node_data_N$Node[i], node_data_N$cluster[i], offset = 0.083, offset.text = 0)
}
for (i in c(20)) {
  plot_tree_N <-plot_tree_N +
    geom_cladelabel(node_data_N$Node[i], node_data_N$cluster[i], offset = 0.04, offset.text = 0)
}
for (i in c(40)) {
  plot_tree_N <-plot_tree_N +
    geom_cladelabel(node_data_N$Node[i], node_data_N$cluster[i], offset = 0.05, offset.text = 0)
}
for (i in c(43)) {
  plot_tree_N <-plot_tree_N +
    geom_cladelabel(node_data_N$Node[i], node_data_N$cluster[i], offset = 0.093, offset.text = 0)
}
# GROUP 5
for (i in c(13, 28, 30)) {
  plot_tree_N <-plot_tree_N +
    geom_cladelabel(node_data_N$Node[i], node_data_N$cluster[i], offset = 0.1, offset.text = 0)
}
for (i in c(9, 19, 16)) {
  plot_tree_N <-plot_tree_N +
    geom_cladelabel(node_data_N$Node[i], node_data_N$cluster[i], offset = 0.097, offset.text = 0)
}
for (i in c(29)) {
  plot_tree_N <-plot_tree_N +
    geom_cladelabel(node_data_N$Node[i], node_data_N$cluster[i], offset = 0.105, offset.text = 0)
}
for (i in c(38)) {
  plot_tree_N <-plot_tree_N +
    geom_cladelabel(node_data_N$Node[i], node_data_N$cluster[i], offset = 0.08, offset.text = 0)
}
for (i in c(39)) {
  plot_tree_N <-plot_tree_N +
    geom_cladelabel(node_data_N$Node[i], node_data_N$cluster[i], offset = 0.07, offset.text = 0)
}
for (i in c(41)) {
  plot_tree_N <-plot_tree_N +
    geom_cladelabel(node_data_N$Node[i], node_data_N$cluster[i], offset = 0.11, offset.text = 0)
}
# GROUP 6
for (i in c(18, 23, 31)) {
  plot_tree_N <-plot_tree_N +
    geom_cladelabel(node_data_N$Node[i], node_data_N$cluster[i], offset = 0.12, offset.text = 0)
}
for (i in c(12)) {
  plot_tree_N <-plot_tree_N +
    geom_cladelabel(node_data_N$Node[i], node_data_N$cluster[i], offset = 0.115, offset.text = 0)
}
# GROUP 7
for (i in c(17, 22, 34)) {
  plot_tree_N <-plot_tree_N +
    geom_cladelabel(node_data_N$Node[i], node_data_N$cluster[i], offset = 0.14, offset.text = 0)
}
# GROUP 8
for (i in c(21, 33)) {
  plot_tree_N <-plot_tree_N +
    geom_cladelabel(node_data_N$Node[i], node_data_N$cluster[i], offset = 0.155, offset.text = 0)
}
# GROUP 9
for (i in c(32)) {
  plot_tree_N <-plot_tree_N +
    geom_cladelabel(node_data_N$Node[i], node_data_N$cluster[i], offset = 0.17, offset.text = 0)
}

#-------------------------------------------------------------------------------

plot_tree_N <- flip(plot_tree_N, 872, 571)
plot_tree_N
# Plot with everything on it!

ggsave("figures/Lineageplot_tree_N10.png", 
       plot = last_plot(),
       height = 15, width = 30)
# Save it

combined <- grid.arrange(plot_tree, plot_tree_N, nrow = 2)

ggsave("figures/Lineageplot_tree_combined.png", 
       plot = combined,
       height = 20, width = 30)
# Save it

for (i in 1:567) {
  sequence_data$"N_cluster"[i] <- sequence_data_N$cluster[which(sequence_data_N$ID == sequence_data$ID[i])]
}

sequence_data$"N_cluster"[1] 

sequence_data_N$cluster <- as.character(sequence_data_N$cluster)
