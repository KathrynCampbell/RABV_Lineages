#'---------------------------------------------------------
#'title: Cosmopolitan Lineage Assignment - Total Automation N Gene
#'author: Kathryn Campbell
#'date: 05/10/2020
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
tree <- read_annotated(file="Trees/051020_GLUE_CosmoSeqs_N_align.fasta.nex") # GLUE sequences

# Need to edit node.comment so it's in the format required 
tree$node.comment<-gsub(".*&bootstrap=", "", tree$node.comment)


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

sequence_data <- lineage_assignment(tree, min.support = 95, max.support = 100, alignment, metadata, sequences = 10)[[2]]
node_data <- lineage_assignment(tree, min.support = 95, max.support = 100, alignment, metadata, sequences = 10)[[1]]

#---------------------------------------------------------------------------------------
# Compare with WGS assignment
# Different number of lineages - definitely different to WGS by 7 - not yet clear which 

sequence_WGS_10s <- read.csv(file = "Outputs/sequence_data_cosmo10.csv")

for (i in 1:567) {
  sequence_WGS_10s$"N_cluster"[i] <- sequence_data$cluster[which(sequence_data$ID == sequence_WGS_10s$ID[i])]
}

sequence_data$cluster <- as.factor(sequence_data$cluster)
# Make a nice figure to save! Similar to before but with text sizes, legend etc
plot_tree<-ggtree(tree) %<+% sequence_data +
  geom_tippoint(na.rm = T, aes(colour = (cluster))) +
  theme(legend.position = c(0.05, 0.83),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8)) +
  guides(colour=guide_legend(override.aes=list(alpha=1, size=5)))
# Plot each clade bar
# ---------------------------------------------------------------------------------------------
# GROUP 1
for (i in c(1, 44)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0, offset.text = 0)
}
# GROUP 2
for (i in c(2, 5)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.055, offset.text = 0)
}
for (i in c(8)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.012, offset.text = 0)
}
for (i in c(25)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.074, offset.text = 0)
}
# GROUP 3
for (i in c(3, 6, 24)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.065, offset.text = 0)
}
for (i in c(27)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.078, offset.text = 0)
}
for (i in c(7, 15)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.069, offset.text = 0)
}
for (i in c(42)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.066, offset.text = 0)
}
for (i in c(11)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.022, offset.text = 0)
}
for (i in c(37)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.09, offset.text = 0)
}
# GROUP 4
for (i in c(4)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.08, offset.text = 0)
}
for (i in c(10)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.085, offset.text = 0)
}
for (i in c(14, 26, 35, 36)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.083, offset.text = 0)
}
for (i in c(20)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.04, offset.text = 0)
}
for (i in c(40)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.05, offset.text = 0)
}
for (i in c(43)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.093, offset.text = 0)
}
# GROUP 5
for (i in c(13, 28, 30)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.1, offset.text = 0)
}
for (i in c(9, 19, 16)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.097, offset.text = 0)
}
for (i in c(29)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.105, offset.text = 0)
}
for (i in c(38)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.08, offset.text = 0)
}
for (i in c(39)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.07, offset.text = 0)
}
for (i in c(41)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.11, offset.text = 0)
}
# GROUP 6
for (i in c(18, 23, 31)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.12, offset.text = 0)
}
for (i in c(12)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.115, offset.text = 0)
}
# GROUP 7
for (i in c(17, 22, 34)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.14, offset.text = 0)
}
# GROUP 8
for (i in c(21, 33)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.155, offset.text = 0)
}
# GROUP 9
for (i in c(32)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.17, offset.text = 0)
}

#-------------------------------------------------------------------------------

plot_tree <- flip(plot_tree, 872, 571)
plot_tree
# Plot with everything on it!

ggsave("figures/Lineageplot_tree_N10.png", 
       plot = last_plot(),
       height = 15, width = 30)
# Save it

write.csv(sequence_data, "Outputs/sequence_data_cosmoN10.csv", row.names=F)
write.csv(node_data, "Outputs/node_data_cosmoN10.csv", row.names=F)
