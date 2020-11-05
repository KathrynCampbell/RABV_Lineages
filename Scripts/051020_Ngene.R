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

sequence_data <- lineage_assignment(tree, min.support = 70, max.support = 100, alignment, metadata, sequences = 10)[[2]]
node_data <- lineage_assignment(tree, min.support = 70, max.support = 100, alignment, metadata, sequences = 10)[[1]]

#---------------------------------------------------------------------------------------
# Compare with WGS assignment
# Different number of lineages - definitely different to WGS by 7 - not yet clear which 

sequence_WGS_10s <- read.csv(file = "Outputs/sequence_data_cosmo10.csv")

sequence_WGS_70 <- read.csv(file = "Outputs/sequence_data_cosmo.csv")

sequence_N_80 <- read.csv(file = "Outputs/sequence_data_cosmo_N_80.csv")
sequence_WGS_80 <- read.csv(file = "Outputs/sequence_data_cosmo_80.csv")

sequence_N_90 <- read.csv(file = "Outputs/sequence_data_cosmo_N_90.csv")
sequence_WGS_90 <- read.csv(file = "Outputs/sequence_data_cosmo_90.csv")

sequence_N_100 <- read.csv(file = "Outputs/sequence_data_cosmo_N_100.csv")
sequence_WGS_100 <- read.csv(file = "Outputs/sequence_data_cosmo_100.csv")

sequence_N_65 <- read.csv(file = "Outputs/sequence_data_cosmo_N_65.csv")
sequence_WGS_65 <- read.csv(file = "Outputs/sequence_data_cosmo_65.csv")

sequence_N_60 <- read.csv(file = "Outputs/sequence_data_cosmo_N_60.csv")
sequence_WGS_60 <- read.csv(file = "Outputs/sequence_data_cosmo_60.csv")

sequence_N_50 <- read.csv(file = "Outputs/sequence_data_cosmo_N_50.csv")
sequence_WGS_50 <- read.csv(file = "Outputs/sequence_data_cosmo_50.csv")

node_data_all <- read.csv(file = "Outputs/node_data_cosmo.csv")

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
# # GROUP 1
# for (i in c(1, 17)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.01, offset.text = 0)
# }
# # GROUP 2
# for (i in c(2, 4)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.06, offset.text = 0)
# }
# for (i in c(6)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.017, offset.text = 0)
# }
# for (i in c(22)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.079, offset.text = 0)
# }
# for (i in c(96)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.075, offset.text = 0)
# }
# for (i in c(97)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.086, offset.text = 0)
# }
# for (i in c(98)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.016, offset.text = 0)
# }
# # GROUP 3
# for (i in c(3, 5, 8)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.07, offset.text = 0)
# }
# for (i in c(12, 20)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.075, offset.text = 0)
# }
# for (i in c(7)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.027, offset.text = 0)
# }
# for (i in c(27)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.095, offset.text = 0)
# }
# for (i in c(41, 47)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.082, offset.text = 0)
# }
# for (i in c(68)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.073, offset.text = 0)
# }
# # GROUP 4
# for (i in c(9, 10)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.08, offset.text = 0)
# }
# for (i in c(11)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.037, offset.text = 0)
# }
# for (i in c(23, 28, 18)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.082, offset.text = 0)
# }
# for (i in c(29, 26)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.085, offset.text = 0)
# }
# for (i in c(31)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.08, offset.text = 0)
# }
# for (i in c(36)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.082, offset.text = 0)
# }
# for (i in c(37)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.095, offset.text = 0)
# }
# for (i in c(38)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.105, offset.text = 0)
# }
# for (i in c(39, 40)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.049, offset.text = 0)
# }
# for (i in c(42)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.085, offset.text = 0)
# }
# for (i in c(59, 60, 53, 54, 57)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.092, offset.text = 0)
# }
# for (i in c(69)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.091, offset.text = 0)
# }
# for (i in c(87, 88, 90)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.091, offset.text = 0)
# }
# for (i in c(76)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.091, offset.text = 0)
# }
# for (i in c(75)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.081, offset.text = 0)
# }
# # GROUP 5
# for (i in c(13, 14)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.095, offset.text = 0)
# }
# for (i in c(15)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.055, offset.text = 0)
# }
# for (i in c(24)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.096, offset.text = 0)
# }
# for (i in c(34, 32)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.1, offset.text = 0)
# }
# for (i in c(43)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.1, offset.text = 0)
# }
# for (i in c(48)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.121, offset.text = 0)
# }
# for (i in c(58, 55)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.1, offset.text = 0)
# }
# for (i in c(51, 52)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.065, offset.text = 0)
# }
# for (i in c(74, 77)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.1, offset.text = 0)
# }
# for (i in c(78)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.096, offset.text = 0)
# }
# for (i in c(79)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.11, offset.text = 0)
# }
# for (i in c(86, 89, 91, 95)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.11, offset.text = 0)
# }
# # GROUP 6
# for (i in c(19, 16, 33)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.12, offset.text = 0)
# }
# for (i in c(46)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.12, offset.text = 0)
# }
# for (i in c(49)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.1, offset.text = 0)
# }
# for (i in c(50)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.09, offset.text = 0)
# }
# for (i in c(56)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.12, offset.text = 0)
# }
# for (i in c(61)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.125, offset.text = 0)
# }
# for (i in c(62)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.12, offset.text = 0)
# }
# for (i in c(65)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.115, offset.text = 0)
# }
# for (i in c(80)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.14, offset.text = 0)
# }
# for (i in c(92, 93)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.125, offset.text = 0)
# }
# for (i in c(84, 85)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.097, offset.text = 0)
# }
# # GROUP 7
# for (i in c(25, 21)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.135, offset.text = 0)
# }
# for (i in c(44)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.135, offset.text = 0)
# }
# for (i in c(66)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.135, offset.text = 0)
# }
# for (i in c(94)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.137, offset.text = 0)
# }
# for (i in c(83)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.114, offset.text = 0)
# }
# for (i in c(81)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.118, offset.text = 0)
# }
# for (i in c(82)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.125, offset.text = 0)
# }
# # GROUP 8
# for (i in c(30, 35)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.15, offset.text = 0)
# }
# for (i in c(63)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.15, offset.text = 0)
# }
# for (i in c(72, 73)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.15, offset.text = 0)
# }
# # GROUP 9
# for (i in c(45)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.17, offset.text = 0)
# }
# for (i in c(70, 71)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.17, offset.text = 0)
# }
# # GROUP 10
# for (i in c(64)) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.19, offset.text = 0)
# }
#-------------------------------------------------------------------------------
for (i in c(1:59)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.01*i, offset.text = 0)
}
plot_tree
flip(plot_tree, 872, 571)
# Plot with everything on it!

ggsave("figures/Lineageplot_tree_N.png", 
       plot = last_plot(),
       height = 15, width = 30)
# Save it

identify(plot_tree)


