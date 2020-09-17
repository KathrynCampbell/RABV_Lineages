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

sequence_data$cluster <- as.factor(sequence_data$cluster)

attach(tree)

plot_tree<-ggtree(tree) %<+% sequence_data +
  geom_tippoint(na.rm = T, aes(colour = (cluster))) +
  theme(legend.position = c(0.12, 0.83),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8)) +
  guides(colour=guide_legend(override.aes=list(alpha=1, size=5)))


# for (i in c(1:length(sequence_data$ID))) {
#   plot_tree <-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.01*i)
# }


plot_tree


