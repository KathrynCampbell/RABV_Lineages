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

sequence_data <- lineage_assignment(tree, min.support = 70, max.support = 100, alignment, metadata)[[2]]
node_data <- lineage_assignment(tree, min.support = 70, max.support = 100, alignment, metadata)[[1]]

#---------------------------------------------------------------------------------------
# Compare with WGS assignment
# Different number of lineages - definitely different to WGS by 7 - not yet clear which 

sequence_data_all <- read.csv(file = "Outputs/sequence_data_cosmo.csv")
node_data_all <- read.csv(file = "Outputs/node_data_cosmo.csv")

for (i in 1:567) {
  sequence_data$WGS_cluster[i] <- sequence_data_all$cluster[which(sequence_data_all$ID == sequence_data$ID[i])]
}

