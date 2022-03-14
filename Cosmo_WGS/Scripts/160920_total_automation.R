#!/usr/bin/env Rscript

rm(list=ls())

args = "Cosmo_WGS"

#'---------------------------------------------------------
#'title: Cosmopolitan WGS Lineage Assignment
#'author: Kathryn Campbell
#'date: 16/09/2020
#'---------------------------------------------------------

#############################################
#            INSTALL PACKAGES               #
#############################################
library(seqinr)
library(ape)
library(MADDOG)
library(ggplot2)

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
tree <- read.tree(file = paste(args, "/Trees/", args, "_aligned.fasta.contree", sep = ""))
# Sequence names got messed up in MAFFT, need to fix these so they match metadata and alignment
# Also node comment is sometimes weird, fix it
# #KB- can replace above 2 lines with this:
tree$tip.label <- gsub("\\..*", "", tree$tip.label, perl = T)
tree$node.comment<- gsub(".*=", "", tree$node.label, perl = T)

#'**METADATA**
#'========================================================================================================
#' The metadata must contain the element 'year' which lists the collection year for each sequence 
#' And the element 'ID' which lists all the sequence ID's
#' These sequence ID's must match the sequence ID's in the tree and alignment
#'=========================================================================================================
metadata <- read.csv(file = paste(args, "/", args, "_metadata.csv", sep = ""))

#'**ALIGNMENT**
#'========================================================================================================
#' The alignment must contain the element 'seq' which contains the sequences
#' And the element 'nam' which lists all the sequence ID's
#' These sequence ID's must match the sequence ID's in the metadata and tree
#'=========================================================================================================
alignment <- read.alignment(file = (paste(args, "/Alignment/", args, "_aligned.fasta", sep = "")), format = "fasta")

# Sequence names got messed up in MAFFT, need to fix these so they match metadata and alignment
# #KB- can replace above 2 lines with this:
alignment$nam <- gsub("\\..*", "", alignment$nam, perl = T)


#'**TIMETREE**
#'========================================================================================================
#' 
#'=========================================================================================================
ancestral <- read.alignment(file = (paste(args, "/Timetree/ancestral_sequences.fasta", sep = "")), format = "fasta")
ancestral$nam <- gsub("\\..*", "", ancestral$nam, perl = T)

#############################################
#           RUN THE ASSIGNMENT              #
#############################################

node_data<-MADDOG::node_info(tree, 70, alignment, metadata, ancestral)
sequence_data<-MADDOG::seq_designation(tree, 70, alignment, metadata, ancestral)
lineage_info<-MADDOG::lineage_info(sequence_data, metadata)
#---------------------------------------------------------------------------------------

write.csv(sequence_data, file = (paste(args, "/Outputs/", args, "_sequence_data.csv", sep = "")), row.names=F)
write.csv(node_data, file = (paste(args, "/Outputs/", args, "_node_data.csv", sep = "")), row.names=F)
write.csv(lineage_info, file = (paste(args, "/Outputs/", args, "_lineage_info.csv", sep = "")), row.names=F)

#############################################
#                 WGS PLOT                  #
#############################################
plot_tree<-MADDOG::lineage_tree(lineage_info, node_data, tree, metadata, sequence_data)

ggsave(paste(args, "/Figures/", args, "_lineage_tree.png", sep = ""),
       plot = plot_tree)


sunburst<-MADDOG::sunburst(lineage_info, node_data, tree, metadata, sequence_data)

htmlwidgets::saveWidget(plotly::as_widget(sunburst), paste(args, "/Figures/", args, "_sunburst.html", sep = ""))

map<-MADDOG::lineage_map(lineage_info, node_data, tree, metadata, sequence_data)

ggsave(paste(args, "/Figures/", args, "_lineage_map.png", sep = ""),
       plot = map)
