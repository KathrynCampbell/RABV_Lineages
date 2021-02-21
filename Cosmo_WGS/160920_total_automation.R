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
library(dplyr)
library(phangorn)
library(caper)
library(stringr)
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

# Need to edit column names so they match what is required
metadata<- metadata %>%
  rename(ID = sequence.sequenceID,
         year = sequence.latest_collection_year,
         country = sequence.m49_country.display_name)

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
#'========================================================================================================
#' Function in 'R/lineage_assignment.r'
#' This function will split the tree into lineages, and assign a lineage number to each sequence
#' It may take a little while!
#' 
#' The lineages are defined according to Rambaut et al. (2020) in which there must be:
#'    - 70% support at the defining node
#'    - At least 10 genomes with 95% coverage
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
#'    - If using ultrafast bootstrap support, use min.support = 95 and max.support = 100
#'    - If using posterior support, use min.support = 0.7 and max.support = 1.0
#' The minimum number of sequences to define a lineage can also be changed in the sequences argument
#' 
#'=========================================================================================================

sequence_data <- lineage_assignment(tree, min.support = 95, max.support = 100, alignment, metadata, ancestral)[[2]]
node_data <- lineage_assignment(tree, min.support = 95, max.support = 100, alignment, metadata, ancestral)[[1]]

#---------------------------------------------------------------------------------------
#
# Everything above this is part of the lineage assignment script
# Everything below is added extras specific to the cosmopolitan lineages - 
# Makes a plot with informative naming
#
#-------------------------------------------------------------------------------

# #############################################
#           RENAME THE LINEAGES               #
# #############################################

# Add a column of previous assignment to the table for comparison
# sequence_data$previous <- NA
# for (i in 1:length(sequence_data$ID)) {
#   sequence_data$previous[i]<-
#     metadata$alignment.displayName[which(metadata$ID == sequence_data$ID[i])]
# }
# 
# # Plot a nice figure to save
# plot_tree<-ggtree(tree) %<+% sequence_data +
#   geom_tippoint(na.rm = T, aes(colour = (cluster))) +
#   ggtitle(paste(args, "Lineage Tree", sep = ""))
# 
# # Plot each clade bar
# for (i in c(1:length(node_data$Node))) {
#   plot_tree<-plot_tree +
#     geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.5+0.15*i, offset.text = 0)
# }
# 
# plot_tree
# 
# # Just doing this by eye currently, must be a better way
# # ------------------------------
node_data$cluster[1]<-"A1"
node_data$cluster[2]<-"A1.1"
node_data$cluster[3]<-"A1.1.1"
node_data$cluster[4]<-"B1"
node_data$cluster[5]<-"C1"
node_data$cluster[6]<-"AF1b_A1"
node_data$cluster[7]<-"AF1b_A1.1"
node_data$cluster[8]<-"AF1b_A1.1.1"
node_data$cluster[9]<-"B1.1"
node_data$cluster[10]<-"B1.2"
node_data$cluster[11]<-"ME1a_A1"
node_data$cluster[12]<-"ME1a_A1.1"
node_data$cluster[13]<-"AF1b_B1"
node_data$cluster[14]<-"ME1a_A1.1.1"
node_data$cluster[15]<-"ME1a_B1"
node_data$cluster[16]<-"ME2_A1"
node_data$cluster[17]<-"D1"
node_data$cluster[18]<-"ME1a_B1.1"
node_data$cluster[19]<-"ME2_A1.1"
node_data$cluster[20]<-"D1.1"
node_data$cluster[21]<-"B1.1.1"
node_data$cluster[22]<-"ME1a_B1.1.1"
node_data$cluster[23]<-"ME2_A1.1.1"
node_data$cluster[24]<-"E1"
node_data$cluster[25]<-"ME1a_C1"
node_data$cluster[26]<-"ME2_B1"
node_data$cluster[27]<-"AF1b_C1"
node_data$cluster[28]<-"AF1b_B1.1"
node_data$cluster[29]<-"NEE_A1"
node_data$cluster[30]<-"CA1_A1"
node_data$cluster[31]<-"ME1a_C1.1"
node_data$cluster[32]<-"CA2_A1"
node_data$cluster[33]<-"AF1b_B1.1.1"
node_data$cluster[34]<-"NEE_A1.1"
node_data$cluster[35]<-"CA1_A1.1"
node_data$cluster[36]<-"CA2_A1.1"
node_data$cluster[37]<-"AF1b_A1.2"
node_data$cluster[38]<-"AF1b_D1"
node_data$cluster[39]<-"AF1b_B1.2"
node_data$cluster[40]<-"AF1a_A1"
node_data$cluster[41]<-"D1.1.1"
node_data$cluster[42]<-"CA1_A1.1.1"
node_data$cluster[43]<-"ME1a_C1.1.1"
node_data$cluster[44]<-"ME2_B1.1"
node_data$cluster[45]<-"AF1b_D1.1"
node_data$cluster[46]<-"AF1b_B1.2.1"
node_data$cluster[47]<-"F1"
node_data$cluster[48]<-"AM2a_A1"
node_data$cluster[49]<-"D1.1.2"
node_data$cluster[50]<-"WE_A1"
node_data$cluster[51]<-"CE_A1"
node_data$cluster[52]<-"EE_A1"
node_data$cluster[53]<-"NEE_A1.1.1"
node_data$cluster[54]<-"CA1_B1"
node_data$cluster[55]<-"ME1a_A1.1.2"
node_data$cluster[56]<-"ME1a_D1"
node_data$cluster[57]<-"ME1a_C1.1.2"
node_data$cluster[58]<-"ME1a_A1.1.3"
node_data$cluster[59]<-"ME1a_A1.1.4"
node_data$cluster[60]<-"ME2_B1.1.1"
node_data$cluster[61]<-"CA2_A1.1.1"
node_data$cluster[62]<-"AF1b_A1.2.1"
node_data$cluster[63]<-"AF1b_C1.1"
node_data$cluster[64]<-"AF1b_D1.1.1"
node_data$cluster[65]<-"AF1b_B1.3"
node_data$cluster[66]<-"AF1b_E1"
node_data$cluster[67]<-"AF1b_A1.2"
node_data$cluster[68]<-"AF1a_A1.1"

#-------------------------------------------------------------------------------
# # Renamed according to Rambaut et al (2020) with A1.1.1.1 becoming a new letter (e.g. C1)
# 
for (i in 1:length(node_data$cluster)) {
  sequence_data$cluster[which(sequence_data$cluster == i)] <- node_data$cluster[i]
}
# Rename the lineages in the sequence assignment table

#############################################
#                 WGS PLOT                  #
#############################################
attach(tree)
sequence_data$cluster <- as.factor(sequence_data$cluster)

# Plot a nice figure to save
plot_tree<-ggtree(tree) %<+% sequence_data +
  geom_tippoint(na.rm = T, aes(colour = (cluster))) +
  theme(legend.position = c(0.1, 0.83),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8)) +
  guides(colour=guide_legend(override.aes=list(alpha=1, size=5))) +
  ggtitle(paste(args, "Lineage Tree", sep = " "))

# Plot each clade bar
# ---------------------------------------------------------------------------------------------
# GROUP 1
for (i in c(1)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0, offset.text = 0)
}
# GROUP 2
for (i in c(2)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.01*i, offset.text = 0)
}
# GROUP 3
for (i in c(3)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.04, offset.text = 0)
}
# GROUP 4
for (i in c(17)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.055, offset.text = 0)
}
for (i in c(4)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.12, offset.text = 0)
}
for (i in c(32,5)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.13, offset.text = 0)
}
for (i in c(40)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.14, offset.text = 0)
}

# GROUP 5
for (i in c(20)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.07, offset.text = 0)
}
for (i in c(10,9,16,6)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.14, offset.text = 0)
}
for (i in c(36,68)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.15, offset.text = 0)
}
# GROUP 6
for (i in c(11,50,51,52,29,21,19,61)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.165, offset.text = 0)
}
for (i in c(41,48,49)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.12, offset.text = 0)
}
for (i in c(7,37,67)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.16, offset.text = 0)
}
# GROUP 7
for (i in c(12,34,24,23,8,62)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.2, offset.text = 0)
}
for (i in c(47)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.16, offset.text = 0)
}
# GROUP 8
for (i in c(55,58,59,14,53,30,26,13,27)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.22, offset.text = 0)
}
# GROUP 9
for (i in c(15,35,44,28,39,63,65)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.25, offset.text = 0)
}
# GROUP 10
for (i in c(18,42,60,33,46)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.27, offset.text = 0)
}
# GROUP 11
for (i in c(22,54,38,66)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.29, offset.text = 0)
}
# GROUP 12
for (i in c(25,45)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.31, offset.text = 0)
}
# GROUP 13
for (i in c(31,64)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.33, offset.text = 0)
}
# GROUP 14
for (i in c(57,43)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.35, offset.text = 0)
}
# GROUP 15
for (i in c(56)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.37, offset.text = 0)
}
plot_tree
# Plot with everything on it!


ggsave(paste(args, "/Figures/", args, "_lineage_tree.png", sep = ""), 
       plot = plot_tree,
       height = 20, width = 30)
# Save it

#KB added row.names=F to avoid a column of row numbers
write.csv(sequence_data, file = (paste(args, "/Outputs/", args, "_sequence_data.csv", sep = "")), row.names=F)
write.csv(node_data, file = (paste(args, "/Outputs/", args, "_node_data.csv", sep = "")), row.names=F)


#############################################
#                TIMESERIES                 #
#############################################
# Plot a timeseries of sequences per year split by cluster
timeseries<- ggplot(sequence_data, aes(x=factor(Year), fill=cluster))+
  geom_histogram(stat="count") +
  ggtitle(paste(args, "Time Series", sep = " "));timeseries

# Save it
ggsave(paste(args, "/Figures/", args, "_timeseries.png", sep = ""),
       plot = timeseries,
       height = 20, width = 30)


#############################################
#         LINEAGE INFORMATION TABLE         #
#############################################

# Extract the place information from the metadata
for (i in 1:length(sequence_data$ID)) {
  sequence_data$Country[i] <- metadata$country[(which(metadata$ID == sequence_data$ID[i]))]
}

# Create a data frame ready to fill in information about each cluster
clusters <- sequence_data %>%
  group_by(cluster) %>%
  summarise()
clusters$country <- NA
clusters$year_first <- NA
clusters$year_last <- NA
clusters$max_distance <- NA
clusters$mean_distance <- NA

# For each cluster, find and list the earliest collection year, the latest collection year and all the places
# that cluster has been found
for (i in 1:length(clusters$cluster)) {
  clusters$year_first[i] <- sequence_data %>%
    filter(cluster == clusters$cluster[i])%>%
    group_by(Year)%>%
    summarise()%>%
    min()
  clusters$year_last[i] <- sequence_data %>%
    filter(cluster == clusters$cluster[i])%>%
    group_by(Year)%>%
    summarise()%>%
    max()
  clusters$country[i]<-paste((sequence_data %>%
                                filter(cluster == clusters$cluster[i] & Country !="-") %>%
                                group_by(Country) %>%
                                summarise()))
}

# Add another column listing the number of sequences assigned to each cluster
clusters$n_seqs<-(sequence_data %>%
                    group_by(cluster)%>%
                    summarise(n=n()))$n

# For each lineage, calculate the pairwise distance for all the sequences allocated to each lineage
# Extract the mean and max distance for each lineage
for (i in 1:length(clusters$cluster)) {
  numbers<-which(alignment$nam %in% (sequence_data$ID[which(sequence_data$cluster == clusters$cluster[i])]))
  
  test_align <- as.alignment(alignment$seq[c(numbers)])
  test_align$nb <- length(numbers)
  test_align$nam <- alignment$nam[c(numbers)]
  test_align$seq <- alignment$seq[c(numbers)]
  test_align$com <- "NA"
  
  test_align <- as.matrix(test_align)
  distances <- as.matrix(dist.gene(test_align, method = "pairwise", pairwise.deletion = F, variance = F))
  clusters$max_distance[i] <- max(distances)
  clusters$mean_distance[i] <- mean(distances)
}

write.csv(clusters, file = (paste(args, "/Outputs/", args, "_lineage_info.csv", sep = "")), row.names=F)
