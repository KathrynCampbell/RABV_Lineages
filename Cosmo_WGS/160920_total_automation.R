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
library(ips)
library(adephylo)


#############################################
#          SOURCE THE FUNCTION              #
#############################################

source("R/lineage_assignment.r")
source("R/lineage_naming.R")


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

sequence_data <- lineage_assignment(tree, min.support = 70, max.support = 100, alignment, metadata, ancestral)[[2]]
node_data <- lineage_assignment(tree, min.support = 70, max.support = 100, alignment, metadata, ancestral)[[1]]

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

node_data<-lineage_naming(sequence_data, metadata, node_data, tree)

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
plot_tree<-ggtree(tree, colour = "grey50", ladderize = T) %<+% sequence_data +
  geom_tippoint(aes(color=cluster), size=3)  +
  theme(legend.position = c(0.9, 0.15),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 20)) +
  guides(colour=guide_legend(override.aes=list(alpha=1, size=7))) +
  ggtitle(paste(args, "Lineage Tree", sep = " "))+
  theme(plot.title = element_text(size = 40, face = "bold"))

# + scale_color_manual(values=c("yellow","yellow2", "yellow3",
#                               "steelblue1", "steelblue3",
#                               "lawngreen", "green", "green2", "green3","green4",
#                               "olivedrab1","olivedrab2","olivedrab3","olivedrab4","darkolivegreen","darkgreen",
#                               "palegreen","palegreen3","chartreuse","chartreuse3","chartreuse4","springgreen",
#                               "blue","blue2","blue3","blue4",
#                               "lightcoral",
#                               "purple","purple3","purple4","darkorchid","darkorchid1",
#                               "gold","gold3","gold4",
#                               "seagreen1",
#                               "cyan","cyan3","cyan4","darkgreen",
#                               "firebrick4",
#                               "tan2",
#                               "navy",
#                               "violetred","violetred1","violetred2","violetred3","violetred4","deeppink4",
#                               "hotpink","hotpink3","hotpink4",
#                               "magenta","magenta2","magenta3","magenta4","deeppink",
#                               "red","red3","red4","tomato","tomato3","tomato4",
#                               "chocolate1","chocolate3","chocolate4",
#                               "khaki3",
#                               "grey68"))

# Plot each clade bar
# ---------------------------------------------------------------------------------------------
node_data$group<-NA

for (i in 1:length(node_data$Node)) {
  node_data$temp[i]<-length(which(node_data$Node %in% (descendants(tree, node_data$Node[i], type = "all"))))
}

node_data$group[1]<-1

for (i in 1:50) {
  group<-which(node_data$group == i)
  if (length(group != 0)) {
    for (j in 1:length(group)) {
      children<-which(node_data$Node %in% descendants(tree, (node_data$Node[group[j]]), type = "all", ignore.tip = T))
      node_data$group[children]<-i+1
    }
  }
}

for (i in 1:length(node_data$Node)) {
  plot_tree <-plot_tree +
    geom_cladelabel(node_data$Node[i], node_data$cluster[i], offset = 0.2*node_data$group[i], offset.text = 0, fontsize = 5)
}

plot_tree
# Plot with everything on it!


ggsave(paste(args, "/Figures/", args, "_lineage_tree.png", sep = ""), 
       plot = plot_tree,
       height = 20, width = 40)
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
# Add another column listing the number of sequences assigned to each cluster
clusters$n_seqs<-(sequence_data %>%
                    group_by(cluster)%>%
                    summarise(n=n()))$n

clusters<-clusters[-c(which(is.na(clusters$cluster))),]
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

# For each lineage, calculate the pairwise distance for all the sequences allocated to each lineage
# Extract the mean and max distance for each lineage

distances<-as.matrix(distTips(tree, tips = "all", method = "patristic"))

for (i in 1:length(node_data$Node)) {
  numbers<-Descendants(tree, node_data$Node[i], type = "tips")
  subset<-distances[,which(rownames(distances) %in% alignment$nam[numbers[[1]]])]
  
  clusters$max_distance[i] <- max(subset)
  clusters$mean_distance[i] <- mean(subset)
}

write.csv(clusters, file = (paste(args, "/Outputs/", args, "_lineage_info.csv", sep = "")), row.names=F)

