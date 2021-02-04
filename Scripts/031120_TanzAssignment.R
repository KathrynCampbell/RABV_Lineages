#'---------------------------------------------------------
#'title: Tanzanian Lineage Assignment
#'author: Kathryn Campbell
#'date: 02/11/2020
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
tree_N <- read_annotated(file="Trees/021120_GLUE_TanzSeqs_N_align.nex") # GLUE sequences
tree_WGS <- read_annotated(file="Trees/021120_GLUE_TanzSeqs_align.nex") # GLUE sequences

# Need to edit node.comment so it's in the format required 
tree_N$node.comment<-gsub(".*&bootstrap=", "", tree_N$node.comment)
tree_WGS$node.comment<-gsub(".*&bootstrap=", "", tree_WGS$node.comment)


#'**METADATA**
#'========================================================================================================
#' The metadata must contain the element 'year' which lists the collection year for each sequence 
#' And the element 'ID' which lists all the sequence ID's
#' These sequence ID's must match the sequence ID's in the tree and alignment
#'=========================================================================================================
metadata <- read.csv("Sequences/021120_GLUE_TanzMeta.csv")

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
alignment_N <- read.alignment("Sequences/021120_GLUE_TanzSeqs_N_align.fasta", format = "fasta")
alignment_WGS <- read.alignment("Sequences/021120_GLUE_TanzSeqs_align.fasta", format = "fasta")



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

sequence_data_N <- lineage_assignment(tree_N, min.support = 95, max.support = 100, alignment_N, metadata, sequences = 10)[[2]]
node_data_N <- lineage_assignment(tree_N, min.support = 95, max.support = 100, alignment_N, metadata, sequences = 10)[[1]]

sequence_data_WGS <- lineage_assignment(tree_WGS, min.support = 95, max.support = 100, alignment_WGS, metadata, sequences = 10)[[2]]
node_data_WGS <- lineage_assignment(tree_WGS, min.support = 95, max.support = 100, alignment_WGS, metadata, sequences = 10)[[1]]

#---------------------------------------------------------------------------------------
#
# Everything above this is part of the lineage assignment script
# Everything below is added extras specific to the Tanzanian lineages - 
# Makes a plot with informative naming
#
#---------------------------------------------------------------------------------------

#############################################
#            PLOT WGS LINEAGES              #
#############################################

node_data_WGS$cluster[1]<-"A1"
node_data_WGS$cluster[2]<-"A1.1"
node_data_WGS$cluster[3]<-"A1.1.1"
node_data_WGS$cluster[4]<-"B1"
node_data_WGS$cluster[5]<-"B1.1"
node_data_WGS$cluster[6]<-"B1.1.1"
node_data_WGS$cluster[7]<-"A1.1.2"
node_data_WGS$cluster[8]<-"C1"
node_data_WGS$cluster[9]<-"C1.1"
node_data_WGS$cluster[10]<-"C1.1.1"
node_data_WGS$cluster[11]<-"D1"
node_data_WGS$cluster[12]<-"E1"
node_data_WGS$cluster[13]<-"A1.1.3"

for (i in 1:length(node_data_WGS$cluster)) {
  sequence_data_WGS$cluster[which(sequence_data_WGS$cluster == i)] <- node_data_WGS$cluster[i]
}
# Rename the lineages in the sequence assignment table 
# Renamed according to Rambaut et al (2020) with A1.1.1.1 becoming a new letter (e.g. C1)

sequence_data_WGS$cluster <- as.factor(sequence_data_WGS$cluster)

# Plot a nice figure to save
plot_tree_WGS<-ggtree(tree_WGS) %<+% sequence_data_WGS +
  geom_tippoint(na.rm = T, aes(colour = (cluster))) +
  theme(legend.position = c(0.05, 0.83),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8)) +
  guides(colour=guide_legend(override.aes=list(alpha=1, size=5))) +
  ggtitle("WGS")

# Plot each clade bar
for (i in c(1:13)) {
    plot_tree_WGS <-plot_tree_WGS +
      geom_cladelabel(node_data_WGS$Node[i], node_data_WGS$cluster[i], offset = 0.01*i, offset.text = 0)
}

plot_tree_WGS


#############################################
#          PLOT N GENE LINEAGES             #
#############################################

sequence_data_N$cluster <- as.factor(sequence_data_N$cluster)

# Plot a nice figure to save
plot_tree_N<-ggtree(tree_N)  %<+% sequence_data_N +
  geom_tippoint(na.rm = T, aes(colour = (cluster))) +
  theme(legend.position = c(0.05, 0.83),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8)) +
  guides(colour=guide_legend(override.aes=list(alpha=1, size=5))) +
  ggtitle("N gene")


# For some reason the node_data_N table doesn't work, but can see from the tree the assignments
# Do these manually:

# identify(plot_tree_N)

# 1 = 207
# 2 = 208
# 3 = 325
# 4 = 345
# 5 = 210
# 6 = 346

node_data_N <- data.frame(Node = c(207, 209, 274, 345, 275, 352), cluster = c(1:6))

# Rename the lineages in the sequence assignment table 
# Renamed according to Rambaut et al (2020) with A1.1.1.1 becoming a new letter (e.g. C1)

node_data_N$cluster[1]<-"A1"
node_data_N$cluster[2]<-"A1.1"
node_data_N$cluster[3]<-"A1.1.1"
node_data_N$cluster[4]<-"A1.2"
node_data_N$cluster[5]<-"B1"
node_data_N$cluster[6]<-"A1.2.1"

node_data_N$cluster[1]

for (i in (as.factor(1:length(node_data_N$cluster)))) {
  sequence_data_N$cluster[which(sequence_data_N$cluster == i)] <- node_data_N$cluster[i]
}

# Plot each clade bar
for (i in c(1:6)) {
  plot_tree_N <-plot_tree_N +
    geom_cladelabel(node_data_N$Node[i], node_data_N$cluster[i], offset = 0.01*i, offset.text = 0)
}

plot_tree_N


#############################################
#         WGS AND N GENE COMPARISON         #
#############################################

combined <- grid.arrange(plot_tree_N, plot_tree_WGS, ncol = 2)

ggsave("figures/Lineageplot_tree_combined_Tanz.png", 
       plot = combined,
       height = 20, width = 30)
# Save it

for (i in 1:205) {
  sequence_data_WGS$"N_cluster"[i] <- sequence_data_N$cluster[which(sequence_data_N$ID == sequence_data_WGS$ID[i])]
}
# Add a column in the sequence data to compare N gene assignment - not working currently
sequence_data_N$cluster


#############################################
#                TIMESERIES                 #
#############################################

# Plot a timeseries of sequences per year split by cluster
timeseries_WGS<- ggplot(sequence_data_WGS, aes(x=factor(Year), fill=cluster))+
  geom_histogram(stat="count") +
  ggtitle("Tanz Time Series WGS"); timeseries_WGS

# Save it
ggsave("figures/TimeSeries_Tanz_WGS.png", 
       plot = timeseries_WGS,
       height = 20, width = 30)

# Currently not working for N gene - need to fix

#############################################
#         LINEAGE INFORMATION TABLE         #
#############################################

# Extract the place information from the metadata
for (i in 1:length(sequence_data_WGS$ID)) {
  sequence_data_WGS$Place[i] <- metadata$sequence.gb_place_sampled[(which(metadata$ID == sequence_data_WGS$ID[i]))]
}

# Create a data frame ready to fill in information about each cluster
clusters_WGS <- sequence_data_WGS %>%
  group_by(cluster) %>%
  summarise()
clusters_WGS$place <- NA

# For each cluster, find and list the earliest colleciton year, the latest collection year and all the places
# that cluster has been found
for (i in 1:length(clusters_WGS$cluster)) {
  clusters_WGS$year_first[i] <- sequence_data_WGS %>%
    filter(cluster == clusters_WGS$cluster[i])%>%
    group_by(Year)%>%
    summarise()%>%
    min()
  clusters_WGS$year_last[i] <- sequence_data_WGS %>%
    filter(cluster == clusters_WGS$cluster[i])%>%
    group_by(Year)%>%
    summarise()%>%
    max()
  clusters_WGS$place[i]<-paste((sequence_data_WGS %>%
                                  filter(cluster == clusters_WGS$cluster[i] & Place !="-") %>%
                                  group_by(Place) %>%
                                  summarise()))
}

# Add another column listing the number of sequences assigned to each cluster
clusters_WGS$n_seqs<-(sequence_data_WGS %>%
                    group_by(cluster)%>%
                    summarise(n=n()))$n

# For each lineage, calculate the pairwise distance for all the sequences allocated to each lineage
# Extract the mean and max distance for each lineage
for (i in 1:length(clusters_WGS$cluster)) {
  numbers<-which(alignment_WGS$nam %in% (sequence_data_WGS$ID[which(sequence_data_WGS$cluster == clusters_WGS$cluster[i])]))
  
  test_align <- as.alignment(alignment_WGS$seq[c(numbers)])
  test_align$nb <- length(numbers)
  test_align$nam <- alignment_WGS$nam[c(numbers)]
  test_align$seq <- alignment_WGS$seq[c(numbers)]
  test_align$com <- "NA"
  
  test_align <- as.matrix(test_align)
  distances <- as.matrix(dist.gene(test_align, method = "pairwise", pairwise.deletion = F, variance = F))
  clusters_WGS$max_distance[i] <- max(distances)
  clusters_WGS$mean_distance[i] <- mean(distances)
}