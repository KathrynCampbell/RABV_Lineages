#'---
#'title: Cosmopolitan Lineage Assignment
#'author: Kathryn Campbell
#'date: 22/07/2020
#'---

#KB added suggestions 08/10/20, 22/10/20
#KB: check out https://stackoverflow.com/questions/50604055/ggplot2-geom-bar-position-failure/50609038
## Fixed!
rm(list=ls())
library(ggplot2)
library(caper)
library(phylotate)
library(dplyr)
library(ggpubr)
#############################################
#            IMPORT THE DATA                #
#############################################
#KB added ignore row names (remove if you change csv export)
sequence_data <- read.csv("Outputs/sequence_data_cosmo.csv")
node_data <- read.csv("Outputs/node_data_cosmo.csv")
# These were created in the lineage assignment script (Scripts/160920_total_automation.R)

tree <- read_annotated(file="Trees/230720_Cosmo_copy.nex.txt") # GLUE sequences
# #KB- can replace above 2 lines with this:
tree$tip.label <- gsub("\\..*", "", tree$tip.label, perl = T)

attach(sequence_data)
#############################################
#            IMPORT THE DATA                #
#############################################
# Test some different time series resolutions
# 
# All lineages included
# -------------------------------------------

all<- ggplot(sequence_data, aes(x=factor(Year), fill=cluster))+
  geom_histogram(stat="count") +
  ggtitle("All"); all

# Major only
# -------------------------------------------
sequence_data_2 <- sequence_data
sequence_data_2$cluster <- gsub("\\..*","",sequence_data$cluster)
major <- ggplot(sequence_data_2, aes(x=factor(Year), fill=cluster))+
  geom_histogram(stat="count", position="stack") +
  ggtitle("Major"); major

# Split into groups as seen in figures/LineageTreeGroups.png
# -------------------------------------------

# These groups are the node numbers that correspond to the lineages in each group identified 
# -------------------------------------------
group1 <- c(1,46,65,47,115)
group2 <- c(3,84,85,54,41,2,22,64,113,66)
group3 <- c(4,86,55,6,5,20,109,62,26,112,114)
group4 <- c(9,23,10,81,82,83,87,7,58,99,100,101,44,25,18,36,108,37)
group5 <- c(27,68,21,50,73,51,75,11,53,8,98,59,29,19,107,63)
group6 <- c(31,49,70,71,72,73,74,75,12,53,86,87,13,42,43,95,57,102,45,30,105,106,110,111,64)
group7 <- c(38,69,32,33,39,28,80,14,94,97,60,35)
group8 <- c(48,76,52,78,34,15,103,61)
group9 <- c(77,40,16,104)
group10 <- c(67,79,17)
group11 <- c(24,91,92)
group12 <- c(88,89,90,93)
# -------------------------------------------

# RUN FROM HERE EACH TIME TESTING NEW GROUPS - needs this first bit too
# Create duplicate data frames so the original ones can be used again later
sequence_data_2 <- sequence_data
node_data_2 <- node_data
node_data_2 <- node_data[-c(group12, group11, group10, group9, group8, group7),]
# Remove the nodes of the lineages discarded as options for the later assignment

# Order the nodes of interest by the number of times they overlap the other nodes of interest (descending)
node_data_2 <- node_data_2[order(-node_data_2$overlaps),]

for (i in 1:(length(node_data_2$Node))) {
  sequence_data_2$cluster[which(sequence_data_2$ID %in% clade.members(node_data_2$Node[i], tree, include.nodes = F, tip.labels = T))] <- node_data_2$cluster[i]
}
# Reassign the lineages with the reduced number of options

groups6 <- ggplot(sequence_data_2, aes(x=Year, fill=cluster))+
  geom_histogram(stat="count") +
  #KB edit: can add in text angle to make axes easier to read:
  theme(axis.text.x = element_text(angle = 45)) +
  ggtitle("Groups_6"); groups6


ggarrange(all, major, groups6, legend = F, ncol=2, nrow=2)

