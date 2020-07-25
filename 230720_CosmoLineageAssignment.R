#'---
#'title: Cosmopolitan Lineage Assignment
#'author: Kathryn Campbell
#'date: 22/07/2020
#'---

rm(list = ls())
library(seqinr)
library(ape)
library(phytools)
library(treeio)
library(dplyr)
library(TreeTools)
library(adephylo)
library(phangorn)
library(phylotate)
library(caper)
library(ggtree)

#Import the tree
Cosmotree<-read_annotated(file="Trees/230720_Cosmo_copy.nex.txt")

#Sequence names got messed up in either MAFFT or IQTREE, need to fix these
Cosmotree$tip.label<-sub("(?<=\\.).*$", "", Cosmotree$tip.label, perl = T)
Cosmotree$tip.label<-gsub("\\.", "", Cosmotree$tip.label, perl = T)

#Import the metadata
Cosmometa<-read.csv("Sequences/220720_GLUE_CosmoMeta.csv")
#Plot the tree
# NEED TO COME BACK AND MAKE THIS NICER
ggtree(Cosmotree) +
  geom_nodelab()

#Identify nodes with a bootstrap of over 70
nodes_70<-which(Cosmotree$node.comment > 70 | Cosmotree$node.comment == 100); nodes_70

m<-matrix(ncol=2, nrow=500)
node_data<-data.frame(m)
node_data[,1]<-nodes_70
names(node_data)<-c("Node", "n tips")
# Make a dataframe ready for values to be put in
# Fill the first column with the numbers of the nodes identified in the previous steps

for(i in 1:500) {
  node_data[i,2]<-length(Descendants(Cosmotree, (nodes_70[i]), type = c("tips"))[[1]])
}

View(node_data)
# For each node identified in the previous step, count the number of tips descended from that node

nodes_5<-node_data[(which(node_data[,2]>=5)),]
# Only carry forwards nodes which have more than 5 tips descended from it
# This has been identified as the definition for a cluster in previous studies

Cosmometa$sequence.sequenceID[which(Cosmometa$sequence.gb_length < 11400)]
# In the metadata, do an initial sweep for any with less than 95% coverage at whole genome level 
# (whole genome c12kb, therefore anything less than 11.4kb can't have 95% coverage)
# Seqs with less than 95% coverage  = "AB839169" "MN726813" "MN726834" "MN726806" "MN726828" "MN726803" "MN726811" "MN726822" "MN726823" "MN726832"
which(Cosmotree$tip.label == "MN726832")
# These correspond to tip numbers 544, 76, 69, 73, 79, 114, 43, 45, 68, 44

r1<-Ancestors(Cosmotree, c(544), 'all')
r2<-Ancestors(Cosmotree, c(76), 'all')
r3<-Ancestors(Cosmotree, c(69), 'all')
r4<-Ancestors(Cosmotree, c(73), 'all')
r5<-Ancestors(Cosmotree, c(79), 'all')
r6<-Ancestors(Cosmotree, c(114), 'all')
r7<-Ancestors(Cosmotree, c(43), 'all')
r8<-Ancestors(Cosmotree, c(45), 'all')
r9<-Ancestors(Cosmotree, c(68), 'all')
r10<-Ancestors(Cosmotree, c(44), 'all')
# List the ancestor nodes for each of these tip numbers

removes<-c(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10)
# Turn the lists of nodes into one long list 

test<-data.frame(table(removes))
# Make a table to count the number the removed sequences descended from each node (e.g. for the deeper nodes, all 10 are descended)

names(test)<-c('Node', 'freq')
# Change the names

which(test[,1] %in% nodes_5[,1])
# Not all nodes are included in the nodes_5 data (some are already excluded) 
#which ones identified as ancestors of the sequences to re

newtest<-test[-c(1, 27, 28, 32, 33, 34, 35, 36, 39, 40, 41, 46, 47, 48),]; newtest
# Get rid of the nodes not in the nodes_5 data 

newtest[,1]
numbers<-c(569,  570,  571,  572,  573,  574,  575,  576,  577,  578,  579,  580,  581,  582,  583,  584,  585,  586,  587,  588,  589,  590,  591,  592,  593,  643,644,  659,  666,  667,  1105, 1106, 1107, 1115)
# Get a list of the nodes needing to be changed; this seemed the easiest way to do it as everything else caused errors!

nodes_new<-nodes_5

for (i in numbers) {
  nodes_new[which(nodes_new == i), 2]<-(nodes_5[which(nodes_5 == i), 2] - (newtest[which(newtest == i), 2]))
}
# Take away the number of removed tips from the previous total number of tips calculated for each node

nodes_new
# Works!!! Checked a few manually 

nodes_5<-nodes_new[(which(nodes_new[,2]>=5)),]
# Redo this to remove any that now have less than 5, and write over the old nodes_5 so this is updated with the new tip numbers

