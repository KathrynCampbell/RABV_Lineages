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
library(stringr)





#############################################
#            IMPORT THE DATA                #
#############################################

#Import the tree
Cosmotree<-read_annotated(file="Trees/230720_Cosmo_copy.nex.txt")

#Sequence names got messed up in MAFFT, need to fix these
Cosmotree$tip.label<-sub("(?<=\\.).*$", "", Cosmotree$tip.label, perl = T)
Cosmotree$tip.label<-gsub("\\.", "", Cosmotree$tip.label, perl = T)

#Import the metadata
Cosmometa<-read.csv("Sequences/220720_GLUE_CosmoMeta.csv")

#Import the alignment
Cosmoalign<-read.alignment("Sequences/220720_GLUE_CosmoSeqs_align.fasta", format = "fasta")

#Sequence names got messed up in MAFFT, need to fix these
Cosmoalign$nam<-sub("(?<=\\.).*$", "", Cosmotree$tip.label, perl = T)
Cosmoalign$nam<-gsub("\\.", "", Cosmotree$tip.label, perl = T)




#############################################
#            PLOT THE TREE                  #
#############################################

# NEED TO COME BACK AND MAKE THIS NICER
ggtree(Cosmotree) +
  geom_nodelab()





#############################################
#            BOOTSTRAP SUPPORT              #
#############################################

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





#############################################
#            95% COVERAGE WGS               #
#############################################

# Make a dataframe ready to fill with info about number of gaps and N bases
m<-matrix(ncol=5, nrow=567)
seq_data<-data.frame(m)
names(seq_data)<-c("ID", "N", "-", "Length_before", "Length_after")
seq_data$ID<-Cosmoalign$nam
seq_data$Length_before<-nchar(Cosmoalign$seq[[1]])
# Add a column with the length of the alignment

for (i in 1:567) {
  seq_data$N[i]<-str_count(Cosmoalign$seq[[i]], pattern = 'n')
  seq_data$`-`[i]<-str_count(Cosmoalign$seq[[i]], pattern = '-')
  seq_data$Length_after[i]<-(seq_data$Length_before[i] - seq_data$N[i] - seq_data$`-`[i])
}
# For each sequence, count the number of n bases and the number of gaps
# Calculate the length after removing these 

seq_data$ID[which(seq_data$Length_after < (seq_data$Length_before * 0.95))]
# Seqs with less than 95% coverage  =  "KY210300" "KX148205" "MK760768" "KF154998" "DQ875050" "LC325820" "DQ099525" "LT839616" 
# "JQ944709" "AB839170" "GU565703" "KX148110" "KX148111" "KX148102" "HQ450386" "KX708502" "JQ685970" "KX708499" "KX708501" "KX708500"
which(Cosmotree$tip.label == "KX708500")
# These correspond to tip numbers 9, 141, 534, 535, 538, 539, 541, 542, 543, 545, 549, 553, 554, 555, 560, 561, 564, 565, 566, 567

r1<-Ancestors(Cosmotree, c(9), 'all')
r2<-Ancestors(Cosmotree, c(141), 'all')
r3<-Ancestors(Cosmotree, c(534), 'all')
r4<-Ancestors(Cosmotree, c(535), 'all')
r5<-Ancestors(Cosmotree, c(538), 'all')
r6<-Ancestors(Cosmotree, c(539), 'all')
r7<-Ancestors(Cosmotree, c(541), 'all')
r8<-Ancestors(Cosmotree, c(542), 'all')
r9<-Ancestors(Cosmotree, c(543), 'all')
r10<-Ancestors(Cosmotree, c(545), 'all')
r11<-Ancestors(Cosmotree, c(549), 'all')
r12<-Ancestors(Cosmotree, c(553), 'all')
r13<-Ancestors(Cosmotree, c(554), 'all')
r14<-Ancestors(Cosmotree, c(555), 'all')
r15<-Ancestors(Cosmotree, c(560), 'all')
r16<-Ancestors(Cosmotree, c(561), 'all')
r17<-Ancestors(Cosmotree, c(564), 'all')
r18<-Ancestors(Cosmotree, c(565), 'all')
r19<-Ancestors(Cosmotree, c(566), 'all')
r20<-Ancestors(Cosmotree, c(567), 'all')

# List the ancestor nodes for each of these tip numbers

removes<-c(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, r16, r17, r18, r19, r20)
# Turn the lists of nodes into one long list 

remove_counts<-data.frame(table(removes))
# Make a table to count the number the removed sequences descended from each node (e.g. for the deeper nodes, all 10 are descended)

names(remove_counts)<-c('Node', 'freq')
# Change the names

which(remove_counts[,1] %in% nodes_5[,1])
# Not all nodes are included in the nodes_5 data (some are already excluded) 
#which ones identified as ancestors of the sequences to re

new_remove<-remove_counts[-c(1, 29, 30, 31, 33, 34, 35, 36, 37, 38, 39, 42, 43, 49, 50, 51, 52, 54, 55, 56, 57, 58, 62, 63, 64, 65, 67, 68),]; new_remove
# Get rid of the nodes not in the nodes_5 data 

more(new_remove[,1])
numbers<-c(569, 570,  571,  572,  573,  574,  575,  576,  577,  578,  579,  580,  581,  582,  583,  584,  585,  586,  587,  588,
           589,  590,  591, 592,  593,  594, 595,  599,  715,  1098, 1105, 1106, 1107, 1108, 1109, 1115, 1122, 1123, 1124, 1131)
# Get a list of the nodes needing to be changed; this seemed the easiest way to do it as everything else caused errors!

nodes_new<-nodes_5

for (i in numbers) {
  nodes_new[which(nodes_new == i), 2]<-(nodes_5[which(nodes_5 == i), 2] - (new_remove[which(new_remove == i), 2]))
}
# Take away the number of removed tips from the previous total number of tips calculated for each node

nodes_new
# Works!!! Checked a few manually 

nodes_5<-nodes_new[(which(nodes_new[,2]>=5)),]
# Redo this to remove any that now have less than 5, and write over the old nodes_5 so this is updated with the new tip numbers





#############################################
#         DIFFERENCE FROM ANCESTOR          #
#############################################

rm(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, r16, r17, r18, r19, r20)

nodes_5<-nodes_new[(which(nodes_new[,2]<=250)),]

help("geom_cladelabel")

ggtree(Cosmotree) + 
  geom_cladelabel(647, "647") +
  geom_cladelabel(651, "651", offset = 0.005) +
  geom_cladelabel(655, "655", offset = 0.005) +
  geom_cladelabel(667, "667", offset = 0.005) +
  geom_cladelabel(675, "675", offset = 0.005) +
  geom_cladelabel(717, "717", offset = 0.005) +
  geom_cladelabel(736, "736", offset = 0.005) +
  geom_cladelabel(756, "756", offset = 0.005)

tree<-ggtree(Cosmotree) +
  xlim(0, 0.1)

tree %<+% Cosmometa +
  geom_tippoint(aes(colour = (alignment.name)))

nodes_5<-nodes_5[order(nodes_5$`n tips`),]
nodes_5<-nodes_5[(which(nodes_5[,2]>=5)),]

m<-matrix(ncol=2, nrow=567)
tip_nodes<-data.frame(m)
names(tip_nodes)<-c("ID", "n_clades")
tip_nodes$ID<-Cosmotree$tip.label

for (i in (1:567)) {
  tip_nodes[i,2]<-
    length(
      intersect(
        Ancestors(Cosmotree, c(which(Cosmotree$tip.label == tip_nodes[i,1]))),
        nodes_5[,1]
      )
    )
}

tip_nodes<-tip_nodes[order(tip_nodes$n_clades),]

tip_nodes[which(tip_nodes$n_clades == 1),1]

intersect(Ancestors(Cosmotree, c(which(Cosmotree$tip.label ==  'KF154998'))), nodes_5[,1])

 # 576 (base), 1088, 1098