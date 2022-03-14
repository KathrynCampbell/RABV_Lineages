rm(list=ls())

tree<-ape::read.tree("Cosmo_N/Trees/Cosmo_N_outgroup_aligned.fasta.contree")
lineage_info<-read.csv("Cosmo_N/Outputs/Cosmo_N_updated_lineage_info.csv")
sequence_data<-read.csv("Cosmo_N/Outputs/Cosmo_N_sequence_data.csv")
metadata<-read.csv("Cosmo_WGS/Cosmo_WGS_metadata_outgroup.csv")
alignment<-seqinr::read.alignment("Cosmo_N/Alignment/Cosmo_N_aligned.fasta", format = "fasta")
ancestral <- seqinr::read.alignment(file = "Cosmo_N/Timetree/ancestral_sequences.fasta", format = "fasta")
ancestral$nam <- gsub("\\..*", "", ancestral$nam, perl = T)

alignment_matrix <- seqinr::as.matrix.alignment(alignment)
ancestral_matrix <- seqinr::as.matrix.alignment(ancestral)

tree$tip.label <- gsub("\\..*", "", tree$tip.label, perl = T)
tree$node.comment<- gsub(".*=", "", tree$node.label, perl = T)

nodes_70 <- which(tree$node.comment > 70 | tree$node.comment == 100)
nodes_70 <- nodes_70 + length(tree$tip.label)

node_data <- data.frame(Node = nodes_70, n_tips = NA)
# Make a dataframe ready for values to be put in
# Fill the first column with the numbers of the nodes identified in the previous steps

for(i in 1:length(nodes_70)) {
  node_data[i,2] <- length(phangorn::Descendants(tree, nodes_70[i], type = "tips")[[1]])
}
# For each node identified in the previous step, count the number of tips descended from that node

nodes_5 <- node_data[(which(node_data[,2]%in%5:9)),]

tests<-data.frame(lineage = lineage_info$lineage[which(lineage_info$n_N_seqs >= 5)], clades = NA)

for (x in 1:length(tests$lineage)) {
  MRCA<-getMRCA(tree, tip = (sequence_data$ID[which(sequence_data$cluster == tests$lineage[x])]))
  if (length(MRCA)!= 0) {
    descendents<-getDescendants(tree, node=MRCA)
    nodes<-which(descendents %in% nodes_5$Node)
    
    if(length(nodes) != 0){
      for (i in 1:length(nodes)) {
        clades<-unique(
          sequence_data$cluster[
            which(sequence_data$ID %in% tree$tip.label[getDescendants(tree, node=descendents[nodes[i]])])])
        if (length(clades)==1) {
          if (clades == tests$lineage[x]){
            tests$clades[x]<-paste(tests$clades[x], descendents[nodes[i]], sep = ",")
          }
        }
      }
    }
  }
}

tests<-tests[-c(which(is.na(tests$clades))),]

tests$clades<-gsub("NA,", "", tests$clades)

tests$n_clades<-NA

for (i in 1:length(tests$lineage)) {
  tests$n_clades[i]<-length((strsplit(tests$clades, ","))[[i]])
}

noi<-tests$clades[1]

for (i in 2:length(tests$lineage)) {
  noi<-paste(noi, tests$clades[i], sep = ",")
}

noi<-strsplit(noi, ",")[[1]]

noi<-data.frame(node = noi, lineage = NA, tips = NA)

for (i in 1:length(noi$node)) {
  noi$lineage[i]<-tests$lineage[grep(noi$node[i], tests$clades)]
}

noi$node<-as.integer(noi$node)

for (i in 1:length(noi$node)) {
  noi$tips[i]<-length(Descendants(tree, noi$node[i], "tips")[[1]])
}

noi$tips<-as.integer(noi$tips)

int<-data.frame(lineage = unique(noi$lineage), count = NA)

for (i in 1:length(int$lineage)) {
  int$count[i]<-length(which(noi$lineage == int$lineage[i]))
}

int$test<-NA

numbers<-which(int$count == 1)

for (i in 1:length(numbers)) {
  if ((lineage_info$n_N_seqs[which(lineage_info$lineage == int$lineage[numbers[i]])] -
       noi$tips[which(noi$lineage == int$lineage[numbers[i]])]) == 0) {
    int$test[numbers[i]]<-"N"
  }
}

nodes_remove <- phangorn::Ancestors(tree,
                                    (which(tree$tip.label
                                           %in% (sequence_data$ID[which(sequence_data$Length_after < (sequence_data$Length_before * 0.95))])
                                    )),
                                    'all')


if (length(nodes_remove)>0) {
  removes <- nodes_remove[[1]]
  for (i in 2:(length(nodes_remove))) {
    removes <- c(removes, nodes_remove[[i]])
  }
  remove_counts <- data.frame(table(removes))
  # Make a table to count the number the removed sequences descended from each node (e.g. for the deeper nodes, all 10 are descended)
  
  names(remove_counts) <-c('Node', 'freq')
  # Change the names
  
  remove_counts$Node <- as.integer(levels(remove_counts$Node))
  # Need to change this, or it creates many levels and causes errors
  
  new_remove <- remove_counts[which(remove_counts[,1] %in% noi[,1]),]; new_remove
  # Not all nodes are included in the nodes_5 data (some are already excluded)
  # Get rid of the nodes not in the nodes_5 data
  
  nodes_new<-noi
  
  for (i in 1:length(new_remove$Node)) {
    nodes_new$tips[which(nodes_new$node == new_remove$Node[i])]<-
      (nodes_new$tips[which(nodes_new$node == new_remove$Node[i])] - new_remove$freq[i])
  }
  # Take away the number of removed tips from the previous total number of tips calculated for each node
  
  noi <- nodes_new[(which(nodes_new$tips %in% 5:9)),] # Redo this to remove any that now have less than 5, and write over the old nodes_5 so this is updated with the new tip numbers
}

noi$diff <- NA
nodes_reduced <- data.frame(Nodes = (noi$node - (1+length(tree$tip.label))))

for (i in 1:length(noi$node)) {
  cm <- caper::clade.members(noi$node[i], tree, include.nodes = F, tip.labels = T)
  seq_cm <- which(sequence_data$ID %in% cm)
  old <- which(row.names(ancestral_matrix) == paste("NODE_", (sprintf("%07d", nodes_reduced$Nodes[i])), sep=""))
  
  tips <- which(row.names(ancestral_matrix) %in% cm)
  x <- which(ancestral_matrix[old,] != ancestral_matrix[(tips[1]),])
  
  for (j in tips[-c(1)]) {
    x <- x[which(x %in% (which(ancestral_matrix[old,] != ancestral_matrix[j,])))]
    print(x)
    noi$diff[i] <- length(x)
  }
}

nodes_diff <- noi[(which(noi$diff!=0)),] # Get rid of the ones with no differences straight away

#############################################
#         OVERLAPPING TIPS REMOVAL          #
#############################################
# Add a column to nodes_diff and for each node, count how many of the other nodes of interest are descended from it
nodes_diff$overlaps <- NA
for (i in 1:length(nodes_diff$node)) {
  nodes_diff$overlaps[i] <- length(which((phangorn::allDescendants(tree)[[(nodes_diff[i,1])]]) %in% nodes_diff[,1]))
}

nodes_diff<-nodes_diff[(which(nodes_diff$overlaps == 0)),]

nodes_diff$distance<-NA

for (i in 1:length(nodes_diff$node)) {
  parent<-Ancestors(tree, nodes_diff$node[i], "parent")
  nodes_diff$distance[i]<-castor::get_pairwise_distances(tree, nodes_diff$node[1], parent)
}

distances<-as.matrix(adephylo::distTips(tree, tips = "all", method = "patristic"))

nodes_diff$top<-NA

for (i in 1:length(nodes_diff$node)) {
  subset<-distances[which(rownames(distances)
                          %in% sequence_data$ID[which(sequence_data$cluster == nodes_diff$lineage[i])]),
                    which(colnames(distances)
                          %in% sequence_data$ID[which(sequence_data$cluster == nodes_diff$lineage[i])])]
  
  nodes_diff$top[i]<-quantile(subset, 0.95)
}

nodes_diff<-nodes_diff[which(nodes_diff$distance >= nodes_diff$top),]

sequence_data$Country<-gsub("-", NA, sequence_data$Country)

nodes_diff$country<-NA

for (i in 1:length(nodes_diff$node)) {
  countries<-unique(sequence_data$Country[
    which(sequence_data$ID %in% tree$tip.label[Descendants(tree, nodes_diff$node[i], "tips")[[1]]])])
  
  if(length(which(is.na(countries))) != 0) {
    countries<-countries[-c(which(is.na(countries)))]
  } 
  if (length(countries) != 1){
    nodes_diff$country[i]<-list(c(countries))
  } else {
    nodes_diff$country[i]<-countries
  }
}

nodes_diff$year_first<-NA
nodes_diff$year_last<-NA

for (i in 1:length(nodes_diff$node)) {
  nodes_diff$year_first[i]<-min(sequence_data$Year[
    which(sequence_data$ID %in% tree$tip.label[Descendants(tree, nodes_diff$node[i], "tips")[[1]]])])
  
  nodes_diff$year_last[i]<-max(sequence_data$Year[
    which(sequence_data$ID %in% tree$tip.label[Descendants(tree, nodes_diff$node[i], "tips")[[1]]])])
}

nodes_diff<-nodes_diff[which((nodes_diff$year_last - nodes_diff$year_first) < 5),]

if(length(grep(",", nodes_diff$country)) != 0) {
  manual<-nodes_diff[grep(",", nodes_diff$country)]
}

emerging_lineages<-nodes_diff[c(1:3, 6, 8:10)]

lengths<-data.frame(setNames(tree$edge.length[sapply(1:length(tree$tip.label),
                                                     function(x,y) which (y==x),y=tree$edge[,2])],
                             tree$tip.label))

colnames(lengths)<-"lengths"

seqs<-rownames(lengths)[which(lengths$lengths >= quantile(lengths$lengths, .95))]

seqs<-seqs[which(seqs %in% alignment$nam)]

seqs<-seqs[which(seqs %in% sequence_data$ID)]

lengths<-data.frame(ID = rownames(lengths), length = lengths$lengths)

lengths<-lengths[which(lengths$ID %in% seqs),]

lengths$close<-NA
lengths$lineage<-NA
lengths$cutoff<-NA

for (i in 1:length(lengths$ID)) {
  lengths$lineage[i]<-sequence_data$cluster[which(sequence_data$ID == lengths$ID[i])]
}

for (i in 1:length(lengths$ID)) {
  sequences<-sequence_data$ID[which(sequence_data$cluster == lengths$lineage[i])]
  subset<-distances[which(rownames(distances)%in% sequences),
                    which(colnames(distances)%in% sequences)]
  lengths$cutoff[i]<-quantile(subset, 0.95)
  
}

for (i in 1:length(lengths$ID)) {
  up<-alignment$nam[(which(alignment$nam == lengths$ID[i]))+1][[1]]
  down<-alignment$nam[(which(alignment$nam == lengths$ID[i]))-1][[1]]
  
  test<-which(tree$tip.label == lengths$ID[i])
  up<-which(tree$tip.label == up)
  down<-which(tree$tip.label == down)
  up<-castor::get_pairwise_distances(tree, test, up)
  down<-castor::get_pairwise_distances(tree, test, down)

  if(down > up && down >= lengths$cutoff[i]) {
    lengths$close[i]<-alignment$nam[(which(alignment$nam == lengths$ID[i]))-1]
  } else {
    if (up >= lengths$cutoff[i]) {
      lengths$close[i]<-alignment$nam[(which(alignment$nam == lengths$ID[i]))+1]
    } else {
      lengths$close[i] <- NA
    }
  }
}

lengths<-lengths[-c(which(is.na(lengths$close))),]

undersampled<-data.frame(lineage = lineage_info$lineage, n_singletons = NA, singleton_countries = NA, singleton_years = NA)

undersampled<-undersampled[which(undersampled$lineage %in% lengths$lineage),]

for (i in 1:length(undersampled$lineage)) {
  undersampled$n_singletons[i]<-length(which(lengths$lineage == undersampled$lineage[i]))
}

metadata<-read.csv("Cosmo_N/Cosmo_N_metadata.csv")

for (i in 1:length(lengths$lineage)) {
  lengths$year[i]<-metadata$year[which(metadata$ID == lengths$ID[i])]
  lengths$country[i]<-metadata$country[which(metadata$ID == lengths$ID[i])]
}

for (i in 1:length(undersampled$lineage)) {
  undersampled$singleton_countries[i]<-list(unique(lengths$country[which(lengths$lineage == undersampled$lineage[i])]))
  undersampled$singleton_years[i]<-list(unique(lengths$year[which(lengths$lineage == undersampled$lineage[i])]))
}

undersampled <- apply(undersampled,2,as.character)

write.csv(undersampled, "Emerg_under_singleton/singletons_of_interest.csv", row.names = F)

emerging_lineages$lineage<-paste(emerging_lineages$lineage, "E1", sep = "_")

emerging_lineages$lineage[
  which(duplicated(emerging_lineages$lineage))]<-
  gsub("_E1", "_E2", emerging_lineages$lineage[which(duplicated(emerging_lineages$lineage))])

emerging_lineages$lineage[
  which(duplicated(emerging_lineages$lineage))]<-
  gsub("_E2", "_E3", emerging_lineages$lineage[which(duplicated(emerging_lineages$lineage))])

emerging_lineages<-emerging_lineages[,2:7]

emerging_lineages <- apply(emerging_lineages,2,as.character)

write.csv(emerging_lineages, "Emerg_under_singleton/emerging_or_undersampled.csv", row.names = F)

